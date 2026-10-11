# MATPL4qNEP_CPU

纯 CPU 版 **MatPL qnep**（LAMMPS `pair_style matpl/nep`），用于在**无 GPU** 集群上运行
`nep4_charge2` 势函数，长程走 LAMMPS Ewald，**MPI 并行（大体系近线性）**，并**修正了总能量（含长程）**，
可用于 TI / 自由能等依赖能量绝对值的计算。

---

## 1. 背景：上游 CPU 版有什么问题

上游 [MatPL-2026.3](https://github.com/...) 提供了 CPU 的 NEP pair style（`matpl/nep` + 倒空间 Ewald），
能跑 qnep，但有两处对"能量"不利的地方：

1. **`pe`（总能）漏了长程 Ewald 能**——代码里 `eng_vdwl += total_potential`，而 `total_potential`
   只由描述符（短程 NEP）累加；倒空间 Ewald 那段**只算力/virial，没算能量**。
   做自由能/TI 时势能绝对值会错（实测比 calorine 少 ~4 eV/胞，且随构型变）。
2. 长程能**无法单独打印**。

本仓库给出：CPU pair style 源码 + **两处补丁** + **一键构建脚本** + 运行模板。
修复后与 **calorine**（同一个 qnep 势的另一套 CPU 实现）**逐位一致**。

## 2. 改动细节（相对上游 MatPL CPU 源码）

### 补丁 1 — 把长程(Ewald)能计入总能 `pe`（`build/patch_pe.py`）
对照 GPUMD `src/main_nep/nep_charge.cu`：`nep4_charge2` 是 **charge_mode==2（纯倒空间 Ewald）**，
其能核为 `g_pe[n] += K_C · Σ_k G(k)|S(k)|² / N`（self/real-space 只给 charge_mode==1）。
CPU 版 `D_real` 已编码结构因子（`Σ_n q_n·D_real[n] = 2·K_C·Σ_k G(k)|S(k)|²`），故：

- `nep_cpu.cpp::find_force_charge_ewald_cpu` 增加形参 `double& total_ewald`；
  末尾插入 `total_ewald += 0.5 * Σ_n charge[n]·D_real[n];`
- 调用点补上实参。
- **不要**额外加 self 项（`-0.5·(2α/√π)Σq²`）——那是 charge1 的，加了会过冲 ~20 eV。

### 补丁 2 — 长程能单独输出到 `ecoul`（`build/patch_pe2.py`）
- `nep_cpu.h`：`NEP3_CPU` 加成员 `double ewald_energy`；
- ewald 把倒空间能累加到该成员；`pair_nep.cpp` 里 `eng_coul += nep_cpu_models[0].ewald_energy;`。
- 结果：`pe = eng_vdwl(短程) + eng_coul(长程)` 仍正确，且可分别打印。

> 实测（324 原子，VASP 轨迹一帧）：`pe=-915.59464`，`evdwl=-918.55784`（短程），
> `ecoul=+2.9631925`（长程）。

细节见 [`BUILD_NOTES.md`](BUILD_NOTES.md)。

## 3. 编译

### 依赖
| 需要 | 不需要 |
|---|---|
| gcc/g++ ≥ 8（带 OpenMP 更好，但 OpenMP 对速度无帮助） | **CUDA** |
| MPI（`mpirun` / `mpicxx`） | **libtorch** |
| CMake ≥ 3.21 | |
| FFTW3（LAMMPS 自带或系统） | |
| 一份**干净**的 LAMMPS 源码（`stable_29Aug2024` 已验证） | |

> 只依赖 LAMMPS 源码；MatPL 的 CPU pair style 源码已随本仓库提供（`matpl_cpu_src/`）。

### 一键构建
```bash
LAMMPS_SRC=/path/to/lammps-stable_29Aug2024 bash build/build_nep_cpu.sh
# 产物: ./build_nep_cpu/lammps/build/lmp
```

### 手动构建（等价流程）
```bash
# 1. 复制干净 LAMMPS
rsync -a --exclude 'src/Obj_*' --exclude 'build' $LAMMPS_SRC/ ./build_nep_cpu/lammps/
cd build_nep_cpu/lammps/src
# 2. 删掉会撞名的旧 USER-NEP (类名 PairNEP) + make 残留
rm -f pair_NEP.cpp pair_NEP.h nep.cpp nep.h neighbor_nep.* ewald_nep.* nep_utilities.h
rm -rf USER-NEP && make no-all purge    # 别删 lmptype.h / lmpfftsettings.h!
# 3. 拷入 CPU pair style
cp <repo>/matpl_cpu_src/nep_cpu.cpp  <repo>/matpl_cpu_src/nep_cpu.h  .
cp <repo>/matpl_cpu_src/pair_nep.cpp <repo>/matpl_cpu_src/pair_nep.h .
# 4. 打补丁
python3 <repo>/build/patch_pe.py  ./nep_cpu.cpp
python3 <repo>/build/patch_pe2.py .
# 5. 编译
cd .. && mkdir build && cd build
cmake -C ../cmake/presets/basic.cmake ../cmake
cmake --build . -j 32
```

## 4. 运行

```
units        metal
atom_style   charge            # 电荷模型必须
pair_style   matpl/nep nep.txt kspace ewald
pair_coeff   * * O H           # 顺序与 nep.txt 头一致
```

要点：
- **必须 `mpirun -np N lmp`**（**不要 `srun`**：这份 LAMMPS 是 OpenMPI 版，srun 会退化成 N 个互相独立的进程）。
- 长程用 **`kspace ewald`**（`nep4_charge2` 是纯 Ewald；比 `pppm` 快且对口）；**不要再写 `kspace_style`**。
- 能量：`thermo_style custom step temp pe evdwl ecoul ...`
  → `pe`=总能，`evdwl`=短程，`ecoul`=长程(Ewald)。
- 纯 MPI 并行，OpenMP 无加速（别设 `OMP_NUM_THREADS`）。

运行模板见 [`template/`](template/)（`in.lammps` / `submit.sh` / `make_supercell.py`）。

## 5. 附加: 自定义 fix `spring/self/scale`（TI 连续缩放弹簧用）

做自由能 / TI 时常需要把"拴到初始位置"的弹簧按 `(1-λ)` **连续缩放**。LAMMPS 的
`fix spring/self` 的 `k` **只能是常数**；`fix adapt` 不支持缩放 fix（实测报
`Nothing to adapt`）；`fix ti/spring` 的参考态又只有弹簧，无法保留额外的软势（UF）。
因此本仓库附带一个很小的自定义 fix：

```
fix ID group spring/self/scale k_or_varname [x|y|z|xy|xz|yz|xyz]
```

- 与 `fix spring/self` 完全一致，只是 `arg[3]` 可以是 **equal-style 变量**（写成 `v_<name>`），
  **每步重新求值**，于是 `k` 可随 `λ` 连续变化：
  ```
  variable  kappa_2 equal 10.6889*(1.0-v_lambda)
  fix       f_spring O spring/self/scale v_kappa_2
  fix_modify f_spring energy no
  ```
- 配合 `pair_style hybrid/scaled` 可整体缩放多体势（NEP 无需支持 `extract`）：
  ```
  pair_style hybrid/scaled v_lambda matpl/nep nep.txt kspace ewald v_one_minus_lambda ufm 5.0
  pair_coeff * * matpl/nep O H
  ```
- **验证**：`λ=0.5`（`k=5.34445`）变量版与 `fix spring/self 5.34445` 常数版
  200 步**逐行完全一致**；`λ=1` 时 `f_spring=0`，`λ=0` 时弹簧全开。
- 源文件：[`fixes/fix_spring_self_scale.{cpp,h}`](fixes/)（基于 LAMMPS `fix_spring_self`，
  改动仅限"k 读变量"）。

## 6. 性能（单节点，`kspace ewald`）

| 体系 | 配置 | 速度 |
|---|---|---|
| 10,368 原子 | 1 核 | 0.36 步/s（2.7 s/步） |
| 10,368 原子 | 40 核 | **11.1 步/s**（0.09 s/步，3 ps ≈ 9 min） |
| 46,656 原子 | 96 核 | 5.7 步/s（**70.8×**，效率 74%） |

体系统大 MPI 扩展越好（10k 以下太小的会因 ghost/通信占比高而变差）。
补长程能对速度**无可见影响**（±1% 噪声）。

## 7. 许可 / 致谢
- NEP 势函数与描述符的实现源自 **GPUMD**（GPL）；LAMMPS 与 **MatPL** 亦为 GPL。
  本仓库是它们的 CPU 派生物与补丁，**遵循上游 GPL**，版权归原作者。
- 交叉验证使用了 **calorine**（GPUMD 的 Python 接口）。

## 目录结构
```
README.md             本文件
BUILD_NOTES.md        详细构建记录 + 改动细节（含踩坑）
build/
  build_nep_cpu.sh    一键构建
  patch_pe.py         补丁1: 长程能 -> pe
  patch_pe2.py        补丁2: 长程能 -> ecoul
matpl_cpu_src/        上游 MatPL-2026.3 的 CPU pair style 源码（未改）
patched/              打完两处补丁的最终文件（可直接替换/参考）
fixes/                自定义 fix spring/self/scale（TI 连续缩放弹簧用）
template/             运行模板 (in.lammps / submit.sh / make_supercell.py / README)
```
