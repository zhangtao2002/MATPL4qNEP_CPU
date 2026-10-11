# qnep (MatPL) 纯 CPU 模拟模板

用 LAMMPS + `pair_style matpl/nep`（CPU-only 版）跑你的 `nep4_charge2` 势函数。
描述符、力、Ewald/电荷**全部 MPI 并行**，大体系近线性扩展。

## 1. 环境

```bash
source ~/deepmd-kit/bin/activate                 # 提供 OpenMPI + gcc13
export LMP=~/build_nep_cpu/lammps/build/lmp      # 纯 CPU 的 qnep-MatPL LAMMPS
```

这个 `lmp` 是 2026-10 在 hpc3 上编的：MatPL-2026.3 的 `lmp_nepkokkos_cmake`
打 `kknep-patch-cpu.sh` 补丁 + `cmake -C basic.cmake`，无 CUDA、无 torch。
源码树 `~/build_nep_cpu/lammps/`。

## 2. 目录内容

| 文件 | 说明 |
|---|---|
| `nep.txt` | 势函数（`nep4_charge2`） |
| `in.lammps` | LAMMPS 输入模板（改温度/步数/结构名即可） |
| `data.lmp` | 示例结构（324 原子 O+H，`atom_style charge`） |
| `submit.sh` | SLURM 提交模板 |
| `make_supercell.py` | 由小胞造超胞（可选） |

## 3. 关键点（照抄）

- **`atom_style charge` 必须**（电荷模型）。
- **`pair_style matpl/nep nep.txt kspace ewald`** + `pair_coeff * * O H`。
- 长程用 **`kspace ewald`**：你的 `nep4_charge2` 就是**纯 Ewald**，比 `pppm`
  **快 ~1.4×**，结果一致（能量相同、压强差 ~0.4%）。别再写 `kspace pppm`。
  * 也**不要**另外再写 `kspace_style pppm`（会报 KSpace/Pair 不兼容）；长程由
    pair style 的 `kspace ewald` 参数接管。
- **只有 MPI 有用，OpenMP 无效**（实测）。用 `srun`/`mpirun -np N` 铺核。
- **体系别太小**：10k 原子以下 ghst/通信占比高、扩展差；几万原子起才近线性。

## 4. 性能参考（单节点 nudt96c，46,656 原子）

| np | timesteps/s | 加速比 |
|---|---|---|
| 1 | 0.081 | 1.0× |
| 16 | 1.12 | 13.8× |
| 40 | 2.48 | 30.6× |
| 96 | 5.74 | 70.8× |

（0.5 fs 步长 → 96 核 ~0.17 s/步，3 ps ≈ 17 分钟）

## 5. 跑

```bash
sbatch submit.sh          # 用模板里的分区/核数
# 或交互测试：
mpirun -np 16 $LMP -in in.lammps
```

## 6. 换自己的体系

1. 准备 `data.lmp`（`atom_style charge`，Atoms 段格式 `id type q x y z`，q 先填 0）。
   可由 `atomsk` 造 O 晶格 + `add_h_to_both.py` 加 H 后转成 LAMMPS data。
2. `pair_coeff * * <元素列表>` 的顺序要和 `nep.txt` 头一致（这里是 `O H`）。
3. 改 `in.lammps` 里的 `TEBEG`/温度、`run` 步数、`read_data` 名字。
