#!/bin/bash
# ============================================================================
# 一键构建「纯 CPU 的 qnep-MatPL LAMMPS」(pair_style matpl/nep, 含长程能)
#
#   LAMMPS_SRC=/path/to/lammps-stable_29Aug2024 bash build/build_nep_cpu.sh
#
# 需要 : gcc/g++ >= 8, MPI(mpirun/mpicxx), cmake >= 3.21, FFTW3(LAMMPS 自带/系统)
# 不需要: CUDA / libtorch
# 产物 : $BUILD_ROOT/lammps/build/lmp   (默认 $PWD/build_nep_cpu)
# ============================================================================
set -e

REPO="$(cd "$(dirname "$0")/.." && pwd)"
LAMMPS_SRC="${LAMMPS_SRC:?请设置 LAMMPS_SRC=干净的 LAMMPS 源码目录}"
BUILD_ROOT="${BUILD_ROOT:-$PWD/build_nep_cpu}"
JOBS="${JOBS:-$(nproc)}"

echo ">> [1/6] 复制 LAMMPS 源 -> $BUILD_ROOT/lammps"
mkdir -p "$BUILD_ROOT"
if command -v rsync >/dev/null 2>&1; then
    rsync -a --exclude 'src/Obj_*' --exclude 'build' "$LAMMPS_SRC/" "$BUILD_ROOT/lammps/"
else
    cp -a "$LAMMPS_SRC" "$BUILD_ROOT/lammps"
    rm -rf "$BUILD_ROOT"/lammps/src/Obj_*
fi
cd "$BUILD_ROOT/lammps/src"

echo ">> [2/6] 清掉会冲突的旧 USER-NEP + make 构建残留"
# USER-NEP 的类名也叫 PairNEP, 会和 MatPL 的 pair_nep.cpp 撞名; 必须删
rm -f pair_NEP.cpp pair_NEP.h nep.cpp nep.h \
      neighbor_nep.cpp neighbor_nep.h ewald_nep.cpp ewald_nep.h nep_utilities.h
rm -rf USER-NEP
make no-all purge >/dev/null 2>&1 || true
rm -f style_*.h lmpinstalledpkgs.h lmpgitversion.h
# ⚠ 千万不要删 lmptype.h / lmpfftsettings.h —— 它们是真实源文件!

echo ">> [3/6] 拷入 MatPL CPU pair style 源码 (nep_cpu.* + pair_nep.*)"
cp "$REPO"/matpl_cpu_src/nep_cpu.cpp "$REPO"/matpl_cpu_src/nep_cpu.h \
   "$REPO"/matpl_cpu_src/pair_nep.cpp "$REPO"/matpl_cpu_src/pair_nep.h .

echo ">> [3b/6] 拷入自定义 fix spring/self/scale (k 可读变量, 用于 TI 连续缩放弹簧)"
cp "$REPO"/fixes/fix_spring_self_scale.cpp "$REPO"/fixes/fix_spring_self_scale.h .

echo ">> [4/6] 打补丁: 长程(Ewald)能进 pe, 并单独输出到 ecoul"
python3 "$REPO"/build/patch_pe.py  "$PWD/nep_cpu.cpp"
python3 "$REPO"/build/patch_pe2.py "$PWD"

echo ">> [5/6] cmake 配置 (basic.cmake 自带 KSPACE/FFTW3, 纯 CPU)"
cd "$BUILD_ROOT/lammps"
mkdir -p build && cd build
cmake -C ../cmake/presets/basic.cmake ../cmake

echo ">> [6/6] 编译 (-j$JOBS)"
cmake --build . -j "$JOBS"

echo
echo "=========================================================="
echo " 完成: $BUILD_ROOT/lammps/build/lmp"
echo " 用法: pair_style matpl/nep nep.txt kspace ewald"
echo "       pair_coeff * * <元素列表>"
echo "       (运行时必须用 mpirun, 不要用 srun)"
echo "=========================================================="
