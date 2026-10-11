#!/bin/bash
#SBATCH -N 1
#SBATCH -n 96
#SBATCH --ntasks-per-node=96
#SBATCH --partition=nudt96c          # 或 nudt40c (40 核/节点)
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err
#SBATCH --job-name=nepmd
#
# 用法:  sbatch submit.sh
# 纯 CPU, 全部 MPI 并行 (描述符/力/Ewald)。OpenMP 无效, 别开。
# 体系建议 >= 几万原子才近线性扩展 (见 README)。
#
# 注意: 这个 LAMMPS 是 OpenMPI 版, 必须用 mpirun (srun 会退化成 N 个独立进程!)。
#   多节点: -N>1 时照用 mpirun, SLURM 会提供 hostfile。

source ~/deepmd-kit/bin/activate
LMP=~/build_nep_cpu/lammps/build/lmp

mpirun -np "$SLURM_NTASKS" "$LMP" -in in.lammps
