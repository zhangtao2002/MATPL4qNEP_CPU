#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""给 MatPL 的 ewald 路径补上长程(Ewald)能量，使其进入 LAMMPS 的 pe。

改动 nep_cpu.cpp:
  1. find_force_charge_ewald_cpu 增加 double& total_potential 形参;
  2. 函数末尾把 E_recip=(1/2)Σ q_n D_real[n]  (= K_C Σ_k G(k)|S(k)|²) 和
     E_self = -0.5 * (2α/√π) Σ q_n²  累加进 total_potential;
  3. 调用点补上 total_potential 实参。
（D_real 已编码倒空间结构因子，能量可由它直接得到。）
"""
import re
import sys

p = sys.argv[1] if len(sys.argv) > 1 else \
    "/data/home/nudtzhangtao/build_nep_cpu/lammps/src/nep_cpu.cpp"
src = open(p).read()
assert "#include" in src

# 1. 签名加 total_potential
m = re.search(
    r"(void find_force_charge_ewald_cpu\(.*?int virial_components,\n)(  const int model_index\))",
    src, re.S)
assert m, "找不到 ewald 签名"
src = src[:m.end(1)] + "  double& total_potential,\n" + src[m.end(1):]

# 2. 函数末尾累加能量 (插到该函数里最后那句 zero_global_mean(D_real) 之前)
i0 = src.index("void find_force_charge_ewald_cpu(")
k = src.index("zero_global_mean(nlocal, natoms_global, mpi_comm, D_real);", i0)
energy = (
    "  {\n"
    "    // qnep 的 charge2 模型: 只有倒空间 Ewald 能 (GPUMD main_nep/nep_charge.cu:\n"
    "    //   g_pe[n] += K_C * sum_k G(k)|S(k)|^2 / N ; self/real-space 只给 charge1)\n"
    "    // D_real[n] 已编码倒空间结构因子: sum_n q_n D_real[n] = 2 K_C sum_k G|S|^2\n"
    "    double e_recip = 0.0;\n"
    "    for (int n = 0; n < nlocal; ++n) e_recip += 0.5 * charge[n] * D_real[n];\n"
    "    total_potential += e_recip;\n"
    "  }\n")
src = src[:k] + energy + src[k:]

# 3. 调用点补实参
old_call = ("      find_force_charge_ewald_cpu(\n"
            "        paramb, nlocal, natoms_global, pos, box, mpi_comm, charge, D_real, force, total_virial, virial,\n"
            "        virial_components, model_index);")
new_call = ("      find_force_charge_ewald_cpu(\n"
            "        paramb, nlocal, natoms_global, pos, box, mpi_comm, charge, D_real, force, total_virial, virial,\n"
            "        virial_components, total_potential, model_index);")
assert old_call in src, "找不到调用点"
src = src.replace(old_call, new_call)

open(p, "w").write(src)
print("patched:", p)
