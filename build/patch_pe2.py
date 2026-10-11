#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""让长程(Ewald)能单独进入 LAMMPS 的 eng_coul，便于用 thermo 打印 (ecoul)。

前提：已应用 patch_pe.py（倒空间能已算好，现加到 total_potential）。
本补丁把它改到 ewald_energy 成员 -> pair 里 eng_coul += ewald_energy。
于是:  pe = eng_vdwl(短程) + eng_coul(长程) 仍正确,  且 thermo 可单独打 ecoul。
"""
import sys

SRC = sys.argv[1] if len(sys.argv) > 1 else "/data/home/nudtzhangtao/build_nep_cpu/lammps/src"
cpp = SRC + "/nep_cpu.cpp"
h = SRC + "/nep_cpu.h"
pair = SRC + "/pair_nep.cpp"

# --- nep_cpu.h: 加成员 ---
s = open(h).read()
assert "double charge_alpha_factor" in s
if "ewald_energy" not in s:
    s = s.replace("    double charge_alpha_factor = 0.0;\n",
                  "    double charge_alpha_factor = 0.0;\n"
                  "    double ewald_energy = 0.0;   // 最近一次算出的长程(Ewald)能\n", 1)
    open(h, "w").write(s)
    print("h: + ewald_energy")

# --- nep_cpu.cpp ---
s = open(cpp).read()
old_sig = ("  int virial_components,\n  double& total_potential,\n  const int model_index)\n"
           "{\n  zero_global_mean(nlocal, natoms_global, mpi_comm, charge);")
new_sig = ("  int virial_components,\n  double& total_ewald,\n  const int model_index)\n"
           "{\n  zero_global_mean(nlocal, natoms_global, mpi_comm, charge);")
assert old_sig in s, "ewald 签名没找到"
s = s.replace(old_sig, new_sig, 1)
assert "    total_potential += e_recip;" in s, "e_recip 累加行没找到"
s = s.replace("    total_potential += e_recip;", "    total_ewald += e_recip;", 1)
# compute_for_lammps: 每次清零成员 + 传成员给 ewald
assert "    num_atoms = nlocal;\n  }\n  find_descriptor_for_lammps(" in s
s = s.replace("    num_atoms = nlocal;\n  }\n  find_descriptor_for_lammps(",
              "    num_atoms = nlocal;\n  }\n  ewald_energy = 0.0;\n  find_descriptor_for_lammps(", 1)
assert "        virial_components, total_potential, model_index);" in s
s = s.replace("        virial_components, total_potential, model_index);",
              "        virial_components, ewald_energy, model_index);", 1)
open(cpp, "w").write(s)
print("cpp: patched")

# --- pair_nep.cpp: eng_coul ---
s = open(pair).read()
a = ("        if (eflag) {\n"
     "            eng_vdwl += total_potential;\n"
     "        }")
b = ("        if (eflag) {\n"
     "            eng_vdwl += total_potential;\n"
     "            eng_coul += nep_cpu_models[0].ewald_energy;   // 长程(Ewald)能\n"
     "        }")
assert a in s, "单模型 eflag 段没找到"
s = s.replace(a, b, 1)
c = ("                if (eflag) {\n"
     "                    eng_vdwl = total_potential;\n"
     "                }")
d = ("                if (eflag) {\n"
     "                    eng_vdwl = total_potential;\n"
     "                    eng_coul = nep_cpu_models[0].ewald_energy;   // 长程(Ewald)能\n"
     "                }")
assert c in s, "多模型 eflag 段没找到"
s = s.replace(c, d, 1)
open(pair, "w").write(s)
print("pair: patched eng_coul")
