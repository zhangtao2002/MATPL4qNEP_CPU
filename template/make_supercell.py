#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""由小的 LAMMPS data（atom_style charge）沿 a/b/c 复制成超胞。

用法: python3 make_supercell.py data.lmp data_big.lmp 4 4 2
（4 4 2 = 在 a/b/c 方向各复制几次）
"""
import sys


def main():
    if len(sys.argv) != 6:
        raise SystemExit("用法: make_supercell.py <in.data> <out.data> <n1> <n2> <n3>")
    src, dst = sys.argv[1], sys.argv[2]
    n1, n2, n3 = (int(x) for x in sys.argv[3:6])

    L = open(src).read().splitlines()
    box = {}
    start = None
    for i, ln in enumerate(L):
        s = ln.split()
        if ln.strip().endswith("xlo xhi"):
            box["x"] = (float(s[0]), float(s[1]))
        elif ln.strip().endswith("ylo yhi"):
            box["y"] = (float(s[0]), float(s[1]))
        elif ln.strip().endswith("zlo zhi"):
            box["z"] = (float(s[0]), float(s[1]))
        elif ln.strip().endswith("xy xz yz"):
            box["tilt"] = tuple(float(v) for v in s[:3])
        elif ln.strip().startswith("Atoms"):
            start = i + 2
    lx = box["x"][1] - box["x"][0]
    ly = box["y"][1] - box["y"][0]
    lz = box["z"][1] - box["z"][0]
    xy, xz, yz = box.get("tilt", (0.0, 0.0, 0.0))

    atoms = []
    for ln in L[start:]:
        w = ln.split()
        if len(w) != 6:
            break
        atoms.append((int(w[0]), int(w[1]), float(w[2]), float(w[3]), float(w[4]), float(w[5])))

    out = []
    for a in range(n1):
        for b in range(n2):
            for c in range(n3):
                dx, dy, dz = a * lx, b * ly, c * lz
                for (idx, tp, q, x, y, z) in atoms:
                    out.append((tp, q, x + dx + (xy / ly) * dy + (xz / lz) * dz,
                                y + dy + (yz / lz) * dz, z + dz))

    with open(dst, "w") as f:
        f.write("supercell %dx%dx%d\n\n%d atoms\n2 atom types\n\n" % (n1, n2, n3, len(out)))
        f.write("0.0 %.8f xlo xhi\n0.0 %.8f ylo yhi\n0.0 %.8f zlo zhi\n" %
                (n1 * lx, n2 * ly, n3 * lz))
        f.write("%.8f %.8f %.8f xy xz yz\n\nMasses\n\n1 15.9994\n2 1.008\n\nAtoms # charge\n\n" %
                (xy, xz, yz))
        for k, (tp, q, x, y, z) in enumerate(out):
            f.write("%d %d %.6f %.8f %.8f %.8f\n" % (k + 1, tp, q, x, y, z))
    print("wrote %s: %d atoms" % (dst, len(out)))


if __name__ == "__main__":
    main()
