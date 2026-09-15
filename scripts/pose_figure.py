# pose_figure.py — 对接位姿 3D 图（matplotlib）：蛋白 Cα 灰线 + 近距残基棒 + PFOS 球棍
# 用法: python pose_figure.py TMEM45A ; python pose_figure.py FLT1
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3D

BASE = "E:/Zcode/0911"

def parse_pdb(path):
    atoms = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            atoms.append({
                "name": line[12:16].strip(),
                "res": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "resi": int(line[22:26]),
                "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                "elem": line[76:78].strip(),
            })
    return atoms

def render(gene):
    rec = parse_pdb(BASE + "/data/receptors/%s_clean.pdb" % gene)
    lig = parse_pdb(BASE + "/results/md/%s_PFOS_ligand_openmm.pdb" % gene)
    rec_xyz = np.array([a["xyz"] for a in rec])
    lig_xyz = np.array([a["xyz"] for a in lig])
    d = np.linalg.norm(rec_xyz[:, None, :] - lig_xyz[None, :, :], axis=2)
    near = np.where(d.min(axis=1) < 4.5)[0]

    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection="3d")

    # 蛋白：Cα 骨架灰色走线（按残基序取首个原子）
    seen, ca_idx = set(), []
    for i, a in enumerate(rec):
        if a["resi"] not in seen:
            seen.add(a["resi"]); ca_idx.append(i)
    ca = rec_xyz[ca_idx]
    ax.plot(ca[:, 0], ca[:, 1], ca[:, 2], color="0.78", lw=1.0, alpha=0.9, zorder=1)

    # 近距残基（<4.5 Å）：侧链棒状橙色
    for i in near:
        a = rec[i]
        ax.plot(*zip(rec_xyz[i], rec_xyz[i]), color="none")  # 占位防 zorder 退化
        ax.scatter(*rec_xyz[i], s=14, c="#f4a582", depthshade=False, zorder=3)

    # 配体：PFOS 球棍（C 灰、F 绿、O 红、S 黄）
    elem_color = {"C": "#4d4d4d", "F": "#7fc97f", "O": "#e41a1c", "S": "#ffdf33"}
    for i, a in enumerate(lig):
        c = elem_color.get(a["elem"], "#999999")
        ax.scatter(*lig_xyz[i], s=120, c=c, edgecolors="k", linewidths=0.5,
                   depthshade=False, zorder=5)
    # 配体内连接（距离 <2.0 Å 的重原子对）
    n = len(lig_xyz)
    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(lig_xyz[i] - lig_xyz[j]) < 2.0:
                ax.plot(*zip(lig_xyz[i], lig_xyz[j]), color="#4d4d4d", lw=2.0, zorder=4)
    # 配体-蛋白近距接触虚线
    for j in near:
        i = int(np.argmin(d[j]))
        if d[j, i] < 3.6:
            ax.plot(*zip(lig_xyz[i], rec_xyz[j]), color="#2b6cb0", lw=0.8, ls="--", zorder=2)

    ax.set_title("%s – PFOS docked pose (best Vina mode)" % gene, fontsize=13)
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)"); ax.set_zlabel("z (Å)")
    ax.view_init(elev=18, azim=35)
    from matplotlib.lines import Line2D
    legend_elems = [
        Line2D([0], [0], color="0.78", lw=2, label="Protein Cα trace"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#f4a582", markersize=8, label="Contact residues (<4.5 Å)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#7fc97f", markersize=8, label="PFOS (F)"),
        Line2D([0], [0], color="#2b6cb0", ls="--", lw=1, label="Ligand–protein contact"),
    ]
    ax.legend(handles=legend_elems, loc="upper left", fontsize=9)
    plt.tight_layout()
    out = BASE + "/results/figures/fig_dock_pose_%s.png" % gene
    plt.savefig(out, dpi=300)
    plt.close()
    print("saved", out, "| near residues:", len(near))

if __name__ == "__main__":
    for g in (sys.argv[1:] or ["TMEM45A", "FLT1"]):
        render(g)
