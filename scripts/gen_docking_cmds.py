# gen_docking_cmds.py — 只读受体坐标、只打印 Vina 命令（不做任何文件写入）
# 盒子：坐标范围 + 6 A 填充，单边上限 100 A（超限以质心为中心截取）
import glob, os

REC_DIR = "data/receptors"
LIGS = {"PFOA": "data/compound/PFOA.pdbqt", "PFOS": "data/compound/PFOS.pdbqt"}
PAD, CAP, EXH, CPU = 6.0, 100.0, 16, 6

def bbox(pdb):
    xs, ys, zs = [], [], []
    with open(pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
    return min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)

for rec in sorted(glob.glob(os.path.join(REC_DIR, "*_receptor.pdbqt"))):
    rec = rec.replace(os.sep, "/")  # bash 消费，必须正斜杠
    gene = os.path.basename(rec).replace("_receptor.pdbqt", "")
    clean = REC_DIR + "/" + gene + "_clean.pdb"
    x0, x1, y0, y1, z0, z1 = bbox(clean)
    cx, cy, cz = (x0+x1)/2, (y0+y1)/2, (z0+z1)/2
    sx = round(min(x1-x0 + 2*PAD, CAP), 2)
    sy = round(min(y1-y0 + 2*PAD, CAP), 2)
    sz = round(min(z1-z0 + 2*PAD, CAP), 2)
    for lig_name, lig in LIGS.items():
        lig = lig  # 字面量路径，已是正斜杠
        out = "results/docking/%s_%s_out.pdbqt" % (gene, lig_name)
        log = "results/docking/%s_%s.vina.log" % (gene, lig_name)
        if os.path.exists(out) and os.path.getsize(out) > 0:
            print('echo "[skip] %s-%s exists" >&2' % (gene, lig_name))
            continue
        print('echo "[run] %s-%s" >&2' % (gene, lig_name))
        print('tools/vina_1.2.7_win.exe --receptor %s --ligand %s '
              '--center_x %.3f --center_y %.3f --center_z %.3f '
              '--size_x %s --size_y %s --size_z %s '
              '--exhaustiveness %d --cpu %d --num_modes 9 '
              '--out %s > %s 2>&1' % (rec, lig, cx, cy, cz, sx, sy, sz,
                                      EXH, CPU, out, log))
