# md_verify_local.py — 本地复核三体系 MD 分析 CSV
# 严格级：TMEM45A（有落盘拓扑）全量重算对比；ENG/FLT1：DCD 头核验 + CSV 统计特征
import os, sys
import numpy as np
import pandas as pd
import mdtraj as mt

BASE = "E:/Zcode/0911"
MD = os.path.join(BASE, "results", "md_gpu")

def dcd_info(path):
    from mdtraj.formats.dcd import DCDTrajectoryFile
    with DCDTrajectoryFile(path) as f:
        xyz, cell, _ = f.read()
    return xyz.shape  # (frames, atoms, 3)

def csv_stats(path):
    d = pd.read_csv(path)
    out = {"rows": len(d), "nan": int(d.isna().sum().sum())}
    for c in d.columns:
        if c == "time_ns":
            continue
        v = d[c].astype(float)
        out[c] = (round(float(v.iloc[0]), 3), round(float(v.iloc[-1]), 3),
                  round(float(v.mean()), 3), round(float(v.std()), 3))
    return out

def verify_light(gene):
    dcd = os.path.join(MD, gene + "_PFOS_prod.dcd")
    frm, nat, _ = dcd_info(dcd)
    m = csv_stats(os.path.join(MD, gene + "_PFOS_md_metrics.csv"))
    c = csv_stats(os.path.join(MD, gene + "_PFOS_md_contacts.csv"))
    print("[%s] DCD: %d frames x %d atoms" % (gene, frm, nat))
    print("  metrics rows=%d nan=%d" % (m["rows"], m["nan"]))
    for k, v in m.items():
        if k not in ("rows", "nan"):
            print("    %-16s first/last/mean/std = %s" % (k, v))
    print("  contacts rows=%d nan=%d" % (c["rows"], c["nan"]))
    for k, v in c.items():
        if k not in ("rows", "nan"):
            print("    %-24s first/last/mean/std = %s" % (k, v))
    ok = (m["rows"] == 500 and c["rows"] == 500 and m["nan"] == 0 and c["nan"] == 0)
    print("  [%s] light verify: %s" % (gene, "PASS" if ok else "FAIL"))
    return ok

def verify_tmem45a_full():
    dcd = os.path.join(MD, "TMEM45A_PFOS_prod.dcd")
    top = os.path.join(MD, "last_system_topology.pdb")
    mdtop = mt.Topology.from_openmm(PDBFile(top).topology) if False else mt.load_pdb(top).topology
    t = mt.load(dcd, top=mdtop)
    print("[TMEM45A-full] frames=%d atoms=%d" % (t.n_frames, t.n_atoms))
    protein = t.topology.select("protein")
    lig = t.topology.select("resname LIG")
    rmsd_prot = mt.rmsd(t, t, 0, atom_indices=protein) * 10
    t.superpose(t, 0, atom_indices=protein)
    rmsd_lig = mt.rmsd(t, t, 0, atom_indices=lig) * 10
    rg = mt.compute_rg(t.atom_slice(protein)) * 10
    ref = pd.read_csv(os.path.join(MD, "TMEM45A_PFOS_md_metrics.csv"))
    diffs = {
        "rmsd_protein_A": float(np.max(np.abs(rmsd_prot - ref["rmsd_protein_A"].values))),
        "rmsd_ligand_A": float(np.max(np.abs(rmsd_lig - ref["rmsd_ligand_A"].values))),
        "rg_protein_A": float(np.max(np.abs(rg - ref["rg_protein_A"].values))),
    }
    for k, v in diffs.items():
        print("    max|delta| %-16s = %.4f A" % (k, v))
    ok = all(v < 0.05 for v in diffs.values()) and len(ref) == t.n_frames
    print("  [TMEM45A-full] strict verify: %s" % ("PASS" if ok else "FAIL"))
    return ok

if __name__ == "__main__":
    from openmm.app import PDBFile  # noqa (unused alternative)
    ok = True
    ok &= verify_light("ENG")
    ok &= verify_light("FLT1")
    ok &= verify_light("TMEM45A")
    ok &= verify_tmem45a_full()
    print("OVERALL:", "PASS" if ok else "FAIL")
