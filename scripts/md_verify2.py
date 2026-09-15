# md_verify2.py — 补查：FLT1/TMEM45A CSV 统计 + TMEM45A 构建态配体-蛋白最小距离（拓扑PDB内含坐标）
import os
import numpy as np
import pandas as pd
import mdtraj as mt

BASE = "E:/Zcode/0911"
MD = os.path.join(BASE, "results", "md_gpu")

def csv_stats(gene):
    m = pd.read_csv(os.path.join(MD, gene + "_PFOS_md_metrics.csv"))
    c = pd.read_csv(os.path.join(MD, gene + "_PFOS_md_contacts.csv"))
    print("[%s] metrics rows=%d | RMSD_prot mean=%.2f last=%.2f | Rg mean=%.2f | RMSD_lig mean=%.2f last=%.2f"
          % (gene, len(m), m.rmsd_protein_A.mean(), m.rmsd_protein_A.iloc[-1],
             m.rg_protein_A.mean(), m.rmsd_ligand_A.mean(), m.rmsd_ligand_A.iloc[-1]))
    print("      hbonds mean=%.2f | contacts mean=%.2f max=%d | contacts>0 frames=%d/%d"
          % (c.hbonds_lig_prot.mean(), c.heavy_atom_contacts.mean(),
             c.heavy_atom_contacts.max(), (c.heavy_atom_contacts > 0).sum(), len(c)))

def build_state_distance():
    top = os.path.join(MD, "last_system_topology.pdb")
    t = mt.load_pdb(top)
    lig = t.topology.select("resname LIG")
    prot = t.topology.select("protein")
    print("[TMEM45A build state] atoms=%d lig=%d prot=%d" % (t.n_atoms, len(lig), len(prot)))
    xl = t.xyz[0][lig]; xp = t.xyz[0][prot]
    d = np.linalg.norm(xl[:, None, :] - xp[None, :, :], axis=2)
    print("  frame-0 min ligand-protein distance = %.3f nm (%.2f A)" % (d.min(), d.min() * 10))
    print("  contact pairs <0.4nm at build: %d" % int((d < 0.4).sum()))
    i, j = np.unravel_index(np.argmin(d), d.shape)
    a1 = t.topology.atom(int(lig[i])); a2 = t.topology.atom(int(prot[j]))
    print("  closest pair: %s(%s) -- %s(%s)" % (a1.name, a1.residue, a2.name, a2.residue))

def dcd_header(gene):
    from mdtraj.formats.dcd import DCDTrajectoryFile
    p = os.path.join(MD, gene + "_PFOS_prod.dcd")
    with DCDTrajectoryFile(p) as f:
        xyz, cell, _ = f.read(n_frames=1)
    print("[%s] DCD first-frame: %d atoms; file OK" % (gene, xyz.shape[1]))

if __name__ == "__main__":
    csv_stats("FLT1")
    csv_stats("TMEM45A")
    build_state_distance()
    dcd_header("FLT1")
    dcd_header("TMEM45A")
