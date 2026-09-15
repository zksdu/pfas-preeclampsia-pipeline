# md_diag_eng.py — ENG 诊断：重建拓扑 → frame-0 配体距离 + 按残基 RMSF 找无序尾 → 核心域 RMSD
import os
import numpy as np
import mdtraj as mt
from openmm.app import PDBFile, Modeller, ForceField
from openmm import unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit.topology import Molecule
import io

BASE = "E:/Zcode/0911"
MD = os.path.join(BASE, "results", "md_gpu")

def rebuild_eng():
    rec = PDBFile(os.path.join(BASE, "data/receptors/ENG_clean.pdb"))
    lig = PDBFile(os.path.join(BASE, "results/md/ENG_PFOS_ligand_openmm.pdb"))
    mol = Molecule.from_file(os.path.join(BASE, "data/compound/PFOS_deprot_3d.sdf"))
    mol.name = "LIG"
    mol.assign_partial_charges("mmff94")
    sgen = SMIRNOFFTemplateGenerator(molecules=[mol])
    xml = sgen.generate_residue_template(mol)
    mod = Modeller(rec.topology, rec.positions)
    ff = ForceField("amber14-all.xml")
    mod.addHydrogens(ff)
    mod.add(lig.topology, lig.positions)
    ffs = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    ffs.loadFile(io.StringIO(xml))
    mod.addSolvent(ffs, model="tip3p", padding=0.8 * unit.nanometer, ionicStrength=0.15 * unit.molar)
    return mod

def main():
    mod = rebuild_eng()
    n = mod.topology.getNumAtoms()
    print("ENG rebuild atoms:", n)
    dcd = os.path.join(MD, "ENG_PFOS_prod.dcd")
    from mdtraj.formats.dcd import DCDTrajectoryFile
    with DCDTrajectoryFile(dcd) as f:
        xyz0, cell, _ = f.read(n_frames=1)
    print("ENG DCD atoms:", xyz0.shape[1])
    if xyz0.shape[1] != n:
        print("!! topology mismatch (%d vs %d) - 用 DCD 帧坐标 + 拓扑残基归属不可行，改用仅重算方案" % (n, xyz0.shape[1]))
        return
    mdtop = mt.Topology.from_openmm(mod.topology)
    t = mt.Trajectory(xyz=xyz0.astype(np.float32), topology=mdtop)
    lig = mdtop.select("resname LIG")
    prot = mdtop.select("protein")
    xl = t.xyz[0][lig]; xp = t.xyz[0][prot]
    d = np.linalg.norm(xl[:, None, :] - xp[None, :, :], axis=2)
    print("ENG production frame-0 min ligand-protein distance: %.3f nm" % d.min())

    # 载入全轨迹算 RMSF（426k 原子 x 500 帧 ≈ 2.5 GB 内存，可接受）
    t = mt.load(dcd, top=mdtop)
    print("loaded full ENG traj:", t.n_frames, "frames")
    t.superpose(t, 0, atom_indices=prot)
    ca = mdtop.select("protein and name CA")
    rmsf = mt.rmsf(t, 0, atom_indices=ca) * 10
    # 找低波动核心：RMSF < 3 A 的 Cα
    core_ca = ca[rmsf < 3.0]
    print("core CAs (RMSF<3A): %d / %d" % (len(core_ca), len(ca)))
    rmsd_core = mt.rmsd(t, t, 0, atom_indices=core_ca) * 10
    print("ENG core-RMSD: mean=%.2f last=%.2f (all-atom was mean=10.6 last=13.3)"
          % (rmsd_core.mean(), rmsd_core[-1]))
    top_flex = np.argsort(-rmsf)[:10]
    print("top-10 flexible residues (idx):", sorted(core_ca_idx if False else [int(ca[i]) for i in top_flex]))
    # 保存核心 RMSD 供画图
    times = np.arange(t.n_frames) * 0.02
    import pandas as pd
    pd.DataFrame({"time_ns": times, "rmsd_core_A": rmsd_core}).to_csv(
        os.path.join(MD, "ENG_PFOS_md_core_rmsd.csv"), index=False)
    print("saved ENG_PFOS_md_core_rmsd.csv")

if __name__ == "__main__":
    main()
