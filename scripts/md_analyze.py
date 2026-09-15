# md_analyze.py — MD 轨迹分析：RMSD / Rg / 氢键 / 配体-蛋白接触（conda-md 环境，mdtraj）
# 用法: python md_analyze.py <gene> <lig>
# 拓扑用同参数确定性重建（仅内存传递，不落盘 PDB）
import os, sys
import numpy as np
import pandas as pd
import mdtraj as mt
from openmm.app import PDBFile

BASE = os.environ.get("MD_BASE")
if not BASE:
    BASE = "E:/Zcode/0911" if os.path.isdir("E:/Zcode/0911") else os.getcwd()

def in_base(path):
    p = os.path.normpath(os.path.abspath(path))
    if os.path.commonpath([BASE, p]) != BASE:
        raise ValueError("path escapes project dir")
    return p

def build_openmm_topology(gene, lig):
    """按 md_build_run.py 相同参数确定性重建体系拓扑（内存中）。"""
    from openmm.app import PDBFile, Modeller, ForceField
    from openmm import unit
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator
    from openff.toolkit.topology import Molecule
    import io
    rec = PDBFile(in_base(BASE + "/data/receptors/" + gene + "_clean.pdb"))
    ligf = PDBFile(in_base(BASE + "/results/md/" + gene + "_" + lig + "_ligand_openmm.pdb"))
    mol = Molecule.from_file(in_base(BASE + "/data/compound/" + lig + "_deprot_3d.sdf"))
    mol.name = "LIG"
    mol.assign_partial_charges("mmff94")
    sgen = SMIRNOFFTemplateGenerator(molecules=[mol])
    xml = sgen.generate_residue_template(mol)
    mod = Modeller(rec.topology, rec.positions)
    ff = ForceField("amber14-all.xml")
    mod.addHydrogens(ff)
    mod.add(ligf.topology, ligf.positions)
    ffs = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    ffs.loadFile(io.StringIO(xml))
    mod.addSolvent(ffs, model="tip3p", padding=0.8 * unit.nanometer,
                   ionicStrength=0.15 * unit.molar)
    return mod.topology

def main(gene, lig):
    outdir = in_base(BASE + "/results/md")
    dcd = in_base(outdir + "/" + gene + "_" + lig + "_prod.dcd")
    top_pdb = in_base(outdir + "/last_system_topology.pdb")
    if not (os.path.exists(top_pdb) and os.path.getsize(top_pdb) > 10000):
        raise RuntimeError("last_system_topology.pdb missing/invalid - 请先跑对应的 md_build_run.py")
    mdtop = mt.Topology.from_openmm(PDBFile(top_pdb).topology)
    t = mt.load(dcd, top=mdtop)
    print("frames:", t.n_frames, "atoms:", t.n_atoms, flush=True)
    protein = t.topology.select("protein")
    ligatoms = t.topology.select("resname LIG")

    rmsd_prot = mt.rmsd(t, t, 0, atom_indices=protein) * 10  # nm -> A
    t.superpose(t, 0, atom_indices=protein)
    rmsd_lig = mt.rmsd(t, t, 0, atom_indices=ligatoms) * 10
    rg = mt.compute_rg(t.atom_slice(protein)) * 10
    times = np.arange(t.n_frames) * 0.02  # 每 20 ps 一帧（DT_FS=4fs x 5000 步）
    pd.DataFrame({"time_ns": times, "rmsd_protein_A": rmsd_prot,
                  "rmsd_ligand_A": rmsd_lig, "rg_protein_A": rg}).to_csv(
        in_base(outdir + "/" + gene + "_" + lig + "_md_metrics.csv"), index=False)

    # 配体-蛋白氢键（Wernet-Nilsson 判据；兼容不同 mdtraj 版本 API 位置）
    try:
        wn = mt.geometry.hydrogen_bonds.wernet_nilsson
    except AttributeError:
        from mdtraj.geometry import hbond
        wn = hbond.wernet_nilsson
    hb = wn(t, exclude_water=True)
    res_of = t.topology.atom
    hb_lig = []
    for frame_hb in hb:
        c = 0
        for d, h, a in frame_hb:
            rd = res_of(d).residue.name
            ra = res_of(a).residue.name
            if (rd == "LIG") != (ra == "LIG"):
                c += 1
        hb_lig.append(c)

    # 配体-蛋白重原子接触（<4A）
    contacts = []
    for fi in range(t.n_frames):
        xyz_l = t.xyz[fi][ligatoms]
        xyz_p = t.xyz[fi][protein]
        d = np.linalg.norm(xyz_l[:, None, :] - xyz_p[None, :, :], axis=2)
        contacts.append(int((d < 0.4).sum()))
    pd.DataFrame({"time_ns": times, "hbonds_lig_prot": hb_lig,
                  "heavy_atom_contacts": contacts}).to_csv(
        in_base(outdir + "/" + gene + "_" + lig + "_md_contacts.csv"), index=False)
    print("ANALYZE DONE", flush=True)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
