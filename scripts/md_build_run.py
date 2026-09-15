# md_build_run.py — OpenMM 全原子 MD：蛋白-PFOS 复合物构建 + 平衡 + 产出（conda-md 环境）
# 用法: python md_build_run.py <gene> <lig> <ns>
# 依赖: openmm, openmmforcefields, openff-toolkit, rdkit（conda-forge env: .conda-md）
import os, sys, time
import numpy as np
import openmm
from openmm import unit
from openmm.app import PDBFile, Modeller, ForceField, DCDReporter, StateDataReporter, CheckpointReporter, HBonds
from openmmforcefields.generators import SMIRNOFFTemplateGenerator, SystemGenerator
from openff.toolkit.topology import Molecule

DT_FS = 4.0  # 氢质量重分配（HMR 3 Da）后允许 4 fs 步长；三体系统一协议

# 可移植根目录：优先 MD_BASE 环境变量，本地默认 E:/Zcode/0911，云端用当前目录
BASE = os.environ.get("MD_BASE")
if not BASE:
    BASE = "E:/Zcode/0911" if os.path.isdir("E:/Zcode/0911") else os.getcwd()

def in_base(path):
    p = os.path.normpath(os.path.abspath(path))
    if os.path.commonpath([BASE, p]) != BASE:
        raise ValueError("path escapes project dir")
    return p

def main(gene, lig, ns, eq_ps=200):
    outdir = in_base(BASE + "/results/md")
    os.makedirs(outdir, exist_ok=True)
    rec_pdb = in_base(BASE + "/data/receptors/" + gene + "_clean.pdb")
    lig_pdb = in_base(BASE + "/results/md/" + gene + "_" + lig + "_ligand_openmm.pdb")
    lig_sdf = in_base(BASE + "/data/compound/" + lig + "_deprot_3d.sdf")
    dcd_out = outdir + "/" + gene + "_" + lig + "_prod.dcd"
    log_out = outdir + "/" + gene + "_" + lig + "_md.log"
    chk_out = outdir + "/" + gene + "_" + lig + "_checkpoint.chk"

    # 1. 配体 openff Molecule（去质子化阴离子；MMFF94 电荷——Windows 无 AmberTools/AM1-BCC 的可复现替代方案）
    import io
    mol = Molecule.from_file(lig_sdf)
    mol.name = "LIG"
    mol.assign_partial_charges("mmff94")   # 电荷总和守恒 = 配体形式电荷
    lig_template = SMIRNOFFTemplateGenerator(molecules=[mol])
    lig_ffxml = lig_template.generate_residue_template(mol)

    # 2. 蛋白先加氢（LIG 无模板会令 addHydrogens 抛错，故在合并配体之前执行），再合并配体
    rec = PDBFile(rec_pdb)
    lig = PDBFile(lig_pdb)
    modeller = Modeller(rec.topology, rec.positions)
    ff = ForceField("amber14-all.xml")
    print("stage: protein loaded", flush=True)
    modeller.addHydrogens(ff)
    print("stage: protein H added", flush=True)
    modeller.add(lig.topology, lig.positions)
    print("stage: ligand merged", flush=True)

    # 3. 显式溶剂（ForceField 需预先注册配体模板，addSolvent 内部 createSystem 才能匹配 LIG）
    ff_solvent = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    ff_solvent.loadFile(io.StringIO(lig_ffxml))
    modeller.addSolvent(ff_solvent,
                        model="tip3p", padding=0.8 * unit.nanometer, ionicStrength=0.15 * unit.molar)
    print("system built:", modeller.topology.getNumAtoms(), "atoms", flush=True)

    # 保存体系拓扑+初始坐标（addHydrogens/addSolvent 有跨运行不确定性，
    # 轨迹分析必须复用本次运行的真实拓扑；文件名固定，链条顺序执行互不覆盖）
    import mdtraj as mt
    mdtop = mt.Topology.from_openmm(modeller.topology)
    xyz_nm = np.array(modeller.positions.value_in_unit(unit.nanometer))
    mt.Trajectory(xyz=xyz_nm, topology=mdtop).save_pdb("results/md/last_system_topology.pdb")
    print("topology saved: results/md/last_system_topology.pdb", flush=True)

    # 4. SystemGenerator（蛋白 ff14SB + TIP3P + 配体 Sage/SMIRNOFF + 预存 MMFF94 电荷）
    #    HMR：氢质量加重到 3 Da + HBonds 刚性约束 → 4 fs 步长稳定
    gen = SystemGenerator(forcefields=["amber14-all.xml", "amber14/tip3p.xml"],
                          small_molecule_forcefield="openff-2.2.1",
                          molecules=[mol],
                          forcefield_kwargs={"constraints": HBonds, "rigidWater": True,
                                             "hydrogenMass": 3 * unit.dalton})
    system = gen.create_system(modeller.topology)
    print("system created; forces:", [system.getForce(i).getName() for i in range(system.getNumForces())], flush=True)

    # 5. 平台自动选择：CUDA > OpenCL > CPU
    integrator = openmm.LangevinMiddleIntegrator(310 * unit.kelvin, 1 / unit.picosecond, DT_FS * unit.femtoseconds)
    platform, prop = None, {}
    for name in ("CUDA", "OpenCL", "CPU"):
        try:
            platform = openmm.Platform.getPlatformByName(name)
            if name == "CUDA":
                prop = {"Precision": "mixed"}
            elif name == "CPU":
                prop = {"Threads": "10"}
            print("platform:", name, flush=True)
            break
        except Exception:
            continue
    if platform is None:
        raise RuntimeError("no OpenMM platform available")
    simulation = openmm.app.Simulation(modeller.topology, system, integrator, platform, prop)
    simulation.context.setPositions(modeller.positions)
    print("minimizing...", flush=True)
    simulation.minimizeEnergy(maxIterations=20000)
    simulation.context.setVelocitiesToTemperature(310 * unit.kelvin, 123)

    # 6. NPT 平衡 200 ps + 产出 NPT
    simulation.system.addForce(openmm.MonteCarloBarostat(1 * unit.atmosphere, 310 * unit.kelvin, 25))
    simulation.context.reinitialize(preserveState=True)
    simulation.reporters.append(StateDataReporter(log_out, 500, step=True, time=True, potentialEnergy=True,
                                                  temperature=True, volume=True, density=True, speed=True))
    simulation.reporters.append(CheckpointReporter(chk_out, 5000))
    eq_steps = int(eq_ps / (DT_FS / 1000))  # 平衡时长（ps）
    print("equilibrating %g ps..." % eq_ps, flush=True)
    simulation.step(eq_steps)
    ns_steps = int(ns * 1e6 / DT_FS)  # ns -> 步数
    simulation.reporters.append(DCDReporter(dcd_out, 5000))  # 每 20 ps 一帧
    print("production %g ns..." % ns, flush=True)
    t0 = time.time()
    simulation.step(ns_steps)
    print("MD DONE %.1f min" % ((time.time() - t0) / 60), flush=True)

if __name__ == "__main__":
    eq = float(sys.argv[4]) if len(sys.argv) > 4 else 200.0
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]), eq)
