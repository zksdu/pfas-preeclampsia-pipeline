# receptor_convert.py — *_clean.pdb -> *_receptor.pdbqt（Vina 刚性受体：极性氢、合并非极性氢）
import glob, os
from openbabel import openbabel as obm

obm.obErrorLog.SetOutputLevel(obm.obError)
os.chdir(r"E:/Zcode/0911/data/receptors")

ok, fail = [], []
for clean in sorted(glob.glob("*_clean.pdb")):
    gene = clean[:-10]  # strip _clean.pdb
    out = f"{gene}_receptor.pdbqt"
    conv = obm.OBConversion()
    conv.SetInAndOutFormats("pdb", "pdbqt")
    conv.AddOption("r", conv.OUTOPTIONS)  # rigid receptor
    conv.AddOption("c", conv.OUTOPTIONS)  # single combined molecule
    mol = obm.OBMol()
    if not conv.ReadFile(mol, clean):
        fail.append((gene, "read")); continue
    mol.AddHydrogens(True, True, 7.4)     # polarOnly, correct pH
    if not conv.WriteFile(mol, out):
        fail.append((gene, "write")); continue
    ok.append((gene, os.path.getsize(out)))

for g, s in ok:
    print(f"[OK] {g}: {s} bytes")
for g, why in fail:
    print(f"[FAIL] {g}: {why}")
print(f"DONE ok={len(ok)} fail={len(fail)}")
