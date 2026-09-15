#!/bin/bash
# 解析全部 _out.pdbqt 的最优打分 -> docking_summary.csv
cd /e/Zcode/0911
echo "gene,ligand,affinity_kcal_mol" > results/docking/docking_summary.csv
for f in results/docking/*_out.pdbqt; do
  [ -s "$f" ] || continue
  b=$(basename "$f" _out.pdbqt); g="${b%_*}"; l="${b##*_}"
  a=$(grep -m1 "REMARK VINA RESULT:" "$f" | awk '{print $4}')
  echo "$g,$l,$a" >> results/docking/docking_summary.csv
done
echo "PARSED"
cat results/docking/docking_summary.csv
