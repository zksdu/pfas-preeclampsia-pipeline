#!/bin/bash
# 受体制备 v2：AlphaFold v6 结构下载 -> ATOM 过滤 -> obabel 加氢转 PDBQT（刚性受体）
cd /e/Zcode/0911
mkdir -p data/receptors logs

prep_one() {  # prep_one <基因> <uniprot>
  local gene="$1" up="$2"
  local raw="data/receptors/AF_${gene}.pdb"
  local clean="data/receptors/${gene}_clean.pdb"
  local rec="data/receptors/${gene}_receptor.pdbqt"
  if [ ! -s "$clean" ]; then
    curl -sSL --ssl-no-revoke --http1.1 --retry 4 --max-time 300 -o "$raw" \
      "https://alphafold.ebi.ac.uk/files/AF-${up}-F1-model_v6.pdb"
    local sz=$(wc -c < "$raw")
    if [ "$sz" -lt 10000 ]; then echo "[FAIL] $gene ($up): only $sz bytes"; return 1; fi
    grep -E "^(ATOM|TER)" "$raw" > "$clean"
    echo "[OK] $gene ($up): $sz bytes"
  fi
  obabel -ipdb "$clean" -opdbqt -O "$rec" -xr -p 7.4 2>/dev/null
  echo "    $gene receptor: $(wc -c < "$rec") bytes"
}

prep_one FLT1     P17948
prep_one ENG      P17813
prep_one TREM1    Q9NP99
prep_one LEP      P41159
prep_one NDRG1    Q92597
prep_one P4HA1    P13674
prep_one COL17A1  Q9UMD9
prep_one CP       P00450
prep_one CRH      P06850
prep_one IL1R2    P27930
prep_one NRIP1    P48552
prep_one BHLHE40  O14503
prep_one HILPDA   Q9Y5L2
prep_one SERPINA3 P01011
prep_one SFXN3    Q9BWM7
prep_one SPX      Q9BT56
prep_one TMEM45A  Q9NWC5
prep_one C12ORF75 Q8TAD7
echo "RECEPTOR PREP DONE"
