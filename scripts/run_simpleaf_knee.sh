#!/bin/bash
# Run simpleaf with knee-calling on all 48 lanes, then translate barcodes and merge.
# No GEX whitelist used — simpleaf's internal knee detection selects cells.
# Post-quant: Feature barcodes → GEX format via 10x translation table.
#
# Usage: bash run_simpleaf_knee.sh

set -uo pipefail

# ── Config ──
CORES=12                                    # per-lane simpleaf threads
INDEX=/data/yunzliu/simpleaf_index_k13       # piscem index
T2G=/data/yunzliu/references/t2g_2col_guide.tsv
FASTQ_BASE=/data/yunzliu/raw_fastq/full/by_lane
OUT_BASE=/data/yunzliu/results/guide_extraction/simpleaf_knee_final
LOG_BASE=/data/yunzliu/logs/guide_extraction/simpleaf_knee_final
TRANS_TABLE=/data/yunzliu/scdata/cellranger_ref/cellranger_whitelist_translation_3v3.txt
CHEMISTRY=10xv3
RESOLUTION=parsimony-gene
CONDA_ENV=simpleaf
CONDA_BASE=$HOME/software/miniconda3

source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"
export ALEVIN_FRY_HOME="$CONDA_BASE/envs/$CONDA_ENV/opt/alevin-fry"

mkdir -p "$OUT_BASE" "$LOG_BASE"

echo "=== simpleaf knee mode: 48 lanes ==="
echo "Chemistry: $CHEMISTRY | Resolution: $RESOLUTION | Cores/lane: $CORES"
echo ""

for lane_num in $(seq -w 1 48); do
    LANE="lane_${lane_num}"
    OUT_DIR="$OUT_BASE/$LANE/simpleaf_quant"
    LOG_FILE="$LOG_BASE/${LANE}.log"
    
    # Skip if already done
    if [ -f "$OUT_DIR/af_quant/alevin/quants_mat.mtx" ]; then
        echo "[$LANE] Already done, skipping"
        continue
    fi
    
    echo "[$LANE] Running simpleaf quant (knee)..."
    mkdir -p "$OUT_DIR"
    
    R1=$(ls $FASTQ_BASE/$LANE/lane${lane_num}_sgRNA_*_R1_001.fastq.gz 2>/dev/null | sort | tr '\n' ',' | sed 's/,$//')
    R2=$(ls $FASTQ_BASE/$LANE/lane${lane_num}_sgRNA_*_R2_001.fastq.gz 2>/dev/null | sort | tr '\n' ',' | sed 's/,$//')
    
    if [ -z "$R1" ]; then
        echo "[$LANE] No FASTQ found, skipping"
        continue
    fi
    
    simpleaf quant \
        --chemistry "$CHEMISTRY" \
        --output "$OUT_DIR" \
        --threads $CORES \
        --index "$INDEX" \
        --reads1 "$R1" \
        --reads2 "$R2" \
        --t2g-map "$T2G" \
        --resolution "$RESOLUTION" \
        --knee \
        >> "$LOG_FILE" 2>&1
    
    echo "[$LANE] Done."
done

echo ""
echo "=== All lanes quantified. Now translating barcodes and merging... ==="

# Translate and merge via Python
python3 << 'PYEOF'
import gzip, os, re, sys
import scipy.sparse as sp
from collections import OrderedDict

OUT_BASE = "/data/yunzliu/results/guide_extraction/simpleaf_knee_final"
TRANS_TABLE = "/data/yunzliu/scdata/cellranger_ref/cellranger_whitelist_translation_3v3.txt"

# Load translation Feature→GEX
feat_to_gex = {}
with open(TRANS_TABLE) as f:
    for line in f:
        feat, gex = line.strip().split()
        feat_to_gex[feat] = gex

# Collect all features across lanes
all_f = OrderedDict()
lane_data = []

for lane_num in range(1, 49):
    lane = f"lane_{lane_num:02d}"
    quant_dir = f"{OUT_BASE}/{lane}/simpleaf_quant/af_quant/alevin"
    mtx_f = f"{quant_dir}/quants_mat.mtx"
    rows_f = f"{quant_dir}/quants_mat_rows.txt"
    cols_f = f"{quant_dir}/quants_mat_cols.txt"
    
    if not os.path.exists(mtx_f):
        print(f"  {lane}: SKIP (no output)")
        continue
    
    # Load features
    fids = []
    with open(rows_f) as f:
        for line in f:
            fid = line.strip()
            fids.append(fid)
            if fid not in all_f:
                all_f[fid] = len(all_f)
    
    # Load barcodes and translate Feature→GEX
    bcs = []
    untrans = 0
    with open(cols_f) as f:
        for line in f:
            bc = line.strip()
            bc_16 = bc.split('-')[0] if '-' in bc else bc
            gex_bc = feat_to_gex.get(bc_16)
            if gex_bc:
                bcs.append(f"{gex_bc}-L{lane_num:02d}")
            else:
                bcs.append(f"{bc_16}-L{lane_num:02d}")
                untrans += 1
    
    # Load matrix
    mtx = sp.mmread(mtx_f).tocsc()
    
    lane_data.append((lane, mtx, bcs, fids))
    print(f"  {lane}: {len(bcs)} cells, {mtx.nnz} nnz, {untrans} untranslated")

print(f"\nFeatures: {len(all_f)}")
fl = list(all_f.keys())
g_idx = {f: i for i, f in enumerate(fl)}

# Stack matrices
blocks = []
all_bc = []
for lane, mtx, bcs, fids in lane_data:
    local_to_global = [g_idx[fid] for fid in fids]
    coo = mtx.tocoo()
    rg = [local_to_global[r] for r in coo.row]
    aligned = sp.coo_matrix((coo.data, (rg, coo.col)), shape=(len(fl), len(bcs))).tocsc()
    blocks.append(aligned)
    all_bc.extend(bcs)

merged = sp.hstack(blocks, format='csc')
print(f"Merged: {merged.shape[0]}f x {merged.shape[1]}c, {merged.nnz:,} nnz")

# Write MEX
MERGED = f"{OUT_BASE}/merged"
os.makedirs(MERGED, exist_ok=True)
mt = merged.T.tocsc()

with gzip.open(f"{MERGED}/merged_matrix.mtx.gz", 'wt') as f:
    f.write("%%MatrixMarket matrix coordinate integer general\n")
    f.write(f"{mt.shape[0]} {mt.shape[1]} {mt.nnz}\n")
    for r, c, v in zip(*sp.find(mt)):
        f.write(f"{r+1} {c+1} {int(v)}\n")

with gzip.open(f"{MERGED}/merged_barcodes.tsv.gz", 'wt') as f:
    for bc in all_bc: f.write(bc + "\n")

with gzip.open(f"{MERGED}/merged_features.tsv.gz", 'wt') as f:
    for fid in fl: f.write(f"{fid}\t{fid}\tCRISPR Guide Capture\n")

for fn in ['merged_matrix.mtx.gz', 'merged_barcodes.tsv.gz', 'merged_features.tsv.gz']:
    print(f"  {fn}: {os.path.getsize(MERGED+'/'+fn)//1024:,} KB")
print(f"\nDone: {MERGED}/")
PYEOF

echo "=== Complete ==="
