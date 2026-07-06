#!/usr/bin/env python3
"""
Convert standardized assignment CSV to per-cell perturbation call.

Reads the unified-schema assignment CSV and produces perturbation_obs.csv —
a (cell -> perturbation identity) mapping ready for multimodal integration.

Supports three guide_design modes:

    single — 1 guide per perturbation (CRISPRko). top-1 -> gene symbol.
    dual   — 2 guides per construct (CRISPRi). top-2 -> construct pair_id.
    multi  — N guides per construct. top-N -> construct_id.

The guide_csv format depends on guide_design (see --help for schemas).

Usage:
    python make_perturbation_obs.py \\
        --input assignments.csv --output perturbation_obs.csv \\
        --guide-design dual --guide-csv guide_library.csv --method pgmm_em
"""

import argparse, csv, time
from collections import defaultdict


CONFIDENCE_TIERS = {
    "pgmm_em": {
        "score_type": "prob_gaussian",
        "high": 0.90,    # >= 0.90: high confidence
        "mid":  0.75,    # >= 0.75: assignment threshold
    },
    "umi_threshold": {
        "score_type": "umi_count",
        "high": 10,      # >= 10 UMI: high
        "mid":  5,       # >= 5: medium
    },
    "fishash": {
        "score_type": "neg_log_pval",
        "high": 15,      # more negative = more significant
        "mid":  10,
    },
}

# How many top guides to consider per cell for construct matching,
# keyed by guide_design.  "multi" uses all available guides.
TOP_GUIDES = {"single": 1, "dual": 2}


def classify_confidence(method, score):
    tiers = CONFIDENCE_TIERS.get(method)
    if tiers is None:
        return "unknown"
    if score >= tiers["high"]:
        return "high"
    elif score >= tiers["mid"]:
        return "mid"
    else:
        return "low"


# ---------------------------------------------------------------------------
# Guide mapping loaders — unified internal format: {guide_id: (construct_id, gene)}
# ---------------------------------------------------------------------------

def _expand_dual_csv(csv_path):
    """dual-guide CSV (sgID_A, sgID_B, gene, pair_id) -> {guide_id: (pid, gene)}."""
    guide_map = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            gene = row.get("gene", "").strip()
            pid = (row.get("pair_id") or row.get("unique sgRNA pair ID", "")).strip()
            for col in ("sgID_A", "sgID_B"):
                sg = row.get(col, "").strip()
                if sg:
                    guide_map[sg] = (pid, gene)
    return guide_map


def _load_single_or_multi_csv(csv_path):
    """single/multi CSV (guide_id, gene [, construct_id]) -> {guide_id: (cid, gene)}.

    - single mode: construct_id column is optional/absent -> cid = "".
    - multi  mode: construct_id column is required.
    """
    guide_map = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            gid = row.get("guide_id", "").strip()
            gene = row.get("gene", "").strip()
            if not gid:
                continue
            cid = row.get("construct_id", "").strip() if "construct_id" in row else ""
            guide_map[gid] = (cid, gene)
    return guide_map


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Standardized assignment -> per-cell perturbation call")
    parser.add_argument("--input", required=True,
                        help="Standardized assignment CSV")
    parser.add_argument("--output", required=True,
                        help="Output perturbation_obs.csv")
    parser.add_argument("--guide-design", required=True,
                        choices=["single", "dual", "multi"],
                        help="single | dual | multi")
    parser.add_argument("--guide-csv", default="",
                        help="Guide mapping CSV. "
                             "single/multi: columns guide_id,gene[,construct_id]. "
                             "dual: columns sgID_A,sgID_B,gene,pair_id. "
                             "single mode: omit if guide_id == gene name.")
    parser.add_argument("--method", required=True,
                        help="Assignment method name")
    args = parser.parse_args()

    t0 = time.time()

    # ---- Load guide mapping ----
    guide_map = {}   # guide_id -> (construct_id, gene)

    if args.guide_design == "dual":
        if not args.guide_csv:
            parser.error("--guide-csv required for dual guide_design")
        guide_map = _expand_dual_csv(args.guide_csv)
        n_constructs = len(set(cid for cid, _ in guide_map.values() if cid))
        print(f"Dual-guide: {len(guide_map):,} sgIDs -> {n_constructs} constructs")

    elif args.guide_design == "multi":
        if not args.guide_csv:
            parser.error("--guide-csv required for multi guide_design")
        guide_map = _load_single_or_multi_csv(args.guide_csv)
        has_cid = any(cid for cid, _ in guide_map.values())
        if not has_cid:
            parser.error("multi guide_design requires construct_id column in guide_csv")
        n_constructs = len(set(cid for cid, _ in guide_map.values() if cid))
        print(f"Multi-guide: {len(guide_map)} guides -> {n_constructs} constructs")

    else:  # single
        if args.guide_csv:
            guide_map = _load_single_or_multi_csv(args.guide_csv)
            print(f"Single-guide mapping: {len(guide_map)} guide_id -> gene entries")
        else:
            print("Single-guide: no guide_csv (guide_id used as gene name)")

    # ---- Load assignments ----
    top_n = TOP_GUIDES.get(args.guide_design, 1)
    print(f"Loading {args.input} (top-{top_n} per cell) ...", end=' ', flush=True)
    cell_guides = defaultdict(list)
    with open(args.input) as f:
        for row in csv.DictReader(f):
            cell = row.get("cell_barcode", "").strip()
            guide = row.get("guide_id", "").strip()
            rank = int(row.get("rank", 0))
            score = float(row.get("score", 0))
            if cell and guide and rank <= top_n:
                cell_guides[cell].append((rank, guide, score))

    for cell in cell_guides:
        cell_guides[cell].sort(key=lambda x: x[0])

    n_cells = len(cell_guides)
    n_total = sum(len(v) for v in cell_guides.values())
    print(f"{n_cells:,} cells, {n_total:,} guide entries  [{time.time()-t0:.1f}s]")

    # ---- Build per-cell perturbation call ----
    rows = []
    n_ambiguous = 0
    n_na = 0

    for cell, guides in cell_guides.items():
        top1_guide = guides[0][1]
        top1_score = guides[0][2]
        confidence = classify_confidence(args.method, top1_score)

        if args.guide_design == "single":
            # ---- single: top-1 guide -> gene ----
            if guide_map:
                _, gene = guide_map.get(top1_guide, ("", top1_guide))
                if not gene:
                    gene = top1_guide
            else:
                gene = top1_guide  # guide_id IS the gene name

            perturbation = gene if gene else "NA"
            if perturbation == "NA":
                n_na += 1

            rows.append({
                "cell_barcode": cell,
                "perturbation": perturbation,
                "n_guides_assigned": 1,
                "assignment_score": f"{top1_score:.6f}",
                "assignment_confidence": confidence,
                "assignment_method": args.method,
            })

        elif args.guide_design == "dual":
            # ---- dual: top-2 guides -> same construct? ----
            cids = []
            for _, g, _ in guides:
                entry = guide_map.get(g)
                if entry:
                    cid, _gene = entry
                    if cid:
                        cids.append(cid)

            if len(cids) >= 2 and len(set(cids)) == 1:
                perturbation = cids[0]
                n_guides_used = len(cids)
            elif len(cids) >= 1 and len(set(cids)) == 1:
                perturbation = cids[0]
                n_guides_used = len(cids)
            elif len(cids) == 0:
                perturbation = "NA"
                n_guides_used = 0
                n_na += 1
            else:
                perturbation = "ambiguous_pair"
                n_guides_used = len(guides)
                n_ambiguous += 1

            rows.append({
                "cell_barcode": cell,
                "perturbation": perturbation,
                "n_guides_assigned": n_guides_used,
                "assignment_score": f"{top1_score:.6f}",
                "assignment_confidence": confidence,
                "assignment_method": args.method,
            })

        else:  # multi
            # ---- multi: top-N guides -> same construct? ----
            cids = []
            n_with_cid = 0
            for _, g, _ in guides:
                entry = guide_map.get(g)
                if entry:
                    cid, _gene = entry
                    if cid:
                        cids.append(cid)
                        n_with_cid += 1

            if len(cids) >= 1 and len(set(cids)) == 1:
                perturbation = cids[0]
                n_guides_used = n_with_cid
            elif len(cids) == 0:
                perturbation = "NA"
                n_guides_used = 0
                n_na += 1
            else:
                perturbation = "ambiguous_construct"
                n_guides_used = n_with_cid
                n_ambiguous += 1

            rows.append({
                "cell_barcode": cell,
                "perturbation": perturbation,
                "n_guides_assigned": n_guides_used,
                "assignment_score": f"{top1_score:.6f}",
                "assignment_confidence": confidence,
                "assignment_method": args.method,
            })

    # ---- Write ----
    with open(args.output, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=[
            "cell_barcode", "perturbation", "n_guides_assigned",
            "assignment_score", "assignment_confidence", "assignment_method",
        ])
        w.writeheader()
        w.writerows(rows)

    print(f"  -> {args.output}  ({len(rows):,} cells)  [{time.time()-t0:.1f}s]")
    if n_ambiguous:
        pct = n_ambiguous / max(n_cells, 1) * 100
        tag = "pair" if args.guide_design == "dual" else "construct"
        print(f"  Note: {n_ambiguous} cells with ambiguous_{tag} ({pct:.1f}%)")
    if n_na:
        pct = n_na / max(n_cells, 1) * 100
        print(f"  Note: {n_na} cells with NA perturbation ({pct:.1f}%)")


if __name__ == "__main__":
    main()
