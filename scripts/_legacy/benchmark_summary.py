#!/usr/bin/env python3
"""
benchmark_summary.py — 收集 Snakemake benchmark 文件生成时间汇总

Snakemake 内置的 `benchmark:` 指令已为每个 rule 生成 TSV 文件，包含:
  - s (seconds): 任务 wall time
  - h:m:s: 任务 CPU 时间
  - max_rss, max_vms, max_uss, max_pss: 内存指标
  - io_in, io_out: 磁盘 IO
  - cpu_load: CPU 负载

本脚本收集所有 benchmark 文件，按 stage + lane 汇总。
"""

import os, sys, glob, csv, re
from collections import defaultdict
from datetime import datetime


def parse_benchmark_dir(log_dir: str) -> dict:
    """解析 Snakemake log 目录下的所有 benchmark TSV 文件."""
    benchmark_dir = os.path.join(log_dir, "benchmark")
    if not os.path.isdir(benchmark_dir):
        print(f"Benchmark dir not found: {benchmark_dir}", file=sys.stderr)
        return {}

    data = defaultdict(list)  # stage → [(lane, wall_seconds, cpu_seconds, max_rss)]

    for fpath in sorted(glob.glob(os.path.join(benchmark_dir, "*.tsv"))):
        fname = os.path.basename(fpath)
        # Parse stage and lane from filename: e.g. "guide_quant_lane_01.tsv"
        parts = fname.replace(".tsv", "").split("_")
        lane = None
        stage = None
        if "lane" in parts:
            idx = parts.index("lane")
            lane = f"lane_{parts[idx+1]}"
            stage = "_".join(parts[:idx])
        elif "whitelist" in fname:
            stage = "whitelist"
            lane = parts[-1] if parts[-1].startswith("lane") else "all"
        else:
            continue

        with open(fpath) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                wall_s = float(row.get("s", 0))
                cpu_s = float(row.get("h:m:s", 0))
                rss = float(row.get("max_rss", 0))
                io_in = float(row.get("io_in", 0))
                io_out = float(row.get("io_out", 0))
                data[stage].append({
                    "lane": lane,
                    "wall_s": wall_s,
                    "cpu_s": cpu_s,
                    "max_rss_mb": rss,
                    "io_in_mb": io_in,
                    "io_out_mb": io_out,
                })
    return data


def print_summary(data: dict) -> None:
    """打印时间汇总表."""
    import numpy as np

    print(f"\n{'='*85}")
    print(f"  BENCHMARK SUMMARY — {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print(f"{'='*85}")

    for stage, entries in sorted(data.items()):
        if not entries:
            continue
        walls = [e["wall_s"] for e in entries]
        rss_values = [e["max_rss_mb"] for e in entries if e["max_rss_mb"] > 0]

        print(f"\n── {stage} ({len(entries)} tasks) ──")
        print(f"  Wall time:  mean={np.mean(walls)/60:.1f} min, "
              f"median={np.median(walls)/60:.1f} min, "
              f"min={np.min(walls)/60:.1f} min, "
              f"max={np.max(walls)/60:.1f} min")
        print(f"  Sum wall:   {sum(walls)/60:.1f} min "
              f"({sum(walls)/3600:.1f} h)")
        if rss_values:
            print(f"  Memory RSS: mean={np.mean(rss_values):.0f} MB, "
                  f"max={np.max(rss_values):.0f} MB")

        # Show top 3 slowest
        sorted_entries = sorted(entries, key=lambda x: -x["wall_s"])
        print(f"  Slowest:    ", end="")
        for e in sorted_entries[:3]:
            print(f"{e['lane']}({e['wall_s']/60:.1f}m)", end="  ")
        print()

    # Total
    total_wall = sum(e["wall_s"] for entries in data.values() for e in entries)
    total_cpu = sum(e["cpu_s"] for entries in data.values() for e in entries)
    total_io_in = sum(e["io_in_mb"] for entries in data.values() for e in entries)
    total_io_out = sum(e["io_out_mb"] for entries in data.values() for e in entries)
    print(f"\n{'='*85}")
    print(f"  TOTAL wall time: {total_wall/60:.1f} min ({total_wall/3600:.1f} h)")
    print(f"  TOTAL CPU time:  {total_cpu/3600:.1f} h")
    print(f"  TOTAL IO read:   {total_io_in/1024:.1f} GB")
    print(f"  TOTAL IO write:  {total_io_out/1024:.1f} GB")
    print(f"{'='*85}")


def export_gantt_csv(data: dict, output_path: str) -> None:
    """导出甘特图友好格式 CSV."""
    rows = []
    for stage, entries in sorted(data.items()):
        for e in entries:
            rows.append({
                "stage": stage,
                "lane": e["lane"],
                "wall_s": e["wall_s"],
                "cpu_s": e["cpu_s"],
                "max_rss_mb": e["max_rss_mb"],
            })
    if not rows:
        return
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nGantt CSV exported: {output_path}")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Collect Snakemake benchmark timing")
    p.add_argument("log_dir", help="Snakemake log directory (contains benchmark/)")
    p.add_argument("--gantt", help="Export Gantt-ready CSV", default=None)
    args = p.parse_args()

    data = parse_benchmark_dir(args.log_dir)
    if not data:
        print("No benchmark data found. Is the run complete?", file=sys.stderr)
        sys.exit(1)
    print_summary(data)
    if args.gantt:
        export_gantt_csv(data, args.gantt)
