#!/usr/bin/env python3
"""Parse samples.yaml and output TSV for bash consumption."""
import yaml, sys

with open(sys.argv[1]) as f:
    cfg = yaml.safe_load(f)

groups = cfg.get('groups', {})
if not groups:
    sys.exit(1)

for gname in sorted(groups.keys()):
    g = groups[gname]
    d = g.get('sgRNA_fastq_dir', '')
    p1 = g.get('sgRNA_r1_pattern', '')
    p2 = g.get('sgRNA_r2_pattern', '')
    print(f'{gname}\t{d}\t{p1}\t{p2}')
