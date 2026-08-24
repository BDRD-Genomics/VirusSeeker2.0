#!/usr/bin/env python3
from collections import defaultdict
from pathlib import Path
import sys

nodes_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
root = int(sys.argv[3])

children = defaultdict(list)

with nodes_path.open("r", encoding="utf-8") as fh:
    for line in fh:
        fields = [x.strip() for x in line.split("|")]
        if len(fields) < 2 or not fields[0] or not fields[1]:
            continue
        taxid = int(fields[0])
        parent = int(fields[1])
        if taxid != parent:
            children[parent].append(taxid)

stack = [root]
seen = set()

while stack:
    taxid = stack.pop()
    if taxid in seen:
        continue
    seen.add(taxid)
    stack.extend(children.get(taxid, ()))

with output_path.open("w", encoding="utf-8") as out:
    for taxid in sorted(seen):
        out.write(f"{taxid}\n")

print(f"Wrote {len(seen):,} viral taxids to {output_path}")
