#!/usr/bin/env python3
import argparse
from pathlib import Path

from ete3 import NCBITaxa


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxdump", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    try:
        ncbi = NCBITaxa(dbfile=args.output, taxdump_file=args.taxdump)
    except TypeError:
        from ete3.ncbi_taxonomy.ncbiquery import update_db
        update_db(args.output, args.taxdump)
        ncbi = NCBITaxa(dbfile=args.output)

    names = ncbi.get_taxid_translator([9606])
    if 9606 not in names:
        raise SystemExit("ETE3 validation failed: taxid 9606 was not resolved")

    print("ETE3 taxa.sqlite created:", args.output)
    print("ETE3 validation:", names[9606])

    sidecar = Path(args.output + ".traverse.pkl")
    if sidecar.exists():
        print("ETE3 traversal cache created:", sidecar)


if __name__ == "__main__":
    main()
