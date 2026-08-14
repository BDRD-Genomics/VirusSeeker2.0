#!/usr/bin/env python3
import argparse
import csv
import sqlite3
from pathlib import Path


def load_table(cur, source: Path, table: str) -> None:
    cur.execute(
        f"CREATE TABLE {table} ("
        "accession TEXT, "
        "accession_version TEXT PRIMARY KEY, "
        "tax_id INTEGER, "
        "gi INTEGER)"
    )

    with source.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        rows = []
        for row in reader:
            if len(row) < 4:
                continue
            try:
                tax_id = int(row[2])
            except ValueError:
                continue
            try:
                gi = int(row[3])
            except ValueError:
                gi = 0
            rows.append((row[0], row[1], tax_id, gi))
            if len(rows) >= 100000:
                cur.executemany(f"INSERT OR REPLACE INTO {table} VALUES (?,?,?,?)", rows)
                rows.clear()
        if rows:
            cur.executemany(f"INSERT OR REPLACE INTO {table} VALUES (?,?,?,?)", rows)

    cur.execute(f"CREATE INDEX idx_{table}_acc ON {table}(accession)")
    cur.execute(f"CREATE INDEX idx_{table}_taxid ON {table}(tax_id)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nucl", required=True, type=Path)
    ap.add_argument("--prot", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    if args.output.exists():
        args.output.unlink()

    con = sqlite3.connect(args.output)
    try:
        cur = con.cursor()
        load_table(cur, args.nucl, "acc_taxid_nucl")
        load_table(cur, args.prot, "acc_taxid_prot")
        con.commit()
    finally:
        con.close()


if __name__ == "__main__":
    main()
