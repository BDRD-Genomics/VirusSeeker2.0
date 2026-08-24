#!/usr/bin/env python3
import argparse
import sqlite3


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("db")
    args = ap.parse_args()

    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    try:
        cur = con.cursor()
        expected = [
            ("acc_taxid_nucl", "NC_045512.2", 2697049),
            ("acc_taxid_nucl", "NC_012920.1", 9606),
            ("acc_taxid_prot", "YP_009724389.1", 2697049),
        ]
        for table, accession, taxid in expected:
            row = cur.execute(
                f"SELECT tax_id FROM {table} WHERE accession_version=?", (accession,)
            ).fetchone()
            if not row or row[0] != taxid:
                raise SystemExit(f"ERROR: {accession} -> taxid {taxid} missing from {table}")
        nprot = cur.execute(
            "SELECT COUNT(*) FROM acc_taxid_prot WHERE tax_id=2697049"
        ).fetchone()[0]
        if nprot < 1:
            raise SystemExit("ERROR: no SARS-CoV-2 protein mappings found")
        print(f"Public SARS taxonomy mappings OK; {nprot} SARS-CoV-2 protein accession(s)")
    finally:
        con.close()


if __name__ == "__main__":
    main()
