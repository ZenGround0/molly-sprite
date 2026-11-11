"""Render cube-face views for every enzyme downloaded by download_enzymes.py.

Usage:
    python render_all.py [--cif glycolysis_cif] [--out glycolysis_views]
"""
import argparse
import csv
from pathlib import Path

from cube_views import make_cube_views


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cif", default="glycolysis_cif", help="Directory written by download_enzymes.py (default: %(default)s)")
    parser.add_argument("--out", default="glycolysis_views", help="Directory for the PNGs (default: %(default)s)")
    args = parser.parse_args()

    cif_dir, out_dir = Path(args.cif), Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(cif_dir / "manifest.csv") as f:
        rows = list(csv.DictReader(f))

    for row in rows:
        step = int(row["step"])
        print(f"[{step}/10] {row['enzyme']} ({row['pdb_id']})")
        make_cube_views(cif_dir / row["file"],
                        out=out_dir / f"{step:02d}_{row['pdb_id'].lower()}.png",
                        title=f"{step}. {row['enzyme']} ({row['pdb_id']})")


if __name__ == "__main__":
    main()
