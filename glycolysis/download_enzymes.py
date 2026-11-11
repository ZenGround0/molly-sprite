"""Download mmCIF structures for the ten enzymes of glycolysis from the RCSB PDB.

By default each enzyme's first biological assembly is downloaded (the oligomer
the enzyme works as, e.g. the aldolase tetramer). Use --asymmetric-unit to get
the deposited asymmetric unit instead. Files are saved gzipped, which
cellscape reads directly, and a manifest.csv is written alongside them for
render_all.py.

Usage:
    python download_enzymes.py [--out glycolysis_cif] [--asymmetric-unit]
"""
import argparse
import csv
import sys
import urllib.error
import urllib.request
from pathlib import Path

# (step, enzyme, PDB ID, what the entry is). All human except where noted.
ENZYMES = [
    (1, "Hexokinase", "1HKB", "human brain hexokinase I"),
    (2, "Glucose-6-phosphate isomerase", "1IAT", "human phosphoglucose isomerase"),
    (3, "Phosphofructokinase-1", "4XYJ", "human platelet PFK-1 with ATP"),
    (4, "Aldolase", "4ALD", "human muscle aldolase A with FBP"),
    (5, "Triosephosphate isomerase", "1HTI", "human TIM"),
    (6, "Glyceraldehyde-3-phosphate dehydrogenase", "1U8F", "human placental GAPDH"),
    (7, "Phosphoglycerate kinase", "2WZB", "human PGK1, closed conformation"),
    (8, "Phosphoglycerate mutase", "1YFK", "human PGAM (B type)"),
    (9, "Enolase", "2PSN", "human enolase 1"),
    (10, "Pyruvate kinase", "3BJF", "human PKM2"),
]

DOWNLOAD_URL = "https://files.rcsb.org/download/{name}.cif.gz"
USER_AGENT = "cellscape-example/1.0"


def fetch(url, dest):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        data = response.read()
    dest.write_bytes(data)
    return len(data)


def download(pdb_id, out_dir, assembly=True):
    """Download one entry, preferring biological assembly 1. Returns the saved path."""
    names = [f"{pdb_id}-assembly1", pdb_id] if assembly else [pdb_id]
    for name in names:
        dest = out_dir / f"{name.lower()}.cif.gz"
        if dest.exists():
            print(f"  {dest.name} already present, skipping")
            return dest
        try:
            size = fetch(DOWNLOAD_URL.format(name=name), dest)
        except urllib.error.HTTPError as e:
            if e.code == 404 and name != names[-1]:
                print(f"  no {name} on RCSB, falling back to the asymmetric unit")
                continue
            raise
        print(f"  saved {dest.name} ({size / 1024:.0f} KB)")
        return dest


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="glycolysis_cif", help="Directory to save files into (default: %(default)s)")
    parser.add_argument("--asymmetric-unit", action="store_true",
                        help="Download the deposited asymmetric unit instead of biological assembly 1")
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    failed = []
    for step, enzyme, pdb_id, description in ENZYMES:
        print(f"[{step}/10] {enzyme} ({pdb_id}: {description})")
        try:
            path = download(pdb_id, out_dir, assembly=not args.asymmetric_unit)
        except (urllib.error.URLError, OSError) as e:
            print(f"  FAILED: {e}", file=sys.stderr)
            failed.append(pdb_id)
            continue
        rows.append({"step": step, "enzyme": enzyme, "pdb_id": pdb_id,
                     "description": description, "file": path.name})

    manifest = out_dir / "manifest.csv"
    with open(manifest, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["step", "enzyme", "pdb_id", "description", "file"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {manifest} ({len(rows)} entries)")

    if failed:
        sys.exit(f"Failed to download: {', '.join(failed)}")


if __name__ == "__main__":
    main()
