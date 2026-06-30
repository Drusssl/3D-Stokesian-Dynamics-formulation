#!/usr/bin/env python3
import csv
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="surge_mpl_"))

import matplotlib.pyplot as plt


def main():
    if len(sys.argv) < 2:
        print("usage: plot_trajectory.py two_sphere_trajectory.csv [out.png]")
        return 2
    csv_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else csv_path.with_suffix(".png")

    rows = list(csv.DictReader(open(csv_path, newline="")))
    x1 = [float(r["x1"]) for r in rows]
    y1 = [float(r["y1"]) for r in rows]
    x2 = [float(r["x2"]) for r in rows]
    y2 = [float(r["y2"]) for r in rows]

    plt.figure(figsize=(6, 6))
    plt.plot(x1, y1, label="sphere 1")
    plt.plot(x2, y2, label="sphere 2")
    plt.axis("equal")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=180)
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
