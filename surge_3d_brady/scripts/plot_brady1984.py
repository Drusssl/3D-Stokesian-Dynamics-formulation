#!/usr/bin/env python3
import csv
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="surge_mpl_"))

import matplotlib.pyplot as plt


def read_xy(path, x_key, y_key):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return [float(r[x_key]) for r in rows], [float(r[y_key]) for r in rows]


def main():
    if len(sys.argv) < 2:
        print("usage: plot_brady1984.py results/brady_A_1000 [out.png]")
        return 2
    run_dir = Path(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else run_dir / "summary.png"
    angular = next(run_dir.glob("*_angular.csv"))
    radial = next(run_dir.glob("*_radial.csv"))
    ydist = next(run_dir.glob("*_y.csv"))

    theta, gtheta = read_xy(angular, "theta_deg", "count_per_sample")
    r, gr = read_xy(radial, "r", "count_per_sample")
    y, gy = read_xy(ydist, "y", "count_per_sample")

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    axes[0].plot(theta, gtheta, marker="o")
    axes[0].set_xlabel("theta [deg]")
    axes[0].set_ylabel("near-contact count/sample")
    axes[0].set_title("Angular")

    axes[1].plot(r, gr)
    axes[1].set_xlim(1.8, 7.0)
    axes[1].set_xlabel("r/a")
    axes[1].set_ylabel("count/sample")
    axes[1].set_title("Radial")

    axes[2].plot(y, gy)
    axes[2].set_xlabel("y/a")
    axes[2].set_ylabel("count/sample")
    axes[2].set_title("Layering")

    fig.tight_layout()
    fig.savefig(out, dpi=180)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

