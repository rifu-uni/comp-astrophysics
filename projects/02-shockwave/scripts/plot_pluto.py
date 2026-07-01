"""
pyPLUTO plots for a single Shock_Cloud setup.

Usage:
    python3 scripts/plot_pluto.py --setup 01 --out plots/01

Reads data.0001.dbl from output/pluto/<setup>/ and produces one PNG per
physical variable (rho, vx, vy, vz, Bx, By, Bz, p) plus a 2x4 summary panel.
"""

import argparse
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pyPLUTO as pp


VARS = (
    ("rho", r"$\rho$", "viridis"),
    ("vx1", r"$v_x$", "RdBu_r"),
    ("vx2", r"$v_y$", "RdBu_r"),
    ("vx3", r"$v_z$", "RdBu_r"),
    ("Bx1", r"$B_x$", "RdBu_r"),
    ("Bx2", r"$B_y$", "RdBu_r"),
    ("Bx3", r"$B_z$", "RdBu_r"),
    ("prs", r"$p$", "inferno"),
)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--setup", required=True,
                    help="setup number, e.g. 01, 02, 09")
    ap.add_argument("--out", required=True,
                    help="output directory for PNGs")
    ap.add_argument("--root", default=".",
                    help="project root containing output/pluto/<setup>/")
    ap.add_argument("--nout", type=int, default=-1,
                    help="output number to plot (default: -1 = final)")
    args = ap.parse_args()

    path = os.path.join(args.root, "output", "pluto", args.setup)
    if not os.path.isdir(path):
        print(f"error: no PLUTO output at {path}", file=sys.stderr)
        return 1

    os.makedirs(args.out, exist_ok=True)
    print(f"Loading {path} (nout={args.nout})")
    D = pp.Load(path=path, nout=args.nout, datatype="dbl",
                var=[v for v, _, _ in VARS], text=False)

    x1, x2 = D.x1, D.x2
    t = float(D.ntime)

    for var, label, cmap in VARS:
        arr = getattr(D, var)
        vmin, vmax = np.percentile(arr, [1, 99])
        if cmap.endswith("_r"):
            vmin, vmax = -max(abs(vmin), abs(vmax)), max(abs(vmin), abs(vmax))
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.pcolormesh(x1, x2, arr.T, vmin=vmin, vmax=vmax,
                           cmap=cmap, shading="auto")
        ax.set_aspect("equal")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel("x"); ax.set_ylabel("y")
        ax.set_title(f"Shock_Cloud setup {args.setup}: {label} at t={t:.3f}")
        fig.colorbar(im, ax=ax, label=label)
        fig.tight_layout()
        p = os.path.join(args.out, f"{var}.png")
        fig.savefig(p, dpi=140); plt.close(fig)
        print(f"  wrote {p}")

    # Summary panel in the report order: rho, velocity components,
    # magnetic components, pressure.
    fig, axes = plt.subplots(2, 4, figsize=(18, 8))
    for ax, var in zip(axes.flat, [v for v, _, _ in VARS]):
        arr = getattr(D, var)
        vmin, vmax = np.percentile(arr, [1, 99])
        cmap = next(c for v, _, c in VARS if v == var)
        if cmap.endswith("_r"):
            vmin, vmax = -max(abs(vmin), abs(vmax)), max(abs(vmin), abs(vmax))
        im = ax.pcolormesh(x1, x2, arr.T, vmin=vmin, vmax=vmax,
                           cmap=cmap, shading="auto")
        ax.set_aspect("equal")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_title(next(l for v, l, _ in VARS if v == var))
        ax.set_xlabel("x")
    axes[0, 0].set_ylabel("y"); axes[1, 0].set_ylabel("y")
    fig.suptitle(f"Shock_Cloud setup {args.setup} at t={t:.3f}")
    fig.tight_layout()
    p = os.path.join(args.out, "summary.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    print(f"  wrote {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
