"""
Postprocess the Python solver output (output/python/data.npz) and produce
the same 8-variable panel as scripts/plot_pluto.py.

Usage:
    python -m python.postprocess
    python -m python.postprocess --in output/python/data.npz --out plots/python
"""

from __future__ import annotations

import argparse
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


VARS = (
    ("rho", r"$\rho$", "viridis"),
    ("vx",  r"$v_x$", "RdBu_r"),
    ("vy",  r"$v_y$", "RdBu_r"),
    ("vz",  r"$v_z$", "RdBu_r"),
    ("Bx",  r"$B_x$", "RdBu_r"),
    ("By",  r"$B_y$", "RdBu_r"),
    ("Bz",  r"$B_z$", "RdBu_r"),
    ("prs", r"$p$",   "inferno"),
)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="inp", default="output/python/data.npz")
    ap.add_argument("--out", default="plots/python")
    args = ap.parse_args()

    if not os.path.exists(args.inp):
        raise SystemExit(f"error: {args.inp} not found; run `make python` first.")

    D = np.load(args.inp)
    x1 = D["x1"]; x2 = D["x2"]; t = float(D["t"])
    os.makedirs(args.out, exist_ok=True)

    # Per-variable PNGs
    for name, label, cmap in VARS:
        arr = D[name]
        # For diverging colormaps centered at zero
        if cmap.startswith("RdBu"):
            vmax = np.percentile(np.abs(arr), 99.0)
            vmin = -vmax
        else:
            vmin = float(np.percentile(arr, 1.0))
            vmax = float(np.percentile(arr, 99.0))
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.pcolormesh(x1, x2, arr.T, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$y$")
        ax.set_title(f"Python solver: {label}  (t = {t:.4f})")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=label)
        out_png = os.path.join(args.out, f"{name}.png")
        fig.savefig(out_png, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {out_png}")

    # 2x4 summary panel in the report order.
    fig, axes = plt.subplots(2, 4, figsize=(18, 8))
    for ax, (name, label, cmap) in zip(axes.flat, VARS):
        arr = D[name]
        if cmap.startswith("RdBu"):
            vmax = np.percentile(np.abs(arr), 99.0)
            vmin = -vmax
        else:
            vmin = float(np.percentile(arr, 1.0))
            vmax = float(np.percentile(arr, 99.0))
        im = ax.pcolormesh(x1, x2, arr.T, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
        ax.set_aspect("equal")
        ax.set_title(label)
        ax.set_xlabel(r"$x$"); ax.set_ylabel(r"$y$")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(f"Python solver: 2.5D MHD shock-cloud (t = {t:.4f})", fontsize=14)
    summary = os.path.join(args.out, "summary.png")
    fig.tight_layout()
    fig.savefig(summary, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {summary}")


if __name__ == "__main__":
    main()
