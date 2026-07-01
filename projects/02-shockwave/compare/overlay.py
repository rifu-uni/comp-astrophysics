"""
Compare PLUTO and Python outputs side-by-side.

PLUTO reads data.*.dbl from output/pluto/<setup>/.
Python reads output/python/data.npz by default, or --python-in if supplied.

Outputs go to output/compare/:
  - overlay_rho.png, overlay_Bx.png, overlay_By.png, overlay_p.png  (2x1)
  - diff_rho.png, diff_Bx.png, diff_By.png, diff_p.png  (|PLUTO - Python|)

Usage:
    python -m compare.overlay
    python -m compare.overlay --setup 01 --out output/compare
    python -m compare.overlay --setup 02+09 --out output/compare-02+09
"""

from __future__ import annotations

import argparse
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pyPLUTO as pp


VARS = (
    ("rho", r"$\rho$", "viridis"),
    ("Bx1", r"$B_x$", "RdBu_r"),
    ("Bx2", r"$B_y$", "RdBu_r"),
    ("prs", r"$p$",   "inferno"),
)
PY_KEY = {
    "rho": "rho",
    "vx1": "vx",  "vx2": "vy",  "vx3": "vz",
    "Bx1": "Bx",  "Bx2": "By",  "Bx3": "Bz",
    "prs": "prs",
}


def load_pluto(setup: str, root: str):
    path = os.path.join(root, "output", "pluto", setup)
    D = pp.Load(path=path, nout=-1, datatype="dbl",
                var=[v for v, _, _ in VARS], text=False)
    return D


def load_python(root: str, python_in: str | None = None):
    path = python_in or os.path.join(root, "output", "python", "data.npz")
    if not os.path.exists(path):
        raise SystemExit(f"error: Python output not found at {path}")
    D = np.load(path)
    return D


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--setup", default="01",
                    help="PLUTO setup to compare against (default: 01)")
    ap.add_argument("--root", default=".",
                    help="project root containing output/pluto/NN/ and output/python/")
    ap.add_argument("--out", default="output/compare",
                    help="output directory for figures")
    ap.add_argument("--python-in", default=None,
                    help="explicit Python .npz output path")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    P = load_pluto(args.setup, args.root)
    Y = load_python(args.root, args.python_in)

    # Build the matching grids.  PLUTO has x1, x2 (cell-centered).
    x1p = np.asarray(P.x1); x2p = np.asarray(P.x2)
    x1y = Y["x1"]; x2y = Y["x2"]

    for name, label, cmap in VARS:
        pl_var = P.__getattr__(name)
        py_var = Y[PY_KEY[name]]

        # Both should be on the same (256, 256) grid; verify and use directly.
        if pl_var.shape != py_var.shape:
            print(f"  WARN: {name} shape mismatch PLUTO={pl_var.shape} "
                  f"Python={py_var.shape} -- skipping")
            continue

        # If grids differ, interpolate Python -> PLUTO grid; here we
        # just check that they match closely enough.
        if x1p.shape != x1y.shape or not np.allclose(x1p, x1y, atol=1e-10):
            print(f"  WARN: x1 grids differ for {name}; using Python's grid for plots")
            x1, x2 = x1y, x2y
        else:
            x1, x2 = x1p, x2p

        if cmap.startswith("RdBu"):
            vmax = np.percentile(np.abs(pl_var), 99.0)
            vmin = -vmax
        else:
            vmin = float(np.percentile(pl_var, 1.0))
            vmax = float(np.percentile(pl_var, 99.0))

        # 2x1 overlay: PLUTO | Python
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        for ax, data, ttl in zip(axes, (pl_var, py_var), ("PLUTO", "Python")):
            im = ax.pcolormesh(x1, x2, data.T, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
            ax.set_aspect("equal")
            ax.set_title(f"{ttl}: {label}")
            ax.set_xlabel(r"$x$"); ax.set_ylabel(r"$y$")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        out_png = os.path.join(args.out, f"overlay_{name}.png")
        fig.tight_layout()
        fig.savefig(out_png, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {out_png}")

        # 1x1 difference map
        diff = np.abs(pl_var - py_var)
        dmax = float(np.percentile(diff, 99.5))
        if dmax <= 0:
            dmax = 1.0
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.pcolormesh(x1, x2, diff.T, cmap="magma", vmin=0, vmax=dmax, shading="auto")
        ax.set_aspect("equal")
        ax.set_title(f"|PLUTO - Python|: {label}   (max={diff.max():.2g})")
        ax.set_xlabel(r"$x$"); ax.set_ylabel(r"$y$")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=f"|Δ{label}|")
        out_png = os.path.join(args.out, f"diff_{name}.png")
        fig.tight_layout()
        fig.savefig(out_png, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {out_png}")


if __name__ == "__main__":
    main()
