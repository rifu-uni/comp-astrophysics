"""
Compute L1, L2, Linf norms between PLUTO and Python outputs for each
physical variable, plus a few derived quantities (|B|, |v|, plasma beta,
magnetization).

Writes output/compare/norms.txt and prints a table to stdout. Python reads
output/python/data.npz by default, or --python-in if supplied.

Usage:
    python -m compare.norms
    python -m compare.norms --setup 01 --out output/compare
    python -m compare.norms --setup 02+09 --out output/compare-02+09
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

import pyPLUTO as pp


VARS = (
    "rho", "vx1", "vx2", "vx3", "Bx1", "Bx2", "Bx3", "prs",
)
PY_KEY = {
    "rho": "rho",
    "vx1": "vx",  "vx2": "vy",  "vx3": "vz",
    "Bx1": "Bx",  "Bx2": "By",  "Bx3": "Bz",
    "prs": "prs",
}
DERIVED = (
    ("|v|",   lambda P, Y: np.sqrt(P.vx1**2 + P.vx2**2)),
    ("|B|",   lambda P, Y: np.sqrt(P.Bx1**2 + P.Bx2**2 + P.Bx3**2)),
    ("beta",  lambda P, Y: 2.0 * P.prs / np.maximum(P.Bx1**2 + P.Bx2**2 + P.Bx3**2, 1e-20)),
)


def load_pluto(setup: str, root: str):
    path = os.path.join(root, "output", "pluto", setup)
    D = pp.Load(path=path, nout=-1, datatype="dbl",
                var=list(VARS), text=False)
    return D


def load_python(root: str, python_in: str | None = None):
    path = python_in or os.path.join(root, "output", "python", "data.npz")
    if not os.path.exists(path):
        raise SystemExit(f"error: Python output not found at {path}")
    D = np.load(path)
    return D


def l_norms(a: np.ndarray, b: np.ndarray):
    diff = a - b
    return (
        float(np.mean(np.abs(diff))),                              # L1 (cell-averaged)
        float(np.sqrt(np.mean(diff * diff))),                       # L2
        float(np.max(np.abs(diff))),                                # Linf
    )


def fmt(x: float) -> str:
    if not np.isfinite(x):
        return "  NaN  "
    if x == 0:
        return "0.0e+00"
    return f"{x:.3e}"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--setup", default="01",
                    help="PLUTO setup to compare against (default: 01)")
    ap.add_argument("--root", default=".",
                    help="project root containing output/pluto/NN/ and output/python/")
    ap.add_argument("--out", default="output/compare",
                    help="output directory for norms.txt")
    ap.add_argument("--python-in", default=None,
                    help="explicit Python .npz output path")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    P = load_pluto(args.setup, args.root)
    Y = load_python(args.root, args.python_in)

    t_pluto = float(P.ntime)
    t_python = float(Y["t"])
    print(f"  t_Python = {t_python:.5f}, t_PLUTO = {t_pluto:.5f} "
          f"(difference = {t_python - t_pluto:+.2e})")

    rows = []
    rows.append(("var", "L1", "L2", "Linf", "Linf_normalized_by_Python_max"))
    for v in VARS:
        pl = P.__getattr__(v)
        py = Y[PY_KEY[v]]
        l1, l2, linf = l_norms(pl, py)
        py_max = float(np.max(np.abs(py)))
        rows.append((v, l1, l2, linf, linf / py_max if py_max > 0 else float("nan")))

    for label, fn in DERIVED:
        pl = fn(P, None)
        # Build the corresponding Python field
        if label == "|v|":
            py = np.sqrt(Y["vx"]**2 + Y["vy"]**2)
        elif label == "|B|":
            py = np.sqrt(Y["Bx"]**2 + Y["By"]**2 + Y["Bz"]**2)
        elif label == "beta":
            B2 = Y["Bx"]**2 + Y["By"]**2 + Y["Bz"]**2
            py = 2.0 * Y["prs"] / np.maximum(B2, 1e-20)
        l1, l2, linf = l_norms(pl, py)
        py_max = float(np.max(np.abs(py)))
        rows.append((label, l1, l2, linf, linf / py_max if py_max > 0 else float("nan")))

    # Pretty-print
    print()
    print(f"{'var':<8} {'L1':>14} {'L2':>14} {'Linf':>14} {'Linf/py_max':>14}")
    print("-" * 66)
    for name, l1, l2, linf, nrm in rows[1:]:
        print(f"{name:<8} {fmt(l1):>14} {fmt(l2):>14} {fmt(linf):>14} {fmt(nrm):>14}")

    # Write to file
    out_path = os.path.join(args.out, "norms.txt")
    with open(out_path, "w") as f:
        f.write(f"# PLUTO ({args.setup}) vs Python solver\n")
        f.write(f"# t_PLUTO = {t_pluto}\n# t_Python = {t_python}\n\n")
        f.write(f"{'var':<8} {'L1':>14} {'L2':>14} {'Linf':>14} {'Linf_normalized':>20}\n")
        for name, l1, l2, linf, nrm in rows[1:]:
            f.write(f"{name:<8} {l1:14.6e} {l2:14.6e} {linf:14.6e} {nrm:20.6e}\n")
    print(f"\n  wrote {out_path}")


if __name__ == "__main__":
    main()
