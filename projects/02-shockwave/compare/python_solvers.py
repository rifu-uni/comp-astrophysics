"""
Compare the baseline Python solver against the CT Python solver.

The script consumes two solver .npz files with the shared output contract
and writes a compact norm table plus direct difference plots.

Usage:
    python -m compare.python_solvers
    python -m compare.python_solvers --baseline output/python/data.npz \
        --candidate output/python-ct/data.npz --out output/compare-python-solvers
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
    ("vx", r"$v_x$", "RdBu_r"),
    ("vy", r"$v_y$", "RdBu_r"),
    ("vz", r"$v_z$", "RdBu_r"),
    ("Bx", r"$B_x$", "RdBu_r"),
    ("By", r"$B_y$", "RdBu_r"),
    ("Bz", r"$B_z$", "RdBu_r"),
    ("prs", r"$p$", "inferno"),
)
DIFF_PLOTS = ("rho", "prs", "Bx", "By")
DIAGNOSTICS = (
    "mass",
    "pressure_floor_cells",
    "cell_divB_max",
    "cell_divB_mean",
    "face_divB_max",
    "face_divB_mean",
    "pressure_floor_repairs",
    "nsteps",
    "wall_s",
)
PRS_FLOOR = 1e-5


def _load(path: str):
    if not os.path.exists(path):
        raise SystemExit(f"error: solver output not found at {path}")
    return np.load(path)


def _norms(a: np.ndarray, b: np.ndarray):
    diff = b - a
    return (
        float(np.mean(np.abs(diff))),
        float(np.sqrt(np.mean(diff * diff))),
        float(np.max(np.abs(diff))),
    )


def _fmt(x: float) -> str:
    if not np.isfinite(x):
        return "NaN"
    if x == 0:
        return "0.000e+00"
    return f"{x:.3e}"


def _cell_divergence(data) -> tuple[float, float]:
    x = np.asarray(data["x1"])
    y = np.asarray(data["x2"])
    dBx = np.gradient(np.asarray(data["Bx"]), x, axis=0, edge_order=2)
    dBy = np.gradient(np.asarray(data["By"]), y, axis=1, edge_order=2)
    div = dBx + dBy
    return float(np.max(np.abs(div))), float(np.mean(np.abs(div)))


def _diagnostics(data) -> dict[str, float | int | None]:
    dx = float(np.mean(np.diff(data["x1"])))
    dy = float(np.mean(np.diff(data["x2"])))
    cc_max, cc_mean = _cell_divergence(data)
    return {
        "mass": float(np.sum(data["rho"]) * dx * dy),
        "pressure_floor_cells": int(np.count_nonzero(data["prs"] <= PRS_FLOOR * (1.0 + 1e-12))),
        "cell_divB_max": cc_max,
        "cell_divB_mean": cc_mean,
        "face_divB_max": float(data["divB_max"]) if "divB_max" in data else None,
        "face_divB_mean": float(data["divB_mean"]) if "divB_mean" in data else None,
        "pressure_floor_repairs": int(data["pressure_floor_repairs"]) if "pressure_floor_repairs" in data else None,
        "nsteps": int(data["nsteps"]) if "nsteps" in data else None,
        "wall_s": float(data["wall_s"]) if "wall_s" in data else None,
    }


def _diagnostic_value(values: dict[str, float | int | None], key: str) -> str:
    value = values.get(key)
    if value is None:
        return "n/a"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return _fmt(float(value))


def _plot_pair_and_difference(base, cand, out_dir: str) -> None:
    x1 = base["x1"]
    x2 = base["x2"]
    for name, label, cmap in VARS:
        a = base[name]
        b = cand[name]
        if cmap.startswith("RdBu"):
            vmax = float(np.percentile(np.abs(np.concatenate([a.ravel(), b.ravel()])), 99.0))
            vmin = -vmax
        else:
            both = np.concatenate([a.ravel(), b.ravel()])
            vmin = float(np.percentile(both, 1.0))
            vmax = float(np.percentile(both, 99.0))
        if vmax <= vmin:
            vmax = vmin + 1.0

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        for ax, arr, title in zip(axes, (a, b), ("Baseline", "CT")):
            im = ax.pcolormesh(x1, x2, arr.T, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
            ax.set_aspect("equal")
            ax.set_title(f"{title}: {label}")
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(r"$y$")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.tight_layout()
        path = os.path.join(out_dir, f"overlay_{name}.png")
        fig.savefig(path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {path}")

        if name not in DIFF_PLOTS:
            continue
        diff = b - a
        dmax = float(np.percentile(np.abs(diff), 99.5))
        if dmax <= 0.0:
            dmax = 1.0
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.pcolormesh(x1, x2, diff.T, cmap="RdBu_r", vmin=-dmax, vmax=dmax, shading="auto")
        ax.set_aspect("equal")
        ax.set_title(f"CT - baseline: {label}   (max |diff|={np.max(np.abs(diff)):.2g})")
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$y$")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=f"Delta {label}")
        path = os.path.join(out_dir, f"diff_{name}.png")
        fig.tight_layout()
        fig.savefig(path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--baseline", default="output/python/data.npz")
    ap.add_argument("--candidate", default="output/python-ct/data.npz")
    ap.add_argument("--out", default="output/compare-python-solvers")
    args = ap.parse_args()

    base = _load(args.baseline)
    cand = _load(args.candidate)
    os.makedirs(args.out, exist_ok=True)

    for axis in ("x1", "x2"):
        if base[axis].shape != cand[axis].shape or not np.allclose(base[axis], cand[axis]):
            raise SystemExit(f"error: {axis} grid mismatch between Python outputs")
    for name, _label, _cmap in VARS:
        if base[name].shape != cand[name].shape:
            raise SystemExit(f"error: {name} shape mismatch: {base[name].shape} vs {cand[name].shape}")

    rows = [("var", "L1", "L2", "Linf", "Linf_normalized_by_baseline_max")]
    for name, _label, _cmap in VARS:
        l1, l2, linf = _norms(base[name], cand[name])
        base_max = float(np.max(np.abs(base[name])))
        rows.append((name, l1, l2, linf, linf / base_max if base_max > 0 else float("nan")))
    base_diag = _diagnostics(base)
    cand_diag = _diagnostics(cand)

    print()
    print(f"{'var':<8} {'L1':>14} {'L2':>14} {'Linf':>14} {'Linf/base_max':>16}")
    print("-" * 72)
    for name, l1, l2, linf, normed in rows[1:]:
        print(f"{name:<8} {_fmt(l1):>14} {_fmt(l2):>14} {_fmt(linf):>14} {_fmt(normed):>16}")

    out_path = os.path.join(args.out, "norms.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("# Python baseline vs Python CT solver\n")
        f.write(f"# baseline = {args.baseline}\n")
        f.write(f"# candidate = {args.candidate}\n")
        f.write(f"# t_baseline = {float(base['t'])}\n")
        f.write(f"# t_candidate = {float(cand['t'])}\n\n")
        f.write(f"{'var':<8} {'L1':>14} {'L2':>14} {'Linf':>14} {'Linf_normalized':>20}\n")
        for name, l1, l2, linf, normed in rows[1:]:
            f.write(f"{name:<8} {l1:14.6e} {l2:14.6e} {linf:14.6e} {normed:20.6e}\n")
        f.write("\n# diagnostics\n")
        f.write(f"{'diagnostic':<24} {'baseline':>18} {'ct':>18}\n")
        for key in DIAGNOSTICS:
            f.write(f"{key:<24} {_diagnostic_value(base_diag, key):>18} {_diagnostic_value(cand_diag, key):>18}\n")
    print(f"\n  wrote {out_path}")

    _plot_pair_and_difference(base, cand, args.out)


if __name__ == "__main__":
    main()
