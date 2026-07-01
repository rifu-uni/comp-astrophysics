"""
Generate reproducible GIF animations for the Shock_Cloud simulations.

Examples:
    python scripts/animate.py --source pluto --setup 01 --out plots/animations/pluto_01_rho.gif
    python scripts/animate.py --source python --out plots/animations/python_rho.gif
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


VAR_STYLE = {
    "rho": (r"$\rho$", "viridis"),
    "vx": (r"$v_x$", "RdBu_r"),
    "vy": (r"$v_y$", "RdBu_r"),
    "vz": (r"$v_z$", "RdBu_r"),
    "vx1": (r"$v_x$", "RdBu_r"),
    "vx2": (r"$v_y$", "RdBu_r"),
    "vx3": (r"$v_z$", "RdBu_r"),
    "Bx": (r"$B_x$", "RdBu_r"),
    "By": (r"$B_y$", "RdBu_r"),
    "Bz": (r"$B_z$", "RdBu_r"),
    "Bx1": (r"$B_x$", "RdBu_r"),
    "Bx2": (r"$B_y$", "RdBu_r"),
    "Bx3": (r"$B_z$", "RdBu_r"),
    "prs": (r"$p$", "inferno"),
}

PLUTO_ALIASES = {
    "vx": "vx1",
    "vy": "vx2",
    "vz": "vx3",
    "Bx": "Bx1",
    "By": "Bx2",
    "Bz": "Bx3",
}
PYTHON_ALIASES = {
    "vx1": "vx",
    "vx2": "vy",
    "vx3": "vz",
    "Bx1": "Bx",
    "Bx2": "By",
    "Bx3": "Bz",
}


def _frame_numbers(path: Path) -> list[int]:
    frames = []
    for file_path in path.glob("data.*.dbl"):
        match = re.fullmatch(r"data\.(\d+)\.dbl", file_path.name)
        if match:
            frames.append(int(match.group(1)))
    return sorted(frames)


def _load_pluto(root: Path, setup: str, var: str):
    try:
        import pyPLUTO as pp
    except ImportError as exc:
        raise SystemExit("pyPLUTO is required for PLUTO animations.") from exc

    pluto_var = PLUTO_ALIASES.get(var, var)
    path = root / "output" / "pluto" / setup
    if not path.is_dir():
        raise SystemExit(f"error: no PLUTO output at {path}; run `make pluto-{setup}`.")

    nouts = _frame_numbers(path)
    if len(nouts) < 2:
        raise SystemExit(
            f"error: found {len(nouts)} PLUTO frame(s) at {path}; "
            "rerun after enabling snapshot output."
        )

    frames = []
    times = []
    x1 = x2 = None
    for nout in nouts:
        data = pp.Load(path=str(path), nout=nout, datatype="dbl",
                       var=[pluto_var], text=False)
        if x1 is None:
            x1 = np.asarray(data.x1)
            x2 = np.asarray(data.x2)
        frames.append(np.asarray(getattr(data, pluto_var), dtype=np.float64))
        times.append(float(data.ntime))
    return x1, x2, np.asarray(times), frames, pluto_var


def _load_python(root: Path, var: str, python_in: str | None = None):
    py_var = PYTHON_ALIASES.get(var, var)
    path = Path(python_in) if python_in else root / "output" / "python" / "data.npz"
    if not path.exists():
        raise SystemExit(f"error: no Python output at {path}; run `make python`.")

    data = np.load(path)
    series_key = f"{py_var}_series"
    if "snapshot_t" not in data or series_key not in data:
        raise SystemExit(
            f"error: {path} has no {series_key}/snapshot_t; "
            f"rerun `make python` with --snapshot-vars {py_var}."
        )

    x1 = np.asarray(data["x1"])
    x2 = np.asarray(data["x2"])
    times = np.asarray(data["snapshot_t"], dtype=np.float64)
    frames = [np.asarray(frame, dtype=np.float64) for frame in data[series_key]]
    if len(frames) < 2:
        raise SystemExit("error: Python output has fewer than two snapshots.")
    return x1, x2, times, frames, py_var


def _extent_from_centers(x1: np.ndarray, x2: np.ndarray) -> tuple[float, float, float, float]:
    dx = float(np.median(np.diff(x1))) if x1.size > 1 else 1.0
    dy = float(np.median(np.diff(x2))) if x2.size > 1 else 1.0
    return (
        float(x1[0] - 0.5 * dx),
        float(x1[-1] + 0.5 * dx),
        float(x2[0] - 0.5 * dy),
        float(x2[-1] + 0.5 * dy),
    )


def _limits(frames: list[np.ndarray], cmap: str) -> tuple[float, float]:
    values = np.concatenate([frame.ravel() for frame in frames])
    if cmap.startswith("RdBu"):
        vmax = float(np.percentile(np.abs(values), 99.0))
        return -vmax, vmax
    return tuple(float(v) for v in np.percentile(values, [1.0, 99.0]))


def _write_gif(
    x1: np.ndarray,
    x2: np.ndarray,
    times: np.ndarray,
    frames: list[np.ndarray],
    var: str,
    title_prefix: str,
    out: Path,
    fps: int,
    dpi: int,
) -> None:
    label, cmap = VAR_STYLE.get(var, (var, "viridis"))
    vmin, vmax = _limits(frames, cmap)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6, 5))
    image = ax.imshow(
        frames[0].T,
        origin="lower",
        extent=_extent_from_centers(x1, x2),
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        interpolation="nearest",
        animated=True,
    )
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    title = ax.set_title("")
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(label)

    def update(idx: int):
        image.set_data(frames[idx].T)
        title.set_text(f"{title_prefix}: {label} at t={times[idx]:.4f}")
        return image, title

    update(0)
    fig.tight_layout()
    anim = FuncAnimation(fig, update, frames=len(frames), interval=1000 / fps)
    anim.save(out, writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=("pluto", "python"), required=True)
    parser.add_argument("--setup", default="01",
                        help="PLUTO setup number for --source pluto")
    parser.add_argument("--var", default="rho",
                        help="variable to animate; default: rho")
    parser.add_argument("--root", default=".",
                        help="project root containing output/")
    parser.add_argument("--out", required=True,
                        help="output GIF path")
    parser.add_argument("--python-in", default=None,
                        help="explicit Python .npz output path")
    parser.add_argument("--title", default=None,
                        help="optional title prefix for the animation")
    parser.add_argument("--fps", type=int, default=8)
    parser.add_argument("--dpi", type=int, default=110)
    args = parser.parse_args()

    root = Path(args.root)
    if args.source == "pluto":
        x1, x2, times, frames, var = _load_pluto(root, args.setup, args.var)
        title = args.title or f"PLUTO setup {args.setup}"
    else:
        x1, x2, times, frames, var = _load_python(root, args.var, args.python_in)
        title = args.title or "Python solver"

    _write_gif(
        x1=x1,
        x2=x2,
        times=times,
        frames=frames,
        var=var,
        title_prefix=title,
        out=Path(args.out),
        fps=args.fps,
        dpi=args.dpi,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
