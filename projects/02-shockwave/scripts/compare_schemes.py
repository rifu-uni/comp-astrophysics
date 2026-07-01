"""
Cross-scheme comparison for Shock_Cloud PLUTO setups (01 / 02 / 09 / 02+09).

Reads the data.*.dbl outputs from output/pluto/{01,02,09,02+09}/, computes L1/L2
differences between every pair of setups for each MHD variable, and produces:

- output/compare-schemes/norms.txt        # table of L1/L2 norms per variable
- output/compare-schemes/density.png       # side-by-side density maps
- output/compare-schemes/diff_01_09.png   # |01 - 09| density map (entropy-fix)
- output/compare-schemes/diff_01_02.png   # |01 - 02| density map (ChTr vs RK2)
- output/compare-schemes/diff_02_02_09.png  # |02 - 02+09| density map

Setups 01 / 02 / 09 / 02+09 share ICs, BCs, grid, and runtime controls:

  01: RK2, linear recon, constrained transport, no entropy fix
  02: ChTr, linear recon, constrained transport, no entropy fix
  09: RK2, linear recon, constrained transport, SELECTIVE entropy fix
  02+09: ChTr, linear recon, constrained transport, SELECTIVE entropy fix

So |01 - 02| isolates the time-stepping effect, and |01 - 09| isolates
the entropy-fix effect. |02 - 02+09| applies the entropy-fix probe on the
Characteristic Tracing branch, while |09 - 02+09| applies the time-stepping
probe with the entropy fix enabled.
"""

import argparse
import itertools
import os
import sys
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pyPLUTO as pp


SETUPS = ("01", "02", "09", "02+09")
SETUP_LABELS = {
    "01": "01 (RK2, no entropy)",
    "02": "02 (ChTr, no entropy)",
    "09": "09 (RK2, SELECTIVE entropy)",
    "02+09": "02+09 (ChTr, SELECTIVE entropy)",
}
VARS = ("rho", "vx1", "vx2", "Bx1", "Bx2", "Bx3", "prs")


def load_setup(setup: str, root: str, nout: int = -1):
    path = os.path.join(root, "output", "pluto", setup)
    if not os.path.isdir(path):
        raise FileNotFoundError(f"No PLUTO output for setup {setup} at {path}")
    D = pp.Load(path=path, nout=nout, datatype="dbl",
                var=list(VARS), text=False)
    return D


def l1_l2(a: np.ndarray, b: np.ndarray):
    diff = a - b
    L1 = float(np.mean(np.abs(diff)))
    L2 = float(np.sqrt(np.mean(diff ** 2)))
    return L1, L2, diff


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="output/compare-schemes",
                    help="output directory")
    ap.add_argument("--root", default=".",
                    help="project root containing output/pluto/")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    print(f"Output: {args.out}")

    print(f"Loading setups {SETUPS} from {args.root}/output/pluto/...")
    D = {s: load_setup(s, args.root) for s in SETUPS}

    # ---- L1 / L2 norms ----
    norm_path = os.path.join(args.out, "norms.txt")
    with open(norm_path, "w") as f:
        f.write(f"Shock_Cloud setup comparison (01 vs 02 vs 09 vs 02+09)\n")
        f.write(f"All setups share the same grid (256x256), tstop=0.06, and ICs.\n")
        f.write(f"Differences are L1 = mean(|a-b|), L2 = sqrt(mean((a-b)^2)).\n\n")
        pairs = list(itertools.combinations(SETUPS, 2))
        f.write(f"{'pair':<10} {'var':<6} {'L1':>14} {'L2':>14}\n")
        f.write("-" * 46 + "\n")
        for a, b in pairs:
            for v in VARS:
                L1, L2, _ = l1_l2(getattr(D[a], v), getattr(D[b], v))
                f.write(f"{a}-{b:<6} {v:<6} {L1:>14.6e} {L2:>14.6e}\n")
            f.write("\n")
    print(f"  wrote {norm_path}")

    # ---- Side-by-side density ----
    rho = {s: getattr(D[s], "rho") for s in SETUPS}
    x1, x2 = D["01"].x1, D["01"].x2
    vmin, vmax = float(min(rho[s].min() for s in SETUPS)), float(max(rho[s].max() for s in SETUPS))

    fig, axes = plt.subplots(2, 2, figsize=(11, 9), sharey=True,
                             constrained_layout=True)
    for ax, s in zip(axes.flat, SETUPS):
        im = ax.pcolormesh(x1, x2, rho[s].T, vmin=vmin, vmax=vmax, shading="auto")
        ax.set_aspect("equal"); ax.set_title(SETUP_LABELS[s])
        ax.set_xlabel("x"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    axes[0, 0].set_ylabel("y")
    axes[1, 0].set_ylabel("y")
    fig.suptitle(r"Shock_Cloud density $\rho$ at $t=0.06$ (PLUTO 01 / 02 / 09 / 02+09)")
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.8, label=r"$\rho$")
    p = os.path.join(args.out, "density.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    print(f"  wrote {p}")

    # ---- Difference maps ----
    for a, b in itertools.combinations(SETUPS, 2):
        label = f"diff_{a.replace('+', '_')}_{b.replace('+', '_')}"
        _, _, diff = l1_l2(rho[a], rho[b])
        vmax_d = float(np.percentile(np.abs(diff), 99.5)) or 1.0
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.pcolormesh(x1, x2, np.abs(diff).T,
                           vmin=0.0, vmax=vmax_d, cmap="magma", shading="auto")
        ax.set_aspect("equal")
        ax.set_title(rf"$|\rho_{{{a}}} - \rho_{{{b}}}|$  ({SETUP_LABELS[a]} vs {SETUP_LABELS[b]})")
        ax.set_xlabel("x"); ax.set_ylabel("y")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        fig.colorbar(im, ax=ax, label=r"$|\Delta\rho|$")
        fig.tight_layout()
        p = os.path.join(args.out, f"{label}.png")
        fig.savefig(p, dpi=140); plt.close(fig)
        print(f"  wrote {p}")

    print("Done.")


if __name__ == "__main__":
    sys.exit(main())
