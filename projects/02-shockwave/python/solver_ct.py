"""
Independent 2.5D ideal-MHD solver with constrained transport.

This module keeps the baseline solver's finite-volume ingredients
(MUSCL/minmod reconstruction, HLL/Rusanov fluxes, RK2 stepping, positivity
fallback) and changes only the magnetic-field representation:

  - conserved fluid variables stay cell-centered;
  - Bx is stored on x-faces with shape (NX+1, NY);
  - By is stored on y-faces with shape (NX, NY+1);
  - cell-centered Bx/By are reconstructed from face averages for fluxes and
    for the public .npz output;
  - face fields are advanced with an edge-centered Ez = vx*By - vy*Bx.

The output contract intentionally matches python.solver so plotting,
animation, and PLUTO comparison scripts can consume either solver.

Usage:
    python -m python.solver_ct
    python -m python.solver_ct --nx 128 --tstop 0.03
    python -m python.solver_ct --out output/python-ct/data.npz
    python -m python.solver_ct --snapshot-count 21
"""

from __future__ import annotations

import argparse
import json
import os
import time
from dataclasses import dataclass

import numpy as np

from .solver import (
    BZ_R,
    GAMMA,
    GM1,
    IBX,
    IBY,
    IEN,
    IRHO,
    NVAR,
    PRIMITIVE_NAMES,
    PRS_FLOOR,
    RHO_FLOOR,
    cons_to_prim,
    max_wave_speed,
    muscl_fluxes,
    prim_to_cons,
    setup_ic,
)


@dataclass
class SolverConfig:
    nx: int = 256
    ny: int = 256
    tstop: float = 0.06
    cfl: float = 0.3
    first_dt: float = 1e-5
    out: str = "output/python-ct/data.npz"
    progress_every: int = 50
    s_max_abort: float = 1e10
    snapshot_count: int = 0
    snapshot_vars: tuple[str, ...] = ("rho",)


def _faces_from_cells(Bx: np.ndarray, By: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Initialize staggered face fields from cell-centered fields."""
    nx, ny = Bx.shape
    Bx_face = np.empty((nx + 1, ny), dtype=Bx.dtype)
    By_face = np.empty((nx, ny + 1), dtype=By.dtype)

    Bx_face[1:-1, :] = 0.5 * (Bx[:-1, :] + Bx[1:, :])
    Bx_face[0, :] = Bx[0, :]
    Bx_face[-1, :] = Bx[-1, :]

    By_face[:, 1:-1] = 0.5 * (By[:, :-1] + By[:, 1:])
    By_face[:, 0] = By[:, 0]
    By_face[:, -1] = By[:, -1]
    return Bx_face, By_face


def _cell_centered_b(
    Bx_face: np.ndarray, By_face: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Average staggered fields back to cell centers."""
    Bx = 0.5 * (Bx_face[:-1, :] + Bx_face[1:, :])
    By = 0.5 * (By_face[:, :-1] + By_face[:, 1:])
    return Bx, By


def _raw_pressure(U: np.ndarray) -> np.ndarray:
    """Gas pressure before pressure flooring."""
    rho = U[IRHO]
    rho_safe = np.where(rho > 0.0, rho, np.nan)
    vx = U[1] / rho_safe
    vy = U[2] / rho_safe
    vz = U[3] / rho_safe
    Bx = U[4]
    By = U[5]
    Bz = U[6]
    kinetic = 0.5 * rho_safe * (vx * vx + vy * vy + vz * vz)
    magnetic = 0.5 * (Bx * Bx + By * By + Bz * Bz)
    return GM1 * (U[IEN] - kinetic - magnetic)


def _enforce_floors(U: np.ndarray) -> tuple[np.ndarray, int]:
    """Repair density/pressure floors and return pressure-repair count."""
    U = U.copy()
    U[IRHO] = np.maximum(U[IRHO], RHO_FLOOR)

    rho = U[IRHO]
    vx = U[1] / rho
    vy = U[2] / rho
    vz = U[3] / rho
    Bx = U[4]
    By = U[5]
    Bz = U[6]
    kinetic = 0.5 * rho * (vx * vx + vy * vy + vz * vz)
    magnetic = 0.5 * (Bx * Bx + By * By + Bz * Bz)
    p = GM1 * (U[IEN] - kinetic - magnetic)
    bad_p = (~np.isfinite(p)) | (p < PRS_FLOOR)
    if np.any(bad_p):
        U[IEN, bad_p] = PRS_FLOOR / GM1 + kinetic[bad_p] + magnetic[bad_p]
    return U, int(np.count_nonzero(bad_p))


def _sync_cell_b(
    U: np.ndarray, Bx_face: np.ndarray, By_face: np.ndarray
) -> tuple[np.ndarray, int]:
    """Copy face-averaged Bx/By into U and repair pressure if needed.

    The finite-volume energy update was computed with the cell-centered
    magnetic field present at the start of the stage.  When CT replaces Bx/By
    with face-averaged values, adjust total energy by the magnetic-energy
    difference so kinetic and thermal energy are preserved by the projection.
    """
    U = U.copy()
    Bx, By = _cell_centered_b(Bx_face, By_face)
    old_b2 = U[IBX] * U[IBX] + U[IBY] * U[IBY] + U[6] * U[6]
    new_b2 = Bx * Bx + By * By + U[6] * U[6]
    U[IBX] = Bx
    U[IBY] = By
    U[IEN] += 0.5 * (new_b2 - old_b2)
    return _enforce_floors(U)


def _primitive_dict_ct(
    U: np.ndarray, Bx_face: np.ndarray, By_face: np.ndarray
) -> dict[str, np.ndarray]:
    U_synced, _ = _sync_cell_b(U, Bx_face, By_face)
    return dict(zip(PRIMITIVE_NAMES, cons_to_prim(U_synced)))


def _edge_emf(
    U: np.ndarray, Bx_face: np.ndarray, By_face: np.ndarray
) -> np.ndarray:
    """Return upwind-stabilized edge-centered Ez with shape (NX+1, NY+1)."""
    del Bx_face, By_face
    rho, vx, vy, _vz, p, Bx, By, Bz = cons_to_prim(U)
    Ez_cell = vx * By - vy * Bx
    Ez_pad = np.pad(Ez_cell, ((1, 1), (1, 1)), mode="edge")
    Ez_center = 0.25 * (
        Ez_pad[:-1, :-1]
        + Ez_pad[1:, :-1]
        + Ez_pad[:-1, 1:]
        + Ez_pad[1:, 1:]
    )

    # A centered EMF preserves the discrete CT divergence but is too sharp for
    # the deliberately simple HLL/RK2 fluid update. Add local LF dissipation to
    # the edge EMF: +alpha_x dBy/2 damps By jumps seen by the x-update and
    # -alpha_y dBx/2 damps Bx jumps seen by the y-update.
    cs2 = GAMMA * p / rho
    ca2 = (Bx * Bx + By * By + Bz * Bz) / rho
    cf = np.sqrt(np.maximum(cs2 + ca2, 0.0))
    speed_x = np.abs(vx) + cf
    speed_y = np.abs(vy) + cf

    By_pad = np.pad(By, ((1, 1), (1, 1)), mode="edge")
    Bx_pad = np.pad(Bx, ((1, 1), (1, 1)), mode="edge")
    sx_pad = np.pad(speed_x, ((1, 1), (1, 1)), mode="edge")
    sy_pad = np.pad(speed_y, ((1, 1), (1, 1)), mode="edge")

    By_left = 0.5 * (By_pad[:-1, :-1] + By_pad[:-1, 1:])
    By_right = 0.5 * (By_pad[1:, :-1] + By_pad[1:, 1:])
    Bx_bottom = 0.5 * (Bx_pad[:-1, :-1] + Bx_pad[1:, :-1])
    Bx_top = 0.5 * (Bx_pad[:-1, 1:] + Bx_pad[1:, 1:])

    ax_left = np.maximum(sx_pad[:-1, :-1], sx_pad[:-1, 1:])
    ax_right = np.maximum(sx_pad[1:, :-1], sx_pad[1:, 1:])
    ay_bottom = np.maximum(sy_pad[:-1, :-1], sy_pad[1:, :-1])
    ay_top = np.maximum(sy_pad[:-1, 1:], sy_pad[1:, 1:])
    alpha_x = np.maximum(ax_left, ax_right)
    alpha_y = np.maximum(ay_bottom, ay_top)

    return (
        Ez_center
        + 0.5 * alpha_x * (By_right - By_left)
        - 0.5 * alpha_y * (Bx_top - Bx_bottom)
    )


def _ct_rhs(
    U: np.ndarray, Bx_face: np.ndarray, By_face: np.ndarray, dx: float, dy: float
) -> tuple[np.ndarray, np.ndarray]:
    """Face-field RHS from the discrete curl of edge-centered Ez."""
    Ez = _edge_emf(U, Bx_face, By_face)
    dBx = -(Ez[:, 1:] - Ez[:, :-1]) / dy
    dBy = +(Ez[1:, :] - Ez[:-1, :]) / dx
    return dBx, dBy


def _discrete_div_b(Bx_face: np.ndarray, By_face: np.ndarray, dx: float, dy: float):
    """Discrete face divergence on cell centers."""
    return (
        (Bx_face[1:, :] - Bx_face[:-1, :]) / dx
        + (By_face[:, 1:] - By_face[:, :-1]) / dy
    )


def _stage_rhs(
    U: np.ndarray, Bx_face: np.ndarray, By_face: np.ndarray, dx: float, dy: float
):
    """Fluid and CT right-hand sides for one RK stage."""
    U_synced, floor_count = _sync_cell_b(U, Bx_face, By_face)
    dFx, dGy, src = muscl_fluxes(U_synced, dx, dy)
    k = -(dFx + dGy) + src
    # Bx/By are advanced by CT, not by the cell-centered finite-volume update.
    k[IBX] = 0.0
    k[IBY] = 0.0
    dBx, dBy = _ct_rhs(U_synced, Bx_face, By_face, dx, dy)
    return U_synced, k, dBx, dBy, floor_count


def run(cfg: SolverConfig) -> dict:
    """Integrate to tstop and save the final CT state."""
    nx, ny = cfg.nx, cfg.ny
    dx = 1.0 / nx
    dy = 1.0 / ny

    rho, vx, vy, vz, p, Bx, By, Bz = setup_ic(nx, ny)
    Bx_face, By_face = _faces_from_cells(Bx, By)
    U = prim_to_cons(rho, vx, vy, vz, p, *_cell_centered_b(Bx_face, By_face), Bz)
    U, initial_floor_repairs = _sync_cell_b(U, Bx_face, By_face)

    t = 0.0
    step = 0
    t0 = time.time()
    pressure_floor_repairs = initial_floor_repairs
    snapshot_targets = (
        np.linspace(0.0, cfg.tstop, cfg.snapshot_count)
        if cfg.snapshot_count > 0
        else np.array([], dtype=np.float64)
    )
    snapshot_next = 0
    snapshot_t: list[float] = []
    snapshots: dict[str, list[np.ndarray]] = {
        name: [] for name in cfg.snapshot_vars
    }

    def save_snapshot() -> None:
        prim = _primitive_dict_ct(U, Bx_face, By_face)
        for name in cfg.snapshot_vars:
            if name not in prim:
                valid = ", ".join(PRIMITIVE_NAMES)
                raise ValueError(f"unknown snapshot variable {name!r}; valid: {valid}")
            snapshots[name].append(np.asarray(prim[name], dtype=np.float32).copy())
        snapshot_t.append(float(t))

    Smax = max_wave_speed(U)
    dt = min(cfg.first_dt, cfg.cfl * min(dx, dy) / Smax)

    if snapshot_targets.size and snapshot_targets[0] <= 0.0:
        save_snapshot()
        snapshot_next = 1

    while t < cfg.tstop:
        U0, k1, dBx1, dBy1, floor_count = _stage_rhs(U, Bx_face, By_face, dx, dy)
        pressure_floor_repairs += floor_count

        U1, stage_floor = _enforce_floors(U0 + dt * k1)
        pressure_floor_repairs += stage_floor
        Bx1 = Bx_face + dt * dBx1
        By1 = By_face + dt * dBy1
        U1, stage_floor = _sync_cell_b(U1, Bx1, By1)
        pressure_floor_repairs += stage_floor

        U1_synced, k2, dBx2, dBy2, floor_count = _stage_rhs(U1, Bx1, By1, dx, dy)
        del U1_synced
        pressure_floor_repairs += floor_count

        U, stage_floor = _enforce_floors(U0 + 0.5 * dt * (k1 + k2))
        pressure_floor_repairs += stage_floor
        Bx_face = Bx_face + 0.5 * dt * (dBx1 + dBx2)
        By_face = By_face + 0.5 * dt * (dBy1 + dBy2)
        U, stage_floor = _sync_cell_b(U, Bx_face, By_face)
        pressure_floor_repairs += stage_floor

        t += dt
        step += 1

        while (
            snapshot_next < snapshot_targets.size
            and t + 1e-12 >= snapshot_targets[snapshot_next]
        ):
            save_snapshot()
            snapshot_next += 1

        Smax = max_wave_speed(U)
        if not np.isfinite(Smax) or Smax > cfg.s_max_abort:
            print(
                f"  ABORT at step {step}: Smax = {Smax} (numerical instability). "
                f"Saving snapshot.",
                flush=True,
            )
            break
        dt_new = cfg.cfl * min(dx, dy) / Smax
        dt = min(dt_new, cfg.tstop - t)

        if step % cfg.progress_every == 0 or t >= cfg.tstop:
            now = time.time()
            divB = _discrete_div_b(Bx_face, By_face, dx, dy)
            print(
                f"  step {step:5d}  t={t:.5f}  dt={dt:.3e}  "
                f"Smax={Smax:.2f}  divBmax={np.max(np.abs(divB)):.2e}  "
                f"wall={now - t0:.1f}s",
                flush=True,
            )

    wall = time.time() - t0
    print(f"  done: {step} steps in {wall:.2f} s")

    U, _ = _sync_cell_b(U, Bx_face, By_face)
    rho, vx, vy, vz, p, Bx, By, Bz = cons_to_prim(U)
    divB = _discrete_div_b(Bx_face, By_face, dx, dy)
    pressure_floor_cells = int(np.count_nonzero(p <= PRS_FLOOR * (1.0 + 1e-12)))

    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    os.makedirs(os.path.dirname(cfg.out) or ".", exist_ok=True)
    payload = {
        "x1": x,
        "x2": y,
        "t": t,
        "rho": rho,
        "vx": vx,
        "vy": vy,
        "vz": vz,
        "prs": p,
        "Bx": Bx,
        "By": By,
        "Bz": Bz,
        "divB_max": float(np.max(np.abs(divB))),
        "divB_mean": float(np.mean(np.abs(divB))),
        "pressure_floor_cells": pressure_floor_cells,
        "pressure_floor_repairs": int(pressure_floor_repairs),
        "nsteps": int(step),
        "wall_s": float(wall),
    }
    if snapshot_t:
        payload["snapshot_t"] = np.asarray(snapshot_t, dtype=np.float64)
        for name, series in snapshots.items():
            payload[f"{name}_series"] = np.stack(series, axis=0)
    np.savez(cfg.out, **payload)

    return {
        "nsteps": step,
        "wall_s": wall,
        "t_final": t,
        "out": cfg.out,
        "nx": nx,
        "ny": ny,
        "snapshot_count": len(snapshot_t),
        "divB_max": payload["divB_max"],
        "divB_mean": payload["divB_mean"],
        "pressure_floor_cells": pressure_floor_cells,
        "pressure_floor_repairs": int(pressure_floor_repairs),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--nx", type=int, default=256)
    ap.add_argument("--ny", type=int, default=256)
    ap.add_argument("--tstop", type=float, default=0.06)
    ap.add_argument("--cfl", type=float, default=0.3)
    ap.add_argument("--first-dt", type=float, default=1e-5)
    ap.add_argument("--out", type=str, default="output/python-ct/data.npz")
    ap.add_argument("--progress-every", type=int, default=50)
    ap.add_argument("--snapshot-count", type=int, default=0,
                    help="number of snapshots to save for animations")
    ap.add_argument("--snapshot-vars", type=str, default="rho",
                    help="comma-separated primitive variables to snapshot")
    args = ap.parse_args()

    snapshot_vars = tuple(
        name.strip() for name in args.snapshot_vars.split(",") if name.strip()
    )
    cfg = SolverConfig(
        nx=args.nx, ny=args.ny, tstop=args.tstop, cfl=args.cfl,
        first_dt=args.first_dt, out=args.out,
        progress_every=args.progress_every,
        snapshot_count=args.snapshot_count,
        snapshot_vars=snapshot_vars,
    )
    print(
        f"python/solver_ct.py:  NX={cfg.nx} NY={cfg.ny} tstop={cfg.tstop} "
        f"CFL={cfg.cfl} -> {cfg.out}"
    )
    stats = run(cfg)
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
