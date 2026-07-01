"""
Independent 2.5D ideal-MHD solver for the MHD shock-cloud benchmark.

Scheme: 2nd-order MUSCL (minmod) + HLL Riemann solver + RK2 time stepping,
with a first-order fallback for nonphysical reconstructed interface states.
Grid: cell-centered, 256 x 256, [0,1] x [0,1].
Boundaries:
  - x = 0  (X1-beg): outflow, zero-gradient
  - x = 1  (X1-end): userdef, pre-shock state
  - y = 0  (X2-beg): outflow, zero-gradient
  - y = 1  (X2-end): outflow, zero-gradient

Same initial conditions as PLUTO Test_Problems/MHD/Shock_Cloud/setup 01
(see ../../inputs/init.c).  Positivity floors repair density/pressure after
RK stages; this keeps the run finite but is not strictly conservative.

This is intentionally a *simple* reference solver: no constrained
transport, no 8-wave source, no entropy fix, no AMR.  Divergence of B
is therefore not machine-zero; the report discusses the consequence.

Usage:
    python -m python.solver                       # defaults: 256, 0.06
    python -m python.solver --nx 128 --tstop 0.03
    python -m python.solver --out output/python/data.npz
    python -m python.solver --snapshot-count 21   # add rho time-series
"""

from __future__ import annotations

import argparse
import json
import os
import time
from dataclasses import dataclass

import numpy as np


# ============================================================
# Problem constants (from inputs/init.c and inputs/pluto.ini)
# ============================================================
GAMMA = 5.0 / 3.0
GM1 = GAMMA - 1.0

B_POST = 2.1826182
B_PRE = 0.56418958
RADIUS = 0.15
CLOUD_RHO = 10.0
CLOUD_X0 = 0.8
CLOUD_Y0 = 0.5

# Post-shock (left, x < 0.6)
RHO_L = 3.86859
VX_L = 0.0
P_L = 167.345
BY_L = +B_POST
BZ_L = -B_POST

# Pre-shock (right, x >= 0.6)
RHO_R = 1.0
VX_R = -11.2536
P_R = 1.0
BY_R = +B_PRE
BZ_R = +B_PRE

# Conservative variable indices
NVAR = 8
IRHO, IVX, IVY, IVZ, IBX, IBY, IBZ, IEN = range(NVAR)
RHO_FLOOR = 1e-3
PRS_FLOOR = 1e-5

# Primitive-to-primitive (and friends) axis conventions:
#  - ux = vx,  uy = vy
#  - For the HLL Riemann solver the "normal" direction is either x or y,
#    picked by a scalar index dir_idx (0 = x, 1 = y).
NORMAL = {0: "x", 1: "y"}


# ============================================================
# Initial condition
# ============================================================
def setup_ic(NX: int, NY: int, x1max: float = 1.0, x2max: float = 1.0):
    """Return primitive variables on cell centers, each shape (NX, NY)."""
    dx = x1max / NX
    dy = x2max / NY
    x = (np.arange(NX) + 0.5) * dx
    y = (np.arange(NY) + 0.5) * dy
    X, Y = np.meshgrid(x, y, indexing="ij")

    post = X < 0.6
    rho = np.where(post, RHO_L, RHO_R).astype(np.float64)
    vx = np.where(post, VX_L, VX_R).astype(np.float64)
    p = np.where(post, P_L, P_R).astype(np.float64)
    By = np.where(post, BY_L, BY_R).astype(np.float64)
    Bz = np.where(post, BZ_L, BZ_R).astype(np.float64)

    vy = np.zeros_like(X)
    vz = np.zeros_like(X)
    Bx = np.zeros_like(X)

    r = np.sqrt((X - CLOUD_X0) ** 2 + (Y - CLOUD_Y0) ** 2)
    rho = np.where(r < RADIUS, CLOUD_RHO, rho)

    return rho, vx, vy, vz, p, Bx, By, Bz


def prim_to_cons(rho, vx, vy, vz, p, Bx, By, Bz):
    """Primitive (each NX,NY) -> conservative (NVAR, NX, NY)."""
    B2 = Bx * Bx + By * By + Bz * Bz
    v2 = vx * vx + vy * vy + vz * vz
    E = p / GM1 + 0.5 * rho * v2 + 0.5 * B2
    return np.stack(
        [rho, rho * vx, rho * vy, rho * vz, Bx, By, Bz, E], axis=0
    )


def cons_to_prim(U):
    """Conservative (NVAR, NX, NY) -> tuple of 8 primitive (NX, NY) arrays.

    Applies a density and pressure floor for stability in shocks.
    """
    rho = np.maximum(U[IRHO], RHO_FLOOR)
    inv = 1.0 / rho
    vx = U[IVX] * inv
    vy = U[IVY] * inv
    vz = U[IVZ] * inv
    Bx = U[IBX]
    By = U[IBY]
    Bz = U[IBZ]
    E = U[IEN]
    B2 = Bx * Bx + By * By + Bz * Bz
    v2 = vx * vx + vy * vy + vz * vz
    p = GM1 * (E - 0.5 * rho * v2 - 0.5 * B2)
    p = np.maximum(p, PRS_FLOOR)
    return rho, vx, vy, vz, p, Bx, By, Bz


def _pressure_from_cons(U: np.ndarray) -> np.ndarray:
    """Raw gas pressure from conservative state, without floors."""
    rho = U[IRHO]
    rho_safe = np.where(rho > 0.0, rho, np.nan)
    vx = U[IVX] / rho_safe
    vy = U[IVY] / rho_safe
    vz = U[IVZ] / rho_safe
    Bx = U[IBX]
    By = U[IBY]
    Bz = U[IBZ]
    v2 = vx * vx + vy * vy + vz * vz
    B2 = Bx * Bx + By * By + Bz * Bz
    return GM1 * (U[IEN] - 0.5 * rho_safe * v2 - 0.5 * B2)


def _physical_state_mask(U: np.ndarray) -> np.ndarray:
    """True where conservative states are usable without floor repair."""
    rho = U[IRHO]
    p = _pressure_from_cons(U)
    return (
        np.isfinite(rho)
        & np.isfinite(p)
        & (rho > RHO_FLOOR)
        & (p > PRS_FLOOR)
    )


def _fallback_bad_interfaces(UL: np.ndarray, UR: np.ndarray,
                             UL0: np.ndarray, UR0: np.ndarray):
    """Replace nonphysical MUSCL states by first-order states interface-wise."""
    bad = ~(_physical_state_mask(UL) & _physical_state_mask(UR))
    if np.any(bad):
        UL = UL.copy()
        UR = UR.copy()
        UL[:, bad] = UL0[:, bad]
        UR[:, bad] = UR0[:, bad]
    return UL, UR


def _enforce_physical_floors(U: np.ndarray) -> np.ndarray:
    """Repair cell states after an update by enforcing rho and pressure floors."""
    U = U.copy()
    U[IRHO] = np.maximum(U[IRHO], RHO_FLOOR)

    rho = U[IRHO]
    vx = U[IVX] / rho
    vy = U[IVY] / rho
    vz = U[IVZ] / rho
    Bx = U[IBX]
    By = U[IBY]
    Bz = U[IBZ]
    kinetic = 0.5 * rho * (vx * vx + vy * vy + vz * vz)
    magnetic = 0.5 * (Bx * Bx + By * By + Bz * Bz)
    p = GM1 * (U[IEN] - kinetic - magnetic)
    bad_p = (~np.isfinite(p)) | (p < PRS_FLOOR)
    if np.any(bad_p):
        U[IEN, bad_p] = PRS_FLOOR / GM1 + kinetic[bad_p] + magnetic[bad_p]
    return U


# ============================================================
# HLL Riemann solver and fluxes
# ============================================================
def _fast_magnetosonic(cs2: np.ndarray, ca2: np.ndarray,
                       Bn: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Fast magnetosonic speed squared (Toro 2009, eq. 4.66).

    c_f^2 = 0.5 * (c_s^2 + c_A^2 + sqrt((c_s^2 + c_A^2)^2 - 4 c_s^2 Bn^2/rho))
    """
    s = cs2 + ca2
    disc = np.maximum(s * s - 4.0 * cs2 * Bn * Bn / rho, 0.0)
    return 0.5 * (s + np.sqrt(disc))


def flux_dir(prim: np.ndarray, dir_idx: int) -> np.ndarray:
    """Compute flux F (or G) in direction dir_idx.  prim is (NVAR, ...)."""
    rho = prim[0]
    vx = prim[1]
    vy = prim[2]
    vz = prim[3]
    p = prim[4]
    Bx = prim[5]
    By = prim[6]
    Bz = prim[7]

    B2 = Bx * Bx + By * By + Bz * Bz
    pt = p + 0.5 * B2
    vdotB = vx * Bx + vy * By + vz * Bz
    E = p / GM1 + 0.5 * rho * (vx * vx + vy * vy + vz * vz) + 0.5 * B2

    if dir_idx == 0:
        F0 = rho * vx
        F1 = rho * vx * vx + pt - Bx * Bx
        F2 = rho * vx * vy - Bx * By
        F3 = rho * vx * vz - Bx * Bz
        F4 = Bx * vx
        F5 = By * vx - Bx * vy
        F6 = Bz * vx - Bx * vz
        F7 = (E + pt) * vx - vdotB * Bx
    else:
        F0 = rho * vy
        F1 = rho * vy * vx - By * Bx
        F2 = rho * vy * vy + pt - By * By
        F3 = rho * vy * vz - By * Bz
        F4 = Bx * vy - By * vx
        F5 = By * vy
        F6 = Bz * vy - By * vz
        F7 = (E + pt) * vy - vdotB * By

    return np.stack([F0, F1, F2, F3, F4, F5, F6, F7], axis=0)


def _prim_to_cons_batch(prim: np.ndarray) -> np.ndarray:
    """prim (..., 8) -> cons (..., 8)."""
    rho = prim[..., 0]
    vx = prim[..., 1]
    vy = prim[..., 2]
    vz = prim[..., 3]
    p = prim[..., 4]
    Bx = prim[..., 5]
    By = prim[..., 6]
    Bz = prim[..., 7]
    B2 = Bx * Bx + By * By + Bz * Bz
    v2 = vx * vx + vy * vy + vz * vz
    E = p / GM1 + 0.5 * rho * v2 + 0.5 * B2
    return np.stack(
        [rho, rho * vx, rho * vy, rho * vz, Bx, By, Bz, E], axis=-1
    )


def _cons_to_prim_batch(U: np.ndarray) -> np.ndarray:
    """cons (..., 8) -> prim (..., 8)."""
    rho = np.maximum(U[..., 0], RHO_FLOOR)
    inv = 1.0 / rho
    vx = U[..., 1] * inv
    vy = U[..., 2] * inv
    vz = U[..., 3] * inv
    Bx = U[..., 4]
    By = U[..., 5]
    Bz = U[..., 6]
    E = U[..., 7]
    B2 = Bx * Bx + By * By + Bz * Bz
    v2 = vx * vx + vy * vy + vz * vz
    p = GM1 * (E - 0.5 * rho * v2 - 0.5 * B2)
    p = np.maximum(p, PRS_FLOOR)
    return np.stack([rho, vx, vy, vz, p, Bx, By, Bz], axis=-1)


def hll_flux(UL: np.ndarray, UR: np.ndarray,
             dir_idx: int) -> np.ndarray:
    """HLL/Einfeldt (Rusanov) flux at each interface.

    UL, UR are conservative with NVAR as axis 0.  Uses Einfeldt wave-speed
    estimates
        S_L = -S_max,  S_R = +S_max
        S_max = max(|u_n|_L + c_fL, |u_n|_R + c_fR)
    which is the most diffusive variant of the HLL family.  It is less
    accurate than Roe or HLLD at contact discontinuities but is the most
    stable for strong MHD shocks at moderate resolution, which is what we
    need to make a run at 256^2 with simple numerics.
    """
    rhoL, vxL, vyL, vzL, prL, BxL, ByL, BzL = cons_to_prim(UL)
    rhoR, vxR, vyR, vzR, prR, BxR, ByR, BzR = cons_to_prim(UR)

    unL = vxL if dir_idx == 0 else vyL
    unR = vxR if dir_idx == 0 else vyR

    cs2L = GAMMA * prL / rhoL
    cs2R = GAMMA * prR / rhoR
    ca2L = (BxL * BxL + ByL * ByL + BzL * BzL) / rhoL
    ca2R = (BxR * BxR + ByR * ByR + BzR * BzR) / rhoR
    BnL = BxL if dir_idx == 0 else ByL
    BnR = BxR if dir_idx == 0 else ByR
    cf2L = _fast_magnetosonic(cs2L, ca2L, BnL, rhoL)
    cf2R = _fast_magnetosonic(cs2R, ca2R, BnR, rhoR)
    cfL = np.sqrt(np.maximum(cf2L, 0.0))
    cfR = np.sqrt(np.maximum(cf2R, 0.0))

    # Einfeldt (Rusanov) wave speeds: S_L = -S_max, S_R = +S_max
    Smax = np.maximum(np.abs(unL) + cfL, np.abs(unR) + cfR)
    SL = -Smax
    SR = +Smax

    primL = np.stack([rhoL, vxL, vyL, vzL, prL, BxL, ByL, BzL], axis=0)
    primR = np.stack([rhoR, vxR, vyR, vzR, prR, BxR, ByR, BzR], axis=0)
    FL = flux_dir(primL, dir_idx)
    FR = flux_dir(primR, dir_idx)

    dSR = SR - SL  # = 2 * Smax, always positive
    Fhll = (SR * FL - SL * FR + SL * SR * (UR - UL)) / dSR
    return Fhll


# ============================================================
# MUSCL reconstruction with minmod limiter
# ============================================================
def _minmod(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.where(a * b > 0, np.where(np.abs(a) < np.abs(b), a, b), 0.0)


def _apply_bcs(U: np.ndarray) -> np.ndarray:
    """Pad U with 1 ghost cell on each side and fill boundary values.

    U has shape (NVAR, NX, NY).  Returns (NVAR, NX+2, NY+2).
    """
    NX, NY = U.shape[1], U.shape[2]
    Up = np.empty((NVAR, NX + 2, NY + 2), dtype=U.dtype)
    Up[:, 1:-1, 1:-1] = U

    pre_shock = prim_to_cons_pre_shock()  # (NVAR,)

    # x-boundaries
    # X1-beg (x=0): outflow, zero-gradient => ghost = first interior cell
    Up[:, 0, 1:-1] = U[:, 0, :]
    # X1-end (x=1): userdef, pre-shock state
    Up[:, -1, 1:-1] = pre_shock[:, None]

    # y-boundaries
    # X2-beg/X2-end: outflow, zero-gradient
    Up[:, 1:-1, 0] = U[:, :, 0]
    Up[:, 1:-1, -1] = U[:, :, -1]

    # corners: x=0 sides use interior; x=1 sides use pre-shock
    Up[:, 0, 0] = U[:, 0, 0]
    Up[:, 0, -1] = U[:, 0, -1]
    Up[:, -1, 0] = pre_shock
    Up[:, -1, -1] = pre_shock
    return Up


def prim_to_cons_pre_shock() -> np.ndarray:
    """Pre-shock state as a conservative vector of length NVAR."""
    rho = RHO_R
    vx = VX_R
    vy = 0.0
    vz = 0.0
    Bx = 0.0
    By = BY_R
    Bz = BZ_R
    p = P_R
    B2 = Bx * Bx + By * By + Bz * Bz
    v2 = vx * vx + vy * vy + vz * vz
    E = p / GM1 + 0.5 * rho * v2 + 0.5 * B2
    return np.array([rho, rho * vx, rho * vy, rho * vz, Bx, By, Bz, E])


def muscl_fluxes(U: np.ndarray, dx: float, dy: float):
    """Compute the flux divergence dF/dx + dG/dy at each interior cell.

    Returns (dFx, dGy, src) of shape (NVAR, NX, NY) such that the explicit
    update is

        U_new = U - dt * (dFx + dGy) + dt * src
    """
    Up = _apply_bcs(U)
    NX, NY = U.shape[1], U.shape[2]

    # ---- x-direction ----
    # One-sided differences over the padded array, then minmod.
    # d_x[j] = U[j+1] - U[j] for padded index j = 0..NX+1 (NX+2 values).
    d_x = Up[:, 1:, 1:-1] - Up[:, :-1, 1:-1]            # (NVAR, NX+1, NY)
    dL_x = d_x[:, :-1, :]                                # U[j] - U[j-1], j=1..NX
    dR_x = d_x[:, 1:, :]                                 # U[j+1] - U[j], j=1..NX
    slope_x = _minmod(dL_x, dR_x)                        # (NVAR, NX, NY), slope at interior cell j (padded j+1)

    # Slope at padded indices 2..NX+1 (NX values), for UR.
    # Take slope_x[:, 1:, :] (padded 2..NX) and append a zero-slope at the
    # right ghost cell, since the HLL flux at the rightmost interface is
    # already controlled by the pre-shock state in the ghost cell.
    slope_x_right = np.concatenate(
        [slope_x[:, 1:, :], np.zeros((NVAR, 1, NY), dtype=slope_x.dtype)],
        axis=1,
    )                                                    # (NVAR, NX, NY)

    # Interface i+1/2 for i=0..NX-1 (padded index i+1 in [1, NX])
    ULx = Up[:, 1:NX+1, 1:-1] + 0.5 * slope_x
    URx = Up[:, 2:NX+2, 1:-1] - 0.5 * slope_x_right
    ULx, URx = _fallback_bad_interfaces(
        ULx, URx, Up[:, 1:NX+1, 1:-1], Up[:, 2:NX+2, 1:-1]
    )
    Fx = hll_flux(ULx, URx, dir_idx=0)                   # (NVAR, NX, NY)

    # dF_x at interior cell i = Fx[i+1/2] - Fx[i-1/2].
    # Pad Fx on the LEFT ONLY with the outflow boundary flux (= Fx[0]).
    # The right boundary contribution uses Fx[NX-1] (already in Fx).
    Fx_padded = np.concatenate([Fx[:, :1, :], Fx], axis=1)  # (NVAR, NX+1, NY)
    dFx = (Fx_padded[:, 1:, :] - Fx_padded[:, :-1, :]) / dx

    # ---- y-direction (analogous) ----
    d_y = Up[:, 1:-1, 1:] - Up[:, 1:-1, :-1]            # (NVAR, NX, NY+1)
    dL_y = d_y[:, :, :-1]                                # U[j] - U[j-1]
    dR_y = d_y[:, :, 1:]                                 # U[j+1] - U[j]
    slope_y = _minmod(dL_y, dR_y)                        # (NVAR, NX, NY)

    slope_y_right = np.concatenate(
        [slope_y[:, :, 1:], np.zeros((NVAR, NX, 1), dtype=slope_y.dtype)],
        axis=2,
    )                                                    # (NVAR, NX, NY)

    ULy = Up[:, 1:-1, 1:NY+1] + 0.5 * slope_y
    URy = Up[:, 1:-1, 2:NY+2] - 0.5 * slope_y_right
    ULy, URy = _fallback_bad_interfaces(
        ULy, URy, Up[:, 1:-1, 1:NY+1], Up[:, 1:-1, 2:NY+2]
    )
    Gy = hll_flux(ULy, URy, dir_idx=1)                   # (NVAR, NX, NY)

    Gy_padded = np.concatenate([Gy[:, :, :1], Gy], axis=2)  # (NVAR, NX, NY+1)
    dGy = (Gy_padded[:, :, 1:] - Gy_padded[:, :, :-1]) / dy

    return dFx, dGy, np.zeros_like(U)  # 8-wave source term disabled (unstable for this problem)


# ============================================================
# Time stepping
# ============================================================
def max_wave_speed(U: np.ndarray) -> float:
    """Global max fast-magnetosonic speed + flow speed across the domain."""
    rho, vx, vy, vz, p, Bx, By, Bz = cons_to_prim(U)
    cs2 = GAMMA * p / rho
    ca2 = (Bx * Bx + By * By + Bz * Bz) / rho
    Bn_x = Bx
    Bn_y = By
    cf2_x = _fast_magnetosonic(cs2, ca2, Bn_x, rho)
    cf2_y = _fast_magnetosonic(cs2, ca2, Bn_y, rho)
    cf = np.sqrt(np.maximum(np.maximum(cf2_x, cf2_y), 0.0))
    speed = np.maximum(np.abs(vx) + cf, np.abs(vy) + cf)
    return float(np.max(speed))


@dataclass
class SolverConfig:
    nx: int = 256
    ny: int = 256
    tstop: float = 0.06
    cfl: float = 0.3
    first_dt: float = 1e-5
    out: str = "output/python/data.npz"
    progress_every: int = 50
    s_max_abort: float = 1e10
    snapshot_count: int = 0
    snapshot_vars: tuple[str, ...] = ("rho",)


PRIMITIVE_NAMES = ("rho", "vx", "vy", "vz", "prs", "Bx", "By", "Bz")


def _primitive_dict(U: np.ndarray) -> dict[str, np.ndarray]:
    return dict(zip(PRIMITIVE_NAMES, cons_to_prim(U)))


def run(cfg: SolverConfig) -> dict:
    """Integrate to tstop and save the final state.  Returns a stats dict."""
    NX, NY = cfg.nx, cfg.ny
    dx = 1.0 / NX
    dy = 1.0 / NY

    rho, vx, vy, vz, p, Bx, By, Bz = setup_ic(NX, NY)
    U = prim_to_cons(rho, vx, vy, vz, p, Bx, By, Bz)

    t = 0.0
    step = 0
    dt = cfg.first_dt
    t0 = time.time()
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
        prim = _primitive_dict(U)
        for name in cfg.snapshot_vars:
            if name not in prim:
                valid = ", ".join(PRIMITIVE_NAMES)
                raise ValueError(f"unknown snapshot variable {name!r}; valid: {valid}")
            snapshots[name].append(np.asarray(prim[name], dtype=np.float32).copy())
        snapshot_t.append(float(t))

    # Estimate initial wave speed and adjust first dt if needed
    Smax = max_wave_speed(U)
    dt = min(cfg.first_dt, cfg.cfl * min(dx, dy) / Smax)
    t_last_print = t0

    if snapshot_targets.size and snapshot_targets[0] <= 0.0:
        save_snapshot()
        snapshot_next = 1

    while t < cfg.tstop:
        # RK2:
        #   k1 = -RHS(U)
        #   U1  = U + dt*k1
        #   k2  = -RHS(U1)
        #   U   = U + 0.5*dt*(k1 + k2)
        dFx, dGy, src = muscl_fluxes(U, dx, dy)
        k1 = -(dFx + dGy) + src
        U1 = _enforce_physical_floors(U + dt * k1)

        dFx1, dGy1, src1 = muscl_fluxes(U1, dx, dy)
        k2 = -(dFx1 + dGy1) + src1
        U = _enforce_physical_floors(U + 0.5 * dt * (k1 + k2))

        t += dt
        step += 1

        while (
            snapshot_next < snapshot_targets.size
            and t + 1e-12 >= snapshot_targets[snapshot_next]
        ):
            save_snapshot()
            snapshot_next += 1

        # Adapt dt
        Smax = max_wave_speed(U)
        if not np.isfinite(Smax) or Smax > cfg.s_max_abort:
            print(
                f"  ABORT at step {step}: Smax = {Smax} (numerical instability). "
                f"Saving snapshot.",
                flush=True,
            )
            break
        dt_new = cfg.cfl * min(dx, dy) / Smax
        dt = min(dt_new, cfg.tstop - t)  # don't overshoot

        if step % cfg.progress_every == 0 or t >= cfg.tstop:
            now = time.time()
            print(
                f"  step {step:5d}  t={t:.5f}  dt={dt:.3e}  "
                f"Smax={Smax:.2f}  wall={now - t0:.1f}s",
                flush=True,
            )
            t_last_print = now

    wall = time.time() - t0
    print(f"  done: {step} steps in {wall:.2f} s")

    # Save
    rho, vx, vy, vz, p, Bx, By, Bz = cons_to_prim(U)
    x = (np.arange(NX) + 0.5) * dx
    y = (np.arange(NY) + 0.5) * dy
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
        "nx": NX, "ny": NY,
        "snapshot_count": len(snapshot_t),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--nx", type=int, default=256)
    ap.add_argument("--ny", type=int, default=256)
    ap.add_argument("--tstop", type=float, default=0.06)
    ap.add_argument("--cfl", type=float, default=0.3)
    ap.add_argument("--first-dt", type=float, default=1e-5)
    ap.add_argument("--out", type=str, default="output/python/data.npz")
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
    print(f"python/solver.py:  NX={cfg.nx} NY={cfg.ny} tstop={cfg.tstop} "
          f"CFL={cfg.cfl} -> {cfg.out}")
    stats = run(cfg)
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
