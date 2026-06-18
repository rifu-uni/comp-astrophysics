import numpy as np
from src.integrators import rk2_heun_step

def ddx(f, dx):
    return (np.roll(f, -1) - np.roll(f, 1)) / (2.0 * dx)

def alfven_ic(x, eps=0.1):
    phi = 2.0 * np.pi * x
    vy = eps * np.sin(phi)
    vz = eps * np.cos(phi)
    By = eps * np.sin(phi)
    Bz = eps * np.cos(phi)
    return vy, vz, By, Bz

def _make_rhs(B0, rho0, mu0, dx, N):
    def rhs_flat(t, y):
        vy = y[0:N]; vz = y[N:2*N]
        By = y[2*N:3*N]; Bz = y[3*N:]
        coeff = B0 / (mu0 * rho0)
        dvy = coeff * ddx(By, dx)
        dvz = coeff * ddx(Bz, dx)
        dBy = B0 * ddx(vy, dx)
        dBz = B0 * ddx(vz, dx)
        return np.concatenate([dvy, dvz, dBy, dBz])
    return rhs_flat

def alfven_analytic(x, t, eps=0.1):
    phi = 2.0 * np.pi * (x - t)
    return eps * np.sin(phi), eps * np.cos(phi)

def solve_alfven(N=128, L=1.0, rho0=1.0, P0=0.1, B0=1.0,
                 eps=0.1, mu0=1.0, CFL=0.4, t_stop=5.0,
                 save_every=8):
    dx = L / N
    x = np.linspace(dx/2, L - dx/2, N)
    vy0, vz0, By0, Bz0 = alfven_ic(x, eps)
    v_A = B0 / np.sqrt(mu0 * rho0)
    dt = CFL * dx / v_A

    f = _make_rhs(B0, rho0, mu0, dx, N)
    y = np.concatenate([vy0, vz0, By0, Bz0])

    snapshots = []
    snap_times = []
    t = 0.0
    step = 0

    while t < t_stop:
        dt_step = min(dt, t_stop - t)
        y = rk2_heun_step(f, t, y, dt_step)
        t += dt_step
        step += 1
        if step % save_every == 0:
            snapshots.append((y[2*N:3*N].copy(), y[3*N:].copy(), t))
            snap_times.append(t)

    vy, vz, By, Bz = y[0:N], y[N:2*N], y[2*N:3*N], y[3*N:]

    return {
        'x': x, 'dx': dx, 'N': N, 'L': L,
        'vy': vy, 'vz': vz, 'By': By, 'Bz': Bz,
        'rho0': rho0, 'P0': P0, 'B0': B0, 'eps': eps, 'mu0': mu0,
        'v_A': v_A, 'dt': dt, 'CFL': CFL, 't': t, 't_stop': t_stop,
        'step': step, 'snapshots': snapshots, 'snap_times': snap_times,
    }

def convergence_test(resolutions=(32, 64, 128, 256, 512),
                     B0=1.0, rho0=1.0, mu0=1.0, eps=0.1,
                     CFL=0.4, t_stop=5.0):
    L1_vals = []
    L2_vals = []
    for N in resolutions:
        dx = 1.0 / N
        x = np.linspace(dx/2, 1.0 - dx/2, N)
        vy0, vz0, By0, Bz0 = alfven_ic(x, eps)
        v_A = B0 / np.sqrt(mu0 * rho0)
        dt = CFL * dx / v_A

        f = _make_rhs(B0, rho0, mu0, dx, N)
        y = np.concatenate([vy0, vz0, By0, Bz0])
        t = 0.0

        while t < t_stop:
            dt_step = min(dt, t_stop - t)
            y = rk2_heun_step(f, t, y, dt_step)
            t += dt_step

        By_a, _ = alfven_analytic(x, t_stop, eps)
        err = np.abs(y[2*N:3*N] - By_a)
        L1_vals.append(np.mean(err))
        L2_vals.append(np.sqrt(np.mean(err**2)))
    return list(resolutions), L1_vals, L2_vals
