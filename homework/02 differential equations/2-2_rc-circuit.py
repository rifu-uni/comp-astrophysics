#!/usr/bin/env python
# coding: utf-8

# # RC circuit — finding the half-charge time
# Un capacitor de capacitancia $C$ se carga a través de una resistencia $R$ alimentada por un voltaje constante $V$. Partiendo de una carga cero, solucione la ode de carga hacia adelante en el tiempo con RK4. Luego, aplique el método de la secante para encontrar el tiempo exacto $t_m$ en el que la carga alcanza exactamente la mitad de su valor máximo $Q_{\max} = CV$.
# 
# $$
# \frac{dq}{dt} = \frac{(V - q/C)}{R}
# $$
# 
# $q(0) = 0$, $V = 10V$, $R = 1000\Omega$, $C = 1 \times 10^{-3}F$,
# Time span: $[0 - 5 \cdot RC]$, step $h = 0.01s$

# In[6]:


import sys
sys.path.append('../../')
import src.integrators as ode
import src.roots as roots
import numpy as np
import matplotlib.pyplot as plt


# In[7]:


# Parameters
V = 10.0
R = 1000.0
C = 1e-3
Q_max = C * V
h = 0.01
t0 = 0.0
tf = 5.0 * R * C

# ODE definition
def rc_ode(t: float, q: np.ndarray) -> np.ndarray:
    return np.array([(V - q[0]/C) / R])

# Analytical solution
def q_exact(t: float) -> float:
    return C * V * (1 - np.exp(-t / (R * C)))


# In[8]:


# Solve using all four methods
methods = {
    'Euler': ode.euler_step,
    'RK2 (Midpoint)': ode.rk2_midpoint_step,
    'RK4 (Runge)': ode.rk4_runge_step,
    'RK4 (Kutta)': ode.rk4_kutta_step
}

results = {}
for name, step_fn in methods.items():
    t_vals, q_vals = ode.integrate(
        f=rc_ode,
        t_span=(t0, tf),
        y0=np.array([0.0]),
        h=h,
        step_fn=step_fn
    )
    results[name] = (t_vals, q_vals[:, 0])


# In[9]:


# Root finding with Secant Method
def q_residual(t: float) -> float:
    """
    Residual function f(t) = q(t) - Q_max/2.
    Evaluates q(t) by integrating up to t.
    """
    if t <= 0:
        return -Q_max / 2.0
    _, q = ode.integrate(rc_ode, (0.0, t), np.array([0.0]), h, ode.rk4_runge_step)
    return q[-1, 0] - Q_max / 2.0

# Initial guesses for secant method: t=0 and t=RC
t_m = roots.find_root(
    f=q_residual,
    df=None,
    state0=(0.0, R*C),
    step_fn=roots.secant_step,
    tol=1e-8
)

print(f"Exact half-charge time (t_m): {t_m:.5f} s")
print(f"Analytical t_m = -R*C*ln(0.5): {-R*C*np.log(0.5):.5f} s")


# In[10]:


# Plot the result
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

t_ref = results['RK4 (Runge)'][0]
q_ana = q_exact(t_ref)
ax1.plot(t_ref, q_ana, 'k--', label='Analytical', lw=2)

for name, (t, q) in results.items():
    ax1.plot(t, q, label=name)

ax1.axhline(Q_max, linestyle='-', color='gray', alpha=0.5, label='Q_max')
ax1.axhline(Q_max / 2.0, linestyle=':', color='gray', label='Q_max / 2')
ax1.axvline(t_m, linestyle=':', color='red', label=f't_m (RK4 Runge) = {t_m:.4f} s')
ax1.plot([t_m], [Q_max / 2.0], 'ko', zorder=5)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Charge (C)')
ax1.set_title('RC Circuit Charging')
ax1.legend()
ax1.grid(True)

# Error plot
for name, (t, q) in results.items():
    error = np.abs(q - q_exact(t))
    ax2.plot(t[1:], error[1:], label=f'{name} Error')

ax2.set_yscale('log')
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Absolute Error (C)')
ax2.set_title('Absolute Error Comparison')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.show()

