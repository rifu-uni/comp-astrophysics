import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.mhd import alfven_analytic, convergence_test

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
data = np.load(os.path.join(OUT, 'alfven_solution.npz'))
x, By_py, Bz_py, t = data['x'], data['By'], data['Bz'], float(data['t'])
eps = 0.1

By_a, Bz_a = alfven_analytic(x, t, eps)
for name, py, an in [('By', By_py, By_a), ('Bz', Bz_py, Bz_a)]:
    err = np.abs(py - an)
    print(f"  L1({name})={np.mean(err):.2e}  L2({name})={np.sqrt(np.mean(err**2)):.2e}")

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle(f'Python Solver vs Analytic (t={t:.2f})', fontsize=14)
for i, (py, an, name) in enumerate([(By_py, By_a, 'By'), (Bz_py, Bz_a, 'Bz')]):
    axes[0,i].plot(x, py, label=f'Python {name}')
    axes[0,i].plot(x, an, '--', label=f'Analytic {name}')
    axes[0,i].set_ylabel(name); axes[0,i].legend(); axes[0,i].grid(True, alpha=0.3)
    err = np.abs(py - an)
    axes[1,i].plot(x, err)
    axes[1,i].set_ylabel(f'|{name} error|'); axes[1,i].set_xlabel('x')
    axes[1,i].grid(True, alpha=0.3)
plt.tight_layout(); plt.savefig(os.path.join(OUT, 'compare_python.png'), dpi=150)

resolutions, L1_vals, L2_vals = convergence_test()
fig2, ax2 = plt.subplots(figsize=(8, 5))
ax2.loglog(resolutions, L1_vals, 'bo-', label='L1')
ax2.loglog(resolutions, L2_vals, 'rs-', label='L2')
ax2.loglog(resolutions, [1/r**2 for r in resolutions], 'k--', alpha=0.5, label='O(N⁻²)')
ax2.set_xlabel('N'); ax2.set_ylabel('Error'); ax2.set_title('Convergence (RK2)')
ax2.legend(); ax2.grid(True, alpha=0.3)
plt.tight_layout(); plt.savefig(os.path.join(OUT, 'convergence.png'), dpi=150)

print("\nConvergence:")
for N, l1, l2 in zip(resolutions, L1_vals, L2_vals):
    print(f"  N={N:4d}: L1={l1:.2e}  L2={l2:.2e}")
print("Saved: compare_python.png, convergence.png")
