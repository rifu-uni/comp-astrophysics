import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from src.mhd import solve_alfven

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sol = solve_alfven(save_every=8)
x, By, Bz = sol['x'], sol['By'], sol['Bz']
eps, L, phi0 = sol['eps'], sol['L'], 2.0 * np.pi * sol['x']
snapshots, snap_times = sol['snapshots'], sol['snap_times']

print(f"Done: {sol['step']} steps, {len(snapshots)} frames, t={sol['t']:.4f}")

fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
fig.suptitle('1D CP Alfvén Wave — Python (RK2)', fontsize=13)
for ax, var, c, label in [(axes[0], By, 'b', 'By'), (axes[1], Bz, 'r', 'Bz')]:
    ax.plot(x, eps * np.sin(phi0 if label == 'By' else phi0), 'k--', alpha=0.3, label='t=0')
    ax.plot(x, var, f'{c}-', lw=1.5, label=f'{label} t={snap_times[-1]:.2f}')
    ax.set_ylabel(label); ax.legend(); ax.grid(True, alpha=0.3)
axes[1].set_xlabel('x')
plt.tight_layout(); plt.savefig(os.path.join(OUT, 'alfven_final.png'), dpi=150)

fig2, axes2 = plt.subplots(2, 1, figsize=(10, 6))
fig2.suptitle('CP Alfvén Wave — Time Evolution (RK2)', fontsize=13)
for ax, c in [(axes2[0], 'b'), (axes2[1], 'r')]:
    ax.set_xlim(0, L); ax.set_ylim(-eps*1.5, eps*1.5); ax.grid(True, alpha=0.3)
    ax.plot(x, eps*np.sin(phi0), 'k--', alpha=0.3, label='t=0')
axes2[0].set_ylabel('By'); axes2[1].set_ylabel('Bz'); axes2[1].set_xlabel('x')
line_By, = axes2[0].plot([], [], 'b-', lw=2, label='By(t)')
line_Bz, = axes2[1].plot([], [], 'r-', lw=2, label='Bz(t)')
time_text = axes2[0].text(0.02, 0.88, '', transform=axes2[0].transAxes, fontsize=10,
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
for ax in axes2: ax.legend()

def animate(i):
    By_i, Bz_i, t_i = snapshots[i]
    line_By.set_data(x, By_i); line_Bz.set_data(x, Bz_i)
    time_text.set_text(f't = {t_i:.2f}')
    return line_By, line_Bz, time_text

ani = animation.FuncAnimation(fig2, animate, frames=len(snapshots), interval=30, blit=True)
ani.save(os.path.join(OUT, 'alfven_animation.gif'), writer='pillow', fps=30, dpi=100)
plt.close('all')
np.savez(os.path.join(OUT, 'alfven_solution.npz'),
         x=x, By=By, Bz=Bz, vy=sol['vy'], vz=sol['vz'],
         rho=np.ones(sol['N'])*sol['rho0'], P=np.ones(sol['N'])*sol['P0'], t=sol['t'])
print("Saved: alfven_final.png, alfven_animation.gif, alfven_solution.npz")
