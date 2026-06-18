import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from src.mhd import solve_alfven

OUT = os.path.dirname(os.path.abspath(__file__))

sol = solve_alfven(save_every=8)
x      = sol['x']
By     = sol['By']
Bz     = sol['Bz']
eps    = sol['eps']
phi0   = 2.0 * np.pi * x
L      = sol['L']
snapshots = sol['snapshots']
snap_times = sol['snap_times']

print(f"Done: {sol['step']} steps, {len(snapshots)} frames, t={sol['t']:.4f}")

fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
fig.suptitle('1D CP Alfvén Wave — Python (RK2)', fontsize=13)
ax1, ax2 = axes
ax1.plot(x, eps * np.sin(phi0), 'k--', alpha=0.3, label='t=0')
ax1.plot(x, By, 'b-', lw=1.5, label=f't={snap_times[-1]:.2f}')
ax1.set_ylabel('By'); ax1.legend(); ax1.grid(True, alpha=0.3)
ax2.plot(x, eps * np.cos(phi0), 'k--', alpha=0.3, label='t=0')
ax2.plot(x, Bz, 'r-', lw=1.5, label=f't={snap_times[-1]:.2f}')
ax2.set_ylabel('Bz'); ax2.set_xlabel('x')
ax2.legend(); ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'alfven_final.png'), dpi=150)

fig2, axes2 = plt.subplots(2, 1, figsize=(10, 6))
fig2.suptitle('Onda de Alfvén CP - 1D MHD Ideal (RK2)', fontsize=13)
ax1, ax2 = axes2
ax1.set_xlim(0, L); ax1.set_ylim(-eps*1.5, eps*1.5)
ax1.set_ylabel('By'); ax1.grid(True, alpha=0.3)
ax1.plot(x, eps*np.sin(phi0), 'k--', alpha=0.3, label='t=0')
ax2.set_xlim(0, L); ax2.set_ylim(-eps*1.5, eps*1.5)
ax2.set_ylabel('Bz'); ax2.set_xlabel('x'); ax2.grid(True, alpha=0.3)
ax2.plot(x, eps*np.cos(phi0), 'k--', alpha=0.3, label='t=0')
line_By, = ax1.plot([], [], 'b-', lw=2, label='By(t)')
line_Bz, = ax2.plot([], [], 'r-', lw=2, label='Bz(t)')
time_text = ax1.text(0.02, 0.88, '', transform=ax1.transAxes, fontsize=10,
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
ax1.legend(); ax2.legend()

def init():
    line_By.set_data([], []); line_Bz.set_data([], [])
    time_text.set_text(''); return line_By, line_Bz, time_text

def animate(i):
    By_i, Bz_i, t_i = snapshots[i]
    line_By.set_data(x, By_i); line_Bz.set_data(x, Bz_i)
    time_text.set_text(f't = {t_i:.2f}')
    return line_By, line_Bz, time_text

ani = animation.FuncAnimation(fig2, animate, init_func=init,
                               frames=len(snapshots), interval=30, blit=True)
ani.save(os.path.join(OUT, 'alfven_animation.gif'), writer='pillow', fps=30, dpi=100)
plt.close('all')

np.savez(os.path.join(OUT, 'alfven_solution.npz'),
         x=x, By=By, Bz=Bz, vy=sol['vy'], vz=sol['vz'],
         rho=np.ones(sol['N'])*sol['rho0'], P=np.ones(sol['N'])*sol['P0'], t=sol['t'])
print("Saved: alfven_final.png, alfven_animation.gif, alfven_solution.npz")
