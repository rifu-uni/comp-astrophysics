import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.vtk_reader import read_pluto_vtk_coarse
from src.mhd import alfven_analytic

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUILD = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'build')
vtk_file = os.path.join(BUILD, 'data.0001.vtk')

if not os.path.isfile(vtk_file):
    print(f"VTK file not found: {vtk_file}"); sys.exit(1)

fields = read_pluto_vtk_coarse(vtk_file)
t_pluto = 1.0
x = np.linspace(1/256, 1 - 1/256, 128)
eps = 0.1
By_analytic, Bz_analytic = alfven_analytic(x, t_pluto, eps)

error_By = np.abs(fields['Bx2'] - By_analytic)
error_Bz = np.abs(fields['Bx3'] - Bz_analytic)
L1_By = np.mean(error_By); L1_Bz = np.mean(error_Bz)
L2_By = np.sqrt(np.mean(error_By**2)); L2_Bz = np.sqrt(np.mean(error_Bz**2))

print(f"PLUTO vs Analytic at t={t_pluto}:")
print(f"  L1( By ) = {L1_By:.2e}")
print(f"  L1( Bz ) = {L1_Bz:.2e}")
print(f"  L2( By ) = {L2_By:.2e}")
print(f"  L2( Bz ) = {L2_Bz:.2e}")

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle('PLUTO CP_Alfven — Comparison with Analytic Solution', fontsize=14)
axes[0,0].plot(x, fields['Bx2'], label='PLUTO By')
axes[0,0].plot(x, By_analytic, '--', label='Analytic By')
axes[0,0].set_ylabel('By'); axes[0,0].legend(); axes[0,0].grid(True, alpha=0.3)
axes[0,1].plot(x, fields['Bx3'], label='PLUTO Bz')
axes[0,1].plot(x, Bz_analytic, '--', label='Analytic Bz')
axes[0,1].set_ylabel('Bz'); axes[0,1].legend(); axes[0,1].grid(True, alpha=0.3)
axes[1,0].plot(x, error_By)
axes[1,0].set_ylabel('|By error|'); axes[1,0].set_xlabel('x')
axes[1,0].set_title(f'L1={L1_By:.2e}, L2={L2_By:.2e}')
axes[1,0].grid(True, alpha=0.3)
axes[1,1].plot(x, error_Bz)
axes[1,1].set_ylabel('|Bz error|'); axes[1,1].set_xlabel('x')
axes[1,1].set_title(f'L1={L1_Bz:.2e}, L2={L2_Bz:.2e}')
axes[1,1].grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'alfven_error.png'), dpi=150)
print("Saved: alfven_error.png")
