"""
Diagnostic script to understand heating speed and curing rate.
"""
import numpy as np
import matplotlib.pyplot as plt
from MODEL_PARAMS import *

# === Time scales ===
print("=== Time scales ===")
print(f"t_diff = {tdiff:.1f} s = {tdiff/60:.1f} min")
print(f"Dimensionless t*=1 => physical t = {tdiff:.0f} s")
print()

gamma = h*L/k
print(f"gamma = hL/k = {gamma:.1f}  (Biot number)")
print()

# === Reaction rates ===
print("=== Reaction rate at T=Tm (u=0), alpha=0.1 ===")
rate_phys = A * np.exp(-E/(R*Tm)) * 0.1**m * (1-0.1)**n
rate_dimless = rate_phys * tdiff
print(f"  Physical rate:      da/dt  = {rate_phys:.6e} s^-1")
print(f"  Dimensionless rate: da/dt* = {rate_dimless:.4f}")
print(f"  Over t*=1 (1 diffusion time): delta_alpha ~ {rate_dimless:.4f}")
print()

# === Rough curing time ===
peak_rate_phys = A * np.exp(-E/(R*Tm)) * 0.4**m * (1-0.4)**n
time_to_cure_phys = 0.8 / peak_rate_phys
time_to_cure_dimless = time_to_cure_phys / tdiff
print(f"=== Rough curing time estimate (at T=Tm, alpha~0.4) ===")
print(f"  Peak physical rate: {peak_rate_phys:.6e} s^-1")
print(f"  Rough time 0.1->0.9: {time_to_cure_phys:.0f} s = {time_to_cure_phys/3600:.1f} hours")
print(f"  In dimensionless units: t* = {time_to_cure_dimless:.1f}")
print(f"  Current sim runs to t*=1, need t*~{time_to_cure_dimless:.0f}")
print()

# === Integrate alpha ODE at T=Tm to see S-curve ===
print("=== Integrating alpha ODE at T=Tm (best case) ===")
T_end_star = 15.0
N_ode = 50000
dt_ode = T_end_star / N_ode
alpha_ode = np.zeros(N_ode+1)
alpha_ode[0] = 0.1
t_ode = np.linspace(0, T_end_star, N_ode+1)

for i in range(N_ode):
    rate = tdiff * A * np.exp(-E/(R*Tm)) * alpha_ode[i]**m * (1-alpha_ode[i])**n
    alpha_ode[i+1] = min(alpha_ode[i] + dt_ode * rate, 1.0)

idx_90 = np.argmax(alpha_ode >= 0.9)
if idx_90 > 0:
    t_90 = t_ode[idx_90]
    print(f"  Alpha reaches 0.9 at t* = {t_90:.2f} ({t_90*tdiff:.0f} s = {t_90*tdiff/3600:.1f} h)")
else:
    print(f"  Alpha does NOT reach 0.9 within t*={T_end_star}")
print(f"  Final alpha at t*={T_end_star}: {alpha_ode[-1]:.4f}")
print()

# === Run Heat equation and check center temperature ===
print("=== Running Heat equation (M=50, N=1000) ===")
from Heat import Heat
pde = Heat()
x, y, t_heat, U_heat = pde.solve(M=50, N=1000)
mid = len(x) // 2
T_centre = (Tm - T0) * U_heat[mid, mid, :] + Tm

# Find when center reaches 99% of Tm
idx_99 = np.argmax(T_centre >= 0.99 * Tm)
if idx_99 > 0:
    t99 = t_heat[idx_99] * tdiff
    print(f"  Centre reaches 99% of Tm at t = {t99:.0f} s ({t99/60:.1f} min), t* = {t_heat[idx_99]:.4f}")

idx_Ts = np.argmax(T_centre >= T_s)
if idx_Ts > 0:
    tTs = t_heat[idx_Ts] * tdiff
    print(f"  Centre reaches T_s={T_s} K at t = {tTs:.0f} s ({tTs/60:.1f} min), t* = {t_heat[idx_Ts]:.4f}")

print(f"  Final centre T = {T_centre[-1]:.2f} K")

# === Plots ===
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Plot 1: Centre temperature over time (Heat only)
axes[0].plot(t_heat * tdiff / 60, T_centre, 'b-')
axes[0].axhline(Tm, color='r', ls='--', label=f"Tm = {Tm} K")
axes[0].axhline(T_s, color='orange', ls='--', label=f"Ts = {T_s} K")
axes[0].set_xlabel("Time [min]")
axes[0].set_ylabel("Centre Temperature [K]")
axes[0].set_title("Heat eq: Centre T vs time")
axes[0].legend()

# Plot 2: Alpha S-curve (ODE at T=Tm)
axes[1].plot(t_ode * tdiff / 3600, alpha_ode, 'g-')
axes[1].axhline(0.9, color='r', ls='--', label="90% cure")
axes[1].set_xlabel("Time [hours]")
axes[1].set_ylabel("alpha")
axes[1].set_title("Alpha ODE at T=Tm (best case S-curve)")
axes[1].legend()

# Plot 3: Alpha S-curve zoomed to current sim time
mask = t_ode <= 1.0
axes[2].plot(t_ode[mask] * tdiff / 60, alpha_ode[mask], 'g-', label="ODE at T=Tm")
axes[2].set_xlabel("Time [min]")
axes[2].set_ylabel("alpha")
axes[2].set_title("Alpha during 1 diffusion time (current sim)")
axes[2].legend()

fig.tight_layout()
plt.savefig("diagnose.png", dpi=150)
print("\nSaved diagnose.png")
plt.show()
