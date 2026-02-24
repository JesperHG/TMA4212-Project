"""
Verify the Vulcanization fix: run with T_end=15 and check for S-curve in alpha.
"""
import numpy as np
import matplotlib.pyplot as plt
from MODEL_PARAMS import Tm, T0, T_s, tdiff
from Vulcanization import Vulcanization

# Run with T_end=15 (about 17 hours), N=3000 for reasonable resolution
print("Running Vulcanization with T_end=15, M=30, N=3000 ...")
pde = Vulcanization()
x, y, t, U, alpha = pde.solve(M=30, N=3000, T_end=15)
mid = len(x) // 2

T_centre = (Tm - T0) * U[mid, mid, :] + Tm
alpha_centre = alpha[mid, mid, :]
T_max_over_time = (Tm - T0) * U.max(axis=(0,1)) + Tm

# Find time to 90% cure at centre
idx_90 = np.argmax(alpha_centre >= 0.9)
if idx_90 > 0:
    t_90 = t[idx_90]
    print(f"Centre alpha reaches 0.9 at t* = {t_90:.2f} ({t_90*tdiff/3600:.1f} hours)")
else:
    print(f"Centre alpha does NOT reach 0.9 (final: {alpha_centre[-1]:.4f})")

print(f"Max temperature anywhere: {T_max_over_time.max():.2f} K (mold = {Tm} K)")
print(f"Final centre alpha: {alpha_centre[-1]:.4f}")
print(f"Final centre T: {T_centre[-1]:.2f} K")

# Plots
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

# Centre temperature
axes[0,0].plot(t * tdiff / 3600, T_centre, 'b-')
axes[0,0].axhline(Tm, color='r', ls='--', label=f"Tm = {Tm} K")
axes[0,0].axhline(T_s, color='orange', ls='--', label=f"Ts = {T_s} K")
axes[0,0].set_xlabel("Time [hours]")
axes[0,0].set_ylabel("Temperature [K]")
axes[0,0].set_title("Centre Temperature")
axes[0,0].legend()

# Centre alpha (the S-curve)
axes[0,1].plot(t * tdiff / 3600, alpha_centre, 'g-', linewidth=2)
axes[0,1].axhline(0.9, color='r', ls='--', label="90% cure")
axes[0,1].set_xlabel("Time [hours]")
axes[0,1].set_ylabel("alpha")
axes[0,1].set_title("Centre Degree of Cure (S-curve)")
axes[0,1].legend()

# Max temperature over time
axes[1,0].plot(t * tdiff / 3600, T_max_over_time, 'r-')
axes[1,0].axhline(Tm, color='gray', ls='--', label=f"Tm = {Tm} K")
axes[1,0].set_xlabel("Time [hours]")
axes[1,0].set_ylabel("Max T [K]")
axes[1,0].set_title("Max Temperature (overheating check)")
axes[1,0].legend()

# Alpha snapshots across space at different times
snapshot_indices = [0, len(t)//10, len(t)//4, len(t)//2, 3*len(t)//4, -1]
for idx in snapshot_indices:
    axes[1,1].plot(x, alpha[:, mid, idx], label=f"t={t[idx]*tdiff/3600:.1f}h")
axes[1,1].set_xlabel("x")
axes[1,1].set_ylabel("alpha")
axes[1,1].set_title("Alpha profile (y=0.5) at different times")
axes[1,1].legend(fontsize=7)

fig.tight_layout()
plt.savefig("verify_vulcanization.png", dpi=150)
print("Saved verify_vulcanization.png")
plt.show()
