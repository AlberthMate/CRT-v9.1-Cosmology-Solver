import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ---- Cosmological parameters (Planck 2018) --------------------------------
H0       = 67.97e3 / 3.085677581e22   # [s^-1]
G        = 6.67430e-11                 # [m^3 kg^-1 s^-2]
rho_crit = 3.0 * H0**2 / (8.0 * np.pi * G)
Omega_m0, Omega_r0, xi = 0.315, 9.2e-5, 1.13

# ---- CRT v9.1 parameters (calibrated) ------------------------------------
a_star     = 0.5    # freeze-in transition scale factor
n_freeze   = 1.0    # freeze-in steepness
beta       = 1.5    # holographic fall-off exponent
r_IR       = 0.9    # L_IR = r_IR * c / H0
coeff      = xi / r_IR**2   # Omega_Lambda = coeff / Lambda^2
kappa_eff  = 0.30   # dimensionless calibrated mobility (C_hol~3.2 absorbed)
Lambda_ini = 1.38   # L_ini/L_IR, normalized at present epoch (see Sec. V.A)

# ---- Dimensionless derivative of Sigma_tot(Lambda) -----------------------
def g(Lam):
    """
    g(Lambda) = dSigma/dLambda for the effective entropy ansatz:
    Sigma(Lambda) = Lambda**2 / (1 + Lambda**beta).
    Note: C_hol is absorbed into kappa_eff during calibration.
    """
    u = Lam**beta
    d = 1.0 + u
    return (2.0 * Lam * d - beta * Lam * u) / d**2

# ---- Dimensionless ODE system --------------------------------------------
def system(tau, y):
    a, h, Lam = y
    if a <= 0.0 or h <= 0.0 or Lam <= 0.0:
        return [0.0, 0.0, 0.0]

    Om = Omega_m0 / a**3
    Or = Omega_r0 / a**4
    OL = coeff / Lam**2

    kappa_a = kappa_eff * np.exp(-(a / a_star)**n_freeze)
    dLam    = kappa_a * Lam**2 * g(Lam)    # entropy-ascending flow

    # Dark-energy equation of state (kinematic definition)
    w_L   = np.clip(-1.0 + (2.0/3.0) * (dLam / (Lam * h)), -3.0, 0.5)
    Otot  = Om + Or + OL
    w_eff = (Or / 3.0 + OL * w_L) / Otot if Otot > 1e-30 else 0.0

    return [a * h,  -1.5 * h**2 * (1.0 + w_eff),  dLam]

# ---- Terminal event: stop integration at a = 1 (today) ------------------
def stop_today(tau, y): return y[0] - 1.0
stop_today.terminal  = True
stop_today.direction = +1

# ---- Initial conditions at z_ini = 1000 ----------------------------------
a_ini = 1.0 / 1001.0
h_ini = np.sqrt(Omega_m0 / a_ini**3 + Omega_r0 / a_ini**4)
tau_i = H0 * 1.0e13    # dimensionless time at t ~ 10^13 s
tau_f = H0 * 9.0e17    # upper bound (terminal event stops integration)

sol = solve_ivp(
    system, (tau_i, tau_f),
    [a_ini, h_ini, Lambda_ini],
    method       = 'RK45',
    events       = stop_today,
    rtol         = 1e-9,
    atol         = 1e-11,
    dense_output = True,
    max_step     = (tau_f - tau_i) / 3000,
)

# ---- Post-processing -----------------------------------------------------
a_s, h_s, L_s = sol.y
z_s   = 1.0 / a_s - 1.0
kp    = kappa_eff * np.exp(-(a_s / a_star)**n_freeze)
w_s   = -1.0 + (2.0 / 3.0) * (kp * L_s**2 * g(L_s) / (L_s * h_s))

Lf    = L_s[-1]
w0    = w_s[-1]
rho_L = coeff / Lf**2 * rho_crit

print(f"w_0        = {w0:.4f}")
print(f"rho_Lambda = {rho_L:.4e}  kg/m^3")

# ---- Plot w(z) -----------------------------------------------------------
mask = z_s < 5.0
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(z_s[mask], w_s[mask], 'b-', lw=2, label='CRT v9.1')
ax.axhline(-1.0, color='k', ls='--', lw=1.5, label=r'$\Lambda$CDM')
ax.axhline(-0.76, color='r', ls='-.', lw=1.5,
           label=r'DESI DR2: $w_0 = -0.76\,\pm\,0.13$')
ax.fill_between([0, 5], -0.89, -0.63, color='r', alpha=0.12)
ax.set_xlim(5, 0);  ax.set_ylim(-1.15, -0.60)
ax.set_xlabel('Redshift $z$')
ax.set_ylabel('Equation of state $w(z)$')
ax.legend();  ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('wz_crt_v91.pdf')
\end{lstlisting}
