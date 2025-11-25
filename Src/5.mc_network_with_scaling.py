import seaborn as sns
import numpy as np
import pandas as pd
import os
import math
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.linalg import schur, expm
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

DATA_DIR = "/home/zilani/1.Md_Zilani/Python/Project 2/Data/"
OUT_DIR = "/home/zilani/1.Md_Zilani/Python/Project 2/Result/"
os.makedirs(OUT_DIR, exist_ok=True)

TALYS_FILES = [
    f"{DATA_DIR}/⁵⁶Fe(n,g).csv",
    f"{DATA_DIR}/⁵⁷Fe(n,g).csv",
    f"{DATA_DIR}/⁵⁸Fe(n,g).csv",
    f"{DATA_DIR}/⁵⁹Fe(n,g).csv",
]
MC_MACS_ISOTOPES = ["⁵⁶Fe", "⁵⁷Fe", "⁵⁸Fe", "⁵⁹Fe"]

DECAY_CSV = f"{DATA_DIR}/decay_data.csv"

T = 1e9             
n_n = 1e20         
N = 200            
nt = 800           
time_scale_factor = 100.0   
min_t_end = 1e-12
max_t_end = 1e7

use_smoothing = True
savgol_window = 51
savgol_poly = 3

rng = np.random.default_rng(12345)

EV_TO_J = 1.602176634e-19
MEV_TO_J = EV_TO_J * 1e6
BARN_TO_M2 = 1e-28

def load_sigma_csv(path):
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
  
    E = df['Energy(MeV)'].values * MEV_TO_J
    s = df['CrossSection(b)'].values * BARN_TO_M2

    E, idx = np.unique(E, return_index=True)
    s = s[idx]
    order = np.argsort(E)
    return E[order], s[order]

def sigma_v_from_sigmaE(E_J, sigma_m2, T):
  
    k_B = 1.380649e-23
    kT = k_B * T
    mu = 1.675e-27  
 
    pref = 2.0 / math.sqrt(math.pi) * (1.0 / (kT)**1.5) * math.sqrt(1.0 / (2.0 * mu))
    sigma_interp = interp1d(E_J, sigma_m2, kind="linear", fill_value=0.0, bounds_error=False)

    def integrand(E):
        return sigma_interp(E) * np.sqrt(E) * np.exp(-E / kT)

    E_min = max(E_J.min(), 1e-40)
    E_max = E_J.max()
    val, err = quad(integrand, E_min, E_max, limit=400)
    return pref * val

def build_A_matrix(lam_cap, lam_beta):
    l0, l1, l2, l3 = lam_cap
    b3 = lam_beta[3]
    b4 = lam_beta[4]
    A = np.zeros((5,5), dtype=float)
    A[0,0] = -l0
    A[1,0] = l0; A[1,1] = -l1
    A[2,1] = l1; A[2,2] = -l2
    A[3,2] = l2; A[3,3] = -(l3 + b3)
    A[4,3] = l3; A[4,4] = -b4
    return A

def solve_chain_schur(A, Y0, t_eval):
    Tmat, Q = schur(A, output='real')
    z0 = Q.T @ Y0
    Ys = np.zeros((A.shape[0], len(t_eval)), dtype=float)
    for j, t in enumerate(t_eval):
        expTt = expm(Tmat * t)
        z = expTt @ z0
        Ys[:, j] = (Q @ z).real
    return Ys

if __name__ == "__main__":
    print("1) Computing deterministic <σv> from TALYS sigma(E)...")
    sigma_v_list = []
    for path in TALYS_FILES:
        E_J, s_m2 = load_sigma_csv(path)
        sv = sigma_v_from_sigmaE(E_J, s_m2, T)
        sigma_v_list.append(sv)
        print(f"  {os.path.basename(path)}  ->  <σv> = {sv:.3e} m^3/s")
    sigma_v_list = np.array(sigma_v_list, dtype=float)
    lam_cap_det = n_n * sigma_v_list
    for i, lam in enumerate(lam_cap_det):
        print(f"  Deterministic lam_cap[{i}] = {lam:.3e} s^-1")

  
    decay = pd.read_csv(DECAY_CSV)
    decay_dict = dict(zip(decay['Isotope'], decay['HalfLife_s']))
    lam_beta = [
        0.0, 0.0, 0.0,
        np.log(2)/decay_dict["59Fe"],
        np.log(2)/decay_dict["60Fe"]
    ]

 
    macs_means = []
    macs_stds = []
    for iso in MC_MACS_ISOTOPES:
        path = f"{DATA_DIR}/mc_macs_summary_{iso}_ng_TALYS.csv"
        df = pd.read_csv(path)
        df = df[df["kT_keV"] == 30.0]
        if df.empty:
            raise RuntimeError(f"Missing MACS summary for {iso} at 30 keV: {path}")
        macs_means.append(float(df["MACS_mean(b)"].values[0]))
        macs_stds.append(float(df["MACS_std(b)"].values[0]))
    macs_means = np.array(macs_means, dtype=float)
    macs_stds = np.array(macs_stds, dtype=float)
    print("\nMACS means (barn):", macs_means)
    print("MACS stds  (barn):", macs_stds)


    combined = np.concatenate([lam_cap_det, np.array([lam_beta[3], lam_beta[4]])])
    max_rate = np.max(combined)
    tau_min = 1.0 / max_rate if max_rate > 0 else 1.0
    t_end = float(np.clip(time_scale_factor * tau_min, min_t_end, max_t_end))
    print(f"\nEstimated max_rate = {max_rate:.3e} s^-1, tau_min = {tau_min:.3e} s")
    print(f"Using time span: 0 .. {t_end:.3e} s, nt = {nt}")
    t_eval = np.linspace(0.0, t_end, nt)

    
    all_solutions = np.empty((N, 5, nt), dtype=float)
    Y0 = np.array([1.0, 0.0, 0.0, 0.0, 0.0])

    lam_cap_mc_means = np.zeros(4, dtype=float)  
    for i in range(N):
        
        samples = rng.normal(macs_means, macs_stds)
     
        samples = np.clip(samples, 1e-12, None)

       
        f = samples / macs_means

      
        sigma_v_mc = sigma_v_list * f       
        lam_cap_mc = n_n * sigma_v_mc       

        lam_cap_mc_means += lam_cap_mc / N   

   
        A = build_A_matrix(lam_cap_mc.tolist(), lam_beta)
        Ys = solve_chain_schur(A, Y0, t_eval)
        all_solutions[i] = Ys

    print("\nMean MC lam_cap (s^-1) across samples:")
    for j, m in enumerate(lam_cap_mc_means):
        print(f"  lam_cap_mc_mean[{j}] = {m:.3e} (deterministic lam_cap[{j}] = {lam_cap_det[j]:.3e})")


    p16 = np.percentile(all_solutions, 16, axis=0)
    p50 = np.percentile(all_solutions, 50, axis=0)
    p84 = np.percentile(all_solutions, 84, axis=0)


    if use_smoothing:
        if savgol_window >= nt:
            savgol_window = nt - 1 if (nt - 1) % 2 == 1 else nt - 2
        if savgol_window < 5:
            use_smoothing = False

    labels = ["56Fe","57Fe","58Fe","59Fe","60Fe"]
    for idx in range(5):
        y50 = p50[idx].copy()
        y16 = p16[idx].copy()
        y84 = p84[idx].copy()
        if use_smoothing:
            try:
                y50 = savgol_filter(y50, savgol_window, savgol_poly)
                y16 = savgol_filter(y16, savgol_window, savgol_poly)
                y84 = savgol_filter(y84, savgol_window, savgol_poly)
            except Exception as e:
                print("Smoothing error:", e)
           
        plt.figure(figsize=(9, 5))

        base_color = sns.color_palette("viridis", 8)[5]
        band_color = sns.color_palette("viridis", 8)[3]
      
        plt.plot(
            t_eval, y50,
            color=base_color,
            linewidth=2.5,
            label="Median"
        )
   
        plt.fill_between(
            t_eval, y16, y84,
            color=band_color,
            alpha=0.25,
            label="16th–84th Percentile"
        )

        # Title & labels
        plt.title(
            f"{labels[idx]} Monte-Carlo Abundance Uncertainty ",
            fontsize=15,
            fontweight="bold",
            pad=10
        )
        plt.xlabel("Time (s)", fontsize=12)
        plt.ylabel("Abundance Y", fontsize=12)

        # Light scientific grid
        plt.grid(alpha=0.30, linestyle="--")

        # Legend styling
        plt.legend(
            fontsize=10,
            frameon=True,
            fancybox=True
        )

        plt.tight_layout()

        # High-resolution save
        out = os.path.join(OUT_DIR, f"mc_scaled_uncertainty_{labels[idx]}.png")
        plt.savefig(out, dpi=350)
        print("Saved:", out)

     
        plt.show()

    print("Done.")
