
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.linalg import schur, expm
from scipy.signal import savgol_filter


# User-configurable params

DATA_DIR = "/home/zilani/1.Md_Zilani/Python/Project 2/Data/"   
DECAY_CSV = "/home/zilani/1.Md_Zilani/Python/Project 2/Data/decay_data.csv"
OUTDIR = "/home/zilani/1.Md_Zilani/Python/Project 2/Result"
os.makedirs(OUTDIR, exist_ok=True)

N = 100                # Monte-Carlo samples (analytic method is very fast)
nt = 500              # time resolution for plots (increase for smoother curves)
use_smoothing = True   # apply Savitzky-Golay smoothing to percentile curves
savgol_window = 51     # must be odd and <= nt
savgol_poly = 3

# Physics / test defaults
T = 1e9                # temperature (unused for correct MACS -> sigma conversion, left for bookkeeping)
n_n = 1e20             # neutron density (change if you want slower captures for testing)

# Auto-time parameters
auto_time = True       # use automatic time window determined from rates
time_scale_factor = 100.0  # simulate to (time_scale_factor * tau_min), where tau_min = 1/max(rate)
min_t_end = 1e-12      # avoid zero
max_t_end = 1e7

rng = np.random.default_rng()


def load_mc_summary(isotope, data_dir=DATA_DIR, kT_keV=30.0):
    path = f"{data_dir}/mc_macs_summary_{isotope}_ng_TALYS.csv"
    df = pd.read_csv(path)
    df = df[df["kT_keV"] == kT_keV]
    if df.empty:
        raise RuntimeError(f"No MACS summary found for {isotope} at kT={kT_keV} keV in {path}")
    mean = float(df["MACS_mean(b)"].values[0])
    std  = float(df["MACS_std(b)"].values[0])
    return mean, std

def macs_to_sigma_v(macs_b):
    """
    Convert MACS (in barns) to cross-section <σ> in m^2.
    IMPORTANT: MACS is already velocity-averaged; do NOT multiply by thermal velocity.
    """
    return macs_b * 1e-28   # barns -> m^2

def build_A_matrix(lam_cap, lam_beta):
    """
    Build 5x5 coefficient matrix A for the chain:
    56 -> 57 -> 58 -> 59 -> 60
    Y' = A Y
    lam_cap: list-like [l0, l1, l2, l3]
    lam_beta: list-like [0,0,0,b3,b4]
    """
    l0, l1, l2, l3 = lam_cap
    b3 = lam_beta[3]
    b4 = lam_beta[4]

    A = np.zeros((5,5), dtype=float)
    A[0,0] = -l0

    A[1,0] = l0
    A[1,1] = -l1

    A[2,1] = l1
    A[2,2] = -l2

    A[3,2] = l2
    A[3,3] = -(l3 + b3)

    A[4,3] = l3
    A[4,4] = -b4

    return A

def solve_chain_schur(A, Y0, t_eval):
    """
    Solve Y(t) = exp(A t) Y0 using real Schur decomposition:
    A = Q T Q^T  -> exp(A t) = Q exp(T t) Q^T
    Returns Ys of shape (5, len(t_eval)).
    """
    Tmat, Q = schur(A, output='real')
    z0 = Q.T @ Y0
    nt_local = len(t_eval)
    Ys = np.zeros((A.shape[0], nt_local), dtype=float)
    for j, t in enumerate(t_eval):
        expTt = expm(Tmat * t)
        z = expTt @ z0
        Ys[:, j] = (Q @ z).real
    return Ys

# ------------------------
# MAIN
# ------------------------
if __name__ == "__main__":
    # Load MACS summaries
    mean_56, std_56 = load_mc_summary("56Fe")
    mean_57, std_57 = load_mc_summary("57Fe")
    mean_58, std_58 = load_mc_summary("58Fe")
    mean_59, std_59 = load_mc_summary("59Fe")

    # Load decays
    decay = pd.read_csv(DECAY_CSV)
    decay_dict = dict(zip(decay['Isotope'], decay['HalfLife_s']))
    lam_beta = [
        0.0, 0.0, 0.0,
        np.log(2)/decay_dict["59Fe"],
        np.log(2)/decay_dict["60Fe"]
    ]

    # Pre-allocate results
    all_solutions = np.empty((N, 5, nt), dtype=float)
    Y0 = np.array([1.0, 0.0, 0.0, 0.0, 0.0])

    # We'll first do a single-sample test to compute time window if auto_time=True
    # Sample means for this estimate (deterministic)
    macs_means = [mean_56, mean_57, mean_58, mean_59]
    sigma_v_means = [macs_to_sigma_v(m) for m in macs_means]
    lam_cap_means = np.array([n_n * sv for sv in sigma_v_means], dtype=float)

    # include beta rates in timescale estimate
    combined_rates = np.concatenate([lam_cap_means, np.array([lam_beta[3], lam_beta[4]])])
    max_rate = np.max(combined_rates)
    if max_rate <= 0:
        tau_min = 1.0
    else:
        tau_min = 1.0 / max_rate

    if auto_time:
        t_end = float(np.clip(time_scale_factor * tau_min, min_t_end, max_t_end))
    else:
        t_end = max_t_end  # fallback

    print(f"n_n = {n_n:.3e}, max_rate (est) = {max_rate:.3e} s^-1, tau_min = {tau_min:.3e} s")
    print(f"Auto time span selected: t = 0 .. {t_end:.3e} s (nt = {nt})")

    t_eval = np.linspace(0.0, t_end, nt)

    # Monte-Carlo loop
    for i in range(N):
        # sample MACS
        m56 = float(rng.normal(mean_56, std_56))
        m57 = float(rng.normal(mean_57, std_57))
        m58 = float(rng.normal(mean_58, std_58))
        m59 = float(rng.normal(mean_59, std_59))

        macs_list = [m56, m57, m58, m59]
        sigma_v_list = [macs_to_sigma_v(x) for x in macs_list]   # m^2
        lam_cap = [n_n * sv for sv in sigma_v_list]               # s^-1

        A = build_A_matrix(lam_cap, lam_beta)
        Ys = solve_chain_schur(A, Y0, t_eval)    # shape (5, nt)
        all_solutions[i] = Ys

    # compute percentile bands along the MC axis (axis=0 -> (N,5,nt) -> percentile -> (5,nt))
    p16 = np.percentile(all_solutions, 16, axis=0)
    p50 = np.percentile(all_solutions, 50, axis=0)
    p84 = np.percentile(all_solutions, 84, axis=0)

    # optional smoothing (applied per isotope)
    if use_smoothing:
        # ensure window is valid
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
                print("Smoothing failed:", e)
                # fallback to unsmoothed

        plt.figure(figsize=(8,4.5))
        plt.plot(t_eval, y50, label='Median', lw=2)
        plt.fill_between(t_eval, y16, y84, alpha=0.3, label='16th-84th')
        plt.title(f"{labels[idx]} abundance uncertainty")
        plt.xlabel("Time (s)")
        plt.ylabel("Y")
        plt.grid(True, ls='--', alpha=0.4)
        plt.legend()
        plt.tight_layout()

        outpath = os.path.join(OUTDIR, f"uncertainty_{labels[idx]}.png")
        plt.savefig(outpath, dpi=200)
        print("Saved:", outpath)
        plt.close()

    print("Done.")

