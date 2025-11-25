import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import quad
import math
import time
import matplotlib.pyplot as plt
import seaborn as sns
import os

EV_TO_J = 1.602176634e-19
MEV_TO_J = EV_TO_J * 1e6
BARN_TO_M2 = 1e-28
sqrt_pi = math.sqrt(math.pi)

def macs_from_arrays(energy_MeV, sigma_b, kT_keV):
    E_J = np.array(energy_MeV, dtype=float) * MEV_TO_J
    sigma_m2 = np.array(sigma_b, dtype=float) * BARN_TO_M2
    kT_J = kT_keV * 1e3 * EV_TO_J

    E_J, idx = np.unique(E_J, return_index=True)
    sigma_m2 = sigma_m2[idx]
    sort_idx = np.argsort(E_J)
    E_J = E_J[sort_idx]
    sigma_m2 = sigma_m2[sort_idx]

    interp = interp1d(E_J, sigma_m2, kind='linear', 
                      bounds_error=False, fill_value=0.0)

    def integrand(E):
        return float(interp(E)) * E * np.exp(-E / kT_J)

    E_min = max(E_J.min(), 1e-40)
    E_max = E_J.max()

    integral, _ = quad(integrand, E_min, E_max, limit=300)
    pref = (2.0 / sqrt_pi) / (kT_J * kT_J)
    macs_m2 = pref * integral

    return macs_m2 / BARN_TO_M2  

def monte_carlo_macs(energy_MeV, sigma_b, kT_keV,
                     n_samples=200, sigma_rel=0.10, seed=None):

    rng = np.random.default_rng(seed)
    samples = np.zeros(n_samples)

    E = np.array(energy_MeV, dtype=float)
    sigma = np.array(sigma_b, dtype=float)

    for i in range(n_samples):
     
        factors = rng.normal(1.0, sigma_rel, size=sigma.shape)
        factors = np.maximum(factors, 1e-6)
        sigma_pert = sigma * factors

        samples[i] = macs_from_arrays(E, sigma_pert, kT_keV)

    return samples


if __name__ == "__main__":

    # ===== File Inputs =====
    talys_file = "/home/zilani/1.Md_Zilani/Python/Project 2/Data/⁵⁹Fe(n,g).csv"
    kT_list = [10, 30, 60]  # keV
    n_samples = 300
    sigma_rel = 0.10    # 10% uncertainty

    outdir = "/home/zilani/1.Md_Zilani/Python/Project 2/Result"
    os.makedirs(outdir, exist_ok=True)

    # read TALYS file
    df = pd.read_csv(talys_file)
    E = df['Energy(MeV)'].values
    sigma = df['CrossSection(b)'].values

    results = []


    for kT in kT_list:

        t0 = time.time()

        samples = monte_carlo_macs(
            E, sigma, kT,
            n_samples=n_samples,
            sigma_rel=sigma_rel,
            seed=42 )

        mean = samples.mean()
        std = samples.std()

        results.append({
            'Reaction': talys_file.split('/')[-1].replace('.csv', ''),
            'kT_keV': kT,
            'MACS_mean(b)': mean,
            'MACS_std(b)': std })

        print(f"kT={kT} keV -> mean={mean:.6e} b, std={std:.3e} b "
              f"(time {time.time()-t0:.1f}s)")

        plt.figure(figsize=(7, 4))
   
        color = sns.color_palette("viridis", 8)[5]
  
        sns.histplot(
            samples,
            bins=30,
            kde=True,
            stat="density",
            color=color,
            edgecolor="black",
            linewidth=0.7,
            alpha=0.80  )

        plt.axvline(mean, color="red", linestyle="--", linewidth=1.8,
                    label=f"Mean = {mean:.3e} b")

        plt.axvspan(mean - std, mean + std,
                    color="red", alpha=0.15,
                    label=f"±1σ = {std:.3e} b")

        plt.title(
            f"MACS Distribution — {talys_file.split('/')[-1].replace('.csv','')} (kT = {kT} keV)",
            fontsize=14, fontweight="bold"
        )
        plt.xlabel("MACS (barn)", fontsize=12)
        plt.ylabel("Probability Density", fontsize=12)

        plt.grid(alpha=0.25, linestyle="--")
        plt.legend(fontsize=10, frameon=True, fancybox=True)

        plt.tight_layout()

        out_png = f"{outdir}/mc_hist_{talys_file.split('/')[-1].replace('.csv','')}_kT{kT}.png"
        plt.savefig(out_png, dpi=350)
        print("Saved:", out_png)
        plt.show()

    out_df = pd.DataFrame(results)
    out_csv = f"{outdir}/mc_macs_summary_{talys_file.split('/')[-1].replace('.csv','')}.csv"
    out_df.to_csv(out_csv, index=False)
    print("Summary saved to", out_csv)




