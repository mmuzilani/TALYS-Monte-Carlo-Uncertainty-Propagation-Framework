
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

EV_TO_J = 1.602176634e-19
MEV_TO_J = EV_TO_J * 1e6
BARN_TO_M2 = 1e-28


talys_file = "/home/zilani/1.Md_Zilani/Python/Project 2/Data/⁵⁷Fe(n,g).csv"
kT_keV = 30.0

df = pd.read_csv(talys_file)
df.columns = [c.strip() for c in df.columns]

E_MeV = df['Energy(MeV)'].values
sigma_b = df['CrossSection(b)'].values


mask = (E_MeV > 0)
E_MeV = E_MeV[mask]
sigma_b = sigma_b[mask]

if len(E_MeV) == 0:
    print("ERROR: TALYS file contains no positive energies.")
    exit()


E_J = E_MeV * MEV_TO_J
sigma_m2 = sigma_b * BARN_TO_M2


uniq_idx = np.unique(E_J, return_index=True)[1]
E_J = E_J[uniq_idx]
sigma_m2 = sigma_m2[uniq_idx]

sort_idx = np.argsort(E_J)
E_J = E_J[sort_idx]
sigma_m2 = sigma_m2[sort_idx]


sigma_interp = interp1d(E_J, sigma_m2, kind="linear",
                        bounds_error=False, fill_value=0.0)


E_min = max(E_J.min(), 1e-40)   
E_max = E_J.max()

E_plot = np.logspace(np.log10(E_min), np.log10(E_max), 2000)

kT_J = kT_keV * 1e3 * EV_TO_J
integrand = sigma_interp(E_plot) * E_plot * np.exp(-E_plot / kT_J)


E_plot_keV = E_plot / (1e3 * EV_TO_J)


plt.figure(figsize=(7.5, 4.8))
plt.loglog(E_plot_keV, integrand, lw=2,
           label=f"Integrand (kT={kT_keV} keV)")

plt.axvline(kT_keV, color="red", linestyle="--", label="kT")
plt.xlabel("Energy (keV)")
plt.ylabel("σ(E) E exp(-E/kT)")
plt.title(f"Maxwellian Integrand for ⁵⁷Fe(n,g)")
plt.grid(True, which='both', ls='--', alpha=0.4)
plt.legend()
plt.tight_layout()

out_png = "/home/zilani/1.Md_Zilani/Python/Project 2/Result/integrand_" \
          + talys_file.split('/')[-1].replace('.csv', '') + ".png"

plt.savefig(out_png, dpi=200)
print("Saved:", out_png)

plt.show()
