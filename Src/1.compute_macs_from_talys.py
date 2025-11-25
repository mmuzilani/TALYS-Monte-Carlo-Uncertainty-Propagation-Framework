import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import quad
import os
import math


EV_TO_J = 1.602176634e-19
MEV_TO_J = EV_TO_J * 1e6
BARN_TO_M2 = 1e-28
sqrt_pi = math.sqrt(math.pi)



sup_map = str.maketrans({
    "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
    "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹"
})
def to_superscript(x):
    return str(x).translate(sup_map)



def macs_from_arrays(energy_MeV, sigma_b, kT_keV):
   
  
    E_J = np.array(energy_MeV, dtype=float) * MEV_TO_J
    sigma_m2 = np.array(sigma_b, dtype=float) * BARN_TO_M2
    kT_J = kT_keV * 1e3 * EV_TO_J


    E_J, unique_idx = np.unique(E_J, return_index=True)
    sigma_m2 = sigma_m2[unique_idx]

 
    sort_idx = np.argsort(E_J)
    E_J = E_J[sort_idx]
    sigma_m2 = sigma_m2[sort_idx]

    
    sigma_interp = interp1d(
        E_J, sigma_m2, kind="linear",
        bounds_error=False, fill_value=0.0
    )

   
    def integrand(E):
        return float(sigma_interp(E)) * E * np.exp(-E / kT_J)

    E_min = max(E_J.min(), 1e-40)
    E_max = E_J.max()


    integral, err = quad(integrand, E_min, E_max, limit=400)

    prefactor = (2.0 / sqrt_pi) / (kT_J * kT_J)
    macs_m2 = prefactor * integral
    macs_b = macs_m2 / BARN_TO_M2

    return macs_b


def compute_macs_for_files(talys_files, kT_list_keV, out_csv_path):
    rows = []

    for fpath in talys_files:
        df = pd.read_csv(fpath)
        df.columns = [c.strip() for c in df.columns]

        if 'Energy(MeV)' not in df.columns or 'CrossSection(b)' not in df.columns:
            raise ValueError(f"File {fpath} missing required columns!")

        E = df['Energy(MeV)'].values
        sigma = df['CrossSection(b)'].values

   
        base = os.path.basename(fpath)
        target = base.split('_')[0]     # e.g. "56Fe"

        import re
        m = re.match(r'(\d+)([A-Za-z]+)', target)
        if m:
            A = int(m.group(1))
            Z = m.group(2)

            target_label = f"{to_superscript(A)}{Z}"      
            product_label = f"{to_superscript(A+1)}{Z}"   
        else:
            target_label = target
            product_label = target + "_p"

        reaction_label = f"{target_label}(n,g){product_label}"


        for kT in kT_list_keV:
            macs_b = macs_from_arrays(E, sigma, kT)
            rows.append({
                "Reaction": reaction_label,
                "kT_keV": kT,
                "MACS(b)": macs_b
            })
            print(f"{reaction_label}   kT={kT:>4} keV  ->  MACS = {macs_b:.6e} b")

    out_df = pd.DataFrame(rows)
    out_df.to_csv(out_csv_path, index=False)
    print("\nSaved MACS to:", out_csv_path)
    return out_df


if __name__ == "__main__":

    talys_files = [
        "/home/zilani/1.Md_Zilani/Python/Project 2/Data/56Fe_ng_TALYS.csv",
        "/home/zilani/1.Md_Zilani/Python/Project 2/Data/57Fe_ng_TALYS.csv",
        "/home/zilani/1.Md_Zilani/Python/Project 2/Data/58Fe_ng_TALYS.csv",
        "/home/zilani/1.Md_Zilani/Python/Project 2/Data/59Fe_ng_TALYS.csv"
    ]

    kT_list = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0]  # keV

    out_csv = "/home/zilani/1.Md_Zilani/Python/Project 2/Data/Fe_MACS_TALYS.csv"

    compute_macs_for_files(talys_files, kT_list, out_csv)

