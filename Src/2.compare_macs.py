
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


talys_macs_csv = "/home/zilani/1.Md_Zilani/Python/Project 2/Data/Fe_MACS_TALYS.csv"
expt_macs_csv  = "/home/zilani/1.Md_Zilani/Python/Project 2/Data/Fe_MACS_Kd.csv"  # combined experimental file


talys = pd.read_csv(talys_macs_csv)
expt  = pd.read_csv(expt_macs_csv)


talys.columns = [c.strip() for c in talys.columns]
expt.columns  = [c.strip() for c in expt.columns]


cmp = talys.merge(expt, on=["Reaction","kT_keV"], suffixes=("_talys","_exp"))

if cmp.empty:
    print("No direct matches found. Check Reaction names and kT_keV values in both files.")
else:

    cmp['pct_diff'] = 100.0 * (cmp['MACS(b)_talys'] - cmp['MACS(b)_exp']) / cmp['MACS(b)_exp']

    with pd.option_context('display.max_rows', None, 'display.width', 120):
        print(cmp[['Reaction','kT_keV','MACS(b)_talys','MACS(b)_exp','pct_diff']])

   
    reactions = cmp['Reaction'].unique()
    for r in reactions:
        sub = cmp[cmp['Reaction']==r]
        plt.figure(figsize=(6,4))
        plt.scatter(sub['kT_keV'], sub['MACS(b)_exp'], label='KADoNiS MACS', marker='o')
        plt.scatter(sub['kT_keV'], sub['MACS(b)_talys'], label='TALYS MACS', marker='x')
        for i,row in sub.iterrows():
            plt.annotate(f"{row['pct_diff']:.1f}%", (row['kT_keV'], row['MACS(b)_talys']), textcoords="offset points", xytext=(5,-8))
        plt.xlabel("kT (keV)")
        plt.ylabel("MACS (b)")
        plt.title(r)
        plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()
       
        out_png = os.path.join("/home/zilani/1.Md_Zilani/Python/Project 2/Result/", f"macs_compare_{r.replace('/','_').replace(' ','')}.png")
        plt.savefig(out_png, dpi=150)
        print("Saved", out_png)
        plt.show()



