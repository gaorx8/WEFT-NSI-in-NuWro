#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

mpl.rcParams.update({
    "figure.figsize": (3.45, 2.60),
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "font.size": 9,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.minor.size": 2.0,
    "ytick.minor.size": 2.0,
})

INFILE = "weft_sigma_mu.txt"
OUT_DIR = Path("figs_test")
OUT_DIR.mkdir(exist_ok=True)

EPS_R = 4e-2
EPS_R2 = EPS_R**2

COL_LL = "sigma_case0_p1"   # SM
COL_RR = "sigma_case1_p1"   # quadratic
COL_LR = "sigma_case2_p1"   # interference

USE_MINUS_DLR = False 

with open(INFILE, "r", encoding="utf-8") as f:
    header_line = f.readline().strip()

tokens = header_line.split()

if len(tokens) >= 2 and tokens[0] == "Energy" and tokens[1].startswith("[MeV"):
    tokens = ["Energy"] + tokens[2:]

if len(tokens) >= 1 and tokens[-1].startswith("[cm"):
    tokens = tokens[:-1]

colnames = tokens

df = pd.read_csv(
    INFILE,
    sep=r"\s+",
    engine="python",
    header=None,
    names=colnames,
    skiprows=1,
)

need_cols = ["Energy", COL_LL, COL_LR, COL_RR]
missing = [c for c in need_cols if c not in df.columns]
if missing:
    raise RuntimeError(
        f"Missing columns {missing}. Parsed columns are:\n{list(df.columns)}\n"
        f"First header tokens were:\n{header_line}"
    )

for c in need_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

E_GEV = df["Energy"].to_numpy(dtype=float) / 1000.0
LL = df[COL_LL].to_numpy(dtype=float)
LR = df[COL_LR].to_numpy(dtype=float)
RR = df[COL_RR].to_numpy(dtype=float)

mask = (
    np.isfinite(E_GEV)
    & np.isfinite(LL)
    & np.isfinite(LR)
    & np.isfinite(RR)
    & (LL != 0.0)
)
E = E_GEV[mask]
LL = LL[mask]
LR = LR[mask]
RR = RR[mask]

if E.size < 3:
    raise RuntimeError(f"After masking, only {E.size} points remain.")

order = np.argsort(E)
E, LL, LR, RR = E[order], LL[order], LR[order], RR[order]

sgn = -1.0 if USE_MINUS_DLR else +1.0
term_LR = sgn * EPS_R * LR
term_RR = EPS_R2 * RR
sigma_tot = LL + term_LR + term_RR
ratio = sigma_tot / LL
delta = ratio - 1.0

print("\n=== FULL Sanity check (first 10 points after mask/sort) ===")
print(f"eps = {EPS_R:.6f}, eps^2 = {EPS_R2:.6f}, sgn(LR) = {sgn:+.0f}")
print(
    " E[GeV]   "
    "LL           "
    "LR           "
    "RR           "
    "eps*LR       "
    "eps^2*RR     "
    "sigma_tot    "
    "ratio        "
    "delta"
)
for k in range(E.size):
    print(
        f"{E[k]:6.3f}   "
        f"{LL[k]: .6e}  "
        f"{LR[k]: .6e}  "
        f"{RR[k]: .6e}  "
        f"{term_LR[k]: .6e}  "
        f"{term_RR[k]: .6e}  "
        f"{sigma_tot[k]: .6e}  "
        f"{ratio[k]: .6f}  "
        f"{delta[k]: .6e}"
    )

print("\n=== Ranges over all energies ===")
print(f"E      : [{np.nanmin(E):.3f}, {np.nanmax(E):.3f}] GeV, N={E.size}")
print(f"LL     : [{np.nanmin(LL):.3e}, {np.nanmax(LL):.3e}]")
print(f"LR     : [{np.nanmin(LR):.3e}, {np.nanmax(LR):.3e}]")
print(f"RR     : [{np.nanmin(RR):.3e}, {np.nanmax(RR):.3e}]")
print(f"eps*LR : [{np.nanmin(term_LR):.3e}, {np.nanmax(term_LR):.3e}]")
print(f"eps^2*RR: [{np.nanmin(term_RR):.3e}, {np.nanmax(term_RR):.3e}]")
print(f"ratio  : [{np.nanmin(ratio):.6f}, {np.nanmax(ratio):.6f}]")
print(f"delta  : [{np.nanmin(delta):.3e}, {np.nanmax(delta):.3e}]")

fig, ax = plt.subplots()
ax.plot(E, delta, lw=1.4)

ax.set_xlabel(r"$E_\nu\ \mathrm{[GeV]}$")
ax.set_ylabel(
    r"$(\sigma_{\rm SM+NSI}-\sigma_{\rm SM})/\sigma_{\rm SM}$"
)


ax.minorticks_on()
ax.set_xlim(0.0, 10.0)
ax.margins(x=0.0)

ymin = float(np.nanmin(ratio))
ymax = np.max(np.abs(delta))
pad = 0.15 * ymax if ymax > 0 else 1e-4
ax.set_ylim(-ymax - pad, ymax + pad)

fig.tight_layout(pad=0.2)

tag = "minusdLR" if USE_MINUS_DLR else "plusLR"
out_pdf = OUT_DIR / f"ratio_numu_R_sig_over_SM_{tag}_p1.pdf"
out_png = OUT_DIR / f"ratio_numu_R_sig_over_SM_{tag}_p1.png"
fig.savefig(out_pdf)
fig.savefig(out_png)
plt.close(fig)

print("\n✓ Saved:", out_pdf)
print("✓ Saved:", out_png)
