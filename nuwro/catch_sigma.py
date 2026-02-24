#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Scan NuWro cross‑sections for all WEFT cases (0‑10) and
selected sub‑configurations. Only cases 0,1,2,7,8 will scan
weft_param = 2,3; other cases scan only 0,1.
"""

import subprocess
import re
import pathlib
import sys
import numpy as np

PARAM_FILE = pathlib.Path("params.txt")        
NUWRO_EXEC = pathlib.Path("nuwro")              
OUT_TABLE  = pathlib.Path("weft_sigma_mu.txt")     

# ENERGIES_GEV  = [0.2, 0.25, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2, 2.8, 3.6,
#                  4.4, 5.2, 6, 7, 8, 9, 10]
# ENERGIES_GEV  = [0.2, 0.4, 0.6, 1.0, 1.6, 2, 3, 4, 5, 6, 7, 8, 9, 10]
# ENERGIES_GEV  = [0.15, 0.2, 0.4, 0.6, 1.0, 3, 6, 10]
ENERGIES_GEV  = [0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10]
WEFT_INDIVIDUAL_CASE    = list(range(11))                  
# WEFT_INDIVIDUAL_CASE    = [10]  
SPECIAL_CASE = {0, 1, 2, 7, 8}                 
ENERGIES_MEV  = (np.array(ENERGIES_GEV) * 1000.0).astype(float)
 
SIGMA_RE = re.compile(r"([-+]?\d\.\d+e[-+]\d+)\s*cm2", re.I)

try:
    ORIG_PARAM = PARAM_FILE.read_text().splitlines()
except FileNotFoundError:
    sys.exit(f"Parameter file '{PARAM_FILE}' not found. Abort.")

def params_for(case: int) -> list[int]:
    """Return allowed weft_param values for a given case."""
    return [0, 1, 2, 3] if case in SPECIAL_CASE else [0, 1]
    # return [0] if case in SPECIAL_CASES else [0]

def rewrite_param(lines: list[str], *,
                  energy: float | None = None,
                  weft_individual_case: int | None = None,
                  weft_param: int | None = None) -> list[str]:
    """Return a *new* list with requested fields replaced."""
    out: list[str] = []
    found_case = found_param = False

    for ln in lines:
        s = ln.strip()
        if energy is not None and s.startswith("beam_energy"):
            out.append(f"beam_energy   = {energy:.5f}")
        elif weft_individual_case is not None and s.startswith("weft_individual_case"):
            out.append(f"weft_individual_case = {weft_individual_case}")
            found_case = True
        elif weft_param is not None and s.startswith("weft_param"):
            out.append(f"weft_param = {weft_param}")
            found_param = True
        else:
            out.append(ln)

    if weft_individual_case is not None and not found_case:
        out.append(f"weft_individual_case = {weft_individual_case}")
    if weft_param is not None and not found_param:
        out.append(f"weft_param = {weft_param}")

    return out

col_labels = [f"sigma_case{c}_p{p}"
              for c in WEFT_INDIVIDUAL_CASE
              for p in params_for(c)]
header = "Energy [MeV]\t" + "\t".join(col_labels) + "  [cm^-2]\n"
OUT_TABLE.write_text(header)

for Ev in ENERGIES_MEV:
    row = [f"{Ev:.1f}"] 

    for case in WEFT_INDIVIDUAL_CASE:
        for param in params_for(case):
            new_lines = rewrite_param(ORIG_PARAM,
                                      energy=Ev,
                                      weft_individual_case=case,
                                      weft_param=param)
            PARAM_FILE.write_text("\n".join(new_lines))

            res = subprocess.run(
                [str(NUWRO_EXEC), "-i", str(PARAM_FILE), "-o", "dummy.root"],
                text=True, capture_output=True)

            m = SIGMA_RE.search(res.stdout)
            row.append(m.group(1) if m else "NaN")
            if not m:
                print(f"[warn] Ev={Ev:.1f} MeV  case={case}  param={param} — σ not found!")

    with OUT_TABLE.open("a") as f:
        f.write("\t".join(row) + "\n")

    print(f"✓ {Ev/1000:.2f} GeV done")

