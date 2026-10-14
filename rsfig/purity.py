"""
Red-sequence box purity per cluster and colour (panel e).

  members     : spec-z (DESI/NED), |Z_best - z_cluster| <= dz
  non-members : spec-z (DESI/NED), |Z_best - z_cluster| >  dz
  both with mag <= mag_max

  1. fit colour = a * mag + b to the members (start from the redder half,
     iterate with clip_fit-sigma clipping so the red sequence drives the fit)
  2. half-width  delta_c = tol_bg * std_c / std_(b-g), std = clip_sigma-clipped
     standard deviation of the members
  3. purity = N_mem,in / (N_mem,in + N_non,in)
"""

import numpy as np
from astropy.table import Table

from . import COLORS
from .members import clipped


def fit_red_sequence(mag, color, clip=2.5, maxiter=10):
    """Iteratively sigma-clipped straight-line fit; returns (slope, intercept)."""
    keep = color >= np.median(color)
    for _ in range(maxiter):
        a, b = np.polyfit(mag[keep], color[keep], 1)
        resid = color - (a * mag + b)
        sig = np.std(resid[keep], ddof=2)
        new = np.abs(resid) <= clip * sig
        if new.sum() < 3 or np.array_equal(new, keep):
            break
        keep = new
    return a, b


def purity_table(cat, zmap, sel, pur):
    mag = cat[pur["mag"]]
    rows = []
    for name in cat.targets(sel["exclude_prefix"]):
        zc = zmap[name]
        spec = ((cat.cluster == name) & np.isin(cat.source, sel["spec_sources"])
                & np.isfinite(cat.zbest) & np.isfinite(mag) & (mag <= pur["mag_max"]))
        mem = spec & (np.abs(cat.zbest - zc) <= sel["dz"])
        if mem.sum() < sel["min_members"]:
            continue
        non = spec & (np.abs(cat.zbest - zc) > sel["dz"])

        row = {"CLUSTER": name, "Z_MEMBERS_MEAN": float(cat.zbest[mem].mean()),
               "N_MEMBERS": int(mem.sum()), "N_NONMEMBERS": int(non.sum())}
        std_bg = None
        for col, *_ in COLORS:
            c = cat[col]
            m_ok, n_ok = mem & np.isfinite(c), non & np.isfinite(c)
            a, b = fit_red_sequence(mag[m_ok], c[m_ok], clip=pur["clip_fit"])
            res_m = np.abs(c[m_ok] - (a * mag[m_ok] + b))
            res_n = np.abs(c[n_ok] - (a * mag[n_ok] + b))

            std = float(np.std(clipped(c[m_ok], sel["clip_sigma"]), ddof=1))
            if col == "color_bg":
                std_bg = std
            tol = pur["tol_bg"] * std / std_bg

            n_in_m, n_in_n = int((res_m <= tol).sum()), int((res_n <= tol).sum())
            row.update({
                f"{col}_SLOPE": a, f"{col}_STD": std, f"{col}_TOL": tol,
                f"{col}_FRAC_IN": n_in_m / len(res_m),
                f"{col}_NMEM_IN": n_in_m, f"{col}_NNON_IN": n_in_n,
                f"{col}_PURITY": n_in_m / (n_in_m + n_in_n) if (n_in_m + n_in_n) else np.nan,
            })
        rows.append(row)
    return Table(rows=rows)
