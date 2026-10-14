"""Spectroscopic cluster members and their colour statistics (panels c, d, f)."""

import numpy as np
from astropy.stats import sigma_clip
from astropy.table import Table

from . import COLORS


def member_mask(cat, name, zc, sel):
    """DESI/NED objects of target `name` within +-dz of the cluster redshift."""
    return ((cat.cluster == name) & np.isin(cat.source, sel["spec_sources"])
            & np.isfinite(cat.zbest) & (np.abs(cat.zbest - zc) <= sel["dz"]))


def clipped(values, clip_sigma):
    """Finite values that survive iterative sigma clipping."""
    v = values[np.isfinite(values)]
    if len(v) < 2:
        return v
    return v[~sigma_clip(v, sigma=clip_sigma, maxiters=None).mask]


def member_colors(cat, name, zc, sel):
    """Clipped member colours of one cluster, in COLORS order (violin panel)."""
    mem = member_mask(cat, name, zc, sel)
    return [clipped(cat[col][mem], sel["clip_sigma"]) for col, *_ in COLORS]


def scatter_table(cat, zmap, sel):
    """
    One row per (non-excluded) target: member count, mean member redshift and,
    per colour, the clipped std / mean / median (+ errors) of the members.
    """
    rows = []
    for name in cat.targets(sel["exclude_prefix"]):
        zc = zmap[name]
        mem = member_mask(cat, name, zc, sel)
        row = {"CLUSTER": name, "Z_CLUSTER": zc, "N_MEMBERS": int(mem.sum()),
               "Z_MEMBERS_MEAN": float(np.mean(cat.zbest[mem])) if mem.any() else np.nan}
        for col, *_ in COLORS:
            c = cat[col][mem]
            c = c[np.isfinite(c)]
            if len(c) >= 2:
                kept = clipped(c, sel["clip_sigma"])
                std = float(np.std(kept, ddof=1)) if len(kept) >= 2 else np.nan
                sem = std / np.sqrt(len(kept))
                row.update({
                    f"{col}_N": len(kept), f"{col}_NCLIP": len(c) - len(kept),
                    f"{col}_STD": std, f"{col}_MEAN": float(np.mean(kept)), f"{col}_SEM": sem,
                    f"{col}_MEDIAN": float(np.median(kept)),
                    # standard error of the median, large-N Gaussian approximation
                    f"{col}_MEDIAN_ERR": np.sqrt(np.pi / 2) * sem,
                })
            else:
                mean = float(c[0]) if len(c) else np.nan
                row.update({
                    f"{col}_N": len(c), f"{col}_NCLIP": 0, f"{col}_STD": np.nan,
                    f"{col}_MEAN": mean, f"{col}_SEM": np.nan,
                    f"{col}_MEDIAN": mean, f"{col}_MEDIAN_ERR": np.nan,
                })
        row["FLAG_LOW_N"] = row["N_MEMBERS"] < sel["min_members"]
        rows.append(row)
    return Table(rows=rows)
