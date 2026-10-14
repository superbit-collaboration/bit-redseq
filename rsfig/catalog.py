"""Read the MEGA colour-magnitude catalog and the SuperBIT target redshifts."""

import os

import numpy as np
import pandas as pd
from astropy.table import Table

import superbit_lensing

from . import COLORS

# Shipped with superbit-lensing: Name, ..., z (cluster redshift) for every target
DEFAULT_TARGETS = os.path.join(os.path.dirname(os.path.dirname(superbit_lensing.__file__)),
                               "data", "SuperBIT_target_galactic_coords.csv")

REQUIRED_COLUMNS = ["CLUSTER", "Z_best", "Z_source"] + [c for c, *_ in COLORS]


def _float(col):
    return np.ma.filled(np.ma.asarray(col, dtype=float), np.nan)


class MegaCatalog:
    """
    Plain-numpy view of the columns used by the analysis (read once).
    `extra` adds further columns, e.g. the magnitude used for the purity box.
    """

    def __init__(self, path, extra=()):
        self.path = os.path.abspath(path)
        tab = Table.read(self.path)
        wanted = list(dict.fromkeys(REQUIRED_COLUMNS + list(extra)))
        missing = [c for c in wanted if c not in tab.colnames]
        if missing:
            raise ValueError(f"{path} is missing columns: {missing}")
        self.n = len(tab)
        self.cluster = np.char.strip(np.asarray(tab["CLUSTER"]).astype(str))
        self.source = np.char.strip(np.asarray(tab["Z_source"]).astype(str))
        self.zbest = _float(tab["Z_best"])
        self.columns = {c: _float(tab[c]) for c in wanted
                        if c not in ("CLUSTER", "Z_source", "Z_best")}

    def __getitem__(self, name):
        return self.columns[name]

    def targets(self, exclude_prefix="COSMOS"):
        names = sorted(set(self.cluster))
        if exclude_prefix:
            names = [t for t in names if not t.upper().startswith(exclude_prefix.upper())]
        return names


def load_target_redshifts(path=None):
    """Dict target name -> cluster redshift."""
    df = pd.read_csv(path or DEFAULT_TARGETS)
    return dict(zip(df["Name"].str.strip(), df["z"]))


def check_redshifts(names, zmap):
    missing = [t for t in names if not (t in zmap and np.isfinite(zmap[t]) and zmap[t] > 0)]
    if missing:
        raise ValueError(f"No valid cluster redshift in the target list for: {missing}")
