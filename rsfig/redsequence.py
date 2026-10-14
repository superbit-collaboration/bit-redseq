"""
Red-sequence selection and density map of one cluster (panels a, b), using
superbit_lensing.rs_catcher.ClusterRedSequenceAnalysis without modifying it.

compute_red_sequence() in rs_catcher also makes and saves its own figure, so
the same steps are run here one by one, minus the plotting.
"""

import contextlib
import os

import numpy as np

from superbit_lensing.rs_catcher import ClusterRedSequenceAnalysis


@contextlib.contextmanager
def _cwd(path):
    old = os.getcwd()
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


def run_red_sequence(catalog_path, rs, workdir):
    """
    Returns the fitted ClusterRedSequenceAnalysis. load_data() writes a
    per-cluster catalog and some folders into the working directory, so it
    runs inside `workdir`.
    """
    with _cwd(workdir):
        an = ClusterRedSequenceAnalysis(rs["cluster"], megafilename=os.path.abspath(catalog_path),
                                        delz=rs["delz"], color_index=rs["color_index"],
                                        only_specz=True)
        an.load_data()
        an.m_lower_limit, an.m_upper_limit = rs["m_lower_limit"], rs["m_upper_limit"]
        an.fit_red_sequence_line(fixed_intercept=rs["fixed_intercept"])
        an.tolerance = rs["tolerance"]
        an.ra_max, an.ra_min = np.max(an.cm_cat["ra"]), np.min(an.cm_cat["ra"])
        an.dec_max, an.dec_min = np.max(an.cm_cat["dec"]), np.min(an.cm_cat["dec"])
        pred = an.a * an.m_b + an.b
        an.red_sequence_mask = ((np.abs(an.color_index - pred) < an.tolerance)
                                & (an.m_b < an.m_upper_limit) & (an.m_b > an.m_lower_limit))
        an.ra_red = an.cm_cat["ra"][an.red_sequence_mask]
        an.dec_red = an.cm_cat["dec"][an.red_sequence_mask]
        an.m_b_red = an.cm_cat["m_b"][an.red_sequence_mask]
        an.create_density_map(resolution=rs["resolution"], sigma=rs["sigma"],
                              mode=rs["density_mode"])
    print(f"Red sequence ({rs['cluster']}): color = {an.a:.3f} m_b + {an.b:.3f}, "
          f"{int(an.red_sequence_mask.sum())} objects")
    return an
