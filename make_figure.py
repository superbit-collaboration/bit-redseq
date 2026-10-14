#!/usr/bin/env python
"""
Build the red-sequence summary figure from a MEGA colour-magnitude catalog.

    python make_figure.py --catalog FINAL_mega_color_mag_catalog_XXXX.fits

Writes into --outdir (default ./output):
    fig_redsequence_summary.pdf/.png                    the figure (count labels in (e)
                                                        only if plot.purity_counts: true)
    cluster_color_scatter_specz.fits                    per-cluster member colour stats
    cluster_redsequence_purity.fits                     per-cluster purity
"""

import argparse
import os

import yaml

from rsfig.catalog import MegaCatalog, load_target_redshifts, check_redshifts
from rsfig.config import load_config
from rsfig.members import member_colors, scatter_table
from rsfig.purity import purity_table
from rsfig.redsequence import run_red_sequence
from rsfig.plotting import make_figure


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalog", required=True, help="MEGA colour-magnitude FITS catalog")
    ap.add_argument("--targets", default=None,
                    help="target list CSV with Name and z columns "
                         "(default: superbit-lensing/data/SuperBIT_target_galactic_coords.csv)")
    ap.add_argument("--config", default=None, help="YAML config (default: config/default.yaml)")
    ap.add_argument("--outdir", default="output", help="output directory (default: ./output)")
    ap.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE",
                    help="override config values, e.g. --set selection.dz=0.01 plot.purity_ymin=0.4")
    return ap.parse_args()


def main():
    args = parse_args()
    overrides = dict(kv.split("=", 1) for kv in args.set)
    cfg = load_config(args.config, {k: yaml.safe_load(v) for k, v in overrides.items()})
    sel, pur, p = cfg["selection"], cfg["purity"], cfg["plot"]
    os.makedirs(args.outdir, exist_ok=True)

    cat = MegaCatalog(args.catalog, extra=[pur["mag"]])
    zmap = load_target_redshifts(args.targets)
    names = cat.targets(sel["exclude_prefix"])
    check_redshifts(names, zmap)
    print(f"{cat.path}: {cat.n} objects, {len(names)} targets")

    scatter = scatter_table(cat, zmap, sel)
    purity = purity_table(cat, zmap, sel, pur)
    scatter.write(os.path.join(args.outdir, "cluster_color_scatter_specz.fits"), overwrite=True)
    purity.write(os.path.join(args.outdir, "cluster_redsequence_purity.fits"), overwrite=True)
    print(f"{(~scatter['FLAG_LOW_N']).sum()} clusters with >= {sel['min_members']} members")

    vname = cfg["violin_cluster"]
    violin = member_colors(cat, vname, zmap[vname], sel)
    an = run_red_sequence(args.catalog, cfg["red_sequence"], os.path.join(args.outdir, "rs_work"))

    make_figure(an, violin, scatter, purity, p, os.path.join(args.outdir, "fig_redsequence_summary"),
                show_counts=p["purity_counts"])

    with open(os.path.join(args.outdir, "config_used.yaml"), "w") as f:
        yaml.safe_dump({"catalog": os.path.abspath(args.catalog),
                        "targets": args.targets or "superbit-lensing default", **cfg},
                       f, sort_keys=False)


if __name__ == "__main__":
    main()
