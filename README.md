# redseq-summary-figure

Reproduces the SuperBIT red-sequence summary figure (`fig_redsequence_summary.pdf`)
from a MEGA colour–magnitude catalog.

| Panel | Content |
|---|---|
| (a) | Abell 1689 colour–magnitude diagram with the red-sequence selection |
| (b) | Smoothed number density of the red-sequence galaxies |
| (c) | b−g, u−g, u−b distributions of the Abell 1689 spectroscopic members |
| (d) | Colour scatter of spectroscopic members vs redshift (numbers: N members) |
| (e) | Red-sequence box purity vs redshift |
| (f) | Median member colour vs redshift |

Panels (a) and (b) use `superbit_lensing.rs_catcher.ClusterRedSequenceAnalysis`
(the package is called, not modified).

## Requirements

- [`superbit-lensing`](https://github.com/superbit-collaboration/superbit-lensing) installed
  (it provides `rs_catcher` and the target list with cluster redshifts)
- Python packages: see `requirements.txt`

The catalog must contain at least `CLUSTER`, `ra`, `dec`, `Z_best`, `Z_source`,
`m_b`, `m_g`, `color_bg`, `color_ug`, `color_ub`, `FLUX_AUTO_{u,b,g}` and `R_b`.

## Usage

```bash
python make_figure.py --catalog /path/to/FINAL_mega_color_mag_catalog_XXXX.fits
```

Outputs go to `./output/` (change with `--outdir`):

| File | |
|---|---|
| `fig_redsequence_summary.pdf` / `.png` | the figure |
| `cluster_color_scatter_specz.fits` | per-cluster member colour statistics (panels d, f) |
| `cluster_redsequence_purity.fits` | per-cluster purity, box widths, counts (panel e) |
| `config_used.yaml` | exact settings of the run |
| `rs_work/` | scratch files written by `rs_catcher.load_data()` |

Other options:

```bash
--targets my_targets.csv      # Name, z columns; default: superbit-lensing/data/SuperBIT_target_galactic_coords.csv
--config my_config.yaml       # default: config/default.yaml
--set selection.dz=0.01 plot.purity_ymin=0.4   # override single settings
--set plot.purity_counts=true  # label (e) with members / non-members per cluster
```

By default panel (e) has no count labels and its axis starts at `plot.purity_ymin`
(0.5), so clusters with lower purity (MACSJ0416, ~0.2) fall outside it. With
`plot.purity_counts: true` the labels are drawn and the axis extends to show every
cluster.

## Method

All settings live in [`config/default.yaml`](config/default.yaml).

**Spectroscopic members** of a target: objects with `Z_source` in {DESI, NED} and
|`Z_best` − z_cluster| ≤ 0.02, with z_cluster from the target list. COSMOS fields are
skipped; clusters with fewer than 10 members are left out of panels (d)–(f).

**Colour scatter (d), median colour (f), violins (c):** per colour, members are
5σ-clipped (iteratively); the scatter is the standard deviation of the remaining
members, the median error is √(π/2)·σ/√N.

**Purity (e):** members and non-members (DESI/NED with |Δz| > 0.02), both with
m_g ≤ 22.5. For each colour c a line c = a·m_g + b is fitted to the members (starting
from the redder half, 2.5σ clipping). A galaxy is in the box if
|c − (a·m_g + b)| ≤ δ_c, with δ_c = 0.1 × σ_c / σ_(b−g), so δ = 0.1 mag for b−g.

  P_c = N_mem,in / (N_mem,in + N_non,in)

**Red-sequence selection (a, b):** `rs_catcher` with a flat line fitted to the
Abell 1689 members in b−g vs b, ±0.1 mag, 16 < b < 23; number density on 0.47′ pixels,
Gaussian smoothing σ = 1.51 pixels.

In panels (d)–(f), clusters closer than 0.0115 in redshift are displaced slightly so
their labels don't overlap; the true redshifts are in the tables.

## Tests

```bash
python -m pytest tests      # or, without pytest:
PYTHONPATH=. python tests/test_rsfig.py
```
