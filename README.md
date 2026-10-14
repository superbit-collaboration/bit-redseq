# redseq-summary-figure

Reproduces the red-sequence summary figure of [arXiv:2610.06635](https://arxiv.org/abs/2610.06635)
from a SuperBIT MEGA colour–magnitude catalog.

## Requirements

[`superbit-lensing`](https://github.com/superbit-collaboration/superbit-lensing) (for `rs_catcher`
and the target redshifts), plus `pip install -r requirements.txt`.

## Usage

```bash
python make_figure.py --catalog FINAL_mega_color_mag_catalog_XXXX.fits
```

Writes `fig_redsequence_summary.pdf/.png`, the per-cluster tables and the settings used
to `./output/`. All settings are in [`config/default.yaml`](config/default.yaml); override
any of them with `--set`, e.g. `--set plot.purity_counts=true selection.dz=0.01`.
