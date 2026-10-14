import numpy as np

from rsfig.config import load_config
from rsfig.plotting import spread_positions
from rsfig.purity import fit_red_sequence


def test_spread_positions_keeps_order_mean_and_spacing():
    x = np.array([0.0920, 0.0942, 0.0943, 0.0995, 0.30])
    s = spread_positions(x, 0.0115)
    assert np.all(np.diff(s[np.argsort(x)]) >= 0.0115 - 1e-12)
    assert np.isclose(s[:4].mean(), x[:4].mean())  # crowded group stays centred
    assert s[4] == x[4]                            # isolated point untouched


def test_fit_red_sequence_ignores_blue_cloud():
    rng = np.random.default_rng(1)
    mag = rng.uniform(18, 22, 300)
    red = 1.0 - 0.02 * (mag - 20) + rng.normal(0, 0.03, 300)
    blue = rng.normal(0.4, 0.2, 300)
    a, b = fit_red_sequence(np.r_[mag, mag], np.r_[red, blue])
    assert abs(a + 0.02) < 0.01
    assert abs((a * 20 + b) - 1.0) < 0.02


def test_config_override():
    cfg = load_config(overrides={"selection.dz": 0.01})
    assert cfg["selection"]["dz"] == 0.01
    assert cfg["purity"]["tol_bg"] == 0.1


if __name__ == "__main__":  # run without pytest: python tests/test_rsfig.py
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"ok  {name}")
