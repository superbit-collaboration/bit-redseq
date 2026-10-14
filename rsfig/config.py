"""Load the YAML configuration."""

import copy
import os

import yaml

DEFAULT_CONFIG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                              "config", "default.yaml")


def load_config(path=None, overrides=None):
    """
    Read a YAML config (default: config/default.yaml) and apply dotted
    overrides, e.g. {"selection.dz": 0.01, "plot.purity_ymin": 0.4}.
    """
    with open(path or DEFAULT_CONFIG) as f:
        cfg = copy.deepcopy(yaml.safe_load(f))
    for key, value in (overrides or {}).items():
        node = cfg
        *parents, leaf = key.split(".")
        for part in parents:
            node = node[part]
        if leaf not in node:
            raise KeyError(f"Unknown config key: {key}")
        node[leaf] = value
    return cfg
