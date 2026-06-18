from __future__ import annotations

import argparse
import tomllib
from pathlib import Path


def load_toml(path: Path) -> dict:
    with open(path, "rb") as f:
        return tomllib.load(f)


def apply_toml_defaults(
    parser: argparse.ArgumentParser,
    toml_key_to_dest: dict[str, str],
    path_dests: set[str] | None = None,
) -> None:
    pre, _ = parser.parse_known_args()
    if getattr(pre, "experiment", None) is None:
        return
    cfg = load_toml(pre.experiment)
    path_dests = path_dests or set()
    defaults = {}
    for toml_key, dest in toml_key_to_dest.items():
        parts = toml_key.split(".")
        node = cfg
        for part in parts[:-1]:
            node = node.get(part, {})
        if parts[-1] in node:
            val = node[parts[-1]]
            if dest in path_dests:
                val = Path(val)
            defaults[dest] = val
    if defaults:
        parser.set_defaults(**defaults)


def add_experiment_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "experiment", nargs="?", type=Path, default=None,
        help="TOML experiment file.",
    )
