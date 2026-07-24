"""Record the conditions a result was produced under.

Call `stamp()` at the top of every notebook and commit the JSON alongside the
results. A number without its conditions cannot be reproduced, and a result
nobody can reproduce cannot be contested.
"""

from __future__ import annotations

import json
import platform
import random
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

SEED = 0


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        return "unknown"


def _gpu() -> str:
    try:
        import torch  # noqa: PLC0415

        if torch.cuda.is_available():
            return f"{torch.cuda.get_device_name(0)} ({torch.__version__})"
        return f"cpu (torch {torch.__version__})"
    except ImportError:
        return "cpu (torch not installed)"


def stamp(experiment: str, out_dir: str | Path = "results", **extra) -> dict:
    """Seed everything, capture the environment, write a provenance record."""
    random.seed(SEED)
    try:
        import numpy as np  # noqa: PLC0415

        np.random.seed(SEED)
    except ImportError:
        pass
    try:
        import torch  # noqa: PLC0415

        torch.manual_seed(SEED)
    except ImportError:
        pass

    record = {
        "experiment": experiment,
        "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": _git_commit(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "device": _gpu(),
        "seed": SEED,
        "packages": _versions("pandas", "numpy", "scipy", "statsmodels", "transformers"),
        **extra,
    }

    out = Path(out_dir) / f"provenance_{experiment}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=2))
    return record


def _versions(*names: str) -> dict[str, str]:
    from importlib.metadata import PackageNotFoundError, version  # noqa: PLC0415

    out = {}
    for n in names:
        try:
            out[n] = version(n)
        except PackageNotFoundError:
            out[n] = "absent"
    return out
