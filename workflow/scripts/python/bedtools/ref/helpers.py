from pathlib import Path
import common.config as cfg


def filter_sort_ref_outputs(refkeys: str, prefix: Path) -> dict[str, Path]:
    """Build the paths for the filter_ref* rules."""
    r = cfg.RefkeyConfiguration(refkeys)
    fa = r.value + "filtered_ref.fa"
    return {
        "fa": prefix / fa,
        "index": prefix / (fa + ".fai"),
        "genome": prefix / (r.value + "genome.txt"),
    }
