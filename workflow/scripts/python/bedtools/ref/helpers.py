from pathlib import Path
import common.config as cfg


def filter_sort_ref_outputs(refkeys: str, bgzip: bool, prefix: Path) -> dict[str, Path]:
    """Build the paths for the filter_ref* rules."""
    r = cfg.RefkeyConfiguration(refkeys)
    fa_ext, idx_ext, genome_stem = (".gz", ".gzi", "gz_") if bgzip else ("", ".fai", "")
    fa = r.value + "filtered_ref.fa" + fa_ext
    return {
        "fa": prefix / fa,
        "index": prefix / (fa + idx_ext),
        "genome": prefix / (genome_stem + "genome.txt"),
    }
