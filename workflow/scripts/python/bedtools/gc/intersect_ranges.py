from pathlib import Path
from typing import Any, NamedTuple, Callable
import common.config as cfg
from common.functional import DesignError
from common.io import bgzip_file, gunzip, check_processes
from common.bed import (
    complementBed,
    subtractBed,
    intersectBed,
    multiIntersectBed,
    mergeBed,
)
import json


class GCInput(NamedTuple):
    bed: Path
    fraction: int
    is_range_bound: bool


def write_simple_range_beds(
    final_path: Callable[[str], Path],
    gs: list[GCInput],
    genome: Path,
    is_low: bool,
    log: Path,
) -> list[str]:
    def fmt_out(bigger_frac: int, smaller_frac: int) -> Path:
        lower_frac, upper_frac = (
            (smaller_frac, bigger_frac) if is_low else (bigger_frac, smaller_frac)
        )
        return final_path(f"gc{lower_frac}to{upper_frac}_slop50")

    pairs = zip(gs[1:], gs[:-1]) if is_low else zip(gs[:-1], gs[1:])
    torun = [
        (fmt_out(bigger.fraction, smaller.fraction), bigger, smaller)
        for bigger, smaller in pairs
    ]
    for out, bigger, smaller in torun:
        p1, o1 = gunzip(bigger.bed)
        p2, o2 = subtractBed(o1, smaller.bed, genome)
        bgzip_file(o2, out)
        check_processes([p1, p2], log)
    return [str(t[0]) for t in torun]


def write_middle_range_bed(
    final_path: Callable[[str], Path],
    lower: GCInput,
    upper: GCInput,
    genome: Path,
    gapless: Path,
    log: Path,
) -> str:
    out = final_path(f"gc{lower.fraction}to{upper.fraction}_slop50")
    with open(upper.bed, "rb") as i:
        p1, o0 = complementBed(i, genome)
        p2, o1 = subtractBed(o0, lower.bed, genome)
        p3, o2 = intersectBed(o1, gapless, genome)
        bgzip_file(o2, out)
        check_processes([p1, p2, p3], log)
    return str(out)


def write_intersected_range_beds(
    final_path: Callable[[str], Path],
    low: list[GCInput],
    high: list[GCInput],
    log: Path,
) -> list[str]:
    pairs = zip(
        [x for x in low if x.is_range_bound],
        [x for x in reversed(high) if x.is_range_bound],
    )
    torun = [
        (i, final_path(f"gclt{b1.fraction}orgt{b2.fraction}_slop50"), b1, b2)
        for i, (b1, b2) in enumerate(pairs)
    ]
    for i, bed_out, b1, b2 in torun:
        p1, o1 = multiIntersectBed([b1.bed, b2.bed])
        p2, o2 = mergeBed(o1, [])
        bgzip_file(o2, bed_out)
        check_processes([p1, p2], log)
    return [str(t[1]) for t in torun]


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    path_pattern = cfg.smk_to_param_str(smk, "path_pattern")

    low_paths = cfg.smk_to_inputs_name(smk, "low")
    high_paths = cfg.smk_to_inputs_name(smk, "high")
    genome_path = cfg.smk_to_input_name(smk, "genome")
    gapless_path = cfg.smk_to_input_name(smk, "gapless")
    log_path = cfg.smk_to_log(smk)
    out_path = cfg.smk_to_output(smk)

    rfk = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)
    bd = sconf.to_build_data(cfg.strip_full_refkey(rfk), bk)
    gps = bd.build.include.gc
    if gps is None:
        raise DesignError
    # ASSUME both of these input lists are sorted by GC fraction
    low = [GCInput(p, f, r) for p, (f, r) in zip(low_paths, gps.low)]
    high = [GCInput(p, f, r) for p, (f, r) in zip(high_paths, gps.high)]
    genome = Path(genome_path)
    gapless = Path(gapless_path)

    def final_path(name: str) -> Path:
        p = Path(path_pattern.format(name))
        p.parent.mkdir(exist_ok=True, parents=True)
        return p

    low_strats = write_simple_range_beds(final_path, low, genome, True, log_path)
    high_strats = write_simple_range_beds(final_path, high, genome, False, log_path)
    range_strat = write_middle_range_bed(
        final_path,
        low[-1],
        high[0],
        genome,
        gapless,
        log_path,
    )
    inter_strats = write_intersected_range_beds(
        final_path,
        low,
        high,
        log_path,
    )

    # ASSUME there is one extreme denoted in the config (otherwise we won't have
    # something to put in "widest extreme"
    with open(out_path, "w") as f:
        # put the first low and last high input here since these are already
        # in the final directory
        ranges: list[str] = [
            str(low[0].bed),
            *low_strats,
            range_strat,
            *high_strats,
            str(high[-1].bed),
        ]
        obj = {
            "gc_ranges": ranges,
            "widest_extreme": inter_strats[0],
            "other_extremes": inter_strats[1:],
        }
        json.dump(obj, f, indent=4)


main(snakemake, snakemake.config)  # type: ignore
