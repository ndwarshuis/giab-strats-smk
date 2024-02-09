from typing import Any
from pathlib import Path
import common.config as cfg
from common.bed import (
    write_bed,
    bed_to_stream,
    read_bed_default,
    complementBed,
    mergeBed,
    subtractBed,
)
from common.io import check_processes
from common.functional import match1_unsafe, match2_unsafe
import pandas as pd


def read_genome_bed(p: Path) -> pd.DataFrame:
    df = pd.read_table(
        p,
        header=None,
        names=["chrom", "end"],
        dtype={"chrom": str, "end": int},
    )
    df["start"] = 0
    return df[["chrom", "start", "end"]].copy()


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    inputs = smk.input
    auto_out = Path(smk.output["auto"])
    parY_out = Path(smk.output["parY"])
    log = Path(smk.log[0])
    ws: dict[str, Any] = smk.wildcards

    def go(
        x: cfg.BuildData_[cfg.RefSourceT, cfg.AnyBedT, cfg.AnyBedT_, cfg.AnySrcT]
    ) -> cfg.BedFile[cfg.AnyBedT] | None:
        return x.refdata.strat_inputs.gap

    # convert genome to bed file (where each region is just the length of
    # one chromosome)
    genome_path = Path(inputs["genome"])
    genome_bed = read_genome_bed(genome_path)

    # If we have gap input, make the gapless file, otherwise just symlink to the
    # genome bed file (which just means the entire genome is gapless)
    if not hasattr(inputs, "gaps"):
        write_bed(auto_out, genome_bed)
        parY_out.symlink_to(auto_out.resolve())
    else:
        gap_inputs: list[Path] = inputs["gaps"]
        rk = cfg.wc_to_reffinalkey(ws)
        bk = cfg.wc_to_buildkey(ws)

        bd = sconf.to_build_data(cfg.strip_full_refkey(rk), bk)
        gaps_df = sconf.with_build_data_and_bed_full(
            rk,
            bk,
            go,
            lambda bd, bf: match1_unsafe(
                gap_inputs, lambda i: cfg.read_filter_sort_hap_bed(bd, bf, i)
            ),
            lambda bd, bf: match1_unsafe(
                gap_inputs, lambda i: cfg.read_filter_sort_dip1to1_bed(bd, bf, i)
            ),
            lambda hap, bd, bf: match1_unsafe(
                gap_inputs,
                lambda i: hap.from_either(*cfg.read_filter_sort_dip1to2_bed(bd, bf, i)),
            ),
            lambda bd, bf: match2_unsafe(
                gap_inputs,
                lambda i0, i1: cfg.read_filter_sort_dip2to1_bed(bd, bf, (i0, i1)),
            ),
            lambda hap, bd, bf: match2_unsafe(
                gap_inputs,
                lambda i0, i1: cfg.read_filter_sort_dip2to2_bed(
                    bd,
                    bf,
                    *hap.from_either((i0, cfg.Haplotype.HAP1), (i1, cfg.Haplotype.HAP2))
                ),
            ),
        )

        with bed_to_stream(gaps_df) as s:
            p1, o = mergeBed(s, ["-d", "100"])
            gaps_bed = read_bed_default(o)

        with bed_to_stream(gaps_bed) as s:
            p2, o = complementBed(s, genome_path)
            gaps_with_parY = read_bed_default(o)

        # If we have a parY bed and chrY is included, subtract parY from the
        # gaps bed, otherwise just link them since we have nothing to subtract
        # off
        if hasattr(inputs, "parY") and bd.want_xy_y:
            parY_src = Path(inputs["parY"])
            with bed_to_stream(gaps_with_parY) as s:
                p3, o = subtractBed(s, parY_src, genome_path)
                gaps_no_parY = read_bed_default(o)
            check_processes([p1, p2, p3], log)

            write_bed(parY_out, gaps_with_parY)
            write_bed(auto_out, gaps_no_parY)
        else:
            check_processes([p1, p2], log)

            write_bed(auto_out, gaps_bed)
            parY_out.symlink_to(auto_out.resolve())


main(snakemake, snakemake.config)  # type: ignore
