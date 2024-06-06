from typing import Any
from pathlib import Path
import common.config as cfg
from common.bed import (
    write_bed,
    bed_to_stream,
    complementBed,
    mergeBed,
    subtractBed,
)
from common.io import check_processes, tee, bgzip_file
from common.functional import match1_unsafe, match2_unsafe, DesignError
import pandas as pd


# convert genome to bed file (where each region is just the length of one
# chromosome)
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
    genome_path = cfg.smk_to_input_name(smk, "genome")
    auto_out = cfg.smk_to_output_name(smk, "auto")
    parY_out = cfg.smk_to_output_name(smk, "parY")
    log = cfg.smk_to_log(smk)
    ws: dict[str, Any] = smk.wildcards

    def go(
        x: cfg.BuildData_[cfg.RefSrcT, cfg.AnyBedT, cfg.AnyVcfT, cfg.AnyBedTxtT]
    ) -> cfg.BedFile[cfg.AnyBedT] | cfg.AnyBedTxtT | None:
        return x.refdata.strat_inputs.gap

    # If we have gap input, make the gapless file, otherwise just symlink to the
    # genome bed file (which just means the entire genome is gapless)
    try:
        gap_inputs = cfg.smk_to_inputs_name(smk, "gaps")

        gaps_df = sconf.with_build_data_and_bed_full(
            cfg.wc_to_reffinalkey(ws),
            cfg.wc_to_buildkey(ws),
            go,
            lambda bd, bf: match1_unsafe(
                gap_inputs, lambda i: cfg.read_filter_sort_hap_bed(bd, bf, i)
            ),
            lambda bd, bf: match1_unsafe(
                gap_inputs, lambda i: cfg.read_filter_sort_dip1to1_bed(bd, bf, i)
            ),
            lambda hap, bd, bf: match1_unsafe(
                gap_inputs,
                lambda i: hap.choose(*cfg.read_filter_sort_dip1to2_bed(bd, bf, i)),
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
                    hap.choose(i0, i1),
                    hap,
                ),
            ),
        )

        with bed_to_stream(gaps_df) as s:
            p1, o1 = mergeBed(s, ["-d", "100"])
            p2, o2 = complementBed(o1, genome_path)

            # If we have a parY bed, subtract parY from the gaps bed, otherwise
            # just link them since we have nothing to subtract off
            try:
                parY_path = cfg.smk_to_input_name(smk, "parY")
                p3, o3, o4 = tee(o2)
                p4, o5 = subtractBed(o3, parY_path, genome_path)

                bgzip_file(o4, parY_out)
                bgzip_file(o5, auto_out)

                check_processes([p1, p2, p3, p4], log)
            except DesignError:
                bgzip_file(o2, auto_out)
                parY_out.symlink_to(auto_out.resolve())

                check_processes([p1, p2], log)

    except DesignError:
        genome_bed = read_genome_bed(genome_path)
        write_bed(auto_out, genome_bed)
        parY_out.symlink_to(auto_out.resolve())


main(snakemake, snakemake.config)  # type: ignore
