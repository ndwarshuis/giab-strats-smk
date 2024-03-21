from pathlib import Path
from typing import Any
import common.config as cfg
from common.bed import (
    complementBed,
    read_bed_default,
    filter_sort_bed,
    bed_to_stream,
)
from common.io import bgzip_file, check_processes

# a general bed-sorting script which will puke if the refkey is incorrect


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    bed_in = cfg.smk_to_input_name(smk, "bed")
    genome_in = cfg.smk_to_input_name(smk, "genome")
    log = cfg.smk_to_log(smk)
    out = cfg.smk_to_output(smk)

    m = sconf.with_build_data_split_full_nohap(
        cfg.wc_to_reffinalkey(ws),
        cfg.wc_to_buildkey(ws),
        lambda hap, bd: bd.refdata.ref.noop_conversion(bd.build_chrs).split(hap),
        lambda hap, bd: hap.choose(*bd.refdata.ref.noop_conversion(bd.build_chrs)),
    )
    im = m.init_mapper
    fm = m.final_mapper

    df = read_bed_default(bed_in)
    with bed_to_stream(filter_sort_bed(im, fm, df)) as s:
        p1, o = complementBed(s, genome_in)
        bgzip_file(o, out)
        check_processes([p1], log)


main(snakemake, snakemake.config)  # type: ignore
