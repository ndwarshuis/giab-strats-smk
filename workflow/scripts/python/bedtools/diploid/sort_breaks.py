from pathlib import Path
from typing import Any
import common.config as cfg
from common.bed import (
    complementBed,
    bgzip_file,
    read_bed_default,
    filter_sort_bed,
    bed_to_stream,
)

# a general bed-sorting script which will puke if the refkey is incorrect


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    bed_in = Path(smk.input["bed"][0])
    genome_in = Path(smk.input["genome"][0])

    m = sconf.with_build_data_split_full_nohap(
        cfg.wc_to_reffinalkey(ws),
        cfg.wc_to_buildkey(ws),
        lambda hap, bd: bd.refdata.ref.noop_conversion(bd.chr_indices).split(hap),
        lambda hap, bd: hap.from_either(
            *bd.refdata.ref.noop_conversion(bd.chr_indices)
        ),
    )
    im = m.init_mapper
    fm = m.final_mapper

    df = read_bed_default(bed_in)
    with bed_to_stream(filter_sort_bed(im, fm, df)) as s:
        _, o = complementBed(s, genome_in)
        bgzip_file(o, smk.output[0])


main(snakemake, snakemake.config)  # type: ignore
