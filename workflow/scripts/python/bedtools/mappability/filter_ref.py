import re
from typing import Any
import common.config as cfg
from pathlib import Path
from common.samtools import filter_sort_fasta_nocompress

# This is the one exception to the pattern wherein we filter all input data to
# the chromosomes we actually care about. In the case of mappability, we want to
# run GEM against all random/unplaced contigs, since these are hard to map to
# begin with (otherwise they wouldn't be unplaced). Here, we need to filter
# the fasta to include all the primary contigs we care about + all the other
# unplaced/random contigs. Later we will filter out all the non-primary contigs.


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    fa = Path(smk.input["fa"])
    idx = Path(smk.input["idx"])
    out = Path(smk.output[0])

    # Parse the index file to get a list of all chromosomes in the FASTA
    # (including all the weird unmapped contigs we don't normally care about).
    # This obviously assumes this index was taken before the FASTA was filtered
    # and sorted. It shouldn't matter that these are not sorted since later
    # in the mappability pipeline we will sort again anyways since wig2bed
    # doesn't guarantee any particular order
    with open(idx, "rt") as f:
        all_chrs = [x.split("\t")[0] for x in f]

    def run_samtools(
        bd: cfg.BuildData_[cfg.RefSourceT, cfg.AnyBedT, cfg.AnyBedT_, cfg.AnySrcT],
        pat: cfg.ChrPattern,
    ) -> None:
        def any_map_pattern(c: str) -> bool:
            return any(
                [re.match(p, c) is not None for p in bd.refdata.mappability_patterns]
            )

        main_chrs = pat.to_names(bd.chr_indices)
        chrs = [c for c in all_chrs if c in main_chrs or any_map_pattern(c)]
        filter_sort_fasta_nocompress(fa, out, chrs)

    sconf.with_build_data_full(
        cfg.wc_to_reffinalkey(ws),
        cfg.wc_to_buildkey(ws),
        lambda bd: run_samtools(bd, bd.refdata.ref.chr_pattern),
        lambda bd: run_samtools(bd, bd.refdata.ref.chr_pattern),
        lambda hap, bd: run_samtools(bd, bd.refdata.ref.chr_pattern.from_either(hap)),
    )


main(snakemake, snakemake.config)  # type: ignore
