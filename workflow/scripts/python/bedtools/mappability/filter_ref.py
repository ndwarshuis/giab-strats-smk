import re
from typing import Any
import common.config as cfg
from common.samtools import filter_sort_fasta
from common.bed import ChrName

# This is the one exception to the pattern wherein we filter all input data to
# the chromosomes we actually care about. In the case of mappability, we want to
# run GEM against all random/unplaced contigs, since these are hard to map to
# begin with (otherwise they wouldn't be unplaced). Here, we need to filter
# the fasta to include all the primary contigs we care about + all the other
# unplaced/random contigs. Later we will filter out all the non-primary contigs.


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    fa = cfg.smk_to_input_name(smk, "fa")
    idx = cfg.smk_to_input_name(smk, "idx")
    log = cfg.smk_to_log(smk)
    out = cfg.smk_to_output(smk)

    # Parse the index file to get a list of all chromosomes in the FASTA
    # (including all the weird unmapped contigs we don't normally care about).
    # This obviously assumes this index was taken before the FASTA was filtered
    # and sorted. It shouldn't matter that these are not sorted since later
    # in the mappability pipeline we will sort again anyways since wig2bed
    # doesn't guarantee any particular order
    with open(idx, "rt") as f:
        all_chrs = [x.split("\t")[0] for x in f]

    def run_samtools(
        bd: cfg.BuildData_[cfg.RefSrcT, cfg.AnyBedT, cfg.AnyVcfT, cfg.AnyBedTxtT],
        pat: cfg.ChrPattern,
    ) -> None:
        def any_map_pattern(c: str) -> bool:
            return any(
                [re.match(p, c) is not None for p in bd.refdata.mappability_patterns]
            )

        main_chrs = pat.to_names(bd.build_chrs)
        chrs = cfg.OrderedHapChrNames(
            [ChrName(c) for c in all_chrs if c in main_chrs or any_map_pattern(c)]
        )
        filter_sort_fasta(fa, out, log, chrs)

    sconf.with_build_data_split_full(
        cfg.wc_to_reffinalkey(ws),
        cfg.wc_to_buildkey(ws),
        lambda bd: run_samtools(bd, bd.refdata.ref.chr_pattern),
        lambda hap, bd: run_samtools(
            bd, bd.refdata.ref.chr_pattern.to_hap_pattern(hap)
        ),
        lambda hap, bd: run_samtools(bd, bd.refdata.ref.chr_pattern.choose(hap)),
    )


main(snakemake, snakemake.config)  # type: ignore
