from typing import Any
from pathlib import Path
import subprocess as sp
from common.samtools import chr_filter_sort_fasta
import common.config as cfg
from common.functional import DesignError

# Given a FASTA and a build, filter and sort the FASTA for the chromosomes
# we want, write it either uncompressed or compressed, then index the output
# fasta and write a genome file

Finished = sp.CompletedProcess[bytes]


def parse_refkeys_config(p: Path) -> cfg.RefkeyConfiguration:
    try:
        return next(c for c in cfg.RefkeyConfiguration if p.name.startswith(c.value))
    except StopIteration:
        raise DesignError(f"could not parse refkeys config from path: {p}")


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    fa = Path(smk.input[0])

    fa_out = Path(smk.output["fa"])
    index_out = Path(smk.output["index"])
    genome_out = Path(smk.output["genome"])
    log = Path(smk.log[0])

    def go(cs: set[cfg.ChrIndex], pat: cfg.ChrPattern) -> Finished:
        return chr_filter_sort_fasta(fa, fa_out, log, cs, pat)

    def hap1(bd: cfg.HapBuildData) -> Finished:
        return go(bd.chr_indices, bd.refdata.ref.chr_pattern)

    def dip1(bd: cfg.Dip1BuildData) -> Finished:
        return go(bd.chr_indices, bd.refdata.ref.chr_pattern)

    def _dip1_split(hap: cfg.Haplotype, bd: cfg.Dip1BuildData) -> Finished:
        pat = bd.refdata.ref.chr_pattern.to_hap_pattern(hap)
        return go(bd.chr_indices, pat)

    def dip2(hap: cfg.Haplotype, bd: cfg.Dip2BuildData) -> Finished:
        pat = bd.refdata.ref.chr_pattern.from_either(hap)
        return go(bd.chr_indices, pat)

    def standard(rk: cfg.RefKeyFullS, bk: cfg.BuildKey) -> Finished:
        return sconf.with_build_data_full(rk, bk, hap1, dip1, dip2)

    def dip1_split(rk: cfg.RefKeyFullS, bk: cfg.BuildKey) -> Finished:
        return sconf.with_build_data_split_full(rk, bk, hap1, _dip1_split, dip2)

    def dip1_split_nohap(rk: cfg.RefKeyFullS, bk: cfg.BuildKey) -> Finished:
        return sconf.with_build_data_split_full_nohap(rk, bk, _dip1_split, dip2)

    f = cfg.choose_refkey_configuration(
        parse_refkeys_config(fa_out),
        standard,
        dip1_split,
        dip1_split_nohap,
    )

    proc_fa = f(cfg.wc_to_reffinalkey(ws), cfg.wc_to_buildkey(ws))

    # Write the index and genome. The genome file is simply the first two
    # columns of the .fai (uncompressed) index.
    with open(index_out, "w") as i, open(genome_out, "w") as g:
        cmd = ["samtools", "faidx", str(fa_out), "-o", "-"]
        proc_idx = sp.run(cmd, stdout=sp.PIPE, text=True)
        for x in proc_idx.stdout.splitlines():
            i.write(x + "\n")
            g.write("\t".join(x.split("\t")[0:2]) + "\n")

    # TODO log errors here if the index step fails
    if proc_fa.returncode != 0 or proc_idx.returncode != 0:
        exit(1)


main(snakemake, snakemake.config)  # type: ignore
