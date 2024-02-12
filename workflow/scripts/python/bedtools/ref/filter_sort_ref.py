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

    rk = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)

    cis = sconf.to_build_data(cfg.strip_full_refkey(rk), bk).build_chrs

    def go(pat: cfg.ChrPattern) -> Finished:
        return chr_filter_sort_fasta(fa, fa_out, log, pat.to_names(cis))

    def hap1(rd: cfg.HapRefData) -> Finished:
        return go(rd.ref.chr_pattern)

    def dip1(rd: cfg.Dip1RefData) -> Finished:
        return go(rd.ref.chr_pattern)

    def _dip1_split(hap: cfg.Haplotype, rd: cfg.Dip1RefData) -> Finished:
        return go(rd.ref.chr_pattern.to_hap_pattern(hap))

    def dip2(hap: cfg.Haplotype, rd: cfg.Dip2RefData) -> Finished:
        return go(rd.ref.chr_pattern.choose(hap))

    def standard(rk: cfg.RefKeyFullS) -> Finished:
        return sconf.with_ref_data_full(rk, hap1, dip1, dip2)

    def dip1_split(rk: cfg.RefKeyFullS) -> Finished:
        return sconf.with_ref_data_split_full(rk, hap1, _dip1_split, dip2)

    def dip1_split_nohap(rk: cfg.RefKeyFullS) -> Finished:
        return sconf.with_ref_data_split_full_nohap(rk, _dip1_split, dip2)

    f = cfg.choose_refkey_configuration(
        parse_refkeys_config(fa_out),
        standard,
        dip1_split,
        dip1_split_nohap,
    )

    proc_fa = f(rk)

    if proc_fa.returncode != 0:
        exit(1)

    # Write the index and genome. The genome file is simply the first two
    # columns of the .fai (uncompressed) index.
    with open(index_out, "w") as i, open(genome_out, "w") as g:
        cmd = ["samtools", "faidx", str(fa_out), "-o", "-"]
        proc_idx = sp.run(cmd, stdout=sp.PIPE, text=True)
        for x in proc_idx.stdout.splitlines():
            i.write(x + "\n")
            g.write("\t".join(x.split("\t")[0:2]) + "\n")

    # TODO log errors here if the index step fails
    if proc_idx.returncode != 0:
        exit(1)


main(snakemake, snakemake.config)  # type: ignore
