from typing import Any
from pathlib import Path
import subprocess as sp
from common.samtools import chr_filter_sort_fasta
import common.config as cfg
import common.io as io

# Given a FASTA and a build, filter and sort the FASTA for the chromosomes
# we want, write it either uncompressed or compressed, then index the output
# fasta and write a genome file

Finished = sp.CompletedProcess[bytes]


def parse_refkeys_config(p: Path) -> tuple[bool, bool]:
    return cfg.prefix_to_refkey_config(p.name)


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    fa = cfg.smk_to_input(smk)
    fa_out = cfg.smk_to_output_name(smk, "fa")
    index_out = cfg.smk_to_output_name(smk, "index")
    genome_out = cfg.smk_to_output_name(smk, "genome")
    log = cfg.smk_to_log(smk)

    rk = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)

    split, nohap = parse_refkeys_config(fa_out)

    cis = sconf.to_build_data(cfg.strip_full_refkey(rk), bk).build_chrs

    def go(pat: cfg.ChrPattern) -> Finished:
        return chr_filter_sort_fasta(fa, fa_out, log, pat.to_names(cis))

    def hap1(rd: cfg.HapRefData) -> Finished:
        return go(rd.ref.chr_pattern)

    def dip1(rd: cfg.Dip1RefData) -> Finished:
        return go(rd.ref.chr_pattern)

    def dip1_split(hap: cfg.Haplotype, rd: cfg.Dip1RefData) -> Finished:
        return go(rd.ref.chr_pattern.to_hap_pattern(hap))

    def dip2(hap: cfg.Haplotype, rd: cfg.Dip2RefData) -> Finished:
        return go(rd.ref.chr_pattern.choose(hap))

    proc_fa = sconf.with_ref_data_full_rconf(
        rk,
        split,
        nohap,
        hap1,
        dip1,
        dip1_split,
        dip2,
    )

    if proc_fa.returncode != 0:
        exit(1)

    # Write the index and genome. The genome file is simply the first two
    # columns of the .fai (uncompressed) index.
    with open(index_out, "w") as i, open(genome_out, "w") as g:
        cmd = ["samtools", "faidx", str(fa_out), "-o", "-"]
        proc_idx = sp.run(cmd, stdout=sp.PIPE)
        for x in proc_idx.stdout.splitlines():
            i.write(x.decode() + "\n")
            g.write("\t".join(x.decode().split("\t")[0:2]) + "\n")

    # TODO log errors here if the index step fails
    io.check_processes([proc_idx], log)


main(snakemake, snakemake.config)  # type: ignore
