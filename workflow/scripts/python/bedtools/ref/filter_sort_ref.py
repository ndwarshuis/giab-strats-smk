from typing import Any
from pathlib import Path
import subprocess as sp
from common.samtools import filter_sort_fasta
import common.config as cfg
from common.functional import DesignError

# Given a FASTA and a build, filter and sort the FASTA for the chromosomes
# we want, write it either uncompressed or compressed, then index the output
# fasta and write a genome file

Finished = tuple[bytes, bytes, int]


def parse_refkeys_config(p: Path) -> cfg.RefkeyConfiguration:
    try:
        return next(c for c in cfg.RefkeyConfiguration if p.name.startswith(c.value))
    except StopIteration:
        raise DesignError(f"could not parse refkeys config from path: {p}")


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    fa = Path(smk.input[0])

    fa_out = Path(smk.output["fasta"])
    genome_out = Path(smk.output["genome"])
    # NOTE: index output path not needed directly
    log = Path(smk.log[0])

    bgzip = fa_out.name.endswith(".gz")

    def go(cs: set[cfg.ChrIndex], pat: cfg.ChrPattern) -> Finished:
        return filter_sort_fasta(fa, fa_out, cs, pat, bgzip)

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
    _, err, rc = f(cfg.wc_to_reffinalkey(ws), cfg.wc_to_buildkey(ws))

    # If something goes wrong, write it down and bail
    with open(log, "wb") as lf:
        lf.write(err)

    if rc != 0:
        exit(1)

    # Write the index
    sp.run(["samtools", "faidx", str(fa_out)])

    # Write the genome. Note that we need to do this separately since I don't
    # know how to stream a .gzi file in a format that would make this easy. The
    # genome file is simply the first two columns of the .fai (uncompressed)
    # index.
    with open(genome_out, "w") as g:
        cmd = ["samtools", "faidx", str(fa_out), "-o", "-"]
        p = sp.run(cmd, stdout=sp.PIPE, text=True)
        for x in p.stdout.splitlines():
            xs = x.split("\t")
            g.write("\t".join([xs[0], "0", xs[1]]) + "\n")


main(snakemake, snakemake.config)  # type: ignore
