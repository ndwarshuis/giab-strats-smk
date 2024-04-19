import gzip
import re
from pathlib import Path
from typing import Any
import common.config as cfg
from common.functional import not_none_unsafe
from common.bed import ChrName, InternalChrIndex


def get_chr_paired_mapper(
    sconf: cfg.GiabStrats,
    rk: cfg.RefKeyFullS,
    bk: cfg.BuildKey,
) -> dict[ChrName, tuple[cfg.ShortChrName, cfg.Haplotype, InternalChrIndex]]:
    cis = sconf.to_build_data(cfg.strip_full_refkey(rk), bk).build_chrs
    xs = sconf.with_ref_data_full(
        rk,
        lambda rd: rd.ref.chr_pattern.to_chr_data(cis, cfg.Haplotype.HAP1),
        lambda rd: rd.ref.chr_pattern.to_chr_data(cis),
        lambda hap, rd: rd.ref.chr_pattern.choose(hap).to_chr_data(cis, hap),
    )
    return {x.name: (x.shortname, x.haplotype, x.idx) for x in xs}


def get_distribution(
    ipath: Path,
    opath: Path,
    genome: dict[str, int],
    k: int,
) -> None:
    acc: dict[str, dict[int, int]] = {}

    with gzip.open(ipath, "rt") as f:
        for x in f:
            s = x.split("\t")
            chrom = s[0]
            start = int(s[1])
            end = int(s[2])
            cur = start
            if chrom not in acc:
                acc[chrom] = {}
            while cur < end:
                idx = floor(cur / k) * k
                if idx not in acc[chrom]:
                    acc[chrom][idx] = 0
                next_idx = idx + k
                acc[chrom][idx] += (next_idx if end >= next_idx else end) - cur
                cur = next_idx

    with open(opath, "w") as f:
        for chr, ks in acc.items():
            for i, c in ks.items():
                total = genome[chr]
                ys = [chr, str(i / total), str(c / k), str(k / total)]
                f.write("\t".join(ys) + "\n")


def sum_bed_file(path: Path) -> dict[str, int]:
    # ASSUME bed file has been tested for validity
    acc = {}
    with gzip.open(path, "r") as f:
        for i in f:
            cells = i.split(b"\t")
            chrom = cells[0].decode()
            if chrom not in acc:
                acc[chrom] = 0
            acc[chrom] += int(cells[2]) - int(cells[1])

    return acc


def bedpath_to_strat_names(bp: Path) -> tuple[str, str]:
    strat_name = not_none_unsafe(
        re.match("^[^_]+_(.+).bed.gz", bp.name),
        lambda x: x[1],
    )

    return (bp.parent.name, strat_name)


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    gapless_path = cfg.smk_to_input_name(smk, "gapless")
    bedlist = cfg.smk_to_input_name(smk, "bedlist")
    out = cfg.smk_to_output(smk)

    rk = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)

    chr_hap_mapper = get_chr_paired_mapper(sconf, rk, bk)
    gapless_sum = sum_bed_file(gapless_path)

    with open(bedlist, "r") as i, gzip.open(out, "w") as o:
        for bedpath in i:
            bp = Path(bedpath.strip())
            level_name, strat_name = bedpath_to_strat_names(bp)
            # TODO I could make these names a bit shorter by only appending
            # the haplotype to dip2 references since that's the only place
            # we need to do this
            for chrom, bedtotal in sum_bed_file(bp).items():
                shortChromName, hap, idx = chr_hap_mapper[ChrName(chrom)]
                fraction = bedtotal / gapless_sum[chrom]
                fullkey = f"{rk}@{bk}-{hap.value + 1}"
                newline = "\t".join(
                    [
                        fullkey,
                        level_name,
                        strat_name,
                        shortChromName,
                        str(idx - hap.value * 24),
                        str(fraction),
                    ]
                )
                o.write((newline + "\n").encode())


main(snakemake, snakemake.config)  # type: ignore