from pathlib import Path
import subprocess as sp
from common.config import ChrIndex, ChrPattern


def _faidx_cmd(i: Path, chrs: list[str]) -> list[str]:
    return ["samtools", "faidx", str(i), *chrs]


def _chr_faidx_cmd(i: Path, cs: set[ChrIndex], pat: ChrPattern) -> list[str]:
    return _faidx_cmd(i, pat.to_names(cs))


def chr_filter_sort_fasta(
    i: Path,
    o: Path,
    log: Path,
    cs: set[ChrIndex],
    pat: ChrPattern,
) -> sp.CompletedProcess[bytes]:
    """Select and order chromosomes by index and stream to FASTA."""
    with open(o, "wb") as f, open(log, "wb") as lf:
        return sp.run(_chr_faidx_cmd(i, cs, pat), stdout=f, stderr=lf)


def filter_sort_fasta(
    i: Path,
    o: Path,
    log: Path,
    chrs: list[str],
) -> sp.CompletedProcess[bytes]:
    """Select and order chromosomes according to 'chrs' to FASTA.

    Used for mappability where chromosomes are not merely indexes (ie includes
    all the unmapped contigs) and GEM can't use bgzip'ed files.
    """
    with open(o, "wb") as f, open(log, "wb") as lf:
        return sp.run(_faidx_cmd(i, chrs), stdout=f, stderr=lf)
