from pathlib import Path
import subprocess as sp
from common.functional import not_none_unsafe
from common.config import ChrIndex, ChrPattern
from common.io import bgzip_file


def _faidx_cmd(i: Path, chrs: list[str]) -> list[str]:
    return ["samtools", "faidx", str(i), *chrs]


def _chr_faidx_cmd(i: Path, cs: set[ChrIndex], pat: ChrPattern) -> list[str]:
    return _faidx_cmd(i, pat.to_names(cs))


def filter_sort_fasta(
    i: Path,
    o: Path,
    cs: set[ChrIndex],
    pat: ChrPattern,
    bgzip: bool,
) -> tuple[bytes, bytes, int]:
    """Select and order chromosomes by index and stream to (un)compressed FASTA."""
    if bgzip:
        p = sp.Popen(_chr_faidx_cmd(i, cs, pat), stdout=sp.PIPE)
        not_none_unsafe(p.stdout, lambda h: bgzip_file(h, o))
        p.wait()
        return (*p.communicate(), p.returncode)
    else:
        c = sp.run(_chr_faidx_cmd(i, cs, pat))
        return c.stdout, c.stderr, c.returncode


def filter_sort_fasta_nocompress(
    i: Path,
    o: Path,
    chrs: list[str],
) -> sp.CompletedProcess[bytes]:
    """Select and order chromosomes according to 'chrs' to an uncompressed FASTA.

    Used for mappability where chromosomes are not merely indexes (ie includes
    all the unmapped contigs) and GEM can't use bgzip'ed files.
    """
    with open(o, "w") as f:
        return sp.run(_faidx_cmd(i, chrs), stdout=f)
