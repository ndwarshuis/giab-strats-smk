from pathlib import Path
import subprocess as sp
from common.config import OrderedHapChrNames


def _faidx_cmd(i: Path, chrs: OrderedHapChrNames) -> list[str]:
    return ["samtools", "faidx", str(i), *chrs]


def chr_filter_sort_fasta(
    i: Path,
    o: Path,
    log: Path,
    cs: OrderedHapChrNames,
) -> sp.CompletedProcess[bytes]:
    """Select and order chromosomes by index and stream to FASTA."""
    with open(o, "wb") as f, open(log, "wb") as lf:
        return sp.run(_faidx_cmd(i, cs), stdout=f, stderr=lf)


def filter_sort_fasta(
    i: Path,
    o: Path,
    log: Path,
    chrs: OrderedHapChrNames,
) -> sp.CompletedProcess[bytes]:
    """Select and order chromosomes according to 'chrs' to FASTA.

    Used for mappability where chromosomes are not merely indexes (ie includes
    all the unmapped contigs) and GEM can't use bgzip'ed files.
    """
    with open(o, "wb") as f, open(log, "wb") as lf:
        return sp.run(_faidx_cmd(i, chrs), stdout=f, stderr=lf)
