import gzip
import threading
import os
import contextlib
import pandas as pd
import subprocess as sp
from typing import NewType, IO, Generator
from pathlib import Path
from common.functional import not_none_unsafe, noop, DesignError
from common.io import spawn_stream, bgzip_file

# from Bio import bgzf  # type: ignore
import csv

# A complete chromosome name like "chr1" or "chr21_PATERNAL"
ChrName = NewType("ChrName", str)

# An integer representing the order of a chromosome during the filter/sort
# steps described below. Ranged 0-41 (all chromosomes for both haplotypes)
InternalChrIndex = NewType("InternalChrIndex", int)

# A mapping b/t chromosome names found in a bed-like file or fasta file and the
# internal sort order of each chromosome. Chromosomes that are not included here
# correspond to those that are excluded from the build.
InitMapper = dict[ChrName, InternalChrIndex]

# A mapping b/t internal chromosome sort order and the chromosome names on the
# target reference.
FinalMapper = dict[InternalChrIndex, ChrName]

# A mapping b/t chromosome names and haplotype; used for splitting inputs that
# contain both haplotypes into two files with one haplotype each. True =
# paternal, False = maternal
SplitMapper = dict[str, bool]

BedColumns = tuple[int, int, int]


def make_split_mapper(im: InitMapper, fm: FinalMapper) -> SplitMapper:
    return {n: i in fm for n, i in im.items()}


def read_bed(
    path: Path,
    columns: BedColumns,
    skip_lines: int,
    sep: str,
    more: list[int],
) -> pd.DataFrame:
    """Read a bed file as a pandas dataframe.

    Return a dataframe where the first three columns are numbered 0, 1, 2 and
    typed str, int, int (first is str regardless of how the chr names are
    formated). Columns from 'more' are appended to the end of the dataframe
    in the order given starting from 3.
    """
    # NOTE: the 'comment="#"' parameter in pandas will strip everything after
    # the '#' in the line, which if at the beginning will include the entire
    # line and it will be skipped, and if not will only obliterate the
    # remainder. Either way this is a problem, since I only care about initial
    # lines starting with '#'. Bed files shouldn't have comments in the middle,
    # and some 'bed' files use '#' as a delimiter within a field.
    #
    # This hacky bit will count the number of lines starting with a '#' and add
    # to the original "skip_lines" parameter, thereby skipping all starting
    # comments as well as the number of lines we wish to skip with 'skip_lines'.
    total_skip = skip_lines
    with gzip.open(path, "rt") as f:
        while line := next(f, None):
            if line.startswith("#"):
                total_skip += 1
            else:
                break
    return read_bed_raw(path, columns, total_skip, sep, more)


def read_bed_raw(
    h: IO[bytes] | Path,
    columns: BedColumns,
    skip_lines: int,
    sep: str,
    more: list[int],
) -> pd.DataFrame:
    bedcols = [*columns, *more]
    df = pd.read_table(
        h,
        header=None,
        usecols=bedcols,
        sep=sep,
        skiprows=skip_lines,
        # satisfy type checker :/
        dtype={
            **{columns[0]: str, columns[1]: int, columns[2]: int},
            **{m: str for m in more},
        },
    )
    return df.set_axis(range(len(bedcols)), axis=1)


def read_bed_default(h: IO[bytes] | Path) -> pd.DataFrame:
    return read_bed_raw(h, (0, 1, 2), 0, "\t", [])


def bed_to_text(df: pd.DataFrame) -> Generator[str, None, None]:
    for r in df.itertuples(index=False):
        yield "\t".join(r)


def write_bed_stream(h: IO[str], df: pd.DataFrame) -> None:
    """Stream bed to handle from a dataframe.

    Dataframe is not checked to make sure it is a "real" bed file.
    """
    w = csv.writer(h, delimiter="\t")
    for r in df.itertuples(index=False):
        w.writerow(r)


@contextlib.contextmanager
def bed_to_stream(df: pd.DataFrame) -> Generator[IO[bytes], None, None]:
    """Transform bed-like dataframe into a bytestream.

    This is useful for using a pandas dataframe in a bedtools subprocess.
    """
    r, w = os.pipe()
    _r = os.fdopen(r, "rb")
    _w = os.fdopen(w, "w")

    # NOTE: python threads can't pass exceptions between themselves and the
    # calling thread, so need to hack something to signal when something bad
    # happens. Do this by closing read side of the thread (which will also make
    # any process that depends on the pipe die since it can't be opened after
    # we close it), and then testing if the read end is closed after the context
    # block.
    def read_df() -> None:
        try:
            write_bed_stream(_w, df)
        except Exception:
            _r.close()
        finally:
            _w.close()

    t = threading.Thread(target=read_df)
    t.start()

    try:
        yield _r
        if _r.closed:
            raise DesignError("bed stream thread exited unexpectedly")
    finally:
        _r.close()
        _w.close()


def write_bed(path: Path, df: pd.DataFrame) -> None:
    """Write a bed file in bgzip format from a dataframe.

    Dataframe is not checked to make sure it is a "real" bed file.
    """
    with bed_to_stream(df) as s:
        bgzip_file(s, path)
    # with bgzf.open(path, "w") as f:
    #     write_bed_stream(f, df)


def complementBed(
    i: IO[bytes] | int,
    genome: Path,
) -> tuple[sp.Popen[bytes], IO[bytes]]:
    cmd = ["complementBed", "-i", "stdin", "-g", str(genome)]
    return spawn_stream(cmd, i)


def subtractBed(
    i: IO[bytes] | int,
    b: Path,
    genome: Path,
) -> tuple[sp.Popen[bytes], IO[bytes]]:
    cmd = ["subtractBed", "-a", "stdin", "-b", str(b), "-sorted", "-g", str(genome)]
    return spawn_stream(cmd, i)


def intersectBed(
    i: IO[bytes] | int,
    b: Path,
    genome: Path,
) -> tuple[sp.Popen[bytes], IO[bytes]]:
    cmd = ["intersectBed", "-a", "stdin", "-b", str(b), "-sorted", "-g", str(genome)]
    return spawn_stream(cmd, i)


def mergeBed(
    i: IO[bytes] | int,
    opts: list[str],
) -> tuple[sp.Popen[bytes], IO[bytes]]:
    return spawn_stream(["mergeBed", "-i", "stdin", *opts], i)


def multiIntersectBed(ps: list[Path]) -> tuple[sp.Popen[bytes], IO[bytes]]:
    cmd = ["multiIntersectBed", "-i", *[str(p) for p in ps]]
    p = sp.Popen(cmd, stdout=sp.PIPE)
    return p, not_none_unsafe(p.stdout, noop)


def sort_bed_numerically(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """Sort a bed file encoded by a dataframe.

    Assumes the first three columns correspond to coordinates, and that all are
    integer typed. Use 'n = 2' to sort only by chr/start, and 'n=1' to sort only
    by chr.

    """
    cols = df.columns.tolist()
    bycols = [cols[i] for i in range(0, n)]
    return df.sort_values(
        by=bycols,
        axis=0,
        ignore_index=True,
    )


def filter_sort_bed(
    from_map: InitMapper,
    to_map: FinalMapper,
    df: pd.DataFrame,
    n: int = 3,
) -> pd.DataFrame:
    """Filter and sort a bed file.

    Arguments:
    from_map - dict containing chr name -> int mappings (int = order)
    to_map - dict containing int -> chr name mappings
    df - dataframe to sort

    Assumes the first three columns correspond to the coordinates of a bed
    file.

    Any chr name not specified in 'from_map' will be removed (hence the filter).
    Furthermore, 'to_map' should contain at least all corresponding entries
    from 'from_map', otherwise the final df will have NaNs.
    """
    chr_col = df.columns.tolist()[0]
    df[chr_col] = df[chr_col].map(from_map)
    df = sort_bed_numerically(df.dropna(subset=[chr_col]), n)
    df[chr_col] = df[chr_col].map(to_map)
    return df


def split_bed(
    split_map: SplitMapper,
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split bed-like dataframe according to 'split_map' into two dataframes.

    Assume the first file is paternal and the second is maternal
    """
    chr_col = df.columns.tolist()[0]
    sp = df[chr_col].map(split_map)
    return df[sp], df[~sp]
