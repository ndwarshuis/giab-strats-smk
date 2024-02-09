import hashlib
from logging import Logger
from typing import IO
from pathlib import Path
from Bio import bgzf  # type: ignore
import subprocess as sp
from common.functional import not_none_unsafe, noop
import gzip


def get_md5(path: str, unzip: bool = False) -> str:
    h = hashlib.md5()
    do_unzip = path.endswith(".gz") and unzip is True
    with gzip.open(path, "rb") if do_unzip else open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            h.update(chunk)
    return h.hexdigest()


def setup_logging(path: str, console: bool = False) -> Logger:
    import logging

    logging.basicConfig(filename=path, level=logging.INFO)
    logging.captureWarnings(True)
    logger = logging.getLogger()
    if console:
        logger.addHandler(logging.StreamHandler())
    return logger


def is_gzip(p: Path) -> bool:
    # test if gzip by trying to read first byte
    with gzip.open(p, "r") as f:
        try:
            f.read(1)
            return True
        except gzip.BadGzipFile:
            return False


def is_bgzip(p: Path) -> bool:
    # since bgzip is in blocks (vs gzip), determine if in bgzip by
    # attempting to seek first block
    with open(p, "rb") as f:
        try:
            next(bgzf.BgzfBlocks(f), None)
            return True
        except ValueError:
            return False


def spawn_stream(
    cmd: list[str],
    i: IO[bytes] | int,
) -> tuple[sp.Popen[bytes], IO[bytes]]:
    p = sp.Popen(cmd, stdin=i, stdout=sp.PIPE)
    # ASSUME since we typed the inputs so that the stdin/stdout can only take
    # file descriptors or file streams, the return for each will never be
    # none
    return p, not_none_unsafe(p.stdout, noop)


def bgzip_file(i: IO[bytes], p: Path) -> sp.CompletedProcess[bytes]:
    with open(p, "wb") as f:
        return bgzip(i, f)


def bgzip(i: IO[bytes], o: IO[bytes]) -> sp.CompletedProcess[bytes]:
    """Stream bgzip to endpoint.

    NOTE: this will block since this is almost always going to be the
    final step in a pipeline.
    """
    return sp.run(["bgzip", "-c"], stdin=i, stdout=o)


def gunzip(i: Path) -> tuple[sp.Popen[bytes], IO[bytes]]:
    """Stream bgzip to endpoint.

    NOTE: this will block since this is almost always going to be the
    final step in a pipeline.
    """
    p = sp.Popen(["gunzip", "-c", i], stdout=sp.PIPE)
    return (p, not_none_unsafe(p.stdout, noop))
