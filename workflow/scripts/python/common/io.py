import hashlib
from logging import Logger
from typing import IO
from pathlib import Path
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


def is_gzip_stream(i: IO[bytes]) -> bool:
    # test if gzip by inspecting the first two bytes in the stream
    return i.read(2) == b"\x1f\x8b"


def is_bgzip_stream(i: IO[bytes]) -> bool:
    # since bgzip is in blocks (vs gzip), determine if in bgzip by
    # attempting to seek first block
    return i.read(4) == b"\x1f\x8b\x08\x04"


def is_gzip(p: Path) -> bool:
    # test if gzip by trying to read first byte
    with open(p, "rb") as f:
        return is_gzip_stream(f)


def is_bgzip(p: Path) -> bool:
    # since bgzip is in blocks (vs gzip), determine if in bgzip by
    # attempting to seek first block
    with open(p, "rb") as f:
        return is_bgzip_stream(f)


def spawn_stream(
    cmd: list[str],
    i: IO[bytes] | int | None = None,
) -> tuple[sp.Popen[bytes], IO[bytes]]:
    p = sp.Popen(cmd, stdin=i, stdout=sp.PIPE, stderr=sp.PIPE)
    # ASSUME since we typed the inputs so that the stdin/stdout can only take
    # file descriptors or file streams, the return for each will never be
    # none
    return p, not_none_unsafe(p.stdout, noop)


def bgzip_file(i: IO[bytes], p: Path) -> sp.CompletedProcess[bytes]:
    with open(p, "wb") as f:
        return bgzip(i, f)


def gzip_(i: IO[bytes], o: IO[bytes]) -> sp.CompletedProcess[bytes]:
    """Stream gzip to endpoint.

    NOTE: this will block since this is almost always going to be the
    final step in a pipeline.
    """
    return sp.run(["gzip", "-c"], stdin=i, stdout=o)


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


def gunzip_stream(i: IO[bytes]) -> tuple[sp.Popen[bytes], IO[bytes]]:
    p = sp.Popen(["gunzip", "-c"], stdin=i, stdout=sp.PIPE)
    return (p, not_none_unsafe(p.stdout, noop))


def check_processes(
    ps: list[sp.Popen[bytes] | sp.CompletedProcess[bytes]], log: Path
) -> None:
    some_error = False
    with open(log, "w") as lf:
        for p in ps:
            if isinstance(p, sp.CompletedProcess):
                err = p.stderr
            else:
                _, err = p.communicate()  # break deadlocks if there are any
            if p.returncode != 0:
                some_error = True

                args = p.args
                if isinstance(args, list):
                    cmd = " ".join(args)
                elif isinstance(args, bytes):
                    cmd = args.decode()
                else:
                    cmd = str(args)
                lf.write(f"{cmd}: {err.decode()}\n")

    if some_error:
        exit(1)
