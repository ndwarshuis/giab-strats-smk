from pathlib import Path
import subprocess as sp
from typing import Callable, IO
from typing_extensions import assert_never
import common.io as io
from common.functional import DesignError
import common.config as cfg

# hacky curl/gzip wrapper; this exists because I got tired of writing
# specialized rules to convert gzip/nozip files to bgzip and back :/
# Solution: force bgzip for references and gzip for bed

smk_log = snakemake.log[0]  # type: ignore

log = io.setup_logging(smk_log)

CURL = ["curl", "-f", "-Ss", "-L", "-q"]


def main(bbBin: str, opath: str, src: cfg.BedSrc | cfg.RefSrc | None) -> None:
    if isinstance(src, cfg.FileSrc_):
        # ASSUME these are already tested via the pydantic class for the
        # proper file format
        Path(opath).symlink_to(Path(src.filepath).resolve())

    elif isinstance(src, cfg.HttpSrc_):
        with open(opath, "wb") as f:
            curlcmd = [*CURL, src.url]
            bbcmd = [bbBin, src.url, "stdout"]

            # to test the format of downloaded files, sample the first 65000
            # bytes (which should be enough to get one block of a bgzip file,
            # which will allow us to test for it)
            curltestcmd = [*CURL, "-r", "0-65000", src.url]

            def curl() -> None:
                p = sp.run(curlcmd, stdout=f, stderr=sp.PIPE)
                io.check_processes([p], smk_log)

            def curl_test(testfun: Callable[[IO[bytes]], bool]) -> bool:
                p, o = io.spawn_stream(curltestcmd)
                res = testfun(o)
                io.check_processes([p], smk_log)
                return res

            def curl_gzip(use_bgzip: bool) -> None:
                zipf = io.bgzip if use_bgzip else io.gzip_
                p1, o1 = io.spawn_stream(curlcmd)
                p2 = zipf(o1, f)
                o1.close()
                io.check_processes([p1, p2], smk_log)

            def bb_bed_gzip() -> None:
                p1, o1 = io.spawn_stream(bbcmd)
                p2 = io.gzip_(o1, f)
                o1.close()
                io.check_processes([p1, p2], smk_log)

            # If we are getting a bed(like) file ensure it is in gzip
            # format; if we are getting a bigbed file then force it to a bed
            # then gzip it
            if isinstance(src, cfg.BedHttpSrc):
                # ASSUME that any url that points to a big-bed will end in .bb
                # TODO this may or may not actually be true, in which case
                # it might be sensible to override with a config switch
                if src.url.endswith(".bb"):
                    bb_bed_gzip()
                elif curl_test(io.is_gzip_stream):
                    curl()
                else:
                    curl_gzip(False)

            # if we are getting a fasta file ensure it is in bgzip format
            elif isinstance(src, cfg.RefHttpSrc):
                if curl_test(io.is_bgzip_stream):
                    curl()
                elif curl_test(io.is_gzip_stream):
                    p1, o1 = io.spawn_stream(curlcmd)
                    p2, o2 = io.gunzip_stream(o1)
                    p3 = io.bgzip(o2, f)
                    o1.close()
                    o2.close()
                    io.check_processes([p1, p2, p3], smk_log)
                else:
                    curl_gzip(True)
            else:
                assert_never(src)

    elif src is None:
        raise DesignError("file src is null; this should not happen")
    else:
        assert_never(src)

    if src.md5 is not None and src.md5 != (actual := io.get_md5(opath, True)):
        log.error("md5s don't match; wanted %s, actual %s", src.md5, actual)
        exit(1)


main(snakemake.input[0], snakemake.output[0], snakemake.params.src)  # type: ignore
