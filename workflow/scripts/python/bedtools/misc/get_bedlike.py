from pathlib import Path
from typing_extensions import assert_never
import common.io as io
from common.functional import DesignError
import common.config as cfg

smk_log = snakemake.log[0]  # type: ignore

log = io.setup_logging(smk_log)


def main(bbBin: str, opath: str, src: cfg.BedSrc | None) -> None:
    if isinstance(src, cfg.FileSrc_):
        # ASSUME these are already tested via the pydantic class for the
        # proper file format
        Path(opath).symlink_to(Path(src.filepath).resolve())

    elif isinstance(src, cfg.HttpSrc_):
        with open(opath, "wb") as f:

            # Ensure all bedlike files are in gzip format; if we are getting a
            # bigbed file then force it to a bed then gzip it
            if src.url.endswith(".bb"):
                # ASSUME that any url that points to a big-bed will end in .bb
                #
                # TODO this may or may not actually be true, in which case
                # it might be sensible to override with a config switch
                p1, o1 = io.spawn_stream([bbBin, src.url, "stdout"])
                p2 = io.gzip_(o1, f)
                o1.close()
                io.check_processes([p1, p2], smk_log)
            elif io.curl_test(src.url, io.is_gzip_stream, smk_log):
                io.curl(src.url, f, smk_log)
            else:
                io.curl_gzip(src.url, f, smk_log, False)

    elif src is None:
        raise DesignError("file src is null; this should not happen")
    else:
        assert_never(src)

    if src.md5 is not None and src.md5 != (actual := io.get_md5(opath, True)):
        log.error("md5s don't match; wanted %s, actual %s", src.md5, actual)
        exit(1)


main(snakemake.input[0], snakemake.output[0], snakemake.params.src)  # type: ignore
