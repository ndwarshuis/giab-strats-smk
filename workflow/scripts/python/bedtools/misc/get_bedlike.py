from pathlib import Path
from typing import Any
from typing_extensions import assert_never
import common.io as io
from common.functional import DesignError
import common.config as cfg


def main(smk: Any) -> None:
    bbBin = cfg.smk_to_input(smk)
    opath = cfg.smk_to_output(smk)
    log = cfg.smk_to_log(smk)
    src: cfg.AnyBedSrc | None = smk.params.src

    def check_md5(src: cfg.BedFileSrc) -> None:
        if src.md5 is not None and src.md5 != (actual := io.get_md5(opath, True)):
            with open(log, "a") as f:
                f.write(f"md5s don't match; wanted {src.md5}, actual {actual}\n")
            exit(1)

    if isinstance(src, cfg.FileSrc_):
        # ASSUME these are already tested via the pydantic class for the
        # proper file format
        opath.symlink_to(Path(src.filepath).resolve())
        check_md5(src)

    elif isinstance(src, cfg.HttpSrc_):
        with open(opath, "wb") as f:

            # Ensure all bedlike files are in gzip format; if we are getting a
            # bigbed file then force it to a bed then gzip it
            if src.url.endswith(".bb"):
                # ASSUME that any url that points to a big-bed will end in .bb
                #
                # TODO this may or may not actually be true, in which case
                # it might be sensible to override with a config switch
                p1, o1 = io.spawn_stream([str(bbBin), src.url, "stdout"])
                p2 = io.gzip_(o1, f)
                o1.close()
                io.check_processes([p1, p2], log)
            elif io.curl_test(src.url, io.is_gzip_stream, log):
                io.curl(src.url, f, log)
            else:
                io.curl_gzip(src.url, f, log, False)

        check_md5(src)

    elif isinstance(src, list):
        pass
        # raise DesignError("file src is null; this should not happen")

    elif src is None:
        raise DesignError("file src is null; this should not happen")

    else:
        assert_never(src)


main(snakemake)  # type: ignore
