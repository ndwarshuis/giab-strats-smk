from pathlib import Path
from typing_extensions import assert_never
from typing import Any
import common.io as io
from common.functional import DesignError
import common.config as cfg


def main(smk: Any) -> None:
    opath = cfg.smk_to_output(smk)
    log = cfg.smk_to_log(smk)
    src: cfg.RefSrc | None = smk.params.src
    if isinstance(src, cfg.FileSrc_):
        # ASSUME these are already tested via the pydantic class for the
        # proper file format
        Path(opath).symlink_to(Path(src.filepath).resolve())

    elif isinstance(src, cfg.HttpSrc_):
        with open(opath, "wb") as f:
            # ensure all refs are in bgzip format so we can index them
            if io.curl_test(src.url, io.is_bgzip_stream, log):
                io.curl(src.url, f, log)
            elif io.curl_test(src.url, io.is_gzip_stream, log):
                p1, o1 = io.spawn_stream(io.curl_cmd(src.url))
                p2, o2 = io.gunzip_stream(o1)
                p3 = io.bgzip(o2, f)
                o1.close()
                o2.close()
                io.check_processes([p1, p2, p3], log)
            else:
                io.curl_gzip(src.url, f, log, True)

    elif src is None:
        raise DesignError("file src is null; this should not happen")
    else:
        assert_never(src)

    if src.md5 is not None and src.md5 != (actual := io.get_md5(opath, True)):
        with open(log, "a") as f:
            f.write(f"md5s don't match; wanted {src.md5}, actual {actual}\n")
        exit(1)


main(snakemake)  # type: ignore
