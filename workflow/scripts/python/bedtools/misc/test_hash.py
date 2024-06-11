from typing import Any
from common.io import get_md5
import common.config as cfg


def main(smk: Any) -> None:
    i = cfg.smk_to_input(smk)
    o = cfg.smk_to_output(smk)
    v = cfg.smk_to_param_str(smk, "validhash")
    t = get_md5(i)
    if v == t:
        o.touch()
    else:
        exit(1)


main(snakemake)  # type: ignore
