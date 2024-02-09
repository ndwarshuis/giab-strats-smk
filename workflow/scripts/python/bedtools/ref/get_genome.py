from typing import Any
import pandas as pd
import common.config as cfg
from common.bed import filter_sort_bed
from common.functional import DesignError


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, Any] = smk.wildcards

    allowed_refkeys = str(smk.params["allowed_refkeys"])

    if allowed_refkeys == "any":
        f = sconf.buildkey_to_ref_mappers
    elif allowed_refkeys == "split":
        f = sconf.buildkey_to_ref_mappers_split
    elif allowed_refkeys == "split_nohap":
        f = sconf.buildkey_to_ref_mappers_split_nohap
    else:
        raise DesignError("wrong 'allowed_keys'")

    im, fm = f(cfg.wc_to_reffinalkey(ws), cfg.wc_to_buildkey(ws))

    # ASSUME the input for this is a .fa.fai file (columns = chr, length)
    df = pd.read_table(
        smk.input[0],
        header=None,
        dtype={0: str, 1: int},
        usecols=[0, 1],
    )
    filtered = filter_sort_bed(im, fm, df, 2)
    filtered.to_csv(smk.output[0], sep="\t", header=False, index=False)


main(snakemake, snakemake.config)  # type: ignore
