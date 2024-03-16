from typing import Any
import common.config as cfg
from common.bed import ChrName
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
    with open(smk.input[0], "r") as i, open(smk.output[0], "w") as o:
        for line in i:
            s = line.split("\t")
            chr = ChrName(s[0])
            try:
                chrIdx = im[chr]
            except KeyError:
                continue
            o.write(f"{fm[chrIdx]}\t{s[1]}")


main(snakemake, snakemake.config)  # type: ignore
