from typing import Any
import common.config as cfg
from common.stratdiff.lib.diff import compare_all
from common.io import setup_logging
from common.functional import DesignError


log = setup_logging(snakemake.log[0])  # type: ignore


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, Any] = smk.wildcards
    out = cfg.smk_to_output(smk)
    new_list_path = cfg.smk_to_input_name(smk, "new_list")
    old_path = cfg.smk_to_input_name(smk, "old")
    rkf = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)
    comparison = sconf.get_comparison(rkf, bk)

    fm = sconf.buildkey_to_ref_mappers(
        cfg.wc_to_reffinalkey(ws),
        cfg.wc_to_buildkey(ws),
    )[1]
    chr_names = [str(s) for s in fm.values()]

    if comparison is None:
        raise DesignError("comparison should not be None")

    outdir = out.parent

    logged = compare_all(
        new_list_path.parent,
        old_path,
        outdir,
        comparison.path_mapper,
        comparison.replacements,
        chr_names,
        comparison.ignore_generated,
        comparison.ignore_other,
    )
    for x in logged:
        log.info(x)


main(snakemake, snakemake.config)  # type: ignore
