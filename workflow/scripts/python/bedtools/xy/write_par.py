from typing import Any
import subprocess as sp
import common.config as cfg


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, Any] = smk.wildcards
    out = cfg.smk_to_output(smk)
    i = cfg.ChrIndex.from_name(cfg.wc_lookup(ws, "sex_chr"))

    # TODO this pattern is DRY?
    rfk = cfg.wc_to_reffinalkey(ws)
    rk, hap = cfg.parse_full_refkey(rfk)
    pat = sconf.refkey_to_xy_ref_chr_pattern(rfk, i)
    cxy = sconf.to_ref_data(rk).strat_inputs.xy

    par_fun = cfg.choose_xy_unsafe(i, cxy.fmt_x_par_unsafe, cxy.fmt_y_par_unsafe)

    par = par_fun(pat)
    with open(out, "wb") as f:
        sp.run(["bgzip", "-c"], input=par, stdout=f, text=True)


main(snakemake, snakemake.config)  # type: ignore
