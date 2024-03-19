import pandas as pd
import json
from pathlib import Path
from typing import Any, Callable
import common.config as cfg
import common.bed as write_bed
from common.functional import not_none_unsafe, match1_unsafe, noop


def filter_ct(df: pd.DataFrame) -> pd.DataFrame:
    return df[~df[3].str.startswith("ct_")]


CVOut = tuple[list[Path], list[Path]]


def filter_cds(fps: cfg.FunctionalParams, df: pd.DataFrame) -> pd.DataFrame:
    return df


def filter_vdj(fps: cfg.FunctionalParams, df: pd.DataFrame) -> pd.DataFrame:
    return df


def main(smk: Any) -> None:
    sconf: cfg.GiabStrats = smk.config
    ws: dict[str, str] = smk.wildcards
    ps: dict[str, str] = smk.params
    inputs: list[Path] = [Path(i) for i in smk.input]
    cds_pattern: str = ps["cds_output"]
    vdj_pattern: str = ps["vdj_output"]

    def to_output_pattern(rk: cfg.RefKeyFull) -> tuple[Path, Path]:
        return (
            cfg.sub_output_path(cds_pattern, rk),
            cfg.sub_output_path(vdj_pattern, rk),
        )

    def hap(bd: cfg.HapBuildData) -> CVOut:
        rd = bd.refdata
        rk = rd.ref.src.key(rd.refkey)
        cpath, vpath = to_output_pattern(rk)
        i = match1_unsafe(inputs, noop)
        bf = not_none_unsafe(cfg.bd_to_functional(bd), noop)
        df = cfg.read_filter_sort_hap_bed(bd, bf, i)
        df_ = filter_cds(bf.fparams, df)
        # df = cfg.read_filter_sort_hap_bed(bd, bf, i)
        # write_bed(cpath, filter_cds(bf.fparams, df))
        # if bd.build.include.vdj:
        #     write_bed(vpath, filter_vdj(bf.fparams, df))

    def dip1(i: Path, bd: cfg.Dip1BuildData) -> CVOut:
        pass

    def dip2(i: Path, bd: cfg.Dip2BuildData) -> CVOut:
        pass

    def write_output1(o: Path, ps: list[Path]) -> None:
        with open(o, "w") as f:
            json.dump([str(p) for p in ps], f)

    # cds, vdj = sconf.with_build_data_and_bed_i(
    #     cfg.wc_to_refkey(ws),
    #     cfg.wc_to_buildkey(ws),
    #     inputs,
    #     lambda bd: cfg.bd_to_si(cfg.si_to_functional, bd),
    #     hap,
    #     dip1to1,
    #     dip1to2,
    #     dip2to1,
    #     lambda i, bf, bd: cfg.wrap_dip_2to2_i_f(dip2to2, i, bf, bd),
    # )

    cs, vs = sconf.with_build_data(
        cfg.wc_to_refkey(ws),
        cfg.wc_to_buildkey(ws),
        hap,
        dip1,
        dip2,
    )

    write_output1(smk.output["cds"], cs)
    write_output1(smk.output["vdj"], vs)


main(snakemake)  # type: ignore
