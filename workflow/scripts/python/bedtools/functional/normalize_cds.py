import pandas as pd
import json
from pathlib import Path
from typing import Any, assert_never
import common.config as cfg
from common.bed import write_bed
from common.functional import both, DesignError, fmap_maybe

CVOut = tuple[pd.DataFrame, tuple[Path, Path]]

VDJ_PAT = "^ID=gene-(IGH|IGK|IGL|TRA|TRB|TRG);"


def main(smk: Any) -> None:
    sconf: cfg.GiabStrats = smk.config
    ws: dict[str, str] = smk.wildcards
    ps: dict[str, str] = smk.params
    inputs: list[Path] = [Path(i) for i in smk.input]
    cds_pattern: str = ps["cds_output"]
    vdj_pattern: str = ps["vdj_output"]

    rk = cfg.wc_to_refkey(ws)
    bk = cfg.wc_to_buildkey(ws)

    bd = sconf.to_build_data(rk, bk)
    fnc = bd.refdata.strat_inputs.functional
    # NOTE: this should never happen because the logic below should prevent the
    # bed file object from propagating if it is None, but the below logic isn't
    # polymorphic enough to understand 'Functional's
    if fnc is None:
        raise DesignError()
    fps = fnc.fparams

    def filter_cds(df: pd.DataFrame) -> pd.DataFrame:
        source_mask = fmap_maybe(lambda x: df[3].str.match(x[0]), fps.source_match)
        offset = 0 if source_mask is None else 1
        type_mask = fmap_maybe(lambda x: df[3 + offset].str.match(x[0]), fps.type_match)
        if source_mask is None:
            if type_mask is None:
                return df
            elif type_mask is not None:
                return df[type_mask].copy()
            else:
                assert_never(type_mask)
        elif source_mask is not None:
            if type_mask is None:
                return df[source_mask].copy()
            elif type_mask is not None:
                return df[source_mask & type_mask].copy()
            else:
                assert_never(type_mask)
        else:
            assert_never(source_mask)

    def filter_vdj(df: pd.DataFrame) -> pd.DataFrame:
        source_offset = 0 if fps.source_match is None else 1
        type_offset = 1 if fps.type_match is None else 1
        return df[df[3 + source_offset + type_offset].str.match(VDJ_PAT)].copy()

    def to_output_pattern(rk: cfg.RefKeyFull) -> tuple[Path, Path]:
        return (
            cfg.sub_output_path(cds_pattern, rk),
            cfg.sub_output_path(vdj_pattern, rk),
        )

    def hap(i: Path, bd: cfg.HapBuildData, bf: cfg.HapBedFile) -> list[CVOut]:
        df = cfg.read_filter_sort_hap_bed(bd, bf, i)
        return [(df, to_output_pattern(bd.refdata.ref.src.key(rk)))]

    def dip1to1(i: Path, bd: cfg.Dip1BuildData, bf: cfg.Dip1BedFile) -> list[CVOut]:
        df = cfg.read_filter_sort_dip1to1_bed(bd, bf, i)
        rfk = bd.refdata.ref.src.key(rk)
        return [(df, to_output_pattern(rfk))]

    def dip1to2(i: Path, bd: cfg.Dip2BuildData, bf: cfg.Dip1BedFile) -> list[CVOut]:
        dfs = cfg.read_filter_sort_dip1to2_bed(bd, bf, i)
        paths = both(to_output_pattern, bd.refdata.ref.src.keys(rk))
        return list(zip(dfs, paths))

    def dip2to1(
        i: tuple[Path, Path],
        bd: cfg.Dip1BuildData,
        bf: cfg.Dip2BedFile,
    ) -> list[CVOut]:
        df = cfg.read_filter_sort_dip2to1_bed(bd, bf, i)
        rfk = bd.refdata.ref.src.key(rk)
        return [(df, to_output_pattern(rfk))]

    def dip2to2(
        i: Path,
        hap: cfg.Haplotype,
        bd: cfg.Dip2BuildData,
        bf: cfg.Dip2BedFile,
    ) -> CVOut:
        df = cfg.read_filter_sort_dip2to2_bed(bd, bf, i, hap)
        paths = to_output_pattern(cfg.RefKeyFull(rk, hap))
        return (df, paths)

    res = sconf.with_build_data_and_bed_i(
        rk,
        bk,
        inputs,
        lambda bd: cfg.bd_to_si(cfg.si_to_functional, bd),
        hap,
        dip1to1,
        dip1to2,
        dip2to1,
        lambda i, bf, bd: cfg.wrap_dip_2to2_i_f(dip2to2, i, bf, bd),
    )

    def write_output(o: Path, ps: list[Path]) -> None:
        with open(o, "w") as f:
            json.dump([str(p) for p in ps], f)

    for df, (cds_path, vdj_path) in res:
        write_bed(cds_path, filter_cds(df))
        if bd.build.include.vdj:
            write_bed(vdj_path, filter_vdj(df))

    write_output(smk.output["cds"], [c for _, (c, _) in res])
    write_output(smk.output["vdj"], [v for _, (_, v) in res])


main(snakemake)  # type: ignore
