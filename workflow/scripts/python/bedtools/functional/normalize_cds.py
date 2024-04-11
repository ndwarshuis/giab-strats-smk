import pandas as pd
import json
from pathlib import Path
from typing import Any, assert_never, NamedTuple, Callable
import common.config as cfg
from common.bed import write_bed
from common.functional import DesignError, fmap_maybe, match12_unsafe

GFF2Bed = Callable[[pd.DataFrame], pd.DataFrame]


class IOGroup(NamedTuple):
    inputs: list[Path]
    output: Path
    pattern: str
    gff2bed: GFF2Bed


class GFFOut(NamedTuple):
    df: pd.DataFrame
    rk: cfg.RefKeyFull


VDJ_PAT = "^ID=gene-(IGH|IGK|IGL|TRA|TRB|TRG);"


def main(smk: Any) -> None:
    sconf: cfg.GiabStrats = smk.config
    ws: dict[str, str] = smk.wildcards

    rk = cfg.wc_to_refkey(ws)
    bk = cfg.wc_to_buildkey(ws)

    bd = sconf.to_build_data(rk, bk)
    fnc = bd.refdata.strat_inputs.functional.cds
    # NOTE: this should never happen because the logic below should prevent the
    # bed file object from propagating if it is None, but the below logic isn't
    # polymorphic enough to understand 'Functional's
    if fnc is None:
        raise DesignError()
    fps = fnc.cds_params

    # functions to convert the GFF dataframe to whatever we feel like

    def keep_bed_columns(df: pd.DataFrame) -> pd.DataFrame:
        return df[[0, 1, 2]].copy()

    def filter_cds(df: pd.DataFrame) -> pd.DataFrame:
        # ASSUME: income dataframe has the following format:
        # 0: chrom
        # 1: start
        # 2: end
        # 3: attributes
        # 4: source, or type if source is not given
        # 5: type (if given)
        source_mask = fmap_maybe(lambda x: df[4].str.match(x[0]), fps.source_match)
        type_col = 4 + (0 if source_mask is None else 1)
        type_mask = fmap_maybe(lambda x: df[type_col].str.match(x[0]), fps.type_match)
        if source_mask is None:
            if type_mask is None:
                return df
            elif type_mask is not None:
                return df[type_mask]
            else:
                assert_never(type_mask)
        elif source_mask is not None:
            if type_mask is None:
                return df[source_mask]
            elif type_mask is not None:
                return df[source_mask & type_mask]
            else:
                assert_never(type_mask)
        else:
            assert_never(source_mask)

    def filter_mhc(df: pd.DataFrame) -> pd.DataFrame:
        raise DesignError("I am not living")

    def filter_kir(df: pd.DataFrame) -> pd.DataFrame:
        raise DesignError("I'm asleep")

    def filter_vdj(df: pd.DataFrame) -> pd.DataFrame:
        return df[df[3].str.match(VDJ_PAT)]

    # functions to read the GFF dataframe

    def hap(i: Path, bd: cfg.HapBuildData, bf: cfg.HapBedFile) -> list[GFFOut]:
        df = cfg.read_filter_sort_hap_bed(bd, bf, i)
        return [GFFOut(df, bd.refdata.ref.src.key(rk))]

    def dip1to1(i: Path, bd: cfg.Dip1BuildData, bf: cfg.Dip1BedFile) -> list[GFFOut]:
        df = cfg.read_filter_sort_dip1to1_bed(bd, bf, i)
        rfk = bd.refdata.ref.src.key(rk)
        return [GFFOut(df, rfk)]

    def dip1to2(i: Path, bd: cfg.Dip2BuildData, bf: cfg.Dip1BedFile) -> list[GFFOut]:
        dfs = cfg.read_filter_sort_dip1to2_bed(bd, bf, i)
        rks = bd.refdata.ref.src.keys(rk)
        return [GFFOut(d, k) for d, k in zip(dfs, rks)]

    def dip2to1(
        i: tuple[Path, Path],
        bd: cfg.Dip1BuildData,
        bf: cfg.Dip2BedFile,
    ) -> list[GFFOut]:
        df = cfg.read_filter_sort_dip2to1_bed(bd, bf, i)
        rfk = bd.refdata.ref.src.key(rk)
        return [GFFOut(df, rfk)]

    def dip2to2(
        i: Path,
        hap: cfg.Haplotype,
        bd: cfg.Dip2BuildData,
        bf: cfg.Dip2BedFile,
    ) -> GFFOut:
        df = cfg.read_filter_sort_dip2to2_bed(bd, bf, i, hap)
        return GFFOut(df, cfg.RefKeyFull(rk, hap))

    # other functions

    def write_bedpaths(p: Path, ps: list[Path]) -> None:
        with open(p, "w") as f:
            json.dump([str(p) for p in ps], f)

    def write_from_gff(res: list[GFFOut], test: bool, g: IOGroup) -> None:
        if test:

            def write1(res: GFFOut) -> list[Path]:
                bedpath = cfg.sub_output_path(g.pattern, res.rk)
                write_bed(bedpath, keep_bed_columns(g.gff2bed(res.df)))
                return [bedpath]

            def write2(res1: GFFOut, res2: GFFOut) -> list[Path]:
                return write1(res1) + write1(res2)

            bedpaths = match12_unsafe(res, write1, write2)
        else:
            bedpaths = []

        write_bedpaths(g.output, bedpaths)

    def write_region(
        res: list[GFFOut],
        test: bool,
        g: IOGroup,
        bd2bed: cfg.BuildDataToBed,
    ) -> None:
        if len(g.inputs) == 0:
            write_from_gff(res, test, g)
        else:
            if test:
                bedpaths = cfg.filter_sort_bed_main_inner(
                    sconf, rk, bk, g.inputs, g.output, g.pattern, bd2bed
                )
            else:
                bedpaths = []

            write_bedpaths(g.output, bedpaths)

    def to_io_group(name: str, f: GFF2Bed) -> IOGroup:
        return IOGroup(
            inputs=cfg.smk_to_inputs_name(smk, name),
            output=cfg.smk_to_output_name(smk, name),
            pattern=cfg.smk_to_param_str(smk, name),
            gff2bed=f,
        )

    # now actually do things...

    cds = to_io_group("cds", filter_cds)
    mhc = to_io_group("mhc", filter_mhc)
    kir = to_io_group("kir", filter_kir)
    vdj = to_io_group("vdj", filter_vdj)

    # first, read the GFF from the CDS source path
    gff = sconf.with_build_data_and_bed_i(
        rk,
        bk,
        cds.inputs,
        lambda bd: cfg.bd_to_si(cfg.si_to_cds, bd),
        hap,
        dip1to1,
        dip1to2,
        dip2to1,
        lambda i, bf, bd: list(cfg.wrap_dip_2to2_i_f(dip2to2, i, bf, bd)),
    )

    # second, write the CDS bed file (assuming we want it)
    write_from_gff(gff, bd.want_cds, cds)

    # third, convert the GFF to MHC/KIR/VDJ bed files (assuming we want them)
    write_region(gff, bd.want_mhc, mhc, cfg.bd_to_mhc)
    write_region(gff, bd.want_kir, kir, cfg.bd_to_kir)
    write_region(gff, bd.want_vdj, vdj, cfg.bd_to_vdj)


main(snakemake)  # type: ignore
