import pandas as pd
import json
from pathlib import Path
from more_itertools import flatten
from typing import Any, assert_never, NamedTuple, Callable
import common.config as cfg
from common.bed import write_bed
from common.functional import DesignError, fmap_maybe

GFF2Bed = Callable[[pd.DataFrame], pd.DataFrame]


class IOGroup(NamedTuple):
    inputs: list[Path]
    output: Path
    pattern: str
    gff2bed: GFF2Bed


class GFFOut(NamedTuple):
    df: pd.DataFrame
    rk: cfg.RefKeyFull


# TODO don't hardcode this
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

        # TODO dirty hack to get this to fail error if we try to call for a
        # manually-specific CDS bed file (which is not a GFF and thus cannot be
        # filtered)
        if not isinstance(fnc, cfg.BedFile):
            raise DesignError()
        fps = fnc.cds_params

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

    def to_io_group(name: str, f: GFF2Bed, allow_empty: bool) -> IOGroup:
        return IOGroup(
            inputs=cfg.smk_to_inputs_name(smk, name, allow_empty),
            output=cfg.smk_to_output_name(smk, name),
            pattern=cfg.smk_to_param_str(smk, name),
            gff2bed=f,
        )

    cds = to_io_group("cds", filter_cds, False)
    mhc = to_io_group("mhc", filter_mhc, True)
    kir = to_io_group("kir", filter_kir, True)
    vdj = to_io_group("vdj", filter_vdj, True)

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

            def write2(res: cfg.Double[GFFOut]) -> list[Path]:
                return [*flatten(res.map(write1).as_list)]

            bedpaths = cfg.match12_unsafe(res, write1, write2)
        else:
            bedpaths = []

        write_bedpaths(g.output, bedpaths)

    def write_region_inner(test: bool, g: IOGroup, bd2bed: cfg.BuildDataToBed) -> None:
        if test:
            bedpaths = cfg.filter_sort_bed_main_inner(
                sconf, rk, bk, g.inputs, g.output, g.pattern, bd2bed
            )
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
            write_region_inner(test, g, bd2bed)

    def write_gff_all(gff: list[GFFOut]) -> None:
        write_from_gff(gff, bd.want_cds, cds)
        write_region(gff, bd.want_mhc, mhc, cfg.bd_to_mhc)
        write_region(gff, bd.want_kir, kir, cfg.bd_to_kir)
        write_region(gff, bd.want_vdj, vdj, cfg.bd_to_vdj)

    def write_no_gff_all() -> None:
        write_region_inner(bd.want_cds, cds, cfg.bd_to_cds)
        write_region_inner(bd.want_mhc, mhc, cfg.bd_to_mhc)
        write_region_inner(bd.want_kir, kir, cfg.bd_to_kir)
        write_region_inner(bd.want_vdj, vdj, cfg.bd_to_vdj)

    # functions to read the GFF dataframe

    def hap(bd: cfg.HapBuildData, bf: cfg.HapBedFileOrCoords) -> None:
        if isinstance(bf, cfg.BedFile):
            gff = cfg.match1_unsafe(
                cds.inputs,
                lambda i: [
                    GFFOut(
                        cfg.read_filter_sort_hap_bed(bd, bf, i),
                        bd.refdata.ref.src.key(rk).elem,
                    )
                ],
            )
            write_gff_all(gff)
        else:
            write_no_gff_all()

    def dip1to1(bd: cfg.Dip1BuildData, bf: cfg.Dip1BedFileOrCoords) -> None:
        if isinstance(bf, cfg.BedFile):
            gff = cfg.match1_unsafe(
                cds.inputs,
                lambda i: [
                    GFFOut(
                        cfg.read_filter_sort_dip1to1_bed(bd, bf, i),
                        bd.refdata.ref.src.key(rk).elem,
                    )
                ],
            )
            write_gff_all(gff)
        else:
            write_no_gff_all()

    def dip1to2(bd: cfg.Dip2BuildData, bf: cfg.Dip1BedFileOrCoords) -> None:
        if isinstance(bf, cfg.BedFile):
            gff = cfg.match1_unsafe(
                cds.inputs,
                lambda i: cfg.read_filter_sort_dip1to2_bed(bd, bf, i).both(
                    lambda df, hap: GFFOut(df, bd.refdata.ref.src.keys(rk).choose(hap))
                ),
            )
            write_gff_all(gff.as_list)
        else:
            write_no_gff_all()

    def dip2to1(
        bd: cfg.Dip1BuildData,
        bf: cfg.Dip2BedFileOrCoords,
    ) -> None:
        if isinstance(bf, cfg.BedFile):
            gff = cfg.match2_unsafe(
                cds.inputs,
                lambda i: [
                    GFFOut(
                        cfg.read_filter_sort_dip2to1_bed(bd, bf, i),
                        bd.refdata.ref.src.key(rk).elem,
                    )
                ],
            )
            write_gff_all(gff)
        else:
            write_no_gff_all()

    def dip2to2(bd: cfg.Dip2BuildData, bf: cfg.Dip2BedFileOrCoords) -> None:
        if isinstance(bf, cfg.BedFile):
            gff = cfg.match2_unsafe(
                cds.inputs,
                lambda i2: i2.both(
                    lambda i, hap: GFFOut(
                        cfg.read_filter_sort_dip2to2_bed(bd, bf, i, hap),
                        bd.refdata.ref.src.keys(rk).choose(hap),
                    )
                ),
            )
            write_gff_all(gff.as_list)
        else:
            write_no_gff_all()

    sconf.with_build_data_and_bed(
        rk,
        bk,
        cfg.bd_to_cds,
        hap,
        dip1to1,
        dip1to2,
        dip2to1,
        dip2to2,
    )


main(snakemake)  # type: ignore
