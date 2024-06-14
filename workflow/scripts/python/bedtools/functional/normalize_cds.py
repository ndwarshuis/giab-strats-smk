import pandas as pd
import json
from pathlib import Path
from typing import Any, assert_never, NamedTuple, Callable, TypeVar
import common.config as cfg
from common.bed import write_bed
from common.functional import DesignError, fmap_maybe, raise_inline, with_first

# This is a rather complicated script due to the fact that any of the inputs
# might be blank. Since we can (in theory) make any of the immunological regions
# from the CDS regions, we don't need to supply the immuno regions if we have
# the CDS. Likewise, the immuno regions can be specified directly as coordinates
# (like all other bed files) and thus won't rely on the CDS regions in this
# case.
#
# As a consequence, the logic for this is unique wrt to the other scripts which
# read bed files, due to the fact that it is perfectly normal for a bed getter
# function to return None (ie no bed supplied).

X = TypeVar("X")

GFF2Bed = Callable[[pd.DataFrame], pd.DataFrame]


class IOGroup(NamedTuple):
    inputs: list[Path]
    output: Path
    pattern: str | None
    gff2bed: GFF2Bed
    getbed: cfg.BuildDataToBed


Out = tuple[list[Path], list[Path], list[Path], list[Path]]


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

    def to_io_group(
        name: str,
        f: GFF2Bed,
        g: cfg.BuildDataToBed,
        wanted: bool,
    ) -> IOGroup:
        return IOGroup(
            inputs=cfg.smk_to_inputs_name(smk, name, True),
            output=cfg.smk_to_output_name(smk, name),
            pattern=cfg.smk_to_param_str(smk, name) if wanted else None,
            gff2bed=f,
            getbed=g,
        )

    cds = to_io_group("cds", filter_cds, cfg.bd_to_cds, bd.want_cds)
    mhc = to_io_group("mhc", filter_mhc, cfg.bd_to_mhc, bd.want_mhc)
    kir = to_io_group("kir", filter_kir, cfg.bd_to_kir, bd.want_kir)
    vdj = to_io_group("vdj", filter_vdj, cfg.bd_to_vdj, bd.want_vdj)

    def with_immuno(f: Callable[[IOGroup], X]) -> tuple[X, X, X]:
        return (f(kir), f(mhc), f(vdj))

    def with_pattern(g: IOGroup, f: Callable[[str], list[X]]) -> list[X]:
        if g.pattern is not None:
            return f(g.pattern)
        else:
            return []

    def cds_maybe1(out: Callable[[str], Path], df: pd.DataFrame) -> list[Path]:
        return with_pattern(
            cds, lambda p: [with_first(out(p), lambda o: write_bed(o, df))]
        )

    def cds_maybe2(
        out: Callable[[str], cfg.Double[Path]],
        df: cfg.Double[pd.DataFrame],
    ) -> list[Path]:
        return with_pattern(
            cds,
            lambda p: with_first(
                out(p),
                lambda o: o.both(lambda _o, hap: write_bed(_o, df.choose(hap))),
            ).as_list,
        )

    def filter_bed1(g: IOGroup, df: pd.DataFrame | None, o: Path) -> None:
        write_bed(o, g.gff2bed(df)) if df is not None else raise_inline()

    def filter_bed2(
        g: IOGroup,
        df: cfg.Double[pd.DataFrame] | None,
        o: cfg.Double[Path],
    ) -> None:
        (
            o.both(lambda _o, h: write_bed(_o, g.gff2bed(df.choose(h))))
            if df is not None
            else raise_inline()
        )

    def hap(rd: cfg.HapRefData) -> Out:
        bd = rd.to_build_data_unsafe(bk)
        bf = cfg.bd_to_cds(bd)

        def out(pattern: str) -> Path:
            return cfg.sub_output_path(
                pattern,
                bd.refdata.ref.src.key(bd.refdata.refkey).elem,
            )

        def go(g: IOGroup, df: pd.DataFrame | None) -> list[Path]:
            return with_pattern(
                g,
                lambda p: [
                    with_first(
                        out(p),
                        lambda o: cfg.with_inputs_null_hap(
                            g.getbed(bd),
                            g.inputs,
                            lambda i, bf: cfg.read_write_filter_sort_hap_bed(
                                i, o, bd, bf
                            ),
                            lambda bc: write_bed(o, cfg.build_hap_coords_df(bd, bc)),
                            lambda: filter_bed1(g, df, o),
                        ),
                    )
                ],
            )

        def with_bedfile(i: Path, bf: cfg.HapBedFile) -> Out:
            df = cfg.read_filter_sort_hap_bed(bd, bf, i)
            c = cds_maybe1(out, df)
            x = with_immuno(lambda g: go(g, df))
            return (c, *x)

        def with_bedcoords(bc: cfg.HapBedCoords) -> Out:
            df = cfg.build_hap_coords_df(bd, bc)
            c = cds_maybe1(out, df)
            x = with_immuno(lambda g: go(g, None))
            return (c, *x)

        def only_immuno() -> Out:
            x = with_immuno(lambda g: go(g, None))
            return ([], *x)

        return cfg.with_inputs_null_hap(
            bf,
            cds.inputs,
            with_bedfile,
            with_bedcoords,
            only_immuno,
        )

    def dip1(rd: cfg.Dip1RefData) -> Out:
        bd = rd.to_build_data_unsafe(bk)
        bf = cfg.bd_to_cds(bd)

        def out(pattern: str) -> Path:
            return cfg.sub_output_path(
                pattern,
                bd.refdata.ref.src.key(bd.refdata.refkey).elem,
            )

        def go(g: IOGroup, df: pd.DataFrame | None = None) -> list[Path]:
            return with_pattern(
                g,
                lambda p: [
                    with_first(
                        out(p),
                        lambda o: (
                            filter_bed1(g, df, o)
                            if (bf := g.getbed(bd)) is None
                            else cfg.with_dip_bedfile(
                                bf,
                                lambda bf: cfg.with_inputs_dip1(
                                    bf,
                                    g.inputs,
                                    lambda i, bf: cfg.read_write_filter_sort_dip1to1_bed(
                                        i, o, bd, bf
                                    ),
                                    lambda bc: write_bed(
                                        o, cfg.build_dip1to1_coords_df(bd, bc)
                                    ),
                                ),
                                lambda bf: cfg.with_inputs_dip2(
                                    bf,
                                    g.inputs,
                                    lambda i, bf: cfg.read_write_filter_sort_dip2to1_bed(
                                        i, o, bd, bf
                                    ),
                                    lambda bc: write_bed(
                                        o, cfg.build_dip2to1_coords_df(bd, bc)
                                    ),
                                ),
                            )
                        ),
                    ),
                ],
            )

        def with_bedfile1(i: Path, bf: cfg.Dip1BedFile) -> Out:
            df = cfg.read_filter_sort_dip1to1_bed(bd, bf, i)
            c = cds_maybe1(out, df)
            x = with_immuno(lambda g: go(g, df))
            return (c, *x)

        def with_bedcoords1(bc: cfg.Dip1BedCoords) -> Out:
            df = cfg.build_dip1to1_coords_df(bd, bc)
            c = cds_maybe1(out, df)
            x = with_immuno(lambda g: go(g))
            return (c, *x)

        def with_bedfile2(i: cfg.Double[Path], bf: cfg.Dip2BedFile) -> Out:
            df = cfg.read_filter_sort_dip2to1_bed(bd, bf, i)
            c = cds_maybe1(out, df)
            x = with_immuno(lambda g: go(g, df))
            return (c, *x)

        def with_bedcoords2(bc: cfg.Dip2BedCoords) -> Out:
            df = cfg.build_dip2to1_coords_df(bd, bc)
            c = cds_maybe1(out, df)
            x = with_immuno(lambda g: go(g))
            return (c, *x)

        def only_immuno() -> Out:
            x = with_immuno(lambda g: go(g))
            return ([], *x)

        if bf is None:
            return only_immuno()
        else:
            return cfg.with_dip_bedfile(
                bf,
                lambda bf: cfg.with_inputs_null_dip1(
                    bf,
                    cds.inputs,
                    with_bedfile1,
                    with_bedcoords1,
                    only_immuno,
                ),
                lambda bf: cfg.with_inputs_null_dip2(
                    bf,
                    cds.inputs,
                    with_bedfile2,
                    with_bedcoords2,
                    only_immuno,
                ),
            )

    def dip2(rd: cfg.Dip2RefData) -> Out:
        bd = rd.to_build_data_unsafe(bk)
        bf = cfg.bd_to_cds(bd)

        def out(pattern: str) -> cfg.Double[Path]:
            return bd.refdata.ref.src.keys(bd.refdata.refkey).map(
                lambda k: cfg.sub_output_path(pattern, k)
            )

        def go(g: IOGroup, df: cfg.Double[pd.DataFrame] | None = None) -> list[Path]:
            return with_pattern(
                g,
                lambda p: with_first(
                    out(p),
                    lambda o: (
                        filter_bed2(g, df, o)
                        if (bf := g.getbed(bd)) is None
                        else cfg.with_dip_bedfile(
                            bf,
                            lambda bf: cfg.with_inputs_dip1(
                                bf,
                                g.inputs,
                                lambda i, bf: cfg.read_write_filter_sort_dip1to2_bed(
                                    i, o, bd, bf
                                ),
                                lambda bc: cfg.build_dip1to2_coords_df(bd, bc).both(
                                    lambda df, hap: write_bed(o.choose(hap), df)
                                ),
                            ),
                            lambda bf: cfg.with_inputs_dip2(
                                bf,
                                g.inputs,
                                lambda i, bf: o.both(
                                    lambda _o, hap: cfg.read_write_filter_sort_dip2to2_bed(
                                        i.choose(hap), _o, hap, bd, bf
                                    )
                                ),
                                lambda bc: o.both(
                                    lambda _o, hap: write_bed(
                                        _o,
                                        cfg.build_dip2to2_coords_df(hap, bd, bc),
                                    )
                                ),
                            ),
                        )
                    ),
                ).as_list,
            )

        def with_bedfile1(i: Path, bf: cfg.Dip1BedFile) -> Out:
            df = cfg.read_filter_sort_dip1to2_bed(bd, bf, i)
            c = cds_maybe2(out, df)
            x = with_immuno(lambda g: go(g, df))
            return (c, *x)

        def with_bedcoords1(bc: cfg.Dip1BedCoords) -> Out:
            df = cfg.build_dip1to2_coords_df(bd, bc)
            c = cds_maybe2(out, df)
            x = with_immuno(lambda g: go(g))
            return (c, *x)

        def with_bedfile2(i: cfg.Double[Path], bf: cfg.Dip2BedFile) -> Out:
            df = i.both(
                lambda _i, hap: cfg.read_filter_sort_dip2to2_bed(bd, bf, _i, hap)
            )
            c = cds_maybe2(out, df)
            x = with_immuno(lambda g: go(g, df))
            return (c, *x)

        def with_bedcoords2(bc: cfg.Dip2BedCoords) -> Out:
            df = cfg.make_double(lambda hap: cfg.build_dip2to2_coords_df(hap, bd, bc))
            c = cds_maybe2(out, df)
            x = with_immuno(lambda g: go(g))
            return (c, *x)

        def only_immuno() -> Out:
            x = with_immuno(lambda g: go(g))
            return ([], *x)

        if bf is None:
            return only_immuno()
        else:
            return cfg.with_dip_bedfile(
                bf,
                lambda bf: cfg.with_inputs_null_dip1(
                    bf,
                    cds.inputs,
                    with_bedfile1,
                    with_bedcoords1,
                    only_immuno,
                ),
                lambda bf: cfg.with_inputs_null_dip2(
                    bf,
                    cds.inputs,
                    with_bedfile2,
                    with_bedcoords2,
                    only_immuno,
                ),
            )

    def write_output(p: Path, ps: list[Path]) -> None:
        with open(p, "w") as f:
            json.dump([str(p) for p in ps], f)

    c, k, m, v = sconf.with_ref_data(rk, hap, dip1, dip2)

    write_output(cds.output, c)
    write_output(kir.output, k)
    write_output(mhc.output, m)
    write_output(vdj.output, v)


main(snakemake)  # type: ignore
