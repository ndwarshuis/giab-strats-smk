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


class FGroup(NamedTuple):
    io: IOGroup
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
        # ASSUME: incoming dataframe has the following format:
        # 0: chrom
        # 1: start
        # 2: end
        # 3: attributes
        # 4: source, or type if source is not given
        # 5: type (if given)

        # TODO dirty hack to get this to fail if we try to use a yaml-specified
        # CDS bed file (which is not a GFF and thus cannot be filtered)
        if not isinstance(fnc, cfg.BedFile):
            raise DesignError()
        fps = fnc.cds_params

        source_mask = fmap_maybe(lambda x: df[4].str.match(x[0]), fps.source_match)
        type_col = 4 + (0 if source_mask is None else 1)
        type_mask = fmap_maybe(lambda x: df[type_col].str.match(x[0]), fps.type_match)

        def go(_df: pd.DataFrame) -> pd.DataFrame:
            if source_mask is None:
                if type_mask is None:
                    return _df
                elif type_mask is not None:
                    return _df[type_mask]
                else:
                    assert_never(type_mask)
            elif source_mask is not None:
                if type_mask is None:
                    return _df[source_mask]
                elif type_mask is not None:
                    return _df[source_mask & type_mask]
                else:
                    assert_never(type_mask)
            else:
                assert_never(source_mask)

        return keep_bed_columns(go(df))

    def filter_mhc(df: pd.DataFrame) -> pd.DataFrame:
        raise DesignError("I am not living")

    def filter_kir(df: pd.DataFrame) -> pd.DataFrame:
        raise DesignError("I'm asleep")

    def filter_vdj(df: pd.DataFrame) -> pd.DataFrame:
        return df[df[3].str.match(VDJ_PAT)]

    def to_io_group(name: str, wanted: bool) -> IOGroup:
        return IOGroup(
            inputs=cfg.smk_to_inputs_name(smk, name, True),
            output=cfg.smk_to_output_name(smk, name),
            pattern=cfg.smk_to_param_str(smk, name) if wanted else None,
        )

    def to_f_group(
        name: str,
        f: GFF2Bed,
        g: cfg.BuildDataToBed,
        wanted: bool,
    ) -> FGroup:
        return FGroup(io=to_io_group(name, wanted), gff2bed=f, getbed=g)

    cds = to_io_group("cds", bd.want_cds)
    mhc = to_f_group("mhc", filter_mhc, cfg.bd_to_mhc, bd.want_mhc)
    kir = to_f_group("kir", filter_kir, cfg.bd_to_kir, bd.want_kir)
    vdj = to_f_group("vdj", filter_vdj, cfg.bd_to_vdj, bd.want_vdj)

    want_cds = cds.pattern is not None
    want_any_immuno = (
        kir.io.pattern is not None
        or mhc.io.pattern is not None
        or vdj.io.pattern is not None
    )

    def with_immuno(c: list[Path], f: Callable[[FGroup], list[Path]]) -> Out:
        return (c, f(kir), f(mhc), f(vdj))

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

    def with_cds_bed1(
        out: Callable[[str], Path],
        get_cds: Callable[[], pd.DataFrame],
        get_immuno: Callable[[FGroup, pd.DataFrame | None], list[Path]],
    ) -> Out:
        if want_cds or want_any_immuno:
            df = get_cds()
            c = cds_maybe1(out, filter_cds(df))
            return with_immuno(c, lambda g: get_immuno(g, df))
        else:
            # if I don't want to build any of the outputs, this script shouldn't
            # be called to begin with
            raise DesignError()

    def with_cds_bed2(
        out: Callable[[str], cfg.Double[Path]],
        get_cds: Callable[[], cfg.Double[pd.DataFrame]],
        get_immuno: Callable[[FGroup, cfg.Double[pd.DataFrame] | None], list[Path]],
    ) -> Out:
        if want_cds or want_any_immuno:
            df = get_cds()
            c = cds_maybe2(out, df.map(filter_cds))
            return with_immuno(c, lambda g: get_immuno(g, df))
        else:
            # if I don't want to build any of the outputs, this script shouldn't
            # be called to begin with
            raise DesignError()

    def with_cds_coords1(
        out: Callable[[str], Path],
        df: pd.DataFrame,
        get_immuno: Callable[[FGroup, pd.DataFrame | None], list[Path]],
    ) -> Out:
        c = cds_maybe1(out, df)
        return with_immuno(c, lambda g: get_immuno(g, None))

    def with_cds_coords2(
        out: Callable[[str], cfg.Double[Path]],
        df: cfg.Double[pd.DataFrame],
        get_immuno: Callable[[FGroup, cfg.Double[pd.DataFrame] | None], list[Path]],
    ) -> Out:
        c = cds_maybe2(out, df)
        return with_immuno(c, lambda g: get_immuno(g, None))

    def filter_bed1(g: FGroup, df: pd.DataFrame | None, o: Path) -> None:
        (
            write_bed(o, keep_bed_columns(g.gff2bed(df)))
            if df is not None
            else raise_inline()
        )

    def filter_bed2(
        g: FGroup,
        df: cfg.Double[pd.DataFrame] | None,
        o: cfg.Double[Path],
    ) -> None:
        (
            o.both(
                lambda _o, h: write_bed(_o, keep_bed_columns(g.gff2bed(df.choose(h))))
            )
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

        def go(g: FGroup, df: pd.DataFrame | None) -> list[Path]:
            return with_pattern(
                g.io,
                lambda p: [
                    with_first(
                        out(p),
                        lambda o: cfg.with_inputs_null_hap(
                            g.getbed(bd),
                            g.io.inputs,
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
            return with_cds_bed1(
                out, lambda: cfg.read_filter_sort_hap_bed(bd, bf, i), go
            )

        def with_bedcoords(bc: cfg.HapBedCoords) -> Out:
            return with_cds_coords1(out, cfg.build_hap_coords_df(bd, bc), go)

        return cfg.with_inputs_null_hap(
            bf,
            cds.inputs,
            with_bedfile,
            with_bedcoords,
            lambda: with_immuno([], lambda g: go(g, None)),
        )

    def dip1(rd: cfg.Dip1RefData) -> Out:
        bd = rd.to_build_data_unsafe(bk)
        bf = cfg.bd_to_cds(bd)

        def out(pattern: str) -> Path:
            return cfg.sub_output_path(
                pattern,
                bd.refdata.ref.src.key(bd.refdata.refkey).elem,
            )

        def go(g: FGroup, df: pd.DataFrame | None = None) -> list[Path]:
            return with_pattern(
                g.io,
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
                                    g.io.inputs,
                                    lambda i, bf: cfg.read_write_filter_sort_dip1to1_bed(
                                        i, o, bd, bf
                                    ),
                                    lambda bc: write_bed(
                                        o, cfg.build_dip1to1_coords_df(bd, bc)
                                    ),
                                ),
                                lambda bf: cfg.with_inputs_dip2(
                                    bf,
                                    g.io.inputs,
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
            return with_cds_bed1(
                out, lambda: cfg.read_filter_sort_dip1to1_bed(bd, bf, i), go
            )

        def with_bedcoords1(bc: cfg.Dip1BedCoords) -> Out:
            return with_cds_coords1(out, cfg.build_dip1to1_coords_df(bd, bc), go)

        def with_bedfile2(i: cfg.Double[Path], bf: cfg.Dip2BedFile) -> Out:
            return with_cds_bed1(
                out, lambda: cfg.read_filter_sort_dip2to1_bed(bd, bf, i), go
            )

        def with_bedcoords2(bc: cfg.Dip2BedCoords) -> Out:
            return with_cds_coords1(out, cfg.build_dip2to1_coords_df(bd, bc), go)

        def only_immuno() -> Out:
            return with_immuno([], lambda g: go(g))

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

        def go(g: FGroup, df: cfg.Double[pd.DataFrame] | None = None) -> list[Path]:
            return with_pattern(
                g.io,
                lambda p: with_first(
                    out(p),
                    lambda o: (
                        filter_bed2(g, df, o)
                        if (bf := g.getbed(bd)) is None
                        else cfg.with_dip_bedfile(
                            bf,
                            lambda bf: cfg.with_inputs_dip1(
                                bf,
                                g.io.inputs,
                                lambda i, bf: cfg.read_write_filter_sort_dip1to2_bed(
                                    i, o, bd, bf
                                ),
                                lambda bc: cfg.build_dip1to2_coords_df(bd, bc).both(
                                    lambda df, hap: write_bed(o.choose(hap), df)
                                ),
                            ),
                            lambda bf: cfg.with_inputs_dip2(
                                bf,
                                g.io.inputs,
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
            return with_cds_bed2(
                out, lambda: cfg.read_filter_sort_dip1to2_bed(bd, bf, i), go
            )

        def with_bedcoords1(bc: cfg.Dip1BedCoords) -> Out:
            return with_cds_coords2(out, cfg.build_dip1to2_coords_df(bd, bc), go)

        def with_bedfile2(i: cfg.Double[Path], bf: cfg.Dip2BedFile) -> Out:
            return with_cds_bed2(
                out,
                lambda: i.both(
                    lambda _i, hap: cfg.read_filter_sort_dip2to2_bed(bd, bf, _i, hap)
                ),
                go,
            )

        def with_bedcoords2(bc: cfg.Dip2BedCoords) -> Out:
            df = cfg.make_double(lambda hap: cfg.build_dip2to2_coords_df(hap, bd, bc))
            return with_cds_coords2(out, df, go)

        def only_immuno() -> Out:
            return with_immuno([], lambda g: go(g))

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
    write_output(kir.io.output, k)
    write_output(mhc.io.output, m)
    write_output(vdj.io.output, v)


main(snakemake)  # type: ignore
