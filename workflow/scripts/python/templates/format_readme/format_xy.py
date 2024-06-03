import jinja2 as j2
from typing import Any, assert_never
import common.config as cfg
from common.functional import DesignError, fmap_maybe
import template_utils as tu


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)

    paths = smk.params["paths"]

    if not isinstance(paths, cfg.SexPaths):
        raise DesignError()

    auto_path = paths.auto
    sex = paths.sex

    def feature_src_txt(use_x: bool, p: cfg.XYFeaturePaths) -> str:
        name = "X" if use_x else "Y"

        filter_phrase = " and ".join(
            [f"'{x}'" for x in [p.xtr, p.ampliconic] if x is not None]
        )

        src_para = cfg.format_src(
            p.bed.bed.documentation.elem,
            p.src,
            f"The {name} chromosome feature bed",
        )
        params_para = cfg.format_bed_params(p.bed.params)
        filter_para = (
            "Specific features were isolated by filtering column "
            f"{p.bed.level_col + 1} for {filter_phrase}"
        )

        return "\n\n".join(
            [cfg.readme_fill(p) for p in [src_para, params_para, filter_para]]
        )

    def to_src(
        x: cfg.SubSexPaths | None, y: cfg.SubSexPaths | None
    ) -> dict[str, str | None]:
        def get_par_src(p: cfg.SubSexPaths | None) -> str | None:
            return fmap_maybe(lambda z: fmap_maybe(lambda p: p.doc, z.par), p)

        def get_feature_src(p: cfg.SubSexPaths | None, use_x: bool) -> str | None:
            return fmap_maybe(
                lambda z: fmap_maybe(lambda p: feature_src_txt(use_x, p), z.features), p
            )

        return {
            "x_par_src": get_par_src(x),
            "y_par_src": get_par_src(y),
            "x_features_src": get_feature_src(x, True),
            "y_features_src": get_feature_src(y, False),
        }

    if isinstance(sex, cfg.MaleHapSexPaths):
        src = to_src(sex.x, sex.y)
    elif isinstance(sex, cfg.Dip1SexPaths):
        src = to_src(sex.sex1, sex.sex2)
    elif isinstance(sex, cfg.Dip2SexPaths):
        src = (
            to_src(None, sex.paths)
            if sex.hap is cfg.Haplotype.PAT
            else to_src(sex.paths, None)
        )
    else:
        assert_never(sex)

    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")

    def render_description(t: j2.Template) -> str:
        return t.render(
            all_auto_file=fmap_maybe(lambda x: x.name, auto_path),
            par_files=[p.name for p in sex.par_paths],
            non_par_files=[p.name for p in sex.nonpar_paths],
            ampliconic_files=[p.name for p in sex.ampliconic_paths],
            xtr_files=[p.name for p in sex.xtr_paths],
        )

    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools", "samtools"})

    def render_methods(t: j2.Template) -> str:
        return t.render(
            **src,
            deps=bedtools_deps,
        )

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        "the sex chromosomes",
        cfg.CoreLevel.XY,
        sconf.refkey_haplotypes(rfk),
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
