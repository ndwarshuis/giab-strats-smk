import jinja2 as j2
from typing import Any
import common.config as cfg
from common.functional import DesignError
import template_utils as tu


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)
    rk = cfg.strip_full_refkey(rfk)
    bk = cfg.wc_to_buildkey(ws)
    bd = sconf.to_build_data(rk, bk)

    paths: cfg.SexPaths = smk.params["paths"]

    if not isinstance(paths, cfg.SexPaths):
        return DesignError()

    xy = bd.refdata.strat_inputs.xy

    if xy is None:
        raise DesignError()

    xy_features = xy.features

    if xy_features is None:
        raise DesignError()

    # TODO not DRY (XTR and Ampliconic used elsewhere, define constants)
    filter_phrase = " and ".join(
        [
            f"'{x}'"
            for x, test in [
                (xy_features.xtr, "XTR"),
                (xy_features.ampliconic, "Ampliconic"),
            ]
            if test
        ]
    )

    def feature_src_txt(use_x: bool) -> str:
        x_bed, name = (xy_features.x_bed, "X") if use_x else (xy_features.y_bed, "Y")

        src_para = cfg.format_src(
            x_bed.bed.documentation.elem,
            cfg.smk_to_input_name(smk, f"{name}_features_inputs"),
            f"The {name} chromosome feature bed",
        )
        params_para = cfg.format_bed_params(x_bed.params)
        filter_para = (
            "Specific features were isolated by filtering column "
            f"{x_bed.level_col + 1} for {filter_phrase}"
        )

        return "\n\n".join(
            [cfg.readme_fill(p) for p in [src_para, params_para, filter_para]]
        )

    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")

    def render_description(t: j2.Template) -> str:
        return t.render(
            all_auto_file=cfg.smk_to_param_path(smk, "auto_path").name,
            par_files=[p.name for p in cfg.smk_to_param_paths(smk, "par_paths")],
            nonpar_files=[p.name for p in cfg.smk_to_param_paths(smk, "nonpar_paths")],
            ampliconic_files=[
                p.name for p in cfg.smk_to_param_paths(smk, "ampliconic_paths")
            ],
            xtr_files=[p.name for p in cfg.smk_to_param_paths(smk, "xtr_paths")],
        )

    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools", "samtools"})

    # TODO need a way to control which X and/or Y are included
    def render_methods(t: j2.Template) -> str:
        return t.render(
            deps=bedtools_deps,
        )

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        "segmental duplications",
        cfg.CoreLevel.SEGDUPS,
        sconf.refkey_haplotypes(rfk),
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
