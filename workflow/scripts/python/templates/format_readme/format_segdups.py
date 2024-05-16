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

    segdups_src = bd.refdata.strat_inputs.segdups

    if segdups_src is None:
        raise DesignError()

    src_txt = sconf.with_build_data_and_bed_doc(
        rfk,
        bk,
        cfg.smk_to_inputs_name(smk, "segdups_inputs"),
        cfg.bd_to_superdups,
        "The superdups file",
        None,
    )

    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")

    def render_description(t: j2.Template) -> str:
        return t.render(
            segdups_file=cfg.smk_to_param_path(smk, "segdups_path").name,
            not_segdups_file=cfg.smk_to_param_path(smk, "not_segdups_path").name,
            long_segdups_file=cfg.smk_to_param_path(smk, "long_segdups_path").name,
            not_long_segdups_file=cfg.smk_to_param_path(
                smk, "not_long_segdups_path"
            ).name,
        )

    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools", "samtools"})

    def render_methods(t: j2.Template) -> str:
        return t.render(src_txt=src_txt, deps=bedtools_deps)

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
