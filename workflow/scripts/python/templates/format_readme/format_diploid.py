import jinja2 as j2
from pathlib import Path
from typing import Any
import common.config as cfg
from common.functional import DesignError
import template_utils as tu


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)

    paths = smk.params["paths"]

    if not isinstance(paths, cfg.DiploidPaths):
        raise DesignError()

    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")
    dipcall_env_path = cfg.smk_to_input_name(smk, "diploid_env")

    def fmt_names(ps: list[Path]) -> list[str]:
        return [tu.sub_rk(rfk, p.name) for p in ps]

    def render_description(t: j2.Template) -> str:
        return t.render(
            het_files=fmt_names(paths.hets),
            SNVorSV_het_files=fmt_names(paths.SNVorSV_hets),
            hom_files=fmt_names(paths.homs),
            SNVorSV_hom_files=fmt_names(paths.SNVorSV_homs),
        )

    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools", "samtools"})
    dipcall_deps = tu.env_dependencies(
        dipcall_env_path, {"k8", "minimap2", "htsbox", "samtools", "bedtools"}
    )

    def render_methods(t: j2.Template) -> str:
        return t.render(
            bedtools_deps=bedtools_deps,
            dipcall_deps=dipcall_deps,
        )

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        "regions specific to diploid references",
        cfg.CoreLevel.DIPLOID,
        sconf.refkey_haplotypes(rfk),
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
