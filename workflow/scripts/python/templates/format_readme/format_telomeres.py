import jinja2 as j2
from typing import Any
import common.config as cfg
from common.functional import DesignError
import template_utils as tu


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)

    paths = smk.params["paths"]

    if not isinstance(paths, cfg.TelomerePaths):
        raise DesignError()

    bedtools_env_path = cfg.smk_to_input_name(smk, "seqtk_env")
    seqtk_deps = tu.env_dependencies(bedtools_env_path, {"seqtk"})

    def render_description(t: j2.Template) -> str:
        return t.render(file=tu.sub_rk(rfk, paths.telomeres.name))

    def render_methods(t: j2.Template) -> str:
        return t.render(deps=seqtk_deps)

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        "telomeric regions",
        cfg.CoreLevel.TELOMERES,
        sconf.refkey_haplotypes(rfk),
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
