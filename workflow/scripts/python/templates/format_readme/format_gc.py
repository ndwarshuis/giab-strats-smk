import jinja2 as j2
from typing import Any
import common.config as cfg
from common.functional import DesignError
import template_utils as tu

# ASSUME ranges are non-empty, and dual extremes are non-empty


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)

    paths = smk.params["paths"]

    if not isinstance(paths, cfg.GCPaths):
        raise DesignError()

    def render_description(t: j2.Template) -> str:
        return t.render(
            low_file=paths.lowGC.name,
            high_file=paths.highGC.name,
            range_files=[p.name for p in paths.middleGC],
            dual_files=[p.name for p in paths.extremes],
        )

    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")
    seqtk_env_path = cfg.smk_to_input_name(smk, "seqtk_env")

    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools"})
    seqtk_deps = tu.env_dependencies(seqtk_env_path, {"seqtk"})

    def render_methods(t: j2.Template) -> str:
        return t.render(
            lower_fractions=paths.params.low_fractions,
            upper_fractions=paths.params.high_fractions,
            deps={**bedtools_deps, **seqtk_deps},
        )

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        "GC rich/poor regions",
        cfg.CoreLevel.GC.value,
        sconf.refkey_haplotypes(rfk),
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
