import json
import jinja2 as j2
from os.path import basename
from typing import Any
import common.config as cfg
from common.functional import DesignError
import template_utils as tu

# ASSUME ranges are non-empty, and dual extremes are non-empty


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)
    rk = cfg.strip_full_refkey(rfk)
    bk = cfg.wc_to_buildkey(ws)
    bd = sconf.to_build_data(rk, bk)
    gc_params = bd.build.include.gc

    if gc_params is None:
        raise DesignError()

    gc_inputs_path = cfg.smk_to_input_name(smk, "gc_inputs")

    with open(gc_inputs_path, "r") as f:
        u = json.load(f)
        dual_extremes = [
            basename(p) for p in [u["widest_extreme"]] + u["other_extremes"]
        ]
        ranges = [basename(p) for p in u["gc_ranges"]]
        low_extreme = ranges[0]
        high_extreme = ranges[-1]
        middle_ranges = ranges[1:-1]

    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")
    seqtk_env_path = cfg.smk_to_input_name(smk, "seqtk_env")

    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools", "samtools"})
    seqtk_deps = tu.env_dependencies(seqtk_env_path, {"seqtk"})

    def render_description(t: j2.Template) -> str:
        return t.render(
            low_file=low_extreme,
            high_file=high_extreme,
            range_files=middle_ranges,
            dual_files=dual_extremes,
        )

    def render_methods(t: j2.Template) -> str:
        return t.render(
            lower_fractions=gc_params.low_fractions,
            upper_fractions=gc_params.high_fractions,
            deps={**bedtools_deps, **seqtk_deps},
        )

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        "coding regions",
        cfg.CoreLevel.FUNCTIONAL,
        sconf.refkey_haplotypes(rfk),
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore