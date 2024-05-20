from pathlib import Path
import jinja2 as j2
from typing import Any
import common.config as cfg
import template_utils as tu


def relative_path(p: Path) -> str:
    return str(Path("..") / p.parent.name / p.name)


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)

    # NOTE: in theory I could have made this flexible enough that each of the
    # inputs was optional, since the only thing 'union' implies is a
    # multiintersect/merge. For other reasons, the only thing that is variable
    # here is the xy input, thus only the following is needed:
    xy = cfg.smk_to_inputs_name(smk, "xy")
    xy_txt = "" if xy == [] else "difficult XY regions"

    alldifficult_desc = (
        f"This contains the above regions plus {xy_txt}tandem "
        "repeats, and high/low GC regions."
    )

    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")

    def render_description(t: j2.Template) -> str:
        return t.render(
            segdup_map_file=cfg.smk_to_param_path(smk, "segdup_map_path").name,
            not_segdup_map_file=cfg.smk_to_param_path(smk, "not_segdup_map_path").name,
            alldifficult_file=cfg.smk_to_param_path(smk, "alldifficult_path").name,
            not_alldifficult_file=cfg.smk_to_param_path(
                smk, "not_alldifficult_path"
            ).name,
            alldifficult_desc=alldifficult_desc,
        )

    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools", "samtools"})

    all_difficult_files = [
        relative_path(p)
        for p in [
            cfg.smk_to_input_name(smk, "repeats"),
            *cfg.smk_to_inputs_name(smk, "xy"),
            *([cfg.smk_to_input_name(smk, "gc")] if "gc" in smk.input else []),
        ]
    ]

    def render_methods(t: j2.Template) -> str:
        return t.render(
            map_file=relative_path(cfg.smk_to_input_name(smk, "lowmap")),
            segdup_file=relative_path(cfg.smk_to_input_name(smk, "segdups")),
            all_difficult_files=all_difficult_files,
            deps=bedtools_deps,
        )

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        "unions of other difficult stratifications",
        cfg.CoreLevel.UNION,
        sconf.refkey_haplotypes(rfk),
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
