import jinja2 as j2
from os.path import basename
from typing import Any
import common.config as cfg
import template_utils as tu
from urllib.parse import unquote
import json


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rk = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)
    bd = sconf.to_build_data(cfg.strip_full_refkey(rk), bk)
    map_params = [*bd.build.include.mappability]

    lowmap_path = cfg.smk_to_input_name(smk, "lowmap")
    notlowmap_path = cfg.smk_to_input_name(smk, "notinlowmap")

    map_env_path = cfg.smk_to_input_name(smk, "map_env")
    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")

    out = cfg.smk_to_output(smk)

    with open(lowmap_path, "r") as f:
        u = json.load(f)
        all_lowmap: str = basename(u["all_lowmap"])
        single_lowmap: list[str] = [*map(basename, u["single_lowmap"])]

    map_deps = tu.env_dependencies(map_env_path, {"bedops"})
    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools"})

    def render_description(t: j2.Template) -> str:
        return t.render(
            single_lowmap_files=single_lowmap,
            all_lowmap_file=all_lowmap,
            not_all_lowmap_file=basename(notlowmap_path),
            params=map_params,
        )

    def render_methods(t: j2.Template) -> str:
        return t.render(
            gemurl=unquote(sconf.tools.gemlib),
            params=map_params,
            deps={**map_deps, **bedtools_deps},
        )

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        "regions that are difficult to map for short reads",
        cfg.CoreLevel.MAPPABILITY,
    )

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
