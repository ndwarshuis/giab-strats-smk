import jinja2 as j2
from typing import Any
import common.config as cfg
import template_utils as tu
from urllib.parse import unquote
from common.functional import DesignError


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rk = cfg.wc_to_reffinalkey(ws)

    paths = smk.params["paths"]

    if not isinstance(paths, cfg.LowmapPaths):
        raise DesignError()

    def render_description(t: j2.Template) -> str:
        return t.render(
            single_lowmap_files=[tu.sub_rk(rk, p.name) for p in paths.single],
            all_lowmap_file=tu.sub_rk(rk, paths.union.positive.name),
            not_all_lowmap_file=tu.sub_rk(rk, paths.union.negative.name),
            params=paths.params,
        )

    map_env_path = cfg.smk_to_input_name(smk, "map_env")
    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")

    map_deps = tu.env_dependencies(map_env_path, {"bedops"})
    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools"})

    def render_methods(t: j2.Template) -> str:
        return t.render(
            gemurl=unquote(sconf.tools.gemlib),
            params=paths.params,
            deps={**map_deps, **bedtools_deps},
        )

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        "regions that are difficult to map for short reads",
        cfg.CoreLevel.MAPPABILITY,
        sconf.refkey_haplotypes(rk),
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
