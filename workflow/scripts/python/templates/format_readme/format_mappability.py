import jinja2 as j2
from pathlib import Path
from os.path import basename
from typing import Any, Callable
import common.config as cfg
from common.functional import DesignError, filter_dict_strict
from urllib.parse import unquote
import json
import yaml


def env_dependencies(p: Path, deps: set[str]) -> dict[str, str]:
    with open(p, "r") as f:
        x = yaml.safe_load(f)
        try:
            ds: list[str] = x["dependencies"]
            ys = {s[0].strip(): s[1].strip() for d in ds if len(s := d.split("=")) == 2}
            return filter_dict_strict(ys, deps)
        except KeyError:
            raise DesignError(f"dependencies not found in {p}")


def load_template_path(p: Path) -> j2.Template:
    env = j2.Environment(
        loader=j2.FileSystemLoader(p.parent),
        undefined=j2.StrictUndefined,
    )
    return env.get_template(p.name)


def render_readme(
    smk: Any,
    render_description: Callable[[j2.Template], str],
    render_methods: Callable[[j2.Template], str],
    desc: str,
    level: cfg.CoreLevel,
) -> str:
    common_path = cfg.smk_to_input_name(smk, "common")
    description_path = cfg.smk_to_input_name(smk, "description")
    methods_path = cfg.smk_to_input_name(smk, "methods")

    desc_template = load_template_path(description_path)
    methods_template = load_template_path(methods_path)

    return load_template_path(common_path).render(
        description_text=render_description(desc_template),
        methods_text=render_methods(methods_template),
        group=level.value,
        description=desc,
    )


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

    map_deps = env_dependencies(map_env_path, {"bedops"})
    bedtools_deps = env_dependencies(bedtools_env_path, {"bedtools"})

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

    txt = render_readme(
        smk,
        render_description,
        render_methods,
        "regions that are difficult to map for short reads",
        cfg.CoreLevel.MAPPABILITY,
    )

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
