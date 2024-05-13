import yaml
import jinja2 as j2
from typing import Any, Callable
from pathlib import Path
from common.functional import DesignError, filter_dict_strict
import common.config as cfg


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
    hs: list[cfg.Haplotype],
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
        haplotypes=[h.value for h in hs],
    )
