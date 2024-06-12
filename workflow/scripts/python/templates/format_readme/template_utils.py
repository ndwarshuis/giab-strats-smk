import yaml
import jinja2 as j2
from typing import Any, Callable
from pathlib import Path
from common.functional import DesignError, filter_dict_strict, from_maybe
import common.config as cfg


def env_dependencies(p: Path, deps: set[str]) -> dict[str, str]:
    with open(p, "r") as f:
        x = yaml.safe_load(f)
        try:
            ds: list[str] = x["dependencies"]
            ys = {
                s[0].strip(): s[1].strip()
                for d in ds
                if len(s := d.split("=")) in [2, 3]
            }
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
    level: str,
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
        group=level,
        description=desc,
        haplotypes=[h.name for h in hs],
    )


def sub_rk(rk: cfg.RefKeyFullS, n: str) -> str:
    return cfg.sub_wildcard(n, "ref_final_key", rk)


def fmt_other_srcs(
    sconf: cfg.GiabStrats,
    rk: cfg.RefKeyFullS,
    bk: cfg.BuildKey,
    lk: cfg.OtherLevelKey,
    paths: dict[cfg.OtherStratKey, cfg.SourceOutputPaths],
) -> dict[cfg.OtherStratKey, str]:
    # NOTE sub these three keys because unlike other sources, these are
    # build specific and are also denoted by the other strat/level keys. Other
    # sources are denoted by just the ref_src_key (which should already be
    # subbed).
    return {
        sk: sconf.with_build_data_and_bed_doc(
            rk,
            bk,
            p.source.map(
                lambda x: cfg.sub_wildcards_path(
                    x,
                    {
                        "build_key": bk,
                        "other_level_key": lk,
                        "other_strat_key": sk,
                    },
                ),
            ),
            lambda bd: cfg.bd_to_other_bed(lk, sk, bd),
            None,
            5,
        )
        for sk, p in paths.items()
    }


def fmt_other_descriptions(
    sconf: cfg.GiabStrats,
    rk: cfg.RefKeyFullS,
    bk: cfg.BuildKey,
    lk: cfg.OtherLevelKey,
    paths: dict[cfg.OtherStratKey, cfg.SourceOutputPaths],
) -> dict[str, str]:
    def go(sk: cfg.OtherStratKey) -> str:
        d = sconf.with_build_data_full(
            rk,
            bk,
            lambda bd: (
                b.description
                if (b := cfg.bd_to_other(lk, sk, bd)) is not None
                else None
            ),
            lambda bd: (
                b.description
                if (b := cfg.bd_to_other(lk, sk, bd)) is not None
                else None
            ),
            lambda _, bd: (
                b.description
                if (b := cfg.bd_to_other(lk, sk, bd)) is not None
                else None
            ),
        )
        return from_maybe("No description", d)

    return {sub_rk(rk, p.output.name): go(k) for k, p in paths.items()}
