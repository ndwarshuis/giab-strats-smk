import jinja2 as j2
from typing import Any
import common.config as cfg
from common.functional import DesignError, fmap_maybe_def, fmap_maybe
import template_utils as tu


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)

    paths = smk.params["paths"]

    if not isinstance(paths, cfg.MiscPaths):
        raise DesignError()

    lk = paths.desc.key

    # TODO not DRY
    def format_other(k: cfg.OtherStratKey) -> str:
        d = sconf.with_build_data_full(
            rfk,
            bk,
            lambda bd: fmap_maybe(lambda x: x.description, cfg.bd_to_other(lk, k, bd)),
            lambda bd: fmap_maybe(lambda x: x.description, cfg.bd_to_other(lk, k, bd)),
            lambda _, bd: fmap_maybe(
                lambda x: x.description, cfg.bd_to_other(lk, k, bd)
            ),
        )
        return fmap_maybe_def("No description", lambda x: tu.sub_rk(rfk, x), d)

    desc = {k: (format_other(k), v.output) for k, v in paths.paths.items()}
    src = {
        k: sconf.with_build_data_and_bed_doc(
            rfk,
            bk,
            cfg.map_single_or_double(
                lambda x: cfg.sub_wildcards_path(
                    x,
                    {"build_key": bk, "other_level_key": lk},
                ),
                p.source,
            ),
            lambda bd: cfg.bd_to_other(lk, k, bd),
            None,
            5,
        )
        for k, p in paths.paths.items()
    }

    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")
    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools", "samtools"})

    def render_description(t: j2.Template) -> str:
        return t.render(desc=desc)

    def render_methods(t: j2.Template) -> str:
        return t.render(
            src=src,
            deps=bedtools_deps,
        )

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        paths.desc.desc,
        paths.desc.key,
        sconf.refkey_haplotypes(rfk),
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
