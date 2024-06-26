from pathlib import Path
import jinja2 as j2
from typing import Any
import common.config as cfg
from common.functional import DesignError, fmap_maybe_def, fmap_maybe
import template_utils as tu


# TODO smells quite general ;)
def concat_comma(xs: list[str]) -> str:
    if len(xs) == 0:
        return ""
    elif len(xs) == 1:
        return xs[0]
    elif len(xs) == 2:
        return f"{xs[0]} and {xs[1]}"
    else:
        first = ", ".join(xs[:-1])
        return f"{first}, and {xs[-1]}"


def relative_path(rk: cfg.RefKeyFullS, p: Path) -> str:
    return str(Path("..") / p.parent.name / p.name)


def from_all_diff(
    rk: cfg.RefKeyFullS,
    a: cfg.AllDifficultPaths,
) -> tuple[str, str, str, list[str]]:
    gc_txt = "high/low GC regions" if a.gc_input else None
    repeat_txt = "tandem repeats" if a.repeat_input else None
    xy_txt = "difficult XY regions" if a.xy_inputs else None

    src_txt = concat_comma([x for x in [gc_txt, repeat_txt, xy_txt] if x is not None])

    return (
        a.output.positive.name,
        a.output.negative.name,
        f"This contains the above regions plus {src_txt}.",
        [relative_path(rk, p) for p in a._all_inputs],
    )


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)

    paths = smk.params["paths"]

    if not isinstance(paths, cfg.UnionPaths):
        raise DesignError()

    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")

    segdup_out = paths.segdup_lowmap.output
    empty: tuple[str | None, str | None, str | None, list[str]] = (
        None,
        None,
        None,
        [],
    )
    alldiff_pos, alldiff_neg, alldiff_desc, alldiff_sources = fmap_maybe_def(
        empty,
        lambda x: from_all_diff(rfk, x),
        paths.all_difficult,
    )

    def fmt_name(p: Path) -> str:
        return tu.sub_rk(rfk, p.name)

    def render_description(t: j2.Template) -> str:
        return t.render(
            segdup_map_file=fmt_name(segdup_out.positive),
            not_segdup_map_file=fmt_name(segdup_out.negative),
            alldifficult_file=fmap_maybe(lambda x: tu.sub_rk(rfk, x), alldiff_pos),
            not_alldifficult_file=fmap_maybe(lambda x: tu.sub_rk(rfk, x), alldiff_neg),
            alldifficult_desc=alldiff_desc,
        )

    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools", "samtools"})

    def render_methods(t: j2.Template) -> str:
        return t.render(
            map_file=relative_path(
                rfk,
                paths.segdup_lowmap.lowmap_input,
            ),
            segdup_file=relative_path(
                rfk,
                paths.segdup_lowmap.segdup_input.all_segdups.positive,
            ),
            all_difficult_files=alldiff_sources,
            deps=bedtools_deps,
        )

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        "combinations of other difficult stratifications",
        cfg.CoreLevel.UNION.value,
        sconf.refkey_haplotypes(rfk),
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
