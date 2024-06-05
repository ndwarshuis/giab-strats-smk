import jinja2 as j2
from pathlib import Path
from typing import Any
import common.config as cfg
from common.functional import DesignError
import template_utils as tu


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)

    paths = smk.params["paths"]

    if not isinstance(paths, cfg.SegdupPaths):
        raise DesignError()

    src_txt = sconf.with_build_data_and_bed_doc(
        rfk,
        bk,
        paths.superdups,
        cfg.bd_to_superdups,
        "The superdups file",
        4,
    )

    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")

    def fmt_name(p: Path) -> str:
        return tu.sub_rk(rfk, p.name)

    def render_description(t: j2.Template) -> str:
        return t.render(
            segdups_file=fmt_name(paths.all_segdups.positive),
            not_segdups_file=fmt_name(paths.all_segdups.negative),
            long_segdups_file=fmt_name(paths.long_segdups.positive),
            not_long_segdups_file=fmt_name(paths.long_segdups.negative),
        )

    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools", "samtools"})

    def render_methods(t: j2.Template) -> str:
        return t.render(src_txt=src_txt, deps=bedtools_deps)

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        "segmental duplications",
        cfg.CoreLevel.SEGDUPS,
        sconf.refkey_haplotypes(rfk),
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
