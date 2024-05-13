import jinja2 as j2
from typing import Any, assert_never
import common.config as cfg
from common.functional import DesignError
from common.bed import BedLines
import template_utils as tu


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)
    rk = cfg.strip_full_refkey(rfk)
    bk = cfg.wc_to_buildkey(ws)
    bd = sconf.to_build_data(rk, bk)

    cds_src = bd.refdata.strat_inputs.functional.cds

    if cds_src is None:
        raise DesignError()

    params = cds_src.params
    cds_params = cds_src.params

    # def format_src(src: cfg.BedFileSrc | BedLines) -> str:
    #     if isinstance(src, cfg.HapChrFileSrc):
    #         src.src
    #     elif isinstance(src, cfg.HapChrTxtSrc):
    #         pass
    #     else:
    #         assert_never(src)

    bed = sconf.with_build_data_and_bed(
        rk,
        bk,
        cfg.bd_to_cds,
        lambda _, bf: bf.bed.documentation.as_list,
        lambda _, bf: bf.bed.documentation.as_list,
        lambda _, bf: bf.bed.documentation.as_list,
        lambda _, bf: bf.bed.documentation.as_list,
        lambda _, bf: bf.bed.documentation.as_list,
    )

    # if isinstance(bed, cfg.Single):
    #     pass
    # elif isinstance(bed, cfg.Double):
    #     pass
    # else:
    #     assert_never(bed)

    cds_path = cfg.smk_to_input_name(smk, "cds")
    notcds_path = cfg.smk_to_input_name(smk, "notcds")

    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")

    out = cfg.smk_to_output(smk)

    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools"})

    def render_description(t: j2.Template) -> str:
        return t.render(
            cds_file=cds_path.name,
            notcds_file=notcds_path.name,
        )

    def render_methods(t: j2.Template) -> str:
        return t.render(
            deps=bedtools_deps,
        )

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        "coding regions",
        cfg.CoreLevel.FUNCTIONAL,
        sconf.refkey_haplotypes(rfk),
    )

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
