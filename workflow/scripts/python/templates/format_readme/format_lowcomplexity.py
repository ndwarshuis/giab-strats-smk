import jinja2 as j2
from typing import Any
import common.config as cfg
from common.functional import DesignError, fmap_maybe
import template_utils as tu


def format_cds_params(p: cfg.CDSParams) -> str:
    source_txt = fmap_maybe(
        lambda x: f"column {x[1]+1} matched '{x[0]}'", p.source_match
    )
    type_txt = fmap_maybe(lambda x: f"column {x[1]+1} matched '{x[0]}'", p.type_match)

    return cfg.readme_fill(
        " ".join(
            [
                f"Lines where {source_txt} and {type_txt} were selected.",
                "Coordinates where start == end were removed.",
                "The remainining regions where merged using `mergeBed`.",
            ]
        )
    )


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)
    rk = cfg.strip_full_refkey(rfk)
    bk = cfg.wc_to_buildkey(ws)
    bd = sconf.to_build_data(rk, bk)

    cds_src = bd.refdata.strat_inputs.low_complexity

    inputs: cfg.LowComplexityPaths = smk.params["input_paths"]

    uniform = inputs.uniform_repeats
    sats = inputs.satellites

    if sats is None:
        pass
    else:
        if sats.all_repeats is None:
            pass
        else:
            pass

    src_txt = sconf.with_build_data_and_bed_doc(
        rfk,
        bk,
        cfg.smk_to_inputs_name(smk, "cds_inputs"),
        cfg.bd_to_cds,
        "The GFF file",
        None,
    )

    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")

    cds_path = cfg.smk_to_input_name(smk, "cds")
    notcds_path = cfg.smk_to_input_name(smk, "notcds")

    def render_description(t: j2.Template) -> str:
        return t.render(
            cds_file=cds_path.name,
            notcds_file=notcds_path.name,
        )

    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools", "samtools"})

    def render_methods(t: j2.Template) -> str:
        return t.render(
            # uniform repeat section
            uniform_repeat_paths=uniform,
            # perfect_files=[x.name for x in uniform.perfect],
            # imperfect_files=[x.name for x in uniform.imperfect],
            # homopolymers_file=uniform.homopolymers.name,
            # not_homopolymers_file=uniform.not_homopolymers.name,
            # sat section
            sat_paths=fmap_maybe(lambda i: i.satellites, inputs),
            sat_src=[],
            # sat_file=[],
            # not_sat_file=[],
            # tr section
            rmsk_src=[],
            sat_src=[],
            trf_src=[],
            repeat_paths=fmap_maybe(
                lambda i: fmap_maybe(lambda s: s.all_repeats, i.satellites), inputs
            ),
            # all_filtered_tr_files=[],
            # all_tr_file=[],
            # not_all_tr_file=[],
            # all_repeat_file=[],
            # not_all_repeat_file=[],
            # software section
            repseq_src="TODO",
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

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
