from pathlib import Path
import json
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


def read_paths(p: Path) -> list[Path]:
    with open(p, "r") as f:
        return [Path(p) for p in json.load(f)]


def format_satellite_source(
    sconf: cfg.GiabStrats,
    ps: cfg.SatellitesPaths,
    rfk: cfg.RefKeyFullS,
    bk: cfg.BuildKey,
) -> tuple[str, str | None, str | None]:
    overlap_txt = "Overlapping regions were then merged with `mergeBed`."

    src_paths = read_paths(ps.sat_src)

    bd = sconf.to_build_data_full(rfk, bk)

    if ps.used_censat:

        src_txt = sconf.with_build_data_and_bed_doc(
            rfk,
            bk,
            src_paths,
            cfg.bd_to_satellites,
            None,
            None,
        )

        sat_conf = bd.refdata.strat_inputs.low_complexity.satellites
        if sat_conf is None:
            raise DesignError()

        other_txt = (
            "This bed file was filtered for for values in column "
            f"{sat_conf.sat_col + 1} which started with `ct`."
        )
    else:

        src_txt = sconf.with_build_data_and_bed_doc(
            rfk,
            bk,
            src_paths,
            cfg.bd_to_rmsk,
            None,
            None,
        )

        rmsk_conf = bd.refdata.strat_inputs.low_complexity.rmsk
        if rmsk_conf is None:
            raise DesignError()

        other_txt = (
            "This bed file was then filtered for the `Satellite` class in "
            f"{rmsk_conf.class_col + 1}."
        )

    sat_txt = "\n\n".join(
        [cfg.readme_fill(x) for x in [src_txt, " ".join([other_txt, overlap_txt])]]
    )

    if ps.all_repeats is None:
        return (sat_txt, None, None)
    else:
        trf_txt = sconf.with_build_data_and_bed_doc(
            rfk,
            bk,
            read_paths(ps.all_repeats.trf_src),
            cfg.bd_to_simreps,
            None,
            None,
        )
        if ps.used_censat:
            # TODO not DRY
            rmsk_conf = bd.refdata.strat_inputs.low_complexity.rmsk
            if rmsk_conf is None:
                raise DesignError()
            rmsk_txt = sconf.with_build_data_and_bed_doc(
                rfk,
                bk,
                src_paths,
                cfg.bd_to_rmsk,
                None,
                None,
            )
            other_txt = (
                "This bed file was then filtered for `Low_complexity` and `Simple_Repeat` class in "
                f"{rmsk_conf.class_col + 1}."
            )
            # boooooooo0000000000
            return (sat_txt, rmsk_txt + " " + other_txt, trf_txt)
        else:
            rmsk_conf = bd.refdata.strat_inputs.low_complexity.rmsk
            if rmsk_conf is None:
                raise DesignError()
            rmsk_txt = (
                "The same repeat masker source from above was used here, "
                "except that `Low_complexity` and `Simple_Repeat` classes were "
                f"selected via {rmsk_conf.class_col + 1}."
            )
            return (sat_txt, rmsk_txt, trf_txt)


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)
    rk = cfg.strip_full_refkey(rfk)
    bk = cfg.wc_to_buildkey(ws)
    bd = sconf.to_build_data(rk, bk)

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
