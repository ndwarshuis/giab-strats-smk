from pathlib import Path
import jinja2 as j2
from typing import Any
from urllib.parse import unquote
import common.config as cfg
from common.functional import DesignError, fmap_maybe, fmap_maybe_def
import template_utils as tu


def format_sources(
    sconf: cfg.GiabStrats,
    ps: cfg.SatellitesPaths | None,
    rfk: cfg.RefKeyFullS,
    bk: cfg.BuildKey,
) -> tuple[str | None, str | None, str | None]:
    def go(f: cfg.BuildDataToBed, paths: cfg.Path1or2) -> str:
        return sconf.with_build_data_and_bed_doc(rfk, bk, paths, f, None, 5)

    if ps is None:
        return (None, None, None)

    bd = sconf.to_build_data_full(rfk, bk)
    lc = bd.refdata.strat_inputs.low_complexity

    # fmap_maybe doesn't work here for some reason... :(
    # if lc.rmsk is None:
    #     rmsk_col = None
    # else:
    #     rmsk_col = lc.rmsk.class_col + 1

    overlap_txt = "Overlapping regions were then merged with `mergeBed`."
    if ps.used_censat:
        if lc.satellites is None:
            raise DesignError()

        sat_src_txt = go(cfg.bd_to_satellites, ps.sat_src)
        sat_other_txt = (
            " ".join(
                [
                    "This bed file was filtered for for values in column",
                    f"column {lc.satellites.sat_col + 1} which started with `ct`.",
                    overlap_txt,
                ]
            )
            if isinstance(lc.satellites, cfg.SatFile)
            else None
        )

    else:
        if lc.rmsk is None:
            raise DesignError()

        sat_src_txt = go(cfg.bd_to_rmsk, ps.sat_src)
        sat_other_txt = (
            " ".join(
                [
                    "This bed file was then filtered for the `Satellite` class in",
                    f"column {lc.rmsk.class_col + 1}.",
                    overlap_txt,
                ]
            )
            if isinstance(lc.rmsk, cfg.RMSKFile)
            else None
        )

    empty: list[str] = []

    sat_txt = "\n\n".join(
        [
            sat_src_txt,
            *fmap_maybe_def(empty, lambda x: [cfg.readme_fill(x)], sat_other_txt),
        ]
    )

    if ps.all_repeats is None:
        return (sat_txt, None, None)
    else:
        trf_txt = go(cfg.bd_to_simreps, ps.all_repeats.trf_src)
        if lc.rmsk is None:
            raise DesignError()
        if ps.used_censat:
            rmsk_src_txt = go(cfg.bd_to_rmsk, ps.all_repeats.rmsk_src)
            rmsk_other_txt = (
                cfg.readme_fill(
                    "This bed file was then filtered for `Low_complexity` "
                    f"and `Simple_Repeat` classes in column {lc.rmsk.class_col + 1}."
                )
                if isinstance(lc.rmsk, cfg.RMSKFile)
                else None
            )
            return (
                sat_txt,
                "\n\n".join(
                    [x for x in [rmsk_src_txt, rmsk_other_txt] if x is not None]
                ),
                trf_txt,
            )
        else:
            rmsk_txt = (
                cfg.readme_fill(
                    "The same repeat masker source from above was used here, "
                    "except that `Low_complexity` and `Simple_Repeat` classes were "
                    f"selected via column {lc.rmsk.class_col + 1}."
                )
                if isinstance(lc.rmsk, cfg.RMSKFile)
                else None
            )
            return (sat_txt, rmsk_txt, trf_txt)


def basenames(rk: cfg.RefKeyFullS, ps: list[Path]) -> list[str]:
    return [tu.sub_rk(rk, x.name) for x in ps]


def maybe_filename(p: Path | None) -> str | None:
    return fmap_maybe(lambda x: x.name, p)


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)

    paths = smk.params["paths"]

    if not isinstance(paths, cfg.LowComplexityPaths):
        raise DesignError()

    sat_paths = paths.satellites
    all_reps = fmap_maybe(lambda s: s.all_repeats, sat_paths)

    sat_src, rmsk_src, trf_src = format_sources(sconf, sat_paths, rfk, bk)

    uniform = paths.uniform_repeats

    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")

    _all_paths: list[tuple[str, Path | list[Path] | None]] = [
        # uniform repeat section
        ("perfect_files", uniform.perfect),
        ("imperfect_files", uniform.imperfect),
        ("homopolymers_file", uniform.homopolymers.positive),
        ("not_homopolymers_file", uniform.homopolymers.negative),
        # sat section
        ("sat_file", fmap_maybe(lambda x: x.sats.positive, sat_paths)),
        ("not_sat_file", fmap_maybe(lambda x: x.sats.negative, sat_paths)),
        # tr section
        ("all_filtered_tr_files", fmap_maybe(lambda x: x.filtered_trs, all_reps)),
        ("all_tr_file", fmap_maybe(lambda x: x.all_trs.positive, all_reps)),
        ("not_all_tr_file", fmap_maybe(lambda x: x.all_trs.negative, all_reps)),
        ("all_repeats_file", fmap_maybe(lambda x: x.all_repeats.positive, all_reps)),
        (
            "not_all_repeats_file",
            fmap_maybe(lambda x: x.all_repeats.negative, all_reps),
        ),
    ]

    all_paths = {
        k: (
            None
            if v is None
            else (basenames(rfk, v) if isinstance(v, list) else tu.sub_rk(rfk, v.name))
        )
        for k, v in _all_paths
    }

    def render_description(t: j2.Template) -> str:
        return t.render(**all_paths)

    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools", "samtools"})

    def render_methods(t: j2.Template) -> str:
        return t.render(
            **all_paths,
            sat_src=sat_src,
            rmsk_src=rmsk_src,
            trf_src=trf_src,
            repseq_src=unquote(sconf.tools.repseq),
            deps=bedtools_deps,
        )

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        "homopolymers and tandem repeats",
        cfg.CoreLevel.LOWCOMPLEXITY.value,
        sconf.refkey_haplotypes(rfk),
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
