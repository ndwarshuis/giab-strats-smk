import jinja2 as j2
from typing import Any
import common.config as cfg
from common.functional import (
    fmap_maybe,
    DesignError,
    uncons_maybe,
    fmap_maybe_def,
)
import template_utils as tu

# TODO booooo put this in the main config somewhere
OTHERKEY = cfg.OtherLevelKey("OtherDifficult")


# TODO implementme
def kir_para() -> str:
    raise DesignError("KIR regions were found in Mike Portnoy's last beer bottle.")


# TODO implementme
def mhc_para() -> str:
    raise DesignError(
        "MHC regions were found in Taylor Swift's secret torture chamber."
    )


def vdj_para() -> str:
    return (
        "VDJ regions were found by filtering the GFF for "
        "attributes matching '^ID=gene-(IGH|IGK|IGL|TRA|TRB|TRG);'."
    )


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)

    paths = smk.params["paths"]

    if not isinstance(paths, cfg.OtherDifficultPaths):
        raise DesignError()

    bedtools_env_path = cfg.smk_to_input_name(smk, "bedtools_env")

    def format_other(k: cfg.OtherStratKey) -> str:
        d = sconf.with_build_data_full(
            rfk,
            bk,
            lambda bd: fmap_maybe(
                lambda x: x.description, cfg.bd_to_other(OTHERKEY, k, bd)
            ),
            lambda bd: fmap_maybe(
                lambda x: x.description, cfg.bd_to_other(OTHERKEY, k, bd)
            ),
            lambda _, bd: fmap_maybe(
                lambda x: x.description, cfg.bd_to_other(OTHERKEY, k, bd)
            ),
        )
        return fmap_maybe_def("No description", lambda x: tu.sub_rk(rfk, x), d)

    def fmt_maybe(p: cfg.GapsPaths | None) -> str | None:
        return fmap_maybe(lambda z: tu.sub_rk(rfk, z.output.name), p)

    def render_description(t: j2.Template) -> str:
        return t.render(
            gaps_file=fmt_maybe(paths.gaps),
            vdj_file=fmt_maybe(paths.vdj),
            kir_file=fmt_maybe(paths.kir),
            mhc_file=fmt_maybe(paths.mhc),
            other_files={p: format_other(k) for k, p in paths.other.items()},
        )

    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools", "samtools"})

    # gaps is relatively straightforward, just get the source for the gaps if
    # needed
    gaps_src = (
        sconf.with_build_data_and_bed_doc(
            rfk,
            bk,
            paths.gaps.source,
            cfg.bd_to_gaps,
            "Gaps file",
            4,
        )
        if paths.gaps is not None
        else None
    )

    # The immuno sources (KIR/MHC/VDJ) need special treatment because some of
    # them might be based on the refseq/cds bed file or they may have their
    # own file. Partition the list of these to those from refseq and those not,
    # then process independently. In the end, We will either have a separate
    # section for KIR/MHC/VDJ or a GFF section in place of some/all of these.
    immuno_pairs = [
        (s, n, f, g)
        for s, n, f, g in [
            (paths.kir, "KIR", kir_para, cfg.bd_to_kir),
            (paths.mhc, "MHC", mhc_para, cfg.bd_to_mhc),
            (paths.vdj, "VDJ", vdj_para, cfg.bd_to_vdj),
        ]
        if s is not None
    ]

    refseq = [(s, f) for s, _, f, _ in immuno_pairs if s.source_is_refseq]
    not_refseq = [(s, n, f) for s, n, _, f in immuno_pairs if not s.source_is_refseq]

    if (r := uncons_maybe(refseq)) is None:
        refseq_src = None
    else:
        s0, ss = r

        refseq_para = sconf.with_build_data_and_bed_doc(
            rfk,
            bk,
            s0[0].source,
            cfg.bd_to_cds,
            "Refseq GFF file",
            4,
        )

        refseq_src = "\n\n".join(
            [refseq_para]
            + [cfg.readme_fill(x) for x in [s[1]() for s in [s0, *ss]] if x is not None]
        )

    immuno_srcs = {
        n: sconf.with_build_data_and_bed_doc(rfk, bk, s.source, f, f"{n} file", 4)
        for s, n, f in not_refseq
    }

    other_srcs = {
        k: sconf.with_build_data_and_bed_doc(
            rfk,
            bk,
            p.source,
            lambda bd: cfg.bd_to_other(OTHERKEY, k, bd),
            None,
            5,
        )
        for k, p in paths.other.items()
    }

    def render_methods(t: j2.Template) -> str:
        return t.render(
            gaps_src=gaps_src,
            refseq_src=refseq_src,
            immuno_srcs=immuno_srcs,
            other_srcs=other_srcs,
            deps=bedtools_deps,
        )

    txt = tu.render_readme(
        smk,
        render_description,
        render_methods,
        "miscellaneously difficult regions",
        cfg.CoreLevel.OTHER_DIFFICULT,
        sconf.refkey_haplotypes(rfk),
    )

    out = cfg.smk_to_output(smk)

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
