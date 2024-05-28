import jinja2 as j2
from pathlib import Path
from typing import Any
import common.config as cfg
from common.functional import fmap_maybe, DesignError, from_maybe
import template_utils as tu

# TODO booooo put this in the main config somewhere
OTHERKEY = cfg.OtherLevelKey("OtherDifficult")


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    rfk = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)

    paths: cfg.OtherDifficultPaths = smk.params["paths"]

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
        return from_maybe("No description", d)

    def render_description(t: j2.Template) -> str:
        return t.render(
            gaps_file=fmap_maybe(lambda z: z.name, paths.gaps_output),
            vdj_file=fmap_maybe(lambda z: z.name, paths.vdj_output),
            kir_file=fmap_maybe(lambda z: z.name, paths.kir_output),
            mhc_file=fmap_maybe(lambda z: z.name, paths.mhc_output),
            other_files={p: format_other(k) for k, p in paths.other_outputs.items()},
        )

    bedtools_deps = tu.env_dependencies(bedtools_env_path, {"bedtools", "samtools"})

    gaps_src = (
        sconf.with_build_data_and_bed_doc(
            rfk,
            bk,
            paths.gaps_src,
            cfg.bd_to_gaps,
            "Gaps file",
            None,
        )
        if paths.gaps_src is not None
        else None
    )

    if paths.refseq_src is not None and (
        (paths.vdj_output is not None and paths.vdj_src is None)
        or (paths.mhc_output is not None and paths.mhc_src is None)
        or (paths.kir_output is not None and paths.kir_src is None)
    ):
        refseq_para = sconf.with_build_data_and_bed_doc(
            rfk,
            bk,
            paths.refseq_src,
            cfg.bd_to_cds,
            "Refseq GFF file",
            None,
        )
        if paths.vdj_output is not None and paths.vdj_src is None:
            vdj_para = (
                "VDJ regions were found by filtering the GFF for"
                "attributes matching '^ID=gene-(IGH|IGK|IGL|TRA|TRB|TRG);'."
            )
        else:
            vdj_para = None

        if paths.kir_output is not None and paths.kir_src is None:
            raise DesignError(
                "KIR regions were found in Mike Portnoy's last beer bottle."
            )
        else:
            kir_para = None

        if paths.mhc_output is not None and paths.mhc_src is None:
            raise DesignError(
                "MHC regions were found in Taylor Swift's secret torture chamber."
            )
        else:
            mhc_para = None

        refseq_src = "\n\n".join(
            [refseq_para]
            + [
                cfg.readme_fill(x)
                for x in [vdj_para, kir_para, mhc_para]
                if x is not None
            ]
        )
    else:
        refseq_src = None

    def immuno_src(
        src: cfg.Path1or2 | None,
        output: Path | None,
        f: cfg.BuildDataToBed,
        name: str,
    ) -> str | None:
        if src is not None and output is not None:
            return sconf.with_build_data_and_bed_doc(
                rfk, bk, src, f, f"{name} file", None
            )
        else:
            return None

    vdj_src = immuno_src(paths.vdj_src, paths.vdj_output, cfg.bd_to_vdj, "VDJ")
    kir_src = immuno_src(paths.kir_src, paths.kir_output, cfg.bd_to_kir, "KIR")
    mhc_src = immuno_src(paths.mhc_src, paths.mhc_output, cfg.bd_to_mhc, "MHC")

    other_srcs = {
        k: sconf.with_build_data_and_bed_doc(
            rfk,
            bk,
            paths,
            lambda bd: cfg.bd_to_other(OTHERKEY, k, bd),
            None,
            None,
        )
        for k, paths in paths.other_src.items()
    }

    def render_methods(t: j2.Template) -> str:
        return t.render(
            gaps_src=gaps_src,
            refseq_src=refseq_src,
            vdj_src=vdj_src,
            mhc_src=mhc_src,
            kir_src=kir_src,
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
