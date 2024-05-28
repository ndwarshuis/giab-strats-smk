from common.config import (
    CoreLevel,
    si_to_cds,
    si_to_mhc,
    si_to_kir,
    si_to_vdj,
    strip_full_refkey,
)

# NOTE this contains rules that generate stratifications in Functional,
# OtherDifficult, and others without explicit rules (Ancestry et al)

func = config.to_bed_dirs(CoreLevel.FUNCTIONAL)
odiff = config.to_bed_dirs(CoreLevel.OTHER_DIFFICULT)

# TODO don't hardcode this...
all_region_types = ["cds", "vdj", "kir", "mhc"]


def all_otherdifficult_sources_smk(wildcards):
    return config.all_otherdifficult_sources(
        wildcards["ref_key"],
        wildcards["build_key"],
        Path(rules.download_gaps.output[0]),
        Path(rules.download_cds.output[0]),
        Path(rules.download_kir.output[0]),
        Path(rules.download_mhc.output[0]),
        Path(rules.download_vdj.output[0]),
        {},
    )


def all_otherdifficult_paths_smk(ref_final_key, build_key):
    # TODO actually make other work
    return config.all_otherdifficult_paths(
        ref_final_key,
        build_key,
        Path(rules.download_gaps.output[0]),
        Path(rules.download_cds.output[0]),
        Path(rules.download_kir.output[0]),
        Path(rules.download_mhc.output[0]),
        Path(rules.download_vdj.output[0]),
        {},
        Path(rules.get_gaps.output[0]),
        Path(rules.merge_cds.output[0]),
        Path(rules.invert_cds.output[0]),
        Path(rules.remove_kir_gaps.output[0]),
        Path(rules.remove_mhc_gaps.output[0]),
        Path(rules.remove_vdj_gaps.output[0]),
        {},
    )


use rule download_gaps as download_cds with:
    output:
        func.src.data / "gff.txt.gz",
    params:
        src=lambda w: to_bed_src(si_to_cds, w),
    localrule: True
    log:
        func.src.log / "gff.log",


use rule download_gaps as download_vdj with:
    output:
        func.src.data / "vdj.bed.gz",
    params:
        src=lambda w: to_bed_src(si_to_vdj, w),
    localrule: True
    log:
        func.src.log / "vdj.log",


use rule download_gaps as download_mhc with:
    output:
        func.src.data / "mhc.bed.gz",
    params:
        src=lambda w: to_bed_src(si_to_mhc, w),
    localrule: True
    log:
        func.src.log / "mhc.log",


use rule download_gaps as download_kir with:
    output:
        func.src.data / "kir.bed.gz",
    params:
        src=lambda w: to_bed_src(si_to_kir, w),
    localrule: True
    log:
        func.src.log / "kir.log",


def functional_inputs(wildcards):
    src = all_otherdifficult_sources_smk(wildcards)
    return {
        "cds": src.refseq_paths,
        "mhc": src.mhc_paths,
        "kir": src.kir_paths,
        "vdj": src.vdj_paths,
    }


checkpoint normalize_cds:
    input:
        unpack(functional_inputs),
    output:
        **{k: func.inter.filtersort.data / f"{k}.json" for k in all_region_types},
    params:
        cds=lambda w: to_output_pattern(func, "cds", w),
        vdj=lambda w: to_output_pattern(func, "vdj", w),
        mhc=lambda w: to_output_pattern(func, "mhc", w),
        kir=lambda w: to_output_pattern(func, "kir", w),
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_key, w.build_key, lambda m: m.normalizeCds
        ),
    benchmark:
        func.inter.filtersort.bench / "normalize_cds.txt"
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/functional/normalize_cds.py"


rule merge_cds:
    input:
        bed=lambda w: read_named_checkpoint("normalize_cds", "cds", w),
        genome=rules.filter_sort_ref.output["genome"],
        gapless=rules.get_gapless.output.auto,
    output:
        func.final("refseq_cds"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        mergeBed -i {input.bed} | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


use rule _invert_autosomal_regions as invert_cds with:
    input:
        rules.merge_cds.output,
    output:
        func.final("notinrefseq_cds"),


rule all_cds:
    input:
        rules.merge_cds.output,
        rules.invert_cds.output,
    localrule: True


rule get_gaps:
    input:
        gapless=rules.get_gapless.output.auto,
        genome=rules.filter_sort_ref.output["genome"],
    output:
        odiff.final("gaps_slop15kb"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        complementBed -i {input.gapless} -g {input.genome} | \
        slopBed -i stdin -b 15000 -g {input.genome} | \
        mergeBed -i stdin | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


rule remove_vdj_gaps:
    input:
        bed=lambda w: read_named_checkpoint("normalize_cds", "vdj", w),
        genome=rules.filter_sort_ref.output["genome"],
        gapless=rules.get_gapless.output.auto,
    output:
        odiff.final("VDJ"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        mergeBed -i {input.bed} | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


use rule remove_vdj_gaps as remove_mhc_gaps with:
    input:
        bed=lambda w: read_named_checkpoint("normalize_cds", "mhc", w),
        genome=rules.filter_sort_ref.output["genome"],
        gapless=rules.get_gapless.output.auto,
    output:
        odiff.final("MHC"),


use rule remove_vdj_gaps as remove_kir_gaps with:
    input:
        bed=lambda w: read_named_checkpoint("normalize_cds", "kir", w),
        genome=rules.filter_sort_ref.output["genome"],
        gapless=rules.get_gapless.output.auto,
    output:
        odiff.final("KIR"),


rule functional_readme:
    input:
        common="workflow/templates/common.j2",
        description="workflow/templates/functional_description.j2",
        methods="workflow/templates/functional_methods.j2",
        _sources=lambda w: all_otherdifficult(
            w["ref_final_key"], w["build_key"]
        ).sources.all_functional_inputs,
        bedtools_env="workflow/envs/bedtools.yml",
    params:
        paths=lambda w: all_otherdifficult(w["ref_final_key"], w["build_key"]),
    output:
        func.readme,
    conda:
        "../envs/templates.yml"
    script:
        "../scripts/python/templates/format_readme/format_functional.py"


rule otherdifficult_readme:
    input:
        common="workflow/templates/common.j2",
        description="workflow/templates/otherdifficult_description.j2",
        methods="workflow/templates/otherdifficult_methods.j2",
        bedtools_env="workflow/envs/bedtools.yml",
        _sources=lambda w: all_otherdifficult(
            w["ref_final_key"], w["build_key"]
        ).sources.all_otherdifficult_inputs,
    params:
        paths=lambda w: all_otherdifficult(w["ref_final_key"], w["build_key"]),
    output:
        odiff.readme,
    conda:
        "../envs/templates.yml"
    script:
        "../scripts/python/templates/format_readme/format_otherdifficult.py"
