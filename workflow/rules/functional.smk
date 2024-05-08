from common.config import CoreLevel, si_to_cds, si_to_mhc, si_to_kir, si_to_vdj

func = config.to_bed_dirs(CoreLevel.FUNCTIONAL)

all_region_types = ["cds", "vdj", "kir", "mhc"]


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
    rk = wildcards.ref_key
    bk = wildcards.build_key
    bd = config.to_build_data(rk, bk)

    def go(test, name, f):
        if test:
            rs = getattr(rules, name).output
            return expand(rs, ref_src_key=config.refkey_to_bed_refsrckeys(f, rk))
        else:
            return []

    return {
        "cds": go(bd.want_cds_src, "download_cds", si_to_cds),
        "mhc": go(bd.want_mhc_src, "download_mhc", si_to_mhc),
        "kir": go(bd.want_kir_src, "download_kir", si_to_kir),
        "vdj": go(bd.want_vdj_src, "download_vdj", si_to_vdj),
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
