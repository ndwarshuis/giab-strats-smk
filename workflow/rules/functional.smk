from common.config import CoreLevel, si_to_functional

func = config.to_bed_dirs(CoreLevel.FUNCTIONAL)


# use rule download_ref as download_ftbl with:
#     output:
#         func.src.data / "ftbl.txt.gz",
#     params:
#         src=lambda w: config.refsrckey_to_functional_src(si_to_ftbl, w.ref_src_key),
#     localrule: True
#     log:
#         func.src.log / "ftbl.log",


use rule download_ref as download_functional with:
    output:
        func.src.data / "gff.txt.gz",
    params:
        src=lambda w: to_bed_src(si_to_functional, w),
    localrule: True
    log:
        func.src.log / "gff.log",


checkpoint normalize_cds:
    input:
        lambda w: bed_src_inputs(rules.download_functional.output, si_to_functional, w),
        # unpack(
        #     lambda w: {
        #         k: expand(
        #             p,
        #             ref_src_key=config.refkey_to_functional_refsrckeys(f, w.ref_key),
        #         )
        #         for k, f, p in zip(
        #             ["ftbl", "gff"],
        #             [si_to_ftbl, si_to_gff],
        #             [rules.download_ftbl.output, rules.download_gff.output],
        #         )
        #     }
        # ),
    output:
        **{k: func.inter.filtersort.data / f"{k}.json" for k in ["cds", "vdj"]},
    params:
        cds_output=lambda w: to_output_pattern(func, "cds", w),
        vdj_output=lambda w: to_output_pattern(func, "vdj", w),
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


rule merge_functional:
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


use rule _invert_autosomal_regions as invert_functional with:
    input:
        rules.merge_functional.output,
    output:
        func.final("notinrefseq_cds"),


rule all_functional:
    input:
        rules.merge_functional.output,
        rules.invert_functional.output,
    localrule: True
