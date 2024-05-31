from common.config import CoreLevel, si_to_superdups, strip_full_refkey, MutualPathPair

segdup = config.to_bed_dirs(CoreLevel.SEGDUPS)


def all_segdups_sources(ref_key, build_key):
    return config.all_segdups_sources(
        ref_key,
        build_key,
        Path(rules.download_superdups.output[0]),
    )


def all_segdups(ref_final_key, build_key):
    return config.all_segdups(
        ref_final_key,
        build_key,
        all_segdups_sources(strip_full_refkey(ref_final_key), build_key),
        MutualPathPair(
            Path(rules.merge_superdups.output[0]),
            Path(rules.notin_superdups.output[0]),
        ),
        MutualPathPair(
            Path(rules.filter_long_superdups.output[0]),
            Path(rules.notin_long_superdups.output[0]),
        ),
    )


use rule download_gaps as download_superdups with:
    output:
        segdup.src.data / "superdups.txt.gz",
    log:
        segdup.src.log / "superdups.log",
    params:
        src=lambda w: config.refsrckey_to_bed_src(si_to_superdups, w.ref_src_key),
    localrule: True


checkpoint normalize_superdups:
    input:
        lambda w: all_segdups_sources(w["ref_key"], w["build_key"]).all_sources,
    output:
        segdup.inter.filtersort.data / "segdups.json",
    params:
        output_pattern=lambda w: to_output_pattern(segdup, "segdups", w),
    conda:
        "../envs/bedtools.yml"
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_key, w.build_key, lambda m: m.normalizeSuperdups
        ),
    benchmark:
        segdup.inter.filtersort.bench / "normalize_segdups.txt"
    script:
        "../scripts/python/bedtools/segdups/normalize_superdups.py"


rule merge_superdups:
    input:
        bed=lambda w: read_checkpoint("normalize_superdups", w),
        genome=rules.filter_sort_ref.output["genome"],
        gapless=rules.get_gapless.output.auto,
    output:
        segdup.final("segdups"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        mergeBed -i {input.bed} -d 100 | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


rule filter_long_superdups:
    input:
        rules.merge_superdups.output,
    output:
        segdup.final("segdups_gt10kb"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        gunzip -c {input} | \
        awk '($3-$2 > 10000)' | \
        bgzip -c > {output}
        """


use rule _invert_autosomal_regions as notin_superdups with:
    input:
        rules.merge_superdups.output,
    output:
        segdup.final("notinsegdups"),


use rule notin_superdups as notin_long_superdups with:
    input:
        rules.filter_long_superdups.output,
    output:
        segdup.final("notinsegdups_gt10kb"),


# rule all_segdups:
#     input:
#         rules.merge_superdups.output,
#         rules.filter_long_superdups.output,
#         rules.notin_superdups.output,
#         rules.notin_long_superdups.output,
#     localrule: True


rule segdups_readme:
    input:
        common="workflow/templates/common.j2",
        description="workflow/templates/segdups_description.j2",
        methods="workflow/templates/segdups_methods.j2",
        _sources=lambda w: all_segdups(w["ref_final_key"], w["build_key"]).all_sources,
        bedtools_env="workflow/envs/bedtools.yml",
    output:
        segdup.readme,
    params:
        paths=lambda w: all_segdups(w["ref_final_key"], w["build_key"]),
    conda:
        "../envs/templates.yml"
    script:
        "../scripts/python/templates/format_readme/format_segdups.py"
