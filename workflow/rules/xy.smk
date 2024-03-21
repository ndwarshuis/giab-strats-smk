from common.config import CoreLevel

xy = config.to_bed_dirs(CoreLevel.XY)


use rule download_gaps as download_genome_features_bed with:
    output:
        xy.src.data / "genome_features_{sex_chr}.bed.gz",
    log:
        xy.src.log / "genome_features_{sex_chr}.log",
    params:
        src=lambda w: (
            config.refsrckey_to_y_features_src
            if w.sex_chr == "Y"
            else config.refsrckey_to_x_features_src
        )(w.ref_src_key),
    localrule: True


rule write_PAR_final:
    input:
        bed=rules.write_PAR_intermediate.output,
        gapless=rules.get_gapless.output.parY,
        genome=rules.filter_sort_ref.output["genome"],
    output:
        xy.final("chr{sex_chr}_PAR"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        intersectBed -a {input.bed} -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


rule filter_XTR_features:
    input:
        bed=lambda w: expand_final_to_src(
            rules.download_genome_features_bed.output, w
        )[0],
        gapless=rules.get_gapless.output.parY,
        genome=rules.filter_sort_ref.output["genome"],
    output:
        xy.final("chr{sex_chr}_XTR"),
    log:
        xy.inter.postsort.log / "{sex_chr}_filter_XTR_features.txt",
    conda:
        "../envs/bedtools.yml"
    params:
        level="XTR",
    script:
        "../scripts/python/bedtools/xy/filter_sort_features.py"


use rule filter_XTR_features as filter_ampliconic_features with:
    output:
        xy.final("chr{sex_chr}_ampliconic"),
    log:
        xy.inter.postsort.log / "{sex_chr}_filter_ampliconic_features.txt",
    params:
        level="Ampliconic",


# def par_input(wildcards):
#     test_fun = config.want_xy_y if wildcards.sex_chr == "Y" else config.want_xy_x
#     return (
#         rules.write_PAR_final.output
#         if test_fun(wildcards.ref_key, wildcards.build_key)
#         else rules.write_PAR_intermediate.output
#     )


rule invert_PAR:
    input:
        bed=rules.write_PAR_final.output,
        genome=rules.filter_sort_ref.output["genome"],
        gapless=rules.get_gapless.output.parY,
    output:
        xy.final("chr{sex_chr}_nonPAR"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        complementBed -i {input.bed} -g {input.genome} | \
        grep {wildcards.sex_chr} | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


rule filter_autosomes:
    input:
        rules.get_gapless.output.auto,
    output:
        xy.final("AllAutosomes"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        gunzip -c {input} | \
        grep -v \"X\|Y\" | \
        bgzip -c > {output}
        """


# helper functions for build targets; note that for each of these we need to
# first get a list of the X/Y chromosomes required (in the case of dip2, for
# the haplotype at hand)
def all_xy_features(ref_final_key, build_key):
    bd = config.to_build_data_full(ref_final_key, build_key)
    sex_chrs = config.buildkey_to_wanted_xy_names(ref_final_key, build_key)
    all_targets = [
        (rules.filter_XTR_features.output[0], bd.want_xy_XTR),
        (rules.filter_ampliconic_features.output[0], bd.want_xy_ampliconic),
    ]
    targets = [x for x, y in all_targets if y]
    return expand(
        targets,
        allow_missing=True,
        sex_chr=sex_chrs,
        ref_final_key=ref_final_key,
        build_key=build_key,
    )


def all_xy_PAR(ref_final_key, build_key):
    bd = config.to_build_data_full(ref_final_key, build_key)
    sex_chrs = [
        c
        for c in config.buildkey_to_wanted_xy_names(ref_final_key, build_key)
        if bd.want_xy_PAR(c)
    ]
    return expand(
        rules.invert_PAR.output + rules.write_PAR_final.output,
        allow_missing=True,
        sex_chr=sex_chrs,
        ref_final_key=ref_final_key,
        build_key=build_key,
    )


def all_xy_sex(ref_final_key, build_key):
    return all_xy_PAR(ref_final_key, build_key) + all_xy_features(
        ref_final_key, build_key
    )
