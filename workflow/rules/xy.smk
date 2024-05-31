from common.config import CoreLevel, MutualPathPair

xy = config.to_bed_dirs(CoreLevel.XY)


def all_xy(ref_final_key, build_key):
    return config.all_xy(
        ref_final_key,
        build_key,
        Path(rules.download_genome_features_bed.output[0]),
        Path(rules.filter_XTR_features.output[0]),
        Path(rules.filter_ampliconic_features.output[0]),
        MutualPathPair(
            Path(rules.write_PAR_final.output[0]),
            Path(rules.invert_PAR.output[0]),
        ),
        Path(rules.filter_autosomes.output[0]),
    )


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


# # helper functions for build targets; note that for each of these we need to
# # first get a list of the X/Y chromosomes required (in the case of dip2, for
# # the haplotype at hand)
# def all_xy_features(ref_final_key, build_key):
#     bd = config.to_build_data_full(ref_final_key, build_key)
#     sex_chrs = config.buildkey_to_wanted_xy_names(ref_final_key, build_key)

#     def go(target, test):
#         if test:
#             return expand(
#                 target,
#                 allow_missing=True,
#                 sex_chr=sex_chrs,
#                 ref_final_key=ref_final_key,
#                 build_key=build_key,
#             )
#         else:
#             return []

#     return {
#         "xtr": go(rules.filter_XTR_features.output, bd.have_xy_XTR),
#         "ampliconic": (
#             go(rules.filter_ampliconic_features.output, bd.have_xy_ampliconic)
#         ),
#     }


# # TODO wet
# def all_xy_features_inputs(wildcards):
#     rk = wildcards["ref_final_key"]
#     bk = wildcards["build_key"]
#     bd = config.to_build_data_full(rk, bk)
#     sex_chrs = config.buildkey_to_wanted_xy_names(rk, bk)
#     have_features = bd.have_xy_XTR and bd.have_xy_ampliconic

#     return {
#         f"{c}_features_input": (
#             expand(
#                 rules.filter_XTR_features.input.bed(wildcards),
#                 allow_missing=True,
#                 sex_chr=c,
#                 build_key=build_key,
#             )[0]
#             if (c in sex_chrs and test and have_features)
#             else None
#         )
#         for c in ["X", "Y"]
#     }


# def all_xy_features_flat(ref_final_key, build_key):
#     return [y for x in all_xy_features(ref_final_key, build_key).values() for y in x]


# def _all_xy_PAR_from_rule(rulename, ref_final_key, build_key):
#     bd = config.to_build_data_full(ref_final_key, build_key)
#     sex_chrs = [
#         c
#         for c in config.buildkey_to_wanted_xy_names(ref_final_key, build_key)
#         if bd.have_xy_PAR(c)
#     ]
#     return expand(
#         rulename,
#         allow_missing=True,
#         sex_chr=sex_chrs,
#         ref_final_key=ref_final_key,
#         build_key=build_key,
#     )


# def all_xy_PAR(ref_final_key, build_key):
#     return _all_xy_PAR_from_rule(
#         rules.write_PAR_final.output,
#         ref_final_key,
#         build_key,
#     )


# def all_xy_nonPAR(ref_final_key, build_key):
#     return _all_xy_PAR_from_rule(
#         rules.invert_PAR.output,
#         ref_final_key,
#         build_key,
#     )


# def all_xy_sex(ref_final_key, build_key):
#     return (
#         all_xy_PAR(ref_final_key, build_key)
#         + all_xy_nonPAR(ref_final_key, build_key)
#         + all_xy_features_flat(ref_final_key, build_key)
#     )


# def all_xy_flat(ref_final_key, build_key):
#     return all_xy(ref_final_key, build_key).all_outputs


rule xy_readme:
    input:
        common="workflow/templates/common.j2",
        description="workflow/templates/xy_description.j2",
        methods="workflow/templates/xy_methods.j2",
        bedtools_env="workflow/envs/bedtools.yml",
        _sources=lambda w: all_xy(w["ref_final_key"], w["build_key"]).all_sources,
    params:
        paths=lambda w: all_xy(w["ref_final_key"], w["build_key"]),
    output:
        xy.readme,
    conda:
        "../envs/templates.yml"
    script:
        "../scripts/python/templates/format_readme/format_xy.py"
