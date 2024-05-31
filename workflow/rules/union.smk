from common.config import CoreLevel, MutualPathPair, strip_full_refkey

uni = config.to_bed_dirs(CoreLevel.UNION)


def all_union(ref_final_key, build_key):
    _rk = strip_full_refkey(ref_final_key)
    return config.all_union(
        ref_final_key,
        build_key,
        all_segdups_sources(_rk, bk),
        Path(rules.merge_nonunique.output.all_lowmap),
        Path(rules.intersect_gc_ranges.output.widest_extreme),
        all_low_complexity(ref_final_key, build_key),
        all_xy(ref_final_key, build_key),
        MutualPathPair(
            Path(rules.intersect_segdup_and_map.output[0]),
            Path(rules.invert_segdup_and_map.output[0]),
        ),
        MutualPathPair(
            Path(rules.intersect_alldifficult.output[0]),
            Path(rules.invert_alldifficult.output[0]),
        ),
    )


# TODO this function is silly and wet
# def segdup_and_map_inputs(wildcards):
#     rk = wildcards["ref_final_key"]
#     bk = wildcards["build_key"]
#     bd = config.to_build_data_full(rk, bk)
#     s = {"segdups": rules.merge_superdups.output[0]} if bd.have_and_want_segdups else {}
#     l = (
#         {"lowmap": rules.merge_nonunique.output.all_lowmap}
#         if bd.have_and_want_mappability
#         else {}
#     )
#     return {**s, **l}


rule intersect_segdup_and_map:
    input:
        lambda w: all_union(w["ref_final_key"], w["build_key"]).segdup_lowmap_inputs,
    output:
        uni.final("alllowmapandsegdupregions"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        multiIntersectBed -i {input} | \
        mergeBed -i stdin | \
        bgzip -c \
        > {output}
        """


use rule _invert_autosomal_regions as invert_segdup_and_map with:
    input:
        rules.intersect_segdup_and_map.output,
    output:
        uni.final("notinalllowmapandsegdupregions"),


# TODO ...and so is this one
# def all_difficult_inputs(wildcards):
#     rk = wildcards["ref_final_key"]
#     bk = wildcards["build_key"]
#     bd = config.to_build_data_full(rk, bk)
#     gc = {"gc": rules.intersect_gc_ranges.output.widest_extreme} if bd.want_gc else {}
#     return {
#         **gc,
#         **{
#             "_segdups_lowmap": rules.intersect_segdup_and_map.output[0],
#             "repeats": rules.merge_HPs_and_TRs.output[0],
#             "xy": all_xy_features(rk, bk),
#         },
#     }


# TODO include other difficult in here as well
use rule intersect_segdup_and_map as intersect_alldifficult with:
    input:
        lambda w: all_union(w["ref_final_key"], w["build_key"]).all_difficult_inputs,
        # unpack(all_difficult_inputs),
        # _segdups_lowmap=rules.intersect_segdup_and_map.output[0],
        # repeats=rules.merge_HPs_and_TRs.output[0],
        # xy=lambda w: all_xy_features(w.ref_final_key, w.build_key),
        # gc=lambda w: gc_inputs(w.ref_final_key, w.build_key)["widest_extreme"],
    output:
        uni.final("alldifficultregions"),


use rule invert_segdup_and_map as invert_alldifficult with:
    input:
        bed=rules.intersect_alldifficult.output,
    output:
        uni.final("notinalldifficultregions"),


# rule all_segdup_and_map:
#     input:
#         rules.intersect_segdup_and_map.output,
#         rules.invert_segdup_and_map.output,
#     localrule: True


# rule all_alldifficult:
#     input:
#         rules.intersect_alldifficult.output,
#         rules.invert_alldifficult.output,
#     localrule: True


rule union_readme:
    input:
        common="workflow/templates/common.j2",
        description="workflow/templates/union_description.j2",
        methods="workflow/templates/union_methods.j2",
        bedtools_env="workflow/envs/bedtools.yml",
        _inputs=lambda w: all_union(w["ref_final_key"], w["build_key"]).all_inputs,
    params:
        paths=lambda w: all_union(w["ref_final_key"], w["build_key"]),
    output:
        uni.readme,
    conda:
        "../envs/templates.yml"
    script:
        "../scripts/python/templates/format_readme/format_union.py"
