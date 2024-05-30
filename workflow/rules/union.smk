from common.config import CoreLevel

uni = config.to_bed_dirs(CoreLevel.UNION)


# TODO this function is silly and wet
def segdup_and_map_inputs(wildcards):
    rk = wildcards["ref_final_key"]
    bk = wildcards["build_key"]
    bd = config.to_build_data_full(rk, bk)
    s = {"segdups": rules.merge_superdups.output[0]} if bd.have_and_want_segdups else {}
    l = (
        {"lowmap": rules.merge_nonunique.output.all_lowmap}
        if bd.have_and_want_mappability
        else {}
    )
    return {**s, **l}


rule intersect_segdup_and_map:
    input:
        unpack(segdup_and_map_inputs),
        # segdups=rules.merge_superdups.output,
        # lowmap=lambda w: nonunique_inputs(w.ref_final_key, w.build_key)["all_lowmap"],
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
def all_difficult_inputs(wildcards):
    rk = wildcards["ref_final_key"]
    bk = wildcards["build_key"]
    bd = config.to_build_data_full(rk, bk)
    gc = {"gc": gc_inputs(rk, bk)["widest_extreme"]} if bd.want_gc else {}
    return {
        **gc,
        **{
            "_segdups_lowmap": rules.intersect_segdup_and_map.output[0],
            "repeats": rules.merge_HPs_and_TRs.output[0],
            "xy": all_xy_features(rk, bk),
        },
    }


# TODO include other difficult in here as well
use rule intersect_segdup_and_map as intersect_alldifficult with:
    input:
        unpack(all_difficult_inputs),
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


rule all_segdup_and_map:
    input:
        rules.intersect_segdup_and_map.output,
        rules.invert_segdup_and_map.output,
    localrule: True


rule all_alldifficult:
    input:
        rules.intersect_alldifficult.output,
        rules.invert_alldifficult.output,
    localrule: True


rule union_readme:
    input:
        unpack(segdup_and_map_inputs),
        unpack(all_difficult_inputs),
        common="workflow/templates/common.j2",
        description="workflow/templates/union_description.j2",
        methods="workflow/templates/union_methods.j2",
        bedtools_env="workflow/envs/bedtools.yml",
    params:
        segdup_map_path=rules.intersect_segdup_and_map.output[0],
        not_segdup_map_path=rules.invert_segdup_and_map.output[0],
        alldifficult_path=rules.intersect_alldifficult.output[0],
        not_alldifficult_path=rules.invert_alldifficult.output[0],
    output:
        uni.readme,
    conda:
        "../envs/templates.yml"
    script:
        "../scripts/python/templates/format_readme/format_union.py"
