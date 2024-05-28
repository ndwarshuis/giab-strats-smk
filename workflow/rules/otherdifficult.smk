from common.config import CoreLevel, all_otherdifficult_paths

odiff = config.to_bed_dirs(CoreLevel.OTHER_DIFFICULT)


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


def all_otherdifficult(ref_final_key, build_key):
    # TODO actually make other work
    return all_otherdifficult_paths(
        config,
        ref_final_key,
        build_key,
        Path(rules.download_gaps.output[0]),
        Path(rules.download_cds.output[0]),
        Path(rules.download_vdj.output[0]),
        Path(rules.download_kir.output[0]),
        Path(rules.download_mhc.output[0]),
        {},
        Path(rules.get_gaps.output[0]),
        Path(rules.remove_vdj_gaps.output[0]),
        Path(rules.remove_kir_gaps.output[0]),
        Path(rules.remove_mhc_gaps.output[0]),
        {},
    )


rule otherdifficult_readme:
    input:
        common="workflow/templates/common.j2",
        description="workflow/templates/otherdifficult_description.j2",
        methods="workflow/templates/otherdifficult_methods.j2",
        bedtools_env="workflow/envs/bedtools.yml",
        _sources=lambda w: all_otherdifficult(
            w["ref_final_key"], w["build_key"]
        ).all_inputs,
    params:
        paths=lambda w: all_otherdifficult(w["ref_final_key"], w["build_key"]),
    output:
        odiff.readme,
    conda:
        "../envs/templates.yml"
    script:
        "../scripts/python/templates/format_readme/format_otherdifficult.py"
