from common.config import CoreLevel

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
