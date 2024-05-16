from common.config import CoreLevel

telo = config.to_bed_dirs(CoreLevel.TELOMERES)


rule find_telomeres:
    input:
        ref=rules.filter_sort_ref.output["fa"],
        genome=rules.filter_sort_ref.output["genome"],
        gapless=rules.get_gapless.output.auto,
    output:
        telo.final("telomeres"),
    conda:
        "../envs/seqtk.yml"
    log:
        telo.inter.postsort.log / "seqtk.log",
    shell:
        """
        seqtk telo {input.ref} 2> {log} | \
        cut -f1-3 | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c \
        > {output}
        """


rule telomere_readme:
    input:
        common="workflow/templates/common.j2",
        description="workflow/templates/telomeres_description.j2",
        methods="workflow/templates/telomeres_methods.j2",
        seqtk_env="workflow/envs/seqtk.yml",
    params:
        path=rules.find_telomeres.output[0],
    output:
        telo.readme,
    conda:
        "../envs/templates.yml"
    script:
        "../scripts/python/templates/format_readme/format_telomeres.py"
