from common.config import (
    CoreLevel,
    parse_full_refkey_class,
    strip_full_refkey,
    flip_full_refkey,
)

dip = config.to_bed_dirs(CoreLevel.DIPLOID)

# TODO don't hardcode minimap2 params (which might be changed if we move to a
# different asm)


def minimap_inputs(wildcards):
    rk = wildcards["ref_final_key"]
    other_rk = flip_full_refkey(rk)
    out = if_dip1_else(
        True, True, "filter_sort_split_ref", "filter_sort_ref", wildcards
    )
    fa = out["fa"]
    idx = out["index"]

    def expand_rk(path, rk):
        return expand(path, allow_missing=True, ref_final_key=rk)

    return {
        "this_hap": expand_rk(fa, rk),
        "_this_idx": expand_rk(idx, rk),
        "other_hap": expand_rk(fa, other_rk),
        "_other_idx": expand_rk(idx, other_rk),
    }


# Dipcall normally outputs a bed file that roughly corresponds to regions with
# no breaks. Since we are interested in large structural variation in
# addition to small variants, obtain this bed file using the same process and
# then flip it.
rule cross_align_breaks:
    input:
        unpack(minimap_inputs),
    output:
        dip.inter.postsort.data / "breaks_cross_align.paf.gz",
    threads: lambda w: config.thread_per_chromosome(w.ref_final_key, w.build_key, 8, True, True)
    log:
        dip.inter.postsort.log / "cross_align_breaks.log",
    benchmark:
        dip.inter.postsort.bench / "cross_align_breaks.txt"
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_final_key, w.build_key, lambda m: m.crossAlignBreaks
        ),
    conda:
        "../envs/quasi-dipcall.yml"
    shell:
        """
        minimap2 -c --paf-no-hit -t{threads} --cs -z200000,10000,200 -xasm5 \
          {input.this_hap} \
          {input.other_hap} \
          2> {log} | \
        bgzip -c > {output}
        """


rule breaks_cross_alignment_to_bed:
    input:
        paf=rules.cross_align_breaks.output,
        paftools_bin=rules.download_paftools.output,
    output:
        dip.inter.postsort.data / "breaks_cross_align.bed.gz",
    log:
        dip.inter.postsort.log / "breaks_cross_alignment_to_bed.log",
    conda:
        "../envs/quasi-dipcall.yml"
    shell:
        """
        gunzip -c {input.paf} | \
        sort -k6,6 -k8,8n | \
        k8 {input.paftools_bin} call - 2> {log} | \
        grep ^R | \
        cut -f2,3,4 | \
        gzip -c > {output}
        """


# These need to be sorted since the paftools call script needs a sorted input
# but for convenience we don't sort numerically (doing so would require
# threading complex python logic through the otherwise elegant pipeline above).
# We could also presort (rather than sort twice) but this would require lots of
# space to store the paf rather than a relatively small bed file (ie no sequence
# strings)
#
# TODO misleading name since this also complements after sorting
rule sort_breaks_bed:
    input:
        bed=rules.cross_align_breaks.output,
        genome=lambda w: if_dip1_else(
            True, True, "filter_sort_split_ref", "filter_sort_ref", w
        )["genome"],
    output:
        dip.inter.postsort.data / "sorted_breaks_cross_align.bed.gz",
    log:
        dip.inter.postsort.log / "sort_breaks_bed.txt",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/diploid/sort_breaks.py"


# Get variants (large and small) between the two haplotypes using the same
# process dipcall uses to produce its *.pair.vcf.gz file.
rule cross_align_variants:
    input:
        unpack(minimap_inputs),
    output:
        dip.inter.postsort.data / "variant_cross_align.sam.gz",
    threads: lambda w: config.thread_per_chromosome(w.ref_final_key, w.build_key, 8, True, True)
    log:
        dip.inter.postsort.log / "cross_align_variants.log",
    # this is what dipcall does to produce the pair file in one step
    conda:
        "../envs/quasi-dipcall.yml"
    benchmark:
        dip.inter.postsort.bench / "cross_align_variants.txt"
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_final_key, w.build_key, lambda m: m.crossAlignVariants
        ),
    shell:
        """
        minimap2 -a -t{threads} --cs -z200000,10000,200 -xasm5 \
          {input.this_hap} \
          {input.other_hap} \
          2> {log} | \
        bgzip -c > {output}
        """


rule filter_sort_variant_cross_alignment:
    input:
        aux_bin=rules.download_dipcall_aux.output,
        sam=rules.cross_align_variants.output,
    output:
        dip.inter.postsort.data / "sorted_variant_cross_alignments.bam",
    benchmark:
        dip.inter.postsort.bench / "filter_sort_variant_cross_alignment.txt"
    log:
        dip.inter.postsort.log / "filter_sort_variant_cross_alignment.txt",
    conda:
        "../envs/quasi-dipcall.yml"
    threads: 4
    resources:
        mem_mb=lambda w, threads: config.buildkey_to_malloc(
            w.ref_final_key,
            w.build_key,
            lambda m: m.filterSortVariantCrossAlignment,
        ),
    params:
        mem_per_thread=lambda w, threads, resources: int(resources.mem_mb / threads),
    shell:
        """
        k8 {input.aux_bin} samflt {input.sam} | \
        samtools sort -m{params.mem_per_thread}M --threads {threads} \
        2> {log} > {output}
        """


# TODO not sure if I also need to filter out alignments where one chr is
# clearly not mapped to the other, something like this: awk '{if ((gensub("_PATERNAL", "", 1, $1) == gensub("_MATERNAL", "", 1, $3)) || ("@" == substr($0,0,1))) { print $0 }}'


rule variant_cross_alignment_to_bed:
    input:
        bam=rules.filter_sort_variant_cross_alignment.output,
        hap=lambda w: if_dip1_else(
            True, True, "filter_sort_split_ref", "filter_sort_ref", w
        )["fa"],
    output:
        dip.inter.postsort.data / "variant_cross_align.bed.gz",
    conda:
        "../envs/quasi-dipcall.yml"
    # Convert vcf tob bed using POS and length of REF for the region bounds, and
    # add a fourth column for the absolute value of the indel length (for
    # filtering later)
    shell:
        """
        htsbox pileup -q5 -evcf {input.hap} {input.bam} | \
        grep -v '^#' | \
        awk 'OFS="\t" {{
          ref = length($4);
          alt = length($5);
          print $1, $2-1, $2+length($4)-1, (ref>alt) ? ref-alt : alt-ref
        }}' | \
        gzip -c > {output}
        """


rule filter_SNVorSV:
    input:
        rules.variant_cross_alignment_to_bed.output,
    output:
        dip.inter.postsort.data / "SNVorSV_variants.bed.gz",
    params:
        ## TODO don't hard code this
        sv_cutoff=50,
    conda:
        "../envs/quasi-dipcall.yml"
    shell:
        """
        gunzip -c {input} | \
        awk '$4>={params.sv_cutoff} || ($4==0 && ($3-$2==1))' | \
        cut -f1,2,3 | \
        gzip -c > {output}
        """


# gzip here so we can concatenate later without crossing the beams
rule merge_all_hets:
    input:
        rules.variant_cross_alignment_to_bed.output,
        rules.sort_breaks_bed.output,
    output:
        dip.inter.postsort.data / "het_regions.bed.gz",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        multiIntersectBed -i {input} | \
        mergeBed -i stdin | \
        gzip -c > {output}
        """


use rule merge_all_hets as merge_SNVorSV_hets with:
    input:
        rules.filter_SNVorSV.output,
        rules.sort_breaks_bed.output,
    output:
        dip.inter.postsort.data / "SNV_or_SV_het_regions.bed.gz",


rule combine_dip1_hets:
    input:
        unpack(lambda w: combine_dip_inputs("merge_all_hets", w)),
    output:
        dip.inter.postsort.data / "combined_het_regions.bed.gz",
    shell:
        """
        cat {input.hap1} {input.hap2} > {output}
        """


use rule combine_dip1_hets as combine_SNVorSV_dip1_hets with:
    input:
        unpack(lambda w: combine_dip_inputs("merge_SNVorSV_hets", w)),
    output:
        dip.inter.postsort.data / "combined_SNV_SV_het_regions.bed.gz",


rule merge_het_regions:
    input:
        lambda w: if_dip1_else(False, True, "combine_dip1_hets", "merge_all_hets", w),
    output:
        dip.final("het_regions_{merge_len}k"),
    conda:
        "../envs/bedtools.yml"
    params:
        gapless=rules.get_gapless.output.auto,
        genome=rules.filter_sort_ref.output["genome"],
    wildcard_constraints:
        merge_len=f"\d+",
    shell:
        """
        mergeBed -i {input} -d $(({wildcards.merge_len}*1000)) | \
        intersectBed -a stdin -b {params.gapless} -sorted -g {params.genome} | \
        bgzip -c > {output}
        """


use rule merge_het_regions as merge_het_SNVorSV_regions with:
    input:
        lambda w: if_dip1_else(
            False, True, "combine_SNVorSV_dip1_hets", "merge_SNVorSV_hets", w
        ),
    output:
        dip.final("het_SNVorSV_regions_{merge_len}k"),
    wildcard_constraints:
        merge_len=f"\d+",


use rule _invert_autosomal_regions as invert_het_regions with:
    input:
        rules.merge_het_regions.output,
    output:
        dip.final("hom_regions_{merge_len}k"),
    wildcard_constraints:
        merge_len=f"\d+",


use rule invert_het_regions as invert_het_SNVorSV_regions with:
    input:
        rules.merge_het_SNVorSV_regions.output,
    output:
        dip.final("hom_SNVorSV_regions_{merge_len}k"),
    wildcard_constraints:
        merge_len=f"\d+",


def het_hom_inputs(ref_final_key, build_key):
    bd = config.to_build_data(strip_full_refkey(ref_final_key), build_key)
    return expand(
        rules.invert_het_regions.output
        + rules.invert_het_SNVorSV_regions.output
        + rules.merge_het_regions.output
        + rules.merge_het_SNVorSV_regions.output,
        merge_len=bd.build.include.hets,
        ref_final_key=ref_final_key,
        build_key=build_key,
    )
