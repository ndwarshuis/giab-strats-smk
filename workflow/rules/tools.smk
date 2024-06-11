from pathlib import Path

################################################################################
# download a bunch of stuff to run GEM
#
# NOTE: this is in bioconda, but the bioconda version does not have gem-2-wig
# for some reason


gemlib_bin = Path("GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin")


rule download_gem:
    output:
        config.tools_src_dir / "gemlib.tbz2",
    params:
        url=config.tools.gemlib,
    conda:
        "../envs/utils.yml"
    localrule: True
    shell:
        "curl -sS -L -o {output} {params.url}"


rule unpack_gem:
    input:
        rules.download_gem.output,
    output:
        # called by other binaries
        config.tools_bin_dir / "gem-indexer_fasta2meta+cont",
        config.tools_bin_dir / "gem-indexer_bwt-dna",
        config.tools_bin_dir / "gem-indexer_generate",
        # the things I actually need
        indexer=config.tools_bin_dir / "gem-indexer",
        mappability=config.tools_bin_dir / "gem-mappability",
        gem2wig=config.tools_bin_dir / "gem-2-wig",
    params:
        bins=lambda wildcards, output: " ".join(
            str(gemlib_bin / basename(o)) for o in output
        ),
    shell:
        """
        mkdir -p {config.tools_bin_dir} && \
        tar xjf {input} \
        --directory {config.tools_bin_dir} \
        --strip-components=2 \
        {params.bins}
        """


################################################################################
# repseq


rule download_repseq:
    output:
        config.tools_src_dir / "repseq.tar.gz",
    params:
        url=config.tools.repseq,
    conda:
        "../envs/utils.yml"
    localrule: True
    shell:
        "curl -f -sS -L -o {output} {params.url}"


rule unpack_repseq:
    input:
        rules.download_repseq.output,
    output:
        directory(config.tools_make_dir / "repseq"),
    localrule: True
    shell:
        """
        mkdir {output} && \
        tar xzf {input} --directory {output} --strip-components=1
        """


rule build_repseq:
    input:
        rules.unpack_repseq.output,
    output:
        config.tools_bin_dir / "repseq",
    conda:
        "../envs/build.yml"
    log:
        config.log_tools_dir / "repseq_build.log",
    shell:
        "make -C {input} 2>&1 > {log} && mv {input}/repseq {output}"


################################################################################
# bigBedToBed


use rule download_repseq as download_kent with:
    output:
        config.tools_src_dir / "kent.tar.gz",
    params:
        url=config.tools.kent,


use rule unpack_repseq as unpack_kent with:
    input:
        rules.download_kent.output,
    output:
        directory(config.tools_make_dir / "kent"),


# NOTE this entire thing is simply to get bigBedToBed and bedToBigBed, which I
# can't seem to install using conda without it complaining that the constraints
# for my other packages can't be satisfied (something about openssl and python
# needing different version). This binary is also available on the UCSC download
# portal but it doesn't seem to be version-tagged. "easy" solution: build from
# source. Start by building just the kent libraries and then just build the
# binaries.
rule build_kent:
    input:
        rules.unpack_kent.output[0],
    output:
        bb2bed=config.tools_bin_dir / "bigBedToBed",
        bed2bb=config.tools_bin_dir / "bedToBigBed",
    threads: 8
    log:
        config.log_tools_dir / "tools" / "kent_build.log",
    conda:
        "../envs/kent.yml"
    params:
        bb2bed_dest=lambda _, output: str(Path(output["bb2bed"]).parent),
        bed2bb_dest=lambda _, output: str(Path(output["bed2bb"]).parent),
    shell:
        """
        here=$(pwd)
        libpath="$CONDA_PREFIX/include"

        cd $here/{input}/src
        make clean -j{threads} 2>&1 > /dev/null
        CFLAGS="-I$libpath" CPLUS_INCLUDE_PATH="$libpath" \
          make libs -j{threads} > $here/{log} 2>&1 

        cd $here/{input}/src/utils/bigBedToBed
        make DESTDIR=$here/{params.bb2bed_dest} BINDIR=/ 2>&1 > $here/{log}

        cd $here/{input}/src/utils/bedToBigBed
        make DESTDIR=$here/{params.bed2bb_dest} BINDIR=/ 2>&1 > $here/{log}
        """


################################################################################
# other stuff


use rule download_repseq as download_paftools with:
    output:
        config.tools_src_dir / "paftools.js",
    params:
        url=config.tools.paftools,
    localrule: True


use rule download_repseq as download_dipcall_aux with:
    output:
        config.tools_src_dir / "dipcall_aux.js",
    params:
        url=config.tools.dipcall_aux,
    localrule: True


rule all_tools:
    input:
        rules.unpack_gem.output,
        rules.unpack_gem.output,
        rules.build_kent.output,
        rules.download_paftools.output,
        rules.download_dipcall_aux.output,
