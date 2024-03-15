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
        config.log_root_dir / "tools" / "repseq_build.log",
    shell:
        "make -C {input} 2>&1 > {log} && mv {input}/repseq {output}"


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


# NOTE this entire thing is simply to get bigBedToBed, which I can't seem to
# install using conda without it complaining that the constraints for my other
# packages can't be satisfied (something about openssl and python needing
# different version). This binary is also available on the UCSC download portal
# but it doesn't seem to be version-tagged. "easy" solution: build from source.
# Start by building just the kent libraries and then just build the bigBedToBed
# binary.
rule build_kent:
    input:
        rules.unpack_kent.output[0],
    output:
        config.tools_bin_dir / "bigBedToBed",
    threads: 8
    log:
        config.log_root_dir / "tools" / "kent_build.log",
    conda:
        "../envs/kent.yml"
    params:
        destdir=lambda _, output: str(Path(output[0]).parent),
    shell:
        """
        here=$(pwd)
        libpath="$CONDA_PREFIX/include"

        cd $here/{input}/src
        make clean -j{threads} 2>&1 > /dev/null
        CFLAGS="-I$libpath" make libs -j{threads} 2>&1 > $here/{log}

        cd $here/{input}/src/utils/bigBedToBed
        make DESTDIR=$here/{params.destdir} BINDIR=/ 2>&1 > $here/{log}
        """


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
