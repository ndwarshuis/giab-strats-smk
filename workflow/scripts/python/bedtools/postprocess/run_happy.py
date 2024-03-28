import gzip
import os
import subprocess as sp
from typing import Any
from pathlib import Path
import common.config as cfg


def gzip_empty(path: Path) -> bool:
    with gzip.open(path, "rb") as f:
        return len(f.read(1)) == 0


def main(smk: Any) -> None:
    # ASSUME if filtered benchmark bed is empty then we have nothing to
    # benchmark because the desired chromosomes in this build are not included
    # in the benchmark
    ref_path = cfg.smk_to_input_name(smk, "ref")
    bench_vcf_path = cfg.smk_to_input_name(smk, "bench_vcf")
    bench_bed_path = cfg.smk_to_input_name(smk, "bench_bed")
    query_vcf_path = cfg.smk_to_input_name(smk, "query_vcf")
    strats_path = cfg.smk_to_input_name(smk, "strats")
    log_path = cfg.smk_to_log(smk)
    out_path = cfg.smk_to_output(smk)
    prefix = cfg.smk_to_param_str(smk, "prefix")

    with open(log_path, "w") as f:
        if gzip_empty(bench_bed_path):
            out_path.touch()
            f.write("no overlapping chromosomes, writing empty dummy file\n")
        else:
            cmd: list[str] = [
                "hap.py",
                *["--engine", "vcfeval"],
                "--verbose",
                *["--threads", str(smk.threads)],
                *["--stratification", str(strats_path)],
                *["-f", str(bench_bed_path)],
                *["-o", prefix],
                str(bench_vcf_path),
                str(query_vcf_path),
            ]
            res = sp.run(
                cmd,
                stderr=sp.STDOUT,
                stdout=f,
                env={"HGREF": str(ref_path), **os.environ},
            )
            if res.returncode != 0:
                exit(1)


main(snakemake)  # type: ignore
