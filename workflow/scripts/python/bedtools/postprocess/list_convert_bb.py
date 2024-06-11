from typing import Any
import common.config as cfg
import common.io as io
import subprocess as sp


def main(smk: Any) -> None:
    bed2bb = cfg.smk_to_input_name(smk, "bed2bb")
    genome = cfg.smk_to_input_name(smk, "genome")
    targets = cfg.smk_to_inputs_name(smk, "targets")
    out = cfg.smk_to_output(smk)
    log = cfg.smk_to_log(smk)
    with open(out, "w") as o, open(log, "wa") as g:
        for t in targets:
            if not io.gzip_is_empty(t):
                outfile = str(t).replace(".bed.gz", ".bb")
                sp.run([bed2bb, t, genome, outfile], text=True, stderr=g)
                o.write(outfile + "\n")


main(snakemake)  # type: ignore
