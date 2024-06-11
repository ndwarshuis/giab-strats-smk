from typing import Any
import common.config as cfg
import subprocess as sp


def main(smk: Any) -> None:
    bed2bb = cfg.smk_to_input_name(smk, "bed2bb")
    genome = cfg.smk_to_input_name(smk, "genome")
    targets = cfg.smk_to_input_name(smk, "targets")
    with open(targets, "r") as i, open(cfg.smk_to_output(smk), "w") as o:
        for t in i:
            outfile = t.replace(".bed", ".bb")
            sp.run([bed2bb, t, genome, outfile], text=True)
            o.write(outfile + "\n")


main(snakemake)  # type: ignore
