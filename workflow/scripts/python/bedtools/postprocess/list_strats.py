from typing import Any
import common.config as cfg


def main(smk: Any) -> None:
    o = cfg.smk_to_output(smk)
    with open(o, "w") as f:
        for s in sorted(smk.input):
            f.write(s + "\n")


main(snakemake)  # type: ignore
