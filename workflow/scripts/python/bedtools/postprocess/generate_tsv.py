from pathlib import Path
from typing import Any
import re
import common.config as cfg


def main(smk: Any) -> None:
    ip = cfg.smk_to_input(smk)
    op = cfg.smk_to_output(smk)
    suffix = smk.params["suffix"]
    with open(ip, "r") as i, open(op, "w") as o:
        for f in i:
            fp = Path(f.strip())
            level = re.sub("^[^_]+_", "", fp.name.replace(".bed.gz", ""))
            o.write(f"{level}\t{Path(fp.parent.name) / fp.name}\n")


main(snakemake)  # type: ignore
