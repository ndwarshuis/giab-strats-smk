from typing import Any
from pathlib import Path
from common.io import get_md5
import common.config as cfg


def strat_files(p: Path) -> list[Path]:
    with open(p, "r") as f:
        return [Path(s.strip()) for s in f]


def main(smk: Any) -> None:
    i = cfg.smk_to_input(smk)
    o = cfg.smk_to_output(smk)
    ss = strat_files(i)
    with open(o, "w") as op:
        for s in ss:
            h = get_md5(s)
            p = s.relative_to(s.parent.parent)
            op.write(f"{h}  {p}\n")


main(snakemake)  # type: ignore
