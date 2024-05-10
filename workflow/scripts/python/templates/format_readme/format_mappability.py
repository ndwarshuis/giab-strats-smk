import jinja2 as j2
from os.path import basename
from typing import Any
import common.config as cfg
from urllib.parse import unquote
import json


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    rk = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)
    bd = sconf.to_build_data(cfg.strip_full_refkey(rk), bk)
    map_params = [*bd.build.include.mappability]

    lowmap_path = cfg.smk_to_input_name(smk, "lowmap")
    notlowmap_path = cfg.smk_to_input_name(smk, "notinlowmap")
    out = cfg.smk_to_output(smk)

    with open(lowmap_path, "r") as f:
        u = json.load(f)
        all_lowmap: str = basename(u["all_lowmap"])
        single_lowmap: list[str] = [*map(basename, u["single_lowmap"])]

    env = j2.Environment(
        loader=j2.FileSystemLoader("workflow/templates"),
        undefined=j2.StrictUndefined,
    )
    common = env.get_template("common.j2")
    overview = env.get_template("mappability_overview.j2")
    methods = env.get_template("mappability_methods.j2")

    otxt = overview.render(
        single_lowmap_files=single_lowmap,
        all_lowmap_file=all_lowmap,
        not_all_lowmap_file=basename(notlowmap_path),
        params=map_params,
    )

    mtxt = methods.render(gemurl=unquote(sconf.tools.gemlib))

    txt = common.render(
        overview_text=otxt,
        methods_text=mtxt,
        stratification_group=cfg.CoreLevel.MAPPABILITY.value,
        pi_name=sconf.docs.pi_name,
        pi_email=sconf.docs.pi_email,
        contact_name=sconf.docs.contact_name,
        contact_email=sconf.docs.contact_email,
    )

    with open(out, "w") as f:
        f.write(txt)


main(snakemake, snakemake.config)  # type: ignore
