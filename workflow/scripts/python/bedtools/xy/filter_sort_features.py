from typing import Any
import common.config as cfg
from common.bed import filter_sort_bed, bed_to_stream, intersectBed
from common.io import bgzip_file, check_processes


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    bed_input = cfg.smk_to_input_name(smk, "bed")
    gapless_input = cfg.smk_to_input_name(smk, "gapless")
    genome_input = cfg.smk_to_input_name(smk, "genome")
    bed_output = cfg.smk_to_output(smk)
    log = cfg.smk_to_log(smk)
    level = cfg.smk_to_param_str(smk, "level")

    i = cfg.ChrIndex.from_name_unsafe(ws["sex_chr"])

    rfk = cfg.wc_to_reffinalkey(ws)
    rk = cfg.strip_full_refkey(rfk)

    pat = sconf.refkey_to_xy_ref_chr_pattern(rfk, i)

    bd = sconf.to_build_data(rk, cfg.wc_to_buildkey(ws))
    bf = bd.refdata.strat_inputs.xy_feature_bed_unsafe(i)
    conv = cfg.HapToHapChrConversion(bf.data.chr_pattern, pat, bd.build_chrs)

    df = bf.read(bed_input)
    df_sorted = filter_sort_bed(conv.init_mapper, conv.final_mapper, df)
    level_mask = df_sorted[bf.level_col].str.contains(level)
    df_filtered = df_sorted[level_mask].drop(columns=[bf.level_col])
    # TODO put this in its own rule to simplify script?
    with bed_to_stream(df_filtered) as s:
        p, o = intersectBed(s, gapless_input, genome_input)
        bgzip_file(o, bed_output)
        check_processes([p], log)


main(snakemake, snakemake.config)  # type: ignore
