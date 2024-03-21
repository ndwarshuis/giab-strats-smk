# Functions that I use in lots of rules across rule files. Note that while these
# are (mostly) python functions and I would rather put them in a real python
# file so they can be linted, they make more sense here since they require
# some snakemake-specific functions that don't have types anyways, so they
# might as well go here.


def expand_final_to_src(paths, wildcards):
    """Turn a final ref key into a source refkey.

    This is only used in rules that are based directly on the reference, in
    which case this is trivial because the final and source refkeys are always
    the same. This rule is necessary because we make a distinction at the
    wildcard level b/t source and final.

    """
    return expand(paths, allow_missing=True, ref_src_key=wildcards.ref_final_key)


def expand_split_final_to_src(paths, wildcards):
    """Turn a final ref key into a source refkey.

    This is only used in rules that are based directly on the reference, in
    which case this is trivial because the final and source refkeys are always
    the same. This rule is necessary because we make a distinction at the
    wildcard level b/t source and final.

    """
    rk = wildcards.ref_final_key
    rk_ = rk
    return expand(paths, allow_missing=True, ref_src_key=rk)


def to_output_pattern(dirs, name, wildcards):
    """Create and output pattern for a normalization rule."""
    return expand(
        dirs.inter.filtersort.subbed / f"{name}.bed.gz",
        build_key=wildcards.build_key,
    )[0]


def read_checkpoint(rulename, wildcards, other=[]):
    """Read a normalization checkpoint"""
    r = getattr(checkpoints, rulename)
    rk = wildcards.ref_final_key
    c = r.get(
        ref_key=strip_full_refkey(rk),
        build_key=wildcards.build_key,
        **{k: wildcards[k] for k in other},
    )
    with c.output[0].open() as f:
        return config.refkey_to_normalization_path(rk, f)


def read_named_checkpoint(rulename, name, wildcards):
    """Read a normalization checkpoint (with name)"""
    # TODO not DRY
    r = getattr(checkpoints, rulename)
    rfk = wildcards.ref_final_key
    rk = strip_full_refkey(wildcards.ref_final_key)
    c = r.get(ref_key=rk, build_key=wildcards.build_key)
    # TODO isn't there supposed to be a [0] after [name]?
    with c.output[name].open() as f:
        return config.refkey_to_normalization_path(rfk, f)


def bed_src_inputs(pathlist, f, wildcards):
    """Expand pathlist given a refkey and a function to lookup a source file
    for which refsrckeys can be defined."""
    rk = wildcards.ref_key
    rsks = config.refkey_to_bed_refsrckeys(f, rk)
    return expand(pathlist, ref_src_key=rsks)


def bed_src_bd_inputs(pathlist, f, wildcards):
    """Like 'bed_src_inputs' but for build data. Used for benchmark files."""
    rk = wildcards.ref_key
    bk = wildcards.build_key
    rsks = config.buildkey_to_bed_refsrckeys(f, rk, bk)
    return expand(pathlist, allow_missing=True, ref_src_key=rsks, build_key=bk)


def to_bed_src(f, wildcards):
    """Shorthand for looking up bed file source using a function and wildcards."""
    return config.refsrckey_to_bed_src(f, wildcards.ref_src_key)


def combine_dip_inputs(rule, wildcards):
    """Expand a rule path for both haplotypes.

    This assumes that the ref_final_key wildcard does not already have a
    haplotype. This is useful for dip1 cases where the files sometimes need to
    be split into two haplotypes and recombined.
    """
    # NOTE: the wildcard value of ref_final_key will not have a haplotype (as
    # per logic of dip1 references)
    return {
        r.hap.name: expand(
            getattr(rules, rule).output,
            allow_missing=True,
            ref_final_key=r.name,
        )
        for r in config.refkey_append_if_dip1(wildcards.ref_final_key)
    }


def if_dip1_else(split, nohap, dip1, not_dip1, wildcards):
    """Retrieve a rule if the reference is dip1 or something else.

    split: require dip1 refkey to have a haplotype and error otherwise if True
    nohap: throw and error if refkey is haploid if True
    """
    t = config.refkey_is_dip1(wildcards["ref_final_key"], split, nohap)
    return (getattr(rules, dip1) if t else getattr(rules, not_dip1)).output
