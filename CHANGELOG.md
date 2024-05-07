# 5.0.0

Pipeline now works on diploid genomes. This is a major release and many previous
configuration options will totally break due to the way the yaml file needed to
be represented for the diploid case. However, the output for the haploid case
should not change.

Also added a new stratification group "Diploid" which contains heterozygous and
homozygous regions between two diploid haplotypes.

Other new features:

- validation now includes intra-chromosomal coverage plots showing the coverage
  of each bed file within 1Mbp windows
- bed files can now be specified directly in the yaml config
- bigbeds can now be imported directly
- numerous speed and memory improvements; most rules will now run in constant
  memory, and those that don't can be allocated more memory on a per-build basis
  using the "malloc" directive in yaml

Other breaking changes:

- In the yaml config, there are now three toplevel stratification categories
  corresponding to "haploid" (self explanatory), "diploid1" (two haplotypes in
  one diploid bed file/fasta), or "diploid2" (one haplotype per bed/fasta
  diploid pair). In the latter two cases, the 1 or 2 designates how the final
  stratification beds will be split (ie into one file or two files).
  Furthermore, the two diploid configurations can take either diploid1 or
  diploid2 beds as input (they will be split/combined as needed), which makes
  the configuration very flexible but also more complex. See
  `config/testing.yaml` for examples.
- CDS regions no longer require the "feature table" file; instead use the
  chromosomal pattern directly to map chromosome accession numbers in the GFF to
  chromosomal indices

# 4.1.1

- remove extra tarball parent directories
- fix comment line skipping when reading bed files
- don't make validation directory hidden
- remove extra static config

# 4.1.0

- automatically derive vdj regions from refseq

# 4.0.0

- allow custom chromosome mappings (to deal with the HG2 paternal asms having
'chrX_MATERNAL' and vice versa)

# 3.0.0

- generalize chr prefix into a pattern (to allow recognizing chromosome names
  like "chr1_PATERNAL")
- relax constraints in input files; if not provided, output will not be
  generated; this is useful for cases where the input files do not exist.
  - affected strats: low complexity, xy, functional

# 2.10.0

- add AT/GC to low complexity just for homopolymers
- add AT/GC low complexity to benchmark output

# 2.9.0

- remove AT/GC from low complexity (for now)

# 2.8.1

- fix some random typos and bugs

# 2.8.0

- lower max low complexity length to 150bp

# 2.7.1

- fix typos

# 2.7.0

- make small rules not run on slurm
- fix lots of errors involving the final list of strats (actually involved using
  checkpoints for rules with complex output)
- fix formatting errors in checksum file (which didn't actually allow md5sun to
  run previously)

# 2.6.1

- fix chrom order in coverage plots

# 2.6.0

- add option to remove gaps from imported strat beds

# 2.5.1

- use plaintext and not gzipped file when performing md5 checks for resources

# 2.5.0

- make comparison way faster and more intuitive
- fix missing middle GC content file
- automatically make flanking gaps stratification

# 2.4.1

- fix overlapping level bug

# 2.4.0

- make benchmark subsets configurable

# 2.3.0

- automatically make gaps stratification

# 2.2.0

- add comparison functions to config/pipeline to test how much generated
  strats have changed relative to previous versions
- pipeline now fails on http 404 (or other bad request)

# 2.1.1

- fix typo

# 2.1.0

- make other strat name constraint more permissive

# 2.0.0

- add telomere stratification
- allow external beds to be used as stratifications
- add gc/at homopolymer bed files
- add benchmarking to validation postprocessing

# 1.0.0

- in the beginning there was darkness
