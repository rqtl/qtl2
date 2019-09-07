---
layout: page
title: R/qtl2 functions
description: Annotated/categorized list of functions in R/qtl2
---

### Data import

- `read_cross2` - read data for a cross from a set of files
- `read_csv` - read a csv file, using a particular set of options
- `read_csv_numer` - like `read_csv` but assuming the contents are
  strictly numeric
- `read_pheno` - read phenotype data from a CSV file, plus
  (optionally) phenotype covariate data from a separate CSV file
- `write_control_file` - write the control file for a set of QTL data
- `zip_datafiles` - zip a set of data files (in the format read by `read_cross2`)


### Data subsetting

- `subset.cross2` - subset a cross2 object by individuals or chromosomes
- `"[".cross2` - shorthand for `subset.cross2`, with the form
  `mycross[ind, chr]`
- `subset.calc_genoprob` subset genotype probabilties by indivduals
  or chromosomes
- `"[".calc_genoprob` - shorthand for `subset.calc_genoprob`,
  with the form `probs[ind, chr]`
- `subset.scan1` - subset genome scan results by chromosome or column
- `subset_scan1` - the same as `subset.scan1`
- `subset.sim_geno` - subset genotype imputations by individuals or chromosomes
- `"[".sim_geno` - shorthand for `subset.sim_geno`, with the form
  `simg[ind, chr]`
- `subset.viterbi` - subset inferred genotypes
- `"[".viterbi` - shorthand for `subset.viterbi`, with the form
  `geno[ind, chr]`
- `drop_markers` - drop a vector of markers from a `"cross2"` object
- `drop_nullmarkers` - drop markers with no genotypes from a
  `"cross2"` object
- `pull_markers` - drop all except a vector of markers from `"cross2"` object


### Combining data

- `cbind_expand` - Like `cbind()` but using row names to align the
  rows and expanding with missing values as necessary
- `cbind.calc_genoprob` - combine genotype probabilities for multiple
  chromosomes but on the same set of individuals
- `rbind.calc_genoprob` - combine genotype probabilities for different individuals
- `cbind.scan1` - combine genome scan results for multiple phenotypes/analyses
- `rbind.scan1` - combine genome scan results for different chromosomes
- `c.scan1perm` - combine genome scan permutation results for multiple replicates
- `cbind.scan1perm` - combine genome scan permutation results for
   multiple phenotypes/analyses
- `rbind.scan1perm` - combine genome scan permutation results for
  multiple chromosomes
- `cbind.sim_geno` - combine genotype imputations for multiple chromosomes
   but on the same set of individuals
- `rbind.sim_geno` - combine genotype imputations for different individuals
- `cbind.viterbi` - combine inferred genotypes for multiple chromosomes
   but on the same set of individuals
- `rbind.viterbi` - combine inferred genotypes for different individuals


### Genotype reconstruction

- `calc_genoprob` - calculate conditional genotype probabilities given
  marker data
- `clean_genoprob` - clean up genotype probabilities, setting small
   values to 0
- `genoprob_to_alleleprob` - convert genotype probabilities to allele dosages
- `genoprob_to_snpprob` - convert genotype probabilities to SNP probabilities
- `interp_genoprob` - linear interpolation of genotype probabilities,
  for example to get two sets onto the same map for comparison purposes
- `probs_to_grid` - subset genotype probabilities to a grid of
  pseudomarkers
- `pull_genoprobpos` - pull out the genotype probabilities for a
  particular position
- `pull_genoprobint` - pull out the genotype probabilities for an
  interval


### Genotype imputation

- `maxmarg` - for each individual at each position, find genotype with
  maximum marginal probability
- `guess_phase` - turn imputed genotypes into phased genotypes along chromosomes
- `sim_geno` - multiple imputations of underlying genotypes given
  marker data
- `viterbi` - find mostly likely sequence of true genotypes given
  marker data
- `predict_snpgeno` - predict SNP genotypes in a multiparent
  population from inferred genotypes plus founder strains' SNP alleles.


### Kinship matrix calculations

- `calc_kinship` - calculate genetic similarity among individuals
- `decomp_kinship` - calculate eigen decomposition of a kinship matrix
- `scale_kinship` - scale kinship matrix to be like a correlation matrix


### Marker maps

- `est_map` - re-estimate the inter-marker distances in a genetic map
- `insert_pseudomarkers` - add pseudomarkers into a map of genetic markers
- `calc_grid` - Calculate indicators of which pseudomarker positions
  are along a fixed grid
- `map_to_grid` - subset a map object to the locations on some grid
- `interp_map` - Use interpolate to convert from one map to another
- `reduce_markers` - Reduce marker map to the largest subset that are
  some distance apart


### QTL analysis

- `est_herit` - estimate heritability with linear mixed model
- `fit1` - fit a single-QTL model at a single position
- `scan1` - genome scan with a single-QTL model
- `scan1perm` - permutation test to establish statistical significance
  in genome scan
- `scan1coef` - calculate QTL effects in scan along one chromosome
- `scan1blup` - like `scan1coef`, but calculating treating QTL
  effects as random and calculating BLUPs
- `scan1snps` - single-QTL scan over SNPs in a multi-parent population


### QTL summaries

- `maxlod` - calculate genome-wide maximum LOD score in genome scan results
- `max.scan1` - calculate maximum LOD score in genome scan and the
  position at which it occurred
- `max_scan1` - the same as `max.scan1`
- `find_peaks` - find QTL peaks in genome scan results
- `lod_int` - calculate LOD support intervals from genome scan results
- `bayes_int` - calculate approximate Bayes intervals for QTL position
  from genome scan results
- `summary.scan1perm` - calculate significance thresholds from genome scan
  permutation results
- `summary_scan1perm` - same as `summary.scan1perm`
- `top_snps` - find the top SNPs from a SNP association scan


### Data diagnostics

- `check_cross2` - check for inconsistencies or errors in a `"cross2"` object
- `calc_entropy` - calculate entropy from genotype probabilities, for
  each individual and position
- `calc_errorlod` - calculate genotyping error LOD scores to help
  identify potential genotyping errors and problem markers or individuals
- `calc_geno_freq` - calculate genotype frequencies, by individual or
  marker, from genotype probabilities
- `calc_het` - Calculate heterozygosities, by individual or marker,
  from genotype probabilities
- `chisq_colpairs` - Perform chi-square test for independence for all
  pairs of columns of a matrix
- `convert2cross2` - convert an R/qtl1 `"cross"` object to the R/qtl2
  `"cross2"` format
- `compare_geno` - compare genotypes for all pairs of individuals, to
  look for possible sample duplicates
- `compare_genoprob` - compare two sets of genotype probabilities for
  one individual on a single chromosome
- `summary.compare_geno` - summarize the results of `compare_geno`
- `summary_compare_geno` - same as `summary.compare_geno`
- `max.compare_geno` - from the results of `compare_geno`, show the
  pair with most similar genotypes
- `max_compare_geno` - same as `max.compare_geno`
- `count_xo` - count the number of crossovers in each individual on
  each chromosome, from matrices of inferred genotypes
- `locate_xo` - locate the positions of crossovers in each individual
  on each chromosome, from matrices of inferred genotypes.
- `find_ibd_segments` - in genotypes of a set of inbred lines, find
  genomic segments that are identity-by-descent (IBD)
- `compare_maps` - compare two marker maps, to identify markers
  present in one but not the other, or on different chromosomes or in
  different orders between the maps.
- `find_map_gaps` - find large gaps between markers in a genetic map
- `reduce_map_gaps` - reduce the lengths of gaps in a genetic map


### Data summaries

- `summary.cross2` - summarize a `"cross2"` object
- `chr_names` - names of chromosomes in a `"cross2"` object
- `marker_names` - names of markers in a `"cross2"` object
- `pheno_names` - names of phenotypes in a `"cross2"` object
- `phenocovar_names` - names of "phenotype covariates" (metadata about
  phenotypes) in a `"cross2"` object
- `covar_names` - names of covariates in a `"cross2"` object
- `ind_ids` - return IDs for all individuals in a `"cross2"` object
- `ind_ids_geno` - return IDs for all individuals in a `"cross2"`
  object that have genotype data
- `ind_ids_pheno` - return IDs for all individuals in a `"cross2"`
  object that have phenotype data
- `ind_ids_gnp` - return IDs for all individuals in a `"cross2"`
  object that have both genotype and phenotype data
- `ind_ids_covar` - return IDs for all individuals in a `"cross2"`
  object that have covariate data
- `n_chr` - number of chromosomes in a `"cross2"` object
- `n_ind` - number of individuals in a `"cross2"` object
- `n_ind_geno` - number of individuals in a `"cross2"` object that
  have genotype data
- `n_ind_pheno` - number of individuals in a `"cross2"` object that
  have phenotype data
- `n_ind_gnp` - number of individuals in a `"cross2"` object that
  have both genotype and phenotype data
- `n_ind_covar` - number of individuals in a `"cross2"` object that
  have covariate data
- `n_mar` - number of markers on each chromosome in a `"cross2"`
  object
- `tot_mar` - total number of markers in a `"cross2"` object
- `n_pheno` - number of phenotypes in a `"cross2"` object
- `n_covar`  - number of covariates in a `"cross2"` object
- `n_phenocovar` - number of "phenotype covariates" (metadata on
  phenotypes) in a `"cross2"` object
- `chr_lengths` - calculate chromosome lengths for a map object
- `find_marker` - find marker closest to a particular genomic position
- `find_markerpos` - find the position of a marker
- `n_missing` - number of missing genotypes, by individual or marker
- `n_typed` - number of genotypes, by individual or marker


### QTL plots

- `plot.scan1` - plot genome scan results
- `plot_scan1` - same as `plot.scan1`
- `xpos_scan1` - determine the x-axis location of a particular genomic
  position in a genome scan plot (for adding annotations)
- `add_threshold` - Add horizontal line at a significance threshold to
  a genome scan plot.
- `plot.scan1coef` - plot QTL effects along a chromosome
- `plot_coef` - same as `plot.scan1coef`
- `plot_coefCC` - like `plot_coef` but assuming there are 8 effects
  and using the standard colors for the Collaborative Cross (`CCcolors`)
- `plot_snpasso` - plot SNP association results
- `plot_genes` - plot locations of a set of genes
- `plot_peaks` - plot a summary of QTL positions for multiple
  phenotypes, using the results of `find_peaks`
- `plot_lodpeaks` - scatterplot of LOD scores vs QTL peak locatiosn
  (possibly with intervals) for multiple traits
- `plot_pxg` - plot phenotype versus QTL genotypes


### Diagnostic plots

- `plot.calc_genoprob` - plot the genotype probabilities for one
  individual on one chromosome, as a heat map
- `plot_genoprob` - the same as `plot.calc_genoprob`
- `plot_onegeno` - plot one individual's genome-wide genotypes
- `plot_genoprobcomp` - plot a comparison of two sets of genotype
  probabilities for one individual on one chromosome, as a bivariate
  heat map


### SNP/gene databases

- `create_variant_query_func` - create a function to connect to a SQLite
  database of founder variant information and return a data frame with
  variants for a selected region
- `create_gene_query_func` - create a function to connect to a SQLite
  database of gene annotations and return a data frame with genes in a
  selected region
- `calc_sdp` - convert founder SNP genotypes to a numeric code for the
  strain distribution pattern
- `invert_sdp` - the inverse of `calc_sdp`
- `index_snps` - partition SNPs into groups that are contained within
  common marker intervals and have the same strain distribution
  pattern, and create an index to a set of distinct SNPs, one per
  partition
- `find_index_snp` - For a particular SNP, find the corresponding
  indexed SNP.
- `create_snpinfo` - Create a table of SNP information from a cross2 object.



### Utility functions

- `batch_cols` - identify batches of columns of a matrix that have the
  same pattern of missing values
- `batch_vec` - split a vector into batches, for help in balancing
  parallel code
- `get_common_ids` - find IDs that are present in all of the input objects
- `get_x_covar` - from a `"cross2"` object, get the matrix of
  covariates to be used for the null hypothesis when performing QTL
  analysis on the X chromosome
- `mat2strata` - use the rows of a matrix to define a set of strata
  for a stratified permutation test
- `replace_ids` - Replace the individual IDs in an object
- `replace_ids.calc_genoprob`  - Replace the individual IDs in a `"calc_genoprob"` object
- `replace_ids.cross2` - Replace the individual IDs in a `"cross2"` object
- `replace_ids.sim_geno` - Replace the individual IDs in a `"sim_geno"` object
- `replace_ids.viterbi` - Replace the individual IDs in a `"viterbi"` object
- `align_scan1_map` - aligns the markers/pseudomarkers in a `"scan1"`
  object (output by `scan1()`) and a marker map.
- `clean` - clean an object
- `clean.scan1` - clean a `"scan"` object (replacing negative values
  with `NA` and removing rows were all values are `NA`.
- `clean_scan1` - the same as `clean.scan1`.
- `clean.calc_genoprob` - clean a `"calc_genoprob"` object (setting
  small values to 0)
- `clean_genoprob` - same as `clean.calc_genoprob`
- `qtl2version` - print the installed version of [R/qtl2](https://kbroman.org/qtl2)



### Boring print functions

- `print.cross2` - print method for a `"cross2"` object
- `print.summary.cross2` - print method for the output of `summary.cross2`
- `print.summary.compare_geno` - print method for the output of `summary.compare_geno`
- `print.summary.scan1perm` - print method for the output of `summary.scan1perm`

### Newly added functions (in development version)

- `plot_compare_geno` - Plot histogram of the results of
  `compare_geno()` (_diagnostic plots_)
- `plot.compare_geno` - Same as `plot_compare_geno()` (_diagnostic plots_)
- `recode_snps` - Recode the SNP genotypes so that `1` is for the
  major allele in the founders (_utilities_)
- `calc_raw_het` - Calculate heterozygosity in the raw SNP genotypes (_diagnostics_)
- `calc_raw_maf` - Calculate the minor allele frequency in the raw SNP
  genotypes (_diagnostics_)
- `calc_raw_geno_freq` - Calculate the genotype frequencies in the raw
  SNP data (_diagnostics_)
- `calc_raw_founder_maf` - Calculate the minor allele frequency in the
  founder strains' SNP genotypes (_diagnostics_)
- `founders` - names of the founder strains (_data summaries_)
- `n_founders` - number of founder strains (_data summaries_)
- `scan1max` - genome-wide maximum LOD score from genome scan (_QTL analysis_)
