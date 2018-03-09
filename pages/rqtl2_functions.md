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

- `cbind.calc_genoprob` - combine genotype probabilities for multiple
  sets of individuals
- `rbind.calc_genoprob` - combine genotype probabilities for multiple
  chromosomes but on the same set of individuals
- `cbind.scan1` - combine genome scan results for multiple phenotypes/analyses
- `rbind.scan1` - combine genome scan results for different chromosomes
- `c.scan1perm` - combine genome scan permutation results for multiple replicates
- `cbind.scan1perm` - combine genome scan permutation results for
   multiple phenotypes/analyses
- `rbind.scan1perm` - combine genome scan permutation results for
  multiple chromosomes
- `cbind.sim_geno` - combine genotype imputations for multiple chromosomes
- `rbind.sim_geno` - combine genotype imputations for different individuals
- `cbind.viterbi` - combine inferred genotypes for multiple chromosomes
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


### Genotype imputation

- `maxmarg` - for each individual at each position, find genotype with
  maximum marginal probability
- `guess_phase` - turn imputed genotypes into phased genotypes along chromosomes
- `sim_geno` - multiple imputations of underlying genotypes given
  marker data
- `viterbi` - find mostly likely sequence of true genotypes given
  marker data


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

- `est_herit`
- `fit1`
- `scan1`
- `scan1blup`
- `scan1coef`
- `scan1perm`
- `scan1snps`
- `bayes_int`
- `lod_int`


### QTL summaries

- `summary.scan1perm`
- `max.scan1`
- `find_peaks`
- `max_scan1`
- `maxlod`
- `summary_scan1perm`
- `top_snps`


### Data diagnostics

- `check_cross2`
- `calc_entropy`
- `calc_errorlod`
- `calc_geno_freq`
- `chisq_colpairs`
- `convert2cross2`
- `compare_geno`
- `compare_genoprob`
- `count_xo`
- `find_ibd_segments`
- `find_map_gaps`
- `reduce_map_gaps`
- `locate_xo`
- `max_compare_geno`
- `summary_compare_geno`
- `compare_maps`


### Data summaries

- `max.compare_geno`
- `summary.compare_geno`
- `summary.cross2`
- `chr_lengths`
- `chr_names`
- `covar_names`
- `find_marker`
- `find_markerpos`
- `ind_ids`
- `ind_ids_covar`
- `ind_ids_geno`
- `ind_ids_gnp`
- `ind_ids_pheno`
- `marker_names`
- `n_chr`
- `n_covar`
- `n_ind`
- `n_ind_covar`
- `n_ind_geno`
- `n_ind_gnp`
- `n_ind_pheno`
- `n_mar`
- `n_missing`
- `n_pheno`
- `n_phenocovar`
- `n_typed`
- `pheno_names`
- `phenocovar_names`
- `tot_mar`


### QTL plots

- `plot.scan1`
- `plot.scan1coef`
- `plot_coef`
- `plot_coefCC`
- `plot_genes`
- `plot_peaks`
- `plot_pxg`
- `plot_scan1`
- `plot_snpasso`
- `xpos_scan1`


### Diagnostic plots

- `plot.calc_genoprob`
- `plot_genoprob`
- `plot_genoprobcomp`
- `plot_onegeno`


### SNP/gene databases

- `create_gene_query_func`
- `create_variant_query_func`
- `calc_sdp`
- `invert_sdp`
- `index_snps`


### Utility functions

- `batch_cols`
- `get_common_ids`
- `get_x_covar`
- `mat2strata`
- `qtl2version`


### Boring print functions

- `print.cross2`
- `print.summary.compare_geno`
- `print.summary.cross2`
- `print.summary.scan1perm`


### Other functions
