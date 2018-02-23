---
layout: page
title: R/qtl2 functions
description: Annotated/categorized list of functions in R/qtl2
---

### Data import

- `read_cross2` - read data for a cross from a set of files
- `read_csv`
- `read_csv_numer`
- `read_pheno`
- `write_control_file`
- `zip_datafiles`


### Data subsetting

- `subset.calc_genoprob`
- `"[".calc_genoprob` - subset genotype probabilties
- `subset.cross2`
- `"[".cross2` - subset a cross2 object
- `subset.scan1`
- `subset_scan1`
- `subset.sim_geno`
- `"[".sim_geno` - subset genotype imputations
- `subset.viterbi`
- `"[".viterbi` - subset inferred genotypes
- `drop_markers`
- `drop_nullmarkers`
- `pull_markers`


### Combining data

- `cbind.calc_genoprob` - combine genotype probabilities for multiple chromosomes
- `rbind.calc_genoprob`
- `cbind.scan1` - combine genome scan results for multiple phenotypes/analyses
- `rbind.scan1`
- `c.scan1perm` - combine genome scan permutation results for multiple replicates
- `cbind.scan1perm` - combine genome scan permutation results for
   multiple phenotypes/analyses
- `rbind.scan1perm`
- `cbind.sim_geno` - combine genotype imputations for multiple chromosomes
- `rbind.sim_geno`
- `cbind.viterbi` - combine inferred genotypes for multiple chromosomes
- `rbind.viterbi`


### Genotype reconstruction

- `calc_genoprob`
- `clean_genoprob`
- `genoprob_to_alleleprob`
- `genoprob_to_snpprob`
- `interp_genoprob`
- `probs_to_grid`
- `pull_genoprobpos`


### Genotype imputation

- `maxmarg`
- `guess_phase`
- `sim_geno`
- `viterbi`


### Kinship matrix calculations

- `calc_kinship`
- `decomp_kinship`
- `scale_kinship`


### Marker maps

- `est_map`
- `calc_grid`
- `map_to_grid`
- `insert_pseudomarkers`
- `interp_map`
- `pick_marker_subset`
- `reduce_markers`
- `top_snps`


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
