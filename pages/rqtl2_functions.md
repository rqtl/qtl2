---
layout: page
title: R/qtl2 functions
description: Annotated/categorized list of functions in R/qtl2
---

### Data import

- [`read_cross2`](https://search.r-project.org/CRAN/refmans/qtl2/html/read_cross2.html) - read data for a cross from a set of files
- [`read_csv`](https://search.r-project.org/CRAN/refmans/qtl2/html/read_csv.html) - read a csv file, using a particular set of options
- [`read_csv_numer`](https://search.r-project.org/CRAN/refmans/qtl2/html/read_csv_numer.html) - like `read_csv` but assuming the contents are
  strictly numeric
- [`read_pheno`](https://search.r-project.org/CRAN/refmans/qtl2/html/read_pheno.html) - read phenotype data from a CSV file, plus
  (optionally) phenotype covariate data from a separate CSV file
- [`write_control_file`](https://search.r-project.org/CRAN/refmans/qtl2/html/write_control_file.html) - write the control file for a set of QTL data
- [`zip_datafiles`](https://search.r-project.org/CRAN/refmans/qtl2/html/zip_datafiles.html) - zip a set of data files (in the format read by `read_cross2`)


### Data subsetting

- [`subset.cross2`](https://search.r-project.org/CRAN/refmans/qtl2/html/subset.cross2.html) - subset a cross2 object by individuals or chromosomes
- [`"[".cross2`](https://search.r-project.org/CRAN/refmans/qtl2/html/subset_cross2.html) - shorthand for `subset.cross2`, with the form
  `mycross[ind, chr]`
- [`subset.calc_genoprob`](https://search.r-project.org/CRAN/refmans/qtl2/html/subset.calc_genoprob.html) subset genotype probabilties by indivduals
  or chromosomes
- [`"[".calc_genoprob`](https://search.r-project.org/CRAN/refmans/qtl2/html/subset.calc_genoprob.html) - shorthand for `subset.calc_genoprob`,
  with the form `probs[ind, chr]`
- [`subset.scan1`](https://search.r-project.org/CRAN/refmans/qtl2/html/subset_scan1.html) - subset genome scan results by chromosome or column
- [`subset_scan1`](https://search.r-project.org/CRAN/refmans/qtl2/html/subset_scan1.html) - the same as `subset.scan1`
- [`subset.sim_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/subset.sim_geno.html) - subset genotype imputations by individuals or chromosomes
- [`"[".sim_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/subset.sim_geno.html) - shorthand for `subset.sim_geno`, with the form
  `simg[ind, chr]`
- [`subset.viterbi`](https://search.r-project.org/CRAN/refmans/qtl2/html/subset.viterbi.html) - subset inferred genotypes
- [`"[".viterbi`](https://search.r-project.org/CRAN/refmans/qtl2/html/subset.viterbi.html) - shorthand for `subset.viterbi`, with the form
  `geno[ind, chr]`
- [`drop_markers`](https://search.r-project.org/CRAN/refmans/qtl2/html/drop_markers.html) - drop a vector of markers from a `"cross2"` object
- [`drop_nullmarkers`](https://search.r-project.org/CRAN/refmans/qtl2/html/drop_nullmarkers.html) - drop markers with no genotypes from a
  `"cross2"` object
- [`pull_markers`](https://search.r-project.org/CRAN/refmans/qtl2/html/pull_markers.html) - drop all except a vector of markers from `"cross2"` object


### Combining data

- [`cbind_expand`](https://search.r-project.org/CRAN/refmans/qtl2/html/cbind_expand.html) - Like `cbind()` but using row names to align the
  rows and expanding with missing values as necessary
- [`cbind.calc_genoprob`](https://search.r-project.org/CRAN/refmans/qtl2/html/cbind.calc_genoprob.html) - combine genotype probabilities for multiple
  chromosomes but on the same set of individuals
- [`rbind.calc_genoprob`](https://search.r-project.org/CRAN/refmans/qtl2/html/rbind.calc_genoprob.html) - combine genotype probabilities for different individuals
- [`cbind.scan1`](https://search.r-project.org/CRAN/refmans/qtl2/html/cbind.scan1.html) - combine genome scan results for multiple phenotypes/analyses
- [`rbind.scan1`](https://search.r-project.org/CRAN/refmans/qtl2/html/rbind.scan1.html) - combine genome scan results for different chromosomes
- [`c.scan1perm`](https://search.r-project.org/CRAN/refmans/qtl2/html/c.scan1perm.html) - combine genome scan permutation results for multiple replicates
- [`cbind.scan1perm`](https://search.r-project.org/CRAN/refmans/qtl2/html/cbind.scan1perm.html) - combine genome scan permutation results for
   multiple phenotypes/analyses
- [`rbind.scan1perm`](https://search.r-project.org/CRAN/refmans/qtl2/html/rbind.scan1perm.html) - combine genome scan permutation results for
  multiple chromosomes
- [`cbind.sim_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/cbind.sim_geno.html) - combine genotype imputations for multiple chromosomes
   but on the same set of individuals
- [`rbind.sim_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/rbind.sim_geno.html) - combine genotype imputations for different individuals
- [`cbind.viterbi`](https://search.r-project.org/CRAN/refmans/qtl2/html/cbind.viterbi.html) - combine inferred genotypes for multiple chromosomes
   but on the same set of individuals
- [`rbind.viterbi`](https://search.r-project.org/CRAN/refmans/qtl2/html/rbind.viterbi.html) - combine inferred genotypes for different individuals


### Genotype reconstruction

- [`calc_genoprob`](https://search.r-project.org/CRAN/refmans/qtl2/html/calc_genoprob.html) - calculate conditional genotype probabilities given
  marker data
- [`clean_genoprob`](https://search.r-project.org/CRAN/refmans/qtl2/html/clean_genoprob.html) - clean up genotype probabilities, setting small
   values to 0
- [`genoprob_to_alleleprob`](https://search.r-project.org/CRAN/refmans/qtl2/html/genoprob_to_alleleprob.html) - convert genotype probabilities to allele dosages
- [`genoprob_to_snpprob`](https://search.r-project.org/CRAN/refmans/qtl2/html/genoprob_to_snpprob.html) - convert genotype probabilities to SNP probabilities
- [`interp_genoprob`](https://search.r-project.org/CRAN/refmans/qtl2/html/interp_genoprob.html) - linear interpolation of genotype probabilities,
  for example to get two sets onto the same map for comparison purposes
- [`probs_to_grid`](https://search.r-project.org/CRAN/refmans/qtl2/html/probs_to_grid.html) - subset genotype probabilities to a grid of
  pseudomarkers
- [`pull_genoprobpos`](https://search.r-project.org/CRAN/refmans/qtl2/html/pull_genoprobpos.html) - pull out the genotype probabilities for a
  particular position
- [`pull_genoprobint`](https://search.r-project.org/CRAN/refmans/qtl2/html/pull_genoprobint.html) - pull out the genotype probabilities for an
  interval


### Genotype imputation

- [`maxmarg`](https://search.r-project.org/CRAN/refmans/qtl2/html/maxmarg.html) - for each individual at each position, find genotype with
  maximum marginal probability
- [`guess_phase`](https://search.r-project.org/CRAN/refmans/qtl2/html/guess_phase.html) - turn imputed genotypes into phased genotypes along chromosomes
- [`sim_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/sim_geno.html) - multiple imputations of underlying genotypes given
  marker data
- [`viterbi`](https://search.r-project.org/CRAN/refmans/qtl2/html/viterbi.html) - find mostly likely sequence of true genotypes given
  marker data
- [`predict_snpgeno`](https://search.r-project.org/CRAN/refmans/qtl2/html/predict_snpgeno.html) - predict SNP genotypes in a multiparent
  population from inferred genotypes plus founder strains' SNP alleles.


### Kinship matrix calculations

- [`calc_kinship`](https://search.r-project.org/CRAN/refmans/qtl2/html/calc_kinship.html) - calculate genetic similarity among individuals
- [`decomp_kinship`](https://search.r-project.org/CRAN/refmans/qtl2/html/decomp_kinship.html) - calculate eigen decomposition of a kinship matrix
- [`scale_kinship`](https://search.r-project.org/CRAN/refmans/qtl2/html/scale_kinship.html) - scale kinship matrix to be like a correlation matrix


### Marker maps

- [`est_map`](https://search.r-project.org/CRAN/refmans/qtl2/html/est_map.html) - re-estimate the inter-marker distances in a genetic map
- [`insert_pseudomarkers`](https://search.r-project.org/CRAN/refmans/qtl2/html/insert_pseudomarkers.html) - add pseudomarkers into a map of genetic markers
- [`calc_grid`](https://search.r-project.org/CRAN/refmans/qtl2/html/calc_grid.html) - Calculate indicators of which pseudomarker positions
  are along a fixed grid
- [`map_to_grid`](https://search.r-project.org/CRAN/refmans/qtl2/html/map_to_grid.html) - subset a map object to the locations on some grid
- [`interp_map`](https://search.r-project.org/CRAN/refmans/qtl2/html/interp_map.html) - Use interpolate to convert from one map to another
- [`reduce_markers`](https://search.r-project.org/CRAN/refmans/qtl2/html/reduce_markers.html) - Reduce marker map to the largest subset that are
  some distance apart


### QTL analysis

- [`est_herit`](https://search.r-project.org/CRAN/refmans/qtl2/html/est_herit.html) - estimate heritability with linear mixed model
- [`fit1`](https://search.r-project.org/CRAN/refmans/qtl2/html/fit1.html) - fit a single-QTL model at a single position
- [`scan1`](https://search.r-project.org/CRAN/refmans/qtl2/html/scan1.html) - genome scan with a single-QTL model
- [`scan1perm`](https://search.r-project.org/CRAN/refmans/qtl2/html/scan1perm.html) - permutation test to establish statistical significance
  in genome scan
- [`scan1coef`](https://search.r-project.org/CRAN/refmans/qtl2/html/scan1coef.html) - calculate QTL effects in scan along one chromosome
- [`scan1blup`](https://search.r-project.org/CRAN/refmans/qtl2/html/scan1blup.html) - like `scan1coef`, but calculating treating QTL
  effects as random and calculating BLUPs
- [`scan1snps`](https://search.r-project.org/CRAN/refmans/qtl2/html/scan1snps.html) - single-QTL scan over SNPs in a multi-parent population
- [`scan1max`](https://search.r-project.org/CRAN/refmans/qtl2/html/scan1max.html) - genome-wide maximum LOD score from genome scan


### QTL summaries

- [`maxlod`](https://search.r-project.org/CRAN/refmans/qtl2/html/maxlod.html) - calculate genome-wide maximum LOD score in genome scan results
- [`max.scan1`](https://search.r-project.org/CRAN/refmans/qtl2/html/max_scan1.html) - calculate maximum LOD score in genome scan and the
  position at which it occurred
- [`max_scan1`](https://search.r-project.org/CRAN/refmans/qtl2/html/max_scan1.html) - the same as `max.scan1`
- [`find_peaks`](https://search.r-project.org/CRAN/refmans/qtl2/html/find_peaks.html) - find QTL peaks in genome scan results
- [`lod_int`](https://search.r-project.org/CRAN/refmans/qtl2/html/lod_int.html) - calculate LOD support intervals from genome scan results
- [`bayes_int`](https://search.r-project.org/CRAN/refmans/qtl2/html/bayes_int.html) - calculate approximate Bayes intervals for QTL position
  from genome scan results
- [`summary.scan1perm`](https://search.r-project.org/CRAN/refmans/qtl2/html/summary_scan1perm.html) - calculate significance thresholds from genome scan
  permutation results
- [`summary_scan1perm`](https://search.r-project.org/CRAN/refmans/qtl2/html/summary_scan1perm.html) - same as `summary.scan1perm`
- [`top_snps`](https://search.r-project.org/CRAN/refmans/qtl2/html/top_snps.html) - find the top SNPs from a SNP association scan


### Data diagnostics

- [`check_cross2`](https://search.r-project.org/CRAN/refmans/qtl2/html/check_cross2.html) - check for inconsistencies or errors in a `"cross2"` object
- [`calc_entropy`](https://search.r-project.org/CRAN/refmans/qtl2/html/calc_entropy.html) - calculate entropy from genotype probabilities, for
  each individual and position
- [`calc_errorlod`](https://search.r-project.org/CRAN/refmans/qtl2/html/calc_errorlod.html) - calculate genotyping error LOD scores to help
  identify potential genotyping errors and problem markers or individuals
- [`calc_geno_freq`](https://search.r-project.org/CRAN/refmans/qtl2/html/calc_geno_freq.html) - calculate genotype frequencies, by individual or
  marker, from genotype probabilities
- [`calc_het`](https://search.r-project.org/CRAN/refmans/qtl2/html/calc_het.html) - Calculate heterozygosities, by individual or marker,
  from genotype probabilities
- [`chisq_colpairs`](https://search.r-project.org/CRAN/refmans/qtl2/html/chisq_colpairs.html) - Perform chi-square test for independence for all
  pairs of columns of a matrix
- [`convert2cross2`](https://search.r-project.org/CRAN/refmans/qtl2/html/convert2cross2.html) - convert an R/qtl1 `"cross"` object to the R/qtl2
  `"cross2"` format
- [`compare_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/compare_geno.html) - compare genotypes for all pairs of individuals, to
  look for possible sample duplicates
- [`compare_genoprob`](https://search.r-project.org/CRAN/refmans/qtl2/html/compare_genoprob.html) - compare two sets of genotype probabilities for
  one individual on a single chromosome
- [`summary.compare_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/summary_compare_geno.html) - summarize the results of `compare_geno`
- [`summary_compare_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/summary_compare_geno.html) - same as `summary.compare_geno`
- [`max.compare_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/max_compare_geno.html) - from the results of `compare_geno`, show the
  pair with most similar genotypes
- [`max_compare_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/max_compare_geno.html) - same as `max.compare_geno`
- [`count_xo`](https://search.r-project.org/CRAN/refmans/qtl2/html/count_xo.html) - count the number of crossovers in each individual on
  each chromosome, from matrices of inferred genotypes
- [`locate_xo`](https://search.r-project.org/CRAN/refmans/qtl2/html/locate_xo.html) - locate the positions of crossovers in each individual
  on each chromosome, from matrices of inferred genotypes.
- [`find_ibd_segments`](https://search.r-project.org/CRAN/refmans/qtl2/html/find_ibd_segments.html) - in genotypes of a set of inbred lines, find
  genomic segments that are identity-by-descent (IBD)
- [`compare_maps`](https://search.r-project.org/CRAN/refmans/qtl2/html/compare_maps.html) - compare two marker maps, to identify markers
  present in one but not the other, or on different chromosomes or in
  different orders between the maps.
- [`find_map_gaps`](https://search.r-project.org/CRAN/refmans/qtl2/html/find_map_gaps.html) - find large gaps between markers in a genetic map
- [`reduce_map_gaps`](https://search.r-project.org/CRAN/refmans/qtl2/html/reduce_map_gaps.html) - reduce the lengths of gaps in a genetic map
- [`calc_raw_het`](https://search.r-project.org/CRAN/refmans/qtl2/html/calc_raw_het.html) - Calculate heterozygosity in the raw SNP genotypes
- [`calc_raw_maf`](https://search.r-project.org/CRAN/refmans/qtl2/html/calc_raw_maf.html) - Calculate the minor allele frequency in the raw SNP
  genotypes
- [`calc_raw_geno_freq`](https://search.r-project.org/CRAN/refmans/qtl2/html/calc_raw_geno_freq.html) - Calculate the genotype frequencies in the raw
  SNP data
- [`calc_raw_founder_maf`](https://search.r-project.org/CRAN/refmans/qtl2/html/calc_raw_founder_maf.html) - Calculate the minor allele frequency in the
  founder strains' SNP genotypes


### Data summaries

- [`summary.cross2`](https://search.r-project.org/CRAN/refmans/qtl2/html/summary.cross2.html) - summarize a `"cross2"` object
- [`chr_names`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - names of chromosomes in a `"cross2"` object
- [`marker_names`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - names of markers in a `"cross2"` object
- [`pheno_names`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - names of phenotypes in a `"cross2"` object
- [`phenocovar_names`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - names of "phenotype covariates" (metadata about
  phenotypes) in a `"cross2"` object
- [`covar_names`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - names of covariates in a `"cross2"` object
- [`ind_ids`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - return IDs for all individuals in a `"cross2"` object
- [`ind_ids_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - return IDs for all individuals in a `"cross2"`
  object that have genotype data
- [`ind_ids_pheno`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - return IDs for all individuals in a `"cross2"`
  object that have phenotype data
- [`ind_ids_gnp`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - return IDs for all individuals in a `"cross2"`
  object that have both genotype and phenotype data
- [`ind_ids_covar`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - return IDs for all individuals in a `"cross2"`
  object that have covariate data
- [`n_chr`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - number of chromosomes in a `"cross2"` object
- [`n_ind`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - number of individuals in a `"cross2"` object
- [`n_ind_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - number of individuals in a `"cross2"` object that
  have genotype data
- [`n_ind_pheno`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - number of individuals in a `"cross2"` object that
  have phenotype data
- [`n_ind_gnp`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - number of individuals in a `"cross2"` object that
  have both genotype and phenotype data
- [`n_ind_covar`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - number of individuals in a `"cross2"` object that
  have covariate data
- [`n_mar`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - number of markers on each chromosome in a `"cross2"`
  object
- [`tot_mar`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - total number of markers in a `"cross2"` object
- [`n_pheno`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - number of phenotypes in a `"cross2"` object
- [`n_covar`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html)  - number of covariates in a `"cross2"` object
- [`n_phenocovar`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - number of "phenotype covariates" (metadata on
  phenotypes) in a `"cross2"` object
- [`chr_lengths`](https://search.r-project.org/CRAN/refmans/qtl2/html/chr_lengths.html) - calculate chromosome lengths for a map object
- [`find_marker`](https://search.r-project.org/CRAN/refmans/qtl2/html/find_marker.html) - find marker closest to a particular genomic position
- [`find_markerpos`](https://search.r-project.org/CRAN/refmans/qtl2/html/find_markerpos.html) - find the position of a marker
- [`n_missing`](https://search.r-project.org/CRAN/refmans/qtl2/html/n_missing.html) - number of missing genotypes, by individual or marker
- [`n_typed`](https://search.r-project.org/CRAN/refmans/qtl2/html/n_typed.html) - number of genotypes, by individual or marker
- [`founders`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - names of the founder strains
- [`n_founders`](https://search.r-project.org/CRAN/refmans/qtl2/html/basic_summaries.html) - number of founder strains


### QTL plots

- [`plot.scan1`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_scan1.html) - plot genome scan results
- [`plot_scan1`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_scan1.html) - same as `plot.scan1`
- [`xpos_scan1`](https://search.r-project.org/CRAN/refmans/qtl2/html/xpos_scan1.html) - determine the x-axis location of a particular genomic
  position in a genome scan plot (for adding annotations)
- [`add_threshold`](https://search.r-project.org/CRAN/refmans/qtl2/html/add_threshold.html) - Add horizontal line at a significance threshold to
  a genome scan plot.
- [`plot.scan1coef`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_coef.html) - plot QTL effects along a chromosome
- [`plot_coef`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_coef.html) - same as `plot.scan1coef`
- [`plot_coefCC`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_coef.html) - like `plot_coef` but assuming there are 8 effects
  and using the standard colors for the Collaborative Cross (`CCcolors`)
- [`plot_snpasso`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_snpasso.html) - plot SNP association results
- [`plot_genes`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_genes.html) - plot locations of a set of genes
- [`plot_sdp`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_sdp.html) - plot strain distribution patterns of SNPs in a region
- [`plot_peaks`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_peaks.html) - plot a summary of QTL positions for multiple
  phenotypes, using the results of `find_peaks`
- [`plot_lodpeaks`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_lodpeaks.html) - scatterplot of LOD scores vs QTL peak locatiosn
  (possibly with intervals) for multiple traits
- [`plot_pxg`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_pxg.html) - plot phenotype versus QTL genotypes


### Diagnostic plots

- [`plot.calc_genoprob`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_genoprob.html) - plot the genotype probabilities for one
  individual on one chromosome, as a heat map
- [`plot_genoprob`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_genoprob.html) - the same as `plot.calc_genoprob`
- [`plot_onegeno`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_onegeno.html) - plot one individual's genome-wide genotypes
- [`plot_genoprobcomp`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_genoprobcomp.html) - plot a comparison of two sets of genotype
  probabilities for one individual on one chromosome, as a bivariate
  heat map
- [`plot_compare_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_compare_geno.html) - Plot histogram of the results of
  `compare_geno()`
- [`plot.compare_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/plot_compare_geno.html) - Same as `plot_compare_geno()`


### SNP/gene databases

- [`create_variant_query_func`](https://search.r-project.org/CRAN/refmans/qtl2/html/create_variant_query_func.html) - create a function to connect to a SQLite
  database of founder variant information and return a data frame with
  variants for a selected region
- [`create_gene_query_func`](https://search.r-project.org/CRAN/refmans/qtl2/html/create_gene_query_func.html) - create a function to connect to a SQLite
  database of gene annotations and return a data frame with genes in a
  selected region
- [`calc_sdp`](https://search.r-project.org/CRAN/refmans/qtl2/html/calc_sdp.html) - convert founder SNP genotypes to a numeric code for the
  strain distribution pattern
- [`invert_sdp`](https://search.r-project.org/CRAN/refmans/qtl2/html/invert_sdp.html) - the inverse of `calc_sdp`
- [`index_snps`](https://search.r-project.org/CRAN/refmans/qtl2/html/index_snps.html) - partition SNPs into groups that are contained within
  common marker intervals and have the same strain distribution
  pattern, and create an index to a set of distinct SNPs, one per
  partition
- [`find_index_snp`](https://search.r-project.org/CRAN/refmans/qtl2/html/find_index_snp.html) - For a particular SNP, find the corresponding
  indexed SNP.
- [`create_snpinfo`](https://search.r-project.org/CRAN/refmans/qtl2/html/create_snpinfo.html) - Create a table of SNP information from a cross2 object.
- [`sdp2char`](https://search.r-project.org/CRAN/refmans/qtl2/html/sdp2char.html) - convert strain distribution pattern numeric codes to
  more meaningful character strings



### Utility functions

- [`batch_cols`](https://search.r-project.org/CRAN/refmans/qtl2/html/batch_cols.html) - identify batches of columns of a matrix that have the
  same pattern of missing values
- [`batch_vec`](https://search.r-project.org/CRAN/refmans/qtl2/html/batch_vec.html) - split a vector into batches, for help in balancing
  parallel code
- [`get_common_ids`](https://search.r-project.org/CRAN/refmans/qtl2/html/get_common_ids.html) - find IDs that are present in all of the input objects
- [`get_x_covar`](https://search.r-project.org/CRAN/refmans/qtl2/html/get_x_covar.html) - from a `"cross2"` object, get the matrix of
  covariates to be used for the null hypothesis when performing QTL
  analysis on the X chromosome
- [`mat2strata`](https://search.r-project.org/CRAN/refmans/qtl2/html/mat2strata.html) - use the rows of a matrix to define a set of strata
  for a stratified permutation test
- [`replace_ids`](https://search.r-project.org/CRAN/refmans/qtl2/html/replace_ids.html) - Replace the individual IDs in an object
- [`replace_ids.calc_genoprob`](https://search.r-project.org/CRAN/refmans/qtl2/html/replace_ids.html)  - Replace the individual IDs in a `"calc_genoprob"` object
- [`replace_ids.cross2`](https://search.r-project.org/CRAN/refmans/qtl2/html/replace_ids.html) - Replace the individual IDs in a `"cross2"` object
- [`replace_ids.sim_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/replace_ids.html) - Replace the individual IDs in a `"sim_geno"` object
- [`replace_ids.viterbi`](https://search.r-project.org/CRAN/refmans/qtl2/html/replace_ids.html) - Replace the individual IDs in a `"viterbi"` object
- [`replace_ids.data.frame`](https://search.r-project.org/CRAN/refmans/qtl2/html/replace_ids.html) - Replace the individual IDs (in row names) in a data frame
- [`replace_ids.matrix`](https://search.r-project.org/CRAN/refmans/qtl2/html/replace_ids.html) - Replace the individual IDs (in row names) in a matrix
- [`align_scan1_map`](https://search.r-project.org/CRAN/refmans/qtl2/html/align_scan1_map.html) - aligns the markers/pseudomarkers in a `"scan1"`
  object (output by `scan1()`) and a marker map.
- [`clean`](https://search.r-project.org/CRAN/refmans/qtl2/html/clean.html) - clean an object
- [`clean.scan1`](https://search.r-project.org/CRAN/refmans/qtl2/html/clean_scan1.html) - clean a `"scan"` object (replacing negative values
  with `NA` and removing rows were all values are `NA`.
- [`clean_scan1`](https://search.r-project.org/CRAN/refmans/qtl2/html/clean_scan1.html) - the same as `clean.scan1`.
- [`clean.calc_genoprob`](https://search.r-project.org/CRAN/refmans/qtl2/html/clean_genoprob.html) - clean a `"calc_genoprob"` object (setting
  small values to 0)
- [`clean_genoprob`](https://search.r-project.org/CRAN/refmans/qtl2/html/clean_genoprob.html) - same as `clean.calc_genoprob`
- [`qtl2version`](https://search.r-project.org/CRAN/refmans/qtl2/html/qtl2version.html) - print the installed version of [R/qtl2](https://kbroman.org/qtl2)
- [`recode_snps`](https://search.r-project.org/CRAN/refmans/qtl2/html/recode_snps.html) - Recode the SNP genotypes so that `1` is for the
  major allele in the founders



### Boring print functions

- [`print.cross2`](https://search.r-project.org/CRAN/refmans/qtl2/html/print.cross2.html) - print method for a `"cross2"` object
- `print.summary.cross2` - print method for the output of `summary.cross2`
- [`print.summary.compare_geno`](https://search.r-project.org/CRAN/refmans/qtl2/html/summary_compare_geno.html) - print method for the output of `summary.compare_geno`
- [`print.summary.scan1perm`](https://search.r-project.org/CRAN/refmans/qtl2/html/print.summary.scan1perm.html) - print method for the output of `summary.scan1perm`



### Newly added functions (in development version)
