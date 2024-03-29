# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

arrange_genes <- function(start, end) {
    .Call(`_qtl2_arrange_genes`, start, end)
}

.bayes_int_plain <- function(lod, pos, prob) {
    .Call(`_qtl2_R_bayes_int_plain`, lod, pos, prob)
}

calc_ll_binreg <- function(X, y, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_calc_ll_binreg`, X, y, maxit, tol, qr_tol, eta_max)
}

calc_coef_binreg <- function(X, y, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_calc_coef_binreg`, X, y, maxit, tol, qr_tol, eta_max)
}

calc_coefSE_binreg <- function(X, y, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_calc_coefSE_binreg`, X, y, maxit, tol, qr_tol, eta_max)
}

fit_binreg <- function(X, y, se = TRUE, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_fit_binreg`, X, y, se, maxit, tol, qr_tol, eta_max)
}

calc_ll_binreg_eigenchol <- function(X, y, maxit = 100L, tol = 1e-6, eta_max = 30.0) {
    .Call(`_qtl2_calc_ll_binreg_eigenchol`, X, y, maxit, tol, eta_max)
}

calc_ll_binreg_eigenqr <- function(X, y, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_calc_ll_binreg_eigenqr`, X, y, maxit, tol, qr_tol, eta_max)
}

calc_coef_binreg_eigenqr <- function(X, y, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_calc_coef_binreg_eigenqr`, X, y, maxit, tol, qr_tol, eta_max)
}

calc_coefSE_binreg_eigenqr <- function(X, y, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_calc_coefSE_binreg_eigenqr`, X, y, maxit, tol, qr_tol, eta_max)
}

fit_binreg_eigenqr <- function(X, y, se = TRUE, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_fit_binreg_eigenqr`, X, y, se, maxit, tol, qr_tol, eta_max)
}

calc_ll_binreg_weighted <- function(X, y, weights, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_calc_ll_binreg_weighted`, X, y, weights, maxit, tol, qr_tol, eta_max)
}

calc_coef_binreg_weighted <- function(X, y, weights, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_calc_coef_binreg_weighted`, X, y, weights, maxit, tol, qr_tol, eta_max)
}

calc_coefSE_binreg_weighted <- function(X, y, weights, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_calc_coefSE_binreg_weighted`, X, y, weights, maxit, tol, qr_tol, eta_max)
}

fit_binreg_weighted <- function(X, y, weights, se = TRUE, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_fit_binreg_weighted`, X, y, weights, se, maxit, tol, qr_tol, eta_max)
}

calc_ll_binreg_weighted_eigenchol <- function(X, y, weights, maxit = 100L, tol = 1e-6, eta_max = 30.0) {
    .Call(`_qtl2_calc_ll_binreg_weighted_eigenchol`, X, y, weights, maxit, tol, eta_max)
}

calc_ll_binreg_weighted_eigenqr <- function(X, y, weights, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_calc_ll_binreg_weighted_eigenqr`, X, y, weights, maxit, tol, qr_tol, eta_max)
}

calc_coef_binreg_weighted_eigenqr <- function(X, y, weights, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_calc_coef_binreg_weighted_eigenqr`, X, y, weights, maxit, tol, qr_tol, eta_max)
}

calc_coefSE_binreg_weighted_eigenqr <- function(X, y, weights, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_calc_coefSE_binreg_weighted_eigenqr`, X, y, weights, maxit, tol, qr_tol, eta_max)
}

fit_binreg_weighted_eigenqr <- function(X, y, weights, se = TRUE, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_fit_binreg_weighted_eigenqr`, X, y, weights, se, maxit, tol, qr_tol, eta_max)
}

.calc_kinship <- function(prob_array) {
    .Call(`_qtl2_calc_kinship`, prob_array)
}

.crosstype_supported <- function(crosstype) {
    .Call(`_qtl2_crosstype_supported`, crosstype)
}

.count_invalid_genotypes <- function(crosstype, genotypes, is_X_chr, is_female, cross_info) {
    .Call(`_qtl2_count_invalid_genotypes`, crosstype, genotypes, is_X_chr, is_female, cross_info)
}

.check_crossinfo <- function(crosstype, cross_info, any_x_chr) {
    .Call(`_qtl2_check_crossinfo`, crosstype, cross_info, any_x_chr)
}

.check_is_female_vector <- function(crosstype, is_female, any_x_chr) {
    .Call(`_qtl2_check_is_female_vector`, crosstype, is_female, any_x_chr)
}

check_handle_x_chr <- function(crosstype, any_x_chr) {
    .Call(`_qtl2_check_handle_x_chr`, crosstype, any_x_chr)
}

.chisq_colpairs <- function(input) {
    .Call(`_qtl2_chisq_colpairs`, input)
}

.clean_genoprob <- function(prob_array, value_threshold = 1e-6, column_threshold = 0.01) {
    .Call(`_qtl2_clean_genoprob`, prob_array, value_threshold, column_threshold)
}

.compare_geno <- function(geno) {
    .Call(`_qtl2_compare_geno`, geno)
}

.count_xo <- function(geno, crosstype, is_X_chr) {
    .Call(`_qtl2_count_xo`, geno, crosstype, is_X_chr)
}

.count_xo_3d <- function(geno_array, crosstype, is_X_chr) {
    .Call(`_qtl2_count_xo_3d`, geno_array, crosstype, is_X_chr)
}

mpp_encode_alleles <- function(allele1, allele2, n_alleles, phase_known) {
    .Call(`_qtl2_mpp_encode_alleles`, allele1, allele2, n_alleles, phase_known)
}

mpp_decode_geno <- function(true_gen, n_alleles, phase_known) {
    .Call(`_qtl2_mpp_decode_geno`, true_gen, n_alleles, phase_known)
}

mpp_is_het <- function(true_gen, n_alleles, phase_known) {
    .Call(`_qtl2_mpp_is_het`, true_gen, n_alleles, phase_known)
}

mpp_geno_names <- function(alleles, is_x_chr) {
    .Call(`_qtl2_mpp_geno_names`, alleles, is_x_chr)
}

invert_founder_index <- function(cross_info) {
    .Call(`_qtl2_invert_founder_index`, cross_info)
}

.is_phase_known <- function(crosstype) {
    .Call(`_qtl2_is_phase_known`, crosstype)
}

.find_dup_markers_notexact <- function(Geno, order, markerloc, adjacent_only) {
    .Call(`_qtl2_find_dup_markers_notexact`, Geno, order, markerloc, adjacent_only)
}

.find_ibd_segments <- function(g1, g2, p, error_prob) {
    .Call(`_qtl2_find_ibd_segments`, g1, g2, p, error_prob)
}

.find_peaks <- function(lod, threshold, peakdrop) {
    .Call(`_qtl2_R_find_peaks`, lod, threshold, peakdrop)
}

.find_peaks_and_lodint <- function(lod, threshold, peakdrop, drop) {
    .Call(`_qtl2_R_find_peaks_and_lodint`, lod, threshold, peakdrop, drop)
}

.find_peaks_and_bayesint <- function(lod, pos, threshold, peakdrop, prob) {
    .Call(`_qtl2_R_find_peaks_and_bayesint`, lod, pos, threshold, peakdrop, prob)
}

fit1_binary_addcovar <- function(genoprobs, pheno, addcovar, weights, se = FALSE, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_fit1_binary_addcovar`, genoprobs, pheno, addcovar, weights, se, maxit, tol, qr_tol, eta_max)
}

fit1_binary_intcovar <- function(genoprobs, pheno, addcovar, intcovar, weights, se = TRUE, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_fit1_binary_intcovar`, genoprobs, pheno, addcovar, intcovar, weights, se, maxit, tol, qr_tol, eta_max)
}

fit1_hk_addcovar <- function(genoprobs, pheno, addcovar, weights, se = FALSE, tol = 1e-12) {
    .Call(`_qtl2_fit1_hk_addcovar`, genoprobs, pheno, addcovar, weights, se, tol)
}

fit1_hk_intcovar <- function(genoprobs, pheno, addcovar, intcovar, weights, se = TRUE, tol = 1e-12) {
    .Call(`_qtl2_fit1_hk_intcovar`, genoprobs, pheno, addcovar, intcovar, weights, se, tol)
}

fit1_pg_addcovar <- function(genoprobs, pheno, addcovar, eigenvec, weights, se = FALSE, tol = 1e-12) {
    .Call(`_qtl2_fit1_pg_addcovar`, genoprobs, pheno, addcovar, eigenvec, weights, se, tol)
}

fit1_pg_intcovar <- function(genoprobs, pheno, addcovar, intcovar, eigenvec, weights, se = TRUE, tol = 1e-12) {
    .Call(`_qtl2_fit1_pg_intcovar`, genoprobs, pheno, addcovar, intcovar, eigenvec, weights, se, tol)
}

geno_names <- function(crosstype, alleles, is_x_chr) {
    .Call(`_qtl2_geno_names`, crosstype, alleles, is_x_chr)
}

nalleles <- function(crosstype) {
    .Call(`_qtl2_nalleles`, crosstype)
}

.genoprob_to_alleleprob <- function(crosstype, prob_array, is_x_chr) {
    .Call(`_qtl2_genoprob_to_alleleprob`, crosstype, prob_array, is_x_chr)
}

.get_x_covar <- function(crosstype, is_female, cross_info) {
    .Call(`_qtl2_get_x_covar`, crosstype, is_female, cross_info)
}

.guess_phase_f2A <- function(geno, deterministic) {
    .Call(`_qtl2_guess_phase_f2A`, geno, deterministic)
}

.guess_phase_f2X <- function(geno, deterministic) {
    .Call(`_qtl2_guess_phase_f2X`, geno, deterministic)
}

.guess_phase_A <- function(geno, crosstype, deterministic) {
    .Call(`_qtl2_guess_phase_A`, geno, crosstype, deterministic)
}

.guess_phase_X <- function(geno, crosstype, is_female, deterministic) {
    .Call(`_qtl2_guess_phase_X`, geno, crosstype, is_female, deterministic)
}

.calc_errorlod <- function(crosstype, probs, genotypes, founder_geno, is_X_chr, is_female, cross_info) {
    .Call(`_qtl2_calc_errorlod`, crosstype, probs, genotypes, founder_geno, is_X_chr, is_female, cross_info)
}

.calc_genoprob <- function(crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, rec_frac, marker_index, error_prob) {
    .Call(`_qtl2_calc_genoprob`, crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, rec_frac, marker_index, error_prob)
}

.calc_genoprob2 <- function(crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, rec_frac, marker_index, error_prob) {
    .Call(`_qtl2_calc_genoprob2`, crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, rec_frac, marker_index, error_prob)
}

.est_map <- function(crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, rec_frac, error_prob, max_iterations, tol, verbose) {
    .Call(`_qtl2_est_map`, crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, rec_frac, error_prob, max_iterations, tol, verbose)
}

.est_map2 <- function(crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, cross_group, unique_cross_group, rec_frac, error_prob, max_iterations, tol, verbose) {
    .Call(`_qtl2_est_map2`, crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, cross_group, unique_cross_group, rec_frac, error_prob, max_iterations, tol, verbose)
}

.sim_geno <- function(crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, rec_frac, marker_index, error_prob, n_draws) {
    .Call(`_qtl2_sim_geno`, crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, rec_frac, marker_index, error_prob, n_draws)
}

.sim_geno2 <- function(crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, rec_frac, marker_index, error_prob, n_draws) {
    .Call(`_qtl2_sim_geno2`, crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, rec_frac, marker_index, error_prob, n_draws)
}

addlog <- function(a, b) {
    .Call(`_qtl2_addlog`, a, b)
}

subtractlog <- function(a, b) {
    .Call(`_qtl2_subtractlog`, a, b)
}

.viterbi <- function(crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, rec_frac, marker_index, error_prob) {
    .Call(`_qtl2_viterbi`, crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, rec_frac, marker_index, error_prob)
}

.viterbi2 <- function(crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, rec_frac, marker_index, error_prob) {
    .Call(`_qtl2_viterbi2`, crosstype, genotypes, founder_geno, is_X_chr, is_female, cross_info, rec_frac, marker_index, error_prob)
}

.interp_genoprob_onechr <- function(genoprob, map, pos_index) {
    .Call(`_qtl2_interp_genoprob_onechr`, genoprob, map, pos_index)
}

interpolate_map <- function(oldpos, oldmap, newmap) {
    .Call(`_qtl2_interpolate_map`, oldpos, oldmap, newmap)
}

find_intervals <- function(pos, map, tol = 1e-8) {
    .Call(`_qtl2_find_intervals`, pos, map, tol)
}

calc_rss_linreg <- function(X, Y, tol = 1e-12) {
    .Call(`_qtl2_calc_rss_linreg`, X, Y, tol)
}

calc_coef_linreg <- function(X, y, tol = 1e-12) {
    .Call(`_qtl2_calc_coef_linreg`, X, y, tol)
}

calc_coefSE_linreg <- function(X, y, tol = 1e-12) {
    .Call(`_qtl2_calc_coefSE_linreg`, X, y, tol)
}

calc_resid_linreg <- function(X, Y, tol = 1e-12) {
    .Call(`_qtl2_calc_resid_linreg`, X, Y, tol)
}

calc_resid_linreg_3d <- function(X, P, tol = 1e-12) {
    .Call(`_qtl2_calc_resid_linreg_3d`, X, P, tol)
}

fit_linreg <- function(X, y, se = TRUE, tol = 1e-12) {
    .Call(`_qtl2_fit_linreg`, X, y, se, tol)
}

fit_linreg_eigenchol <- function(X, y, se) {
    .Call(`_qtl2_fit_linreg_eigenchol`, X, y, se)
}

calc_coef_linreg_eigenchol <- function(X, y) {
    .Call(`_qtl2_calc_coef_linreg_eigenchol`, X, y)
}

calc_coefSE_linreg_eigenchol <- function(X, y) {
    .Call(`_qtl2_calc_coefSE_linreg_eigenchol`, X, y)
}

calc_rss_eigenchol <- function(X, y) {
    .Call(`_qtl2_calc_rss_eigenchol`, X, y)
}

calc_fitted_linreg_eigenchol <- function(X, y) {
    .Call(`_qtl2_calc_fitted_linreg_eigenchol`, X, y)
}

fit_linreg_eigenqr <- function(X, y, se, tol = 1e-12) {
    .Call(`_qtl2_fit_linreg_eigenqr`, X, y, se, tol)
}

calc_coef_linreg_eigenqr <- function(X, y, tol = 1e-12) {
    .Call(`_qtl2_calc_coef_linreg_eigenqr`, X, y, tol)
}

calc_coefSE_linreg_eigenqr <- function(X, y, tol = 1e-12) {
    .Call(`_qtl2_calc_coefSE_linreg_eigenqr`, X, y, tol)
}

calc_rss_eigenqr <- function(X, y, tol = 1e-12) {
    .Call(`_qtl2_calc_rss_eigenqr`, X, y, tol)
}

calc_fitted_linreg_eigenqr <- function(X, y, tol = 1e-12) {
    .Call(`_qtl2_calc_fitted_linreg_eigenqr`, X, y, tol)
}

calc_mvrss_eigenchol <- function(X, Y) {
    .Call(`_qtl2_calc_mvrss_eigenchol`, X, Y)
}

calc_mvrss_eigenqr <- function(X, Y, tol = 1e-12) {
    .Call(`_qtl2_calc_mvrss_eigenqr`, X, Y, tol)
}

calc_resid_eigenchol <- function(X, Y) {
    .Call(`_qtl2_calc_resid_eigenchol`, X, Y)
}

calc_resid_eigenqr <- function(X, Y, tol = 1e-12) {
    .Call(`_qtl2_calc_resid_eigenqr`, X, Y, tol)
}

Rcpp_eigen_decomp <- function(A) {
    .Call(`_qtl2_Rcpp_eigen_decomp`, A)
}

Rcpp_eigen_rotation <- function(K, y, X) {
    .Call(`_qtl2_Rcpp_eigen_rotation`, K, y, X)
}

Rcpp_calc_logdetXpX <- function(X) {
    .Call(`_qtl2_Rcpp_calc_logdetXpX`, X)
}

Rcpp_calcLL <- function(hsq, Kva, y, X, reml = TRUE, logdetXpX = NA_real_) {
    .Call(`_qtl2_Rcpp_calcLL`, hsq, Kva, y, X, reml, logdetXpX)
}

Rcpp_calcLL_mat <- function(hsq, Kva, Y, X, reml = TRUE, logdetXpX = NA_real_) {
    .Call(`_qtl2_Rcpp_calcLL_mat`, hsq, Kva, Y, X, reml, logdetXpX)
}

Rcpp_fitLMM <- function(Kva, y, X, reml = TRUE, check_boundary = TRUE, logdetXpX = NA_real_, tol = 1e-4) {
    .Call(`_qtl2_Rcpp_fitLMM`, Kva, y, X, reml, check_boundary, logdetXpX, tol)
}

Rcpp_fitLMM_mat <- function(Kva, Y, X, reml = TRUE, check_boundary = TRUE, logdetXpX = NA_real_, tol = 1e-4) {
    .Call(`_qtl2_Rcpp_fitLMM_mat`, Kva, Y, X, reml, check_boundary, logdetXpX, tol)
}

.locate_xo <- function(geno, map, crosstype, is_X_chr) {
    .Call(`_qtl2_locate_xo`, geno, map, crosstype, is_X_chr)
}

.lod_int_plain <- function(lod, drop) {
    .Call(`_qtl2_R_lod_int_plain`, lod, drop)
}

find_matching_cols <- function(mat, tol = 1e-12) {
    .Call(`_qtl2_find_matching_cols`, mat, tol)
}

find_lin_indep_cols <- function(mat, tol = 1e-12) {
    .Call(`_qtl2_find_lin_indep_cols`, mat, tol)
}

formX_intcovar <- function(probs, addcovar, intcovar, position, has_intercept = TRUE) {
    .Call(`_qtl2_formX_intcovar`, probs, addcovar, intcovar, position, has_intercept)
}

expand_genoprobs_intcovar <- function(probs, intcovar) {
    .Call(`_qtl2_expand_genoprobs_intcovar`, probs, intcovar)
}

weighted_matrix <- function(mat, weights) {
    .Call(`_qtl2_weighted_matrix`, mat, weights)
}

weighted_3darray <- function(array, weights) {
    .Call(`_qtl2_weighted_3darray`, array, weights)
}

matrix_x_matrix <- function(X, Y) {
    .Call(`_qtl2_matrix_x_matrix`, X, Y)
}

matrix_x_vector <- function(X, y) {
    .Call(`_qtl2_matrix_x_vector`, X, y)
}

matrix_x_3darray <- function(X, A) {
    .Call(`_qtl2_matrix_x_3darray`, X, A)
}

.maxmarg <- function(prob_array, minprob, tol) {
    .Call(`_qtl2_maxmarg`, prob_array, minprob, tol)
}

.predict_snpgeno <- function(allele1, allele2, founder_geno) {
    .Call(`_qtl2_predict_snpgeno`, allele1, allele2, founder_geno)
}

random_int <- function(n, low, high) {
    .Call(`_qtl2_random_int`, n, low, high)
}

get_permutation <- function(n) {
    .Call(`_qtl2_get_permutation`, n)
}

permute_nvector <- function(n_perm, x) {
    .Call(`_qtl2_permute_nvector`, n_perm, x)
}

permute_ivector <- function(n_perm, x) {
    .Call(`_qtl2_permute_ivector`, n_perm, x)
}

permute_nvector_stratified <- function(n_perm, x, strata, n_strata = -1L) {
    .Call(`_qtl2_permute_nvector_stratified`, n_perm, x, strata, n_strata)
}

permute_ivector_stratified <- function(n_perm, x, strata, n_strata = -1L) {
    .Call(`_qtl2_permute_ivector_stratified`, n_perm, x, strata, n_strata)
}

.reduce_markers <- function(pos, min_dist, weights) {
    .Call(`_qtl2_reduce_markers`, pos, min_dist, weights)
}

scan_binary_onechr <- function(genoprobs, pheno, addcovar, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_scan_binary_onechr`, genoprobs, pheno, addcovar, maxit, tol, qr_tol, eta_max)
}

scan_binary_onechr_weighted <- function(genoprobs, pheno, addcovar, weights, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_scan_binary_onechr_weighted`, genoprobs, pheno, addcovar, weights, maxit, tol, qr_tol, eta_max)
}

scan_binary_onechr_intcovar_highmem <- function(genoprobs, pheno, addcovar, intcovar, maxit = 100L, tol = 1e-6, qr_tol = 1e-12) {
    .Call(`_qtl2_scan_binary_onechr_intcovar_highmem`, genoprobs, pheno, addcovar, intcovar, maxit, tol, qr_tol)
}

scan_binary_onechr_intcovar_weighted_highmem <- function(genoprobs, pheno, addcovar, intcovar, weights, maxit = 100L, tol = 1e-6, qr_tol = 1e-12) {
    .Call(`_qtl2_scan_binary_onechr_intcovar_weighted_highmem`, genoprobs, pheno, addcovar, intcovar, weights, maxit, tol, qr_tol)
}

scan_binary_onechr_intcovar_lowmem <- function(genoprobs, pheno, addcovar, intcovar, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_scan_binary_onechr_intcovar_lowmem`, genoprobs, pheno, addcovar, intcovar, maxit, tol, qr_tol, eta_max)
}

scan_binary_onechr_intcovar_weighted_lowmem <- function(genoprobs, pheno, addcovar, intcovar, weights, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_scan_binary_onechr_intcovar_weighted_lowmem`, genoprobs, pheno, addcovar, intcovar, weights, maxit, tol, qr_tol, eta_max)
}

scan_hk_onechr_nocovar <- function(genoprobs, pheno, tol = 1e-12) {
    .Call(`_qtl2_scan_hk_onechr_nocovar`, genoprobs, pheno, tol)
}

scan_hk_onechr <- function(genoprobs, pheno, addcovar, tol = 1e-12) {
    .Call(`_qtl2_scan_hk_onechr`, genoprobs, pheno, addcovar, tol)
}

scan_hk_onechr_weighted <- function(genoprobs, pheno, addcovar, weights, tol = 1e-12) {
    .Call(`_qtl2_scan_hk_onechr_weighted`, genoprobs, pheno, addcovar, weights, tol)
}

scan_hk_onechr_intcovar_highmem <- function(genoprobs, pheno, addcovar, intcovar, tol = 1e-12) {
    .Call(`_qtl2_scan_hk_onechr_intcovar_highmem`, genoprobs, pheno, addcovar, intcovar, tol)
}

scan_hk_onechr_intcovar_weighted_highmem <- function(genoprobs, pheno, addcovar, intcovar, weights, tol = 1e-12) {
    .Call(`_qtl2_scan_hk_onechr_intcovar_weighted_highmem`, genoprobs, pheno, addcovar, intcovar, weights, tol)
}

scan_hk_onechr_intcovar_lowmem <- function(genoprobs, pheno, addcovar, intcovar, tol = 1e-12) {
    .Call(`_qtl2_scan_hk_onechr_intcovar_lowmem`, genoprobs, pheno, addcovar, intcovar, tol)
}

scan_hk_onechr_intcovar_weighted_lowmem <- function(genoprobs, pheno, addcovar, intcovar, weights, tol = 1e-12) {
    .Call(`_qtl2_scan_hk_onechr_intcovar_weighted_lowmem`, genoprobs, pheno, addcovar, intcovar, weights, tol)
}

scan_pg_onechr <- function(genoprobs, pheno, addcovar, eigenvec, weights, tol = 1e-12) {
    .Call(`_qtl2_scan_pg_onechr`, genoprobs, pheno, addcovar, eigenvec, weights, tol)
}

scan_pg_onechr_intcovar_highmem <- function(genoprobs, pheno, addcovar, intcovar, eigenvec, weights, tol = 1e-12) {
    .Call(`_qtl2_scan_pg_onechr_intcovar_highmem`, genoprobs, pheno, addcovar, intcovar, eigenvec, weights, tol)
}

scan_pg_onechr_intcovar_lowmem <- function(genoprobs, pheno, addcovar, intcovar, eigenvec, weights, tol = 1e-12) {
    .Call(`_qtl2_scan_pg_onechr_intcovar_lowmem`, genoprobs, pheno, addcovar, intcovar, eigenvec, weights, tol)
}

scanblup <- function(genoprobs, pheno, addcovar, se, reml, tol = 1e-12) {
    .Call(`_qtl2_scanblup`, genoprobs, pheno, addcovar, se, reml, tol)
}

scancoef_binary_addcovar <- function(genoprobs, pheno, addcovar, weights, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_scancoef_binary_addcovar`, genoprobs, pheno, addcovar, weights, maxit, tol, qr_tol, eta_max)
}

scancoef_binary_intcovar <- function(genoprobs, pheno, addcovar, intcovar, weights, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_scancoef_binary_intcovar`, genoprobs, pheno, addcovar, intcovar, weights, maxit, tol, qr_tol, eta_max)
}

scancoefSE_binary_addcovar <- function(genoprobs, pheno, addcovar, weights, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_scancoefSE_binary_addcovar`, genoprobs, pheno, addcovar, weights, maxit, tol, qr_tol, eta_max)
}

scancoefSE_binary_intcovar <- function(genoprobs, pheno, addcovar, intcovar, weights, maxit = 100L, tol = 1e-6, qr_tol = 1e-12, eta_max = 30.0) {
    .Call(`_qtl2_scancoefSE_binary_intcovar`, genoprobs, pheno, addcovar, intcovar, weights, maxit, tol, qr_tol, eta_max)
}

scancoef_hk_addcovar <- function(genoprobs, pheno, addcovar, weights, tol = 1e-12) {
    .Call(`_qtl2_scancoef_hk_addcovar`, genoprobs, pheno, addcovar, weights, tol)
}

scancoef_hk_intcovar <- function(genoprobs, pheno, addcovar, intcovar, weights, tol = 1e-12) {
    .Call(`_qtl2_scancoef_hk_intcovar`, genoprobs, pheno, addcovar, intcovar, weights, tol)
}

scancoefSE_hk_addcovar <- function(genoprobs, pheno, addcovar, weights, tol = 1e-12) {
    .Call(`_qtl2_scancoefSE_hk_addcovar`, genoprobs, pheno, addcovar, weights, tol)
}

scancoefSE_hk_intcovar <- function(genoprobs, pheno, addcovar, intcovar, weights, tol = 1e-12) {
    .Call(`_qtl2_scancoefSE_hk_intcovar`, genoprobs, pheno, addcovar, intcovar, weights, tol)
}

scancoef_pg_addcovar <- function(genoprobs, pheno, addcovar, eigenvec, weights, tol = 1e-12) {
    .Call(`_qtl2_scancoef_pg_addcovar`, genoprobs, pheno, addcovar, eigenvec, weights, tol)
}

scancoef_pg_intcovar <- function(genoprobs, pheno, addcovar, intcovar, eigenvec, weights, tol = 1e-12) {
    .Call(`_qtl2_scancoef_pg_intcovar`, genoprobs, pheno, addcovar, intcovar, eigenvec, weights, tol)
}

scancoefSE_pg_addcovar <- function(genoprobs, pheno, addcovar, eigenvec, weights, tol = 1e-12) {
    .Call(`_qtl2_scancoefSE_pg_addcovar`, genoprobs, pheno, addcovar, eigenvec, weights, tol)
}

scancoefSE_pg_intcovar <- function(genoprobs, pheno, addcovar, intcovar, eigenvec, weights, tol = 1e-12) {
    .Call(`_qtl2_scancoefSE_pg_intcovar`, genoprobs, pheno, addcovar, intcovar, eigenvec, weights, tol)
}

.calc_sdp <- function(geno) {
    .Call(`_qtl2_calc_sdp`, geno)
}

.invert_sdp <- function(sdp, n_str) {
    .Call(`_qtl2_invert_sdp`, sdp, n_str)
}

.alleleprob_to_snpprob <- function(alleleprob, sdp, interval, on_map) {
    .Call(`_qtl2_alleleprob_to_snpprob`, alleleprob, sdp, interval, on_map)
}

genocol_to_snpcol <- function(n_str, sdp) {
    .Call(`_qtl2_genocol_to_snpcol`, n_str, sdp)
}

.genoprob_to_snpprob <- function(genoprob, sdp, interval, on_map) {
    .Call(`_qtl2_genoprob_to_snpprob`, genoprob, sdp, interval, on_map)
}

Xgenocol_to_snpcol <- function(n_str, sdp) {
    .Call(`_qtl2_Xgenocol_to_snpcol`, n_str, sdp)
}

.Xgenoprob_to_snpprob <- function(genoprob, sdp, interval, on_map) {
    .Call(`_qtl2_Xgenoprob_to_snpprob`, genoprob, sdp, interval, on_map)
}

test_init <- function(crosstype, true_gen, is_x_chr, is_female, cross_info) {
    .Call(`_qtl2_test_init`, crosstype, true_gen, is_x_chr, is_female, cross_info)
}

test_emit <- function(crosstype, obs_gen, true_gen, error_prob, founder_geno, is_x_chr, is_female, cross_info) {
    .Call(`_qtl2_test_emit`, crosstype, obs_gen, true_gen, error_prob, founder_geno, is_x_chr, is_female, cross_info)
}

test_step <- function(crosstype, gen_left, gen_right, rec_frac, is_x_chr, is_female, cross_info) {
    .Call(`_qtl2_test_step`, crosstype, gen_left, gen_right, rec_frac, is_x_chr, is_female, cross_info)
}

test_check_geno <- function(crosstype, gen, is_observed_value, is_x_chr, is_female, cross_info) {
    .Call(`_qtl2_test_check_geno`, crosstype, gen, is_observed_value, is_x_chr, is_female, cross_info)
}

test_possible_gen <- function(crosstype, is_x_chr, is_female, cross_info) {
    .Call(`_qtl2_test_possible_gen`, crosstype, is_x_chr, is_female, cross_info)
}

test_ngen <- function(crosstype, is_x_chr) {
    .Call(`_qtl2_test_ngen`, crosstype, is_x_chr)
}

test_nrec <- function(crosstype, gen_left, gen_right, is_x_chr, is_female, cross_info) {
    .Call(`_qtl2_test_nrec`, crosstype, gen_left, gen_right, is_x_chr, is_female, cross_info)
}

test_founder_geno_values <- function(crosstype, founder_geno) {
    .Call(`_qtl2_test_founder_geno_values`, crosstype, founder_geno)
}

need_founder_geno <- function(crosstype) {
    .Call(`_qtl2_need_founder_geno`, crosstype)
}

check_founder_geno_size <- function(crosstype, founder_geno, n_markers) {
    .Call(`_qtl2_check_founder_geno_size`, crosstype, founder_geno, n_markers)
}

test_emitmatrix <- function(crosstype, error_prob, max_obsgeno, founder_geno, is_x_chr, is_female, cross_info) {
    .Call(`_qtl2_test_emitmatrix`, crosstype, error_prob, max_obsgeno, founder_geno, is_x_chr, is_female, cross_info)
}

test_stepmatrix <- function(crosstype, rec_frac, is_x_chr, is_female, cross_info) {
    .Call(`_qtl2_test_stepmatrix`, crosstype, rec_frac, is_x_chr, is_female, cross_info)
}

test_initvector <- function(crosstype, is_x_chr, is_female, cross_info) {
    .Call(`_qtl2_test_initvector`, crosstype, is_x_chr, is_female, cross_info)
}

