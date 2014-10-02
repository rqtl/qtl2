#include <Rcpp.h>
#include <vector>
#include "cross.h"

using namespace Rcpp;

// test init functions from R
// [[Rcpp::export]]
double test_init(String crosstype,
                 int true_geno,
                 bool is_X_chr, bool is_female, IntegerVector cross_info, 
                 bool phase_known=false)
{
    Cross* cross = Cross::Create(crosstype);
    
    std::vector<int> ci(as<std::vector<int> >(cross_info));

    if(phase_known) return cross->initPK(true_geno, is_X_chr, is_female, ci);
    return cross->init(true_geno, is_X_chr, is_female, ci);
}

// test emit functions from R
// [[Rcpp::export]]
double test_emit(String crosstype,
                 int obs_geno, int true_geno, double error_prob,
                 bool is_X_chr, bool is_female, IntegerVector cross_info,
                 bool phase_known=false)
{
    Cross* cross = Cross::Create(crosstype);
    
    std::vector<int> ci(as<std::vector<int> >(cross_info));

    if(phase_known)
        return cross->emitPK(obs_geno, true_geno, error_prob,
                             is_X_chr, is_female, ci);
    
    return cross->emit(obs_geno, true_geno, error_prob,
                       is_X_chr, is_female, ci);
}


// test emit functions from R
// [[Rcpp::export]]
double test_step(String crosstype,
                 int geno_left, int geno_right, double rec_frac,
                 bool is_X_chr, bool is_female, IntegerVector cross_info,
                 bool phase_known=false)
{
    Cross* cross = Cross::Create(crosstype);
    
    std::vector<int> ci(as<std::vector<int> >(cross_info));

    if(phase_known)
        return cross->stepPK(geno_left, geno_right, rec_frac, is_X_chr, is_female, ci);
    return cross->step(geno_left, geno_right, rec_frac, is_X_chr, is_female, ci);
}
