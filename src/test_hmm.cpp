#include <Rcpp.h>
#include <vector>
#include "cross.h"

using namespace Rcpp;

// test init functions from R
// [[Rcpp::export]]
double test_init(String crosstype,
                 int true_gen,
                 bool is_X_chr, bool is_female, IntegerVector cross_info,
                 bool phase_known=false)
{
    Cross* cross = Cross::Create(crosstype);

    std::vector<int> ci(as<std::vector<int> >(cross_info));

    if(phase_known) return cross->initPK(true_gen, is_X_chr, is_female, ci);
    return cross->init(true_gen, is_X_chr, is_female, ci);
}

// test emit functions from R
// [[Rcpp::export]]
double test_emit(String crosstype,
                 int obs_gen, int true_gen, double error_prob,
                 bool is_X_chr, bool is_female, IntegerVector cross_info,
                 bool phase_known=false)
{
    Cross* cross = Cross::Create(crosstype);

    std::vector<int> ci(as<std::vector<int> >(cross_info));

    if(phase_known)
        return cross->emitPK(obs_gen, true_gen, error_prob,
                             is_X_chr, is_female, ci);

    return cross->emit(obs_gen, true_gen, error_prob,
                       is_X_chr, is_female, ci);
}


// test emit functions from R
// [[Rcpp::export]]
double test_step(String crosstype,
                 int gen_left, int gen_right, double rec_frac,
                 bool is_X_chr, bool is_female, IntegerVector cross_info,
                 bool phase_known=false)
{
    Cross* cross = Cross::Create(crosstype);

    std::vector<int> ci(as<std::vector<int> >(cross_info));

    if(phase_known)
        return cross->stepPK(gen_left, gen_right, rec_frac, is_X_chr, is_female, ci);
    return cross->step(gen_left, gen_right, rec_frac, is_X_chr, is_female, ci);
}

// [[Rcpp::export]]
bool test_check_geno(String crosstype, int gen, bool is_observed_value,
                     bool is_X_chr, bool is_female, IntegerVector cross_info,
                     bool phase_known=false)
{
    Cross* cross = Cross::Create(crosstype);

    std::vector<int> ci(as<std::vector<int> >(cross_info));

    if(phase_known)
        return cross->check_genoPK(gen, is_observed_value, is_X_chr, is_female, ci);
    return cross->check_geno(gen, is_observed_value, is_X_chr, is_female, ci);
}

// [[Rcpp::export]]
IntegerVector test_geno(String crosstype,
                        bool is_X_chr, bool is_female, IntegerVector cross_info,
                        bool phase_known=false)
{
    Cross* cross = Cross::Create(crosstype);

    std::vector<int> ci(as<std::vector<int> >(cross_info));

    if(phase_known)
        return wrap(cross->genoPK(is_X_chr, is_female, ci));
    return wrap(cross->geno(is_X_chr, is_female, ci));
}

// [[Rcpp::export]]
IntegerVector test_allgeno(String crosstype, bool is_X_chr, bool phase_known=false)
{
    Cross* cross = Cross::Create(crosstype);

    if(phase_known)
        return wrap(cross->allgenoPK(is_X_chr));
    return wrap(cross->allgeno(is_X_chr));
}

// [[Rcpp::export]]
double test_nrec(String crosstype, int gen_left, int gen_right,
                 bool is_X_chr, bool is_female, IntegerVector cross_info)
{
    Cross* cross = Cross::Create(crosstype);

    std::vector<int> ci(as<std::vector<int> >(cross_info));

    return cross->nrec(gen_left, gen_right, is_X_chr, is_female, ci);
}
