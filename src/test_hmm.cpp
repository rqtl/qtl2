// functions to test basic HMM things from R

#include <Rcpp.h>
#include "cross.h"

using namespace Rcpp;

// test init functions from R
// [[Rcpp::export]]
double test_init(String crosstype,
                 int true_gen,
                 bool is_X_chr, bool is_female, IntegerVector cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    return cross->init(true_gen, is_X_chr, is_female, cross_info);
}

// test emit functions from R
// [[Rcpp::export]]
double test_emit(String crosstype,
                 int obs_gen, int true_gen, double error_prob,
                 bool is_X_chr, bool is_female, IntegerVector cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    return cross->emit(obs_gen, true_gen, error_prob,
                       is_X_chr, is_female, cross_info);
}


// test emit functions from R
// [[Rcpp::export]]
double test_step(String crosstype,
                 int gen_left, int gen_right, double rec_frac,
                 bool is_X_chr, bool is_female, IntegerVector cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    return cross->step(gen_left, gen_right, rec_frac, is_X_chr, is_female, cross_info);
}

// [[Rcpp::export]]
bool test_check_geno(String crosstype, int gen, bool is_observed_value,
                     bool is_X_chr, bool is_female, IntegerVector cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    return cross->check_geno(gen, is_observed_value, is_X_chr, is_female, cross_info);
}

// [[Rcpp::export]]
IntegerVector test_possible_gen(String crosstype,
                                bool is_X_chr, bool is_female, IntegerVector cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    return wrap(cross->possible_gen(is_X_chr, is_female, cross_info));
}

// [[Rcpp::export]]
int test_ngen(String crosstype, bool is_X_chr)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    return cross->ngen(is_X_chr);
}

// [[Rcpp::export]]
double test_nrec(String crosstype, int gen_left, int gen_right,
                 bool is_X_chr, bool is_female, IntegerVector cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    return cross->nrec(gen_left, gen_right, is_X_chr, is_female, cross_info);
}
