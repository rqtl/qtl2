// functions to test basic HMM things from R

#include <Rcpp.h>
#include "cross.h"

using namespace Rcpp;

// test init functions from R
// [[Rcpp::export]]
double test_init(const String& crosstype,
                 const int true_gen,
                 const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    return cross->init(true_gen, is_x_chr, is_female, cross_info);
}

// test emit functions from R
// [[Rcpp::export]]
double test_emit(const String& crosstype,
                 const int obs_gen, const int true_gen, const double error_prob,
                 const IntegerVector& founder_geno, const bool is_x_chr, 
                 const bool is_female, const IntegerVector& cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    return cross->emit(obs_gen, true_gen, error_prob,
                       founder_geno, is_x_chr, is_female, cross_info);
}


// test emit functions from R
// [[Rcpp::export]]
double test_step(const String& crosstype,
                 const int gen_left, const int gen_right, const double rec_frac,
                 const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    return cross->step(gen_left, gen_right, rec_frac, is_x_chr, is_female, cross_info);
}

// [[Rcpp::export]]
bool test_check_geno(const String& crosstype, const int gen, const bool is_observed_value,
                     const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    return cross->check_geno(gen, is_observed_value, is_x_chr, is_female, cross_info);
}

// [[Rcpp::export]]
IntegerVector test_possible_gen(const String& crosstype,
                                const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    return wrap(cross->possible_gen(is_x_chr, is_female, cross_info));
}

// [[Rcpp::export]]
int test_ngen(const String& crosstype, const bool is_x_chr)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    return cross->ngen(is_x_chr);
}

// [[Rcpp::export]]
double test_nrec(const String& crosstype, const int gen_left, const int gen_right,
                 const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    return cross->nrec(gen_left, gen_right, is_x_chr, is_female, cross_info);
}
