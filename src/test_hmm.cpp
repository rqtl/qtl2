// functions to test basic HMM things from R

#include <Rcpp.h>
#include "cross.h"
#include "test_hmm.h"

using namespace Rcpp;

// test init functions from R
// [[Rcpp::export]]
double test_init(const String& crosstype,
                 const int true_gen,
                 const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    double result = cross->init(true_gen, is_x_chr, is_female, cross_info);
    delete cross;
    return result;
}

// test emit functions from R
// [[Rcpp::export]]
double test_emit(const String& crosstype,
                 const int obs_gen, const int true_gen, const double error_prob,
                 const IntegerVector& founder_geno, const bool is_x_chr,
                 const bool is_female, const IntegerVector& cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    double result = cross->emit(obs_gen, true_gen, error_prob,
                                founder_geno, is_x_chr, is_female, cross_info);
    delete cross;
    return result;
}


// test emit functions from R
// [[Rcpp::export]]
double test_step(const String& crosstype,
                 const int gen_left, const int gen_right, const double rec_frac,
                 const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    double result = cross->step(gen_left, gen_right, rec_frac, is_x_chr, is_female, cross_info);
    delete cross;
    return result;
}

// [[Rcpp::export]]
bool test_check_geno(const String& crosstype, const int gen, const bool is_observed_value,
                     const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    bool result = cross->check_geno(gen, is_observed_value, is_x_chr, is_female, cross_info);
    delete cross;
    return result;
}

// [[Rcpp::export]]
IntegerVector test_possible_gen(const String& crosstype,
                                const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    IntegerVector result = wrap(cross->possible_gen(is_x_chr, is_female, cross_info));
    delete cross;
    return result;

}

// [[Rcpp::export]]
int test_ngen(const String& crosstype, const bool is_x_chr)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    int result = cross->ngen(is_x_chr);
    delete cross;
    return result;
}

// [[Rcpp::export]]
double test_nrec(const String& crosstype, const int gen_left, const int gen_right,
                 const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    double result = cross->nrec(gen_left, gen_right, is_x_chr, is_female, cross_info);
    delete cross;
    return result;
}

// [[Rcpp::export]]
bool test_founder_geno_values(const String& crosstype, const IntegerMatrix& founder_geno)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    bool result = cross->check_founder_geno_values(founder_geno);
    delete cross;
    return result;
}

// [[Rcpp::export]]
bool need_founder_geno(const String& crosstype)
{
    QTLCross* cross = QTLCross::Create(crosstype);

    bool result = cross->need_founder_geno();
    delete cross;
    return result;
}
