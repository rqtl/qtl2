// general qtlcross class
//
// see cross.cpp for info on how to add a new cross type

#ifndef CROSS_H
#define CROSS_H

#include <Rcpp.h>

using namespace Rcpp;

class QTLCross
{

public:
    String type;

    String phase_known_type;

    static QTLCross* Create(const String type);

    virtual const bool check_geno(const int gen, const bool is_observed_value,
                                  const bool& is_x_chr, const bool& is_female,
                                  const IntegerVector& cross_info)
    {
        if(is_observed_value && gen==0) return true;
        if(gen==1 || gen==2) return true;

        return false;
    }

    virtual const double init(const int true_gen,
                              const bool& is_x_chr, const bool& is_female,
                              const IntegerVector& cross_info)
    {
        #ifdef DEBUG
        if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
            throw std::range_error("genotype value not allowed");
        #endif

        return -log(2.0);
    }

    virtual const double emit(const int obs_gen, const int true_gen, const double error_prob,
                              const bool& is_x_chr, const bool& is_female,
                              const IntegerVector& cross_info)
    {
        #ifdef DEBUG
        if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
            throw std::range_error("genotype value not allowed");
        #endif


        if(obs_gen==0 || !check_geno(obs_gen, true, is_x_chr, is_female, cross_info))
            return 0.0; // missing or invalid

        if(obs_gen == true_gen) return log(1.0 - error_prob);
        else return log(error_prob);

    }

    virtual const double step(const int gen_left, const int gen_right, const double rec_frac,
                              const bool& is_x_chr, const bool& is_female,
                              const IntegerVector& cross_info)
    {
        #ifdef DEBUG
        if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
           !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
            throw std::range_error("genotype value not allowed");
        #endif

        if(gen_left == gen_right) return log(1.0-rec_frac);
        else return log(rec_frac);
    }

    virtual const int ngen(const bool& is_x_chr)
    {
        return 2;
    }

    virtual const IntegerVector possible_gen(const bool& is_x_chr, const bool& is_female,
                                             const IntegerVector& cross_info)
    {
        int ng = ngen(is_x_chr);
        IntegerVector x(ng);
        for(int i=0; i<ng; i++) x[i] = i+1;
        return x;
    }

    virtual const double nrec(const int gen_left, const int gen_right,
                              const bool& is_x_chr, const bool& is_female,
                              const IntegerVector& cross_info)
    {
        #ifdef DEBUG
        if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
           !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
            throw std::range_error("genotype value not allowed");
        #endif

        if(gen_left == gen_right) return 0.0;
        else return 1.0;
    }

    virtual const double est_rec_frac(const NumericMatrix& gamma, const bool& is_x_chr)
    {
        int n_gen = gamma.rows();
        int n_gen_sq = n_gen*n_gen;

        double denom = 0.0;
        for(int i=0; i<n_gen_sq; i++) denom += gamma[i];

        double diagsum = 0.0;
        for(int i=0; i<n_gen; i++) diagsum += gamma(i,i);

        return 1.0 - diagsum/denom;

    }

    // the following is for checking whether a crosstype is supported from R
    // (some classes, like f2pk, aren't appropriate on the R side
    virtual const bool crosstype_supported()
    {
        return true;
    }

};

#endif // CROSS.H
