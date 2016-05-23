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
    Rcpp::String crosstype;

    Rcpp::String phase_known_crosstype;

    static QTLCross* Create(const Rcpp::String& crosstype);

    virtual ~QTLCross(){};

    virtual const bool check_geno(const int gen, const bool is_observed_value,
                                  const bool is_x_chr, const bool is_female,
                                  const Rcpp::IntegerVector& cross_info)
    {
        if(is_observed_value && gen==0) return true;
        if(gen==1 || gen==2) return true;

        return false;
    }

    virtual const double init(const int true_gen,
                              const bool is_x_chr, const bool is_female,
                              const Rcpp::IntegerVector& cross_info)
    {
        #ifndef NDEBUG
        if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
            throw std::range_error("genotype value not allowed");
        #endif

        return -log(2.0);
    }

    virtual const double emit(const int obs_gen, const int true_gen, const double error_prob,
                              const Rcpp::IntegerVector& founder_geno, const bool is_x_chr,
                              const bool is_female, const Rcpp::IntegerVector& cross_info)
    {
        #ifndef NDEBUG
        if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
            throw std::range_error("genotype value not allowed");
        #endif


        if(obs_gen==0 || !check_geno(obs_gen, true, is_x_chr, is_female, cross_info))
            return 0.0; // missing or invalid

        if(obs_gen == true_gen) return log(1.0 - error_prob);
        else return log(error_prob);

    }

    virtual const double step(const int gen_left, const int gen_right, const double rec_frac,
                              const bool is_x_chr, const bool is_female,
                              const Rcpp::IntegerVector& cross_info)
    {
        #ifndef NDEBUG
        if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
           !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
            throw std::range_error("genotype value not allowed");
        #endif

        if(gen_left == gen_right) return log(1.0-rec_frac);
        else return log(rec_frac);
    }

    virtual const int ngen(const bool is_x_chr)
    {
        return 2;
    }

    virtual const int nalleles()
    {
        return 2;
    }

    virtual const Rcpp::IntegerVector possible_gen(const bool is_x_chr, const bool is_female,
                                                   const Rcpp::IntegerVector& cross_info)
    {
        int ng = ngen(is_x_chr);
        Rcpp::IntegerVector x(ng);
        for(int i=0; i<ng; i++) x[i] = i+1;
        return x;
    }

    virtual const double nrec(const int gen_left, const int gen_right,
                              const bool is_x_chr, const bool is_female,
                              const Rcpp::IntegerVector& cross_info)
    {
        #ifndef NDEBUG
        if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
           !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
            throw std::range_error("genotype value not allowed");
        #endif

        if(gen_left == gen_right) return 0.0;
        else return 1.0;
    }

    virtual const double est_rec_frac(const Rcpp::NumericVector& gamma, const bool is_x_chr,
                                      const Rcpp::IntegerMatrix& cross_info, const int n_gen)
    {
        int n_ind = cross_info.cols();
        int n_gen_sq = n_gen*n_gen;

        double diagsum=0.0;
        for(int ind=0, offset=0; ind<n_ind; ind++, offset += n_gen_sq)
            for(int i=0; i<n_gen; i++) diagsum += gamma[offset+i*n_gen+i];

        double result = 1.0 - diagsum/(double)n_ind;
        if(result < 0.0) result = 0.0;

        return result;
    }

    // the following is for checking whether a crosstype is supported from R
    // (some classes, like f2pk, aren't appropriate on the R side
    virtual const bool crosstype_supported()
    {
        return true;
    }

    // check that founder genotype data has correct no. founders and markers
    // (for crosses with no founder_geno, just return true)
    virtual const bool check_founder_geno_size(const Rcpp::IntegerMatrix& founder_geno, const int n_markers)
    {
        return true;
    }

    // check that founder genotype data contains correct values
    // (for crosses with no founder_geno, just return true)
    virtual const bool check_founder_geno_values(const Rcpp::IntegerMatrix& founder_geno)
    {
        return true;
    }

    // matrix to convert genotype probabilities to allele probabilities
    // if no conversion necessary, it returns a matrix with 0 rows and 0 cols
    virtual const Rcpp::NumericMatrix geno2allele_matrix(const bool is_x_chr)
    {
        return Rcpp::NumericMatrix(0,0);
    }


    // check that cross_info conforms to expectation
    virtual const bool check_crossinfo(const Rcpp::IntegerMatrix& cross_info, const bool any_x_chr)
    {
        //const unsigned int n_col = cross_info.cols();
        //if(n_col > 0)
        //    REprintf("cross_info provided (with %d columns) but ignored for this cross type\n", n_col);

        return true; // don't call it an error
    }

    // check that sex conforms to expectation
    virtual const bool check_is_female_vector(const Rcpp::LogicalVector& is_female, const bool any_x_chr)
    {
        //const unsigned int n = is_female.size();
        //if(n > 0)
        //    REprintf("is_female provided but ignored for this cross type\n");

        return true; // don't call it an error
    }

    // check whether X chr can be handled
    virtual const bool check_handle_x_chr(const bool any_x_chr)
    {
        return true; // most crosses can handle the X chr
    }

    // need founder_geno?
    virtual const bool need_founder_geno()
    {
        return false; // most crosses don't need founder_geno
    }

    // X chromosome covariates
    virtual const Rcpp::NumericMatrix get_x_covar(const Rcpp::LogicalVector& is_female, const Rcpp::IntegerMatrix& cross_info)
    {
        const unsigned int n_ind = is_female.size();
        unsigned int n_female=0;
        for(unsigned int i=0; i<n_ind; i++)
            if(is_female[i]) ++n_female;

        if(n_female==0 || n_female==n_ind) // all male or all female
            return Rcpp::NumericMatrix(n_ind,0);

        // some males and some females; return single-column matrix with sex indicators
        Rcpp::NumericMatrix result(n_ind,1);
        for(unsigned int i=0; i<n_ind; i++) {
            if(is_female[i]) result(i,0) = 0.0;
            else result(i,0) = 1.0;
        }
        colnames(result) = CharacterVector::create("sex");

        return result;
    }

    // genotype names from allele names
    // default version is A,B -> AA,BB
    virtual const std::vector<std::string> geno_names(const std::vector<std::string> alleles,
                                                      const bool is_x_chr)
    {
        if(alleles.size() < 2)
            throw std::range_error("alleles must have length 2");

        std::vector<std::string> result(2);
        for(int i=0; i<2; i++) {
            result[i] = alleles[i] + alleles[i];
        }
        return result;
    }


    // calculate a vector of emission matrices
    virtual const std::vector<Rcpp::NumericMatrix> calc_emitmatrix(const double error_prob,
                                                                   const Rcpp::IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                                                                   const bool is_x_chr, const bool is_female,
                                                                   const Rcpp::IntegerVector& cross_info)
    {
        Rcpp::IntegerVector gen = possible_gen(is_x_chr, is_female, cross_info);
        const unsigned int n_true_gen = gen.size();
        const unsigned int n_obs_gen = 6;

        const unsigned int n_markers = founder_geno.cols();

        std::vector<Rcpp::NumericMatrix> result;
        for(unsigned int i=0; i<n_markers; i++) {
            Rcpp::NumericMatrix emitmatrix(n_obs_gen, n_true_gen);
            for(unsigned int obs_gen=0; obs_gen<n_obs_gen; obs_gen++) {
                for(unsigned int true_gen=0; true_gen<n_true_gen; true_gen++) {
                    emitmatrix(obs_gen,true_gen) = emit(obs_gen, gen[true_gen], error_prob,
                                                        founder_geno(_, i), is_x_chr, is_female, cross_info);
                }
            }
            result.push_back(emitmatrix);
        }

        return result;
    }

    // calculate a vector of transition matrices
    virtual const std::vector<Rcpp::NumericMatrix> calc_stepmatrix(const Rcpp::NumericVector rec_frac,
                                                                   const bool is_x_chr, const bool is_female,
                                                                   const Rcpp::IntegerVector& cross_info)
    {
        Rcpp::IntegerVector gen = possible_gen(is_x_chr, is_female, cross_info);
        const unsigned int n_gen = gen.size();

        const unsigned int n_intervals = rec_frac.size();

        std::vector<Rcpp::NumericMatrix> result;
        for(unsigned int i=0; i<n_intervals; i++) {
            Rcpp::NumericMatrix stepmatrix(n_gen, n_gen);
            for(unsigned int left=0; left<n_gen; left++) {
                for(unsigned int right=0; right<n_gen; right++) {
                    stepmatrix(left,right) = step(gen[left], gen[right], rec_frac[i],
                                                  is_x_chr, is_female, cross_info);
                }
            }
            result.push_back(stepmatrix);
        }

        return result;
    }

    // calculate init probabilities
    virtual const Rcpp::NumericVector calc_initvector(const bool is_x_chr, const bool is_female,
                                        const Rcpp::IntegerVector& cross_info)
    {
        Rcpp::IntegerVector gen = possible_gen(is_x_chr, is_female, cross_info);
        const unsigned int n_gen = gen.size();

        Rcpp::NumericMatrix result(n_gen);
        for(unsigned int g=0; g<n_gen; g++) {
            result[g] = init(gen[g], is_x_chr, is_female, cross_info);
        }

        return result;
    }

};

#endif // CROSS_H
