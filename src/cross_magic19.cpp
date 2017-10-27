// 16-way MAGIC (RIL by selfing) QTLCross class (for HMM)
//
// See Kover et al (2009) PLOS Genet 5: e1000551 doi:10.1371/journal.pgen.1000551
// full diallel (in both directions) among 19 strains
// then 3 generations of random mating (to F4) - 384 families
// then selfed for 6 generations

#include "cross_magic19.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_util.h"
#include "cross_do_util.h"
#include "r_message.h" // defines RQTL2_NODEBUG and r_message()

enum gen {A=1, H=2, B=3, notA=5, notB=4};

const bool MAGIC19::check_geno(const int gen, const bool is_observed_value,
                               const bool is_x_chr, const bool is_female,
                               const IntegerVector& cross_info)
{
    // allow any value 0-5 for observed
    if(is_observed_value) {
        if(gen==0 || gen==A || gen==H || gen==B ||
           gen==notA || gen==notB) return true;
        else return false;
    }

    const int n_geno = 19;

    if(gen>= 1 && gen <= n_geno) return true;

    return false; // otherwise a problem
}

const double MAGIC19::init(const int true_gen,
                           const bool is_x_chr, const bool is_female,
                           const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    return -log(19.0);
}

const double MAGIC19::emit(const int obs_gen, const int true_gen, const double error_prob,
                           const IntegerVector& founder_geno, const bool is_x_chr,
                           const bool is_female, const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(obs_gen==0) return 0.0; // missing

    int f = founder_geno[true_gen-1]; // founder allele
    if(f!=1 && f!=3) return 0.0;      // founder missing -> no information

    if(f == obs_gen) return log(1.0 - error_prob);

    return log(error_prob); // genotyping error
}


const double MAGIC19::step(const int gen_left, const int gen_right, const double rec_frac,
                           const bool is_x_chr, const bool is_female,
                           const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(gen_left == gen_right)
        return log(19.0 - 52.0*rec_frac + 54.0*rec_frac*rec_frac - 18.0*rec_frac*rec_frac*rec_frac) -
            log(1.0+2.0*rec_frac) - log(19.0);
    else
        return log(rec_frac) + log(90.0 - 54.0*rec_frac + 18.0*rec_frac*rec_frac) -
            log(1.0+2.0*rec_frac) - log(19.0) - log(18.0);
}

const IntegerVector MAGIC19::possible_gen(const bool is_x_chr, const bool is_female,
                                          const IntegerVector& cross_info)
{
    int n_geno = 19;
    IntegerVector result(n_geno);

    for(int i=0; i<n_geno; i++) result[i] = i+1;
    return result;
}

const int MAGIC19::ngen(const bool is_x_chr)
{
    return 19;
}

const int MAGIC19::nalleles()
{
    return 19;
}


// check that founder genotype data has correct no. founders and markers
const bool MAGIC19::check_founder_geno_size(const IntegerMatrix& founder_geno, const int n_markers)
{
    bool result=true;

    const int fg_mar = founder_geno.cols();
    const int fg_f   = founder_geno.rows();

    if(fg_mar != n_markers) {
        result = false;
        r_message("founder_geno has incorrect number of markers");
    }

    if(fg_f != 19) {
        result = false;
        r_message("founder_geno should have 19 founders");
    }

    return result;
}

// check that founder genotype data has correct values
const bool MAGIC19::check_founder_geno_values(const IntegerMatrix& founder_geno)
{
    const int fg_mar = founder_geno.cols();
    const int fg_f   = founder_geno.rows();

    for(int f=0; f<fg_f; f++) {
        for(int mar=0; mar<fg_mar; mar++) {
            int fg = founder_geno(f,mar);
            if(fg != 0 && fg != 1 && fg != 3) {
                // at least one invalid value
                r_message("founder_geno contains invalid values; should be in {0, 1, 3}");
                return false;
            }
        }
    }

    return true;
}

const bool MAGIC19::need_founder_geno()
{
    return true;
}

// geno_names from allele names
const std::vector<std::string> MAGIC19::geno_names(const std::vector<std::string> alleles,
                                                   const bool is_x_chr)
{
    if(alleles.size() < 19)
        throw std::range_error("alleles must have length 19");

    const int n_alleles = 19;

    std::vector<std::string> result(n_alleles);

    for(int i=0; i<n_alleles; i++)
        result[i] = alleles[i] + alleles[i];

    return result;
}


const int MAGIC19::nrec(const int gen_left, const int gen_right,
                        const bool is_x_chr, const bool is_female,
                        const Rcpp::IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(gen_left == gen_right) return 0;
    else return 1;
}

const double MAGIC19::est_rec_frac(const Rcpp::NumericVector& gamma, const bool is_x_chr,
                                   const Rcpp::IntegerMatrix& cross_info, const int n_gen)
{
    double R = QTLCross::est_rec_frac(gamma, is_x_chr, cross_info, n_gen);

    // inverse of R = (90-54r+18r^2)/19/(1+2r)
    double A = pow(3.0, -4.5) * sqrt(2475.0-304.0*R) * (18.0-19.0*R)/4.0;
    double B = (18.0-19.0*R)/12.0;
    double C = -pow(A+B, 1.0/3.0);

    double D = 18.0-19.0*R;
    double E = sqrt(2475.0 - 304.0*R)*(18.0-19.0*R)/4.0/pow(3,4.5);
    double F = (18.0-19.0*R)/12.0;
    double G = pow(E+F, 1.0/3.0) * 27.0;

    return C + D/G + 1.0;
}

// check whether X chr can be handled
const bool MAGIC19::check_handle_x_chr(const bool any_x_chr)
{
    if(any_x_chr) {
        r_message("X chr ignored for RIL by selfing.");
        return false;
    }

    return true; // most crosses can handle the X chr
}
