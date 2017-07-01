// phase-known AIL3 (3-way advanced intercross lines) QTLCross class (for HMM, in particular est.map)
// (assuming all F1 hybrids formed followed by random mating with large population)

#include "cross_ail3pk.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_util.h" // mpp_encode_alleles and mpp_decode_geno
#include "r_message.h"

enum gen {AA=1, AB=2, BB=3, notA=5, notB=4,
          A=1, H=2, B=3};

const bool AIL3PK::check_geno(const int gen, const bool is_observed_value,
                              const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    // allow any value 0-5 for observed
    if(is_observed_value) {
        if(gen==0 || gen==AA || gen==AB || gen==BB ||
           gen==notA || gen==notB) return true;
        else return false;
    }

    const int n_genoA = 9;
    const int n_genoX = 12;

    if(!is_x_chr || is_female) { // autosome or female X
        if(gen >= 1 && gen <= n_genoA) return true;
    }
    else { // male X
        if(gen > n_genoA && gen <= n_genoX) return true;
    }

    return false; // invalid
}


const double AIL3PK::init(const int true_gen,
                          const bool is_x_chr, const bool is_female,
                          const IntegerVector& cross_info)
{
    #ifndef NDEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(!is_x_chr || is_female) { // autosome or female X
        return -log(9.0);
    }
    else { // male X
        return -log(3.0);
    }
}

const double AIL3PK::emit(const int obs_gen, const int true_gen, const double error_prob,
                          const IntegerVector& founder_geno, const bool is_x_chr,
                          const bool is_female, const IntegerVector& cross_info)
{

    #ifndef NDEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    const int n_geno = 9;

    if(obs_gen==0) return 0.0; // missing

    if(!is_x_chr || is_female) { // autosome or female X
        const IntegerVector true_alleles = mpp_decode_geno(true_gen, 3, true);
        int f1 = founder_geno[true_alleles[0]-1];
        int f2 = founder_geno[true_alleles[1]-1];

        // treat founder hets as missing
        if(f1==2) f1 = 0;
        if(f2==2) f2 = 0;

        // neither founder alleles observed
        if(f1 == 0 && f2 == 0) return 0.0;

        // one founder allele observed
        if(f1 == 0 || f2 == 0) {

            switch(std::max(f1, f2)) {
            case H: return 0.0; // het compatible with either founder allele
            case A:
                switch(obs_gen) {
                case A: case notB: return log(1.0-error_prob);
                case B: case notA: return log(error_prob);
                case H: return 0.0;
                }
            case B:
                switch(obs_gen) {
                case B: case notA: return log(1.0-error_prob);
                case A: case notB: return log(error_prob);
                case H: return 0.0;
                }
            }
            return 0.0;
        }

        switch((f1+f2)/2) { // values 1, 2, 3
        case A:
            switch(obs_gen) {
            case A: return log(1.0-error_prob);
            case H: return log(error_prob/2.0);
            case B: return log(error_prob/2.0);
            case notA: return log(error_prob);
            case notB: return log(1.0-error_prob/2.0);
            }
        case H:
            switch(obs_gen) {
            case A: return log(error_prob/2.0);
            case H: return log(1.0-error_prob);
            case B: return log(error_prob/2.0);
            case notA: return log(1.0-error_prob/2.0);
            case notB: return log(1.0-error_prob/2.0);
            }
        case B:
            switch(obs_gen) {
            case B: return log(1.0-error_prob);
            case H: return log(error_prob/2.0);
            case A: return log(error_prob/2.0);
            case notB: return log(error_prob);
            case notA: return log(1.0-error_prob/2.0);
            }
        }
        return 0.0;
    }
    else { // male X
        const int founder_allele = founder_geno[(true_gen - n_geno) - 1];

        switch(founder_allele) {
        case A:
            switch(obs_gen) {
            case A: case notB: return log(1.0-error_prob);
            case B: case notA: return log(error_prob);
            }
        case B:
            switch(obs_gen) {
            case B: case notA: return log(1.0-error_prob);
            case A: case notB: return log(error_prob);
            }
        }
        return(0.0);
    }

    return NA_REAL; // shouldn't get here
}


const double AIL3PK::step(const int gen_left, const int gen_right, const double rec_frac,
                          const bool is_x_chr, const bool is_female,
                          const IntegerVector& cross_info)
{
    #ifndef NDEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    const int n_gen = cross_info[0]; // number of generations

    if(is_x_chr && !is_female) { // male X
        const double z = sqrt((1.0-rec_frac)*(9.0-rec_frac));
        const double pAA = (1.0-rec_frac)/3.0 * ( (1.0-rec_frac+z)/(2.0*z) * pow((1.0-rec_frac-z)/4.0, (double)(n_gen-2)) +
                                                  (-1.0+rec_frac+z)/(2.0*z) * pow((1.0-rec_frac+z)/4.0, (double)(n_gen-2))) +
            (2.0-rec_frac)/(2.0*3) * ( (1.0-rec_frac-z)/2.0 * (1.0-rec_frac+z)/(2.0*z) * pow((1.0-rec_frac-z)/4.0,(double)(n_gen-2)) +
                                       (1.0-rec_frac+z)/2.0 * (-1.0+rec_frac+z)/(2.0*z) * pow((1.0-rec_frac+z)/4.0,(double)(n_gen-2))) +
            ( (rec_frac*rec_frac + rec_frac*(z-5.0))/(9*(3.0+rec_frac+z)) * (1.0-rec_frac+z)/(2.0*z) * pow((1.0-rec_frac-z)/4.0,(double)(n_gen-2)) +
              (rec_frac*rec_frac - rec_frac*(z+5.0))/(9*(3.0+rec_frac-z)) * (-1.0+rec_frac+z)/(2.0*z) * pow((1.0-rec_frac+z)/4.0,(double)(n_gen-2)) + 1.0/9.0);
        const double R = (1.0-3.0*pAA);

        if(gen_left == gen_right) return log1p(-R);
        return log(R) - log(2.0);
    }
    else { // autosome or female X
        double pAA;
        if(!is_x_chr) { // autosome
            pAA = (1.0 - (-2.0 + 3.0*rec_frac)*pow(1.0 - rec_frac, (double)(n_gen-2)))/9.0;
        }
        else { // female X
            const double z = sqrt((1.0-rec_frac)*(9.0-rec_frac));

            pAA = (1.0-rec_frac)/3.0 * ( (-1.0/z)*pow((1.0-rec_frac-z)/4.0,(double)(n_gen-2)) +
                                         (1.0/z)*pow((1.0-rec_frac+z)/4.0,(double)(n_gen-2))) +
                (2.0-rec_frac)/6.0 * ( (1.0-rec_frac-z)/2.0 * (-1.0/z)*pow((1.0-rec_frac-z)/4.0,(double)(n_gen-2)) +
                                       (1.0-rec_frac+z)/2.0 *  (1.0/z)*pow((1.0-rec_frac+z)/4.0,(double)(n_gen-2))) +
                ( (rec_frac*rec_frac + rec_frac*(z-5.0))/(9.0*(3.0+rec_frac+z)) * (-1.0/z)*pow((1.0-rec_frac-z)/4.0,(double)(n_gen-2)) +
                  (rec_frac*rec_frac - rec_frac*(z+5.0))/(9.0*(3.0+rec_frac-z)) * (1.0/z)*pow((1.0-rec_frac+z)/4.0,(double)(n_gen-2))  + 1.0/9.0);
        }
        double R = (1.0 - 3.0*pAA);

        const IntegerVector alleles_left = mpp_decode_geno(gen_left, 3, true);
        const IntegerVector alleles_right = mpp_decode_geno(gen_right, 3, true);

        if(alleles_left[0] == alleles_left[1]) { // homozygous
            if(alleles_right[0] == alleles_right[1]) { // homozygous
                if(alleles_left[0] == alleles_right[0]) { // AA -> AA
                    return 2.0*log(1.0-R);
                }
                else { // AA -> BB
                    return 2.0*(log(R) - log(2.0));
                }
            }
            else {
                if(alleles_left[0] == alleles_right[0] ||
                   alleles_left[0] == alleles_right[1]) { // AA -> AB
                    return log(1.0-R) + log(R) - log(2.0);
                }
                else { // AA -> BC
                    return 2.0*(log(R)-log(2.0));
                }
            }
        }
        else {
            if(alleles_right[0] == alleles_right[1]) { // homozygous
                if(alleles_left[0] == alleles_right[0] ||
                   alleles_left[1] == alleles_right[1]) { // AB -> AA
                    return log(1.0-R) + log(R) - log(2.0);
                }
                else { // AB -> CC
                    return 2.0*(log(R) - log(2.0));
                }
            }
            else { // both het
                if(alleles_left[0] == alleles_right[0] &&
                   alleles_left[1] == alleles_right[1]) { // AB -> AB
                    return 2.0*log(1.0-R);
                }
                else if(alleles_left[0] == alleles_right[0] ||
                        alleles_left[1] == alleles_right[1]) { // AB -> AC
                    return log(1.0-R) + log(R) - log(2.0);
                }
                else { // AB -> BC
                    return 2.0*(log(R) - log(2.0));
                }
            }
        }

    }

    return NA_REAL; // shouldn't get here
}

const IntegerVector AIL3PK::possible_gen(const bool is_x_chr, const bool is_female,
                                         const IntegerVector& cross_info)
{
    if(is_x_chr && !is_female) { // male X chromosome
        IntegerVector result = IntegerVector::create(10, 11, 12);
        return result;
    }
    else { // autosome or female X
        IntegerVector result = IntegerVector::create(1,2,3,4,5,6,7,8,9);
        return result;
    }
}

const int AIL3PK::ngen(const bool is_x_chr)
{
    if(is_x_chr) return 12;
    return 9;
}

const int AIL3PK::nalleles()
{
    return 3;
}

const int AIL3PK::nrec(const int gen_left, const int gen_right,
                       const bool is_x_chr, const bool is_female,
                       const IntegerVector& cross_info)
{
    #ifndef NDEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(is_x_chr && gen_left > 9 && gen_right > 9) { // male X
        if(gen_left == gen_right) return(0);
        else return(1);
    }

    // otherwise autosome or female X
    Rcpp::IntegerVector a_left = mpp_decode_geno(gen_left, 3, true);
    Rcpp::IntegerVector a_right = mpp_decode_geno(gen_right, 3, true);

    int result = 0;
    if(a_left[0] != a_right[0]) result++;
    if(a_left[1] != a_right[1]) result++;
    return result;
}

const double AIL3PK::est_rec_frac(const NumericVector& gamma, const bool is_x_chr,
                                  const IntegerMatrix& cross_info, const int n_gen)
{
    Rcpp::stop("est_map not yet available for 3-way AIL");

    return NA_REAL;
}


const NumericMatrix AIL3PK::geno2allele_matrix(const bool is_x_chr)
{
    if(is_x_chr) {
        NumericMatrix result(12,3);
        result(0,0) = 1.0;               // AA female
        result(1,0) = result(1,1) = 0.5; // AB female
        result(2,1) = 1.0;               // BB female
        result(3,0) = result(3,2) = 0.5; // AC female
        result(4,1) = result(4,2) = 0.5; // BC female
        result(5,2) = 1.0;               // CC female
        result(6,1) = result(6,0) = 0.5; // BA female
        result(7,2) = result(7,0) = 0.5; // CA female
        result(8,2) = result(8,1) = 0.5; // CB female

        result(9, 0) = 1.0; // AY male
        result(10,1) = 1.0; // BY male
        result(11,2) = 1.0; // CY male

        return result;
    }
    else {
        NumericMatrix result(9,3);
        result(0,0) = 1.0;               // AA
        result(1,0) = result(1,1) = 0.5; // AB
        result(2,1) = 1.0;               // BB
        result(3,0) = result(3,2) = 0.5; // AC
        result(4,1) = result(4,2) = 0.5; // BC
        result(5,2) = 1.0;               // CC
        result(6,1) = result(6,0) = 0.5; // BA
        result(7,2) = result(7,0) = 0.5; // CA
        result(8,2) = result(8,1) = 0.5; // CB

        return result;
    }
}

// check that sex conforms to expectation
const bool AIL3PK::check_is_female_vector(const LogicalVector& is_female, const bool any_x_chr)
{
    bool result = true;
    const int n = is_female.size();
    if(!any_x_chr) { // all autosomes
        if(n > 0) {
            // not needed here, but don't call this an error
            result = true;
        }
    }
    else { // X chr included
        if(n == 0) {
            result = false;
            r_message("is_female not provided, but needed to handle X chromosome");
        }
        else {
            int n_missing = 0;
            for(int i=0; i<n; i++)
                if(is_female[i] == NA_LOGICAL) ++n_missing;
            if(n_missing > 0) {
                result = false;
                r_message("is_female contains missing values (it shouldn't)");
            }
        }
    }
    return result;
}

// check that cross_info conforms to expectation
const bool AIL3PK::check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr)
{
    bool result = true;
    const int n_row = cross_info.rows();
    const int n_col = cross_info.cols();
    // single column with number of generations (needed no matter what; values should be >= 2)

    if(n_col != 1) {
        result = false;
        r_message("cross_info should have one column, with no. generations");
        return result;
    }

    int n_missing=0;
    int n_invalid=0;
    for(int i=0; i<n_row; i++) {
        if(cross_info[i] == NA_INTEGER) ++n_missing;
        else if(cross_info[i] < 2) ++n_invalid;
    }
    if(n_missing > 0) {
        result = false;
        r_message("cross_info has missing values (it shouldn't)");
    }
    if(n_invalid > 0) {
        result = false;
        r_message("cross_info has invalid values; no. generations should be >= 2");
    }

    return result;
}
