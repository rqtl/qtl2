// RI by sib-mating QTLCross class (for HMM)

#include "cross_risib.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "r_message.h" // defines RQTL2_NODEBUG and r_message()

enum gen_risib {AA=1, BB=2};

const double RISIB::init(const int true_gen,
                         const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    const bool forward_direction = (cross_info.size()<1 || cross_info[0] == 0); // AA female x BB male

    if(is_x_chr) {
        if(forward_direction) {
            if(true_gen == AA) return log(2.0)-log(3.0);
            if(true_gen == BB) return -log(3.0);
        }
        else {
            if(true_gen == BB) return log(2.0)-log(3.0);
            if(true_gen == AA) return -log(3.0);
        }
    }
    else { // autosome
        return -log(2.0);
    }

    return NA_REAL; // can't get here
}

const double RISIB::step(const int gen_left, const int gen_right, const double rec_frac,
                         const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(is_x_chr)  {
        const bool forward_direction = (cross_info.size()<1 || cross_info[0] == 0); // AA female x BB male

        if(forward_direction) {
            if(gen_left == AA) {
                if(gen_right == AA)
                    return log(1.0 + 2.0*rec_frac) - log(1.0 + 4.0*rec_frac);
                if(gen_right == BB)
                    return log(2.0*rec_frac) - log(1.0 + 4.0*rec_frac);
            }

            if(gen_left == BB) {
                if(gen_right == BB)
                    return -log(1.0 + 4.0*rec_frac);
                if(gen_right == AA)
                    return log(4.0*rec_frac) - log(1.0 + 4.0*rec_frac);
            }
        }
        else {
            if(gen_left == AA) {
                if(gen_right == AA)
                    return -log(1.0 + 4.0*rec_frac);
                if(gen_right == BB)
                    return log(4.0*rec_frac) - log(1.0 + 4.0*rec_frac);
            }

            if(gen_left == BB) {
                if(gen_right == BB)
                    return log(1.0 + 2.0*rec_frac) - log(1.0 + 4.0*rec_frac);
                if(gen_right == AA)
                    return log(2.0*rec_frac) - log(1.0 + 4.0*rec_frac);
            }
        }
    }
    else { // autosome
        const double R = 4.0*rec_frac/(1+6.0*rec_frac);

        if(gen_left == gen_right) return log(1.0-R);
        else return log(R);
    }

    return NA_REAL; // can't get here
}

const double RISIB::est_rec_frac(const NumericVector& gamma, const bool is_x_chr,
                                 const IntegerMatrix& cross_info, const int n_gen)
{
    if(is_x_chr) {
        int n_ind = cross_info.cols();
        int n_gen_sq = n_gen*n_gen;

        double denom = 0.0, sum00=0.0, sum01=0.0;
        for(int ind=0, offset=0; ind<n_ind; ind++, offset += n_gen_sq) {
            for(int i=0; i<n_gen_sq; i++) denom += gamma[offset+i];

            if(cross_info.rows()<0 || cross_info[ind]==0) // A x B direction
                sum00 += gamma[offset];
            else                   // B x A direction
                sum00 += gamma[offset+3];

            // recombinants
            sum01 += (gamma[offset+1]+gamma[offset+2]);
        }

        // the MLE is solution to a quadratic equation
        return (2.0*denom - sum00 - 3.0*sum01  -
                (sqrt(sum01*sum01 + (-2.0*sum00-4.0*denom)*sum01 + sum00*sum00 -
                      4.0*denom*sum00 + 4.0*denom*denom)))/8.0/(sum01+sum00-denom);

    }
    else {
        double R = QTLCross::est_rec_frac(gamma, is_x_chr, cross_info, n_gen);

        return R/(4.0-6.0*R);
    }
}

// check that cross_info conforms to expectation
const bool RISIB::check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr)
{
    bool result = true;
    const int n_row = cross_info.rows();
    const int n_col = cross_info.cols();

    if(!any_x_chr) { // all autosomes
        if(n_col > 0) {
            // not needed here, but don't call this an error
            result = true;
        }
    }
    else { // X chr included
        if(n_col == 0) {
            result = false;
            r_message("cross_info not provided, but needed to handle X chromosome");
        }
        else if(n_col > 1) {
            result = false;
            r_message("cross_info has >1 columns, but should have just 1");
        }
        else {
            int n_missing = 0;
            for(int i=0; i<n_row; i++)
                if(cross_info[i] == NA_INTEGER) ++n_missing;
            if(n_missing > 0) {
                result = false;
                r_message("cross_info contains missing values (it shouldn't)");
            }

            int n_invalid = 0;
            for(int i=0; i<n_row; i++)
                if(cross_info[i] != NA_INTEGER && cross_info[i] != 0 && cross_info[i] != 1) ++n_invalid;
            if(n_invalid > 0) {
                result = false;
                r_message("cross_info contains invalid values; should be 0 or 1.");
            }
        }
    }
    return result;
}
