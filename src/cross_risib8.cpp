// 8-way RIL by sib-mating QTLCross class (for HMM)

#include "cross_risib8.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_util.h"
#include "cross_do_util.h"
#include "r_message.h" // defines RQTL2_NODEBUG and r_message()

enum gen {A=1, H=2, B=3, notA=5, notB=4};

const bool RISIB8::check_geno(const int gen, const bool is_observed_value,
                                const bool is_x_chr, const bool is_female,
                                const IntegerVector& cross_info)
{
    // allow any value 0-5 for observed
    if(is_observed_value) {
        if(gen==0 || gen==A || gen==H || gen==B ||
           gen==notA || gen==notB) return true;
        else return false;
    }

    if(!is_x_chr) { // autosome
        const int n_geno = 8;
        if(gen>= 1 && gen <= n_geno) return true;
    }
    else { // X chromosome (can be A, B, C, E, F but not D, G, H)
        if(gen >= 1 && gen <= 8 &&
           gen != cross_info[3] &&
           gen != cross_info[6] &&
           gen != cross_info[7])
            return true;
    }

    return false; // otherwise a problem
}

const double RISIB8::init(const int true_gen,
                            const bool is_x_chr, const bool is_female,
                            const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(!is_x_chr) { // autosome
        return -log(8.0);
    }
    else { // X chromosome Pr(A)=Pr(B)=Pr(C)=1/3
        if(true_gen == cross_info[2]) return -log(3.0);
        else return -log(6.0);
    }
}

const double RISIB8::emit(const int obs_gen, const int true_gen, const double error_prob,
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


const double RISIB8::step(const int gen_left, const int gen_right, const double rec_frac,
                            const bool is_x_chr, const bool is_female,
                            const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(!is_x_chr) {
        // equations are from Broman (2005) Genetics 169:1133-1146
        //    doi:10.1534/genetics.104.035212
        //    see bottom equation in right column on page 1137
        //    (need to multiply by 8 to get conditional probabilities,
        //     and note that there was an error in the i != j case)
        if(gen_left == gen_right)
            return log(1.0 - rec_frac) - log(1.0 + 6.0 * rec_frac);

        return log(rec_frac) - log(1.0 + 6.0 * rec_frac);
    }
    else { // X chr; need to use founder order
        // equations are from Broman (2005) Genetics 169:1133-1146
        //    doi:10.1534/genetics.104.035212
        //    see table 4 page 1137
        //    (need to multiply by the marginal probability, 1/6 or 1/3,
        //     to get these conditional probabilities)
        if(gen_left == gen_right) {
            if(gen_left == cross_info[2])
                return - log(1.0 + 4.0 * rec_frac);
            else
                return log(1.0 - rec_frac) - log(1.0 + 4.0 * rec_frac);
        }

        if(gen_right == cross_info[2])
            return log(2.0) + log(rec_frac) - log(1.0 + 4.0 * rec_frac);

        // otherwise...
        return log(rec_frac) - log(1.0 + 4.0 * rec_frac);
    }
}

const IntegerVector RISIB8::possible_gen(const bool is_x_chr, const bool is_female,
                                       const IntegerVector& cross_info)
{
    if(!is_x_chr) { // autosome
        int n_geno = 8;
        IntegerVector result(n_geno);

        for(int i=0; i<n_geno; i++) result[i] = i+1;
        return result;
    }
    else { // X chromosome
        int n_geno = 5;
        IntegerVector result(n_geno);

        result[0] = cross_info[0];
        result[1] = cross_info[1];
        result[2] = cross_info[2];
        result[3] = cross_info[4];
        result[4] = cross_info[5];

        return result;
    }
}

const int RISIB8::ngen(const bool is_x_chr)
{
    return 8;
}

const int RISIB8::nalleles()
{
    return 8;
}


// check that cross_info conforms to expectation
const bool RISIB8::check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr)
{
    bool result = true;
    const int n_row = cross_info.rows();
    const int n_col = cross_info.cols();
    // 16 columns with order of cross

    if(n_col != 8) {
        result = false;
        r_message("cross_info should have 8 columns, indicating the order of the cross");
        return result;
    }

    int n_missing=0;
    int n_invalid=0;
    for(int i=0; i<n_row; i++) {
        for(int j=0; j<n_col; j++) {
            if(cross_info(i,j) == NA_INTEGER) ++n_missing;
            else if(cross_info(i,j) < 1 || cross_info(i,j)>n_col) ++n_invalid;
        }
        // count values 1..ncol
        IntegerVector counts(n_col);
        for(int j=0; j<n_col; j++) counts[j] = 0; // zero counts
        for(int j=0; j<n_col; j++) ++counts[cross_info(i,j)-1]; // count values
        for(int j=0; j<n_col; j++) {
            if(counts[j] != 1) n_invalid += abs(counts[j] - 1);
        }
    }
    if(n_missing > 0) {
        result = false;
        r_message("cross_info has missing values (it shouldn't)");
    }
    if(n_invalid > 0) {
        result = false;
        r_message("cross_info has invalid values; each row should be permutation of {1, 2, ..., 8}");
    }

    return result;
}


// check that founder genotype data has correct no. founders and markers
const bool RISIB8::check_founder_geno_size(const IntegerMatrix& founder_geno, const int n_markers)
{
    bool result=true;

    const int fg_mar = founder_geno.cols();
    const int fg_f   = founder_geno.rows();

    if(fg_mar != n_markers) {
        result = false;
        r_message("founder_geno has incorrect number of markers");
    }

    if(fg_f != 8) {
        result = false;
        r_message("founder_geno should have 4 founders");
    }

    return result;
}

// check that founder genotype data has correct values
const bool RISIB8::check_founder_geno_values(const IntegerMatrix& founder_geno)
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

const bool RISIB8::need_founder_geno()
{
    return true;
}

// geno_names from allele names
const std::vector<std::string> RISIB8::geno_names(const std::vector<std::string> alleles,
                                                const bool is_x_chr)
{
    if(alleles.size() < 8)
        throw std::range_error("alleles must have length 8");

    const int n_alleles = 8;

    std::vector<std::string> result(n_alleles);

    for(int i=0; i<n_alleles; i++)
        result[i] = alleles[i] + alleles[i];

    return result;
}


const double RISIB8::est_rec_frac(const Rcpp::NumericVector& gamma, const bool is_x_chr,
                                    const Rcpp::IntegerMatrix& cross_info, const int n_gen)
{


    if(!is_x_chr) { // autosome: solve R=7r/(1+6r) for r
        double R = QTLCross::est_rec_frac(gamma, is_x_chr, cross_info, n_gen);
        return R / (7.0 - 6.0 * R);
    }
    else {
        int n_ind = cross_info.cols();
        int n_gen_sq = n_gen*n_gen;

        // three groups of counts
        double a=0.0, b=0.0, c=0.0;

        for(int ind=0, offset=0; ind<n_ind; ind++, offset += n_gen_sq) {
            int founder_c = cross_info(2, ind) - 1; // third founder in cross = "C"

            for(int gl=0; gl<n_gen; gl++) {
                if(gl == founder_c)
                    c += gamma[offset+gl*n_gen+gl];
                else
                    a += gamma[offset+gl*n_gen+gl];

                for(int gr=gl+1; gr<n_gen; gr++) {
                    b += (gamma[offset+gl*n_gen+gr] + gamma[offset+gr*n_gen+gl]);
                }
            }
        }

        // solved via maxima
        return (4*c + b+ 5*a - sqrt(16*c*c + 8*(5*a-b)*c + b*b + 10*a*b + 25*a*a))/8/c;
    }
}

// check whether X chr can be handled
const bool RISIB8::check_handle_x_chr(const bool any_x_chr)
{
    return true; // most crosses can handle the X chr
}

const NumericVector RISIB8::est_map2(const IntegerMatrix& genotypes,
                                     const IntegerMatrix& founder_geno,
                                     const bool is_X_chr,
                                     const LogicalVector& is_female,
                                     const IntegerMatrix& cross_info,
                                     const IntegerVector& cross_group,
                                     const IntegerVector& unique_cross_group,
                                     const NumericVector& rec_frac,
                                     const double error_prob,
                                     const int max_iterations,
                                     const double tol,
                                     const bool verbose)
{
    if(!is_X_chr) { // autosome; can ignore founder order
        const int n_ind = cross_group.size();
        Rcpp::IntegerVector one_group(n_ind);
        for(int i=0; i<n_ind; i++) one_group[i] = 0;
        Rcpp::IntegerVector one_unique_group(1);
        one_unique_group[0] = 0;

        return est_map2_grouped(this->crosstype,
                                genotypes, founder_geno,
                                is_X_chr, is_female, cross_info,
                                one_group, one_unique_group,
                                rec_frac, error_prob, max_iterations,
                                tol, verbose);
    }

    // X chromosome: need to use the lowmem version for now
    return est_map2_lowmem(this->crosstype,
                           genotypes, founder_geno,
                           is_X_chr, is_female, cross_info,
                           cross_group, unique_cross_group,
                           rec_frac, error_prob, max_iterations,
                           tol, verbose);
}
