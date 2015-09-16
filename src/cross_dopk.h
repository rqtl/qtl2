// phase-known Diversity Outcross QTLCross class (for HMM, in particular est.map)

#ifndef CROSS_DOPK_H
#define CROSS_DOPK_H

class DOPK : public QTLCross
{
 public:
    DOPK(){
        crosstype = "dopk";
        phase_known_crosstype = "dopk";
    };
    ~DOPK(){};

    const bool check_geno(const int gen, const bool is_observed_value,
                          const bool is_x_chr, const bool is_female, const IntegerVector& cross_info);

    const double init(const int true_gen,
                      const bool is_x_chr, const bool is_female, const IntegerVector& cross_info);
    const double emit(const int obs_gen, const int true_gen, const double error_prob,
                      const IntegerVector& founder_geno, const bool is_x_chr,
                      const bool is_female, const IntegerVector& cross_info);
    const double step(const int gen_left, const int gen_right, const double rec_frac,
                      const bool is_x_chr, const bool is_female, const IntegerVector& cross_info);

    const IntegerVector possible_gen(const bool is_x_chr, const bool is_female, const IntegerVector& cross_info);

    const int ngen(const bool is_x_chr);

    const double nrec(const int gen_left, const int gen_right,
                      const bool is_x_chr, const bool is_female, const IntegerVector& cross_info);

    // this isn't actually implemented yet; just throws an error
    const double est_rec_frac(const NumericVector& gamma, const bool is_x_chr,
                              const IntegerMatrix& cross_info, const int n_gen);

    // this is to indicate that the f2pk cross type shouldn't exist on the R side
    // (it's strictly a device for using phase-known version for est_map)
    const bool crosstype_supported() {
        return false;
    }

    const NumericMatrix geno2allele_matrix(const bool is_x_chr);

    const bool check_is_female_vector(const LogicalVector& is_female, const bool any_x_chr);

    const bool check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr);

};

#endif // CROSS_DOPK_H
