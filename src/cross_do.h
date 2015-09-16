// Diversity Outcross QTLCross class (for HMM)

#ifndef CROSS_DO_H
#define CROSS_DO_H

class DO : public QTLCross
{
 public:
    DO(){
        crosstype = "do";
        phase_known_crosstype = "dopk";
    };
    ~DO(){};

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

    const NumericMatrix geno2allele_matrix(const bool is_x_chr);

    const bool check_is_female_vector(const LogicalVector& is_female, const bool any_x_chr);

    const bool check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr);

};

#endif // CROSS_DO_H
