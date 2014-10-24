// intercross QTLCross class (for HMM)

#ifndef CROSS_F2_H
#define CROSS_F2_H

class F2 : public QTLCross
{
 public:
    F2(){
        crosstype = "f2";
        phase_known_crosstype = "f2pk";
    };
    ~F2(){};

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

};

#endif // CROSS_F2_H
