#ifndef CROSS_BC_H
#define CROSS_BC_H

class BC : public QTLCross
{
 public:
    BC(){
        type = "bc";
        phase_known_type = "bc";
    };
    ~BC(){};

    const bool check_geno(const int gen, const bool is_observed_value,
                          const bool is_x_chr, const bool is_female, const IntegerVector cross_info);

    const double init(const int true_gen, const bool is_x_chr, const bool is_female,
                      const IntegerVector cross_info);
    const double emit(const int obs_gen, const int true_gen, const double error_prob, 
                      const bool is_x_chr, const bool is_female, const IntegerVector cross_info);
    const double step(const int gen_left, const int gen_right, const double rec_frac,
                      const bool is_x_chr, const bool is_female, const IntegerVector cross_info);

    const IntegerVector possible_gen(const bool is_x_chr, const bool is_female,
                                     const IntegerVector cross_info);
    const int ngen(const bool is_x_chr);

    const double nrec(const int gen_left, const int gen_right,
                      const bool is_x_chr, const bool is_female, const IntegerVector cross_info);

    const double est_rec_frac(NumericMatrix full_gamma, const bool is_x_chr);
};

#endif // CROSS_BC_H
