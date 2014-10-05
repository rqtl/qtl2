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

    bool check_geno(int gen, bool is_observed_value,
                    bool is_X_chr, bool is_female, IntegerVector cross_info);

    double init(int true_gen, bool is_X_chr, bool is_female,
                IntegerVector cross_info);
    double emit(int obs_gen, int true_gen, double error_prob, 
                bool is_X_chr, bool is_female, IntegerVector cross_info);
    double step(int gen_left, int gen_right, double rec_frac,
                bool is_X_chr, bool is_female, IntegerVector cross_info);

    IntegerVector possible_gen(bool is_X_chr, bool is_female,
                             IntegerVector cross_info);
    int ngen(bool is_X_chr);

    double nrec(int gen_left, int gen_right,
                bool is_X_chr, bool is_female, IntegerVector cross_info);

    double est_rec_frac(NumericMatrix full_gamma, bool is_X_chr);
};

#endif // CROSS_BC_H
