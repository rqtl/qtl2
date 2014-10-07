#ifndef CROSS_HAPLOID_H
#define CROSS_HAPLOID_H

class HAPLOID : public QTLCross
{
 public:
    HAPLOID(){
        type = "haploid";
        phase_known_type = "haploid";
    };
    ~HAPLOID(){};

    bool check_geno(int gen, bool is_observed_value,
                    bool ignored1, bool ignored2, IntegerVector ignored3);

    double init(int true_gen, bool is_X_chr, bool is_female,
                IntegerVector cross_info);
    double emit(int obs_gen, int true_gen, double error_prob,
                bool is_X_chr, bool is_female, IntegerVector cross_info);
    double step(int gen_left, int gen_right, double rec_frac,
                bool is_X_chr, bool is_female, IntegerVector cross_info);

    int ngen(bool ignored1);

    double nrec(int gen_left, int gen_right,
                bool ignored1, bool ignored2, IntegerVector ignored3);

    double est_rec_frac(NumericMatrix full_gamma, bool is_X_chr);
};

#endif // CROSS_HAPLOID_H
