#ifndef CROSS_RISELF_H
#define CROSS_RISELF_H

class RISELF : public QTLCross
{
 public:
    RISELF(){
        type = "riself";
        phase_known_type = "riself";
    };
    ~RISELF(){};

    double step(int gen_left, int gen_right, double rec_frac,
                bool is_X_chr, bool is_female, IntegerVector cross_info);

    double est_rec_frac(NumericMatrix full_gamma, bool is_X_chr);
};

#endif // CROSS_RISELF_H
