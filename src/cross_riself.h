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

    const double step(const int gen_left, const int gen_right, const double rec_frac,
                      const bool& is_x_chr, const bool& is_female, const IntegerVector& cross_info);

    const double est_rec_frac(const NumericMatrix& gamma, const bool& is_x_chr);
};

#endif // CROSS_RISELF_H
