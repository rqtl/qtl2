// RI by selfing QTLCross class (for HMM)

#ifndef CROSS_RISELF_H
#define CROSS_RISELF_H

class RISELF : public QTLCross
{
 public:
    RISELF(){
        crosstype = "riself";
        phase_known_crosstype = "riself";
    };
    ~RISELF(){};

    const double step(const int gen_left, const int gen_right, const double rec_frac,
                      const bool is_x_chr, const bool is_female, const IntegerVector& cross_info);

    const double est_rec_frac(const NumericVector& gamma, const bool is_x_chr,
                              const IntegerMatrix& cross_info, const int n_gen);
};

#endif // CROSS_RISELF_H
