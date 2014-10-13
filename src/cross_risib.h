#ifndef CROSS_RISIB_H
#define CROSS_RISIB_H

class RISIB : public QTLCross
{
 public:
    RISIB(){
        type = "risib";
        phase_known_type = "risib";
    };

    ~RISIB(){};

    const double init(const int& true_gen,
                      const bool& is_x_chr, const bool& is_female, const IntegerVector& cross_info);

    const double step(const int& gen_left, const int& gen_right, const double& rec_frac,
                      const bool& is_x_chr, const bool& is_female, const IntegerVector& cross_info);

    const double est_rec_frac(const NumericMatrix& gamma, const bool& is_x_chr);
};

#endif // CROSS_RISIB_H
