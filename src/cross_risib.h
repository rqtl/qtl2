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

    double init(int true_gen,
                bool is_x_chr, bool is_female, IntegerVector cross_info);

    double step(int gen_left, int gen_right, double rec_frac,
                bool is_x_chr, bool is_female, IntegerVector cross_info);

    double est_rec_frac(NumericMatrix full_gamma, bool is_x_chr);
};

#endif // CROSS_RISIB_H
