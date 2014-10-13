#ifndef CROSS_F2PK_H
#define CROSS_F2PK_H

class F2PK : public QTLCross
{
 public:
    F2PK(){
        type = "f2pk";
        phase_known_type = "f2pk";
    };
    ~F2PK(){};

    const bool check_geno(const int& gen, const bool& is_observed_value,
                          const bool& is_x_chr, const bool& is_female, const IntegerVector& cross_info);

    const double init(const int& true_gen,
                      const bool& is_x_chr, const bool& is_female, const IntegerVector& cross_info);
    const double emit(const int& obs_gen, const int& true_gen, const double& error_prob,
                      const bool& is_x_chr, const bool& is_female, const IntegerVector& cross_info);
    const double step(const int& gen_left, const int& gen_right, const double& rec_frac,
                      const bool& is_x_chr, const bool& is_female, const IntegerVector& cross_info);

    const IntegerVector possible_gen(const bool& is_x_chr, const bool& is_female, const IntegerVector& cross_info);

    const int ngen(const bool& is_x_chr);

    const double nrec(const int& gen_left, const int& gen_right,
                      const bool& is_x_chr, const bool& is_female, const IntegerVector& cross_info);

    const double est_rec_frac(const NumericMatrix& gamma, const bool& is_x_chr);

    // this is to indicate that the f2pk cross type shouldn't exist on the R side
    // (it's strictly a device for using phase-known version for est_map)
    const bool crosstype_supported() {
        return false;
    }

};

#endif // CROSS_F2PK_H
