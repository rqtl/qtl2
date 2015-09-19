// phase-known Diversity Outcross QTLCross class (for HMM, in particular est.map)

#ifndef CROSS_DOPK_H
#define CROSS_DOPK_H

class DOPK : public QTLCross
{
 public:
    DOPK(){
        crosstype = "dopk";
        phase_known_crosstype = "dopk";
    };
    ~DOPK(){};

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

    const double nrec(const int gen_left, const int gen_right,
                      const bool is_x_chr, const bool is_female, const IntegerVector& cross_info);

    // this isn't actually implemented yet; just throws an error
    const double est_rec_frac(const NumericVector& gamma, const bool is_x_chr,
                              const IntegerMatrix& cross_info, const int n_gen);

    // this is to indicate that the f2pk cross type shouldn't exist on the R side
    // (it's strictly a device for using phase-known version for est_map)
    const bool crosstype_supported() {
        return false;
    }

    const NumericMatrix geno2allele_matrix(const bool is_x_chr);

    const bool check_is_female_vector(const LogicalVector& is_female, const bool any_x_chr);

    const bool check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr);

    const bool is_het(const int true_gen); // is heterozygous (just for female X or autosome)

    const int encode_alleles(const int allele1, const int allele2); // convert (a1,a2) pair to genotype 1-36
    const IntegerVector decode_geno(const int true_gen);            // convert genotype to (a1,a2) pair

    // helper functions for step()
    const double step_auto(int left, int right, double r, int s,
                           IntegerVector precc_gen, NumericVector precc_alpha);
    const double step_femX(int left, int right, double r, int s,
                           IntegerVector precc_gen, NumericVector precc_alpha);
    const double step_malX(int left, int right, double r, int s,
                           IntegerVector precc_gen, NumericVector precc_alpha);

    const bool check_founder_geno_size(const IntegerMatrix& founder_geno, const int n_markers);
    const bool check_founder_geno_values(const IntegerMatrix& founder_geno);
    const bool need_founder_geno();

};

#endif // CROSS_DOPK_H
