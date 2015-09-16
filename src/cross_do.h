// Diversity Outcross QTLCross class (for HMM)

#ifndef CROSS_DO_H
#define CROSS_DO_H

class DO : public QTLCross
{
 public:
    DO(){
        crosstype = "do";
        phase_known_crosstype = "dopk";
    };
    ~DO(){};

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

    const NumericMatrix geno2allele_matrix(const bool is_x_chr);

    const bool check_is_female_vector(const LogicalVector& is_female, const bool any_x_chr);

    const bool check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr);

    const bool is_het(const int true_gen); // is heterozygous (just for female X or autosome)

    const int encode_alleles(const int allele1, const int allele2); // convert (a1,a2) pair to genotype 1-36
    const IntegerVector decode_geno(const int true_gen);            // convert genotype to (a1,a2) pair

};

#endif // CROSS_DO_H
