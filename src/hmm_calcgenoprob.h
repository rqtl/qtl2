// calculate QTL genotype probabilities
NumericVector calc_genoprob(const String& crosstype,
                            const IntegerMatrix& genotypes, // columns are individuals, rows are markers
                            const IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                            const bool is_X_chr,
                            const LogicalVector& is_female, // length n_ind
                            const IntegerMatrix& cross_info, // columns are individuals
                            const NumericVector& rec_frac,   // length nrow(genotypes)-1
                            const IntegerVector& marker_index, // length nrow(genotypes)
                            const double error_prob);
