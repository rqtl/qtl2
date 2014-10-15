// re-estimate inter-marker recombination fractions
NumericVector est_map(const String& crosstype,
                      const IntegerMatrix& genotypes, // columns are individuals, rows are markers
                      const IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                      const bool is_X_chr,
                      const LogicalVector& is_female,
                      const IntegerMatrix& cross_info,
                      const NumericVector& rec_frac,
                      const double error_prob,
                      const int max_iterations,
                      const double tol,
                      const bool verbose);

