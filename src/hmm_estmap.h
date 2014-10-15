// re-estimate inter-marker recombination fractions
NumericVector est_map(const String& crosstype,
                      const IntegerMatrix& genotypes,
                      const bool is_X_chr,
                      const LogicalVector& is_female,
                      const IntegerMatrix& cross_info,
                      const NumericVector& rec_frac,
                      const double error_prob,
                      const int max_iterations,
                      const double tol,
                      const bool verbose);

