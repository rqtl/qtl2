// forward equations
NumericMatrix forwardEquations(QTLCross* cross,
                               const IntegerVector& genotypes,
                               const bool& is_X_chr,
                               const bool& is_female,
                               const IntegerVector& cross_info,
                               const NumericVector& rec_frac,
                               const IntegerVector& marker_index,
                               const double& error_prob,
                               const IntegerVector& poss_gen);


// backward Equations
NumericMatrix backwardEquations(QTLCross* cross,
                                const IntegerVector& genotypes,
                                const bool& is_X_chr,
                                const bool& is_female,
                                const IntegerVector& cross_info,
                                const NumericVector& rec_frac,
                                const IntegerVector& marker_index,
                                const double& error_prob,
                                const IntegerVector& poss_gen);

// calculate QTL genotype probabilities
NumericVector calc_genoprob(const String& crosstype,
                            const IntegerMatrix& genotypes, // columns are individuals, rows are markers
                            const bool& is_X_chr,
                            const LogicalVector& is_female, // length n_ind
                            const IntegerMatrix& cross_info, // columns are individuals
                            const NumericVector& rec_frac,   // length nrow(genotypes)-1
                            const IntegerVector& marker_index, // length nrow(genotypes)
                            const double& error_prob);


// re-estimate inter-marker recombination fractions
NumericVector est_map(const String& crosstype,
                      const IntegerMatrix& genotypes,
                      const bool& is_X_chr,
                      const LogicalVector& is_female,
                      const IntegerMatrix& cross_info,
                      const NumericVector& rec_frac,
                      const double& error_prob,
                      const int& max_iterations,
                      const double& tol,
                      const bool& verbose);

// Calculate addlog(a,b) = log[exp(a) + exp(b)]
double addlog(const double a, const double b);
