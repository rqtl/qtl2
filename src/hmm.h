// forward equations
NumericMatrix forwardEquations(QTLCross* cross,
                               IntegerVector genotypes,
                               bool is_X_chr,
                               bool is_female,
                               IntegerVector cross_info,
                               NumericVector rec_frac,
                               IntegerVector marker_index,
                               double error_prob,
                               IntegerVector poss_gen);


// backward Equations
NumericMatrix backwardEquations(QTLCross* cross,
                                IntegerVector genotypes,
                                bool is_X_chr,
                                bool is_female,
                                IntegerVector cross_info,
                                NumericVector rec_frac,
                                IntegerVector marker_index,
                                double error_prob,
                                IntegerVector poss_gen);

// calculate QTL genotype probabilities
NumericVector calc_genoprob(String crosstype,
                            IntegerMatrix genotypes, // columns are individuals, rows are markers
                            bool is_X_chr,
                            LogicalVector is_female, // length n_ind
                            IntegerMatrix cross_info, // columns are individuals
                            NumericVector rec_frac,   // length nrow(genotypes)-1
                            IntegerVector marker_index, // length nrow(genotypes)
                            double error_prob);


// re-estimate inter-marker recombination fractions
NumericVector est_map(String crosstype,
                      IntegerMatrix genotypes,
                      bool is_X_chr,
                      LogicalVector is_female,
                      IntegerMatrix cross_info,
                      NumericVector rec_frac,
                      double error_prob,
                      int max_iterations,
                      double tol,
                      bool verbose);

// Calculate addlog(a,b) = log[exp(a) + exp(b)]
double addlog(const double a, const double b);
