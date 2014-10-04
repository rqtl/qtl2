// forward equations
NumericMatrix forwardEquations(Cross* cross,
                               IntegerVector genotypes,
                               bool is_X_chr,
                               bool is_female,
                               IntegerVector cross_info,
                               NumericVector rec_frac,
                               IntegerVector marker_index,
                               double error_prob,
                               IntegerVector poss_gen);


// backward Equations
NumericMatrix backwardEquations(Cross* cross,
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


// Calculate addlog(a,b) = log[exp(a) + exp(b)]
double addlog(const double a, const double b);
