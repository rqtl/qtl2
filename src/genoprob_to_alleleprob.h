// convert genotype probabilities to allele probabilities
NumericVector genoprob_to_alleleprob(const String& crosstype,
                                     const NumericVector& prob_array, // array as n_gen x n_ind x n_pos
                                     const bool is_x_chr);
