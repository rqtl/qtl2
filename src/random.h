// random integer from {low, low+1, ..., high}
int random_int(const int low, const int high);

// vector of random integers from {low, low+1, ..., high}
IntegerVector random_int(const int n, const int low, const int high);

// permute a vector of numbers
NumericVector permute_nvector(const NumericVector x);
IntegerVector permute_ivector(const IntegerVector x);

// permute a vector of numbers in place
void permute_nvector_inplace(NumericVector x);
void permute_ivector_inplace(IntegerVector x);

// get permutation of {0..(n-1)}
IntegerVector get_permutation(const int n);

// get a set of permutations of a vector, as columns of a matrix
NumericMatrix permute_nvector(const int n, const NumericVector x);
IntegerMatrix permute_ivector(const int n, const IntegerVector x);
