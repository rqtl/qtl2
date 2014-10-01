// random integer from {low, low+1, ..., high}
int random_int(const int low, const int high);

// vector of random integers from {low, low+1, ..., high}
IntegerVector random_int(const int n, const int low, const int high);

// permute a vector of numbers
NumericVector permute_nvector(const NumericVector x);

// permute a vector of integers
IntegerVector permute_ivector(const IntegerVector x);

// permute a vector of numbers in place
void permute_nvector_inplace(NumericVector x);

// permute a vector of integers in place
void permute_ivector_inplace(IntegerVector x);

// get permutation of {0..(n-1)}
IntegerVector get_permutation(const int n);

// get a set of permutations of a vector, as columns of a matrix
NumericMatrix get_permutations(const int n, NumericVector x);
