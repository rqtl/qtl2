// random integer from {low, low+1, ..., high}
int random_int(const int low, const int high);

// vector of random integers from {low, low+1, ..., high}
IntegerVector random_int(const int n, const int low, const int high);

// permute a vector of numbers
NumericVector permute_nvector(const NumericVector x);
IntegerVector permute_ivector(const IntegerVector x);
vector<double> permute_nvector(const vector<double> x);
vector<int> permute_ivector(const vector<int> x);

// permute a vector of numbers in place
void permute_nvector_inplace(NumericVector x);
void permute_ivector_inplace(IntegerVector x);
void permute_nvector_inplace(vector<double> x);
void permute_ivector_inplace(vector<int> x);

// get permutation of {0..(n-1)}
IntegerVector get_permutation(const int n);

// get a set of permutations of a vector, as columns of a matrix
NumericMatrix permute_nvector(const int n_perm, const NumericVector x);
IntegerMatrix permute_ivector(const int n_perm, const IntegerVector x);

// stratified permutation
NumericMatrix permute_nvector_stratified(const int n_perm, const NumericVector& x,
                                         const IntegerVector& strata, int n_strata);
IntegerMatrix permute_ivector_stratified(const int n_perm, const IntegerVector& x,
                                         const IntegerVector& strata, int n_strata);
