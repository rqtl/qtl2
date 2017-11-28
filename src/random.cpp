// random number generation

#include "random.h"
<<<<<<< HEAD
#include <Rcpp.h>
using namespace Rcpp;
=======
#include <vector>
#include <map>
#include <Rcpp.h>
using namespace Rcpp;
using std::vector;
using std::map;


// random integer from {low, low+1, ..., high}
int random_int(const int low, const int high)
{
    return (int)R::runif((double)low, double(high+1));
}

// vector of random integers from {low, low+1, ..., high}
// [[Rcpp::export]]
IntegerVector random_int(const int n, const int low, const int high)
{
    IntegerVector result(n);

    for(int i=0; i<n; i++)
        result[i] = random_int(low, high);

    return result;
}

// permute a vector of numbers
NumericVector permute_nvector(const NumericVector x)
{
    const int n = x.size();

    NumericVector result(clone(x));

    for(int i=n-1; i>0; i--)
        std::swap(result[i], result[random_int(0, i)]);

    return result;
}

// permute a vector of integers
IntegerVector permute_ivector(const IntegerVector x)
{
    const int n = x.size();

    IntegerVector result(clone(x));

    for(int i=n-1; i>0; i--)
        std::swap(result[i], result[random_int(0, i)]);

    return result;
}

// permute a vector of numbers
vector<double> permute_nvector(const vector<double> x)
{
    const int n = x.size();

    vector<double> result(x);

    for(int i=n-1; i>0; i--)
        std::swap(result[i], result[random_int(0, i)]);

    return result;
}


// permute a vector of numbers
vector<int> permute_ivector(const vector<int> x)
{
    const int n = x.size();

    vector<int> result(x);

    for(int i=n-1; i>0; i--)
        std::swap(result[i], result[random_int(0, i)]);

    return result;
}


// permute a vector of numbers in place
void permute_nvector_inplace(NumericVector x)
{
    const int n = x.size();

    for(int i=n-1; i>0; i--)
        std::swap(x[i], x[random_int(0, i)]);
}

// permute a vector of integers in place
void permute_ivector_inplace(IntegerVector x)
{
    const int n = x.size();

    for(int i=n-1; i>0; i--)
        std::swap(x[i], x[random_int(0, i)]);
}

// permute a vector of numbers in place
void permute_nvector_inplace(vector<double> x)
{
    const int n = x.size();

    for(int i=n-1; i>0; i--)
        std::swap(x[i], x[random_int(0, i)]);
}

// permute a vector of integers in place
void permute_ivector_inplace(vector<int> x)
{
    const int n = x.size();

    for(int i=n-1; i>0; i--)
        std::swap(x[i], x[random_int(0, i)]);
}

// get permutation of {0..(n-1)}
// [[Rcpp::export]]
IntegerVector get_permutation(const int n)
{
    IntegerVector result(n);

    for(int i=0; i<n; i++) result[i] = i;

    permute_ivector_inplace(result);

    return result;
}

// get a set of permutations of a vector, as columns of a matrix
// [[Rcpp::export]]
NumericMatrix permute_nvector(const int n_perm, const NumericVector x)
{
    const int length = x.size();

    NumericMatrix result(length,n_perm);

    for(int i=0; i<n_perm; i++) {
        NumericVector permx = permute_nvector(x);
        std::copy(permx.begin(), permx.end(), result.begin()+i*length);
    }
>>>>>>> qtl2scan/master


// sample random integer from 0, 1, 2, ..., n-1 with probability p[0], p[1], ...
int sample_int(NumericVector probs)
{
<<<<<<< HEAD
    int n=probs.size();
    int result;
=======
    const int length = x.size();
>>>>>>> qtl2scan/master

    double u = R::runif(0.0, 1.0);

<<<<<<< HEAD
    for(result=0; result < n; result++) {
        if(u <= probs[result]) return result;
        u -= probs[result];
=======
    for(int i=0; i<n_perm; i++) {
        IntegerVector permx = permute_ivector(x);
        std::copy(permx.begin(), permx.end(), result.begin()+i*length);
>>>>>>> qtl2scan/master
    }

    return NA_INTEGER;
}

// sample random integer from 0, 1, 2, ..., n-1 with equal probabilities
int sample_int(int n)
{
<<<<<<< HEAD
    return (int)(unif_rand()*n);
=======
    const int n = x.size();
    NumericMatrix result(n,n_perm);

    if(strata.size() != n)
        throw std::length_error("length(x) != length(strata)");

    if(n_strata < 0) // find maximum strata
        n_strata = max(strata) + 1;

    // map of indices for the strata
    map<int, vector<int> > strata_index;
    for(int i=0; i<n; ++i) {
        if(strata[i] >= n_strata || strata[i] < 0)
            throw std::domain_error("strata should be in [0, n_strata)");
        strata_index[strata[i]].push_back(i);
    }

    for(int perm=0; perm<n_perm; ++perm) {
        // for each stratum:
        for(int stratum=0; stratum < n_strata; ++stratum) {
            // permute indices
            vector<int> index_permuted = permute_ivector(strata_index[stratum]);

            int n = strata_index[stratum].size();
            for(int i=0; i<n; ++i)
                result(strata_index[stratum][i],perm) = x[index_permuted[i]];
        }
    }

    return result;
}

// permute x within strata
//     strata is integer vector {0, 1, 2, ..., n_strata-1}
//     n_strata is the number of strata; if == -1, it is calculated
// [[Rcpp::export]]
IntegerMatrix permute_ivector_stratified(const int n_perm, const IntegerVector& x,
                                         const IntegerVector& strata, int n_strata = -1)
{
    const int n = x.size();
    IntegerMatrix result(n,n_perm);

    if(strata.size() != n)
        throw std::length_error("length(x) != length(strata)");

    if(n_strata < 0) // find maximum strata
        n_strata = max(strata) + 1;

    // map of indices for the strata
    map<int, vector<int> > strata_index;
    for(int i=0; i<n; ++i) {
        if(strata[i] >= n_strata || strata[i] < 0)
            throw std::domain_error("strata should be in [0, n_strata)");
        strata_index[strata[i]].push_back(i);
    }

    for(int perm=0; perm<n_perm; ++perm) {
        // for each stratum:
        for(int stratum=0; stratum < n_strata; ++stratum) {
            // permute indices
            vector<int> index_permuted = permute_ivector(strata_index[stratum]);

            int n = strata_index[stratum].size();
            for(int i=0; i<n; ++i)
                result(strata_index[stratum][i],perm) = x[index_permuted[i]];
        }
    }

    return result;
>>>>>>> qtl2scan/master
}
