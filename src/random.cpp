#include <Rcpp.h>
using namespace Rcpp;

#include "random.h"

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
    
    return(result);
}

// permute a vector of numbers
NumericVector permute_nvector(const NumericVector x)
{
    int n = x.size();

    NumericVector result(clone(x));

    for(int i=n-1; i>0; i--)
        std::swap(result[i], result[random_int(0, i)]);

    return result;
}

// permute a vector of integers
IntegerVector permute_ivector(const IntegerVector x)
{
    int n = x.size();

    IntegerVector result(clone(x));

    for(int i=n-1; i>0; i--)
        std::swap(result[i], result[random_int(0, i)]);

    return result;
}

// permute a vector of numbers in place
void permute_nvector_inplace(NumericVector x)
{
    int n = x.size();

    for(int i=n-1; i>0; i--)
        std::swap(x[i], x[random_int(0, i)]);
}

// permute a vector of integers in place
void permute_ivector_inplace(IntegerVector x)
{
    int n = x.size();

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
NumericMatrix permute_nvector(const int n, const NumericVector x)
{
    int length = x.size();

    NumericMatrix result(length,n);

    for(int i=0; i<n; i++) {
        NumericVector permx = permute_nvector(x);
        std::copy(permx.begin(), permx.end(), result.begin()+i*length);
    }

    return(result);
}

// get a set of permutations of a vector, as columns of a matrix
// [[Rcpp::export]]
IntegerMatrix permute_ivector(const int n, const IntegerVector x)
{
    int length = x.size();

    IntegerMatrix result(length,n);

    for(int i=0; i<n; i++) {
        IntegerVector permx = permute_ivector(x);
        std::copy(permx.begin(), permx.end(), result.begin()+i*length);
    }

    return(result);
}

