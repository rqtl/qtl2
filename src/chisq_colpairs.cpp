// perform chi-square tests on all pairs of columns of a matrix

#include "chisq_colpairs.h"
#include <Rcpp.h>
using namespace Rcpp;


// perform chi-square test on all pairs of columns
// columns assumed to have values {1,2,3,...,k} for some k
// 0's and NA's are ignored
//
// [[Rcpp::export(".chisq_colpairs")]]
NumericMatrix chisq_colpairs(const IntegerMatrix& input) // matrix of integers; should be contiguous
{
    const int n_row = input.rows();
    const int n_col = input.cols();
    if(n_col < 2)
        throw std::invalid_argument("Need at least two columns.");

    NumericMatrix result(n_col,n_col);

    // find max value in each column
    IntegerVector max_value(n_col);
    for(int j=0; j<n_col; j++) {
        result(j,j) = 0.0;
        max_value[j] = 0;
        for(int i=0; i<n_row; i++) {
            if(!IntegerVector::is_na(input(i,j)) && input(i,j) > max_value[j])
                max_value[j] = input(i,j);
        }
    }

    // now do the chi-square tests
    for(int col1=0; col1<n_col-1; col1++) {
        for(int col2=col1+1; col2<n_col; col2++) {

            Rcpp::checkUserInterrupt();  // check for ^C from user

            IntegerVector sums1(max_value[col1]);
            IntegerVector sums2(max_value[col2]);
            for(int k=0; k<max_value[col1]; k++) sums1[k] = 0;
            for(int k=0; k<max_value[col2]; k++) sums2[k] = 0;

            IntegerMatrix counts(max_value[col1], max_value[col2]);
            for(int k1=0; k1<max_value[col1]; k1++)
                for(int k2=0; k2<max_value[col2]; k2++)
                    counts(k1,k2) = 0;
            int total=0;

            for(int k=0; k<n_row; k++) {
                if(!IntegerVector::is_na(input(k,col1)) &&
                   !IntegerVector::is_na(input(k,col2)) &&
                   input(k,col1)>0 && input(k,col2)>0) {
                    sums1[input(k,col1)-1]++;
                    sums2[input(k,col2)-1]++;
                    counts(input(k,col1)-1, input(k,col2)-1)++;
                    total++;
                }
            }

            if(total==0) {
                result(col1,col2) = result(col2,col1) = NA_REAL;
                continue;
            }

            result(col1,col2) = 0.0;
            for(int k1=0; k1<max_value[col1]; k1++) {
                for(int k2=0; k2<max_value[col2]; k2++) {
                    double expected = (double)(sums1[k1]*sums2[k2])/(double)total;
                    if(expected > 0) {
                        double numerator = (expected - (double)counts(k1,k2));
                        result(col1,col2) += numerator*numerator/expected;
                    }
                }
            }
            result(col2,col1) = result(col1,col2);
        }
    }

    return result;
}
