// debugging utilities

#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

#include "debug_util.h"

// print a vector
void print_vector(const NumericVector& x)
{
    int n=x.size();
    for(int i=0; i<n; i++)
        Rcout << x[i] << ' ';
    Rcout << '\n';
}

void print_vector(const IntegerVector& x)
{
    int n=x.size();
    for(int i=0; i<n; i++)
        Rcout << x[i] << ' ';
    Rcout << '\n';
}

void print_vector(const VectorXd& x)
{
    int n=x.size();
    for(int i=0; i<n; i++)
        Rcout << x[i] << ' ';
    Rcout << '\n';
}

void print_vector(const VectorXi& x)
{
    int n=x.size();
    for(int i=0; i<n; i++)
        Rcout << x[i] << ' ';
    Rcout << '\n';
}

void print_matdim(const NumericMatrix& x)
{
    Rcout << x.rows() << " x " << x.cols() << '\n';
}

void print_matdim(const IntegerMatrix& x)
{
    Rcout << x.rows() << " x " << x.cols() << '\n';
}

void print_matdim(const MatrixXd& x)
{
    Rcout << x.rows() << " x " << x.cols() << '\n';
}

void print_matdim(const MatrixXi& x)
{
    Rcout << x.rows() << " x " << x.cols() << '\n';
}
