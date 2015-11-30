// debugging utilities
#ifndef DEBUG_UTIL_H
#define DEBUG_UTIL_H

// print a vector
void print_vector(const Rcpp::NumericVector& x);
void print_vector(const Rcpp::IntegerVector& x);
void print_vector(const Eigen::VectorXd& x);
void print_vector(const Eigen::VectorXi& x);
// print matrix dimension
void print_matdim(const Rcpp::NumericMatrix& x);
void print_matdim(const Rcpp::IntegerMatrix& x);
void print_matdim(const Eigen::MatrixXd& x);
void print_matdim(const Eigen::MatrixXi& x);

#endif // DEBUG_UTIL_H
