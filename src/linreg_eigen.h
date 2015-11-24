// linear regression via RcppEigen
#ifndef LINREG_EIGEN_H
#define LINREG_EIGEN_H

// calc X'X
MatrixXd calc_XpX_eigen(const MatrixXd& A);

// least squares by "LLt" Cholesky decomposition
// needs to be full rank
List fit_linreg_eigenchol(const NumericMatrix& X, const NumericVector& y);

// least squares by "LLt" Cholesky decomposition
// return just the residual sum of squares
// needs to be full rank
double calc_rss_eigenchol(const NumericMatrix& X, const NumericVector& y);

// least squares by QR decomposition with column pivoting
List fit_linreg_eigenqr(const NumericMatrix& X, const NumericVector& y);

// least squares by QR decomposition with column pivoting
// return just the residual sum of squares
double calc_rss_eigenqr(const NumericMatrix& X, const NumericVector& y);

// least squares by "LLt" Cholesky decomposition, with matrix Y
// return vector of RSS
NumericVector calc_mvrss_eigenchol(const NumericMatrix& X, const NumericMatrix& Y);

// least squares by QR decomposition with column pivoting, with matrix Y
// return vector of RSS
NumericVector calc_mvrss_eigenqr(const NumericMatrix& X, const NumericMatrix& Y);

// least squares by "LLt" Cholesky decomposition, with matrix Y
// return matrix of residuals
NumericMatrix calc_resid_eigenchol(const NumericMatrix& X, const NumericMatrix& Y);

// least squares by QR decomposition with column pivoting, with matrix Y
// return matrix of residuals
NumericMatrix calc_resid_eigenqr(const NumericMatrix& X, const NumericMatrix& Y);

#endif // LINREG_EIGEN_H
