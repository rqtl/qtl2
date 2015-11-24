// linear regression via RcppEigen

// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

#include "linreg_eigen.h"
#include "debug_util.h"

// calc X'X
MatrixXd calc_XpX_eigen(const MatrixXd& X)
{
    int n = X.cols();

    return MatrixXd(n,n).setZero().selfadjointView<Lower>()
        .rankUpdate(X.adjoint());
}

// least squares by "LLt" Cholesky decomposition
// [[Rcpp::export]]
List fit_linreg_eigenchol(const NumericMatrix& X, const NumericVector& y)
{
    MatrixXd XX(as<Map<MatrixXd> >(X));
    VectorXd yy(as<Map<VectorXd> >(y));

    int n = XX.rows(), p=XX.cols();
    LLT<MatrixXd> llt = calc_XpX_eigen(XX);

    VectorXd betahat = llt.solve(XX.adjoint() * yy);
    VectorXd fitted = XX * betahat;
    VectorXd resid = yy - fitted;
    int df = n-p;
    double s = resid.norm() / std::sqrt((double)df);
    VectorXd se = s * llt.matrixL().solve(MatrixXd::Identity(p,p)).colwise().norm();
    double rss = resid.squaredNorm();

    return List::create(Named("coef") = betahat,
                        Named("fitted") = fitted,
                        Named("resid") = resid,
                        Named("rss") = rss,
                        Named("sigma") = s,
                        Named("rank") = p,
                        Named("df") = df,
                        Named("SE") = se);
}

// least squares by "LLt" Cholesky decomposition
// return just the residual sum of squares
// [[Rcpp::export]]
double calc_rss_eigenchol(const NumericMatrix& X, const NumericVector& y)
{
    MatrixXd XX(as<Map<MatrixXd> >(X));
    VectorXd yy(as<Map<VectorXd> >(y));

    LLT<MatrixXd> llt = calc_XpX_eigen(XX);

    VectorXd betahat = llt.solve(XX.adjoint() * yy);
    VectorXd fitted = XX * betahat;
    VectorXd resid = yy - fitted;
    return resid.squaredNorm();
}

// least squares by QR decomposition with column pivoting
// [[Rcpp::export]]
List fit_linreg_eigenqr(const NumericMatrix& X, const NumericVector& y)
{
    MatrixXd XX(as<Map<MatrixXd> >(X));
    VectorXd yy(as<Map<VectorXd> >(y));

    typedef ColPivHouseholderQR<MatrixXd> CPivQR;
    typedef CPivQR::PermutationType Permutation;

    int n = XX.rows(), p = XX.cols();

    CPivQR PQR( XX );
    Permutation Pmat( PQR.colsPermutation() );
    int r = PQR.rank();

    VectorXd betahat(p), fitted(n), se(p);

    if(r == p) { // full rank
        betahat = PQR.solve(yy);
        fitted = XX * betahat;

        se = Pmat * PQR.matrixQR().topRows(p).
            triangularView<Upper>().
            solve(MatrixXd::Identity(p, p)).
            rowwise().norm();

    } else {
        MatrixXd Rinv( PQR.matrixQR().topLeftCorner(r,r).
                       triangularView<Upper>().
                       solve(MatrixXd::Identity(r,r)) );
        VectorXd effects( PQR.householderQ().adjoint() * yy );

        betahat.fill(::NA_REAL);
        betahat.head(r) = Rinv * effects.head(r);
        betahat = Pmat*betahat;

        se.fill(::NA_REAL);
        se.head(r) = Rinv.rowwise().norm();
        se = Pmat * se;

        effects.tail(n - r).setZero();
        fitted = PQR.householderQ() * effects;
    }

    VectorXd resid = yy - fitted;
    double rss = resid.squaredNorm();
    int df = n - r;
    double sigma = std::sqrt(rss/(double)df);

    return List::create(Named("coef") = betahat,
                        Named("fitted") = fitted,
                        Named("resid") = resid,
                        Named("rss") = rss,
                        Named("sigma") = sigma,
                        Named("rank") = r,
                        Named("df") = df,
                        Named("SE") = sigma*se);
}

// least squares by QR decomposition with column pivoting
// return just the residual sum of squares
// [[Rcpp::export]]
double calc_rss_eigenqr(const NumericMatrix& X, const NumericVector& y)
{
    MatrixXd XX(as<Map<MatrixXd> >(X));
    VectorXd yy(as<Map<VectorXd> >(y));

    typedef Eigen::ColPivHouseholderQR<MatrixXd> CPivQR;
    typedef CPivQR::PermutationType Permutation;

    int n = XX.rows(), p = XX.cols();

    CPivQR PQR = XX;
    Permutation Pmat = PQR.colsPermutation();
    int r = PQR.rank();

    VectorXd fitted(n);
    if(r == p) { // full rank
        VectorXd betahat = PQR.solve(yy);
        fitted = XX * betahat;
    } else {
        MatrixXd Rinv = PQR.matrixQR().topLeftCorner(r,r)
            .triangularView<Upper>().solve(MatrixXd::Identity(r,r));
        VectorXd effects = PQR.householderQ().adjoint() * yy;
        effects.tail(n - r).setZero();
        fitted = PQR.householderQ() * effects;
    }

    VectorXd resid = yy - fitted;
    return resid.squaredNorm();
}


// least squares by "LLt" Cholesky decomposition, with matrix Y
// return vector of RSS
// [[Rcpp::export]]
NumericVector calc_mvrss_eigenchol(const NumericMatrix& X, const NumericMatrix& Y)
{
    int ncolY = Y.cols();
    int ncolX = X.cols();

    MatrixXd XX(as<Map<MatrixXd> >(X));
    MatrixXd YY(as<Map<MatrixXd> >(Y));

    LLT<MatrixXd> llt = calc_XpX_eigen(XX);

    MatrixXd XXpY(XX.adjoint() * YY);

    MatrixXd betahat(ncolX,ncolY);
    for(int i=0; i<ncolY; i++)
        betahat.col(i) = llt.solve(XXpY.col(i));

    MatrixXd fitted = XX * betahat;
    MatrixXd resid = YY - fitted;

    NumericVector rss(wrap(resid.colwise().squaredNorm().transpose()));
    return rss;
}

// least squares by QR decomposition with column pivoting, with matrix Y
// return vector of RSS
// [[Rcpp::export]]
NumericVector calc_mvrss_eigenqr(const NumericMatrix& X, const NumericMatrix& Y)
{
    MatrixXd XX(as<Map<MatrixXd> >(X));
    MatrixXd YY(as<Map<MatrixXd> >(Y));

    typedef Eigen::ColPivHouseholderQR<MatrixXd> CPivQR;
    typedef CPivQR::PermutationType Permutation;

    int n = XX.rows(), p = XX.cols();
    int k = YY.cols();

    CPivQR PQR = XX;
    Permutation Pmat = PQR.colsPermutation();
    int r = PQR.rank();

    MatrixXd fitted(n,k);

    if(r == p) { // full rank
        MatrixXd betahat(p,k);

        for(int i=0; i<k; i++)
            betahat.col(i) = PQR.solve(YY.col(i));

        fitted = XX * betahat;

    } else {

        MatrixXd Rinv = PQR.matrixQR().topLeftCorner(r,r)
            .triangularView<Upper>().solve(MatrixXd::Identity(r,r));

        for(int i=0; i<k; i++) {
            VectorXd effects = PQR.householderQ().adjoint() * YY.col(i);
            effects.tail(n - r).setZero();

            fitted.col(i) = PQR.householderQ() * effects;
        }
    }

    MatrixXd resid = YY - fitted;

    NumericVector rss(wrap(resid.colwise().squaredNorm().transpose()));
    return rss;
}

// least squares by "LLt" Cholesky decomposition, with matrix Y
// return matrix of residuals
// [[Rcpp::export]]
NumericMatrix calc_resid_eigenchol(const NumericMatrix& X, const NumericMatrix& Y)
{
    int ncolY = Y.cols();
    int ncolX = X.cols();

    MatrixXd XX(as<Map<MatrixXd> >(X));
    MatrixXd YY(as<Map<MatrixXd> >(Y));

    LLT<MatrixXd> llt = calc_XpX_eigen(XX);

    MatrixXd XXpY(XX.adjoint() * YY);

    MatrixXd betahat(ncolX,ncolY);
    for(int i=0; i<ncolY; i++)
        betahat.col(i) = llt.solve(XXpY.col(i));

    MatrixXd fitted = XX * betahat;
    MatrixXd resid = YY - fitted;

    NumericMatrix result(wrap(resid));

    return result;
}

// least squares by QR decomposition with column pivoting, with matrix Y
// return matrix of residuals
// [[Rcpp::export]]
NumericMatrix calc_resid_eigenqr(const NumericMatrix& X, const NumericMatrix& Y)
{
    MatrixXd XX(as<Map<MatrixXd> >(X));
    MatrixXd YY(as<Map<MatrixXd> >(Y));

    typedef Eigen::ColPivHouseholderQR<MatrixXd> CPivQR;
    typedef CPivQR::PermutationType Permutation;

    int n = XX.rows(), p = XX.cols();
    int k = YY.cols();

    CPivQR PQR = XX;
    Permutation Pmat = PQR.colsPermutation();
    int r = PQR.rank();

    MatrixXd fitted(n,k);

    if(r == p) { // full rank
        MatrixXd betahat(p,k);

        for(int i=0; i<k; i++)
            betahat.col(i) = PQR.solve(YY.col(i));

        fitted = XX * betahat;

    } else {

        MatrixXd Rinv = PQR.matrixQR().topLeftCorner(r,r)
            .triangularView<Upper>().solve(MatrixXd::Identity(r,r));

        for(int i=0; i<k; i++) {
            VectorXd effects = PQR.householderQ().adjoint() * YY.col(i);
            effects.tail(n - r).setZero();

            fitted.col(i) = PQR.householderQ() * effects;
        }
    }

    MatrixXd resid = YY - fitted;

    NumericMatrix result(wrap(resid));

    return result;
}
