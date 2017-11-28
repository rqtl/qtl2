// linear regression via RcppEigen

// [[Rcpp::depends(RcppEigen)]]

#include "linreg_eigen.h"
#include <RcppEigen.h>
#include "r_message.h" // defines RQTL2_NODEBUG

using namespace Rcpp;
using namespace Eigen;

// calc X'X
MatrixXd calc_XpX(const MatrixXd& X)
{
    const int n = X.cols();

    return MatrixXd(n,n).setZero().selfadjointView<Lower>()
        .rankUpdate(X.adjoint());
}

// least squares by "LLt" Cholesky decomposition
// [[Rcpp::export]]
List fit_linreg_eigenchol(const NumericMatrix& X, const NumericVector& y, const bool se)
{
    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const VectorXd yy(as<Map<VectorXd> >(y));

    const int n = XX.rows(), p=XX.cols();
    #ifndef RQTL2_NODEBUG
    if(n != y.size())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    LLT<MatrixXd> llt ( calc_XpX(XX) );

    VectorXd betahat = llt.solve(XX.adjoint() * yy);
    VectorXd fitted = XX * betahat;
    VectorXd resid = yy - fitted;
    const int df = n-p;
    const double s = resid.norm() / std::sqrt((double)df);
    const double rss = resid.squaredNorm();
    if(se) {
        VectorXd SE = s * llt.matrixL().solve(MatrixXd::Identity(p,p)).colwise().norm();

        return List::create(Named("coef") = betahat,
                            Named("fitted") = fitted,
                            Named("resid") = resid,
                            Named("rss") = rss,
                            Named("sigma") = s,
                            Named("rank") = p,
                            Named("df") = df,
                            Named("SE") = SE);
    } else {
        return List::create(Named("coef") = betahat,
                            Named("fitted") = fitted,
                            Named("resid") = resid,
                            Named("rss") = rss,
                            Named("sigma") = s,
                            Named("rank") = p,
                            Named("df") = df);
    }
}

// least squares by "LLt" Cholesky decomposition
// return just the coefficients
// [[Rcpp::export]]
NumericVector calc_coef_linreg_eigenchol(const NumericMatrix& X, const NumericVector& y)
{
    #ifndef RQTL2_NODEBUG
    if(X.rows() != y.size())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const VectorXd yy(as<Map<VectorXd> >(y));

    LLT<MatrixXd> llt ( calc_XpX(XX) );

    NumericVector betahat(wrap(llt.solve(XX.adjoint() * yy)));
    return betahat;
}


// least squares by "LLt" Cholesky decomposition
// return the coefficients and SEs
// [[Rcpp::export]]
List calc_coefSE_linreg_eigenchol(const NumericMatrix& X, const NumericVector& y)
{
    #ifndef RQTL2_NODEBUG
    if(X.rows() != y.size())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const VectorXd yy(as<Map<VectorXd> >(y));

    const int n = XX.rows(), p=XX.cols();
    LLT<MatrixXd> llt ( calc_XpX(XX) );

    VectorXd betahat = llt.solve(XX.adjoint() * yy);
    VectorXd fitted = XX * betahat;
    VectorXd resid = yy - fitted;
    const int df = n-p;
    const double s = resid.norm() / std::sqrt((double)df);
    VectorXd se = s * llt.matrixL().solve(MatrixXd::Identity(p,p)).colwise().norm();

    return List::create(Named("coef") = betahat,
                        Named("SE") = se);
}

// least squares by "LLt" Cholesky decomposition
// return just the residual sum of squares
// [[Rcpp::export]]
double calc_rss_eigenchol(const NumericMatrix& X, const NumericVector& y)
{
    #ifndef RQTL2_NODEBUG
    if(X.rows() != y.size())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const VectorXd yy(as<Map<VectorXd> >(y));

    LLT<MatrixXd> llt ( calc_XpX(XX) );

    VectorXd betahat = llt.solve(XX.adjoint() * yy);
    VectorXd fitted = XX * betahat;
    VectorXd resid = yy - fitted;
    return resid.squaredNorm();
}

// least squares by "LLt" Cholesky decomposition
// return just the fitted values
// [[Rcpp::export]]
NumericVector calc_fitted_linreg_eigenchol(const NumericMatrix& X, const NumericVector& y)
{
    #ifndef RQTL2_NODEBUG
    if(X.rows() != y.size())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const VectorXd yy(as<Map<VectorXd> >(y));

    LLT<MatrixXd> llt ( calc_XpX(XX) );

    VectorXd betahat = llt.solve(XX.adjoint() * yy);
    VectorXd fitted = XX * betahat;
    return wrap(fitted);
}

// least squares by QR decomposition with column pivoting
// [[Rcpp::export]]
List fit_linreg_eigenqr(const NumericMatrix& X, const NumericVector& y,
                        const bool se, const double tol=1e-12)
{
    #ifndef RQTL2_NODEBUG
    if(X.rows() != y.size())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const VectorXd yy(as<Map<VectorXd> >(y));

    typedef ColPivHouseholderQR<MatrixXd> CPivQR;
    typedef CPivQR::PermutationType Permutation;

    const int n = XX.rows(), p = XX.cols();

    CPivQR PQR( XX );
    PQR.setThreshold(tol); // set tolerance
    Permutation Pmat( PQR.colsPermutation() );
    const int r = PQR.rank();

    VectorXd betahat(p), fitted(n), SE(p);

    if(r == p) { // full rank
        betahat = PQR.solve(yy);
        fitted = XX * betahat;

        if(se) {
            SE = Pmat * PQR.matrixQR().topRows(p).
                triangularView<Upper>().
                solve(MatrixXd::Identity(p, p)).
                rowwise().norm();
        }

    } else {
        MatrixXd Rinv( PQR.matrixQR().topLeftCorner(r,r).
                       triangularView<Upper>().
                       solve(MatrixXd::Identity(r,r)) );
        VectorXd effects( PQR.householderQ().adjoint() * yy );

        betahat.fill(::NA_REAL);
        betahat.head(r) = Rinv * effects.head(r);
        betahat = Pmat*betahat;

        if(se) {
            SE.fill(::NA_REAL);
            SE.head(r) = Rinv.rowwise().norm();
            SE = Pmat * SE;
        }

        effects.tail(n - r).setZero();
        fitted = PQR.householderQ() * effects;
    }

    VectorXd resid = yy - fitted;
    const double rss = resid.squaredNorm();
    const int df = n - r;
    const double sigma = std::sqrt(rss/(double)df);

    if(se)
        return List::create(Named("coef") = betahat,
                            Named("fitted") = fitted,
                            Named("resid") = resid,
                            Named("rss") = rss,
                            Named("sigma") = sigma,
                            Named("rank") = r,
                            Named("df") = df,
                            Named("SE") = sigma*SE);
    else
        return List::create(Named("coef") = betahat,
                            Named("fitted") = fitted,
                            Named("resid") = resid,
                            Named("rss") = rss,
                            Named("sigma") = sigma,
                            Named("rank") = r,
                            Named("df") = df);
}

// least squares by QR decomposition with column pivoting
// this just returns the coefficients
// [[Rcpp::export]]
NumericVector calc_coef_linreg_eigenqr(const NumericMatrix& X, const NumericVector& y,
                                       const double tol=1e-12)
{
    #ifndef RQTL2_NODEBUG
    if(X.rows() != y.size())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const VectorXd yy(as<Map<VectorXd> >(y));

    typedef ColPivHouseholderQR<MatrixXd> CPivQR;
    typedef CPivQR::PermutationType Permutation;

    const int p = XX.cols();

    CPivQR PQR( XX );
    PQR.setThreshold(tol); // set tolerance
    Permutation Pmat( PQR.colsPermutation() );
    const int r = PQR.rank();

    VectorXd betahat(p);

    if(r == p) { // full rank
        betahat = PQR.solve(yy);
    } else {
        MatrixXd Rinv( PQR.matrixQR().topLeftCorner(r,r).
                       triangularView<Upper>().
                       solve(MatrixXd::Identity(r,r)) );
        VectorXd effects( PQR.householderQ().adjoint() * yy );

        betahat.fill(::NA_REAL);
        betahat.head(r) = Rinv * effects.head(r);
        betahat = Pmat*betahat;
    }

    NumericVector result(wrap(betahat));
    return result;
}

// least squares by QR decomposition with column pivoting
// return the coefficients and SEs
// [[Rcpp::export]]
List calc_coefSE_linreg_eigenqr(const NumericMatrix& X, const NumericVector& y,
                                const double tol=1e-12)
{
    #ifndef RQTL2_NODEBUG
    if(X.rows() != y.size())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const VectorXd yy(as<Map<VectorXd> >(y));

    typedef ColPivHouseholderQR<MatrixXd> CPivQR;
    typedef CPivQR::PermutationType Permutation;

    const int n = XX.rows(), p = XX.cols();

    CPivQR PQR( XX );
    PQR.setThreshold(tol); // set tolerance
    Permutation Pmat( PQR.colsPermutation() );
    const int r = PQR.rank();

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
    const double rss = resid.squaredNorm();
    const int df = n - r;
    const double sigma = std::sqrt(rss/(double)df);

    return List::create(Named("coef") = betahat,
                        Named("SE") = sigma*se);
}

// least squares by QR decomposition with column pivoting
// return just the residual sum of squares
// [[Rcpp::export]]
double calc_rss_eigenqr(const NumericMatrix& X, const NumericVector& y,
                        const double tol=1e-12)
{
    #ifndef RQTL2_NODEBUG
    if(X.rows() != y.size())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const VectorXd yy(as<Map<VectorXd> >(y));

    typedef Eigen::ColPivHouseholderQR<MatrixXd> CPivQR;
    typedef CPivQR::PermutationType Permutation;

    const int n = XX.rows(), p = XX.cols();

    CPivQR PQR ( XX );
    PQR.setThreshold(tol); // set tolerance
    Permutation Pmat = PQR.colsPermutation();
    const int r = PQR.rank();

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


// least squares by QR decomposition with column pivoting
// return just the fitted values
// [[Rcpp::export]]
NumericVector calc_fitted_linreg_eigenqr(const NumericMatrix& X, const NumericVector& y,
                                         const double tol=1e-12)
{
    #ifndef RQTL2_NODEBUG
    if(X.rows() != y.size())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const VectorXd yy(as<Map<VectorXd> >(y));

    typedef Eigen::ColPivHouseholderQR<MatrixXd> CPivQR;
    typedef CPivQR::PermutationType Permutation;

    const int n = XX.rows(), p = XX.cols();

    CPivQR PQR ( XX );
    PQR.setThreshold(tol); // set tolerance
    Permutation Pmat = PQR.colsPermutation();
    const int r = PQR.rank();

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

    return wrap(fitted);
}


// least squares by "LLt" Cholesky decomposition, with matrix Y
// return vector of RSS
// [[Rcpp::export]]
NumericVector calc_mvrss_eigenchol(const NumericMatrix& X, const NumericMatrix& Y)
{
    #ifndef RQTL2_NODEBUG
    if(X.rows() != Y.rows())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    const int ncolY = Y.cols();
    const int ncolX = X.cols();

    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const MatrixXd YY(as<Map<MatrixXd> >(Y));

    LLT<MatrixXd> llt ( calc_XpX(XX) );

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
NumericVector calc_mvrss_eigenqr(const NumericMatrix& X, const NumericMatrix& Y,
                                 const double tol=1e-12)
{
    #ifndef RQTL2_NODEBUG
    if(X.rows() != Y.rows())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const MatrixXd YY(as<Map<MatrixXd> >(Y));

    typedef Eigen::ColPivHouseholderQR<MatrixXd> CPivQR;
    typedef CPivQR::PermutationType Permutation;

    const int n = XX.rows(), p = XX.cols();
    const int k = YY.cols();

    CPivQR PQR ( XX );
    PQR.setThreshold(tol); // set tolerance
    Permutation Pmat = PQR.colsPermutation();
    const int r = PQR.rank();

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
    #ifndef RQTL2_NODEBUG
    if(X.rows() != Y.rows())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    const int ncolY = Y.cols();
    const int ncolX = X.cols();

    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const MatrixXd YY(as<Map<MatrixXd> >(Y));

    LLT<MatrixXd> llt ( calc_XpX(XX) );

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
NumericMatrix calc_resid_eigenqr(const NumericMatrix& X, const NumericMatrix& Y,
                                 const double tol=1e-12)
{
    #ifndef RQTL2_NODEBUG
    if(X.rows() != Y.rows())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const MatrixXd YY(as<Map<MatrixXd> >(Y));

    typedef Eigen::ColPivHouseholderQR<MatrixXd> CPivQR;
    typedef CPivQR::PermutationType Permutation;

    const int n = XX.rows(), p = XX.cols();
    const int k = YY.cols();

    CPivQR PQR ( XX );
    PQR.setThreshold(tol); // set tolerance
    Permutation Pmat = PQR.colsPermutation();
    const int r = PQR.rank();

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
