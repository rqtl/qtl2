// scan chromosome to get BLUPs of coefficients

#include "scan1blup.h"
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

#include "linreg.h"
#include "linreg_eigen.h"  // contains calc_XpX()
#include "lmm.h"
#include "matrix.h"

// Scan a single chromosome to get BLUPs of coefficients
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (must include intercept)
// se        = If TRUE, calculate SEs
// reml      = If TRUE, use REML to estimate variance components; otherwise use maximum
//             likelihood
// preserve_intercept = If FALSE, add the intercept to the BLUPs and remove that column
// tol       = Numeric tolerance
//
// output    = List with twomatrix of coefficients (genotypes x positions)
//
// [[Rcpp::export]]
List scanblup(const NumericVector& genoprobs,
              const NumericVector& pheno,
              const NumericMatrix& addcovar,
              const bool se,
              const bool reml,
              const bool preserve_intercept,
              const double tol=1e-12)
{
    const int n_ind = pheno.size();
    const Dimension d = genoprobs.attr("dim");
    const int n_pos = d[2];
    const int n_gen = d[1];
    const int n_addcovar = addcovar.cols();
    const int x_size = n_ind * n_gen;
    int n_coef = n_gen + n_addcovar;
    if(!preserve_intercept) n_coef--;
    NumericMatrix coef(n_coef, n_pos); // to contain the estimated coefficients
    NumericMatrix SE(n_coef, n_pos); // to contain the estimated SEs

    for(int pos=0; pos < n_pos; pos++) { // loop over positions
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // calculate ZZ' with Z = genoprobs for this position
        MatrixXd Z(n_ind, n_gen);
        std::copy(genoprobs.begin() + x_size*pos, genoprobs.begin() + x_size*(pos+1),
                  Z.data());
        MatrixXd ZZp = calc_XpX(Z.transpose());

        // get eigen decomposition
        std::pair<VectorXd, MatrixXd> ZZp_decomp = eigen_decomp(ZZp);

        // rotate pheno and X
        VectorXd ph(as <Map<VectorXd> >(pheno));
        MatrixXd ac(as <Map<MatrixXd> >(addcovar));
        VectorXd y = ZZp_decomp.second * ph; // D'y
        MatrixXd X = ZZp_decomp.second * ac; // D'X

        // fit LMM
        struct lmm_fit lmm_out = fitLMM(ZZp_decomp.first, y, X, reml, true, NA_REAL, tol);

        // calculate BLUPs
        VectorXd resid = y - X * lmm_out.beta;
        for(int i=0; i<n_ind; i++) // multiply by weights
            resid[i] *= lmm_out.hsq/(lmm_out.hsq * ZZp_decomp.first[i] + (1.0-lmm_out.hsq));
        Z = ZZp_decomp.second * Z;
        VectorXd blup = Z.transpose() * resid;

        // insert estimated coefficients
        if(!preserve_intercept) { // add intercept to BLUPs and drop the intercept column
            for(int i=0; i<n_gen; i++) coef(i,pos) = blup[i] + lmm_out.beta[0];
            for(int i=0; i<n_addcovar-1; i++) coef(n_gen+i,pos) = lmm_out.beta[i+1];
        }
        else {
            for(int i=0; i<n_gen; i++) coef(i,pos) = blup[i];
            for(int i=0; i<n_addcovar; i++) coef(n_gen+i,pos) = lmm_out.beta[i];
        }

        if(se) { // get SEs
            // construct variance matrix
            MatrixXd ZX(n_ind, n_gen+n_addcovar);
            for(int i=0; i<n_gen; i++) ZX.col(i) = Z.col(i);
            for(int i=0; i<n_addcovar; i++) ZX.col(i+n_gen) = ac.col(i);
            MatrixXd Vinv = calc_XpX(ZX)/(1.0-lmm_out.hsq);
            for(int i=0; i<n_gen; i++) Vinv(i,i) += 1.0/lmm_out.hsq;

            // sigma^2 * (Vinv)^(-1)
            MatrixXd V = Vinv.inverse()*lmm_out.sigmasq;

            // insert estimated SEs
            if(!preserve_intercept) { // add intercept to BLUPs and drop the intercept column
                for(int i=0; i<n_gen; i++) SE(i,pos) = sqrt(V(i,i) + V(n_gen,n_gen) + 2.0*V(i,n_gen));
                for(int i=0; i<n_addcovar-1; i++) SE(n_gen+i,pos) = sqrt(V(n_gen+i+1,n_gen+i+1));
            }
            else {
                for(int i=0; i<n_coef; i++) SE(i,pos) = sqrt(V(i,i));
            }

        }
    }

    return List::create(Named("coef") = coef,
                        Named("SE") = SE);
}
