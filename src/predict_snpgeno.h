// get predicted SNP genotypes from inferred genotypes + founder genotypes
#ifndef QTL2_PREDICT_SNPGENO
#define QTL2_PREDICT_SNPGENO


#include <Rcpp.h>

Rcpp::IntegerMatrix predict_snpgeno(const Rcpp::IntegerMatrix& allele1,
                                    const Rcpp::IntegerMatrix& allele2,
                                    const Rcpp::IntegerMatrix& founder_geno);


#endif // QTL2_PREDICT_SNPGENO
