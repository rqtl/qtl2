// Matrix utilities
#ifdef MATRIX_H
#define MATRIX_H

// cbind two or three matrices
IntegerMatrix cbind_imatrix(const IntegerMatrix& mat1, const IntegerMatrix& mat2);
IntegerMatrix cbind_3imatrix(const IntegerMatrix& mat1, const IntegerMatrix& mat2,
                             const IntegerMatrix& mat3);
NumericMatrix cbind_nmatrix(const NumericMatrix& mat1, const NumericMatrix& mat2);
NumericMatrix cbind_3nmatrix(const NumericMatrix& mat1, const NumericMatrix& mat2,
                             const NumericMatrix& mat3);

// rbind two or three matrices
IntegerMatrix rbind_imatrix(const IntegerMatrix& mat1, const IntegerMatrix& mat2);
IntegerMatrix rbind_3imatrix(const IntegerMatrix& mat1, const IntegerMatrix& mat2,
                             const IntegerMatrix& mat3);
NumericMatrix rbind_nmatrix(const NumericMatrix& mat1, const NumericMatrix& mat2);
NumericMatrix rbind_3nmatrix(const NumericMatrix& mat1, const NumericMatrix& mat2,
                             const NumericMatrix& mat3);

#endif // MATRIX_H
