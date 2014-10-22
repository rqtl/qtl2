// Matrix utilities

// cbind two or three matrices
IntegerMatrix cbind_imatrix(const IntegerMatrix& mat1, const IntegerMatrix& mat2);
IntegerMatrix cbind_3imatrix(const IntegerMatrix& mat1, const IntegerMatrix& mat2,
                             const IntegerMatrix& mat3);
NumericMatrix cbind_nmatrix(const NumericMatrix& mat1, const NumericMatrix& mat2);
NumericMatrix cbind_3nmatrix(const NumericMatrix& mat1, const NumericMatrix& mat2,
                             const NumericMatrix& mat3);
