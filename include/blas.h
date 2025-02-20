#ifndef _SPACK_BLAS_H_
#define _SPACK_BLAS_H_
#include<matrix.h>
#include<vector>
void gemv(const csr_matrix&  A, const std::vector<double>& x, std::vector<double>& y);
void gemv(const csr_matrix_sym&  A, const std::vector<double>& x, std::vector<double>& y);
#endif