#ifndef _SPACK_LAPACK_H_
#define _SPACK_LAPACK_H_
#include <matrix.h>
#include <vector>
#include <tuple>

void backward_sobstitute(const csr_matrix &L_upper,
                         const std::vector<double> &b, std::vector<double> &x);
void forward_sobsitute(const csr_matrix &L_lower, const std::vector<double> &b,
                       std::vector<double> &x);
std::vector<size_t> build_etree(const csr_matrix_sym &A);
std::vector<size_t> build_etree_pc(const csr_matrix_sym &A);
std::tuple<std::vector<size_t>, std::vector<size_t>> build_L_layout(const csr_matrix_sym &A, const std::vector<size_t> &etree);
void cholesky_decompose( const csr_matrix_sym  &A, csr_matrix &L);
#endif