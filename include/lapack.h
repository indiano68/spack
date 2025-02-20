#ifndef _SPACK_LAPACK_H_
#define _SPACK_LAPACK_H_
#include <matrix.h>
#include <vector>

void backward_sobstitute(const csr_matrix &L_upper, const std::vector<double> &b, std::vector<double> &x);
void forward_sobsitute(const csr_matrix &L_lower, const std::vector<double> &b, std::vector<double> &x);

#endif