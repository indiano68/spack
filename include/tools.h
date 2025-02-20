/********************************************************************
    File: tools.h
    Author: Erik Fabrizzi
    Contains tools for the manipulation of sparse matrices
**********************************************************************/
#ifndef _SPACK_TOOLS_H_
#define _SPACK_TOOLS_H_
#include <matrix.h>
csr_matrix build_upper_triangular_sparsity(size_t n);
csr_matrix build_lower_triangular_sparsity(size_t n);
csr_matrix build_block_diagonal_sparsity(size_t block_n, size_t n_blocks,
                                         size_t overlap_length);
csr_matrix_sym build_block_diagonal_sparsity_sym(size_t block_n, size_t n_blocks,
                                         size_t overlap_length);
#endif