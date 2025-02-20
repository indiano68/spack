#include <blas.h>
#include <commons.hpp>
#include <fstream>
#include <iostream>
#include <matrix.h>
#include <tools.h>

int main(int argc, char **argv)
{
    std::string prefix = "";
    if(argc>1) prefix = std::string(argv[1]) + "/";
    csr_matrix A = build_block_diagonal_sparsity(3, 1000, 1);
    A.set_data(create_random_vector(A.r_size()));
    A.get_data()[0]=1;
    std::vector<double> x = create_random_vector(A.ncols());
    std::vector<double> y(A.nrows(), 0);
    gemv(A, x, y);
    matrix_to_file(A, prefix + "test_gemv_A.txt");
    vector_to_file(x, prefix + "test_gemv_x.txt");
    vector_to_file(y, prefix + "test_gemv_y.txt");
}

