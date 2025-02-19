#include <fstream>
#include <iostream>

#include <blas.h>
#include <matrix.h>
#include <tools.h>
std::vector<double> create_random_vector(size_t n)
{
    std::vector<double> vec(n);
    for (size_t i = 0; i < n; i++) {
        vec[i] = (double)rand() / RAND_MAX;
    }
    return vec;
}
std::vector<double> create_seq_vector(size_t n)
{
    std::vector<double> vec(n);
    for (size_t i = 0; i < n; i++) {
        vec[i] = (double)i;
    }
    return vec;
}

int main()
{
    csr_matrix A = build_block_diagonal_sparsity(3, 10, 1);
    A.set_data(create_random_vector(A.r_size()));
    std::vector<double> x = create_random_vector(A.ncols());
    std::vector<double> y(A.nrows(), 0);
    gemv(A, x, y);
    std::ofstream file_A("test_gemv_A.txt");
    for (size_t i = 0; i < A.nrows(); i++) {
        for (size_t j = 0; j < A.ncols(); j++) {
            file_A << A(i, j);
            if (j < A.ncols() - 1) {
                file_A << " ";
            }
        }
        file_A << std::endl;
    }
    file_A.close();

    std::ofstream file_x("test_gemv_x.txt");
    for (size_t i = 0; i < y.size(); i++) {
        file_x << x[i] << std::endl;
    }
    file_x.close();

    std::ofstream file_y("test_gemv_y.txt");
    for (size_t i = 0; i < y.size(); i++) {
        file_y << y[i] << std::endl;
    }
    file_y.close();
}