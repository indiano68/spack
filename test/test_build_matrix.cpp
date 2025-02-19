#include <iomanip>
#include <iostream>
#include <vector>

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
    csr_matrix csr = build_upper_triangular_sparsity(10);
    // csr.set_data(create_random_vector(csr.r_size()));
    csr.set_data(create_seq_vector(csr.r_size()));
    for (size_t i = 0; i < csr.nrows(); i++) {
        for (size_t j = 0; j < csr.ncols(); j++)
            std::cout << std::fixed << std::setprecision(3) << csr(i, j) << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
    csr_matrix csr2 = build_lower_triangular_sparsity(10);
    csr2.set_data(create_seq_vector(csr2.r_size()));
    for (size_t i = 0; i < csr2.nrows(); i++) {
        for (size_t j = 0; j < csr2.ncols(); j++)
            std::cout << std::fixed << std::setprecision(3) << csr2(i, j) << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
    csr_matrix csr3 = build_block_diagonal_sparsity(4,5,3);
    csr3.set_data(create_seq_vector(csr3.r_size()));
    for (size_t i = 0; i < csr3.nrows(); i++) {
        for (size_t j = 0; j < csr3.ncols(); j++)
            std::cout << std::fixed << std::setprecision(2) << std::setw(2) << std::setfill(' ') << (int)csr3(i, j) << " ";
        std::cout << std::endl;
    }
}
