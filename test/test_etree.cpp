#include <commons.hpp>
#include <iomanip>
#include <iostream>
#include <lapack.h>
int main()
{
    csr_matrix_sym matrix(8, {0, 1, 2, 3, 6, 9, 11, 12},
                          {0, 1, 2, 1, 2, 3, 0, 1, 4, 0, 5, 6, 1, 2, 6, 7},
                          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});

    for (size_t i = 0; i < matrix.nrows(); i++) {
        for (size_t j = 0; j < matrix.ncols(); j++)
            std::cout << std::setw(3) << std::setprecision(3) << matrix(i, j)
                      << " ";
        std::cout << std::endl;
    }
    auto etree = build_etree_pc(matrix);
    auto L_layout = build_L_layout(matrix, etree);
    csr_matrix L =
        csr_matrix(8, 8, std::get<0>(L_layout), std::get<1>(L_layout),
                   std::vector<double>(std::get<1>(L_layout).size(), 1));
    std::cout << "\n\n" << std::endl;
    for (size_t i = 0; i < L.nrows(); i++) {
        for (size_t j = 0; j < L.ncols(); j++)
            std::cout << std::setw(3) << std::setprecision(3) << L(i, j) << " ";
        std::cout << std::endl;
    }
}