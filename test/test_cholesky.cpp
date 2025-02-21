#include <commons.hpp>
#include <iomanip>
#include <iostream>
#include <lapack.h>

int main(int argc, char **argv)
{
    std::string fiename;
    if (argc > 1)
        fiename = argv[1];
    else
        std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;

    csr_matrix_sym matrix = load_csr_matrix_sym(fiename);

    for (size_t i = 0; i < matrix.nrows(); i++) {
        for (size_t j = 0; j < matrix.ncols(); j++)
            std::cout << std::setw(3) << std::setprecision(3) << matrix(i, j)
                      << " ";
        std::cout << std::endl;
    }
    auto etree = build_etree_pc(matrix);
    auto L_layout = build_L_layout(matrix, etree);
    auto L = csr_matrix(matrix.nrows(), matrix.ncols(), std::get<0>(L_layout),
                   std::get<1>(L_layout), std::vector<double>(std::get<1>(L_layout).size(),1));
    std::cout << "L:" << std::endl;
    cholesky_decompose(matrix, L);
    for (size_t i = 0; i < L.nrows(); i++) {
        for (size_t j = 0; j < L.ncols(); j++)
            std::cout << std::setw(3) << std::setprecision(3) << L(i, j)
                      << " ";
        std::cout << std::endl;
    }
}