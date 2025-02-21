#include <iomanip>
#include <iostream>

#include <commons.hpp>
#include <matrix.h>
#include <tools.h>

int main()
{
    csr_matrix_sym test_matrix(build_block_diagonal_sparsity_sym(4, 5, 2));
    test_matrix.set_data(create_seq_vector(test_matrix.r_size()));
    for (size_t i = 0; i < test_matrix.ncols(); i++) {
        for (size_t j = 0; j < test_matrix.ncols(); j++) {
            std::cout << std::fixed << std::setprecision(2) << test_matrix(i, j)
                      << " ";
        }
        std::cout << std::endl;
    }
}