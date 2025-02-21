#include <commons.hpp>
#include <iostream>
#include <lapack.h>
#include <iomanip>
int main()
{
    csr_matrix_sym matrix(8,{0,1,2,3,6,9,11,12},
    {0,1,2,1,2,3,0,1,4,0,5,6,1,2,6,7},{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1});

    for(size_t i = 0; i < matrix.nrows(); i++)
    {
        for(size_t j = 0; j < matrix.ncols(); j++)
            std::cout << std::setw(3) << std::setprecision(3) << matrix(i, j) << " ";
        std::cout << std::endl;
    }
    print_etree_ascii(build_etree(matrix));
}