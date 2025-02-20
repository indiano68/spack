#include <commons.hpp>
#include <iostream>
#include <iomanip>

int main(int argc, char **argv)
{
    std::string fiename;
    if (argc > 1)
        fiename = argv[1];
    else
        std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;

    csr_matrix_sym matrix = load_csr_matrix_sym(fiename);
    
    for(size_t i = 0; i < matrix.nrows(); i++)
    {
        for(size_t j = 0; j < matrix.ncols(); j++)
            std::cout << std::setw(3) << std::setprecision(3) << matrix(i, j) << " ";
        std::cout << std::endl;
    }
}
