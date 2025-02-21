#include <commons.hpp>
#include <iomanip>
#include <iostream>
#include <lapack.h>

int main(int argc, char **argv)
{
    std::string fiename;
    std::string prefix;
    if (argc > 2)
    {
        fiename = argv[1];
        prefix = argv[2];
        prefix+="/";
    }
    else
        std::cout << "Usage: " << argv[0] << " <filename>" << "<prefix>" << std::endl;
    csr_matrix_sym matrix = load_csr_matrix_sym(fiename);
    auto etree = build_etree_pc(matrix);
    auto L_layout = build_L_layout(matrix, etree);
    auto L = csr_matrix(matrix.nrows(), matrix.ncols(), std::get<0>(L_layout),
                   std::get<1>(L_layout), std::vector<double>(std::get<1>(L_layout).size(),1));
    std::cout << "L:" << std::endl;
    cholesky_decompose(matrix, L);
    matrix_to_file(L, prefix + "L.txt");
}