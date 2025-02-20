#include <iostream>
#include <matrix.h>
#include <commons.hpp>
#include <tools.h>
#include <lapack.h>


int main(int argc, char **argv)
{
    std::string prefix = "";
    if(argc>1) prefix = std::string(argv[1]) + "/";
    csr_matrix upper_triangular = build_upper_triangular_sparsity(1000);
    std::vector<double> b = create_seq_vector(upper_triangular.nrows());
    std::vector<double> x(upper_triangular.ncols(),0);
    upper_triangular.set_data(create_random_vector(upper_triangular.r_size()));
    backward_sobstitute(upper_triangular, b, x);
    matrix_to_file(upper_triangular, prefix + "test_upper_triangular.txt");
    vector_to_file(b, prefix + "test_b.txt");
    vector_to_file(x, prefix + "test_x.txt");
    return 0;
}