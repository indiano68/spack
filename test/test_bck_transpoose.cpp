#include <commons.hpp>
#include <tools.h>
#include <lapack.h>


int main(int argc, char **argv)
{
    std::string prefix = "";
    if(argc>1) prefix = std::string(argv[1]) + "/";
    csr_matrix lower_triangular = build_lower_triangular_sparsity(1000);
    std::vector<double> b = create_random_vector(lower_triangular.nrows());
    std::vector<double> x(lower_triangular.ncols(),0);
    lower_triangular.set_data(create_random_vector(lower_triangular.r_size()));
    backward_sobstitute_transpose(lower_triangular, b, x);
    std::cout << "exit_success" << std::endl;
    matrix_to_file(lower_triangular, prefix + "test_upper_triangular.txt");
    vector_to_file(b, prefix + "test_b.txt");
    vector_to_file(x, prefix + "test_x.txt");
    return 0;
}