#include <lapack.h>
#include <string>
extern void error_handler(std::string message);

void backward_sobstitute(const csr_matrix &L_upper, const std::vector<double> &b,
                         std::vector<double> &x)
{
    const size_t * row_pointers = L_upper.get_row_pointers();
    const size_t * column_indexes = L_upper.get_column_indexes();
    const double * data = L_upper.get_data();

    for (size_t i = L_upper.nrows(); i-- > 0; )
    {
        x[i] = b[i];
        for(size_t j_cursor =row_pointers[i]+1; j_cursor < row_pointers[i+1]; ++j_cursor)
        {
            x[i] -= data[j_cursor] * x[column_indexes[j_cursor]];
        }
        x[i]/= L_upper(i,i);
    }
}

void forward_sobsitute(const csr_matrix &L_lower, const std::vector<double> &b,
                       std::vector<double> &x)
{
    const size_t * row_pointers = L_lower.get_row_pointers();
    const size_t * column_indexes = L_lower.get_column_indexes();
    const double * data = L_lower.get_data();

    for (size_t i = 0; i < L_lower.nrows(); ++i)
    {
        x[i] = b[i];
        for (size_t j_cursor = row_pointers[i]; j_cursor < row_pointers[i+1] - 1; ++j_cursor)
        {
            x[i] -= data[j_cursor] * x[column_indexes[j_cursor]];
        }
        x[i] /= L_lower(i,i);
    }
}

void __symbolic_phase();
csr_matrix cholesky_decompose(const csr_matrix_sym &A)
{
    
}
