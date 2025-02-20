#include <blas.h>

extern void error_handler(std::string error_message);

void gemv(const csr_matrix &A, const std::vector<double> &x,
          std::vector<double> &y)
{

    if (A.nrows() != y.size()) {
        error_handler("Matrix and vector dimensions do not match");
    }
    if (A.ncols() != x.size()) {
        error_handler("Matrix and vector dimensions do not match");
    }
    const size_t *row_pointers = A.get_row_pointers();
    const size_t *column_indexes = A.get_column_indexes();
    const double *A_data = A.get_data();
    const double *x_data = x.data();
    double *y_data = y.data();
    for (size_t row_idx = 0; row_idx < A.nrows(); row_idx++) {
        y_data[row_idx] = 0;
        for (size_t i = row_pointers[row_idx]; i < row_pointers[row_idx + 1];
             i++) {
            y_data[row_idx] += A_data[i] * x_data[column_indexes[i]];
        }
    }
}

