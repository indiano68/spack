#include <lapack.h>
#include <string>
extern void error_handler(std::string message);

void backward_sobstitute(const csr_matrix &L_upper,
                         const std::vector<double> &b, std::vector<double> &x)
{
    const size_t *row_pointers = L_upper.get_row_pointers();
    const size_t *column_indexes = L_upper.get_column_indexes();
    const double *data = L_upper.get_data();

    for (size_t i = L_upper.nrows(); i-- > 0;) {
        x[i] = b[i];
        for (size_t j_cursor = row_pointers[i] + 1;
             j_cursor < row_pointers[i + 1]; ++j_cursor) {
            x[i] -= data[j_cursor] * x[column_indexes[j_cursor]];
        }
        x[i] /= L_upper(i, i);
    }
}

void forward_sobsitute(const csr_matrix &L_lower, const std::vector<double> &b,
                       std::vector<double> &x)
{
    const size_t *row_pointers = L_lower.get_row_pointers();
    const size_t *column_indexes = L_lower.get_column_indexes();
    const double *data = L_lower.get_data();

    for (size_t i = 0; i < L_lower.nrows(); ++i) {
        x[i] = b[i];
        for (size_t j_cursor = row_pointers[i];
             j_cursor < row_pointers[i + 1] - 1; ++j_cursor) {
            x[i] -= data[j_cursor] * x[column_indexes[j_cursor]];
        }
        x[i] /= L_lower(i, i);
    }
}

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

std::vector<size_t> build_etree(const csr_matrix_sym &A)
{
    std::vector<size_t> etree(A.ncols());
    const size_t *row_pointers = A.get_row_pointers();
    const size_t *column_indexes = A.get_column_indexes();
    for (size_t i = 0; i < A.ncols(); ++i) {
        etree[i] = static_cast<size_t>(-1);
        for (size_t j_cursor = row_pointers[i];
             j_cursor < row_pointers[i + 1] - 1; ++j_cursor) {
            uint64_t j_root = column_indexes[j_cursor];
            while (etree[j_root] != static_cast<size_t>(-1) &&
                   etree[j_root] < i) {
                j_root = etree[j_root];
            }
            if (etree[j_root] == static_cast<size_t>(-1)) {
                etree[j_root] = i;
            }
        }
    }
    return etree;
}
