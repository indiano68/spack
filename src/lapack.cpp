#include <cmath>
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

std::vector<size_t> build_etree_pc(const csr_matrix_sym &A)
{
    std::vector<size_t> etree(A.ncols());
    std::vector<size_t> anc_etree(A.ncols());
    const size_t *row_pointers = A.get_row_pointers();
    const size_t *column_indexes = A.get_column_indexes();
    for (size_t i = 0; i < A.ncols(); ++i) {
        etree[i] = static_cast<size_t>(-1);
        anc_etree[i] = static_cast<size_t>(-1);
        for (size_t j_cursor = row_pointers[i];
             j_cursor < row_pointers[i + 1] - 1; ++j_cursor) {
            uint64_t j_root = column_indexes[j_cursor];
            while (anc_etree[j_root] != static_cast<size_t>(-1) &&
                   anc_etree[j_root] < i) {
                size_t buffer = anc_etree[j_root];
                anc_etree[j_root] = i;
                j_root = buffer;
            }
            if (etree[j_root] == static_cast<size_t>(-1)) {
                etree[j_root] = i;
                anc_etree[j_root] = i;
            }
        }
    }
    return etree;
}

std::tuple<std::vector<size_t>, std::vector<size_t>>
build_L_layout(const csr_matrix_sym &A, const std::vector<size_t> &etree)
{
    std::vector<size_t> L_row_pointers(A.nrows() + 1);
    std::vector<size_t> L_column_indexes;
    const size_t *row_pointers = A.get_row_pointers();
    const size_t *column_indexes = A.get_column_indexes();
    std::vector<size_t> mark(A.nrows());
    for (size_t i = 0; i < A.nrows(); ++i) {
        size_t added = 1;
        mark[i] = i;
        for (size_t j_cursor = row_pointers[i];
             j_cursor < row_pointers[i + 1] - 1; ++j_cursor) {
            size_t j = column_indexes[j_cursor];
            while (mark[j] != i) {
                L_column_indexes.push_back(j);
                mark[j] = i;
                j = etree[j];
                ++added;
            }
        }
        L_column_indexes.push_back(i);
        L_row_pointers[i + 1] = L_row_pointers[i] + added;
    }
    L_row_pointers.pop_back();
    return std::make_tuple(L_row_pointers, L_column_indexes);
}

void cholesky_decompose(const csr_matrix_sym &A, csr_matrix &L)
{
    double *L_data = L.get_data();
    const size_t *L_row_pointers = L.get_row_pointers();
    const size_t *L_column_indexes = L.get_column_indexes();
    const double * A_data = A.get_data();
    const size_t * A_row_pointers = A.get_row_pointers();
    const size_t * A_column_indexes = A.get_column_indexes();
    for (size_t i = 0; i < L.nrows(); i++) {
        for (size_t j_cursor = L_row_pointers[i];
             j_cursor < L_row_pointers[i + 1]; ++j_cursor) {
            size_t j = L_column_indexes[j_cursor];
            double sum = 0;
            size_t k_p_cursor = L_row_pointers[j];
            size_t j_a_cursor = A_row_pointers[i];
            for (size_t k_cursor = L_row_pointers[i]; k_cursor < j_cursor; k_cursor++)
            {
                for (; k_p_cursor < L_row_pointers[j+1]-1; k_p_cursor++)
                {
                   if (L_column_indexes[k_p_cursor] == L_column_indexes[k_cursor]) {
                        sum += L_data[k_cursor] * L_data[k_p_cursor]; 
                        break;
                   }
                   else if (L_column_indexes[k_p_cursor] > L_column_indexes[k_cursor])
                       break;
                }
                //    sum += L_data[k_cursor] * L(j, L_column_indexes[k_cursor]); 
            }
            if (i == j)
                L_data[j_cursor] = sqrt(A_data[A_row_pointers[i+1]-1] - sum);
            else
            for (; j_a_cursor < A_row_pointers[i+1]; j_a_cursor++)
            {
               if (L_column_indexes[j_cursor] == A_column_indexes[j_a_cursor]) {
                    L_data[j_cursor] =
                        ((1.0 / L[L_row_pointers[j+1]-1]) * (A_data[j_a_cursor] - sum)); 
                    break;
               }
               
               if (L_column_indexes[j_cursor] < A_column_indexes[j_a_cursor])
               {
                L_data[j_cursor] =
                    ((1.0 / L[L_row_pointers[j+1]-1]) * (-1*sum)); 
                   break;
               }
                
            }
                // L_data[j_cursor] =
                //     ((1.0 / L[L_row_pointers[j+1]-1]) * (A(i, j) - sum)); 
        }
    }
}
