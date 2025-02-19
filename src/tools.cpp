#include <cmath>
#include <matrix.h>
#include <tools.h>
#include <vector>

extern void error_handler(std::string error_message);
csr_matrix build_upper_triangular_sparsity(size_t n)
{
    csr_matrix csr = csr_matrix(n, n, (size_t)((n * (n + 1)) / 2));
    std::vector<size_t> row_pointers(n);
    std::vector<size_t> col_indexes((size_t)((n * (n + 1)) / 2));
    std::vector<double> data((size_t)((n * (n + 1)) / 2), 1.0);
    row_pointers[0] = 0;

    for (size_t i = 0; i < n - 1; i++) {
        row_pointers[i + 1] = row_pointers[i] + n - i;
    }

    for (size_t i = 0; i < n; i++) {
        for (size_t j = i; j < csr.ncols(); j++) {
            col_indexes[row_pointers[i] + (j - i)] = j;
        }
    }
    csr.set_row_pointers(row_pointers);
    csr.set_column_indexes(col_indexes);
    csr.set_data(data);
    return csr;
}
csr_matrix build_lower_triangular_sparsity(size_t n)
{
    csr_matrix csr = csr_matrix(n, n, (size_t)((n * (n + 1)) / 2));
    std::vector<size_t> row_pointers(n);
    std::vector<size_t> col_indexes((size_t)((n * (n + 1)) / 2));
    std::vector<double> data((size_t)((n * (n + 1)) / 2), 1.0);
    row_pointers[0] = 0;
    for (size_t i = 0; i < n - 1; i++) {
        row_pointers[i + 1] = row_pointers[i] + i + 1;
    }

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < i + 1; j++) {
            col_indexes[row_pointers[i] + j] = j;
        }
    }

    csr.set_row_pointers(row_pointers);
    csr.set_column_indexes(col_indexes);
    csr.set_data(data);
    return csr;
}

csr_matrix build_block_diagonal_sparsity(size_t block_n, size_t n_blocks,
                                         size_t overlap_length)
/*TODO: this is overly complicated and must be revised*/
{
    if (overlap_length >= block_n) {
        error_handler("Overlap length must be smaller than block size");
    }

    if (overlap_length == 0) {
        size_t total_rows = block_n * n_blocks;
        size_t total_nnz = total_rows * block_n;

        std::vector<size_t> row_pointers(total_rows, 0);
        std::vector<size_t> col_indexes(total_nnz, 0);
        std::vector<double> data(total_nnz, 1.0);

        for (size_t r = 0; r < total_rows; r++) {
            row_pointers[r] = r * block_n;

            size_t block = r / block_n;
            size_t block_start_col = block * block_n;
            for (size_t j = 0; j < block_n; j++) {
                col_indexes[r * block_n + j] = block_start_col + j;
            }
        }

        return csr_matrix(total_rows, total_rows, row_pointers, col_indexes,
                          data);
    }

    else {
        size_t total_rows =
            block_n + (n_blocks - 1) * (block_n - overlap_length);

        std::vector<size_t> row_pointers(total_rows, 0);
        size_t total_nnz = 0;
        const size_t block_step = block_n - overlap_length;

        for (size_t r = 0; r < total_rows; r++) {
            size_t b_min = 0;
            if (r >= block_n) {
                double tmp = (double)(r - (block_n - 1)) / block_step;
                b_min = static_cast<size_t>(std::ceil(tmp));
            }

            double tmp2 = (double)r / block_step;
            size_t b_max = static_cast<size_t>(std::floor(tmp2));
            if (b_max >= n_blocks)
                b_max = n_blocks - 1;
            size_t union_length =
                (b_max > b_min ? (b_max - b_min) * block_step : 0) + block_n;
            row_pointers[r] = total_nnz;
            total_nnz += union_length;
        }

        std::vector<size_t> col_indexes(total_nnz, 0);
        std::vector<double> data(total_nnz, 1.0);

        size_t nnz_counter = 0;
        for (size_t r = 0; r < total_rows; r++) {
            size_t b_min = 0;
            if (r >= block_n) {
                double tmp = (double)(r - (block_n - 1)) / block_step;
                b_min = static_cast<size_t>(std::ceil(tmp));
            }
            double tmp2 = (double)r / block_step;
            size_t b_max = static_cast<size_t>(std::floor(tmp2));
            if (b_max >= n_blocks)
                b_max = n_blocks - 1;

            size_t union_length =
                (b_max > b_min ? (b_max - b_min) * block_step : 0) + block_n;
            size_t union_start = b_min * block_step;

            for (size_t j = 0; j < union_length; j++) {
                col_indexes[nnz_counter + j] = union_start + j;
            }
            nnz_counter += union_length;
        }

        return csr_matrix(total_rows, total_rows, row_pointers, col_indexes,
                          data);
    }
}
