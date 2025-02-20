#include <cmath>
#include <matrix.h>
#include <tools.h>
#include <vector>

extern void error_handler(std::string error_message);
csr_matrix build_upper_triangular_sparsity(size_t n)
{
    std::vector<size_t> row_pointers(n);
    std::vector<size_t> col_indexes((size_t)((n * (n + 1)) / 2));
    std::vector<double> data((size_t)((n * (n + 1)) / 2), 1.0);
    row_pointers[0] = 0;

    for (size_t i = 0; i < n - 1; i++) {
        row_pointers[i + 1] = row_pointers[i] + n - i;
    }

    for (size_t i = 0; i < n; i++) {
        for (size_t j = i; j < n; j++) {
            col_indexes[row_pointers[i] + (j - i)] = j;
        }
    }
    return csr_matrix(n, n, row_pointers, col_indexes, data);
}
csr_matrix build_lower_triangular_sparsity(size_t n)
{
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
    return csr_matrix(n, n, row_pointers, col_indexes, data);
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
    } else {
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

csr_matrix_sym build_block_diagonal_sparsity_sym(size_t block_n, size_t n_blocks,
                                             size_t overlap_length)
{
    // Overlap must be strictly smaller than the block size.
    if (overlap_length >= block_n) {
        error_handler("Overlap length must be smaller than block size");
    }

    // --- Case 1: No overlap ---
    if (overlap_length == 0) {
        // Global matrix is simply n_blocks independent dense blocks.
        size_t total_rows = block_n * n_blocks;
        // For each block, row r in that block (with local index r_local) stores
        // only the upper-triangular part, i.e. block_n - r_local entries.
        size_t total_nnz = 0;
        for (size_t r = 0; r < total_rows; r++) {
            size_t local = r % block_n;
            total_nnz += block_n - local;
        }

        std::vector<size_t> row_pointers(total_rows, 0);
        std::vector<size_t> col_indexes(total_nnz, 0);
        std::vector<double> data(total_nnz, 1.0);

        size_t nnz_counter = 0;
        for (size_t r = 0; r < total_rows; r++) {
            row_pointers[r] = nnz_counter;
            size_t block = r / block_n;
            size_t local = r % block_n;
            size_t block_start = block * block_n;
            // In the dense block, store only j>=local.
            for (size_t j = local; j < block_n; j++) {
                col_indexes[nnz_counter + (j - local)] = block_start + j;
            }
            nnz_counter += block_n - local;
        }

        return csr_matrix_sym(total_rows, row_pointers, col_indexes,
                          data);
    }
    // --- Case 2: Nonzero overlap ---
    else {
        // Global rows: the first block contributes block_n rows;
        // each subsequent block adds (block_n - overlap_length) new rows.
        size_t total_rows =
            block_n + (n_blocks - 1) * (block_n - overlap_length);

        std::vector<size_t> row_pointers(total_rows, 0);
        size_t total_nnz = 0;
        const size_t block_step =
            block_n - overlap_length; // vertical shift per block

        // First pass: determine, for each row, how many columns will be stored.
        // In the fully stored matrix, the union of nonzeros for row r is
        // computed by:
        //   b_min = smallest block index that covers row r
        //   b_max = largest block index that covers row r (capped at
        //   n_blocks-1)
        // and the full union covers columns
        //   [union_start, union_end) where union_start = b_min * block_step,
        //   and union_end = b_max * block_step + block_n.
        // For symmetric storage we keep only columns j >= r.
        for (size_t r = 0; r < total_rows; r++) {
            size_t b_min = 0;
            if (r >= block_n) {
                double tmp =
                    static_cast<double>(r - (block_n - 1)) / block_step;
                b_min = static_cast<size_t>(std::ceil(tmp));
            }
            double tmp2 = static_cast<double>(r) / block_step;
            size_t b_max = static_cast<size_t>(std::floor(tmp2));
            if (b_max >= n_blocks)
                b_max = n_blocks - 1;

            size_t union_start = b_min * block_step;
            size_t union_end =
                b_max * block_step + block_n; // one past the last col index

            // Only store entries with col index >= r.
            size_t effective_start = (r > union_start ? r : union_start);
            size_t sym_length =
                (union_end > effective_start ? union_end - effective_start : 0);

            row_pointers[r] = total_nnz;
            total_nnz += sym_length;
        }

        std::vector<size_t> col_indexes(total_nnz, 0);
        std::vector<double> data(total_nnz, 1.0);

        size_t nnz_counter = 0;
        // Second pass: fill in the column indices.
        for (size_t r = 0; r < total_rows; r++) {
            size_t b_min = 0;
            if (r >= block_n) {
                double tmp =
                    static_cast<double>(r - (block_n - 1)) / block_step;
                b_min = static_cast<size_t>(std::ceil(tmp));
            }
            double tmp2 = static_cast<double>(r) / block_step;
            size_t b_max = static_cast<size_t>(std::floor(tmp2));
            if (b_max >= n_blocks)
                b_max = n_blocks - 1;

            size_t union_start = b_min * block_step;
            size_t union_end = b_max * block_step + block_n;
            size_t effective_start = (r > union_start ? r : union_start);
            size_t sym_length =
                (union_end > effective_start ? union_end - effective_start : 0);

            for (size_t j = 0; j < sym_length; j++) {
                col_indexes[nnz_counter + j] = effective_start + j;
            }
            nnz_counter += sym_length;
        }

        return csr_matrix_sym(total_rows, row_pointers, col_indexes,
                          data);
    }
}