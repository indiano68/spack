#include <algorithm>
#include <cstring>
#include <iostream>
#include <matrix.h>
#include <string>

void error_handler(std::string error_message)
{
    std::cerr << "ERROR:" << error_message << std::endl;
    std::exit(EXIT_FAILURE);
}
static void warn_handler(std::string warning_message)
{
    std::cerr << "WARNING:" << warning_message << std::endl;
}

/******************** Implement Matrix ********************/

size_t matrix_interface::size() const
{
    return __ncols * __nrows;
}
size_t matrix_interface::nrows() const
{
    return __nrows;
}
size_t matrix_interface::ncols() const
{
    return __ncols;
}

/******************* Implement Dense Matrix ********************/

dense_matrix::dense_matrix(size_t row_number, size_t col_number) noexcept
    : matrix_interface(row_number, col_number)
{
    __allocate();
}
void dense_matrix::__allocate()
{
    __data = std::make_unique<double[]>(__nrows * __ncols);
};
void dense_matrix::set_data(const std::vector<double> &data)
{
    if (__ncols * __nrows != data.size())
        error_handler("Data dimention incorrect in call to set_data (Dense)");
    std::memcpy(__data.get(), data.data(), data.size());
}
double *dense_matrix::get_data()
{
    return __data.get();
}
double &dense_matrix::operator()(size_t row_idx, size_t col_idx)
{
    return __data[__ncols * row_idx + col_idx];
}
double &dense_matrix::operator[](size_t idx)
{
    return __data[idx];
}
const double &dense_matrix::operator()(size_t row_idx, size_t col_idx) const
{
    return __data[__ncols * row_idx + col_idx];
}
const double &dense_matrix::operator[](size_t idx) const
{
    return __data[idx];
}

/****************** Implement Sparse Matrix ******************/

sparse_matrix_interface::sparse_matrix_interface(
    size_t row_number, size_t col_number, size_t non_zero_number) noexcept
    : matrix_interface(row_number, col_number), __n_non_zero(non_zero_number)
{
    __allocate();
}
size_t sparse_matrix_interface::r_size()
{
    return __n_non_zero;
}
void sparse_matrix_interface::__allocate()
{
    __data = std::make_unique<double[]>(__n_non_zero);
};

/******************* Implement CRF *******************/

coo_matrix::coo_matrix(std::size_t n_rows, std::size_t n_cols,
                       size_t n_non_zero) noexcept
    : sparse_matrix_interface(n_rows, n_cols, n_non_zero),
      __column_indexes(std::make_unique<size_t[]>(n_non_zero)),
      __row_indexes(std::make_unique<size_t[]>(n_non_zero)), __populated{false}
{
}
coo_matrix::coo_matrix(size_t n_rows, size_t n_cols,
                       std::vector<size_t> r_indexes,
                       std::vector<size_t> c_indexes,
                       std::vector<double> data) noexcept
    : sparse_matrix_interface(n_rows, n_cols, data.size())
{

    if (r_indexes.size() == c_indexes.size() &&
        c_indexes.size() == data.size()) {
        __column_indexes = std::make_unique<size_t[]>(__n_non_zero);
        __row_indexes = std::make_unique<size_t[]>(__n_non_zero);
        set_indexes(r_indexes, c_indexes);
        set_data(data);
    } else
        error_handler(
            "Incorrectly sized objects passed to crf_matrix_constructor");
}
static bool __check_order(const std::vector<size_t> &row_idx,
                          const std::vector<size_t> &col_idx)
{
    bool result = true;
    for (std::size_t i = 1; i < row_idx.size(); ++i) {
        if (row_idx[i] < row_idx[i - 1]) {
            result = false;
            break;
        }
        if (row_idx[i] == row_idx[i - 1] && col_idx[i] < col_idx[i - 1]) {
            result = false;
            break;
        }
    }
    return result;
}
void coo_matrix::set_indexes(const std::vector<size_t> &row_indexes,
                             const std::vector<size_t> &col_indexes)
{
    if (row_indexes.empty() || col_indexes.empty())
        error_handler("Some of the indexes are empty");
    if (row_indexes.size() != col_indexes.size() ||
        row_indexes.size() != __n_non_zero ||
        col_indexes.size() != __n_non_zero)
        error_handler(
            "One ar all index sizes not matching non-zero element number");
    if (!__check_order(row_indexes, col_indexes))
        error_handler(
            "Row and/or column indexes are not respecting the required "
            "forma (in-order)");
    std::memcpy(__row_indexes.get(), row_indexes.data(),
                __n_non_zero * sizeof(size_t));
    std::memcpy(__column_indexes.get(), col_indexes.data(),
                __n_non_zero * sizeof(size_t));
    __populated[0] = true;
    __populated[1] = true;
    __populated[3] = __populated[0] && __populated[1] && __populated[2];
}
void coo_matrix::set_data(const std::vector<double> &data)
{
    if (data.size() != __n_non_zero)
        error_handler(
            "Number of data entries does not match non-zero element number");
    ;
    std::memcpy(__data.get(), data.data(), __n_non_zero * sizeof(double));
    __populated[2] = true;
    __populated[3] = __populated[0] && __populated[1] && __populated[2];
}
size_t *coo_matrix::get_column_indexes()
{
    return __column_indexes.get();
}
size_t *coo_matrix::get_row_indexes()
{
    return __row_indexes.get();
}
double *coo_matrix::get_data()
{
    return __data.get();
}
bool coo_matrix::is_populated() const
{
    return __populated[3];
}
const double &coo_matrix::operator()(size_t row_idx, size_t col_idx) const
{
    static const double zero_dummy = 0;
    const double *data_ptr = &zero_dummy;
    if (row_idx >= __nrows)
        error_handler("Row index  out of bound");
    if (col_idx >= __ncols)
        error_handler("Col index  out of bound");
    if (!__populated[2])
        warn_handler("Accessing Uninitialized Data");
    for (size_t data_idx = 0; data_idx < __n_non_zero; data_idx++)
        if (__row_indexes[data_idx] == row_idx)
            if (__column_indexes[data_idx] == col_idx) {
                data_ptr = &__data[data_idx];
                break;
            }
    return *data_ptr;
}
const double &coo_matrix::operator[](size_t idx) const
{

    if (idx >= __n_non_zero)
        error_handler("Index  out of bound");
    if (!__populated[2])
        warn_handler("Accessing Uninitialized Data");
    return __data[idx];
}

/******************** Implement CRS ********************/

csr_matrix::csr_matrix(std::size_t n_rows, std::size_t n_cols,
                       size_t n_non_zero) noexcept
    : sparse_matrix_interface(n_rows, n_cols, n_non_zero),
      __row_pointers(std::make_unique<size_t[]>(n_rows + 1)),
      __col_indexes(std::make_unique<size_t[]>(n_non_zero))
{
    __row_pointers[n_rows] = n_non_zero;
}
csr_matrix::csr_matrix(size_t n_rows, size_t n_cols,
                       std::vector<size_t> r_pointers,
                       std::vector<size_t> c_indexes,
                       std::vector<double> data) noexcept
    : sparse_matrix_interface(n_rows, n_cols, data.size())
{
    if (r_pointers.size() == n_rows && c_indexes.size() == data.size()) {
        __row_pointers = std::make_unique<size_t[]>(n_rows+1);
        __col_indexes = std::make_unique<size_t[]>(c_indexes.size());
        __row_pointers[n_rows] = data.size();
        set_row_pointers(r_pointers);
        set_column_indexes(c_indexes);
        set_data(data);
    } else
        error_handler(
            "Incorrectly sized objects passed to csr_matrix_constructor");
}
void csr_matrix::set_column_indexes(const std::vector<size_t> &col_indexes)
{
    if (col_indexes.size() != __n_non_zero)
        error_handler("Incorrect number of column indexes");
    __populated[0] = true;
    __populated[3] = __populated[0] && __populated[1] && __populated[2];
    std::memcpy(__col_indexes.get(), col_indexes.data(),
                sizeof(size_t) * col_indexes.size());
}
void csr_matrix::set_row_pointers(const std::vector<size_t> &row_pointers)
{
    if (row_pointers.size() != __nrows)
        error_handler("Incorrect number of row pointers");
    __populated[1] = true;
    __populated[3] = __populated[0] && __populated[1] && __populated[2];
    std::memcpy(__row_pointers.get(), row_pointers.data(),
                sizeof(size_t) * row_pointers.size());
}
void csr_matrix::set_data(const std::vector<double> &data)
{
    if (data.size() != __n_non_zero)
        error_handler("Incorrect number of row pointers");
    __populated[2] = true;
    __populated[3] = __populated[0] && __populated[1] && __populated[2];

    std::memcpy(__data.get(), data.data(), sizeof(double) * data.size());
}
size_t *csr_matrix::get_column_indexes()
{
    return __col_indexes.get();
}
const size_t *csr_matrix::get_column_indexes() const
{
    return __col_indexes.get();
}
size_t *csr_matrix::get_row_pointers()
{
    return __row_pointers.get();
}
const size_t *csr_matrix::get_row_pointers() const 
{
    return __row_pointers.get();
}
double *csr_matrix::get_data()
{
    return __data.get();
}
const double *csr_matrix::get_data() const
{
    return __data.get();
}
bool csr_matrix::is_populated() const
{
    return __populated[3];
}
const double &csr_matrix::operator()(size_t row_idx, size_t col_idx) const
{
    static const double zero_dummy = 0;
    const double *data_ptr = &zero_dummy;
    if (row_idx >= __nrows)
        error_handler("Row index  out of bound");
    if (col_idx >= __ncols)
        error_handler("Col index  out of bound");
    if (!__populated[2])
        warn_handler("Accessing Uninitialized Data");
    for (size_t data_idx = __row_pointers[row_idx];
         data_idx < __row_pointers[row_idx + 1]; data_idx++)
        if (__col_indexes[data_idx] == col_idx) {
            data_ptr = &__data[data_idx];
            break;
        }
    return *data_ptr;
}
const double &csr_matrix::operator[](size_t idx) const
{
    if (idx >= __n_non_zero)
        error_handler("Out of bound access in csr_matrix::operator[]");
    if (!__populated[2])
        warn_handler("Accessing Uninitialized Data");
    return __data[idx];
}
