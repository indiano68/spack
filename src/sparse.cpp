#include <linalg.h>
#include <cstring>
#include <algorithm>
#include <string>
#include <iostream>

static void error_handler(std::string error_message)
{
    std::cerr <<"ERROR:" << error_message << std::endl;
    std::exit(EXIT_FAILURE);
}


/******************** Implement Matrix ********************/

size_t matrix_interface::size()
{
    return __ncols*__nrows;
}
size_t matrix_interface::nrows()
{
    return __nrows;
}
size_t matrix_interface::ncols()
{
    return __ncols;
}

/******************* Implement Dense Matrix ********************/

dense_matrix::dense_matrix(    
    size_t row_number, 
    size_t col_number
) noexcept: matrix_interface(row_number,col_number)
{
    __allocate();
}
void dense_matrix::__allocate()
{
    __data = std::make_unique<double[]>(__nrows*__ncols);
};
void dense_matrix::set_data(const std::vector<double> & data    )
{
    if(__ncols*__nrows!=data.size()) 
        error_handler("Data dimention incorrect in call to set_data (Dense)");
    std::memcpy(__data.get(), data.data(), data.size());
}
double * dense_matrix::get_data()
{
    return __data.get();
}
double & dense_matrix::operator()(size_t row_idx, size_t col_idx) 
{
    return __data[__ncols*row_idx+col_idx];
}
double & dense_matrix::operator[](size_t idx)
{
    return __data[idx];
}
const double & dense_matrix::operator()(size_t row_idx, size_t col_idx) const 
{
    return __data[__ncols*row_idx+col_idx];
}
const double & dense_matrix::operator[](size_t idx) const 
{
    return __data[idx];
}

/****************** Implement Sparse Matrix ******************/

sparse_matrix_interface::sparse_matrix_interface(
    size_t row_number, 
    size_t col_number, 
    size_t non_zero_number) noexcept: 
    matrix_interface(row_number,col_number), 
    __non_zero_elements(non_zero_number)
{
    __allocate();
}
size_t  sparse_matrix_interface::r_size()
{
    return __non_zero_elements;
}
void sparse_matrix_interface::__allocate()
{
    __data = std::make_unique<double[]>(__non_zero_elements);
};

/******************* Implement CRF *******************/

crf_matrix::crf_matrix(
    std::size_t row_count, 
    std::size_t col_count, 
    size_t non_zero_count) noexcept: 
    sparse_matrix_interface(row_count,col_count,non_zero_count),
    __column_indexes(std::make_unique<size_t[]>(non_zero_count)),
    __row_indexes(std::make_unique<size_t[]>(non_zero_count)),
    __populated{false}{}
crf_matrix::crf_matrix( 
    size_t n_rows,
    size_t n_cols,
    std::vector<size_t> r_indexes, 
    std::vector<size_t> c_indexes, 
    std::vector<double> data) noexcept:
    sparse_matrix_interface(n_rows,n_cols,data.size()) 
{

    if(    r_indexes.size() == c_indexes.size() 
        && c_indexes.size() == data.size())
    {
        __column_indexes = std::make_unique<size_t[]>(__non_zero_elements);
        __row_indexes = std::make_unique<size_t[]>(__non_zero_elements);
        set_indexes(r_indexes,c_indexes);
        set_data(data);
    }
    else 
        error_handler("Incorrectly sized objects passed to crf_matrix_constructor");
}

static bool __check_order(const std::vector<size_t>& row_idx, const std::vector<size_t>& col_idx) 
{
    bool result = true;
    for (std::size_t i = 1; i < row_idx.size(); ++i) {
        if (row_idx[i] < row_idx[i - 1]) {
            result= false;
            break;
        }        
        if (row_idx[i] == row_idx[i - 1] && col_idx[i] < col_idx[i - 1]) {
            result = false;
            break;
        }
    }
    return result;
}


void crf_matrix::set_indexes(
    const std::vector<size_t> & row_indexes, 
    const std::vector<size_t> & col_indexes
)
{
    if(row_indexes.empty() || col_indexes.empty())
        error_handler("Some of the indexes are empty");
    if(row_indexes.size()!= col_indexes.size()  || 
       row_indexes.size()!= __non_zero_elements ||
       col_indexes.size()!= __non_zero_elements)
        error_handler("One ar all index sizes not matching non-zero element number");
    if(!__check_order(row_indexes,col_indexes)) 
        error_handler("Row and/or column indexes are not respecting the required forma (in-order)");
    std::memcpy(__row_indexes.get(),row_indexes.data(),__non_zero_elements); 
    std::memcpy(__column_indexes.get(),col_indexes.data(),__non_zero_elements); 
    __populated[0] =true;
    __populated[1] =true;
    __populated[3] = __populated[0] && __populated[1] && __populated[2];
}
void crf_matrix::set_data(const std::vector<double> &  data)
{
    if(data.size()!=__non_zero_elements) 
        error_handler("Number of data entries does not match non-zero element number");;
    std::memcpy(__data.get(),data.data(),__non_zero_elements); 
    __populated[2] =true;
    __populated[3] = __populated[0] && __populated[1] && __populated[2];
}
size_t * crf_matrix::get_column_indexes()
{
    return __column_indexes.get();
}
size_t * crf_matrix::get_row_indexes()
{
    return __row_indexes.get();
}
double * crf_matrix::get_data()
{
    return __data.get();
}
bool crf_matrix::is_populated() const
{
    return __populated[3];
}

/******************** Implement CRS ********************/