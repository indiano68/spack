#ifndef _LINALG_H_
#define _LINALG_H_

#include <cstddef>
#include <memory>
#include <iostream>
#include <vector>
class matrix_interface
{
public:
    virtual double &operator()(size_t row_idx, size_t col_idx) = 0;
    virtual double &operator[](size_t idx) = 0;
    virtual const double &operator()(size_t row_idx, size_t col_idx) const = 0;
    virtual const double &operator[](size_t idx) const = 0;
    virtual void set_data(const std::vector<double> & data) = 0;
    virtual double * get_data() = 0;
    size_t size();
    size_t nrows();
    size_t ncols();


protected:
    std::size_t __nrows;
    std::size_t __ncols;
    std::unique_ptr<double[]> __data;
    virtual void __allocate() = 0;
    matrix_interface(size_t row_number, size_t col_number) noexcept : __nrows(row_number),
                                                                      __ncols(col_number) {}

private:
    matrix_interface();
};

class dense_matrix : public matrix_interface
{
private:
    inline void __allocate() override;
public:
    dense_matrix(    
        size_t row_number, 
        size_t col_number
    ) noexcept;
    double &operator()(size_t row_idx, size_t col_idx) override;
    double &operator[](size_t idx) override;
    const double &operator()(size_t row_idx, size_t col_idx) const override;
    const double &operator[](size_t idx) const override;
    void set_data(const std::vector<double> & data) override;
    double * get_data() override;
};



class sparse_matrix_interface : public matrix_interface
{
public:
    size_t r_size();
protected:
    size_t __non_zero_elements;
    sparse_matrix_interface(
        size_t row_number,
        size_t col_number,
        size_t non_zero_number) noexcept;
    inline void __allocate() override;
};

class crf_matrix : public sparse_matrix_interface
{
private:
    std::unique_ptr<std::size_t[]> __column_indexes;
    std::unique_ptr<std::size_t[]> __row_indexes;
    bool __populated[4];

public:
    crf_matrix(
        std::size_t row_count,
        std::size_t col_count,
        size_t non_zero_count) noexcept;

    crf_matrix(
        size_t n_rows,
        size_t n_cols,
        std::vector<size_t> r_indexes,
        std::vector<size_t> c_indexes,
        std::vector<double> data) noexcept;

    void set_indexes(const std::vector<size_t> & row_indexes, const std::vector<size_t> & col_indexes);
    size_t *get_column_indexes();
    size_t *get_row_indexes();

    void set_data(const std::vector<double> & data ) override;
    double * get_data() override;

    // double &operator()(size_t row_idx, size_t col_idx);
    // double &operator[](size_t idx);
    // const double &operator()(size_t row_idx, size_t col_idx) const;
    // const double &operator[](size_t idx) const;
    bool is_populated() const;
    void print()
    {
        std::cout << __populated[0]
                  << __populated[1]
                  << __populated[2]
                  << __populated[3] << std::endl;
    }
};

// class csr_matrix: public sparse_matrix_interface
// {
//     private:
//     std::unique_ptr<std::size_t[]> __column_indexes;
//     std::unique_ptr<std::size_t[]> row_pointer;
//     csr_matrix();
//     public:
//     csr_matrix(
//         std::size_t row_count,
//         std::size_t col_count,
//         size_t non_zero_count) noexcept:
//         sparse_matrix_interface(row_count,col_count,non_zero_count),
//         __column_indexes(std::make_unique<size_t[]>(non_zero_count)),
//         row_pointer(std::make_unique<size_t[]>(row_count+1)){}
//     double & operator()(size_t row_idx, size_t col_idx){ return __data[0];}
//     double & operator[](size_t idx){ return __data[0];}
//     const double & operator()(size_t row_idx, size_t col_idx) const
//     { return __data[0];}
//     const double & operator[](size_t idx)                     const
//     { return __data[0];}
//     void print_csr()
//     {
//         std::cout
//         << __nrows
//         << " "
//         << __ncols
//         << " "
//         << __non_zero_elements << std::endl;
//     }
// };

#endif