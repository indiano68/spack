#ifndef _SPACK_LINALG_H_
#define _SPACK_LINALG_H_

#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

/***************** Interfaces *****************/
class matrix_interface
{
  public:
    virtual const double &operator()(size_t row_idx, size_t col_idx) const = 0;
    virtual const double &operator[](size_t idx) const = 0;
    virtual void set_data(const std::vector<double> &data) = 0;
    virtual double *get_data() = 0;
    virtual const double *get_data() const = 0;
    size_t size() const;
    size_t nrows() const;
    size_t ncols() const;

  protected:
    std::size_t __nrows;
    std::size_t __ncols;
    std::unique_ptr<double[]> __data;
    virtual void __allocate() = 0;
    matrix_interface(size_t row_number, size_t col_number) noexcept
        : __nrows(row_number), __ncols(col_number)
    {
    }

  private:
    matrix_interface();
};

class sparse_matrix_interface : public matrix_interface
{
  public:
    size_t r_size();

  protected:
    size_t __n_non_zero;
    sparse_matrix_interface(size_t row_number, size_t col_number,
                            size_t non_zero_number) noexcept;
    void __allocate() override;
};

/***************** DENSE FORMAT  *****************/
class dense_matrix : public matrix_interface
{
  private:
    inline void __allocate() override;

  public:
    dense_matrix(size_t row_number, size_t col_number) noexcept;
    double &operator()(size_t row_idx, size_t col_idx);
    double &operator[](size_t idx);
    const double &operator()(size_t row_idx, size_t col_idx) const override;
    const double &operator[](size_t idx) const override;
    void set_data(const std::vector<double> &data) override;
    double *get_data() override;
};

/***************** SPARSE FORMAT  *****************/
class coo_matrix : public sparse_matrix_interface
{
  private:
    std::unique_ptr<std::size_t[]> __column_indexes;
    std::unique_ptr<std::size_t[]> __row_indexes;
    bool __populated[4];

  public:
    coo_matrix(std::size_t n_rows, std::size_t n_cols,
               size_t n_non_zero) noexcept;

    coo_matrix(size_t n_rows, size_t n_cols, std::vector<size_t> r_indexes,
               std::vector<size_t> c_indexes,
               std::vector<double> data) noexcept;

    void set_indexes(const std::vector<size_t> &row_indexes,
                     const std::vector<size_t> &col_indexes);
    size_t *get_column_indexes();
    size_t *get_row_indexes();

    void set_data(const std::vector<double> &data) override;
    double *get_data() override;

    const double &operator()(size_t row_idx, size_t col_idx) const override;
    const double &operator[](size_t idx) const override;

    bool is_populated() const;
};

class csr_matrix : public sparse_matrix_interface
{
  private:
    std::unique_ptr<size_t[]> __row_pointers;
    std::unique_ptr<size_t[]> __col_indexes;
    bool __populated[4];

  public:
    csr_matrix(std::size_t n_rows, std::size_t n_cols,
               size_t n_non_zero) noexcept;

    csr_matrix(size_t n_rows, size_t n_cols, std::vector<size_t> r_pointers,
               std::vector<size_t> c_indexes,
               std::vector<double> data) noexcept;

    void set_column_indexes(const std::vector<size_t> &col_indexes);
    size_t *get_column_indexes();
    const size_t *get_column_indexes() const;
    void set_row_pointers(const std::vector<size_t> &row_pointers);
    size_t *get_row_pointers();
    const size_t *get_row_pointers() const;

    void set_data(const std::vector<double> &data) override;
    double *get_data() override;
    const double *get_data() const override;

    const double &operator()(size_t row_idx, size_t col_idx) const override;
    const double &operator[](size_t idx) const override;
    bool is_populated() const;
};

class csr_matrix_sym : public sparse_matrix_interface
{
  private:
    std::unique_ptr<size_t[]> __row_pointers;
    std::unique_ptr<size_t[]> __col_indexes;
    bool __populated[4];

  public:
    csr_matrix_sym(std::size_t n, size_t n_non_zero) noexcept;

    csr_matrix_sym(size_t n, std::vector<size_t> r_pointers,
                   std::vector<size_t> c_indexes,
                   std::vector<double> data) noexcept;

    void set_column_indexes(const std::vector<size_t> &col_indexes);
    size_t *get_column_indexes();
    const size_t *get_column_indexes() const;
    void set_row_pointers(const std::vector<size_t> &row_pointers);
    size_t *get_row_pointers();
    const size_t *get_row_pointers() const;

    void set_data(const std::vector<double> &data) override;
    double *get_data() override;
    const double *get_data() const override;

    const double &operator()(size_t row_idx, size_t col_idx) const override;
    const double &operator[](size_t idx) const override;
    bool is_populated() const;
};
#endif