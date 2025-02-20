#ifndef COMMONS_HPP
#define COMMONS_HPP
#include <vector>
#include <fstream>
#include <matrix.h>
#include <limits>

void vector_to_file (const std::vector<double> &x, std::string filename)
{
    std::ofstream file_x(filename);
    file_x.precision(std::numeric_limits<double>::max_digits10);
    for (size_t i = 0; i < x.size(); i++) {
        file_x << x[i] << std::endl;
    }
    file_x.close();
}

void matrix_to_file(const csr_matrix &A, std::string filename)
{
    std::ofstream file_A(filename);
    file_A.precision(std::numeric_limits<double>::max_digits10);
    for (size_t i = 0; i < A.nrows(); i++) {
        for (size_t j = 0; j < A.ncols(); j++) {
            file_A << A(i, j);
            if (j < A.ncols() - 1) {
                file_A << " ";
            }
        }
        file_A << std::endl;
    }
    file_A.close();
}

#include <random>

std::vector<double> create_random_vector(size_t n)
{
    std::vector<double> vec(n);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (size_t i = 0; i < n; i++) {
        vec[i] = dis(gen);
    }
    return vec;
}

std::vector<double> create_seq_vector(size_t n)
{
    std::vector<double> vec(n);
    for (size_t i = 0; i < n; i++) {
        vec[i] = (double)i;
    }
    return vec;
}
#endif