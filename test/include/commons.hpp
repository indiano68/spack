#ifndef COMMONS_HPP
#define COMMONS_HPP
#include <fstream>
#include <limits>
#include <matrix.h>
#include <sstream>
#include <vector>

void vector_to_file(const std::vector<double> &x, std::string filename)
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

csr_matrix_sym load_csr_matrix_sym(const std::string &filename)
{
    std::ifstream file(filename);
    std::vector<std::vector<double>> dense;
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<double> row;
        double val;
        while (ss >> val) {
            row.push_back(val);
        }
        dense.push_back(row);
    }

    size_t n = dense.size();
    std::vector<size_t> row_ptr(n, 0);
    std::vector<size_t> col_idx;
    std::vector<double> vals;

    for (size_t i = 0; i < n; i++) {
        size_t added = 1;
        for (size_t j = 0; j < i; j++) {
            double v = dense[i][j];
            if (v != 0.0) {
                col_idx.push_back(j);
                vals.push_back(v);
                added++;
            }
        }
        col_idx.push_back(i);
        vals.push_back(dense[i][i]);
        if (i != n - 1)
            row_ptr[i + 1] = row_ptr[i] + added;
    }
    return csr_matrix_sym(n, row_ptr, col_idx, vals);
    ;
}

void print_etree_node_ascii(const std::vector<std::vector<size_t>> &children,
                            size_t node, const std::string &prefix, bool isLast)
{
    std::cout << prefix;
    if (isLast)
        std::cout << "\\-- ";
    else
        std::cout << "|-- ";

    std::cout << node << "\n";
    std::string childPrefix = prefix;
    if (isLast)
        childPrefix += "    ";
    else
        childPrefix += "|   ";

    const auto &childList = children[node];
    for (size_t i = 0; i < childList.size(); ++i) {
        bool childIsLast = (i == childList.size() - 1);
        print_etree_node_ascii(children, childList[i], childPrefix,
                               childIsLast);
    }
}

void print_etree_ascii(const std::vector<uint64_t> &etree)
{
    size_t n = etree.size();
    std::vector<std::vector<size_t>> children(n);
    std::vector<size_t> roots;

    for (size_t i = 0; i < n; ++i) {
        if (etree[i] == static_cast<uint64_t>(-1)) {
            roots.push_back(i);
        } else {
            children[etree[i]].push_back(i);
        }
    }
    std::cout << "Elimination Tree:\n";
    for (size_t i = 0; i < roots.size(); ++i) {
        std::cout << roots[i] << "\n";
        const auto &childList = children[roots[i]];
        for (size_t j = 0; j < childList.size(); ++j) {
            bool isLast = (j == childList.size() - 1);
            print_etree_node_ascii(children, childList[j], "", isLast);
        }
        if (i != roots.size() - 1) {
            std::cout << "\n";
        }
    }
}
#endif