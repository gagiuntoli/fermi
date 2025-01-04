#ifndef ELLPACK_H
#define ELLPACK_H

#include <vector>

typedef struct {
    size_t nrows;
    size_t ncols;
    size_t non_zeros_per_row;
    std::vector<size_t> cols;
    std::vector<double> vals;
} ellpack_t;

int ellpack_mvp(std::vector<double> y, ellpack_t matrix, std::vector<double> x);
int ellpack_solve_cg(std::vector<double> x, ellpack_t matrix, std::vector<double> b);
double dot(std::vector<double> y, std::vector<double> x, size_t n);
double norm(std::vector<double> x, size_t n);

#endif
