/*
 *  This source code is part of Fermi: a finite element code
 *  to solve the neutron diffusion problem for nuclear reactor
 *  designs.
 *
 *  Copyright (C) - 2019 - Guido Giuntoli <gagiuntoli@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef ELLPACK_H
#define ELLPACK_H

#include <vector>

struct Ellpack {
  size_t nrows;
  size_t ncols;
  size_t non_zeros_per_row;
  std::vector<int> cols;
  std::vector<double> vals;

  Ellpack() = default;

  Ellpack(size_t nrows_, size_t ncols_, size_t non_zeros_per_row_)
      : nrows(nrows_),
        ncols(ncols_),
        non_zeros_per_row(non_zeros_per_row_),
        cols(std::vector<int>(ncols * non_zeros_per_row, -1)),
        vals(std::vector<double>(ncols * non_zeros_per_row, 0.0)) {}

  size_t getIndex(size_t rowTarget, size_t colTarget) {
    size_t index = rowTarget * non_zeros_per_row;
    for (; index < (rowTarget + 1) * non_zeros_per_row; index++) {
      if (cols[index] == -1 || cols[index] == colTarget) {
        return index;
      }
    }
    return index;
  }

  int insert(size_t row, size_t col, double value) {
    int index = getIndex(row, col);
    if (index < (row + 1) * non_zeros_per_row) {  // entry exists or it is empty
      cols[index] = col;
      vals[index] = value;
      return 0;
    }
    return 1;
  }

  bool get(double &value, size_t row, size_t col) {
    int index = getIndex(row, col);
    if (index < (row + 1) * non_zeros_per_row && cols[index] == col) {  // value exists
      value = vals[index];
      return true;
    }
    return false;
  }
};

int ellpack_mvp(std::vector<double> y, Ellpack matrix, std::vector<double> x);
int ellpack_solve_cg(std::vector<double> x, Ellpack matrix, std::vector<double> b);
double dot(std::vector<double> y, std::vector<double> x, size_t n);
double norm(std::vector<double> x, size_t n);

#endif
