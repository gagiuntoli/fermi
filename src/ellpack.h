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

#include <iomanip>
#include <sstream>
#include <string>
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

  size_t getIndex(size_t rowTarget, size_t colTarget) const;
  int deleteRow(size_t row);
  int insert(size_t row, size_t col, double value);
  int add(size_t row, size_t col, double value);
  bool get(double &value, size_t row, size_t col) const;

  int mvp(std::vector<double> &y, std::vector<double> x) const;
  int solve_cg(std::vector<double> &x, std::vector<double> b) const;

  std::string toString() const;
};

#endif
