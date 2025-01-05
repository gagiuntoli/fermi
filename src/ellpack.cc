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

#include "ellpack.h"

#include <cmath>
#include <cstdlib>

#include "algebra.h"

int Ellpack::mvp(std::vector<double> &y, std::vector<double> x) const {
  for (size_t row = 0; row < nrows; row++) {
    double tmp = 0;
    for (size_t j = 0; j < non_zeros_per_row; j++) {
      int index = cols[row * non_zeros_per_row + j];
      if (index < 0) break;

      tmp += vals[row * non_zeros_per_row + j] * x[index];
    }
    y[row] = tmp;
  }
  return 0;
}

int Ellpack::solve_cg(std::vector<double> &x, std::vector<double> b) const {
  int max_iters = 10000000;
  int n = nrows;
  std::vector<double> r(n);
  std::vector<double> p(n);
  std::vector<double> Ap(n);

  std::fill(x.begin(), x.end(), 0.0);

  mvp(r, x);
  for (int i = 0; i < n; i++) {
    r[i] = b[i] - r[i];
    p[i] = r[i];
  }

  double residual = norm(r, n);
  int iters = 0;
  while (residual / norm(b, n) > 1.0e-3 && iters < max_iters) {
    mvp(Ap, p);

    double rr_old = dot(r, r, n);
    double alpha = rr_old / dot(p, Ap, n);

    for (int i = 0; i < n; i++) {
      x[i] += alpha * p[i];
      r[i] -= alpha * Ap[i];
    }

    double beta = dot(r, r, n) / rr_old;

    for (int i = 0; i < n; i++) {
      p[i] = r[i] + beta * p[i];
    }

    residual = norm(r, n);
    iters++;
  }
  return 0;
}

size_t Ellpack::getIndex(size_t rowTarget, size_t colTarget) const {
  size_t index = rowTarget * non_zeros_per_row;
  for (; index < (rowTarget + 1) * non_zeros_per_row; index++) {
    if (cols[index] == -1 || cols[index] == colTarget) {
      return index;
    }
  }
  return index;
}

int Ellpack::deleteRow(size_t row) {
  size_t start = row * non_zeros_per_row;
  size_t end = (row + 1) * non_zeros_per_row;
  std::fill(cols.begin() + start, cols.begin() + end, -1);
  std::fill(vals.begin() + start, vals.begin() + end, 0.0);
  return 0;
}

int Ellpack::insert(size_t row, size_t col, double value) {
  int index = getIndex(row, col);
  if (index < (row + 1) * non_zeros_per_row) {  // entry exists or it is empty
    cols[index] = col;
    vals[index] = value;
    return 0;
  }
  return 1;
}

int Ellpack::add(size_t row, size_t col, double value) {
  int index = getIndex(row, col);
  if (index < (row + 1) * non_zeros_per_row) {  // entry exists or it is empty
    cols[index] = col;
    vals[index] += value;
    return 0;
  }
  return 1;
}

bool Ellpack::get(double &value, size_t row, size_t col) const {
  int index = getIndex(row, col);
  if (index < (row + 1) * non_zeros_per_row && cols[index] == col) {  // value exists
    value = vals[index];
    return true;
  }
  return false;
}

std::string Ellpack::toString() const {
  std::ostringstream oss;
  oss << "Row: (col, val)" << std::endl;
  for (int i = 0; i < nrows; i++) {
    oss << i << " : ";
    for (int j = 0; j < non_zeros_per_row; j++) {
      int col = cols[i * non_zeros_per_row + j];
      if (col < 0) break;
      oss << "(" << col << "," << std::scientific << std::setprecision(4) << vals[i * non_zeros_per_row + j] << "), ";
    }
    oss << std::endl;
  }
  return oss.str();
}
