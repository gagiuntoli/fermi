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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

int ellpack_mvp(std::vector<double> y, Ellpack matrix, std::vector<double> x) {
  for (size_t row = 0; row < matrix.nrows; row++) {
    double tmp = 0;
    for (size_t j = 0; j < matrix.non_zeros_per_row; j++) {
      size_t index = matrix.cols[row * matrix.non_zeros_per_row + j];

      tmp += matrix.vals[row * matrix.non_zeros_per_row + j] * x[index];
    }
    y[row] = tmp;
  }
  return 0;
}

int ellpack_solve_cg(std::vector<double> x, Ellpack matrix, std::vector<double> b) {
  int max_iters = 10000000;
  int n = matrix.nrows;
  std::vector<double> r(n);
  std::vector<double> p(n);
  std::vector<double> Ap(n);

  std::fill(x.begin(), x.end(), 0.0);

  ellpack_mvp(r, matrix, x);
  for (int i = 0; i < n; i++) {
    r[i] = b[i] - r[i];
    p[i] = r[i];
  }

  double residual = norm(r, n);
  int iters = 0;
  while (residual / norm(b, n) > 1.0e-3 && iters < max_iters) {
    ellpack_mvp(Ap, matrix, p);

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

    // printf("iter: %d residual: %lf alpha: %lf beta: %lf\n", iters, residual, alpha, beta);
    residual = norm(r, n);
    iters++;
  }
  return 0;
}

double dot(std::vector<double> y, std::vector<double> x, size_t n) {
  double result = 0.0;
  for (size_t i = 0; i < n; i++) {
    result += x[i] * y[i];
  }
  return result;
}

double norm(std::vector<double> x, size_t n) { return sqrt(dot(x, x, n)); }
