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
