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

#include "algebra.h"

#include <cmath>

template <>
double Matrix<1>::determinant() {
  return data[0][0];
}

template <>
int Matrix<1>::inverse(Matrix& inverse, double& det) {
  det = determinant();
  if (std::abs(det) < 1e-10) return 1;
  inverse = {{{1.0 / det}}};
  return 0;
}

template <>
double Matrix<2>::determinant() {
  return data[0][0] * data[1][1] - data[0][1] * data[1][0];
}

template <>
int Matrix<2>::inverse(Matrix& inverse, double& det) {
  det = determinant();
  if (std::abs(det) < 1e-10) return 1;

  inverse.data[0][0] = data[1][1] / det;
  inverse.data[0][1] = -data[0][1] / det;
  inverse.data[1][0] = -data[1][0] / det;
  inverse.data[1][1] = data[0][0] / det;
  return 0;
}

template <>
double Matrix<3>::determinant() {
  return data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1]) -
         data[0][1] * (data[1][0] * data[2][2] - data[1][2] * data[2][0]) +
         data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]);
}

template <>
int Matrix<3>::inverse(Matrix& inverse, double& det) {
  det = determinant();
  if (std::abs(det) < 1e-10) return 1;

  std::array<std::array<double, 3>, 3> cofactor;
  cofactor[0][0] = data[1][1] * data[2][2] - data[2][1] * data[1][2];
  cofactor[0][1] = -(data[1][0] * data[2][2] - data[2][0] * data[1][2]);
  cofactor[0][2] = data[1][0] * data[2][1] - data[2][0] * data[1][1];
  cofactor[1][0] = -(data[0][1] * data[2][2] - data[2][1] * data[0][2]);
  cofactor[1][1] = data[0][0] * data[2][2] - data[2][0] * data[0][2];
  cofactor[1][2] = -(data[0][0] * data[2][1] - data[2][0] * data[0][1]);
  cofactor[2][0] = data[0][1] * data[1][2] - data[1][1] * data[0][2];
  cofactor[2][1] = -(data[0][0] * data[1][2] - data[1][0] * data[0][2]);
  cofactor[2][2] = data[0][0] * data[1][1] - data[1][0] * data[0][1];

  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      inverse.data[j][i] = cofactor[i][j] / det;
    }
  }
  return 0;
}

double dot(const std::vector<double>& y, const std::vector<double>& x, size_t n) {
  double result = 0.0;
  for (size_t i = 0; i < n; i++) {
    result += x[i] * y[i];
  }
  return result;
}

double norm(const std::vector<double>& x, size_t n) { return sqrt(dot(x, x, n)); }
