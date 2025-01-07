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

#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <array>
#include <cstdlib>
#include <vector>

template <size_t N>
class Matrix {
 public:
  std::array<std::array<double, N>, N> data;

  Matrix(std::array<std::array<double, N>, N> data_) : data(data_) {}

  double determinant();
  int inverse(Matrix& inverse, double& det);

  std::array<double, N> mvp(const std::array<double, N>& x) const {
    std::array<double, N> y = {};
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        y[i] += data[i][j] * x[j];
      }
    }
    return y;
  }
};

template <>
class Matrix<1> {
 public:
  std::array<std::array<double, 1>, 1> data;

  double determinant() { return data[0][0]; }

  int inverse(Matrix& inverse, double& det) {
    det = determinant();
    if (std::abs(det) < 1e-10) return 1;
    inverse = {{{1.0 / det}}};
    return 0;
  }
};

template <>
class Matrix<2> {
 public:
  std::array<std::array<double, 2>, 2> data;

  double determinant() { return data[0][0] * data[1][1] - data[0][1] * data[1][0]; }

  int inverse(Matrix& inverse, double& det) {
    det = determinant();
    if (std::abs(det) < 1e-10) return 1;

    inverse.data[0][0] = data[1][1] / det;
    inverse.data[0][1] = -data[0][1] / det;
    inverse.data[1][0] = -data[1][0] / det;
    inverse.data[1][1] = data[0][0] / det;
    return 0;
  }

  std::array<double, 2> mvp(const std::array<double, 2>& x) const {
    std::array<double, 2> y = {};
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        y[i] += data[i][j] * x[j];
      }
    }
    return y;
  }
};

template <>
class Matrix<3> {
 public:
  std::array<std::array<double, 3>, 3> data;

  double determinant() {
    return data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1]) -
           data[0][1] * (data[1][0] * data[2][2] - data[1][2] * data[2][0]) +
           data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]);
  }

  int inverse(Matrix& inverse, double& det) {
    det = determinant();
    if (std::abs(det) < 1e-10) return 1;

    Matrix<3> cofactor;
    cofactor.data[0][0] = data[1][1] * data[2][2] - data[2][1] * data[1][2];
    cofactor.data[0][1] = -(data[1][0] * data[2][2] - data[2][0] * data[1][2]);
    cofactor.data[0][2] = data[1][0] * data[2][1] - data[2][0] * data[1][1];
    cofactor.data[1][0] = -(data[0][1] * data[2][2] - data[2][1] * data[0][2]);
    cofactor.data[1][1] = data[0][0] * data[2][2] - data[2][0] * data[0][2];
    cofactor.data[1][2] = -(data[0][0] * data[2][1] - data[2][0] * data[0][1]);
    cofactor.data[2][0] = data[0][1] * data[1][2] - data[1][1] * data[0][2];
    cofactor.data[2][1] = -(data[0][0] * data[1][2] - data[1][0] * data[0][2]);
    cofactor.data[2][2] = data[0][0] * data[1][1] - data[1][0] * data[0][1];

    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        inverse.data[j][i] = cofactor.data[i][j] / det;
      }
    }
    return 0;
  }
};

double dot(const std::vector<double>& y, const std::vector<double>& x, size_t n);
double norm(const std::vector<double>& x, size_t n);

#endif
