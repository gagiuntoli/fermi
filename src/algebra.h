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

  Matrix() : data{} {}
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
double Matrix<1>::determinant();
template <>
int Matrix<1>::inverse(Matrix& inverse, double& det);

template <>
double Matrix<2>::determinant();
template <>
int Matrix<2>::inverse(Matrix& inverse, double& det);

template <>
double Matrix<3>::determinant();
template <>
int Matrix<3>::inverse(Matrix& inverse, double& det);

double dot(const std::vector<double>& y, const std::vector<double>& x, size_t n);
double norm(const std::vector<double>& x, size_t n);

#endif
