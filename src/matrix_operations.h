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

#ifndef MATRIX_OPERATION_H
#define MATRIX_OPERATION_H

#include <array>
#include <cstdlib>

template <size_t N>
class MatrixOperations {
 public:
  using Matrix = std::array<std::array<double, N>, N>;

  static double determinant(const Matrix& matrix);
  static int inverse(Matrix& inverse, const Matrix& matrix, double& det);
};

template <>
class MatrixOperations<1> {
 public:
  using Matrix = std::array<std::array<double, 1>, 1>;

  static double determinant(const Matrix& matrix) { return matrix[0][0]; }

  static int inverse(Matrix& inverse, const Matrix& matrix, double& det) {
    det = determinant(matrix);
    if (std::abs(det) < 1e-10) return 1;
    inverse = {{{1.0 / det}}};
    return 0;
  }
};

template <>
class MatrixOperations<2> {
 public:
  using Matrix = std::array<std::array<double, 2>, 2>;

  static double determinant(const Matrix& matrix) { return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]; }

  static int inverse(Matrix& inverse, const Matrix& matrix, double& det) {
    det = determinant(matrix);
    if (std::abs(det) < 1e-10) return 1;

    inverse = {{
        {{matrix[1][1] / det, -matrix[0][1] / det}},
        {{-matrix[1][0] / det, matrix[0][0] / det}},
    }};
    return 0;
  }
};

template <>
class MatrixOperations<3> {
 public:
  using Matrix = std::array<std::array<double, 3>, 3>;

  static double determinant(const Matrix& matrix) {
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
           matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
           matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
  }

  static int inverse(Matrix& inverse, const Matrix& matrix, double& det) {
    det = determinant(matrix);
    if (std::abs(det) < 1e-10) return 1;

    Matrix cofactor;
    cofactor[0][0] = matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2];
    cofactor[0][1] = -(matrix[1][0] * matrix[2][2] - matrix[2][0] * matrix[1][2]);
    cofactor[0][2] = matrix[1][0] * matrix[2][1] - matrix[2][0] * matrix[1][1];
    cofactor[1][0] = -(matrix[0][1] * matrix[2][2] - matrix[2][1] * matrix[0][2]);
    cofactor[1][1] = matrix[0][0] * matrix[2][2] - matrix[2][0] * matrix[0][2];
    cofactor[1][2] = -(matrix[0][0] * matrix[2][1] - matrix[2][0] * matrix[0][1]);
    cofactor[2][0] = matrix[0][1] * matrix[1][2] - matrix[1][1] * matrix[0][2];
    cofactor[2][1] = -(matrix[0][0] * matrix[1][2] - matrix[1][0] * matrix[0][2]);
    cofactor[2][2] = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];

    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        inverse[j][i] = cofactor[i][j] / det;
      }
    }
    return 0;
  }
};

#endif
