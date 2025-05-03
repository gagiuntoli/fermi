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

#ifndef FERMI_H
#define FERMI_H

#include <array>
#include <cstdlib>

#include "algebra.h"
#include "node.h"

template <size_t N, size_t DIM>
class ShapeBase {
 public:
  using ShapeArray = std::array<std::array<double, N>, N>;
  using DShapeArray = std::array<std::array<std::array<double, N>, DIM>, N>;

  virtual constexpr std::array<Node, N> gaussPoints() const = 0;
  virtual constexpr std::array<double, N> weights() const = 0;
  virtual constexpr ShapeArray sh() const = 0;
  virtual constexpr DShapeArray dsh() const = 0;

  const std::string toString() const {
    std::ostringstream oss;
    oss << "Gauss Points:" << std::endl;
    for (const auto& gp : gaussPoints()) {
      oss << "  " << gp.toString() << "" << std::endl;
    }
    oss << "Weights:" << std::endl;
    for (const auto& w : weights()) {
      oss << "  " << w << "" << std::endl;
    }
    oss << "Shape Functions:" << std::endl;
    const auto& shape = sh();
    for (size_t i = 0; i < shape.size(); ++i) {
      oss << "  ";
      for (size_t j = 0; j < shape[i].size(); ++j) {
        oss << shape[i][j] << " ";
      }
      oss << "" << std::endl;
    }
    oss << "DShape Functions:" << std::endl;
    const auto& dShape = dsh();
    for (size_t i = 0; i < dShape.size(); ++i) {
      for (size_t j = 0; j < dShape[i].size(); ++j) {
        oss << "  ";
        for (size_t k = 0; k < dShape[i][j].size(); ++k) {
          oss << dShape[i][j][k] << " ";
        }
        oss << std::endl;
      }
    }
    return oss.str();
  }

  virtual ~ShapeBase() = default;

  int computeInverseJacobian(Matrix<DIM>& ijac, double& det, std::vector<Node> nodes, size_t gp) const {
    DShapeArray dsh_ = dsh();

    size_t num_nodes = nodes.size();
    Matrix<DIM> jac;
    for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < DIM; j++) {
        jac.data[i][j] = 0.0;
        for (int n = 0; n < num_nodes; n++) {
          jac.data[i][j] += dsh_[n][i][gp] * nodes[n].getCoor(j);
        }
      }
    }
    jac.inverse(ijac, det);
    if (det < 0.0) det *= -1;
    return 0;
  }
};

class ShapeSegment2 : public ShapeBase<2, 1> {
 public:
  static constexpr size_t N = 2;
  static constexpr size_t DIM = 1;

  constexpr std::array<Node, N> gaussPoints() const override {
    return {{{-0.577350269189626, 0, 0}, {+0.577350269189626, 0, 0}}};
  }

  constexpr std::array<double, N> weights() const override { return {1.0, 1.0}; }

  constexpr ShapeArray sh() const override {
    std::array<std::array<double, N>, N> sh{};
    for (size_t gp = 0; gp < N; ++gp) {
      sh[0][gp] = +0.5 * (1.0 - gaussPoints()[gp].x);
      sh[1][gp] = +0.5 * (1.0 + gaussPoints()[gp].x);
    }
    return sh;
  }

  constexpr DShapeArray dsh() const override {
    std::array<std::array<std::array<double, N>, DIM>, N> ds{};
    for (size_t gp = 0; gp < N; ++gp) {
      ds[0][0][gp] = -0.5;
      ds[1][0][gp] = +0.5;
    }
    return ds;
  }
};

class ShapeTria3 : public ShapeBase<3, 2> {
 public:
  static constexpr size_t NGP = 3;
  static constexpr size_t N = 3;
  static constexpr size_t DIM = 2;

  constexpr std::array<Node, N> gaussPoints() const override {
    return {{
        {+0.166666667, +0.166666667, 0},
        {+0.666666667, +0.166666667, 0},
        {+0.166666667, +0.666666667, 0},
    }};
  }

  constexpr std::array<double, N> weights() const override { return {0.166666667, 0.166666667, 0.166666667}; }

  constexpr ShapeArray sh() const override {
    std::array<std::array<double, N>, N> sh{};
    for (size_t gp = 0; gp < NGP; ++gp) {
      sh[0][gp] = 1.0 - gaussPoints()[gp].x - gaussPoints()[gp].y;
      sh[1][gp] = gaussPoints()[gp].x;
      sh[2][gp] = gaussPoints()[gp].y;
    }
    return sh;
  }

  constexpr DShapeArray dsh() const override {
    std::array<std::array<std::array<double, N>, DIM>, N> ds{};
    for (size_t gp = 0; gp < NGP; ++gp) {
      ds[0][0][gp] = -1.0;
      ds[0][1][gp] = -1.0;

      ds[1][0][gp] = +1.0;
      ds[1][1][gp] = +0.0;

      ds[2][0][gp] = +0.0;
      ds[2][1][gp] = +1.0;
    }
    return ds;
  }
};

class ShapeQuad4 : public ShapeBase<4, 2> {
 public:
  static constexpr size_t NGP = 4;
  static constexpr size_t N = 4;
  static constexpr size_t DIM = 2;

  constexpr std::array<Node, N> gaussPoints() const override {
    return {{{-0.577350269189626, -0.577350269189626, 0},
             {-0.577350269189626, +0.577350269189626, 0},
             {+0.577350269189626, -0.577350269189626, 0},
             {+0.577350269189626, +0.577350269189626, 0}}};
  }

  constexpr std::array<double, N> weights() const override { return {1.0, 1.0, 1.0, 1.0}; }

  constexpr ShapeArray sh() const override {
    std::array<std::array<double, N>, N> sh{};
    for (size_t gp = 0; gp < NGP; ++gp) {
      sh[0][gp] = (1.0 - gaussPoints()[gp].x) * (1.0 - gaussPoints()[gp].y) * 0.25;
      sh[1][gp] = (1.0 - gaussPoints()[gp].x) * (1.0 + gaussPoints()[gp].y) * 0.25;
      sh[2][gp] = (1.0 + gaussPoints()[gp].x) * (1.0 - gaussPoints()[gp].y) * 0.25;
      sh[3][gp] = (1.0 + gaussPoints()[gp].x) * (1.0 + gaussPoints()[gp].y) * 0.25;
    }
    return sh;
  }

  constexpr DShapeArray dsh() const override {
    std::array<std::array<std::array<double, N>, DIM>, N> ds{};
    for (size_t gp = 0; gp < NGP; ++gp) {
      ds[0][0][gp] = -1.0 * (1.0 - gaussPoints()[gp].y) * 0.25;
      ds[0][1][gp] = -1.0 * (1.0 - gaussPoints()[gp].x) * 0.25;

      ds[1][0][gp] = -1.0 * (1.0 + gaussPoints()[gp].y) * 0.25;
      ds[1][1][gp] = +1.0 * (1.0 - gaussPoints()[gp].x) * 0.25;

      ds[2][0][gp] = +1.0 * (1.0 - gaussPoints()[gp].y) * 0.25;
      ds[2][1][gp] = -1.0 * (1.0 + gaussPoints()[gp].x) * 0.25;

      ds[3][0][gp] = +1.0 * (1.0 + gaussPoints()[gp].y) * 0.25;
      ds[3][1][gp] = +1.0 * (1.0 + gaussPoints()[gp].x) * 0.25;
    }
    return ds;
  }
};

class ShapeHexa8 : public ShapeBase<8, 3> {
 public:
  static constexpr size_t NGP = 8;
  static constexpr size_t N = 8;
  static constexpr size_t DIM = 3;

  constexpr std::array<Node, N> gaussPoints() const override {
    return {{{-0.577350269189626, -0.577350269189626, -0.577350269189626},
             {-0.577350269189626, +0.577350269189626, -0.577350269189626},
             {+0.577350269189626, -0.577350269189626, -0.577350269189626},
             {+0.577350269189626, +0.577350269189626, -0.577350269189626},
             {-0.577350269189626, -0.577350269189626, +0.577350269189626},
             {-0.577350269189626, +0.577350269189626, +0.577350269189626},
             {+0.577350269189626, -0.577350269189626, +0.577350269189626},
             {+0.577350269189626, +0.577350269189626, +0.577350269189626}}};
  }

  constexpr std::array<double, N> weights() const override { return {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; }

  constexpr ShapeArray sh() const override {
    std::array<std::array<double, N>, N> sh{};
    for (size_t gp = 0; gp < NGP; ++gp) {
      sh[0][gp] = (1.0 - gaussPoints()[gp].x) * (1.0 - gaussPoints()[gp].y) * (1.0 - gaussPoints()[gp].z) * 0.125;
      sh[1][gp] = (1.0 - gaussPoints()[gp].x) * (1.0 + gaussPoints()[gp].y) * (1.0 - gaussPoints()[gp].z) * 0.125;
      sh[2][gp] = (1.0 + gaussPoints()[gp].x) * (1.0 - gaussPoints()[gp].y) * (1.0 - gaussPoints()[gp].z) * 0.125;
      sh[3][gp] = (1.0 + gaussPoints()[gp].x) * (1.0 + gaussPoints()[gp].y) * (1.0 - gaussPoints()[gp].z) * 0.125;
      sh[4][gp] = (1.0 - gaussPoints()[gp].x) * (1.0 - gaussPoints()[gp].y) * (1.0 + gaussPoints()[gp].z) * 0.125;
      sh[5][gp] = (1.0 - gaussPoints()[gp].x) * (1.0 + gaussPoints()[gp].y) * (1.0 + gaussPoints()[gp].z) * 0.125;
      sh[6][gp] = (1.0 + gaussPoints()[gp].x) * (1.0 - gaussPoints()[gp].y) * (1.0 + gaussPoints()[gp].z) * 0.125;
      sh[7][gp] = (1.0 + gaussPoints()[gp].x) * (1.0 + gaussPoints()[gp].y) * (1.0 + gaussPoints()[gp].z) * 0.125;
    }
    return sh;
  }

  constexpr DShapeArray dsh() const override {
    std::array<std::array<std::array<double, N>, DIM>, N> ds{};
    for (size_t gp = 0; gp < NGP; ++gp) {
      ds[0][0][gp] = -1.0 * (1.0 - gaussPoints()[gp].y) * (1.0 - gaussPoints()[gp].z) * 0.125;
      ds[0][1][gp] = (1.0 - gaussPoints()[gp].x) * -1.0 * (1.0 - gaussPoints()[gp].z) * 0.125;
      ds[0][2][gp] = (1.0 - gaussPoints()[gp].x) * (1.0 - gaussPoints()[gp].y) * -1.0 * 0.125;

      ds[1][0][gp] = -1.0 * (1.0 + gaussPoints()[gp].y) * (1.0 - gaussPoints()[gp].z) * 0.125;
      ds[1][1][gp] = (1.0 - gaussPoints()[gp].x) * +1.0 * (1.0 - gaussPoints()[gp].z) * 0.125;
      ds[1][2][gp] = (1.0 - gaussPoints()[gp].x) * (1.0 + gaussPoints()[gp].y) * -1.0 * 0.125;

      ds[2][0][gp] = +1.0 * (1.0 - gaussPoints()[gp].y) * (1.0 - gaussPoints()[gp].z) * 0.125;
      ds[2][1][gp] = (1.0 + gaussPoints()[gp].x) * -1.0 * (1.0 - gaussPoints()[gp].z) * 0.125;
      ds[2][2][gp] = (1.0 + gaussPoints()[gp].x) * (1.0 - gaussPoints()[gp].y) * -1.0 * 0.125;

      ds[3][0][gp] = +1.0 * (1.0 + gaussPoints()[gp].y) * (1.0 - gaussPoints()[gp].z) * 0.125;
      ds[3][1][gp] = (1.0 + gaussPoints()[gp].x) * +1.0 * (1.0 - gaussPoints()[gp].z) * 0.125;
      ds[3][2][gp] = (1.0 + gaussPoints()[gp].x) * (1.0 + gaussPoints()[gp].y) * -1.0 * 0.125;

      ds[4][0][gp] = -1.0 * (1.0 - gaussPoints()[gp].y) * (1.0 + gaussPoints()[gp].z) * 0.125;
      ds[4][1][gp] = (1.0 - gaussPoints()[gp].x) * -1.0 * (1.0 + gaussPoints()[gp].z) * 0.125;
      ds[4][2][gp] = (1.0 - gaussPoints()[gp].x) * (1.0 - gaussPoints()[gp].y) * +1.0 * 0.125;

      ds[5][0][gp] = -1.0 * (1.0 + gaussPoints()[gp].y) * (1.0 + gaussPoints()[gp].z) * 0.125;
      ds[5][1][gp] = (1.0 - gaussPoints()[gp].x) * +1.0 * (1.0 + gaussPoints()[gp].z) * 0.125;
      ds[5][2][gp] = (1.0 - gaussPoints()[gp].x) * (1.0 + gaussPoints()[gp].y) * +1.0 * 0.125;

      ds[6][0][gp] = +1.0 * (1.0 - gaussPoints()[gp].y) * (1.0 + gaussPoints()[gp].z) * 0.125;
      ds[6][1][gp] = (1.0 + gaussPoints()[gp].x) * -1.0 * (1.0 + gaussPoints()[gp].z) * 0.125;
      ds[6][2][gp] = (1.0 + gaussPoints()[gp].x) * (1.0 - gaussPoints()[gp].y) * +1.0 * 0.125;

      ds[7][0][gp] = +1.0 * (1.0 + gaussPoints()[gp].y) * (1.0 + gaussPoints()[gp].z) * 0.125;
      ds[7][1][gp] = (1.0 + gaussPoints()[gp].x) * +1.0 * (1.0 + gaussPoints()[gp].z) * 0.125;
      ds[7][2][gp] = (1.0 + gaussPoints()[gp].x) * (1.0 + gaussPoints()[gp].y) * +1.0 * 0.125;
    }
    return ds;
  }
};

#endif
