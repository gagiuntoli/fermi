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

#include "node.h"

template <size_t N, size_t DIM>
class ShapeBase {
 public:
  using ShapeArray = std::array<std::array<double, N>, N>;
  using DShapeArray = std::array<std::array<std::array<double, N>, DIM>, N>;

  virtual constexpr std::array<Node, N> getGaussPoints() const = 0;
  virtual constexpr std::array<double, N> getWeights() const = 0;
  virtual constexpr ShapeArray getShapeFunctions() const = 0;
  virtual constexpr DShapeArray getDShapeFunctions() const = 0;

  const std::string toString() const {
    std::ostringstream oss;
    oss << "Gauss Points:\n";
    for (const auto& gp : getGaussPoints()) {
      oss << "  " << gp.toString() << "\n";
    }
    oss << "Weights:\n";
    for (const auto& w : getWeights()) {
      oss << "  " << w << "\n";
    }
    oss << "Shape Functions:\n";
    const auto& shape = getShapeFunctions();
    for (size_t i = 0; i < shape.size(); ++i) {
      oss << "  ";
      for (size_t j = 0; j < shape[i].size(); ++j) {
        oss << shape[i][j] << " ";
      }
      oss << "\n";
    }
    oss << "DShape Functions:\n";
    const auto& dShape = getDShapeFunctions();
    for (size_t i = 0; i < dShape.size(); ++i) {
      for (size_t j = 0; j < dShape[i].size(); ++j) {
        oss << "  ";
        for (size_t k = 0; k < dShape[i][j].size(); ++k) {
          oss << dShape[i][j][k] << " ";
        }
        oss << "\n";
      }
    }
    return oss.str();
  }

  virtual ~ShapeBase() = default;
};

class Segment2 : public ShapeBase<2, 1> {
 public:
  static constexpr size_t N = 2;
  static constexpr size_t DIM = 1;

  constexpr std::array<Node, N> getGaussPoints() const override {
    return {{{-0.577350269189626, 0, 0}, {+0.577350269189626, 0, 0}}};
  }

  constexpr std::array<double, N> getWeights() const override { return {1.0, 1.0}; }

  constexpr ShapeArray getShapeFunctions() const override {
    std::array<std::array<double, N>, N> sh{};
    for (size_t gp = 0; gp < N; ++gp) {
      sh[0][gp] = +0.5 * (1.0 - getGaussPoints()[gp].x);
      sh[1][gp] = +0.5 * (1.0 + getGaussPoints()[gp].x);
    }
    return sh;
  }

  constexpr DShapeArray getDShapeFunctions() const override {
    std::array<std::array<std::array<double, N>, DIM>, N> ds{};
    for (size_t gp = 0; gp < N; ++gp) {
      ds[0][0][gp] = -0.5;
      ds[1][0][gp] = +0.5;
    }
    return ds;
  }
};

#endif
