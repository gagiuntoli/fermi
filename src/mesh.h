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

#ifndef MESH_H
#define MESH_H

#include <memory>
#include <string>
#include <vector>

#include "node.h"
#include "shape.h"

struct ElementBase {
  std::vector<size_t> nodeIndexes;
  std::vector<Node> nodes;
  std::shared_ptr<Shape> shape;

  ElementBase(std::shared_ptr<Shape> shape, std::vector<Node> nodes_, std::vector<size_t> nodeIndexes_)
      : shape(shape), nodes(nodes_), nodeIndexes(nodeIndexes_) {}

  virtual ~ElementBase() = default;

  const std::string toString() const {
    std::ostringstream oss;
    for (const auto &node : nodeIndexes) {
      oss << node << ",";
    }
    auto result = oss.str();
    return std::string(result.begin(), result.end() - 1);
  }

  virtual std::vector<double> computeAe() const = 0;
  virtual std::vector<double> computeBe() const = 0;
};

struct Mesh {
  std::vector<Node> nodes;
  std::vector<std::shared_ptr<ElementBase>> elements;

  static inline size_t getIndexStructured(size_t i, size_t j, size_t k, size_t NX, size_t NY, size_t NZ) {
    return i + j * NX + k * (NX * NY);
  }

  static Mesh create1Dlinear(size_t nx, double length);
  static Mesh create2DlinearQuad4(size_t nx, size_t ny, double lx, double ly);
  static Mesh create2DlinearTria3(size_t nx, size_t ny, double lx, double ly);
  static Mesh create3DlinearHexa8(size_t nx, size_t ny, size_t nz, double lx, double ly, double lz);

  const std::string toString() const {
    std::ostringstream oss;
    oss << "Nodes:" << std::endl;
    for (const Node &node : nodes) {
      oss << node.toString() << std::endl;
    }
    oss << std::endl;
    oss << "Elements:" << std::endl;
    for (const auto &element : elements) {
      oss << element->toString() << std::endl;
    }
    return oss.str();
  }
};

#endif
