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

#include "matrix_operations.h"
#include "node.h"

template <size_t DIM>
struct ElementBase {
  std::vector<size_t> nodeIndexes;
  std::vector<Node> nodes;

  ElementBase(std::vector<Node> nodes_, std::vector<size_t> nodeIndexes_) : nodeIndexes(nodeIndexes_), nodes(nodes_) {}

  const std::string toString() const {
    std::ostringstream oss;
    for (const auto &node : nodeIndexes) {
      oss << node << ",";
    }
    auto result = oss.str();
    return std::string(result.begin(), result.end() - 1);
  }
};

template <size_t DIM>
struct Mesh {
  std::vector<Node> nodes;
  std::vector<std::shared_ptr<ElementBase<DIM>>> elements;

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
