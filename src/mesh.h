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

struct ElementBase {
  std::vector<size_t> nodes;

  ElementBase(std::vector<size_t> nodes_) : nodes(nodes_) {}

  virtual std::vector<double> computeElementMatrix() const = 0;

  const std::string toString() const;
};

enum class MeshType { Dim1, Dim2, Dim3 };

struct Mesh {
  std::vector<Node> nodes;
  std::vector<std::shared_ptr<ElementBase>> elements;
  MeshType type;

  const std::string toString() const;
};

#endif
