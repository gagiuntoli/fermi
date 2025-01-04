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

#include "mesh_fermi.h"

#include "element.h"

Mesh mesh_create_structured_1d(size_t npoints, double length) {
  std::vector<Node> nodes;
  double h = length / (npoints - 1);
  for (size_t i = 0; i < npoints; i++) {
    nodes.push_back({i * h, 0, 0});
  }

  std::vector<std::shared_ptr<ElementBase>> elements;
  for (size_t i = 0; i < npoints - 1; i++) {
    elements.push_back(std::make_shared<ElementSegment2>(std::vector<size_t>{i, i + 1}, 1.0, 1.0, 1.0, 1.0));
  }

  return Mesh{.nodes = std::move(nodes), .elements = std::move(elements), .type = MeshType::Dim1};
}
