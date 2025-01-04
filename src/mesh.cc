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

#include "mesh.h"

#include <sstream>

const std::string ElementBase::toString() const {
  std::ostringstream oss;
  for (const auto &node : nodes) {
    oss << node << ",";
  }
  auto result = oss.str();
  return std::string(result.begin(), result.end() - 1);
}

const std::string Mesh::toString() const {
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