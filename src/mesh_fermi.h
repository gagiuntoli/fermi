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

#ifndef MESH_FERMI_H
#define MESH_FERMI_H

#include "element.h"
#include "mesh.h"

inline Mesh mesh_create_structured_1d(size_t npoints, double length) {
  std::vector<Node> nodes;
  double h = length / (npoints - 1);
  for (size_t i = 0; i < npoints; i++) {
    nodes.push_back({i * h, 0, 0});
  }

  std::vector<std::shared_ptr<ElementBase>> elements;
  for (size_t i = 0; i < npoints - 1; i++) {
    std::vector<size_t> nodeIndexes = {i, i + 1};
    std::vector<Node> elemNodes;
    for (const size_t& index : nodeIndexes) {
      elemNodes.push_back(nodes[index]);
    }
    elements.push_back(std::make_shared<ElementSegment2>(elemNodes, nodeIndexes, 0.03, 0.04, 1.0, 1.2));
  }

  return Mesh{.nodes = std::move(nodes), .elements = std::move(elements)};
}

inline Mesh mesh_create_structured_2d_tria3(size_t nx, size_t ny, double lx, double ly) {
  std::vector<Node> nodes;
  double hx = lx / (nx - 1);
  double hy = ly / (ny - 1);
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      nodes.push_back({i * hx, j * hy, 0});
    }
  }

  std::vector<std::shared_ptr<ElementBase>> elements;
  for (size_t i = 0; i < (nx - 1); i++) {
    for (size_t j = 0; j < (ny - 1); j++) {
      std::vector<size_t> nodeIndexes = {i * ny + j, (i + 1) * ny + j, (i + 1) * ny + j + 1};
      std::vector<Node> elemNodes;
      for (const size_t& index : nodeIndexes) {
        elemNodes.push_back(nodes[index]);
      }
      elements.push_back(std::make_shared<Tria3>(elemNodes, nodeIndexes, 0.03, 0.04, 1.0, 1.2));

      nodeIndexes = {i * ny + j, i * ny + j + 1, (i + 1) * ny + j + 1};
      elemNodes.clear();
      for (const size_t& index : nodeIndexes) {
        elemNodes.push_back(nodes[index]);
      }
      elements.push_back(std::make_shared<Tria3>(elemNodes, nodeIndexes, 0.03, 0.04, 1.0, 1.2));
    }
  }

  return Mesh{.nodes = std::move(nodes), .elements = std::move(elements)};
}

inline Mesh mesh_create_structured_2d_quad4(size_t nx, size_t ny, double lx, double ly) {
  std::vector<Node> nodes;
  double hx = lx / (nx - 1);
  double hy = ly / (ny - 1);
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      nodes.push_back({i * hx, j * hy, 0});
    }
  }

  std::vector<std::shared_ptr<ElementBase>> elements;
  for (size_t i = 0; i < (nx - 1); i++) {
    for (size_t j = 0; j < (ny - 1); j++) {
      std::vector<size_t> nodeIndexes = {i * ny + j, i * ny + j + 1, (i + 1) * ny + j, (i + 1) * ny + j + 1};
      std::vector<Node> elemNodes;
      for (const size_t& index : nodeIndexes) {
        elemNodes.push_back(nodes[index]);
      }
      elements.push_back(std::make_shared<Quad4>(elemNodes, nodeIndexes, 0.03, 0.04, 1.0, 1.2));
    }
  }

  return Mesh{.nodes = std::move(nodes), .elements = std::move(elements)};
}

inline Mesh mesh_create_structured_3d_hexa8(size_t nx, size_t ny, size_t nz, double lx, double ly, double lz) {
  std::vector<Node> nodes;
  double hx = lx / (nx - 1);
  double hy = ly / (ny - 1);
  double hz = lz / (nz - 1);
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      for (size_t k = 0; k < nz; k++) {
        nodes.push_back({i * hx, j * hy, k * hz});
      }
    }
  }

  std::vector<std::shared_ptr<ElementBase>> elements;
  for (size_t i = 0; i < (nx - 1); i++) {
    for (size_t j = 0; j < (ny - 1); j++) {
      for (size_t k = 0; k < (nz - 1); k++) {
        std::vector<size_t> nodeIndexes = {i * ny + j + k * (nx * ny),
                                           i * ny + j + 1 + k * (nx * ny),
                                           (i + 1) * ny + j + k * (nx * ny),
                                           (i + 1) * ny + j + 1 + k * (nx * ny),  //
                                           i * ny + j + (k + 1) * (nx * ny),
                                           i * ny + j + 1 + (k + 1) * (nx * ny),
                                           (i + 1) * ny + j + (k + 1) * (nx * ny),
                                           (i + 1) * ny + j + 1 + (k + 1) * (nx * ny)};
        std::vector<Node> elemNodes;
        for (const size_t& index : nodeIndexes) {
          elemNodes.push_back(nodes[index]);
        }
        elements.push_back(std::make_shared<Hexa8>(elemNodes, nodeIndexes, 0.03, 0.04, 1.0, 1.2));
      }
    }
  }

  return Mesh{.nodes = std::move(nodes), .elements = std::move(elements)};
}

#endif
