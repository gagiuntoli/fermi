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

#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <cassert>

#include "element.h"
#include "mesh.h"
#include "sparsex/src/ellpack.h"

inline int assemblyA(Ellpack &A, Mesh &mesh) {
  std::fill(A.vals.begin(), A.vals.end(), 0.0);
  for (const auto &elem_ : mesh.elements) {
    if (auto elem = std::dynamic_pointer_cast<ElementDiffusion>(elem_)) {
      auto Ae = elem->computeAe();

      size_t n = elem->nodeIndexes.size();
      assert(Ae.size() == n * n);

      for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
          size_t row = elem->nodeIndexes[i];
          size_t col = elem->nodeIndexes[j];
          A.add(row, col, Ae[i * n + j]);
        }
      }
    } else {
      return 1;
    }
  }
  return 0;
}

inline int assemblyB(Ellpack &B, Mesh &mesh) {
  std::fill(B.vals.begin(), B.vals.end(), 0.0);
  for (const auto &elem_ : mesh.elements) {
    if (auto elem = std::dynamic_pointer_cast<ElementDiffusion>(elem_)) {
      auto Be = elem->computeBe();

      size_t n = elem->nodeIndexes.size();
      assert(Be.size() == n * n);

      for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
          size_t row = elem->nodeIndexes[i];
          size_t col = elem->nodeIndexes[j];
          B.add(row, col, Be[i * n + j]);
        }
      }
    } else {
      assert(false);
      exit(1);
    }
  }
  return 0;
}

#endif
