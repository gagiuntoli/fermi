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

#ifndef ELEMENT_H
#define ELEMENT_H

#include "fem.h"
#include "mesh.h"

struct ElementDiffusion : public ElementBase {
  double xs_a;
  double xs_f;
  double nu;
  double d;

  ElementDiffusion(std::vector<size_t> nodes_, double xs_a_, double xs_f_, double nu_, double d_)
      : ElementBase(nodes_), xs_a(xs_a_), xs_f(xs_f_), nu(nu_), d(d_) {}
};

struct ElementSegment2 : public ElementDiffusion {
  using ElementDiffusion::ElementDiffusion;

  std::vector<double> computeElementMatrix() const override {
    size_t n = nodes.size();
    std::vector<double> matrix(n * n, 0.0);
    Segment2 segment2;

    for (const auto& gp : segment2.getGaussPoints()) {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
        }
      }
    }

    return matrix;
  }
};

#endif
