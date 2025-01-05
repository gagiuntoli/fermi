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

template <size_t DIM>
struct ElementDiffusion : public ElementBase<DIM> {
  double xs_a;
  double xs_f;
  double nu;
  double d;

  ElementDiffusion(std::vector<Node> nodes_, std::vector<size_t> nodeIndexes_, double xs_a_, double xs_f_, double nu_,
                   double d_)
      : ElementBase<DIM>(nodes_, nodeIndexes_), xs_a(xs_a_), xs_f(xs_f_), nu(nu_), d(d_) {}

  virtual std::vector<double> computeAe() const = 0;
  virtual std::vector<double> computeBe() const = 0;
  virtual int computeInverseJacobian(MatrixOperations<DIM>::Matrix &jacobian, double &det, size_t gp) const = 0;
};

struct ElementSegment2 : public ElementDiffusion<1> {
  using ElementDiffusion::ElementDiffusion;

  int computeInverseJacobian(MatrixOperations<1>::Matrix &inverseJacobian, double &det, size_t gp) const {
    Segment2 segment2;
    auto dshapes = segment2.getDShapeFunctions();
    size_t num_nodes = nodes.size();
    MatrixOperations<1>::Matrix jacobian;
    for (int i = 0; i < 1; i++) {
      for (int j = 0; j < 1; j++) {
        jacobian[i][j] = 0.0;
        for (int n = 0; n < num_nodes; n++) {
          jacobian[i][j] += dshapes[n][i][gp] * nodes[n].x;
        }
      }
    }
    MatrixOperations<1>::inverse(inverseJacobian, jacobian, det);
    return 0;
  }

  std::vector<double> computeAe() const override {
    size_t n = nodes.size();
    std::vector<double> Ae(n * n, 0.0);
    Segment2 segment2;
    double det;

    auto shapes = segment2.getShapeFunctions();
    auto dsh = segment2.getDShapeFunctions();
    auto wgp = segment2.getWeights();
    for (size_t gp = 0; gp < wgp.size(); gp++) {
      MatrixOperations<1>::Matrix ijac;
      computeInverseJacobian(ijac, det, gp);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          auto dsh_i = dsh[i][0][gp] * ijac[0][0];
          auto dsh_j = dsh[j][0][gp] * ijac[0][0];
          Ae[n * i + j] += (+d * dsh_i * dsh_j + xs_a * shapes[i][gp] * shapes[j][gp]) * wgp[gp] * det;
        }
      }
    }
    return Ae;
  }

  std::vector<double> computeBe() const override {
    size_t n = nodes.size();
    std::vector<double> Be(n * n, 0.0);
    Segment2 segment2;

    auto shapes = segment2.getShapeFunctions();
    auto wgp = segment2.getWeights();
    for (size_t gp = 0; gp < wgp.size(); gp++) {
      double det;
      MatrixOperations<1>::Matrix ijac;
      computeInverseJacobian(ijac, det, gp);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          Be[n * i + j] += nu * xs_f * shapes[i][gp] * shapes[j][gp] * wgp[gp] * det;
        }
      }
    }
    return Be;
  }
};

#endif
