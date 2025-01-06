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
};

struct ElementSegment2 : public ElementDiffusion<1> {
  using ElementDiffusion::ElementDiffusion;

  std::vector<double> computeAe() const override {
    size_t n = nodes.size();
    std::vector<double> Ae(n * n, 0.0);
    Segment2 segment2;
    double det;

    auto shapes = segment2.sh();
    auto dsh = segment2.dsh();
    auto wgp = segment2.weights();
    for (size_t gp = 0; gp < wgp.size(); gp++) {
      Matrix<1> ijac;
      segment2.computeInverseJacobian(ijac, det, nodes, gp);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          auto dsh_i = dsh[i][0][gp] * ijac.data[0][0];
          auto dsh_j = dsh[j][0][gp] * ijac.data[0][0];
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

    auto shapes = segment2.sh();
    auto wgp = segment2.weights();
    for (size_t gp = 0; gp < wgp.size(); gp++) {
      double det;
      Matrix<1> ijac;
      segment2.computeInverseJacobian(ijac, det, nodes, gp);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          Be[n * i + j] += nu * xs_f * shapes[i][gp] * shapes[j][gp] * wgp[gp] * det;
        }
      }
    }
    return Be;
  }
};

struct Quad2D : public ElementDiffusion<2> {
  using ElementDiffusion::ElementDiffusion;

  std::vector<double> computeAe() const override {
    size_t n = nodes.size();
    std::vector<double> Ae(n * n, 0.0);
    Quad2 quad4;
    double det;

    auto shapes = quad4.sh();
    auto dsh = quad4.dsh();
    auto wgp = quad4.weights();
    for (size_t gp = 0; gp < wgp.size(); gp++) {
      Matrix<2> ijac;
      quad4.computeInverseJacobian(ijac, det, nodes, gp);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          auto dsh_i = dsh[i][0][gp] * ijac.data[0][0];
          auto dsh_j = dsh[j][0][gp] * ijac.data[0][0];
          Ae[n * i + j] += (+d * dsh_i * dsh_j + xs_a * shapes[i][gp] * shapes[j][gp]) * wgp[gp] * det;
        }
      }
    }
    return Ae;
  }

  std::vector<double> computeBe() const override {
    size_t n = nodes.size();
    std::vector<double> Be(n * n, 0.0);
    Quad2 quad4;

    auto shapes = quad4.sh();
    auto wgp = quad4.weights();
    for (size_t gp = 0; gp < wgp.size(); gp++) {
      double det;
      Matrix<2> ijac;
      quad4.computeInverseJacobian(ijac, det, nodes, gp);
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
