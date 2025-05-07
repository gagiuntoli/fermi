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

#include "algebra.h"
#include "mesh.h"

template <size_t DIM>
struct ElementDiffusion : public ElementBase {
  double xs_a;
  double xs_f;
  double nu;
  double d;

  ElementDiffusion(std::shared_ptr<Shape> shape_, std::vector<Node> nodes_, std::vector<size_t> nodeIndexes_,
                   double xs_a_, double xs_f_, double nu_, double d_)
      : ElementBase(shape_, nodes_, nodeIndexes_), xs_a(xs_a_), xs_f(xs_f_), nu(nu_), d(d_) {}

  int computeInverseJacobian(Matrix<DIM>& ijac, double& det, size_t gp) const {
    size_t num_nodes = nodes.size();
    Matrix<DIM> jac;
    for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < DIM; j++) {
        jac.data[i][j] = 0.0;
        for (int n = 0; n < num_nodes; n++) {
          jac.data[i][j] += shape->dsh[n][i][gp] * nodes[n].getCoor(j);
        }
      }
    }
    jac.inverse(ijac, det);
    if (det < 0.0) det *= -1;
    return 0;
  }

  std::vector<double> computeAe() const override {
    const size_t n = nodes.size();
    std::vector<double> Ae(n * n, 0.0);
    auto sh = shape->sh;
    auto dsh = shape->dsh;
    auto wgp = shape->wgs;
    for (size_t gp = 0; gp < wgp.size(); gp++) {
      Matrix<DIM> ijac;
      double det;
      computeInverseJacobian(ijac, det, gp);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          std::array<double, DIM> dsh_gp_i, dsh_gp_j;
          for (int d = 0; d < DIM; d++) {
            dsh_gp_i[d] = dsh[i][d][gp];
            dsh_gp_j[d] = dsh[j][d][gp];
          }
          std::array<double, DIM> dsh_gp_it = ijac.mvp(dsh_gp_i);
          std::array<double, DIM> dsh_gp_jt = ijac.mvp(dsh_gp_j);
          double dotProd = dot<DIM>(dsh_gp_it, dsh_gp_jt);
          Ae[n * i + j] += (+d * dotProd + xs_a * sh[i][gp] * sh[j][gp]) * wgp[gp] * det;
        }
      }
    }
    return Ae;
  }

  std::vector<double> computeBe() const override {
    size_t n = nodes.size();
    std::vector<double> Be(n * n, 0.0);
    auto sh = shape->sh;
    auto wgp = shape->wgs;
    for (size_t gp = 0; gp < wgp.size(); gp++) {
      double det;
      Matrix<DIM> ijac;
      computeInverseJacobian(ijac, det, gp);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          Be[n * i + j] += nu * xs_f * sh[i][gp] * sh[j][gp] * wgp[gp] * det;
        }
      }
    }
    return Be;
  }
};

#endif
