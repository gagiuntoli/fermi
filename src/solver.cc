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

#include "solver.h"

#include <cmath>

#include "algebra.h"

double solver_keff(std::vector<double> &phi, const Ellpack &A, const Ellpack &B) {
  const size_t MAX_ITERS = 15;
  size_t n = A.nrows;
  std::vector<double> source(n);
  std::vector<double> source_new(n);

  double keff = 1.0;
  size_t iters = 0;
  while (true) {
    B.mvp(source, phi);
    double norm_source = 0;
    for (size_t i = 0; i < n; i++) {
      norm_source += source[i];
      source[i] /= keff;
    }

    A.solve_cg(phi, source);

    B.mvp(source_new, phi);

    double norm_source_new = 0;
    for (size_t i = 0; i < n; i++) {
      norm_source_new += source_new[i];
    }
    double keff_new = keff * norm_source_new / norm_source;

    if (iters > MAX_ITERS) return keff_new;
    keff = keff_new;

    double power = norm(phi, n);
    for (size_t i = 0; i < n; i++) {
      phi[i] /= power;
    }

    iters++;
  }
}
