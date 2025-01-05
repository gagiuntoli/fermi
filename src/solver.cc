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
#include <iostream>

double solver_keff(std::vector<double> &phi, const Ellpack &A, const Ellpack &B) {
  const size_t MAX_ITERS = 15;
  size_t n = A.nrows;
  std::vector<double> source(n);
  std::vector<double> source_new(n);

  double keff = 1.0;
  size_t iters = 0;
  while (true) {
    for (auto p : phi) {
      std::cout << p << std::endl;
    }
    std::cout << std::endl;
    ellpack_mvp(source, B, phi);
    double norm_source = 0;
    for (size_t i = 0; i < n; i++) {
      norm_source += source[i];
      source[i] /= keff;
    }

    ellpack_solve_cg(phi, A, source);

    ellpack_mvp(source_new, B, phi);

    double norm_source_new = 0;
    for (size_t i = 0; i < n; i++) {
      norm_source_new += source_new[i];
    }
    double keff_new = keff * norm_source_new / norm_source;
    std::cout << "i:" << iters << " keff: " << keff << std::endl;

    if (iters > MAX_ITERS) return keff_new;
    // if (iters > MAX_ITERS || (std::abs(keff_new - keff) / keff) < 0.000001) return keff_new;
    keff = keff_new;

    double power = norm(phi, n);
    for (size_t i = 0; i < n; i++) {
      phi[i] /= power;
    }

    iters++;
  }
}
