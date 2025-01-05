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

#include <iostream>

#include "assembly.h"
#include "fem.h"
#include "mesh_fermi.h"
#include "solver.h"

int main(int argc, char **argv) {
  Mesh<1> mesh = mesh_create_structured_1d(10, 10.0);

  Ellpack A;
  A.nrows = mesh.nodes.size();
  A.ncols = A.nrows;
  A.non_zeros_per_row = 3;
  A.cols = std::vector<int>(A.nrows * A.non_zeros_per_row);
  A.vals = std::vector<double>(A.nrows * A.non_zeros_per_row);

  assemblyA(A, mesh);

  std::cout << mesh.toString() << std::endl;

  Segment2 segment;
  std::cout << segment.toString() << std::endl;

  // solver_keff();
  return 0;
}
