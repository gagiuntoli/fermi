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
#include "mesh_fermi.h"
#include "solver.h"

int main(int argc, char **argv) {
  Mesh<1> mesh = mesh_create_structured_1d(100, 50.0);

  size_t nnodes = mesh.nodes.size();
  Ellpack A(nnodes, nnodes, 3);
  Ellpack B(nnodes, nnodes, 3);

  std::vector<double> phi(nnodes, 1.0);
  phi[0] = 0.0;
  phi[nnodes - 1] = 0.0;

  assemblyA(A, mesh);
  A.deleteRow(0);
  A.deleteRow(nnodes - 1);
  A.insert(0, 0, 1.0);
  A.insert(nnodes - 1, nnodes - 1, 1.0);

  assemblyB(B, mesh);
  B.deleteRow(0);
  B.deleteRow(nnodes - 1);
  B.insert(0, 0, 1.0);
  B.insert(nnodes - 1, nnodes - 1, 1.0);

  double keff = solver_keff(phi, A, B);
  std::cout << "keff: " << keff << std::endl;
  return 0;
}
