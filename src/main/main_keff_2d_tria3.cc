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
  size_t NX = 100, NY = 100;
  Mesh<2> mesh = mesh_create_structured_2d_tria3(NX, NY, 1.0, 1.0);

  size_t nnodes = mesh.nodes.size();
  Ellpack A(nnodes, nnodes, 7);
  Ellpack B(nnodes, nnodes, 7);

  assemblyA(A, mesh);
  assemblyB(B, mesh);

  std::vector<double> phi(nnodes, 1.0);
  for (int i = 0; i < NY; i++) {
    phi[i] = 0.0;
    A.deleteRow(i);
    A.insert(0, i, 1.0);
    B.deleteRow(i);
    B.insert(0, i, 1.0);
  }

  for (int i = 0; i < NX; i++) {
    phi[i * NY] = 0.0;
    A.deleteRow(i * NY);
    A.insert(0, i * NY, 1.0);
    B.deleteRow(i * NY);
    B.insert(0, i * NY, 1.0);
  }

  for (int i = 0; i < NY; i++) {
    phi[nnodes - NY + i] = 0.0;
    A.deleteRow(nnodes - NY + i);
    A.insert(nnodes - NY + i, nnodes - NY + 1, 1.0);
    B.deleteRow(nnodes - NY + i);
    B.insert(nnodes - NY + i, nnodes - NY + 1, 1.0);
  }

  for (int i = 0; i < NX; i++) {
    phi[i * NY + (NY - 1)] = 0.0;
    A.deleteRow(i * NY + (NY - 1));
    A.insert(0, i * NY + (NY - 1), 1.0);
    B.deleteRow(i * NY + (NY - 1));
    B.insert(0, i * NY + (NY - 1), 1.0);
  }

  double keff = solver_keff(phi, A, B);
  std::cout << "keff: " << keff << std::endl;
  return 0;
}
