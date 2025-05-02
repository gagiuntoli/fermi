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
  size_t NX = 5, NY = 5, NZ = 5;
  Mesh<3> mesh = mesh_create_structured_3d_hexa8(NX, NY, NZ, 1.0, 1.0, 1.0);

  size_t nnodes = mesh.nodes.size();
  Ellpack A(nnodes, nnodes, 27);
  Ellpack B(nnodes, nnodes, 27);

  assemblyA(A, mesh);
  assemblyB(B, mesh);

  std::vector<double> phi(nnodes, 1.0);
  // z = 0
  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      int index = i * NY + j;
      phi[index] = 0.0;
      A.deleteRow(index);
      A.insert(index, index, 1.0);
      B.deleteRow(index);
      B.insert(index, index, 1.0);
    }
  }

  // z = NZ - 1
  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      int index = i * NY + j + (NZ - 1) * (NX * NY);
      phi[index] = 0.0;
      A.deleteRow(index);
      A.insert(index, index, 1.0);
      B.deleteRow(index);
      B.insert(index, index, 1.0);
    }
  }

  // x = 0
  for (int k = 0; k < NZ; k++) {
    for (int j = 0; j < NY; j++) {
      int index = 0 * NY + j + k * (NX * NY);
      phi[index] = 0.0;
      A.deleteRow(index);
      A.insert(index, index, 1.0);
      B.deleteRow(index);
      B.insert(index, index, 1.0);
    }
  }

  // x = NX - 1
  for (int k = 0; k < NZ; k++) {
    for (int j = 0; j < NY; j++) {
      int index = (NX - 1) * NY + j + k * (NX * NY);
      phi[index] = 0.0;
      A.deleteRow(index);
      A.insert(index, index, 1.0);
      B.deleteRow(index);
      B.insert(index, index, 1.0);
    }
  }

  // y = 0
  for (int i = 0; i < NX; i++) {
    for (int k = 0; k < NZ; k++) {
      int index = i * NY + 0 + k * (NX * NY);
      phi[index] = 0.0;
      A.deleteRow(index);
      A.insert(index, index, 1.0);
      B.deleteRow(index);
      B.insert(index, index, 1.0);
    }
  }

  // y = NY - 1
  for (int i = 0; i < NX; i++) {
    for (int k = 0; k < NZ; k++) {
      int index = i * NY + (NY - 1) + k * (NX * NY);
      phi[index] = 0.0;
      A.deleteRow(index);
      A.insert(index, index, 1.0);
      B.deleteRow(index);
      B.insert(index, index, 1.0);
    }
  }

  double keff = solver_keff(phi, A, B);
  std::cout << "keff: " << keff << std::endl;
  return 0;
}
