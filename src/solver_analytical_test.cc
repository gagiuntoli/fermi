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

#include <gtest/gtest.h>

#include "assembly.h"
#include "mesh.h"
#include "solver.h"

TEST(AnalyticalSolution, 1d_lineal_segment) {
  const size_t NX = 200;
  Mesh mesh = Mesh::create1Dlinear(NX, 50.0);

  Ellpack A(NX, NX, 3);
  Ellpack B(NX, NX, 3);

  std::vector<double> phi(NX, 1.0);

  assemblyA(A, mesh);
  assemblyB(B, mesh);

  int index = Mesh::getIndexStructured(0, 0, 0, NX, 0, 0);
  phi[index] = 0.0;
  A.deleteRow(index);
  A.insert(index, index, 1.0);
  B.deleteRow(index);
  B.insert(index, index, 1.0);

  index = Mesh::getIndexStructured(NX - 1, 0, 0, NX, 0, 0);
  phi[index] = 0.0;
  A.deleteRow(index);
  A.insert(index, index, 1.0);
  B.deleteRow(index);
  B.insert(index, index, 1.0);

  double keff = solver_keff(phi, A, B);
  const double keffAnalytical = 1.151496323;
  std::cout << std::abs(keff - keffAnalytical) << std::endl;

  EXPECT_TRUE(std::abs(keff - keffAnalytical) < 1.0e-5);
}

TEST(AnalyticalSolution, 2d_lineal_quad) {
  size_t NX = 200, NY = 200;
  Mesh mesh = Mesh::create2DlinearQuad4(NX, NY, 50.0, 50.0);

  size_t nnodes = mesh.nodes.size();
  Ellpack A(nnodes, nnodes, 9);
  Ellpack B(nnodes, nnodes, 9);

  assemblyA(A, mesh);
  assemblyB(B, mesh);

  std::vector<double> phi(nnodes, 1.0);

  // Y = 0
  for (int i = 0; i < NX; i++) {
    size_t index = Mesh::getIndexStructured(i, 0, 0, NX, NY, 0);
    phi[index] = 0.0;
    A.deleteRow(index);
    A.insert(index, index, 1.0);
    B.deleteRow(index);
    B.insert(index, index, 1.0);
  }

  // Y = NY - 1
  for (int i = 0; i < NX; i++) {
    size_t index = Mesh::getIndexStructured(i, NY - 1, 0, NX, NY, 0);
    phi[index] = 0.0;
    A.deleteRow(index);
    A.insert(index, index, 1.0);
    B.deleteRow(index);
    B.insert(index, index, 1.0);
  }

  // X = 0
  for (int j = 0; j < NY; j++) {
    size_t index = Mesh::getIndexStructured(0, j, 0, NX, NY, 0);
    phi[index] = 0.0;
    A.deleteRow(index);
    A.insert(index, index, 1.0);
    B.deleteRow(index);
    B.insert(index, index, 1.0);
  }

  // X = NX - 1
  for (int j = 0; j < NY; j++) {
    size_t index = Mesh::getIndexStructured(NX - 1, j, 0, NX, NY, 0);
    phi[index] = 0.0;
    A.deleteRow(index);
    A.insert(index, index, 1.0);
    B.deleteRow(index);
    B.insert(index, index, 1.0);
  }

  double keff = solver_keff(phi, A, B);
  const double keffAnalytical = 1.013304171;

  EXPECT_TRUE(std::abs(keff - keffAnalytical) < 1.0e-4);
}

TEST(AnalyticalSolution, 3d_lineal_hexagon) {
  size_t NX = 30, NY = 30, NZ = 30;
  Mesh mesh = Mesh::create3DlinearHexa8(NX, NY, NZ, 50.0, 50.0, 50.0);

  size_t nnodes = mesh.nodes.size();
  Ellpack A(nnodes, nnodes, 27);
  Ellpack B(nnodes, nnodes, 27);

  assemblyA(A, mesh);
  assemblyB(B, mesh);

  std::vector<double> phi(nnodes, 1.0);

  // z = 0
  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      int index = Mesh::getIndexStructured(i, j, 0, NX, NY, NZ);
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
      int index = Mesh::getIndexStructured(i, j, NZ - 1, NX, NY, NZ);
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
      int index = Mesh::getIndexStructured(0, j, k, NX, NY, NZ);
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
      int index = Mesh::getIndexStructured(NX - 1, j, k, NX, NY, NZ);
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
      int index = Mesh::getIndexStructured(i, 0, k, NX, NY, NZ);
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
      int index = Mesh::getIndexStructured(i, NY - 1, k, NX, NY, NZ);
      phi[index] = 0.0;
      A.deleteRow(index);
      A.insert(index, index, 1.0);
      B.deleteRow(index);
      B.insert(index, index, 1.0);
    }
  }

  double keff = solver_keff(phi, A, B);
  const double keffAnalytical = 0.904727034;

  EXPECT_TRUE(std::abs(keff - keffAnalytical) < 1.0e-3);
}
