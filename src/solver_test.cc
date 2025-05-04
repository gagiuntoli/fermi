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

#include <gtest/gtest.h>

#include "assembly.h"
#include "mesh.h"

TEST(Solver, quad4_tria3) {
  size_t NX = 100, NY = 100;
  Mesh meshQuad4 = Mesh::create2DlinearQuad4(NX, NY, 1.0, 1.0);
  Mesh meshTria3 = Mesh::create2DlinearTria3(NX, NY, 1.0, 1.0);
  double keffQuad4, keffTria3;

  {
    Ellpack A(NX * NY, NX * NY, 9);
    Ellpack B(NX * NY, NX * NY, 9);

    assemblyA(A, meshQuad4);
    assemblyB(B, meshQuad4);

    std::vector<double> phi(NX * NY, 1.0);

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

    keffQuad4 = solver_keff(phi, A, B);
  }

  {
    Ellpack A(NX * NY, NX * NY, 7);
    Ellpack B(NX * NY, NX * NY, 7);

    assemblyA(A, meshTria3);
    assemblyB(B, meshTria3);

    std::vector<double> phi(NX * NY, 1.0);

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

    keffTria3 = solver_keff(phi, A, B);
  }

  EXPECT_TRUE(std::abs(keffTria3 - keffQuad4) < 1.0e-4);
}
