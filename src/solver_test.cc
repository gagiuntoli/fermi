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
#include "mesh_fermi.h"

TEST(Solver, quad4_tria3) {
  size_t NX = 100, NY = 100;
  Mesh<2> meshQuad4 = mesh_create_structured_2d_quad4(NX, NY, 1.0, 1.0);
  Mesh<2> meshTria3 = mesh_create_structured_2d_tria3(NX, NY, 1.0, 1.0);

  size_t nnodesQuad4 = meshQuad4.nodes.size();
  Ellpack AQuad4(nnodesQuad4, nnodesQuad4, 9);
  Ellpack BQuad4(nnodesQuad4, nnodesQuad4, 9);

  assemblyA(AQuad4, meshQuad4);
  assemblyB(BQuad4, meshQuad4);

  std::vector<double> phi(nnodesQuad4, 1.0);
  for (int i = 0; i < NY; i++) {
    phi[i] = 0.0;
    AQuad4.deleteRow(i);
    AQuad4.insert(0, i, 1.0);
    BQuad4.deleteRow(i);
    BQuad4.insert(0, i, 1.0);
  }

  for (int i = 0; i < NX; i++) {
    phi[i * NY] = 0.0;
    AQuad4.deleteRow(i * NY);
    AQuad4.insert(0, i * NY, 1.0);
    BQuad4.deleteRow(i * NY);
    BQuad4.insert(0, i * NY, 1.0);
  }

  for (int i = 0; i < NY; i++) {
    phi[nnodesQuad4 - NY + i] = 0.0;
    AQuad4.deleteRow(nnodesQuad4 - NY + i);
    AQuad4.insert(nnodesQuad4 - NY + i, nnodesQuad4 - NY + 1, 1.0);
    BQuad4.deleteRow(nnodesQuad4 - NY + i);
    BQuad4.insert(nnodesQuad4 - NY + i, nnodesQuad4 - NY + 1, 1.0);
  }

  for (int i = 0; i < NX; i++) {
    phi[i * NY + (NY - 1)] = 0.0;
    AQuad4.deleteRow(i * NY + (NY - 1));
    AQuad4.insert(0, i * NY + (NY - 1), 1.0);
    BQuad4.deleteRow(i * NY + (NY - 1));
    BQuad4.insert(0, i * NY + (NY - 1), 1.0);
  }

  double keffQuad4 = solver_keff(phi, AQuad4, BQuad4);

  size_t nnodesTria3 = meshTria3.nodes.size();
  Ellpack ATria3(nnodesTria3, nnodesTria3, 9);
  Ellpack BTria3(nnodesTria3, nnodesTria3, 9);

  assemblyA(ATria3, meshTria3);
  assemblyB(BTria3, meshTria3);

  std::vector<double> phiTria3(nnodesTria3, 1.0);
  for (int i = 0; i < NY; i++) {
    phiTria3[i] = 0.0;
    ATria3.deleteRow(i);
    ATria3.insert(0, i, 1.0);
    BTria3.deleteRow(i);
    BTria3.insert(0, i, 1.0);
  }

  for (int i = 0; i < NX; i++) {
    phiTria3[i * NY] = 0.0;
    ATria3.deleteRow(i * NY);
    ATria3.insert(0, i * NY, 1.0);
    BTria3.deleteRow(i * NY);
    BTria3.insert(0, i * NY, 1.0);
  }

  for (int i = 0; i < NY; i++) {
    phiTria3[nnodesTria3 - NY + i] = 0.0;
    ATria3.deleteRow(nnodesTria3 - NY + i);
    ATria3.insert(nnodesTria3 - NY + i, nnodesTria3 - NY + 1, 1.0);
    BTria3.deleteRow(nnodesTria3 - NY + i);
    BTria3.insert(nnodesTria3 - NY + i, nnodesTria3 - NY + 1, 1.0);
  }

  for (int i = 0; i < NX; i++) {
    phiTria3[i * NY + (NY - 1)] = 0.0;
    ATria3.deleteRow(i * NY + (NY - 1));
    ATria3.insert(0, i * NY + (NY - 1), 1.0);
    BTria3.deleteRow(i * NY + (NY - 1));
    BTria3.insert(0, i * NY + (NY - 1), 1.0);
  }

  double keffTria3 = solver_keff(phiTria3, ATria3, BTria3);

  EXPECT_TRUE(std::abs(keffTria3 - keffQuad4) < 1.0e-4);
}
