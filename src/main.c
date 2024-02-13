/*
 *  This source code is part of Fermi: a finite element code
 *  to solve the neutron diffusion problem for nuclear reactor
 *  designs.
 *
 *  Copyright (C) - 2019 - Guido Giuntoli <gagiuntoli@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "fermi.h"

MPI_Comm WORLD_Comm; // global communicator
MPI_Comm FERMI_Comm; // local  communicator

int globa_rank; // rank in WORLD_Comm
int globa_size; // size in WORLD_Comm
int local_rank; // rank in FERMI_Comm
int local_size; // size in FERMI_Comm

bool couple_fl;
coupling_t coupling;

PetscViewer kspview;
PetscViewer viewer;

int rank;
int nproc;
int nummat;

int *loc2gold, *loc2gnew;
int ndir;
int *npp;
int ntot;
int *dirIndex;
double *dirValue;
double *dirZeros;

int memory;

list_t list_nodes;
list_t list_ghost;
list_t list_elemv;
list_t list_elems;
list_t list_physe;
list_t list_mater;
list_t list_bound;
list_t list_fun1d; /* list of functions */
list_t list_ctrlr; /* list of control rods */
list_t list_outpu;
list_t list_comms;

mesh_t mesh;

int nke, nbe, nev, its;
double *Ae, *Be, *Me, *be;
Vec phi_n, phi_o, b, b_a;
Vec xlocal;
double ikeff, keff;
Mat A, B, M, K;
KSP ksp;
PC pc; /* preconditioner context */
int Istart, Iend;

char inputfile[32];
char meshfile[32];
char epartfile[32];
char npartfile[32];

calcu_t calcu;

int *idxm;
double **der, ***ode, **sh, **jac, **ijac, **coor, *wp;

int egn, pgn;
int nxs_mat; // number of xs values per material
double **phi;
double power;
double vol;

/* Precursors constants */
double *veloc, *lambda, *beta, *chi, beta_tot;
double dtn;

double **xp_segm_2;
double *wp_segm_2;
double **sh_segm_2;
double ***ds_segm_2;

double **xp_tria_3;
double *wp_tria_3;
double **sh_tria_3;
double ***ds_tria_3;

double **xp_quad_4;
double *wp_quad_4;
double **sh_quad_4;
double ***ds_quad_4;

double **xp_tetra_4;
double *wp_tetra_4;
double **sh_tetra_4;
double ***ds_tetra_4;

double **xp_prism_6;
double *wp_prism_6;
double **sh_prism_6;
double ***ds_prism_6;

double **xp_hexa_8;
double *wp_hexa_8;
double **sh_hexa_8;
double ***ds_hexa_8;

int main(int argc, char **argv) {

  int step = 0;
  char nam[32];
  node_list_t *pNod;

  if (ferinit(argc, argv)) {
    goto error;
  }

  calcu.t = calcu.t0;

  if (calcu.timedep== QS) {

    pNod = calcu.time.head;

    while (pNod) {

      dtn = ((tcontrol_t *)pNod->data)->dt;

      fersrods(calcu.t);

      ferass_ST();

      fersolv_ST();

      fer_norm();

      sprintf(nam, "steady_r%d_t%d", rank, step);
      print_vtk(nam);
      print_out(&phi_n, step);

      calcu.t = calcu.t + dtn;
      step++;
      if (calcu.t > ((tcontrol_t *)pNod->data)->tf - 1.0e-8)
        pNod = pNod->next;
    }
  }

  else if (calcu.timedep == TR) {

    //==================================================
    // Transient simulation
    //
    // 1) Calculate steady state solving Ax = (1/k) Bx
    // 2) Recv information for coupling (if it required)
    // 3) Calculates a "dt" increase in the flux solving
    //    Ax = b
    // 4) Send information for coupling (if it required)
    // 5) Repeat from "2)" up to achieving final time "tf"
    //
    PetscPrintf(FERMI_Comm, "calculating stationary state.\n");

    fersrods(calcu.t);

    ferass_ST();

    fersolv_ST();

    fer_norm();

    sprintf(nam, "steady_r%d", rank);
    print_vtk(nam);

    PetscPrintf(FERMI_Comm, "initial power:%e\n", power);
    PetscPrintf(FERMI_Comm, "time    power   its\n");

    pNod = calcu.time.head;
    while (pNod) {

      dtn = ((tcontrol_t *)pNod->data)->dt;

      fersrods(calcu.t);

      ferass_TR(step);

      fersolv_TR();

      fer_pow(&power);

      PetscPrintf(FERMI_Comm, "%lf %e %d \n", calcu.t, power, its);

      sprintf(nam, "tran_rank%d_t%d", rank, step);
      print_vtk(nam);
      print_out(&phi_n, step);

      calcu.t = calcu.t + dtn;
      step++;
      if (calcu.t > ((tcontrol_t *)pNod->data)->tf - 1.0e-8)
        pNod = pNod->next;
    }
  }

error:

  ferfini();

  return 0;
}
