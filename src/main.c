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

int main(int argc, char **argv)
{

    int step = 0;
    char nam[32];
    node_list_t * pNod;

    if(ferinit(argc,argv)) {
        goto error;
    }

    calcu.t = calcu.t0;

    if(calcu.timedep == QS){

        pNod = calcu.time.head;

        while(pNod){

            dtn = ((tcontrol_t*)pNod->data)->dt;

            fersrods(calcu.t);

            ferass_ST();

            fersolv_ST();

            fer_norm();

            sprintf(nam,"steady_r%d_t%d", rank, step);
            print_vtk(nam);
            print_out(&phi_n, step);

            calcu.t=calcu.t + dtn;
            step ++;
            if(calcu.t > ((tcontrol_t*)pNod->data)->tf - 1.0e-8)
            pNod = pNod->next;

        }
    }

    else if(calcu.timedep == TR){

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
        PetscPrintf(FERMI_Comm,"calculating stationary state.\n");

        fersrods(calcu.t);

        ferass_ST();

        fersolv_ST();

        fer_norm();

        sprintf(nam,"steady_r%d",rank);
        print_vtk(nam);

        PetscPrintf(FERMI_Comm,"initial power:%e\n",power);
        PetscPrintf(FERMI_Comm,"time    power   its\n");

        pNod = calcu.time.head;
        while(pNod){

            dtn = ((tcontrol_t*)pNod->data)->dt;

            fersrods(calcu.t);

            ferass_TR(step);

            fersolv_TR();

            fer_pow(&power);

            PetscPrintf(FERMI_Comm,"%lf %e %d \n", calcu.t, power, its);

            sprintf(nam,"tran_rank%d_t%d", rank, step);
            print_vtk(nam);
            print_out(&phi_n, step);

            calcu.t = calcu.t + dtn;
            step++;
            if(calcu.t > ((tcontrol_t*)pNod->data)->tf - 1.0e-8)
            pNod = pNod->next;

        }
    }

    error:

    ferfini();

    return 0;
}
