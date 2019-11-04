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

int fersolv_ST(void)
{
    /* Power Method
    *
    * Iterate with:
    * s = Bx
    * x = A^-1 s
    *
    * At the end:
    * lambda = (Ax . x) / (x . x)
    */

    int its = 0, ksp_its;
    double x_dot_x, Ax_dot_x, lambda;
    double error;
    double power;

    Vec x, s, Ax;

    VecDuplicate(phi_n, &x);
    VecDuplicate(phi_n, &s);
    VecDuplicate(phi_n, &Ax);

    VecSet(x, 1.0);

    while(its < MAX_ITS_POWER) {

        VecSetValues(x, ndir, dirIndex, dirZeros, INSERT_VALUES);
        VecAssemblyBegin(x);
        VecAssemblyEnd(x);

        MatMult(B, x, s);

        KSPSetOperators(ksp, A, A);
        KSPSolve(ksp, s, x);
        KSPGetIterationNumber(ksp, &ksp_its);

        MatMult(A, x, Ax);
        VecDot(x, x, &x_dot_x);
        VecDot(Ax, x, &Ax_dot_x);
        lambda = Ax_dot_x / x_dot_x;
        keff = 1 / lambda;
        PetscPrintf(FERMI_Comm, "Keff: %lf KSP its: %d\n", keff, ksp_its);

        VecCopy(Ax, x);

        VecCopy(x, phi_n);
        fer_pow(&power);
        VecScale(x, 1 / power);

        its++;
    }

    VecCopy(x, phi_n);

    VecDestroy(&x);
    VecDestroy(&s);
    VecDestroy(&Ax);

    return 0;
}

int fersolv_TR(void)
{
//    PetscViewerASCIIOpen(FERMI_Comm,"b.dat",&viewer);
//    VecView(b,viewer);
//    PetscViewerASCIIOpen(FERMI_Comm,"M.dat",&viewer);
//    MatView(M,viewer);
    KSPSetOperators(ksp,M,M);
    KSPSolve(ksp,b,phi_n);
    KSPGetIterationNumber(ksp,&its);
    return 0;
}
