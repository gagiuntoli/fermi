/* Solves the linear sistem of equations*/

#include "fermi.h"

int fersolv_ST(void)
{
    int    i;
    double error;

    EPSSetOperators(eps,A,B);
    EPSSolve(eps);
    EPSGetIterationNumber(eps,&its);
    EPSGetDimensions(eps,&nev,NULL,NULL);

    for(i=0;i<nev;i++)
    {
      EPSGetEigenpair(eps,i,&ikeff,NULL,phi_n,NULL);
      keff=1/ikeff;
      EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);
      PetscPrintf(PETSC_COMM_WORLD,"Keff: %lf Its:%d error:%e\n",keff,its,error); 
    }

    return 0;
}

int fersolv_TR(void)
{
//    PetscViewerASCIIOpen(PETSC_COMM_WORLD,"b.dat",&viewer);
//    VecView(b,viewer);
//    PetscViewerASCIIOpen(PETSC_COMM_WORLD,"M.dat",&viewer);
//    MatView(M,viewer);
    KSPSetOperators(ksp,M,M);
    KSPSolve(ksp,b,phi_n);
    KSPGetIterationNumber(ksp,&its);
    return 0;
}
