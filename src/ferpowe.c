/* frees all the memory */

#include "funct.h"

int fernorm(void)
{
  // Normalize the flux acording to the power specified (or default)
  double fpower;
  if(ferpowe(&fpower))
    return 1;
  PetscPrintf(PETSC_COMM_WORLD,"Eigen power : %lf\n",fpower); 
  VecScale(phi_n,power/fpower);
  return 0;
}

int ferpowe(double *fpower)
{
  // Calculates the fission power
  int e,g,gp,i,d,error,npe,ngp,locind;
  double power_l=0.0,det,phival;
  pv_t *pv;

  VecGhostUpdateBegin(phi_n,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(phi_n,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(phi_n,&xlocal);
  for(e=0;e<mesh.nelemv;e++)
  {
    pv=(pv_t *)mesh.elemv[e].prop;
    npe=mesh.elemv[e].npe;
    ngp=mesh.elemv[e].ngp;
    for(i=0;i<npe;i++)
    {
      for(d=0;d<DIM;d++)
        coor[i][d]=mesh.node[mesh.elemv[e].nodel[i]].coor[d];
    }
    fem_calwei(npe,DIM,&wp);
    fem_calode(npe,DIM,&ode);
    fem_calshp(npe,DIM,&sh);
    for(i=0;i<npe;i++)
    {
      for(g=0;g<egn;g++)
      {
        locind=mesh.elemv[e].nodel[i]*egn+g;
        VecGetValues(xlocal,1,&locind,&phival);
        for(gp=0;gp<ngp;gp++)
        {
          fem_caljac(coor, ode, npe, gp, DIM, jac);
          fem_invjac(jac, DIM, ijac, &det);
          power_l+=pv->exs_f[g]*phival*sh[i][gp]*det*wp[g];
        }
      }
    }
  }
  VecGhostRestoreLocalForm(phi_n,&xlocal);
  error=MPI_Allreduce(&power_l,fpower,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  if(error)
  {
    PetscPrintf(PETSC_COMM_WORLD,"ferpowe.c:error mpi_allreduce.\n"); 
    return 1;
  }
  return 0;
}
