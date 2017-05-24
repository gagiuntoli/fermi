#include "fermi.h"

/*************************************************************/
int fernorm(void){

  /* Normalize the flux acording to the power specified (or default) */

  double fpower;

  if(ferpowe(&fpower))
    return 1;

  PetscPrintf(PETSC_COMM_WORLD,"Eigen power : %lf\n",fpower); 
  VecScale(phi_n,power/fpower);

  return 0;
}

/*************************************************************/

int ferpowe(double *fpower){

  /* Calculates the fission power */

  int      e,g,gp,i,d;
  int      error,npe,ngp,locind;
  double   power_l;                /* potencia local dentro de este proceso */ 
  double   det,phival;
  pv_t   * pv;

  power_l=0.0;

  VecGhostUpdateBegin(phi_n,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(phi_n,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(phi_n,&xlocal);
  
  /* recorremos todos los elementos de este proceso */
  for(e=0;e<mesh.nelemv;e++){

    pv  = (pv_t *)mesh.elemv[e].prop;
    npe = mesh.elemv[e].npe;         /* numero de vertices por elemento */
    ngp = mesh.elemv[e].ngp;         /* numero de puntos de Gauss por elemento */

    /* armamos un vector con las coordenadas de los vertices para calcular el jac */
    for(i=0;i<npe;i++){
      for(d=0;d<DIM;d++)
        coor[i][d]=mesh.node[mesh.elemv[e].nodel[i]].coor[d];
    }

    /* buscamos los pesos, las funciones de forma y las derivadas */
    fem_calwei(npe,DIM,&wp);
    fem_calode(npe,DIM,&ode);
    fem_calshp(npe,DIM,&sh);

    /* Recorremos los vertices del elemento para calcular el aporte en la integral de cada funcion de forma */
    for(i=0;i<npe;i++){

      /* Recorremos cada una de las energias del neutron */
      for(g=0;g<egn;g++){

        /* calculamos el indice donde está el valor del flujo en el vector distribuido y lo pedimos */
        locind=mesh.elemv[e].nodel[i]*egn+g;
        VecGetValues(xlocal,1,&locind,&phival);

        /* sumamos la contribución a la integral total */
        for(gp=0;gp<ngp;gp++){
          fem_caljac(coor, ode, npe, gp, DIM, jac);
          fem_invjac(jac, DIM, ijac, &det);
          power_l+=pv->exs_f[g]*phival*sh[i][gp]*det*wp[g];
        }
      }

    }

  }
  
  /* Hacemos el Allreduce de las power_l de cada proceso y metemos el resultado en fpower */
  error = MPI_Allreduce(&power_l,fpower,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  if(error){
    PetscPrintf(PETSC_COMM_WORLD,"ferpowe.c:error mpi_allreduce.\n"); 
    return 1;
  }

  return 0;
}
/*************************************************************/
