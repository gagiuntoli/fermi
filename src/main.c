/* Bathe main program*/

#include "funct.h"

int main(int argc,char **argv){

  int step;
  char nam[32];
  node_list_t *pNod;

  if(ferinit(argc,argv))
    goto error;

  calcu.t = calcu.t0;
  step=0;
  
  if(calcu.timedep == QS)
  {

    pNod=calcu.time.head;
    while(pNod)
    {

      dtn=((tcontrol_t*)pNod->data)->dt;

      if(ferstep_ST())
        goto error;

      sprintf(nam,"steady_r%d_t%d",rank,step);
      print_vtk(nam);

      calcu.t=calcu.t + dtn;
      step ++;
      if(calcu.t>((tcontrol_t*)pNod->data)->tf-1.0e-8)
        pNod=pNod->next;

    }

  }else if(calcu.timedep == TR){

    PetscPrintf(PETSC_COMM_WORLD,"calculating stationary state.\n");
    if(ferstep_ST())
      goto error;

    sprintf(nam,"steady_r%d",rank);
    print_vtk(nam);

    PetscPrintf(PETSC_COMM_WORLD,"initial power:%e\n",power);
    PetscPrintf(PETSC_COMM_WORLD,"time    power   its\n");

    pNod=calcu.time.head;
    while(pNod)
    {

      dtn=((tcontrol_t*)pNod->data)->dt;

      if(ferstep_TR(step))
        goto error;

      if(ferpowe(&power))
        goto error;

      PetscPrintf(PETSC_COMM_WORLD,"%lf %e %d \n",calcu.t,power,its);

      sprintf(nam,"tran_rank%d_t%d",rank,step);
      print_vtk(nam);

      calcu.t=calcu.t + dtn;
      step ++;
      if(calcu.t>((tcontrol_t*)pNod->data)->tf-1.0e-8)
        pNod=pNod->next;

    }
  }

error:

  ferfini();

  return 0;
}
