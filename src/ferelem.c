/* Assembly functions */

#include "funct.h"

int ferelem_ABM(int e)
{
  /* A,B & M matrices */
  int i,j,d,gp,npe,ngp;
  double det,val,D,xs_a,nxs_f,vel;
  pv_t *pv=(pv_t *)mesh.elemv[e].prop;

  npe=mesh.elemv[e].npe;
  ngp=mesh.elemv[e].ngp;
  nke=npe*egn;

  for(i=0;i<npe;i++)
  {
    for(d=0;d<egn;d++)
      idxm[i*egn+d]=mesh.elemv[e].nodeg[i]*egn+d;
    for(d=0;d<DIM;d++)
      coor[i][d]=mesh.node[mesh.elemv[e].nodel[i]].coor[d];
  }
  memset(Ae,0.0,nke*nke*sizeof(double));
  memset(Be,0.0,nke*nke*sizeof(double));
  memset(Me,0.0,nke*nke*sizeof(double));
  fem_calwei(npe,DIM,&wp);
  fem_calode(npe,DIM,&ode);
  fem_calshp(npe,DIM,&sh);

  for(gp=0;gp<ngp;gp++)
  {
    fem_caljac(coor,ode,npe,gp,DIM,jac);
    fem_invjac(jac,DIM,ijac,&det);
    fem_calder(ijac,npe,DIM,gp,ode,der);
    if(egn==1)
    {
      D=pv->D[0];
      xs_a=pv->xs_a[0];
      nxs_f=pv->nxs_f[0];
      vel=veloc[0];
      for(i=0;i<npe;i++)
      {
        for(j=0;j<npe;j++)
        {
          fem_dotdsh(i,j,der,3,&val);
          Ae[i*npe+j]+= -(D*val + xs_a*sh[i][gp]*sh[j][gp])*wp[gp]*det;
          Be[i*npe+j]+= ikeff*nxs_f*sh[i][gp]*sh[j][gp]*wp[gp]*det;
          Me[i*npe+j]+= sh[i][gp]*sh[j][gp]/(dtn*vel)*wp[gp]*det;
        }
      }
    }
  }/*gp loop*/ 
  return 0;
}

int ferelem_AB(int e)
{
  /* A & B matrices */
  int i,j,d,gp,npe,ngp;
  double det,val,D,xs_a,nxs_f;
  pv_t *pv=(pv_t *)mesh.elemv[e].prop;

  npe=mesh.elemv[e].npe;
  ngp=mesh.elemv[e].ngp;
  nke=npe*egn;
  for(i=0;i<npe;i++)
  {
    for(d=0;d<egn;d++)
      idxm[i*egn+d]=mesh.elemv[e].nodeg[i]*egn+d;
    for(d=0;d<DIM;d++)
      coor[i][d]=mesh.node[mesh.elemv[e].nodel[i]].coor[d];
  }
  memset(Ae,0.0,nke*nke*sizeof(double));
  memset(Be,0.0,nke*nke*sizeof(double));
  fem_calwei(npe,DIM,&wp);
  fem_calode(npe,DIM,&ode);
  fem_calshp(npe,DIM,&sh);

  for(gp=0;gp<ngp;gp++)
  {
    fem_caljac(coor,ode,npe,gp,DIM,jac);
    fem_invjac(jac,DIM,ijac,&det);
    fem_calder(ijac,npe,DIM,gp,ode,der);
    if(egn==1)
    {
      D=pv->D[0];
      xs_a=pv->xs_a[0];
      nxs_f=pv->nxs_f[0];
      for(i=0;i<npe;i++)
      {
        for(j=0;j<npe;j++)
        {
          fem_dotdsh(i,j,der,3,&val);
          Ae[i*npe+j]+=(D*val + xs_a*sh[i][gp]*sh[j][gp])*wp[gp]*det;
          Be[i*npe+j]+=nxs_f*sh[i][gp]*sh[j][gp]*wp[gp]*det;
        }
      }
    }
  }/*gp loop*/ 
  return 0;
}

int ferelem_M(int e)
{
  /* mass matrix */
  int i,j,d,p,gp,npe,ngp;
  double det,val,D,xs_a,nxs_f,vel;
  pv_t *pv=(pv_t *)mesh.elemv[e].prop;

  npe=mesh.elemv[e].npe;
  ngp=mesh.elemv[e].ngp;
  nke=npe*egn;

  for(i=0;i<npe;i++)
  {
    for(d=0;d<egn;d++)
      idxm[i*egn+d]=mesh.elemv[e].nodeg[i]*egn+d;
    for(d=0;d<DIM;d++)
      coor[i][d]=mesh.node[mesh.elemv[e].nodel[i]].coor[d];
  }

  memset(Ae,0.0,nke*nke*sizeof(double));
  memset(Be,0.0,nke*nke*sizeof(double));
  memset(Me,0.0,nke*nke*sizeof(double));
  memset(be,0.0,nke*sizeof(double));
  fem_calwei(npe,DIM,&wp);
  fem_calode(npe,DIM,&ode);
  fem_calshp(npe,DIM,&sh);

  for(gp=0;gp<ngp;gp++)
  {
    fem_caljac(coor,ode,npe,gp,DIM,jac);
    fem_invjac(jac,DIM,ijac,&det);
    fem_calder(ijac,npe,DIM,gp,ode,der);
    if(egn==1)
    {
      D=pv->D[0];
      xs_a=pv->xs_a[0];
      nxs_f=pv->nxs_f[0];
      vel=veloc[0];
      for(i=0;i<npe;i++)
      {
        for(j=0;j<npe;j++)
        {
          fem_dotdsh(i,j,der,3,&val);
          Ae[i*npe+j]+=(D*val + xs_a*sh[i][gp]*sh[j][gp])*wp[gp]*det;
          Be[i*npe+j]+=nxs_f*sh[i][gp]*sh[j][gp]*wp[gp]*det;
          Me[i*npe+j]+= sh[i][gp]*sh[j][gp]/(dtn*vel)*wp[gp]*det;
          if(pv->hasprec)
          {
            for(p=0;p<pgn;p++)
              be[j]+=lambda[p]*pv->conc[p]*wp[gp]*det;
          }
        }
      }
    }
  }/*gp loop*/ 
  return 0;
}

int ferelem_R(int e,double xs_a,double factor)
{
  /* control rod perturbation matrix */
  int i,j,d,gp,npe,ngp;
  double det;

  npe=mesh.elemv[e].npe;
  ngp=mesh.elemv[e].ngp;
  nke=npe*egn;

  for(i=0;i<npe;i++)
  {
    for(d=0;d<egn;d++)
      idxm[i*egn+d]=mesh.elemv[e].nodeg[i]*egn+d;
    for(d=0;d<DIM;d++)
      coor[i][d]=mesh.node[mesh.elemv[e].nodel[i]].coor[d];
  }
  
  memset(Be,0.0,nke*nke*sizeof(double));
  fem_calwei(npe,DIM,&wp);
  fem_calode(npe,DIM,&ode);
  fem_calshp(npe,DIM,&sh);
  
  for(gp=0;gp<ngp;gp++)
  {
    fem_caljac(coor,ode,npe,gp,DIM,jac);
    fem_invjac(jac,DIM,ijac,&det);
    if(egn==1)
    {
      for(i=0;i<npe;i++)
      {
        for(j=0;j<npe;j++)
          Be[i*npe+j] += factor*xs_a*sh[i][gp]*sh[j][gp]*wp[gp]*det;
      }
    }
  }/*gp loop*/ 
  return 0;
}
