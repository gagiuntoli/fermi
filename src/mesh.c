/* Rutines to performe calculations on a mesh
 * 
 * 
 * list_t *list_nodes list of nodes (locals)
 * list_t *list_ghost list of ghost nodes
 * list_t *list_elemv list of volumetric elements
 * 
 * cpyprop_t cpypropv : pointer to funtion to copy properties from list_t *list_elemvprop to elemv
 * cpyprop_t cpyprops : pointer to funtion to copy properties from list_t *list_elemsprop to elems
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mesh.h"

int mesh_alloc(list_t *list_nodes, list_t *list_ghost, cpynode_t cpynode, list_t *list_elemv, cpyelem_t cpyelemv, list_t *list_elems, cpyelem_t cpyelems, mesh_t *mesh){

  int i;

  if(!mesh)
    return 1;
  mesh->nnodes=list_nodes->sizelist;
  mesh->nghost=list_ghost->sizelist;
  mesh->node=(node_t*)calloc((mesh->nnodes+mesh->nghost),sizeof(node_t));
  if(!mesh->node)
    return 1;
  for(i=0;i<mesh->nnodes;i++)
  {
    cpynode(list_nodes->head,&mesh->node[i]);
    if(list_delfirst(list_nodes))
      return 1;
  }
  for(i=0;i<mesh->nghost;i++)
  {
    cpynode(list_ghost->head,&mesh->node[mesh->nnodes+i]);
    if(list_delfirst(list_ghost))
      return 1;
  }
  mesh->nelemv=list_elemv->sizelist;
  mesh->elemv=(elem_t *)calloc(mesh->nelemv,sizeof(elem_t));
  if(!mesh->elemv)
    return 1;
  for(i=0;i<mesh->nelemv;i++)
  {
    if(cpyelemv(list_elemv->head,&mesh->elemv[i]))
      return 1;
    if(list_delfirst(list_elemv))
      return 1;
  }
  mesh->nelems=list_elems->sizelist;
  mesh->elems=(elem_t *)calloc(mesh->nelems,sizeof(elem_t));
  if(!mesh->elems)
    return 1;
  for(i=0;i<mesh->nelems;i++)
  {
    if(cpyelems(list_elems->head,&mesh->elems[i]))
      return 1;
    if(list_delfirst(list_elems))
      return 1;
  }    
  return 0;    
}


int mesh_renum(mesh_t *mesh, int *loc2gold, int *loc2gnew){

  /* Renumbers the nodes of each element acording to the loc2glo vector*/

  int e,n,i,fl;
  for(e=0;e<mesh->nelemv;e++){
    for(n=0;n<mesh->elemv[e].npe;n++){
      fl=0;
      for(i=0;i<(mesh->nnodes + mesh->nghost);i++){
        if(loc2gold[i]==mesh->elemv[e].nodeg[n]){
          mesh->elemv[e].nodeg[n]=loc2gnew[i];
          mesh->elemv[e].nodel[n]=i;
          fl=1;
          break;
        }
      }
      if(!fl)
        return 1;
    } 
  }
  for(e=0;e<mesh->nelems;e++){
    for(n=0;n<mesh->elems[e].npe;n++){
      fl=0;
      for(i=0;i<(mesh->nnodes + mesh->nghost);i++){
        if(loc2gold[i]==mesh->elems[e].nodeg[n]){
          mesh->elems[e].nodeg[n]=loc2gnew[i];
          mesh->elems[e].nodel[n]=i;
          fl=1;
          break;
        }
      }
      if(!fl)
        return 1;
    } 
  }

  return 0;
}

int mesh_neigh(mesh_t *mesh, int *loc2gnew){

  /* completes the mesh.node[i].elemv & mesh.node[i].elems lists
   * for both : local nodes and ghosts nodes. The lists contains 
   * the local element numeration of those which have the node 
   * as a vertex.
   */
  int i,e,n;

  for(i=0;i<(mesh->nnodes + mesh->nghost);i++){
    list_init(&mesh->node[i].elemvL,sizeof(int),elem_cmp);
    for(e=0;e<mesh->nelemv;e++){
      for(n=0;n<mesh->elemv[e].npe;n++){
        if(mesh->elemv[e].nodeg[n]==loc2gnew[i]){
          list_insert_se(&mesh->node[i].elemvL,(void*)&e);
          break;
        }
      } 
    }
  }

  for(i=0;i<(mesh->nnodes + mesh->nghost);i++){
    list_init(&mesh->node[i].elemsL,sizeof(int),elem_cmp);
    for(e=0;e<mesh->nelems;e++){
      for(n=0;n<mesh->elems[e].npe;n++){
        if(mesh->elems[e].nodeg[n]==loc2gnew[i]){
          list_insert_se(&mesh->node[i].elemsL,(void*)&e);
          break;
        }
      } 
    }
  }
  return 0;    
}


int mesh_carea(mesh_t *mesh,elem_t *elem,int dim,double *area){

  int d;
  double v1[3],v2[3],v3[3],vr[3],mod;
  *area=0.0;

  if(elem->npe==2){
    for(d=0;d<3;d++){
      v1[d]=mesh->node[elem->nodel[1]].coor[d]-mesh->node[elem->nodel[0]].coor[d];
    }
    mesh_vnorm(vr,3,&mod);
    *area+=mod;
    return 0;
  }else if(elem->npe==3){
    for(d=0;d<3;d++){
      v1[d]=mesh->node[elem->nodel[1]].coor[d]-mesh->node[elem->nodel[0]].coor[d];
      v2[d]=mesh->node[elem->nodel[2]].coor[d]-mesh->node[elem->nodel[0]].coor[d];
    }
    mesh_vcross(v1,v2,vr);
    mesh_vnorm(vr,3,&mod);
    *area+=mod/2;
    return 0;
  }else if(elem->npe==4){
    for(d=0;d<3;d++){
      v1[d]=mesh->node[elem->nodel[1]].coor[d]-mesh->node[elem->nodel[0]].coor[d];
      v2[d]=mesh->node[elem->nodel[2]].coor[d]-mesh->node[elem->nodel[0]].coor[d];
      v3[d]=mesh->node[elem->nodel[3]].coor[d]-mesh->node[elem->nodel[0]].coor[d];
    }
    mesh_vcross(v1,v2,vr);
    mesh_vnorm(vr,3,&mod);
    *area+=mod/2;
    mesh_vcross(v3,v2,vr);
    mesh_vnorm(vr,3,&mod);
    *area+=mod/2;
    return 0;
  }
  return 1;
}

int mesh_vnorm(double *vec, int n, double *mod){

  int d;
  if(!vec || !mod)
    return 1;
  *mod=0.0; 
  for(d=0;d<n;d++)
    *mod+=pow(vec[d],2);
  *mod=sqrt(*mod);
  return 0;
}

int mesh_vcross(double *v1, double *v2, double *vr){

  if(!v1 || !v2 || !vr)
    return 1;
  vr[0]= v1[1]*v2[2]-v2[1]*v1[2];
  vr[1]= v1[0]*v2[2]-v2[2]*v1[0];
  vr[2]= v1[0]*v2[1]-v2[0]*v1[1];
  return 0;
}

int elem_cmp(void *a, void *b){
  if ( *(int*)a > *(int*)b){
    return 1;
  }else if(*(int*)a == *(int*)b){
    return 0;
  }else{
    return -1;
  }
}
