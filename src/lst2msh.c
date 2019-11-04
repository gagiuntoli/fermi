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

#include "global.h"
#include "types.h"
#include "mesh.h"
#include "list.h"
#include "gmsh.h"
#include "fun.h"

int cpynode (node_list_t *node_nl, node_t *node){

    int d;
    if(!node || !node_nl) return 1;
    for(d=0;d<3;d++)
        node->coor[d]=((gmshN_t *)(node_nl->data))->coor[d];
    return 0;
}


int cpyelemv (node_list_t *elem_nl, elem_t *elemv){

    int gmshid,d,d1;
    char name[32];
    node_list_t *onode;
    pv_t  * pv;
    pvl_t * pvl;

    if(!elem_nl || !elemv)
        return 1;
    elemv->npe = ((gmshE_t *)(elem_nl->data))->npe;
    elemv->ngp = ((gmshE_t *)(elem_nl->data))->npe;
    elemv->nodel=(int *)calloc(elemv->npe,sizeof(int));
    elemv->nodeg=(int *)calloc(elemv->npe,sizeof(int));
    for(d=0;d<elemv->npe;d++)
        elemv->nodeg[d]=((gmshE_t *)(elem_nl->data))->node[d]-1;

    onode = list_physe.head;
    gmshid = ((gmshE_t *)(elem_nl->data))->gmshid;
    while(gmshid != ((gmshP_t*)(onode->data))->gmshid){
        onode=onode->next;
        if(!onode)
            return 2;
    }
    strcpy(name,((gmshP_t*)(onode->data))->name);
    onode = list_mater.head;
    while(strcmp(name,((pvl_t*)(onode->data))->name) != 0){
        onode=onode->next;
        if(!onode)
            return 1;
    }
    pvl = (pvl_t*)onode->data;

    elemv->prop=(pv_t*)calloc(1,sizeof(pv_t));
    pv=(pv_t*)elemv->prop;

    pv->name  = (char*)calloc(sizeof(pvl->name),sizeof(char)); strcpy(pv->name,pvl->name);
    pv->gmshid=gmshid;
    pv->D     = (double*)calloc(egn,sizeof(double));
//    pv->xs_a  = (double*)calloc(egn,sizeof(double));
    pv->xs_s  = (double*)calloc(egn*(egn-1),sizeof(double));
    pv->nxs_f = (double*)calloc(egn,sizeof(double));
    pv->exs_f = (double*)calloc(egn,sizeof(double));
    pv->chi   = (double*)calloc(egn,sizeof(double));

    for(d=0;d<egn;d++){
        pv->D[d]     = pvl->D[d];
        pv->xs_a     = pvl->xs_a;
        for(d1=0;d1<egn-1;d1++)
            pv->xs_s[d*(egn-1)+d1]  = pvl->xs_s[d*(egn-1)+d1];
        pv->nxs_f[d] = pvl->nxs_f[d];
        pv->exs_f[d] = pvl->exs_f[d];
        pv->chi[d]   = pvl->chi[d];
    }

    return 0;
}

int cpyelems (node_list_t *elem_nl, elem_t *elems){

    int d;
    ps_t *ps;

    if(!elem_nl || !elems)
        return 1;
    elems->npe = ((gmshE_t *)(elem_nl->data))->npe;
    elems->ngp = ((gmshE_t *)(elem_nl->data))->npe;
    elems->nodel=(int *)calloc(elems->npe,sizeof(int));
    elems->nodeg=(int *)calloc(elems->npe,sizeof(int));
    for(d=0;d<elems->npe;d++)
        elems->nodeg[d]=((gmshE_t *)(elem_nl->data))->node[d]-1;
    elems->prop=(ps_t*)calloc(1,sizeof(ps_t));
    ps = (ps_t*)elems->prop;
    ps->gmshid = ((gmshE_t*)elem_nl->data)->gmshid;

    return 0;
}
