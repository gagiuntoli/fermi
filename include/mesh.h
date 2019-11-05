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


#ifndef MESH_H
#define MESH_H


#include "list.h"
#include <math.h>


typedef struct _elem_t{

    int  npe;
    int  ngp;
    int  *nodel;
    int  *nodeg;
    void *prop;

}elem_t;

typedef struct _node_t{

    double coor[3];
    list_t elemvL;
    list_t elemsL;

}node_t;

typedef struct _mesh_t{

    int nelemv;
    int nelems;
    int nnodes;
    int nghost;
    node_t *node;
    elem_t *elemv;
    elem_t *elems;

}mesh_t;

typedef int (*cpyelem_t) (node_list_t *elem_nl, elem_t *elem);
typedef int (*cpynode_t) (node_list_t *node_nl, node_t *node);

int mesh_alloc(list_t *list_nodes, list_t *list_ghost, cpynode_t cpynode, list_t *list_elemv, cpyelem_t cpyelemv, list_t *list_elems, cpyelem_t cpyelems, mesh_t *mesh);
int mesh_renum(mesh_t *mesh, int *loc2gold, int *loc2gnew);
int mesh_neigh(mesh_t *mesh, int *loc2gnew);
int mesh_vnorm(double *vec, int n, double *mod);
int mesh_carea(mesh_t *mesh,elem_t *elem,int dim,double *area);
int mesh_vcross(double *v1, double *v2, double *vr);
int elem_cmp(void *a, void *b);

#endif
