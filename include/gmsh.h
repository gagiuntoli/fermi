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

#ifndef GMSH_H
#define GMSH_H

typedef struct _gmshN_t{

    int    n;
    double coor[3];

}gmshN_t;

typedef struct _gmshE_t{

    int    node[8];
    int    npe;
    int    gmshid;
    int    elemv;

}gmshE_t;

typedef struct _gmshP_t{

    char   name[32];
    int    dim;
    int    gmshid;
    list_t elem;

}gmshP_t;


int gmsh_read(char *file,char *efile,char *nfile,int rank,int dim, list_t *list_nodes, list_t *list_ghost, list_t
*list_elemv, list_t *list_elems, list_t *list_phyce,int **loc2gold,int **loc2gnew,int **npp,int nproc);
int gmsh_readnodes(char *file,int nproc,char *nfile,int rank,list_t *list_nodes);
int gmsh_detghosts(list_t *list_nodes,list_t *list_elemv,list_t *list_ghost);
int gmsh_readghosts(char *file,list_t *list_ghost);
int gmsh_readelemv(char *file,int nproc,char *nfile,int rank,int dim,list_t *list_elemv);
int gmsh_readelems(char *file,int dim,list_t *list_elemv,list_t *list_elems);
int gmsh_readphys(char *file,list_t *list_elemv,list_t *list_elems,list_t *list_phyce);
int gmsh_reenumerate(char *nfile,int rank,list_t *list_nodes,list_t *list_ghost,int **loc2gold,int **loc2gnew,int **npp,int nproc);
int gmsh_phys_elmlist(list_t *list_elemv, list_t *list_physe);
int gmsh_elems_belongs(list_t *list_elemv,gmshE_t *elems,int *nelemv);
int gmsh_phys_belongs(list_t *list_elemv,list_t *list_elems,gmshP_t *phys);
int gmsh_isvol(int code,int dim);
int gmsh_npe(int code);
int gmsh_nodcmp(void *a, void *b);

#endif
