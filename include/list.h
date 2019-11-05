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

#ifndef LIST_H
#define LIST_H

typedef int (*fcmp) (void *, void *);

typedef struct _node_list_t {

    void * data;
    struct _node_list_t * next;

}node_list_t;


typedef struct {

    node_list_t *  head;
    node_list_t *  tail;
    int            sizedata;
    int            sizelist;
    fcmp           cmp;

}list_t;

int list_init(list_t * list, int sizedata, fcmp cmp);
int list_insert_se(list_t * list, void *data);
int list_insertlast(list_t * list, void *data);
int list_delfirst(list_t * list);
int list_del(list_t *list, node_list_t* pNod);
int list_free(list_t *list);

#endif
