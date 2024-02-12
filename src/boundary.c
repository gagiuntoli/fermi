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

#include "fermi.h"

int ferbouset(void) {

  /* For every boundary on the list list_bound searchs for
   * the physical entity associated with the same name and assigns
   * the dimS value for each bound.
   * For each boundary search for the nodes that shears the elements
   * that have the same physical entity associated (gmshid). The nodes
   * are saved on the list "nodeL" and are delete in the last part those
   * which are repeated
   *
   * Generates a vector dirIndex with the dirichlet index positions
   */

  int e, n, i, j, k, h, d, gmshid, index;
  node_list_t *pNod, *pNodA, *pba, *pbb, *pna, *pnb;
  list_t dirIndexL;

  pNod = list_bound.head;
  while (pNod) {
    pNodA = list_physe.head;
    while (pNodA) {
      if (!strcmp(((bound_t *)pNod->data)->name,
                  ((gmshP_t *)pNodA->data)->name))
        break;
      pNodA = pNodA->next;
    }
    if (!pNodA) {
      PetscPrintf(FERMI_Comm, "bound.c:boundary %s has no phys entity.\n",
                  ((bound_t *)pNod->data)->name);
      return 1;
    }
    gmshid = ((gmshP_t *)pNodA->data)->gmshid;
    ((bound_t *)pNod->data)->dimS = ((gmshP_t *)pNodA->data)->dim;
    list_init(&((bound_t *)pNod->data)->nodeL, sizeof(int), cmp_nod);
    for (e = 0; e < mesh.nelems; e++) {
      if (gmshid == ((ps_t *)mesh.elems[e].prop)->gmshid) {
        for (n = 0; n < mesh.elems[e].npe; n++)
          list_insert_se(&((bound_t *)pNod->data)->nodeL,
                         (void *)&mesh.elems[e].nodeg[n]);
      }
    }
    pNod = pNod->next;
  }

  for (i = list_bound.sizelist; i > 0; i--) {
    pbb = list_bound.head;
    for (j = 0; j < i - 1; j++) {
      pbb = pbb->next;
    }
    pnb = ((bound_t *)pbb->data)->nodeL.head;
    pba = list_bound.head;
    for (j = 0; j < i - 1; j++) {
      pna = ((bound_t *)pba->data)->nodeL.head;
      for (h = 0; h < ((bound_t *)pba->data)->nodeL.sizelist; h++) {
        pnb = ((bound_t *)pbb->data)->nodeL.head;
        for (k = 0; k < ((bound_t *)pbb->data)->nodeL.sizelist; k++) {
          if (*(int *)pnb->data == *(int *)pna->data) {
            if (list_del(&((bound_t *)pbb->data)->nodeL, pnb)) {
              return 1;
            }
            break;
          }
          pnb = pnb->next;
        }
        pna = pna->next;
      }
      pba = pba->next;
    }
  }

  //
  //==============================
  //   DIRICHLET INDECES
  //==============================
  //
  list_init(&dirIndexL, sizeof(int), cmp_int);
  pba = list_bound.head;
  while (pba) {
    pna = ((bound_t *)pba->data)->nodeL.head;
    while (pna) {
      for (d = 0; d < egn; d++) {
        index = *(int *)pna->data * egn + d;
        if (((bound_t *)pba->data)->kind == 0)
          /* kind = 0 es dirichlet */
          list_insertlast(&dirIndexL, (void *)&index);
      }
      pna = pna->next;
    }
    pba = pba->next;
  }
  ndir = dirIndexL.sizelist;
  dirIndex = (int *)calloc(ndir, sizeof(int));
  dirValue = (double *)calloc(ndir, sizeof(double));
  dirZeros = (double *)calloc(ndir, sizeof(double));
  pna = dirIndexL.head;
  i = 0;
  while (pna) {
    dirIndex[i] = *(int *)pna->data;
    dirZeros[i] = 0.0;
    i++;
    pna = pna->next;
  }
  //
  //==============================
  //   NEUMANN INDECES
  //==============================
  //
  /* por ahora no le hacemos nada*/
  return 0;
}

int cmp_nod(void *a, void *b) {
  if (*(int *)a > *(int *)b) {
    return 1;
  } else if (*(int *)a == *(int *)b) {
    return 0;
  } else {
    return -1;
  }
}
