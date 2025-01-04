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

#include "fermi.hpp"

int ferfini(void) {

  // CLOSING OUTPUT FILES (in those who need it)
  node_list_t *pn;
  output_t *po;
  pn = list_outpu.head;
  while (pn) {
    po = (output_t *)pn->data;
    switch (po->kind) {

    case 1:
      break;

    case 2:
      // power on physical entities as a function of time
      fclose(po->kind_2.fp);
      break;

    default:
      return 1;
    }
    pn = pn->next;
  }

  PetscFinalize();

  return 0;
}
