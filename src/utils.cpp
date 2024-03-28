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

int cmp_int(void *a, void *b) {
  if (*(int *)a > *(int *)b) {
    return 1;
  } else if (*(int *)a == *(int *)b) {
    return 0;
  } else {
    return -1;
  }
}

int cmp_dou(void *a, void *b) {
  if (*(double *)a > *(double *)b) {
    return 1;
  } else if (*(double *)a == *(double *)b) {
    return 0;
  } else {
    return -1;
  }
}

int strBin2Dec(char *str, int *dec) {

  /* Converts the string in "str" that suppose to have a continuue
   * sequence of "0" and "1" to its decimal representation.
   */

  int i;
  *dec = 0;
  for (i = strlen(str) - 1; i >= 0; i--) {
    if (str[i] == '0' || str[i] == '1') {
      if (str[i] == '1')
        *dec += (int)pow(2, strlen(str) - 1 - i);
    } else {
      return 1;
    }
  }
  return 0;
}
