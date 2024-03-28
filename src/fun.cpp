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

int f1d_init(double *x, double *y, int n, int inter, f1d_t *f1d) {

  int i;
  if (n == 0 || !x || !y || !f1d)
    return 1;
  f1d->n = n;
  f1d->inter = inter;
  f1d->fnum = -1;
  f1d->x = (double *)calloc(n, sizeof(double));
  f1d->y = (double *)calloc(n, sizeof(double));
  if (!f1d->x || !f1d->y)
    return 1;
  for (i = 0; i < n; i++) {
    f1d->x[i] = x[i];
    f1d->y[i] = y[i];
  }
  return 0;
}

int f1d_eval(double x, f1d_t *f1d, double *y) {

  int i;

  if (!f1d)
    return 1;
  if (f1d->n < 1)
    return 1;
  if (f1d->n == 1) {
    *y = f1d->y[0];
    return 0;
  }
  if (x < f1d->x[0]) {
    *y = f1d->y[0];
    return 0;
  }
  i = 1;
  while (i < f1d->n) {
    if (f1d->x[i - 1] <= x && x < f1d->x[i])
      break;
    i++;
  }
  if (i == f1d->n) {
    *y = f1d->y[i - 1];
    return 0;
  }
  *y = (f1d->y[i] - f1d->y[i - 1]) * (x - f1d->x[i - 1]) /
           (f1d->x[i] - f1d->x[i - 1]) +
       f1d->y[i - 1];
  return 0;
}

int cmp_f1d(void *a, void *b) {
  if (((f1d_t *)a)->fnum > ((f1d_t *)b)->fnum) {
    return 1;
  } else if (((f1d_t *)a)->fnum == ((f1d_t *)b)->fnum) {
    return 0;
  } else {
    return -1;
  }
}
