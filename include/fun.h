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


#ifndef FUN_H
#define FUN_H

#include "stdlib.h"
#include "list.h"

enum {INTER1,INTER2};

typedef struct _f1d_t{

    int n;
    int inter;
    int fnum;

    double *x;
    double *y;

}f1d_t;

int f1d_init(double *x, double *y, int n, int inter, f1d_t *f1d);
int f1d_eval(double x, f1d_t *f1d, double *y);
int cmp_f1d(void *a, void *b);

#endif
