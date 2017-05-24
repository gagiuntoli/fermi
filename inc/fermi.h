/* structure types and constants
 * 
 */
#ifndef _FUNCT_H_
#define _FUNCT_H_

#include <slepceps.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "global.h"
#include "mesh.h"
#include "list.h"
#include "types.h"
#include "gmsh.h"
#include "fem.h"
#include "fun.h"
#include "utils.h"
//parser.c
int parse_input(void);
int parse_mesh(void);
int parse_mats(void);
int parse_mode(void);
int parse_func(void);
int parse_boun(void);
int parse_crod(void);
int parse_outp(void);
int parse_material(char *buff, pvl_t *mat);
int cmp_mat(void *a, void *b);
int parse_boundary(char *bufcpy, bound_t *bou);
int cmp_bou(void *a, void *b);
int cmp_time(void *a, void *b);

int ferass_TR(int step);
int ferass_ST(void);
int ferinit(int argc,char **argv);
int ferfini(void);
int ferstep_ST(void);
int ferstep_TR(int step);
int fersolv_ST(void);
int fersolv_TR(void);
int ferelem_AB(int e);
int ferelem_ABM(int e);
int ferelem_M(int e);
int ferelem_R(int e,double xsa,double factor);
//lst2msh.c
int cpynode (node_list_t *node_nl, node_t *node);
int cpyelemv (node_list_t *elem_nl, elem_t *elemv);
int cpyelems (node_list_t *elem_nl, elem_t *elems);
//ferboun.c
int ferbouset(void);
int cmp_nod(void *a, void *b);
//ferpower.c
int fernorm(void);
int ferpowe(double *fpower);
//ferrods.c
int ferirods(void);
int fersrods(double t);
//output.c
int print_struct(int step);
int print_vtk(const char *name);
int printMatrixR(char *name,double *A, int m, int n);
int printMatrixRC(char *name,double **A, int m, int n);
int printVector(char *name,double *vec, int m);
int printMat(char *name,Mat *A, int m);
int printVec(char *name,Vec *vec, int m);
int vtkcode(int dim,int npe);
#endif
