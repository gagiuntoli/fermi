/* structure types and constants
 * 
 */

#undef _TYPES_H_
#define _TYPES_H_

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <slepceps.h>
#include <petsctime.h>
#include <petscao.h>
#include <petscbt.h>
#include <list.h>
#include "petscmat.h"
#include "petscao.h"
#include "petscbt.h"

#include "globals.h"

typedef struct shlist_{
  double **sh;                  // shape function array
  double ***dsh;                // derivatives of shape function array
  double **xp;                  // gauss points
  double *wp;                   // gauss weights
  double **coords;              // coordinate of an element 
  int    npe;                   // number of shapes functions to integrate 
  int    gpn;                   // number of gauss points
  struct shlist_ *next;         
}shlist_t;                      
                                

// typedef struct{
//     element_t      * elem;         // local elements of the mesh
//     surf_element_t * sfelem;       // local elements of the mesh
//     int            num_nod;        // total number of nodes from all the process
//     int            num_elem;       // total number of elements
//     int            num_full_elem;  
//     int            num_sfelem;     // total number of surf elements
//     int            dim;            // dimension of the problem 
//     int            num_nod_l;      
//     int            num_ghost;      
// }mesh_t;

 
// mesh_t * mesh;
shlist_t * shlist;

shlist_t * found_shape(int npe);

double dot(double * v1, double * v2, int n);
double built_jacobian(double **coord, double ***dsh, double jac[MAX_DIM][MAX_DIM], int dim, int n, int gp);
double invert_matrix(double mat[MAX_DIM][MAX_DIM], int size);

void print_initial(void);
void initgp( int n, int o, int dim, double **xp, double *wp);
void init_gauss(void);
void initsh(double **sh, double ***dsh, double **xp, int gpn, int npe, int dim, int o);
void init_shape_list(void);
void elemental_matrix(int e, double *A_e, double *B_e);  // for static problem
void elemental_dynamic_mass(int e, double *Me);          // for transient problem
void set_boundary_condition(void); 
void set_boundary_dynamic_mass(void);
void set_boundary_dynamic_rhs(void);
void normalize_flux(void);


int  node_renumbering(void);
int  determine_ghosts(void);
int  build_boundary_vectors(void);
int  save_vtk(char vtkname[]);
int  save_partition(void);
int  exists(const char *fname);
int  parser1(char *filename);
int  solve_generalized_eigenproblem( Mat *A, Mat *B, Vec * xr, double *keff );
int  solve_linear_problem( Mat *M, Vec *b, Vec * xr );
int  evolute_precursors(void);
int  init_precursors_in_equilibrium(void);
int  exists(const char *fname);
int  readmesh_gmsh(char *filename,int mode);
int  readpartition_metis(char *filename_elements, char *filename_nodes);


int  parser_input(void);

int  gmsh_is_volumetric(int code);
int  gmsh_read_npe_vol(void);
int  gmsh_how_many_nodes(int code);
int  gmsh_read_vol_elements(void);
int  gmsh_read_num_physical(void);
int  gmsh_read_physical(void);
int  gmsh_number_surf_elements(void);
int  gmsh_read_npe_surf(void);
int  gmsh_read_surf_elements(void);

int gmsh_read(char *file,char *efile,char *nfile,int mode,int rank,int dim, list_t *list_coord,list_t *list_elemv,list_t *list_elems,list_t *list_phyce);

int  xstable_read_values(void);
int  xstable_gmsh_number(void);

int  alloc_ownerships(void);
int  alloc_physical_entities(void);
int  alloc(int mode);
int  alloc_vol_elements(void);
int  alloc_xstable(void);
int  alloc_precursors(void);
int  alloc_velocities(void);
int  alloc_npe_surf(void);
int  alloc_surfelements(void);

int  ghost_nodes(void);

int cmp1(void *a, void *b);

int  precursors_read_values(void);

int  velocities_read_values(void);

int  assign_matindex_vol_elements(void);

int  determine_local_materials(void);

int  output_print_structures(void);

int  exist(char *filename);