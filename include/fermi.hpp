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

#ifndef _FERMI_H_
#define _FERMI_H_

#include <mpi.h>
#include <petscksp.h>

#include <stdbool.h>
#include <unistd.h>
#include <vector>

struct Material {
  std::vector<double> D;    // diffusion coefficient
  std::vector<double> xs_a; // absorbsion XS
  std::vector<double> xs_f; // fission XS
  std::vector<double> xs_s; // scattering XS
  std::vector<double> chi;  // factor determing with % of all fissions ends in the groups

  bool operator==(const Material& other) const {
        return D == other.D && xs_a == other.xs_a && xs_f == other.xs_f && xs_s == other.xs_s && chi == other.chi;
  }
};

enum Calculation {
  Keff
};

enum BoundaryCondition {
  Dirichlet,
  Neumann,
};

#define DIM 3
#define NPE 8
#define MAX_ITS_POWER 200

enum { INTER1, INTER2 };

typedef struct f1d_t_ {

  int n;
  int inter;
  int fnum;

  double *x;
  double *y;

} f1d_t;

int f1d_init(double *x, double *y, int n, int inter, f1d_t *f1d);
int f1d_eval(double x, f1d_t *f1d, double *y);
int cmp_f1d(void *a, void *b);

typedef int (*fcmp)(void *, void *);

typedef struct _node_list_t {

  void *data;
  struct _node_list_t *next;

} node_list_t;

typedef struct {

  node_list_t *head;
  node_list_t *tail;
  int sizedata;
  int sizelist;
  fcmp cmp;

} list_t;

typedef struct _gmshN_t {

  int n;
  double coor[3];

} gmshN_t;

typedef struct _gmshE_t {

  int node[8];
  int npe;
  int gmshid;
  int elemv;

} gmshE_t;

typedef struct _gmshP_t {

  char name[32];
  int dim;
  int gmshid;
  list_t elem;

} gmshP_t;

int gmsh_read(char *file, char *efile, char *nfile, int rank, int dim,
              list_t *list_nodes, list_t *list_ghost, list_t *list_elemv,
              list_t *list_elems, list_t *list_phyce, int **loc2gold,
              int **loc2gnew, int **npp, int nproc);
int gmsh_readnodes(char *file, int nproc, char *nfile, int rank,
                   list_t *list_nodes);
int gmsh_detghosts(list_t *list_nodes, list_t *list_elemv, list_t *list_ghost);
int gmsh_readghosts(char *file, list_t *list_ghost);
int gmsh_readelemv(char *file, int nproc, char *nfile, int rank, int dim,
                   list_t *list_elemv);
int gmsh_readelems(char *file, int dim, list_t *list_elemv, list_t *list_elems);
int gmsh_readphys(char *file, list_t *list_elemv, list_t *list_elems,
                  list_t *list_phyce);
int gmsh_reenumerate(char *nfile, int rank, list_t *list_nodes,
                     list_t *list_ghost, int **loc2gold, int **loc2gnew,
                     int **npp, int nproc);
int gmsh_phys_elmlist(list_t *list_elemv, list_t *list_physe);
int gmsh_elems_belongs(list_t *list_elemv, gmshE_t *elems, int *nelemv);
int gmsh_phys_belongs(list_t *list_elemv, list_t *list_elems, gmshP_t *phys);
int gmsh_isvol(int code, int dim);
int gmsh_npe(int code);
int gmsh_nodcmp(void *a, void *b);

int list_init(list_t *list, int sizedata, fcmp cmp);
int list_insert_se(list_t *list, void *data);
int list_insertlast(list_t *list, void *data);
int list_delfirst(list_t *list);
int list_del(list_t *list, node_list_t *pNod);
int list_free(list_t *list);

typedef struct _elem_t {

  int npe;
  int ngp;
  int *nodel;
  int *nodeg;
  void *prop;

} elem_t;

typedef struct _node_t {

  double coor[3];
  list_t elemvL;
  list_t elemsL;

} node_t;

typedef struct _mesh_t {

  int nelemv;
  int nelems;
  int nnodes;
  int nghost;
  node_t *node;
  elem_t *elemv;
  elem_t *elems;

} mesh_t;

typedef int (*cpyelem_t)(node_list_t *elem_nl, elem_t *elem);
typedef int (*cpynode_t)(node_list_t *node_nl, node_t *node);

enum { SEQUENCIAL, PARALLEL };

typedef struct bit {

  unsigned int val : 1;

} bit_t;

typedef struct _propES_t {

  int gmshid; /* ID that correspond to evewry element in the gmshfile */
  int elemv;  /* The elemv local numeration at which belongs */

} ps_t;

typedef struct _bound_t {

  char name[16];
  int kind;
  int order;
  int dimS;
  list_t nodeL;
  list_t elemsL;
  list_t elemvL;

} bound_t;

typedef struct _pv_t {

  char *name;
  int gmshid;

  double *D;     // diffusion coeficient
  double *xs_a;  // absortion XS
  double *nxs_f; // nu x fission XS
  double *exs_f; // energy x fission XS
  double *xs_s;  // scattering XS
  double *xs_r;  // remotion XS
  double *chi;   // fission spectrum

  int hasprec;
  double *conc; // Fission precursors concentration  ( I groups )

} pv_t;

typedef struct _pvl_t {

  char name[16];

  double *D;     // diffusion coeficient
  double *xs_a;  // absortion XS
  double *nxs_f; // nu x fission XS
  double *exs_f; // energy x fission XS
  double *xs_s;  // scattering XS
  double *xs_r;  // remotion XS
  double *chi;   // fission spectrum

  int hasprec;

} pvl_t;

/*************************************************************/
typedef struct _kind_1_t {

  char phys[16];

} kind_1_t;

/* this is used to print localized powers on physical entities on a file */
typedef struct _kind_2_t {

  FILE *fp;
  char file[16];
  char **phys; // array of Physical Entities names
  int nphy;
  int *ids;
  double *pow;

} kind_2_t;

typedef struct _output_t {

  /* esto va a volar pronto */
  char file[16];
  char phys[16];
  double norm[3];

  int kind;

  kind_1_t kind_1;
  kind_2_t kind_2;

} output_t;

typedef struct _comm_1_t {

  char friend_name[64];
  int nphy;
  char **phys; // array of Physical Entities names
  int *ids;
  int remote_rank;
  ;

  /* recv */
  double *xs;

  /* send */
  double *pow;

  MPI_Comm *intercomm;

} comm_1_t;

typedef struct _comm_t {

  int kind;

  comm_1_t comm_1;

} comm_t;

/*************************************************************/

typedef struct _tcontrol_t {

  double tf;
  double dt;

} tcontrol_t;

typedef struct _ctrlrod_t {

  char name_ele[16]; /* elem physical entity name */
  char name_nod[16]; /* node physical entity name */
  int nfun;          /* function id number */
  f1d_t *funins;     /* insertion value function of time */
  double n[3];       /* direction of control rod insertion */
  double p[3];       /* reference point from where insertion starts */
  double xsaval;     /* xsa value to perturb  */
  list_t elemv;      /* list of vol elem that the ctrlrod intersects*/
  list_t xsa;        /* list of absortion xs on element list*/

} ctrlrod_t;

typedef struct {

  double t0;
  double t;
  list_t time;

  int timedep;
  int kmode;
  int mode;
  int exec;

} calcu_t;

/*************************************************************/

// parser.c
int parse_input(void);
int parse_mesh(void);
int parse_mats(void);
int parse_mode(void);
int parse_func(void);
int parse_boun(void);
int parse_crod(void);
int parse_outp(void);
int parse_communication(void);
int parse_material(char *buff, pvl_t *mat);
int cmp_mat(void *a, void *b);
int parse_boundary(char *bufcpy, bound_t *bou);
int cmp_bou(void *a, void *b);
int cmp_time(void *a, void *b);
int get_int(char *buf, const char *name, int *a);
int get_char(char *buf, const char *name, char *a);

int ferass_TR(int step);
int ferass_ST(void);
int ferinit(int argc, char **argv);
int ferfini(void);
int ferstep_ST(void);
int ferstep_TR(int step);
int fersolv_ST(void);
int fersolv_TR(void);
int ferelem_AB(int e);
int ferelem_ABM(int e);
int ferelem_M(int e);
int ferelem_R(int e, double xsa, double factor);
// lst2msh.c
int cpynode(node_list_t *node_nl, node_t *node);
int cpyelemv(node_list_t *elem_nl, elem_t *elemv);
int cpyelems(node_list_t *elem_nl, elem_t *elems);
// ferboun.c
int ferbouset(void);
int cmp_nod(void *a, void *b);

// fer_power.c
int fer_norm(void);
int fer_pow(double *fpower);
int fer_pow_phys(int n, int *ids, double *fpower);
int fer_pow_elem(int e, double *fpower);

// ferrods.c
int ferirods(void);
int fersrods(double t);
// output.c
int print_struct(int step);
int print_vtk(const char *name);
int printMatrixR(char *name, double *A, int m, int n);
int printMatrixRC(char *name, double **A, int m, int n);
int printVector(char *name, double *vec, int m);
int printMat(char *name, Mat *A, int m);
int printVec(char *name, Vec *vec, int m);
int vtkcode(int dim, int npe);
int print_out(Vec *phi, int step);

int fer_comm_init(void);
int fer_comm_step(int order);

extern MPI_Comm WORLD_Comm; // global communicator
extern MPI_Comm FERMI_Comm; // local  communicator

extern int globa_rank; // rank in WORLD_Comm
extern int globa_size; // size in WORLD_Comm
extern int local_rank; // rank in FERMI_Comm
extern int local_size; // size in FERMI_Comm

enum { QS, TR }; // Quasi Static, transient
enum { K1 };     // Calculation of elemental matrix by this K modes

extern PetscViewer kspview;
extern PetscViewer viewer;

extern int rank;
extern int nproc;
extern int nummat;

extern int *loc2gold, *loc2gnew;
extern int ndir;
extern int *npp;
extern int ntot;
extern int *dirIndex;
extern double *dirValue;
extern double *dirZeros;

extern int memory;

extern list_t list_nodes;
extern list_t list_ghost;
extern list_t list_elemv;
extern list_t list_elems;
extern list_t list_physe;
extern list_t list_mater;
extern list_t list_bound;
extern list_t list_fun1d; /* list of functions */
extern list_t list_ctrlr; /* list of control rods */
extern list_t list_outpu;
extern list_t list_comms;

extern mesh_t mesh;

extern int nke, nbe, nev, its;
extern double *Ae, *Be, *Me, *be;
extern Vec phi_n, phi_o, b, b_a;
extern Vec xlocal;
extern double ikeff, keff;
extern Mat A, B, M, K;
extern KSP ksp;
extern PC pc; /* preconditioner context */
extern int Istart, Iend;

extern char inputfile[32];
extern char meshfile[32];
extern char epartfile[32];
extern char npartfile[32];

extern calcu_t calcu;

extern int *idxm;
extern double **der, ***ode, **sh, **jac, **ijac, **coor, *wp;

extern int egn, pgn;
extern int nxs_mat; // number of xs values per material
extern double **phi;
extern double power;
extern double vol;

/* Precursors constants */
extern double *veloc, *lambda, *beta, *chi, beta_tot;
extern double dtn;

int cmp_int(void *a, void *b);
int cmp_dou(void *a, void *b);
int strBin2Dec(char *str, int *dec);

int cmp_int(void *a, void *b);
int cmp_dou(void *a, void *b);
int strBin2Dec(char *str, int *dec);

int mesh_alloc(list_t *list_nodes, list_t *list_ghost, cpynode_t cpynode,
               list_t *list_elemv, cpyelem_t cpyelemv, list_t *list_elems,
               cpyelem_t cpyelems, mesh_t *mesh);
int mesh_renum(mesh_t *mesh, int *loc2gold, int *loc2gnew);
int mesh_neigh(mesh_t *mesh, int *loc2gnew);
int mesh_vnorm(double *vec, int n, double *mod);
int mesh_carea(mesh_t *mesh, elem_t *elem, int dim, double *area);
int mesh_vcross(double *v1, double *v2, double *vr);
int elem_cmp(void *a, void *b);

int fem_inigau(void);
int fem_caljac(double **coor, double ***ds, int npe, int gp, int dim,
               double **jac);
int fem_invjac(double **jac, int dim, double **ijac, double *det);
int fem_calode(int npe, int dim, double ****oder);
int fem_calshp(int npe, int dim, double ***sh);
int fem_calwei(int npe, int dim, double **wp);
int fem_calder(double **ijac, int nsh, int dim, int gp, double ***oder,
               double **der);
int fem_calare(double **pts, int npe, int dim, double *area);
int fem_vecmod(double *vec, int n, double *mod);
int fem_dotdsh(int i, int j, double **derivs, int dim, double *p);
int fem_vcross(double *v1, double *v2, double *vr);

// Segment 2 nodes
extern double **xp_segm_2;
extern double *wp_segm_2;
extern double **sh_segm_2;
extern double ***ds_segm_2;

// Triangle 3 nodes
extern double **xp_tria_3;
extern double *wp_tria_3;
extern double **sh_tria_3;
extern double ***ds_tria_3;

// Quadrangle 4 nodes
extern double **xp_quad_4;
extern double *wp_quad_4;
extern double **sh_quad_4;
extern double ***ds_quad_4;

// Tetrahedron 4 nodes
extern double **xp_tetra_4;
extern double *wp_tetra_4;
extern double **sh_tetra_4;
extern double ***ds_tetra_4;

// Prism 6 nodes
extern double **xp_prism_6;
extern double *wp_prism_6;
extern double **sh_prism_6;
extern double ***ds_prism_6;

// Hexahedron 8 nodes
extern double **xp_hexa_8;
extern double *wp_hexa_8;
extern double **sh_hexa_8;
extern double ***ds_hexa_8;

#endif
