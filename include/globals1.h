#include "petscao.h"
#include <list.h>

#ifndef _GLOBALS_H_
#define _GLOBALS_H_

// MACROS definitions

#define COORD_TOL   1.0e-9            // Minimun separation between nodes in each coordinate

#define GPN_1_2     2                 // Gauss point number for segment     linear shapes funtions
#define GPN_2_3     3                 // Gauss point number for triangle    linear shapes funtions
#define GPN_2_4     4                 // Gauss point number for quadrangle  linear shapes funtions
#define GPN_3_4     4                 // Gauss point number for tetrahedron linear shapes funtions
#define GPN_3_8     8                 // Gauss point number for hexahedron  linear shapes funtions

#define MAXL        32                // Maximum length of command line argument

#define MAX_DIM     3                 // Maximum dimension of the problem
#define MAX_EGN     3                 // Maximum number of energy groups
#define MAX_NPE     8                 // Maximum number of nodes per element
#define MAX_NPSE    4                 // Maximum number of nodes per surface element
#define MAX_GPN     8                 // Maximum number of Gauss points used in integration

#define MAX_MAT     8                 // Maximum number of materials defined on input

// Node Type Options

#define INSIDE_VOL  0                 // This node is not in the border
#define DIRICHLET   1                 // This node has a Dirichlet condition
#define NEUMANN     2                 // This node has a Neumann condition
#define ROBIN       3                 // This node has a Robin condition

// Errors

#define COMMAND_ERROR 1

// Global Variables

#define SEQUEN           1
#define PARALL           2
#define SEGMENT_2        3
#define TRIANGLE_3       4
#define QUAD_4           5
#define TETRAHEDRON_4    6
#define HEXAEDRON_8      7
// #define UNKNOWN          1000

#define BUFFER_LENGHT    128
#define BIGBUFFER_LENGHT 1000
#define WORDS_LENGHT     128

enum {GMSH_KIND};
enum {METIS_KIND};
enum {STATIC, TRANSIENT};
enum {SEQUENCIAL, PARALLEL};

Mat     A;
Mat     B;
Mat     M;     //Mass Matrix for dynamic problem

KSP      ksp;

ST       st;

Vec      xr;
Vec      xr_aux;
Vec      xi;
Vec      *Iv;
Vec      *Cv;
Vec      lxr;
Vec      b;     //RHS vector for dynamic problem
Vec      vout;

PetscInt  Istart;
PetscInt  Iend;
PetscInt  its;
PetscInt  lits;
PetscInt  nev;
PetscInt  maxit;
PetscInt  Istartb;
PetscInt  Iendb;
PetscInt  Istartxr;
PetscInt  Iendxr;

PetscReal tol;

PetscLogDouble t1,
        t2,
        t3,
        t4,
        t5,
        t6,
        dt1,
        dt2,
        dt3,
        dt4;

PetscViewer    viewer;

PetscBool      flg;
PetscBool      evecs;
PetscBool      ishermitian;
PetscBool      terse;
PetscBool      set;
PetscBool      mesh_flag;
PetscBool      mat_flag;
PetscBool      vtk_flag;
PetscBool      epmv_flag;
PetscBool      npmv_flag;
PetscBool      epms_flag;

PetscErrorCode ierr;

ISLocalToGlobalMapping ltog;

AO ao;

VecScatter *ctx;


double  **coord_ghosts;
double  *xvalues;
double  keff;
double  power;
double  energy_per_fision;
double  nu;
double  *velocities;  // (egn components)
double  dt;           // time step for dynamic simulations
double  time;         // total time up to now
double  **cross_sections;

int     Istart;
int     Iend;

int     num_local;     // number of local nodes
int     num_ghost;     // number of ghost nodes
int     num_static_calc;
int     *glo2loc;
int     *loc2glo;
int     *ghosts;
int     *are_mine;

int     num_dirichlet;
int     *is_dirichlet;
int     *dirichlet;
int     *dirichlet_rows;
int     num_neumann;
int     *is_neumann;
int     *neumann;
int     *neumann_rows;
int     *global_new;
int     *nnpp;
int     *is_ghost;             // Number of nodes per process*/
int     integration_form;      // Kind of flux averange for precursors calculations
int     num_total_vol_elem;
int     num_total_surf_elem;

/***************************************************************************/

typedef struct bit{

    unsigned int val : 1;

}bit_t;

typedef struct{

    unsigned long    num_global;            // global numeration of the node
    unsigned long    num_local;             // local numeration of the node
    bit_t            mine;                  // if belongs to this process is 1 otherwise is 0

}node_t;

typedef struct{

    unsigned int    npe;                      // Number of nodes of this element
    unsigned long   mat_index;
    unsigned long   gmsh_number;              //number that appear in the Gmsh file

    bit_t           mine;                     // 1 if all its nodes are from this process
    node_t          *node;                    // Array of nodes of the element

}element_vol_t;

typedef struct{

    unsigned int    npe;                    // Number of nodes of this element

    unsigned int    boundary;
    unsigned long   gmsh_number;

    bit_t           mine;
    node_t          *node;                  // Array of nodes of the element

}element_surf_t;

typedef struct{

    char             name[64];
    unsigned long    gmsh_number;               // gmsh code used to relate this with the elements in the mesh

    double           *D;                        // diffusion coeficient
    double           *xs_a;                     // absortion XS
    double           *nxs_f;                    // nu x fission XS
    double           *exs_f;                    // energy x fission XS
    double           *xs_s;                     // scattering XS
    double           *xs_r;                     // remotion XS
    double           *chi;                      // fission spectrum

    bit_t            has_precursors;
    double           *conc;                     // Fission precursors concentration  ( I groups )

}xstable_t;

typedef struct{

    double        *lambda;                   // Fission precursors decay constant ( I groups )
    double        *beta;                     // Fission precursors fission yield  ( I groups )
    double        *chi;           // chi constants of precursors       ( I groups x egn )
    double        beta_tot;                  // Fission precursors total fission yield

}precursors_t;

typedef struct{

    double t0;
    double tf;
    double t;
    double dt;


    int mode;           // STATIC || TRANSIENT
    int exec;
    int static_steps;

}calculation_flow_t;

typedef struct{

    unsigned long  type;                //if belongs to surf or to vol
    unsigned long  gmsh_number;         //number that appear in the Gmsh file
    char           name[64];            //nombre que le pongo yo a esa entidad
    bit_t          mine;                //one if some of the elements of this process has this material

}physical_entity_t;


int                 execution_mode;         // SEQUENCIAL | PARALLEL
int                 rank;
int                 num_process;   // number of processes
int                 flag_coor;
unsigned int        dim;

int                 mesh_kind;
int                 partition_kind;

unsigned long       num_physical_entities;
unsigned long       num_materials;          //number of local materials
unsigned int        pgn;                    //number of energy groups
unsigned int        egn;                    //number of precursors groups

bit_t               *is_my_vol_element;
bit_t               *is_my_surf_element;
bit_t               *is_my_node;

physical_entity_t   *physical_entity;

unsigned long       num_elements_vol;
element_vol_t       *element_vol;
unsigned int        *npe_vol;

unsigned long       num_elements_surf;
element_surf_t      *element_surf;
unsigned int        *npe_surf;

unsigned long       num_nodes;
unsigned long       *nodes_per_process;
Vec                 coord;                  // coord[numumero total de nodos][3]

unsigned long       total_memory;

list_t              list_coord;
list_t              list_elemv;
list_t              list_elems;
list_t              list_physe;

char                mesh_file_name[WORDS_LENGHT];
char                input_file_name[MAXL];
char                epart_file_name[MAXL];
char                npart_file_name[MAXL];
char                xs_file_name[MAXL];

xstable_t           *xstable;
precursors_t        precursors;
double              *velocities;                      // Neutron medium velocity ( I groups )

calculation_flow_t  calcu;

/***************************************************************************/

// Tetrahedron
double xp_3_4[GPN_3_4][3];
double wp_3_4[GPN_3_4];
double shape_3_4[4][4];
double dshape_3_4[4][3];
double integral_3_4[2];
double integral_grad_3_4[4][4];

// Hexahedron
double xp_3_8[GPN_3_8][3];
double wp_3_8[GPN_3_8];

double shape_3_8[4][4];
double dshape_3_8[4][4];

#endif
