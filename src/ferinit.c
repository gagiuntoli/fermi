#include "funct.h"

int ferinit(int argc,char **argv)
{
  /* Reads the imput file
   * Reads the mesh
   * Reads the mesh
   * Allocs mamory for K, x, b
   */
  int error,i,d,*ghost;
  SlepcInitialize(&argc,&argv,(char*)0,NULL);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &nproc);  
  calcu.exec = (nproc>1)?PARALLEL:SEQUENCIAL;
  if(argc == 1)
  {
    PetscPrintf(PETSC_COMM_WORLD,"main.c:input file NS.\n\n"); 
    return 1;
  }
  //  
  //============================== 
  // PARCING INPUT FILE
  //============================== 
  //    
  list_init(&list_mater, sizeof(pvl_t),cmp_mat);
  list_init(&list_bound, sizeof(bound_t),cmp_bou);
  list_init(&list_fun1d, sizeof(bound_t),cmp_f1d);
  list_init(&list_ctrlr, sizeof(ctrlrod_t),NULL);
  strcpy(inputfile,argv[1]);
  PetscPrintf(PETSC_COMM_WORLD,"Parcing input file.\n");
  error=parse_input();
  if(error!=0)
  {
    PetscPrintf(PETSC_COMM_WORLD,"main.c:error parsing input file.\n");
    return 1;
  }
  //
  //============================== 
  // READING MESH 
  //============================== 
  //   
  list_init(&list_nodes, sizeof(gmshN_t), gmsh_nodcmp);
  list_init(&list_ghost, sizeof(gmshN_t), gmsh_nodcmp);
  list_init(&list_elemv, sizeof(gmshE_t), cmp_int);
  list_init(&list_elems, sizeof(gmshE_t), cmp_int);
  list_init(&list_physe, sizeof(gmshP_t), cmp_int);    
  PetscPrintf(PETSC_COMM_WORLD,"Reading mesh.\n");
  error=gmsh_read(meshfile,epartfile,npartfile,rank,DIM,&list_nodes,&list_ghost,&list_elemv,&list_elems,&list_physe,&loc2gold,&loc2gnew,&npp,nproc);    
  if(error!=0)
  {
    PetscPrintf(PETSC_COMM_WORLD,"main.c:error reading mesh.\n"); 
    return 1;
  }
  ntot=0;
  for(i=0;i<nproc;i++)
    ntot+=npp[i];
  //
  //============================== 
  // PRINTING STRUCTURES
  //============================== 
  //     
  PetscPrintf(PETSC_COMM_WORLD,"Printing structures 1.\n");
  error = print_struct(1);   
  if(error!=0){
    PetscPrintf(PETSC_COMM_WORLD,"main.c:error printing structures.\n"); 
    return 1;
  }

  //
  //============================== 
  // CONSTRUCTING MESH
  //============================== 
  //      
  PetscPrintf(PETSC_COMM_WORLD,"Constructing mesh.\n");
  error=mesh_alloc(&list_nodes, &list_ghost, cpynode, &list_elemv, cpyelemv, &list_elems, cpyelems, &mesh);
  if(error) 
  {
    PetscPrintf(PETSC_COMM_WORLD,"main.c:error allocating mesh.\n"); 
    return 1;
  }
  error=mesh_renum(&mesh,loc2gold,loc2gnew);
  if(error)
  {
    PetscPrintf(PETSC_COMM_WORLD,"main.c:error renumbering mesh nodes.\n"); 
    return 1;
  }
//  error=mesh_neigh(&mesh,loc2gnew);
//  if(error)
//  {
//    PetscPrintf(PETSC_COMM_WORLD,"main.c: Error calculating mesh neightbors.\n",error); 
//    return 1;
//  }

  //      
  //==============================      
  // CONTROL RODS INIT
  //==============================      
  //      
  PetscPrintf(PETSC_COMM_WORLD,"Initializing control rods.\n");
  error=ferirods();
  if(error)
  {
    PetscPrintf(PETSC_COMM_WORLD,"main.c:error control rods initialization.\n",error); 
    return 1;
  }
  //      
  //==============================      
  // ASSEMBLY BOUNDARIES 
  //==============================      
  //      
  PetscPrintf(PETSC_COMM_WORLD,"Assemblying BCs.\n");
  error=ferbouset();
  if(error)
  {
    PetscPrintf(PETSC_COMM_WORLD,"main.c:error assembling BCs.\n"); 
    return 1;
  }
  //      
  //==============================      
  // ALLOCATING MATRICES/VECTORS 
  //==============================      
  //      
  PetscPrintf(PETSC_COMM_WORLD,"Allocating Matrices/Vectors.\n");
  ghost=(int*)calloc(mesh.nghost*DIM,sizeof(int));
  for(i=0;i<mesh.nghost;i++)
  {
    for(d=0;d<DIM;d++)
      ghost[i*egn+d]=loc2gnew[mesh.nnodes+i]*egn+d;
  }

  VecCreateGhost(PETSC_COMM_WORLD,mesh.nnodes*egn,ntot*egn,mesh.nghost*egn,(PetscInt*)ghost,&phi_n); 
  VecDuplicate(phi_n,&phi_o);
  VecDuplicate(phi_n,&b);
  VecDuplicate(phi_n,&b_a);
  MatCreateAIJ(PETSC_COMM_WORLD,mesh.nnodes*egn,mesh.nnodes*egn,ntot*egn,ntot*egn,78,NULL,78,NULL,&A);
  MatCreateAIJ(PETSC_COMM_WORLD,mesh.nnodes*egn,mesh.nnodes*egn,ntot*egn,ntot*egn,78,NULL,78,NULL,&B);
  MatCreateAIJ(PETSC_COMM_WORLD,mesh.nnodes*egn,mesh.nnodes*egn,ntot*egn,ntot*egn,78,NULL,78,NULL,&M);
  MatCreateAIJ(PETSC_COMM_WORLD,mesh.nnodes*egn,mesh.nnodes*egn,ntot*egn,ntot*egn,78,NULL,78,NULL,&K);

  Ae=(double*)calloc(NPE*egn*NPE*egn,sizeof(double));
  Be=(double*)calloc(NPE*egn*NPE*egn,sizeof(double));
  Me=(double*)calloc(NPE*egn*NPE*egn,sizeof(double));
  be=(double*)calloc(NPE*egn,sizeof(double));
  idxm=(int*)calloc(NPE*egn,sizeof(int));
  jac=(double**)calloc(DIM,sizeof(double*));
  for(i=0;i<DIM;i++)
    jac[i]=(double*)calloc(DIM,sizeof(double));
  ijac=(double**)calloc(DIM,sizeof(double*));
  for(i=0;i<DIM;i++)
    ijac[i]=(double*)calloc(DIM,sizeof(double));
  der=(double**)calloc(NPE,sizeof(double*));
  for(i=0;i<NPE;i++)
    der[i]=(double*)calloc(DIM,sizeof(double));
  coor=(double**)calloc(NPE,sizeof(double*));
  for(i=0;i<NPE;i++)
    coor[i]=(double*)calloc(DIM,sizeof(double));

  fem_inigau();
  if(error)
  {
    PetscPrintf(PETSC_COMM_WORLD,"main.c:error gps init.\n\n"); 
    return 1;
  }

  //      
  //==============================      
  // SETTING SOLVER
  //==============================      
  //      
  EPSCreate(PETSC_COMM_WORLD,&eps);
  EPSSetProblemType(eps,EPS_GNHEP);    
  EPSSetOperators(eps,A,B);
  EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);
  EPSSetType(eps,EPSJD);
  EPSSetDimensions(eps,1,PETSC_DEFAULT,PETSC_DEFAULT);
  EPSSetFromOptions(eps);

  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetType(ksp,KSPGMRES);
  KSPGetPC(ksp,&pc);
//  PCSetType(pc,PCLU);
  KSPSetFromOptions(ksp);
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"kspinfo.dat",&kspview);
  KSPView(ksp,kspview);

 // EPSCreate(PETSC_COMM_WORLD,&eps);
 // EPSSetProblemType(eps,EPS_GNHEP);

  //      
  //============================== 
  // PRINTING STRUCTURES
  //============================== 
  //     
  PetscPrintf(PETSC_COMM_WORLD,"Printing structures 2.\n");
  error = print_struct(2);   
  if(error)
  {
    PetscPrintf(PETSC_COMM_WORLD,"main.c:error printing structures.\n"); 
    return 1;
  }
  return 0;
}
