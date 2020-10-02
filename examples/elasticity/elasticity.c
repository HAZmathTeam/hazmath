/*! \file examples/elasticity/elasticity.c
*
*  Created by Peter Ohm on 2/5/20.
*  Copyright 2015_HAZMATH__. All rights reserved.
*
* \brief This program solves Elasticity using finite elements
*
*      -2*div(mu*eps(u)+lam*div(u)) + grad(p)    = f
*             div(u)               = 0
*
*        where eps(u) = (grad u + (grad u)^T)/2 is the symmetric gradient,
*
*        in 2D or 3D.
*
*        Along the boundary of the region, Dirichlet conditions are
*        imposed for u and Neumann for p.  P2-P1 or P2-P0 can be used,
*        though others can be implemented.
*
*
*
*/

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "hazmath.h"
/*********************************************************************************/

/******** Data Input *************************************************/
// PDE variables

// Youngs Modulus
void get_young(REAL *val,REAL* x,REAL time,void *param) {
  *val = 3e4;
}
// Poisson Ratio
void get_nu(REAL *val,REAL* x,REAL time,void *param) {
  *val = 0.499999999;
}

// Lame Coefficients
void get_lam(REAL *val,REAL* x,REAL time,void *param) {
  REAL E;
  REAL nu;
  get_nu(&nu,x,time,param);
  get_young(&E,x,time,param);
  *val = E*nu/((1-2*nu)*(1+nu));
}
void  get_mu(REAL *val,REAL* x,REAL time,void *param) {
  REAL E;
  REAL nu;
  get_nu(&nu,x,time,param);
  get_young(&E,x,time,param);
  *val = E/(1+2*nu);
}

// Right-hand side
void source2D(REAL *val,REAL* x,REAL time,void *param) {
  val[0] = 0.0;// ux
  val[1] = 0.0;// uy
  val[2] = 0.0;// p
}

// Boundary Conditions
void bc2D(REAL *val,REAL* x,REAL time,void *param) {
  val[0] = 0.0;// ux
  val[1] = 0.0;// uy
  val[2] = 0.0;// p
}

/************* Bad Programming ***************************************************/
#include "elasticity_system.h"
/*********************************************************************************/

/****** MAIN DRIVER **************************************************************/
int main (int argc, char* argv[])
{

  printf("\n===========================================================================\n");
  printf("Beginning Program to solve Stokes Equation.\n");
  printf("===========================================================================\n");

  /****** INITIALIZE PARAMETERS **************************************************/
  // Loop Indices
  INT i, j;
  bool SOLVE_GMG = true;

  // Overall CPU Timing
  clock_t clk_overall_start = clock();

  // Set Parameters from Reading in Input File
  input_param inparam;
  param_input_init(&inparam);
  param_input("./input.dat", &inparam);

  // Open gridfile for reading
  printf("\nCreating mesh and FEM spaces:\n");
  FILE* gfid = HAZ_fopen(inparam.gridfile,"r");

  i = 0;
  char GRIDDIR[200];// THIS GETS THE LOCATION OF THE GRID DIRECTORY
  while( inparam.gridfile[i] != '\0'){
    GRIDDIR[i] = inparam.gridfile[i];
    if(inparam.gridfile[i] == '\\' || inparam.gridfile[i] == '/'){ j=0;} else { j++;}
    i++;
  }
  GRIDDIR[i-j] = '\0';
  printf("Path to mesh directory: %s\n",GRIDDIR);

  // Create the mesh (now we assume triangles in 2D or tetrahedra in 3D)
  // File types possible are 0 - old format; 1 - vtk format
  INT mesh_type = 0;
  clock_t clk_mesh_start = clock(); // Time mesh generation FE setup
  mesh_struct mesh;
  printf(" --> loading grid from file: %s\n",inparam.gridfile);
  initialize_mesh(&mesh);
  creategrid_fread(gfid,mesh_type,&mesh);
  fclose(gfid);
  INT dim = mesh.dim;

  // Get Quadrature Nodes for the Mesh
  INT nq1d = inparam.nquad; /* Quadrature points per dimension */
  qcoordinates *cq = get_quadrature(&mesh,nq1d);

  // Get info for and create FEM spaces
  // Order of elements: 0 - P0; 1 - P1; 2 - P2; 20 - Nedlec; 30 - Raviart-Thomas
  INT order_b = 61;
  INT order_u = 1;
  INT order_p = 0;

  // Need Spaces for each component of the velocity plus pressure
  fespace FE_bb; // Displacement bubble
  create_fespace(&FE_bb,&mesh,order_b);
  fespace FE_ux; // Displacement in x direction
  create_fespace(&FE_ux,&mesh,order_u);
  fespace FE_uy; // Displacement in y direction
  create_fespace(&FE_uy,&mesh,order_u);
  fespace FE_uz; // Displacement in z direction
  if(dim==3) create_fespace(&FE_uz,&mesh,order_u);
  fespace FE_p;  // Pressure
  create_fespace(&FE_p,&mesh,order_p);

  // Set Dirichlet Boundaries
  set_dirichlet_bdry(&FE_bb,&mesh,1,1);
  set_dirichlet_bdry(&FE_ux,&mesh,1,1);
  set_dirichlet_bdry(&FE_uy,&mesh,1,1);
  if(dim==3) set_dirichlet_bdry(&FE_uz,&mesh,1,1);
  set_dirichlet_bdry(&FE_p,&mesh,-1,-1);

  // Create Block System with ordering (u,p)
  INT ndof = FE_bb.ndof + FE_ux.ndof + FE_uy.ndof + FE_p.ndof;
  if(dim==3) ndof += FE_uz.ndof;
  // Get Global FE Space
  block_fespace FE;
  FE.nun = dim+1;
  FE.ndof = ndof;
  FE.nbdof = FE_bb.nbdof + FE_ux.nbdof + FE_uy.nbdof + FE_p.nbdof;
  if(dim==3) FE.nbdof += FE_uz.nbdof;
  FE.nspaces = dim+2;
  FE.var_spaces = (fespace **) calloc( FE.nspaces ,sizeof(fespace *));
  FE.var_spaces[0] = &FE_bb;
  FE.var_spaces[1] = &FE_ux;
  FE.var_spaces[2] = &FE_uy;
  if(dim==3) FE.var_spaces[3] = &FE_uz;
  FE.var_spaces[dim+1] = &FE_p;

  // Set Dirichlet Boundaries
  set_dirichlet_bdry_block(&FE,&mesh);


  clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
  printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",(REAL) (clk_mesh_end-clk_mesh_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  printf("***********************************************************************************\n");
  printf("Number of Elements = %d\tOrder of Quadrature = %d\n",mesh.nelm,2*nq1d-1);
  printf("\n\t--- Degrees of Freedom ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d", mesh.nv, mesh.nedge, mesh.nface);
  printf("\t--> DOF: %d\n",FE.ndof);
  printf("\n\t--- Boundaries ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nbv,mesh.nbedge,mesh.nbface);
  printf("\t--> Boundary DOF: %d\n",FE.nbdof);
  printf("***********************************************************************************\n\n");

  /******************** Assemble the matrix and right hand side ******************************/
  /* Here we assemble the discrete system:
  *  The weak form is:
  *
  *  2*mu*<eps(u), eps(v)> + lam*<div u, div v> - <p, div v> = <f, v>
  *  - <div u, q>                                            = 0
  */
  printf("Assembling the matrix and right-hand side:\n");
  clock_t clk_assembly_start = clock();

  // Allocate the right-hand and declare the csr matrix
  dvector b;

  // Put into Block Form
  block_dCSRmat A;
  bdcsr_alloc(dim+2,dim+2,&A);

  // Assemble the matricies without BC first
  if(dim==2) assemble_global_block(&A,&b,Elasticity_system,FEM_Block_RHS_Local,&FE,&mesh,cq,source2D,0.0);

  // Merge Displacement into single block without BC
  block_dCSRmat A2_noBC;
  bdcsr_alloc(2, 2, &A2_noBC);
  dCSRmat Ablk11_noBC = bdcsr_subblk_2_dcsr(&A,     0,   dim,     0,   dim);
  dCSRmat Ablk12_noBC = bdcsr_subblk_2_dcsr(&A,     0,   dim, dim+1, dim+1);
  dCSRmat Ablk21_noBC = bdcsr_subblk_2_dcsr(&A, dim+1, dim+1,     0,   dim);
  dCSRmat Ablk22_noBC = bdcsr_subblk_2_dcsr(&A, dim+1, dim+1, dim+1, dim+1);
  A2_noBC.blocks[0] = &Ablk11_noBC;
  A2_noBC.blocks[1] = &Ablk12_noBC;
  A2_noBC.blocks[2] = &Ablk21_noBC;
  A2_noBC.blocks[3] = &Ablk22_noBC;

  // Eliminate boundary conditions in matrix and rhs
  if(dim==2) {
    eliminate_DirichletBC_blockFE_blockA(bc2D,&FE,&mesh,&b,&A,0.0);
  }

  // Merge Displacement into single block
  block_dCSRmat A2;
  bdcsr_alloc(2, 2, &A2);
  dCSRmat Ablk11 = bdcsr_subblk_2_dcsr(&A,     0,   dim,     0,   dim);
  dCSRmat Ablk12 = bdcsr_subblk_2_dcsr(&A,     0,   dim, dim+1, dim+1);
  dCSRmat Ablk21 = bdcsr_subblk_2_dcsr(&A, dim+1, dim+1,     0,   dim);
  dCSRmat Ablk22 = bdcsr_subblk_2_dcsr(&A, dim+1, dim+1, dim+1, dim+1);
  A2.blocks[0] = &Ablk11;
  A2.blocks[1] = &Ablk12;
  A2.blocks[2] = &Ablk21;
  A2.blocks[3] = &Ablk22;

  clock_t clk_assembly_end = clock();
  printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL) (clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  /**************** Prepare Preconditioners **************************************************/

  // Prepare diagonal blocks
  dCSRmat *A_diag;
  A_diag = (dCSRmat *)calloc(2, sizeof(dCSRmat));

  for(i=0;i<1;i++){ // copy block diagonal to A_diag
    dcsr_alloc(A2.blocks[i*(dim+2)]->row, A2.blocks[i*(dim+2)]->col, A2.blocks[i*(dim+2)]->nnz, &A_diag[i]);
    dcsr_cp(A2.blocks[i*(dim+2)], &A_diag[i]);
  }

  // Get Mass Matrix for p
  dCSRmat Mp;
  assemble_global(&Mp,NULL,assemble_mass_local,&FE_p,&mesh,cq,NULL,one_coeff_scal,0.0);
  dcsr_alloc(Mp.row, Mp.col, Mp.nnz, &A_diag[1]);
  dcsr_cp(&Mp, &A_diag[1]);
  /*******************************************************************************************/

  /***************** Solve *******************************************************************/
  printf("Solving the System:\n");
  clock_t clk_solve_start = clock();

  INT solver_flag = -20;

  // Allocate solution
  dvector sol = dvec_create(ndof);

  // Set initial guess to be all zero
  dvec_set(sol.row, &sol, 0.0);
  // SET SOLN and RANDOM INIT
  dvec_set(b.row,&b,0.0);// solve zero
  dvec_rand(sol.row,&sol);
  for(i=0; i<sol.row; i++){ if( FE.dirichlet[i] == 1){ sol.val[i] = 0.0; } }

  // Set parameters for linear iterative methods
  linear_itsolver_param linear_itparam;
  param_linear_solver_init(&linear_itparam);
  param_linear_solver_set(&linear_itparam,&inparam);

  // Set parameters for AMG
  AMG_param amgparam;
  param_amg_init(&amgparam);
  param_amg_set(&amgparam, &inparam);

  // Solve
  if( SOLVE_GMG ){
/*====================================================================================================*/
    solve_stats solve_info;
    solve_info.iteration_count = 0;
    solve_info.time_setup = 0.0;
    solve_info.time_precondition_setup = 0.0;
    solve_info.time_solve = 0.0;
    // POINTERS FOR STUFF
    INT gmg_type[]          = {999,0};
    INT Schwarz_on_blk[]    = {1,0};
    amgparam.Schwarz_on_blk = Schwarz_on_blk;
    amgparam.AMG_type       = -1;
    amgparam.HAZDIR         = GRIDDIR;

    FILE* fid;
    fid = fopen("A_matrix.dat","w");
    bdcsr_print_matlab(fid,&A2);
    fclose(fid);


    solver_flag = linear_solver_bdcsr_gmg(&A2,&b,&sol,&amgparam,gmg_type,&mesh,&FE,NULL,A_diag,&A2_noBC,&linear_itparam,&solve_info);
    //dvector_print(stdout, &sol);
/*====================================================================================================*/
  } else {
    if(dim==2){
      if (linear_itparam.linear_precond_type == PREC_NULL) {
        solver_flag = linear_solver_bdcsr_krylov(&A2, &b, &sol, &linear_itparam);
      } else {
        solver_flag = linear_solver_bdcsr_krylov_block_2(&A2, &b, &sol, &linear_itparam, NULL, A_diag);
      }
    } else if (dim==3) {
      if (linear_itparam.linear_precond_type == PREC_NULL) {
        solver_flag = linear_solver_bdcsr_krylov(&A, &b, &sol, &linear_itparam);
      } else {
        solver_flag = linear_solver_bdcsr_krylov_block_4(&A, &b, &sol, &linear_itparam, NULL, A_diag);
      }
    }
  }

  // Error Check
  if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n",solver_flag);

  clock_t clk_solve_end = clock();
  printf("Elapsed CPU Time for Solve = %f seconds.\n\n",
  (REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  /********************* Compute Errors if you have exact solution ***************************/

  /*******************************************************************************************/

  /************ Free All the Arrays **********************************************************/
  // CSR
  bdcsr_free( &A );
  dcsr_free( &Mp);
  for(i=0;i<2;i++)
  dcsr_free( &A_diag[i] );
  if(A_diag) free(A_diag);

  // Vectors
  dvec_free( &b );
  dvec_free( &sol );

  // FE Spaces
  free_fespace(&FE_bb);
  free_fespace(&FE_ux);
  free_fespace(&FE_uy);
  if(dim==3) free_fespace(&FE_uz);
  free_fespace(&FE_p);
  free_blockfespace(&FE);

  // Quadrature
  if(cq){
    free_qcoords(cq);
    free(cq);
    cq=NULL;
  }

  // Mesh
  free_mesh(&mesh);

  /*******************************************************************************************/
  clock_t clk_overall_end = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",
  (REAL) (clk_overall_end-clk_overall_start)/CLOCKS_PER_SEC);
  return 0;
}  /* End of Program */
/*********************************************************************************************/
