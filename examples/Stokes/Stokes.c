/*! \file examples/Stokes/Stokes.c
 *
 *  Created by Peter Ohm on 2/5/17.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program solves Stokes PDE using finite elements
 *
 *      -laplace(u) + div(p) = f
 *                   grad(u) = 0
 *
 *
 *        in 2D or 3D
 *
 *        Along the boundary of the region, Dirichlet conditions are
 *        imposed for u and Neumann for p.  P2-P1 or P2-P0 can be used.
 *
 * \note This example shows how to build your own bilinear form for a system.
 *       The forms are found in StokesSystem.h and all Problem Data is found in
 *       StokesData.h
 *
 *
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "global.h"
#include "hazmath.h"
#include "StokesData.h"
#include "StokesSystem.h"
/*********************************************************************************/

/****** MAIN DRIVER **************************************************************/
int main (int argc, char* argv[])
{

  printf("\n===========================================================================\n");
  printf("Beginning Program to solve Stokes Equation.\n");
  printf("===========================================================================\n");

  /****** INITIALIZE PARAMETERS **************************************************/
  INT i;

  // Overall CPU Timing
  clock_t clk_overall_start = clock();

  // Set Parameters from Reading in Input File
  input_param inparam;
  param_input_init(&inparam);
  param_input("./input.dat", &inparam);

  // Open gridfile for reading
  printf("\nCreating mesh and FEM spaces:\n");
  FILE* gfid = HAZ_fopen(inparam.gridfile,"r");

  // Create the mesh (now we assume triangles in 2D or tetrahedra in 3D)
  // File types possible are 0 - old format; 1 - vtk format
  INT mesh_type = 0;
  clock_t clk_mesh_start = clock(); // Time mesh generation FE setup
  trimesh mesh;
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
  INT order_u = 2;
  INT order_p = 0;

  // Need Spaces for each component of the velocity plus pressure
  fespace FE_ux; // Velocity in x direction
  create_fespace(&FE_ux,&mesh,order_u);
  fespace FE_uy; // Velocity in y direction
  create_fespace(&FE_uy,&mesh,order_u);
  fespace FE_uz; // Velocity in z direction
  if(dim==3) create_fespace(&FE_uz,&mesh,order_u);
  fespace FE_p; // Pressure
  create_fespace(&FE_p,&mesh,order_p);

  // Set Dirichlet Boundaries
  set_dirichlet_bdry(&FE_ux,&mesh,1,1);
  set_dirichlet_bdry(&FE_uy,&mesh,1,1);
  if(dim==3) set_dirichlet_bdry(&FE_uz,&mesh,1,1);
  set_dirichlet_bdry(&FE_p,&mesh,1,1);
  for(i=0;i<FE_p.ndof;i++) {
    FE_p.dirichlet[i] = 0;
  }

  // Create Block System with ordering (u,p)
  INT ndof = FE_ux.ndof + FE_uy.ndof + FE_p.ndof;
  if(dim==3) ndof += FE_uz.ndof;
  // Get Global FE Space
  block_fespace FE;
  FE.nun = dim+1;
  FE.ndof = ndof;
  FE.nbdof = FE_ux.nbdof + FE_uy.nbdof + FE_p.nbdof;
  if(dim==3) FE.nbdof += FE_uz.nbdof;
  FE.nspaces = dim+1;
  FE.var_spaces = (fespace **) calloc(dim+1,sizeof(fespace *));

  FE.var_spaces[0] = &FE_ux;
  FE.var_spaces[1] = &FE_uy;
  if(dim==3) FE.var_spaces[2] = &FE_uz;
  FE.var_spaces[dim] = &FE_p;
  
  // Set Dirichlet Boundaries
  set_dirichlet_bdry_block(&FE,&mesh);


  clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
  printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
         (REAL) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
  /*******************************************************************************/

  printf("***********************************************************************************\n");
  printf("Number of Elements = %d\tOrder of Quadrature = %d\n",mesh.nelm,2*nq1d-1);
  printf("\n\t--- Element Type ---\n");
  printf("Velocity Element Type = %d\tPressure Element Type = %d\n",order_u,order_p);
  printf("\n\t--- Degrees of Freedom ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nv,mesh.nedge,mesh.nface);
  printf("\t--> DOF: %d\n",FE.ndof);
  printf("\n\t--- Boundaries ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nbv,mesh.nbedge,mesh.nbface);
  printf("\t--> Boundary DOF: %d\n",FE.nbdof);
  printf("***********************************************************************************\n\n");

  /*** Assemble the matrix and right hand side *******************************/
  /* Here we assemble the discrete system:
   *  The weak form is:
   *
   *  <grad u, grad v> - <p, div v> = <f, v>
   *                   - <div u, q> = 0
   */
  printf("Assembling the matrix and right-hand side:\n");
  clock_t clk_assembly_start = clock();

  // Allocate the right-hand and declare the csr matrix
  dvector b;

  // Put into Block Form
  block_dCSRmat A;
  bdcsr_alloc(dim+1,dim+1,&A);

  // Assemble the matricies without BC first
  if(dim==2) assemble_global_block(&A,&b,local_assembly_Stokes,FEM_Block_RHS_Local,&FE,&mesh,cq,source2D,0.0);
  if(dim==3) assemble_global_block(&A,&b,local_assembly_Stokes,FEM_Block_RHS_Local,&FE,&mesh,cq,source3D,0.0);

  // Eliminate boundary conditions in matrix and rhs
  if(dim==2) {
    eliminate_DirichletBC_blockFE_blockA(bc2D,&FE,&mesh,&b,&A,0.0);
  }
  if(dim==3) {
    eliminate_DirichletBC_blockFE_blockA(bc3D,&FE,&mesh,&b,&A,0.0);
  }

  /**************************************************/
  //  Apply Pressure "BCs" (removes singularity)
  REAL pressureval =0.5;
  INT pressureloc = 0;
  /**************************************************/

  clock_t clk_assembly_end = clock();
  printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL)
         (clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/
  
  /************ Prepare Preconditioners **************************************************/
  
  // Shift for HAZMATH
  for(i=0;i<(FE.nspaces)*(FE.nspaces);i++) {
    if(A.blocks[i])
      dcsr_shift(A.blocks[i],-1);
  }

  // Prepare diagonal blocks
  dCSRmat *A_diag;
  A_diag = (dCSRmat *)calloc(dim+1, sizeof(dCSRmat));

  for(i=0;i<dim;i++){ // copy block diagonal to A_diag
    dcsr_alloc(A.blocks[i*(dim+2)]->row, A.blocks[i*(dim+2)]->col, A.blocks[i*(dim+2)]->nnz, &A_diag[i]);
    dcsr_cp(A.blocks[i*(dim+2)], &A_diag[i]);
  }

  // Get Mass Matrix for p
  dCSRmat Mp;
  assemble_global(&Mp,NULL,assemble_mass_local,&FE_p,&mesh,cq,NULL,one_coeff_scal,0.0);
  dcsr_shift(&Mp,-1);
  dcsr_alloc(Mp.row, Mp.col, Mp.nnz, &A_diag[dim]);
  dcsr_cp(&Mp, &A_diag[dim]);
  /*******************************************************************************************/

  /***************** Solve *******************************************************************/
  printf("Solving the System:\n");
  clock_t clk_solve_start = clock();

  INT solver_flag = -20;

  // Allocate solution
  dvector sol = dvec_create(ndof);
  dvector v_ux = dvec_create(FE_ux.ndof);
  dvector v_uy = dvec_create(FE_uy.ndof);
  dvector v_uz;
  if(dim==3) v_uz = dvec_create(FE_uz.ndof);
  dvector v_p = dvec_create(FE_p.ndof);
  
  // Set initial guess to be all zero
  dvec_set(sol.row, &sol, 0.0);
  // Set initial guess for pressure to match the known "boundary condition" for pressure
  if(dim==2) sol.val[FE_ux.ndof + FE_uy.ndof + pressureloc]  = pressureval;
  if(dim==3) sol.val[FE_ux.ndof + FE_uy.ndof + FE_uz.ndof + pressureloc]  = pressureval;

  // Set parameters for linear iterative methods
  linear_itsolver_param linear_itparam;
  param_linear_solver_init(&linear_itparam);
  param_linear_solver_set(&linear_itparam,&inparam);

  // Solve
  if(dim==2){
    if (linear_itparam.linear_precond_type == PREC_NULL) {
      solver_flag = linear_solver_bdcsr_krylov(&A, &b, &sol, &linear_itparam);
    } else {
      solver_flag = linear_solver_bdcsr_krylov_block_3(&A, &b, &sol, &linear_itparam, NULL, A_diag);
    }
  } else if (dim==3) {
    if (linear_itparam.linear_precond_type == PREC_NULL) {
      solver_flag = linear_solver_bdcsr_krylov(&A, &b, &sol, &linear_itparam);
    } else {
      solver_flag = linear_solver_bdcsr_krylov_block_4(&A, &b, &sol, &linear_itparam, NULL, A_diag);
    }
  }

  // Error Check
  if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n",solver_flag);
  
  clock_t clk_solve_end = clock();
  printf("Elapsed CPU Time for Solve = %f seconds.\n\n",
         (REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/
  
  /********************* Compute Errors if you have exact solution ****************************/
  clock_t clk_error_start = clock();
  REAL* solerrL2 = (REAL *) calloc(dim+1, sizeof(REAL));
  REAL* solerrH1 = (REAL *) calloc(dim+1, sizeof(REAL)); // Note: No H1 error for P0 elements
  if(dim==2){
    L2error_block(solerrL2, sol.val, exact_sol2D, &FE, &mesh, cq, 0.0);
    if(order_p > 0) HDerror_block(solerrH1, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
  } else if(dim==3){
    L2error_block(solerrL2, sol.val, exact_sol3D, &FE, &mesh, cq, 0.0);
    if(order_p > 0) HDerror_block(solerrH1, sol.val, exact_sol3D, Dexact_sol3D, &FE, &mesh, cq, 0.0);
  }

  REAL uerrL2 = 0;
  REAL uerrH1 = 0;
  for(i=0;i<dim;i++) uerrL2 += solerrL2[i]*solerrL2[i];
  for(i=0;i<dim;i++) uerrH1 += solerrH1[i]*solerrH1[i];
  uerrL2 = sqrt(uerrL2);
  uerrH1 = sqrt(uerrH1);
  REAL perrL2 = solerrL2[dim];
  REAL perrH1 = solerrH1[dim];

  printf("*******************************************************\n");
  printf("L2 Norm of u error    = %26.13e\n",uerrL2);
  printf("L2 Norm of p error    = %26.13e\n",perrL2);
  printf("H1 Norm of u error    = %26.13e\n",uerrH1);
  printf("H1 Norm of p error    = %26.13e\n",perrH1);
  printf("*******************************************************\n\n");
  clock_t clk_error_end = clock();
  printf("Elapsed CPU time for getting errors = %lf seconds.\n\n",(REAL)
         (clk_error_end-clk_error_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  // Plotting
  get_unknown_component(&v_ux,&sol,&FE,0);
  get_unknown_component(&v_uy,&sol,&FE,1);
  if(dim==3) get_unknown_component(&v_uz,&sol,&FE,2);
  get_unknown_component(&v_p,&sol,&FE,dim);
 
  char** varname; 
  if(inparam.print_level > 3){
    char* soldump = "output/solution.vtu";
    varname = malloc(10*FE.nspaces*sizeof(char *));
    varname[0] = "ux";
    varname[1] = "uy";
    if(dim==3) varname[2] = "uz";
    varname[dim] = "p ";
    dump_blocksol_vtk(soldump,varname,&mesh,&FE,sol.val);

    if(dim==3) print_matlab_vector_field(&v_ux,&v_uy,&v_uz,&FE_ux);
  }
  /************ Free All the Arrays ***********************************************************/
  // CSR
  bdcsr_free( &A );
  dcsr_free( &Mp);
  for(i=0;i<dim+1;i++)
    dcsr_free( &A_diag[i] );
  if(A_diag) free(A_diag);

  // Vectors
  if(solerrL2) free(solerrL2);
  if(solerrH1) free(solerrH1);
  dvec_free( &b );
  dvec_free( &sol );
  dvec_free( &v_ux );
  dvec_free( &v_uy );
  if(dim==3) dvec_free( &v_uz );
  dvec_free( &v_p );

  // FE Spaces
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

  // Strings
  if(varname) free(varname);  

  /*******************************************************************************************/
  clock_t clk_overall_end = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",
         (REAL) (clk_overall_end-clk_overall_start)/CLOCKS_PER_SEC);
  return 0;
}  /* End of Program */
/*********************************************************************************************/
