/*! \file examples/elasticity/elasticity_mixed.c
*
*  Created by James Adler on 7/4/20.
*  Copyright 2015_HAZMATH__. All rights reserved.
*
* \brief This program solves a linear elasticity problem using finite elements:
*
*       Primal: -div(2*mu*eps(u)) - grad(lam*div u) = f
*
*
*       where eps(u) = (grad u + (grad u)^T)/2 is the symmetric gradient,
*
* \note This example shows how to build your own bilinear form for a system.
*       The forms are found in elasticity_system.h and all Problem Data is found in
*       elasticity_data.h for different test problems.
*
*
*/

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "hazmath.h"
#include "elasticity_data.h"
#include "elasticity_system.h"
/*********************************************************************************/

/****** MAIN DRIVER **************************************************************/
int main (int argc, char* argv[])
{

  // Get Parameters
  REAL mu = 0.0;
  get_mu(&mu,NULL,0.0,NULL);
  REAL lam = 0.0;
  get_lam(&lam,NULL,0.0,NULL);

  printf("\n===========================================================================\n");
  printf("Beginning Program to solve Elasticity in Primal Form (mu = %e\tlam = %e)\n",mu,lam);
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
  INT order_u = 2; // Displacements

  // Need Spaces for each component of the velocity
  fespace FE_ux; // Displacement in x direction
  create_fespace(&FE_ux,&mesh,order_u);
  fespace FE_uy; // Displacement in y direction
  create_fespace(&FE_uy,&mesh,order_u);
  fespace FE_uz; // Displacement in z direction
  if(dim==3) create_fespace(&FE_uz,&mesh,order_u);

  // Set Dirichlet Boundaries
  /* If the mesh data has different flags for different boundaries, here
  * you can set different boundaries to have different conditions.
  * Alternatively, we can relabel the dof_flags of each finite-element space
  * to indicate different boundaries, then set various ones to have different
  * conditions.  A sample code in elasticity_data.h is given and can be called
  * via relabel_boundary2D(&FE_ux) for example.
  * For now, we just use the mesh boundaries and Dirichlet all around.
  */
  set_dirichlet_bdry(&FE_ux,&mesh,1,1);
  set_dirichlet_bdry(&FE_uy,&mesh,1,1);
  if(dim==3) set_dirichlet_bdry(&FE_uz,&mesh,1,1);

  // Create Block FE System with ordering (u,p)
  INT uxdof = FE_ux.ndof;
  INT uydof = FE_uy.ndof;
  INT uzdof;
  INT ndof = uxdof+uydof;
  INT nun = dim; // Number of scalar unknowns
  INT nspaces = dim; // Number of FE spaces
  if(dim==3) {
    uzdof = FE_uz.ndof;
    ndof += uzdof;
  }
  block_fespace FE;
  FE.nun = nun; // Number of unknowns (all scalar quantities)
  FE.ndof = ndof; // Total DoF
  FE.nbdof = FE_ux.nbdof + FE_uy.nbdof; // Total Boundary DoF
  if(dim==3) FE.nbdof += FE_uz.nbdof;
  FE.nspaces = nspaces; // Total number of FE spaces
  FE.var_spaces = (fespace **) calloc(nspaces,sizeof(fespace *));

  FE.var_spaces[0] = &FE_ux;
  FE.var_spaces[1] = &FE_uy;
  if(dim==3) FE.var_spaces[2] = &FE_uz;

  // Set Dirichlet Boundaries on the block system (grabs each spaces conditions)
  set_dirichlet_bdry_block(&FE,&mesh);

  clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
  printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
  (REAL) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
  /*******************************************************************************/

  printf("***********************************************************************************\n");
  printf("Number of Elements = %d\tOrder of Quadrature = %d\n",mesh.nelm,2*nq1d-1);
  printf("\n\t--- Element Type ---\n");
  printf("Displacement Element Type = %d\n",order_u);
  printf("\n\t--- Degrees of Freedom ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nv,mesh.nedge,mesh.nface);
  printf("\t--> DOF: %d\n",FE.ndof);
  printf("***********************************************************************************\n\n");

  /*** Assemble the matrix and right hand side *******************************/
  /* Here we assemble the discrete system in primal form:
  *  The weak form is:
  *
  *   2*mu*<eps(u), eps(v)> + lam*<div u, div v>   = <f, v>
  *
  * Note: We make sure the system is symmetric, and we have included the boundary
  *       integral in the case of a Dirichlet condition on p.
  */
  printf("Assembling the matrix and right-hand side:\n");
  clock_t clk_assembly_start = clock();

  // Declare the rhs vector and the stiffness matrix
  dvector b;
  block_dCSRmat A; // Block matrix for each FE space.
  bdcsr_alloc(nspaces,nspaces,&A);

  // Assemble the matrices without BC first
  if(dim==2) assemble_global_block(&A,&b,elasticity_primal_system,FEM_Block_RHS_Local,&FE,&mesh,cq,usource2D,0.0);
  if(dim==3) assemble_global_block(&A,&b,elasticity_primal_system,FEM_Block_RHS_Local,&FE,&mesh,cq,usource3D,0.0);

  // Allocate solution and set to be 0.0;
  dvector sol = dvec_create(ndof);
  dvec_set(sol.row, &sol, 0.0);

  // Eliminate remaining Dirichlet boundary conditions
  if(dim==2) eliminate_DirichletBC_blockFE_blockA(ubc2D,&FE,&mesh,&b,&A,0.0);
  if(dim==3) eliminate_DirichletBC_blockFE_blockA(ubc3D,&FE,&mesh,&b,&A,0.0);

  clock_t clk_assembly_end = clock();
  printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL)
  (clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  /************ Prepare Preconditioners **************************************************/
  // Prepare diagonal blocks
  dCSRmat *A_diag;
  A_diag = (dCSRmat *)calloc(nspaces, sizeof(dCSRmat));
  for(i=0;i<dim;i++){ // copy block diagonal to A_diag for displacements
    dcsr_alloc(A.blocks[i*(dim+1)]->row, A.blocks[i*(dim+1)]->col, A.blocks[i*(dim+1)]->nnz, &A_diag[i]);
    dcsr_cp(A.blocks[i*(dim+1)], &A_diag[i]);
  }
  /*******************************************************************************************/

  /***************** Solve *******************************************************************/
  printf("Solving the System:\n");
  clock_t clk_solve_start = clock();

  INT solver_flag = -20;

  // Set parameters for linear iterative methods
  linear_itsolver_param linear_itparam;
  param_linear_solver_init(&linear_itparam);
  param_linear_solver_set(&linear_itparam,&inparam);
  INT solver_type = linear_itparam.linear_itsolver_type;
  INT solver_printlevel = linear_itparam.linear_print_level;

  // Solve
  if(solver_type==0) { // Direct Solver
    solver_flag = block_directsolve_UMF(&A,&b,&sol,solver_printlevel);
  } else { // Iterative Solver
    if (linear_itparam.linear_precond_type == PREC_NULL) {
      solver_flag = linear_solver_bdcsr_krylov(&A, &b, &sol, &linear_itparam);
    } else {
      if(dim==2) solver_flag = linear_solver_bdcsr_krylov_block_2(&A, &b, &sol, &linear_itparam, NULL, A_diag);
      if(dim==3) solver_flag = linear_solver_bdcsr_krylov_block_3(&A, &b, &sol, &linear_itparam, NULL, A_diag);
    }
  }

  // Solver Error Check
  if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n",solver_flag);

  clock_t clk_solve_end = clock();
  printf("Elapsed CPU Time for Solve = %f seconds.\n\n",
  (REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  /********************* Compute Errors if you have exact solution ****************************/
  clock_t clk_error_start = clock();
  REAL* solerrL2 = (REAL *) calloc(nspaces, sizeof(REAL));
  REAL* solerrH1 = (REAL *) calloc(nspaces, sizeof(REAL)); // Note: Derivatives of P0 approximations are considered 0
  if(dim==2){
    L2error_block(solerrL2, sol.val, uexact_sol2D, &FE, &mesh, cq, 0.0);
    HDerror_block(solerrH1, sol.val, uexact_sol2D, uDexact_sol2D, &FE, &mesh, cq, 0.0);
  } else if(dim==3){
    L2error_block(solerrL2, sol.val, uexact_sol3D, &FE, &mesh, cq, 0.0);
    HDerror_block(solerrH1, sol.val, uexact_sol3D, uDexact_sol3D, &FE, &mesh, cq, 0.0);
  }

  REAL uerrL2 = 0;
  REAL uerrH1 = 0;
  for(i=0;i<dim;i++) uerrL2 += solerrL2[i]*solerrL2[i];
  for(i=0;i<dim;i++) uerrH1 += solerrH1[i]*solerrH1[i];
  uerrL2 = sqrt(uerrL2);
  uerrH1 = sqrt(uerrH1);

  printf("*******************************************************\n");
  printf("L2 Norm of u error    = %26.13e\n",uerrL2);
  printf("H1 Norm of u error    = %26.13e\n",uerrH1);
  printf("*******************************************************\n\n");
  clock_t clk_error_end = clock();
  printf("Elapsed CPU time for getting errors = %lf seconds.\n\n",(REAL)
  (clk_error_end-clk_error_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  // Plotting in vtu format
  char** varname;
  if(inparam.print_level > 3){
    char* soldump = "output/solution.vtu";
    varname = malloc(10*nspaces*sizeof(char *));
    varname[0] = "ux";
    varname[1] = "uy";
    if(dim==3) varname[2] = "uz";
    dump_blocksol_vtk(soldump,varname,&mesh,&FE,sol.val);
  }

  /************ Free All the Arrays ***********************************************************/
  // CSR
  bdcsr_free( &A );
  for(i=0;i<dim;i++) dcsr_free( &A_diag[i] );
  if(A_diag) free(A_diag);

  // Vectors
  if(solerrL2) free(solerrL2);
  if(solerrH1) free(solerrH1);
  dvec_free( &b );
  dvec_free( &sol );
  //dvec_free( &rhs_bdry );

  // FE Spaces
  free_fespace(&FE_ux);
  free_fespace(&FE_uy);
  if(dim==3) free_fespace(&FE_uz);
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
