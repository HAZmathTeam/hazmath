/*! \file examples/stokes/stokes.c
 *
 *  Created by Peter Ohm on 2/5/17.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program solves Stokes PDE using finite elements
 *
 *      -2*div(eps(u)) + grad(p) = f
 *                 div(u)           = 0
 *
 *        where eps(u) = (grad u + (grad u)^T)/2 is the symmetric gradient,
 *
 *        in 2D or 3D.
 *
 *        Along the boundary of the region, Dirichlet conditions are
 *        imposed for u and Neumann for p.  P2-P1 or P2-P0 can be used,
 *        though others can be implemented.
 *
 * \note This example shows how to build your own bilinear form for a system.
 *       The forms are found in stokes_system.h and all Problem Data is found in
 *       stokes_data.h to simplify the code.  This example also illustrates how to
 *       construct block versions of the finite-element spaces, linear systems, and
 *       how to use block solvers.
 *
 *
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "hazmath.h"
#include "problem_data.h"
#include "variational_system.h"
/*********************************************************************************/

/****** MAIN DRIVER **************************************************************/
int main (int argc, char* argv[])
{

  printf("\n===========================================================================\n");
  printf("Beginning Program to solve Symmetric Gradient System.\n");
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
  INT order_u = 1;

  // Need Spaces for each component of the velocity plus pressure
  fespace FE_ux; // Velocity in x direction
  create_fespace(&FE_ux,&mesh,order_u);
  fespace FE_uy; // Velocity in y direction
  create_fespace(&FE_uy,&mesh,order_u);
  fespace FE_uz; // Velocity in z direction
  if(dim==3) create_fespace(&FE_uz,&mesh,order_u);

  // Set Dirichlet Boundaries
  set_dirichlet_bdry(&FE_ux,&mesh,1,1);
  set_dirichlet_bdry(&FE_uy,&mesh,1,1);
  if(dim==3) set_dirichlet_bdry(&FE_uz,&mesh,1,1);

  // Create Block System with ordering (u,p)
  INT ndof = FE_ux.ndof + FE_uy.ndof;
  if(dim==3) ndof += FE_uz.ndof;
  // Get Global FE Space
  block_fespace FE;
  INT nun = dim;
  INT nspaces = dim;
  FE.nun = nun;
  FE.ndof = ndof;
  FE.nspaces = nspaces;
  FE.var_spaces = (fespace **) calloc(nspaces,sizeof(fespace *));

  FE.var_spaces[0] = &FE_ux;
  FE.var_spaces[1] = &FE_uy;
  if(dim==3) FE.var_spaces[2] = &FE_uz;

  // Set Dirichlet Boundaries
  set_dirichlet_bdry_block(&FE,&mesh);


  clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
  printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
         (REAL) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
  /*******************************************************************************/

  printf("***********************************************************************************\n");
  printf("Number of Elements = %d\tOrder of Quadrature = %d\n",mesh.nelm,2*nq1d-1);
  printf("\n\t--- Element Type ---\n");
  printf("Element Type = %d\n",order_u);
  printf("\n\t--- Degrees of Freedom ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nv,mesh.nedge,mesh.nface);
  printf("\t--> DOF: %d\n",FE.ndof);
  printf("\n\t--- Boundaries ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d\n",mesh.nbv,mesh.nbedge,mesh.nbface);
  printf("***********************************************************************************\n\n");

  /*** Assemble the matrix and right hand side *******************************/
  /* Here we assemble the discrete system:
   *  The weak form is:
   *
   *  <2*eps(u), eps(v)> - <p, div v> = <f, v>
   *                   - <div u, q> = 0
   */
  printf("Assembling the matrix and right-hand side:\n");
  clock_t clk_assembly_start = clock();

  // Allocate the right-hand and declare the csr matrix
  dvector b;

  // Put into Block Form
  block_dCSRmat A;
  bdcsr_alloc(nspaces,nspaces,&A);

  // Assemble the matricies without BC first
  if(dim==2) assemble_global_block(&A,&b,assemble_symmetricDuDv_local,FEM_Block_RHS_Local,&FE,&mesh,cq,source2D,0.0);
  if(dim==3) assemble_global_block(&A,&b,assemble_symmetricDuDv_local,FEM_Block_RHS_Local,&FE,&mesh,cq,source3D,0.0);

  // Eliminate boundary conditions in matrix and rhs
  if(dim==2) {
    eliminate_DirichletBC_blockFE_blockA(bc2D,&FE,&mesh,&b,&A,0.0);
  }
  if(dim==3) {
    eliminate_DirichletBC_blockFE_blockA(bc3D,&FE,&mesh,&b,&A,0.0);
  }

  clock_t clk_assembly_end = clock();
  printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL)
         (clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  /************ Prepare Preconditioners **************************************************/

  // Prepare diagonal blocks
  dCSRmat *A_diag;
  A_diag = (dCSRmat *)calloc(nspaces, sizeof(dCSRmat));

  for(i=0;i<dim;i++){ // copy block diagonal to A_diag
    dcsr_alloc(A.blocks[i*(dim+1)]->row, A.blocks[i*(dim+1)]->col, A.blocks[i*(dim+1)]->nnz, &A_diag[i]);
    dcsr_cp(A.blocks[i*(dim+1)], &A_diag[i]);
  }
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

  // Set initial guess to be all zero
  dvec_set(sol.row, &sol, 0.0);

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
      if(dim==2) solver_flag = linear_solver_bdcsr_krylov_block_3(&A, &b, &sol, &linear_itparam, NULL, A_diag);
      if(dim==3) solver_flag = linear_solver_bdcsr_krylov_block_4(&A, &b, &sol, &linear_itparam, NULL, A_diag);
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
  REAL* solerrL2 = (REAL *) calloc(nspaces, sizeof(REAL));
  REAL* solerrH1 = (REAL *) calloc(nspaces, sizeof(REAL)); // Note: No H1 error for P0 elements
  if(dim==2){
    L2error_block(solerrL2, sol.val, exact_sol2D, &FE, &mesh, cq, 0.0);
    HDerror_block(solerrH1, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
  } else if(dim==3){
    L2error_block(solerrL2, sol.val, exact_sol3D, &FE, &mesh, cq, 0.0);
    HDerror_block(solerrH1, sol.val, exact_sol3D, Dexact_sol3D, &FE, &mesh, cq, 0.0);
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

  // Plotting
  get_unknown_component(&v_ux,&sol,&FE,0);
  get_unknown_component(&v_uy,&sol,&FE,1);
  if(dim==3) get_unknown_component(&v_uz,&sol,&FE,2);

  char** varname;
  if(inparam.print_level > 3){
    char* soldump = "output/solution.vtu";
    varname = malloc(10*FE.nspaces*sizeof(char *));
    varname[0] = "ux";
    varname[1] = "uy";
    if(dim==3) varname[2] = "uz";
    dump_blocksol_vtk(soldump,varname,&mesh,&FE,sol.val);

    // Print in Matlab format to show vector field in nice way.
    if(dim==3) print_matlab_vector_field(&v_ux,&v_uy,&v_uz,&FE_ux);
  }
  /************ Free All the Arrays ***********************************************************/
  // CSR
  bdcsr_free( &A );
  for(i=0;i<nspaces;i++)
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
