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
*       stokes_data.h.  This example also illustrates how to construct block
*       versions of the finite-element spaces, linear systems, and
*       how to use block solvers and special preconditioners designed for
*       Stokes.  These are found in stokes_precond.h.
*
*/

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "hazmath.h"
#include "stokes_data.h"
#include "stokes_system.h"
#include "stokes_precond.h"
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

  // Create the mesh
  INT read_mesh_from_file=0;
  printf("\nCreating mesh and FEM spaces:\n");

  // Time the mesh generation and FE setup
  clock_t clk_mesh_start = clock();

  // Use HAZMATH built in functions for a uniform mesh in 2D or 3D
  mesh_struct mesh;
  INT dim = inparam.spatial_dim;                 // dimension of computational domain
  INT mesh_ref_levels=inparam.refinement_levels; // refinement levels
  INT mesh_ref_type=inparam.refinement_type;     // refinement type (>10 uniform or <10 other)
  INT set_bndry_codes=inparam.boundary_codes;    // set flags for the boundary DoF (1-16 are Dirichlet)
  mesh=make_uniform_mesh(dim,mesh_ref_levels,mesh_ref_type,set_bndry_codes);

  // Get Quadrature Nodes for the Mesh
  INT nq1d = inparam.nquad; /* Quadrature points per dimension */
  qcoordinates *cq = get_quadrature(&mesh,nq1d);

  // Get info for and create FEM spaces
  // Order of elements: 0 - P0; 1 - P1; 2 - P2; 20 - Nedlec; 30 - Raviart-Thomas
  INT order_u = 2;
  INT order_p = 1;

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
  // u = g on all boundaries
  // p is Neumann on all boundaries
  // The mesh is set up so that flag values 1-16 are Dirichlet and 17-32 are Neumann
  // For now mesh marks all physical boundaries as set_bndry_codes value
  set_dirichlet_bdry(&FE_ux,&mesh,1,16);
  set_dirichlet_bdry(&FE_uy,&mesh,1,16);
  if(dim==3) set_dirichlet_bdry(&FE_uz,&mesh,1,16);
  set_dirichlet_bdry(&FE_p,&mesh,-1,-1); // No DoF for p should be Dirichlet

  // Create Block System with ordering (u,p)
  INT udof = FE_ux.ndof + FE_uy.ndof; // Total DoF for u
  if(dim==3) udof += FE_uz.ndof;
  INT pdof = FE_p.ndof; // Total DoF for p
  INT ndof = udof + pdof; // Total DoF in System
  INT nspaces = dim+1; // Number of FE spaces in system
  INT nun = dim+1; // Number of unkonwns (scalar quantities)
  // Get Global FE Space
  block_fespace FE;
  initialize_fesystem(&FE,nspaces,nun,ndof,mesh.nelm);
  FE.var_spaces[0] = &FE_ux;
  FE.var_spaces[1] = &FE_uy;
  if(dim==3) FE.var_spaces[2] = &FE_uz;
  FE.var_spaces[dim] = &FE_p;

  // Set Dirichlet Boundaries and DoF flags
  set_dirichlet_bdry_block(&FE,&mesh);

  clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
  printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
  (REAL) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
  /*******************************************************************************/

  // Summarize Setup
  printf("***********************************************************************************\n");
  printf("\t--- %d-dimensional grid ---\n",dim);
  printf("Number of Elements = %d\tOrder of Quadrature = %d\n",mesh.nelm,2*nq1d-1);
  printf("\n\t--- Element Type ---\n");
  printf("Velocity Element Type = %d\tPressure Element Type = %d\n",order_u,order_p);
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

  // Allocate the right-hand and declare the block csr matrix
  dvector b;

  // Put into Block Form
  block_dCSRmat A;
  bdcsr_alloc(nspaces,nspaces,&A);

  // Assemble the matricies without BC first
  if(dim==2) assemble_global_block(&A,&b,local_assembly_Stokes,FEM_Block_RHS_Local,&FE,&mesh,cq,source2D,0.0);
  if(dim==3) assemble_global_block(&A,&b,local_assembly_Stokes,FEM_Block_RHS_Local,&FE,&mesh,cq,source3D,0.0);

  // Eliminate boundary conditions in matrix and rhs
  if(dim==2) eliminate_DirichletBC_blockFE_blockA(bc2D,&FE,&mesh,&b,&A,0.0);
  if(dim==3) eliminate_DirichletBC_blockFE_blockA(bc3D,&FE,&mesh,&b,&A,0.0);

  clock_t clk_assembly_end = clock();
  printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL)
  (clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  /************ Solve **************************************************/
  printf("Solving the System:\n");
  clock_t clk_solve_start = clock();
  INT solver_flag = -20;

  // Prepare diagonal blocks for preconditioners
  dCSRmat *A_diag;
  A_diag = (dCSRmat *)calloc(dim+1, sizeof(dCSRmat));

  // For velcocities, we use the diagonal blocks of the velocity block
  for(i=0;i<dim;i++){
    dcsr_alloc(A.blocks[i*(dim+2)]->row, A.blocks[i*(dim+2)]->col, A.blocks[i*(dim+2)]->nnz, &A_diag[i]);
    dcsr_cp(A.blocks[i*(dim+2)], &A_diag[i]);
  }
  // For pressure, we use the mass matrix
  dCSRmat Mp;
  assemble_global(&Mp,NULL,assemble_mass_local,&FE_p,&mesh,cq,NULL,one_coeff_scal,0.0);
  dcsr_alloc(Mp.row, Mp.col, Mp.nnz, &A_diag[dim]);
  dcsr_cp(&Mp, &A_diag[dim]);

  // Allocate solution and set initial guess to be all zero
  dvector sol = dvec_create(ndof);
  dvec_set(sol.row, &sol, 0.0);

  //  Apply Pressure "BCs" (removes singularity)
  INT pressureloc = 0;
  REAL pressurecoord[dim];
  for(i=0;i<dim;i++) pressurecoord[i] = 0.0;
  REAL solval[dim+1];
  if(dim==2) exact_sol2D(solval,pressurecoord,0.0,NULL);
  if(dim==3) exact_sol3D(solval,pressurecoord,0.0,NULL);
  // Set initial guess for pressure to match the known "boundary condition" for pressure
  sol.val[udof + pressureloc]  = solval[dim];

  // Set parameters for linear iterative methods
  linear_itsolver_param linear_itparam;
  param_linear_solver_init(&linear_itparam);
  param_linear_solver_set(&linear_itparam,&inparam);
  INT solver_type = linear_itparam.linear_itsolver_type;
  INT solver_printlevel = linear_itparam.linear_print_level;

  /* Set parameters for algebriac multigrid methods */
  AMG_param amgparam;
  param_amg_init(&amgparam);
  param_amg_set(&amgparam, &inparam);
  param_amg_print(&amgparam);

  // Solve
  if(solver_type==0) { // Direct Solver

    solver_flag = block_directsolve_HAZ(&A,&b,&sol,solver_printlevel);

  } else { // Iterative Solver

    if (linear_itparam.linear_precond_type == PREC_NULL) { // No Preconditioner
      solver_flag = linear_solver_bdcsr_krylov(&A, &b, &sol, &linear_itparam);

    } else if (linear_itparam.linear_precond_type >= 10 && linear_itparam.linear_precond_type < 40) {  // General Block Preconditioner
      if(dim==2) solver_flag = linear_solver_bdcsr_krylov_block_3(&A, &b, &sol, &linear_itparam, &amgparam, A_diag);
      if(dim==3) solver_flag = linear_solver_bdcsr_krylov_block_4(&A, &b, &sol, &linear_itparam, &amgparam, A_diag);

    } else { // Special Preconditioner for Stokes

      // get Stokes preconditioner data
      precond_block_data *precdata = get_precond_block_data_stokes(&A, &Mp, &linear_itparam, &amgparam);

      // Setup the preconditioner and choose type
      precond prec;
      prec.data = precdata;
      set_precond_type_stokes(&prec,&linear_itparam);

      // Solve
      solver_flag = solver_bdcsr_linear_itsolver(&A, &b, &sol, &prec, &linear_itparam);

      // Free the preconditioner
      precond_block_data_free_stokes(precdata);
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
  REAL perrL2 = solerrL2[dim];
  REAL perrH1 = solerrH1[dim];
  if(order_p==0) perrH1 = 0.0;

  printf("*******************************************************\n");
  printf("L2 Norm of u error    = %26.13e\n",uerrL2);
  printf("L2 Norm of p error    = %26.13e\n",perrL2);
  printf("H1 Norm of u error    = %26.13e\n",uerrH1);
  if(order_p==0) {
    printf("H1 Norm of p error    = %26.13e (not computed for P0 elements)\n",perrH1);
  } else {
    printf("H1 Norm of p error    = %26.13e\n",perrH1);
  }
  printf("*******************************************************\n\n");
  clock_t clk_error_end = clock();
  printf("Elapsed CPU time for getting errors = %lf seconds.\n\n",(REAL)
  (clk_error_end-clk_error_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  // Plotting
  dvector v_ux = dvec_create(FE_ux.ndof);
  dvector v_uy = dvec_create(FE_uy.ndof);
  dvector v_uz;
  if(dim==3) v_uz = dvec_create(FE_uz.ndof);
  dvector v_p = dvec_create(FE_p.ndof);
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

    // Print in Matlab format to show vector field in nice way.
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
