/*! \file navierstokes.c
*
*  Created by James Adler on 1/30/24.
*  Copyright 2015_HAZMATH__. All rights reserved.
*
* \brief This program solves steady-state Navier-Stokes PDE using finite elements
*
*      -2*div(eps(u)) + grad(p) + Re*u*grad u = f
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
* \note The variational forms are found in ns_system.h and all Problem Data is found in
*       ns_data.h.  
*
* \note P2-P0 in 3D is not necessarily stable - use with caution!
*
*/

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "hazmath.h"
#include "ns_data.h"
#include "ns_system.h"
#include "stokes_precond.h"
/*********************************************************************************/

/****** MAIN DRIVER **************************************************************/
int main (int argc, char* argv[])
{

  printf("\n===========================================================================\n");
  printf("Beginning Program to solve steady Navier-Stokes Equation.\n");
  printf("===========================================================================\n");

  /****** INITIALIZE PARAMETERS **************************************************/
  INT i;

  // Overall CPU Timing
  clock_t clk_overall_start = clock();

  // Set Parameters from Reading in Input File
  input_param inparam;
  param_input_init(&inparam);
  param_input("./input.dat", &inparam);

  // Time the mesh generation and FE setup
  clock_t clk_mesh_start = clock();

  // Create the mesh
  INT read_mesh_from_file=0; // Use HAZmath built in mesh generator
  printf("\nCreating mesh and FEM spaces:\n");

  // Use HAZMATH built in functions for a uniform mesh in 2D or 3D
  INT dim = inparam.spatial_dim;                 // dimension of computational domain
  INT mesh_ref_levels=inparam.refinement_levels; // refinement levels
  INT mesh_ref_type=inparam.refinement_type;     // refinement type (>10 uniform or <10 other)
  INT set_bndry_codes=inparam.boundary_codes;    // set flags for the boundary DoF (1-16 are Dirichlet)
  scomplex *sc = make_uniform_mesh(dim,mesh_ref_levels,mesh_ref_type,set_bndry_codes);
  relabel_mesh(sc); // relabel the mesh so that boundary flags are user-defined in ns_data.h

  // Get Quadrature Nodes for the Mesh
  INT nq1d = inparam.nquad; /* Quadrature points per dimension */
  qcoordinates *cq = get_quadrature(sc,nq1d);

  // Get info for and create FEM spaces
  // Order of elements: 0 - P0; 1 - P1; 2 - P2; 20 - Nedlec; 30 - Raviart-Thomas
  INT order_u = 2;
  INT order_p = 1;

  // Need Spaces for each component of the velocity plus pressure
  fespace FE_ux; // Velocity in x direction
  create_fespace(&FE_ux,sc,order_u);
  fespace FE_uy; // Velocity in y direction
  create_fespace(&FE_uy,sc,order_u);
  fespace FE_uz; // Velocity in z direction
  if(dim==3) create_fespace(&FE_uz,sc,order_u);
  fespace FE_p; // Pressure
  create_fespace(&FE_p,sc,order_p);
  // Relabel boundary flags as defined in ns_data.h
  relabel_boundary(&FE_ux,dim);
  relabel_boundary(&FE_uy,dim);
  if(dim==3) relabel_boundary(&FE_uz,dim);
  relabel_boundary(&FE_p,dim);

  // Set Dirichlet Boundaries
  set_dirichlet_bdry(&FE_ux,sc,1,16);
  set_dirichlet_bdry(&FE_uy,sc,1,16);
  if(dim==3) set_dirichlet_bdry(&FE_uz,sc,1,16);
  set_dirichlet_bdry(&FE_p,sc,-1,-1);

  // Create Block System with ordering (u,p)
  INT udof = FE_ux.ndof + FE_uy.ndof; // Total DoF for u
  if(dim==3) udof += FE_uz.ndof;
  INT pdof = FE_p.ndof; // Total DoF for p
  INT ndof = udof + pdof; // Total DoF in System
  INT nspaces = dim+1; // Number of FE spaces in system
  INT nun = dim+1; // Number of unkonwns (scalar quantities)
  // Get Global FE Space
  block_fespace FE;
  initialize_fesystem(&FE,nspaces,nun,ndof,sc->fem->ns_leaf);
  FE.var_spaces[0] = &FE_ux;
  FE.var_spaces[1] = &FE_uy;
  if(dim==3) FE.var_spaces[2] = &FE_uz;
  FE.var_spaces[dim] = &FE_p;

  // Set Dirichlet Boundaries and DoF flags
  set_dirichlet_bdry_block(&FE,sc);

  clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
  printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
  (REAL) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
  /*******************************************************************************/

  // Summarize Setup
  printf("***********************************************************************************\n");
  printf("\t--- %lld-dimensional grid ---\n",(long long )dim);
  printf("Number of Elements = %lld\tOrder of Quadrature = %lld\n",(long long )sc->fem->ns_leaf,2*(long long )nq1d-1);
  printf("\n\t--- Element Type ---\n");
  printf("Velocity Element Type = %lld\tPressure Element Type = %lld\n",(long long )order_u,(long long )order_p);
  printf("\n\t--- Degrees of Freedom ---\n");
  printf("Vertices: %-7lld\tEdges: %-7lld\tFaces: %-7lld",(long long )sc->nv,(long long )sc->fem->nedge,(long long )sc->fem->nface);
  printf("\t--> DOF: %lld\n",(long long )FE.ndof);
  printf("***********************************************************************************\n\n");

  /*** Assemble the matrix and right hand side *******************************/
  /* Here we assemble the discrete system:
  *  The weak form is:
  *
  *  <2*eps(u), eps(v)> + Re<u*grad(u),v> - <p, div v> = <f, v>
  *                   - <div u, q> = 0
  * 
  *  This is nonlinear, so we use Newton's method to solve the system.
  *  
  *  <2*eps(u), eps(v)> + Re<u0*grad(u) + u*grad(u0),v> - <p, div v> = <f, v> - <2*eps(u0), eps(v)> - Re<u0*grad(u0),v> + <p0, div v>
  *              - <div u, q> = <div u0, q>
  * where u0 is the previous Newton iterate.
  */

  // Initialize Newton Stepping
  printf("Performing Newton Steps:\n");
  clock_t clk_newton_start = clock();
  newton n_it;
  n_it.isblock = 1; // All matrices are blocks
  initialize_newton(&n_it,&inparam,ndof,nun);
  // Change Newton weighting if desired
  n_it.step_length = 1.0;

  // Set initial guess for Newton
  blockFE_Evaluate(n_it.sol->val,initial_guess,&FE,sc,0.0);

  // Dump Initial Guess
  char solout[40];
  char solfinalout[80];
  char** varname = malloc(50*(FE.nspaces)*sizeof(char *));
  varname[0] = "u1";
  varname[1] = "u2";
  if(dim==3) varname[2] = "u3";
  varname[dim] = "p";
  if (inparam.print_level > 3) {
    sprintf(solout,"output/solution_newt000.vtu");
    dump_blocksol_vtk(solout,varname,sc,&FE,n_it.sol->val);
  }

  // Perform initial Jacobian assembly
  assemble_global_system(n_it.Jac_block,n_it.rhs,&FE,sc,cq,local_assembly_NS,n_it.sol,sourcerhs,NULL,0.0);

  // Eliminate Dirichlet boundary conditions in matrix and rhs
  eliminate_DirichletBC_blockFE_blockA(bc,&FE,sc,n_it.rhs,n_it.Jac_block,0.0);

 // Compute Initial Nonlinear Residual + scale for size
  get_residual_norm(&n_it);
  REAL res_norm_scaled = n_it.res_norm/sqrt(n_it.rhs->row);
  printf("\nInitial Residuals:\nl2-norm of Nonlinear Residual = %25.16e\n",n_it.res_norm);
  printf("               Scaled Version = %25.16e\n",res_norm_scaled);
  fflush(stdout);

  // FILE* Jac = fopen("Jac.dat","w");
  //   bdcsr_print_matlab(Jac, n_it.Jac_block);
  //   dvec_write("rhs.dat",n_it.rhs);

  // Check Convergence before starting
  INT newton_stop=0;
  if(res_norm_scaled < n_it.tol) {
    printf("The initial nonlinear residual is below the tolerance.  Not doing any stepping.\n\n");
    newton_stop=1;
  }
  /*******************************************************************************************/

  /************ Prepare Preconditioners and Linear Solvers ***********************************/
  // Set parameters for linear iterative methods
  linear_itsolver_param linear_itparam;
  param_linear_solver_init(&linear_itparam);
  param_linear_solver_set(&linear_itparam, &inparam);
  INT solver_flag=-20;

  // Set parameters for algebriac multigrid methods
  AMG_param amgparam;
  param_amg_init(&amgparam);
  param_amg_set(&amgparam, &inparam);

  // Prepare diagonal blocks
  dCSRmat *A_diag;
  A_diag = (dCSRmat *)calloc(nun, sizeof(dCSRmat));

  // Get Mass Matrix for p
  dCSRmat Mp = dcsr_create(0,0,0);
  assemble_global_single(&Mp, NULL, &FE_p, sc, cq,
                         local_assembly_mass, NULL, NULL, one_coeff_scal, 0.0);
  dcsr_alloc(Mp.row, Mp.col, Mp.nnz, &A_diag[dim]);
  dcsr_cp(&Mp, &A_diag[dim]);
  /*******************************************************************************************/

  // Perform Newton Steps
  while(!newton_stop) {

    // Update Newton Step Data
    update_newtonstep(&n_it);
    printf("==================================\n");
    printf("\t\tNewton step %d\t\t\n",n_it.current_step);
    printf("==================================\n");

    // Prepare block preconditioner
    // For velcocities, we use the diagonal blocks of the velocity block
    for(i=0;i<dim;i++){
      dcsr_alloc(n_it.Jac_block->blocks[i*(nun+1)]->row, n_it.Jac_block->blocks[i*(nun+1)]->col, n_it.Jac_block->blocks[i*(nun+1)]->nnz, &A_diag[i]);
      dcsr_cp(n_it.Jac_block->blocks[i*(nun+1)], &A_diag[i]);
    }

    // Solve
    clock_t clk_solve_start = clock();
    if(linear_itparam.linear_itsolver_type == 0) { // Direct Solver
      printf(" --> using UMFPACK's Direct Solver:\n");
      solver_flag = block_directsolve_HAZ(n_it.Jac_block,n_it.rhs,n_it.update,linear_itparam.linear_print_level);
    } else { // Iterative Solver
      if (linear_itparam.linear_precond_type == PREC_NULL) {
        solver_flag = linear_solver_bdcsr_krylov(n_it.Jac_block,n_it.rhs,n_it.update,&linear_itparam);

      } else {
        solver_flag = linear_solver_bdcsr_krylov_block(n_it.Jac_block,n_it.rhs,n_it.update, &linear_itparam, &amgparam, A_diag);
      }
    }

    // Free preconditioner for next Newton step
    for(i=0;i<dim;i++) dcsr_free(&A_diag[i]);

    // Error Check
    if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n", solver_flag);
    clock_t clk_solve_end = clock();
    printf("Elapsed CPU Time for Solve = %f seconds.\n\n",(REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);
    fflush(stdout);

    // Update Newton Solution
    update_sol_newton(&n_it);

    printf("updatenorm = %f\n",n_it.update_norm);
    // Get norm of update
    get_blockupdate_norm(&n_it,&FE,sc,cq);
    printf("updatenorm = %f\n",n_it.update_norm);


    // Update Jacobian and nonlinear residual with new solution
    assemble_global_system(n_it.Jac_block,n_it.rhs,&FE,sc,cq,local_assembly_NS,n_it.sol,sourcerhs,NULL,0.0);

    // Eliminate Dirichlet boundary conditions in matrix and rhs
    eliminate_DirichletBC_blockFE_blockA(bc,&FE,sc,n_it.rhs,n_it.Jac_block,0.0);

    // Compute Nonlinear Residual and scaled version
    get_residual_norm(&n_it);
    res_norm_scaled = n_it.res_norm/sqrt(n_it.rhs->row);
    printf("Residuals:\nl2-norm of Nonlinear Residual = %25.16e\n",n_it.res_norm);
    printf("               Scaled Version = %25.16e\n",res_norm_scaled);
    printf("\nL2-norm of Update             = %25.16e\n",n_it.update_norm);
    fflush(stdout);

    // Check for Convergence
    newton_stop = check_newton_convergence(&n_it);

    // Dump solution
    if (inparam.print_level > 3) {
      // Solution at each timestep
      sprintf(solout,"output/solution_newt%03d.vtu",n_it.current_step);
      dump_blocksol_vtk(solout,varname,sc,&FE,n_it.sol->val);
    }
  }

  clock_t clk_newton_end = clock();
  printf(" --> elapsed CPU time for Newton = %f seconds.\n\n",(REAL)
  (clk_newton_end-clk_newton_start)/CLOCKS_PER_SEC);
  fflush(stdout);


  /********************* Compute Errors if you have exact solution ****************************/
  clock_t clk_error_start = clock();
  REAL* solerrL2 = (REAL *) calloc(nspaces, sizeof(REAL));
  REAL* solerrH1 = (REAL *) calloc(nspaces, sizeof(REAL)); // Note: No H1 error for P0 elements
  if(dim==2){
    L2error_block(solerrL2, n_it.sol->val, exact_sol2D, &FE, sc, cq, 0.0);
    HDerror_block(solerrH1, n_it.sol->val, exact_sol2D, Dexact_sol2D, &FE, sc, cq, 0.0);
  } else if(dim==3){
    L2error_block(solerrL2, n_it.sol->val, exact_sol3D, &FE, sc, cq, 0.0);
    HDerror_block(solerrH1, n_it.sol->val, exact_sol3D, Dexact_sol3D, &FE, sc, cq, 0.0);
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

  // Collect solutions in pvd file
  if(inparam.print_level > 3) {
    create_pvd("output/solution.pvd",n_it.current_step+1,"solution_newt","timestep");
  }
  /************ Free All the Arrays ***********************************************************/

  // CSR
  dcsr_free(&Mp);
  for(i=0;i<dim+1;i++)
    dcsr_free(&A_diag[i]);
  if(A_diag) free(A_diag);

  // Vectors
  if(solerrL2) free(solerrL2);
  if(solerrH1) free(solerrH1);

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
  haz_scomplex_free(sc);

  // Strings
  if(varname) free(varname);

  // Newton
  free_newton(&n_it);

  /*******************************************************************************************/
  clock_t clk_overall_end = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",
         (REAL) (clk_overall_end-clk_overall_start)/CLOCKS_PER_SEC);
  return 0;
}  /* End of Program */
/*********************************************************************************************/