/*! \file examples/navier_stokes/navier_stokes.c
 *
 *  Created by James Adler on 05/20/2020
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program solves Naiver-Stokes PDE using finite elements
 *
 *      -2*div(eps(u)) + Re*u*grad u + grad(p) = f
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
 *       Note the scaling of Reynolds Number here.  We just put this on the
 *       nonlinear inertia term.  We assume the pressure, p, contains the appropriate
 *       scaling.  This allows for the case of Re=0 equaling Stokes' Equations.
 *
 * \note This example shows how to build your own bilinear form for a nonlinear
 *       system and time-dependent problem. The variational forms are found in
 *       ns_system.h and all Problem Data is found in ns_data.h.
 *
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "hazmath.h"
#include "ns_data.h"
#include "ns_system.h"
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

  /*** Assemble the Initial Jacobian matrix and right hand side **************/
  /* Here we assemble the discrete system:
   *  The nonlinear weak form is:
   *
   *  2*<eps(u), eps(v)> + Re*<u*grad u,v> - <p, div v> = <f, v>
   *                   - <div u, q> = 0
   */

   // Initialize Newton Stepping
   printf("Performing Newton Steps:\n");
   clock_t clk_newton_start = clock();
   newton n_it;
   n_it.isblock = 1; // All matrices are blocks
   initialize_newton(&n_it,&inparam,ndof,FE.nspaces);

   // Change Newton weighting if desired
   n_it.step_length = 1.0;

   // Allocate the Initial Guess
   if(dim==2) blockFE_Evaluate(n_it.sol->val,initial_guess2D,&FE,&mesh,0.0);
   if(dim==3) blockFE_Evaluate(n_it.sol->val,initial_guess3D,&FE,&mesh,0.0);

   // Dump Initial Guess
   char solout[40];
   char** varname = malloc(50*FE.nspaces*sizeof(char *));
   varname[0] = "u1";
   varname[1] = "u2";
   if(dim==3) varname[2] = "u3";
   varname[dim] = "p";
   if (inparam.print_level > 3) {
     sprintf(solout,"output/solution_newt000.vtu");
     dump_blocksol_vtk(solout,varname,&mesh,&FE,n_it.sol->val);
   }

   // Perform initial Jacobian assembly
   if(dim==2) assemble_global_Jacobian(n_it.Jac_block,n_it.rhs,n_it.sol,local_assembly_NavierStokes,&FE,&mesh,cq,source2D,0.0);
   if(dim==3) assemble_global_Jacobian(n_it.Jac_block,n_it.rhs,n_it.sol,local_assembly_NavierStokes,&FE,&mesh,cq,source3D,0.0);

   // Eliminate Dirichlet boundary conditions in matrix and rhs
   if(dim==2) eliminate_DirichletBC_blockFE_blockA(bc2D,&FE,&mesh,n_it.rhs,n_it.Jac_block,0.0);
   if(dim==3) eliminate_DirichletBC_blockFE_blockA(bc3D,&FE,&mesh,n_it.rhs,n_it.Jac_block,0.0);

  /**************************************************/
  //  Apply Pressure "BCs" (removes singularity)
  REAL pressureval =0.5;
  INT pressureloc = 0;
  /**************************************************/

  /*******************************************************************************************/

  /************ Prepare Preconditioners **************************************************/

  // // Prepare diagonal blocks
  // dCSRmat *A_diag;
  // A_diag = (dCSRmat *)calloc(dim+1, sizeof(dCSRmat));
  //
  // for(i=0;i<dim;i++){ // copy block diagonal to A_diag
  //   dcsr_alloc(A.blocks[i*(dim+2)]->row, A.blocks[i*(dim+2)]->col, A.blocks[i*(dim+2)]->nnz, &A_diag[i]);
  //   dcsr_cp(A.blocks[i*(dim+2)], &A_diag[i]);
  // }
  //
  // // Get Mass Matrix for p
  // dCSRmat Mp;
  // assemble_global(&Mp,NULL,assemble_mass_local,&FE_p,&mesh,cq,NULL,one_coeff_scal,0.0);
  // dcsr_alloc(Mp.row, Mp.col, Mp.nnz, &A_diag[dim]);
  // dcsr_cp(&Mp, &A_diag[dim]);
  /*******************************************************************************************/

  // Set parameters for linear iterative methods
  linear_itsolver_param linear_itparam;
  param_linear_solver_init(&linear_itparam);
  param_linear_solver_set(&linear_itparam, &inparam);
  INT solver_flag=-20;

  // Compute Initial Nonlinear Residual + scale for size
  get_residual_norm(&n_it);
  REAL res_norm_scaled = n_it.res_norm/sqrt(n_it.rhs->row);
  printf("\nInitial Residuals:\nl2-norm of Nonlinear Residual = %25.16e\n",n_it.res_norm);
  printf("               Scaled Version = %25.16e\n",res_norm_scaled);

  // Check Convergence before starting
  INT newton_stop=0;
  if(res_norm_scaled < n_it.tol) {
    printf("The initial nonlinear residual is below the tolerance.  Not doing any stepping.\n\n");
    newton_stop=1;
  }

  // Perform Newton Steps
  while(!newton_stop) {

    // Update Newton Step Data
    update_newtonstep(&n_it);
    printf("==================================\n");
    printf("\t\tNewton step %d\t\t\n",n_it.current_step);
    printf("==================================\n");

    // Solve
    clock_t clk_solve_start = clock();
    if(linear_itparam.linear_itsolver_type == 0) { // Direct Solver
      printf(" --> using UMFPACK's Direct Solver:\n");
      solver_flag = block_directsolve_UMF(n_it.Jac_block,n_it.rhs,n_it.update,linear_itparam.linear_print_level);
    } else { // Iterative Solver
      solver_flag = linear_solver_bdcsr_krylov(n_it.Jac_block,n_it.rhs,n_it.update,&linear_itparam);
    }

    // Error Check
    if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n", solver_flag);
    clock_t clk_solve_end = clock();
    printf("Elapsed CPU Time for Solve = %f seconds.\n\n",(REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);

    // Update Newton Solution
    update_sol_newton(&n_it);

    // Get norm of update
    get_blockupdate_norm(&n_it,&FE,&mesh,cq);

    // Update Jacobian and nonlinear residual with new solution
    if(dim==2) assemble_global_Jacobian(n_it.Jac_block,n_it.rhs,n_it.sol,local_assembly_NavierStokes,&FE,&mesh,cq,source2D,0.0);
    if(dim==3) assemble_global_Jacobian(n_it.Jac_block,n_it.rhs,n_it.sol,local_assembly_NavierStokes,&FE,&mesh,cq,source3D,0.0);

    // Eliminate Dirichlet boundary conditions in matrix and rhs
    if(dim==2) eliminate_DirichletBC_blockFE_blockA(bc2D,&FE,&mesh,n_it.rhs,n_it.Jac_block,0.0);
    if(dim==3) eliminate_DirichletBC_blockFE_blockA(bc3D,&FE,&mesh,n_it.rhs,n_it.Jac_block,0.0);

    // Compute Nonlinear Residual and scaled version
    get_residual_norm(&n_it);
    res_norm_scaled = n_it.res_norm/sqrt(n_it.rhs->row);

    // Print Data
    printf("Residuals:\nl2-norm of Nonlinear Residual = %25.16e\n",n_it.res_norm);
    printf("               Scaled Version = %25.16e\n",res_norm_scaled);
    printf("\nL2-norm of Update             = %25.16e\n",n_it.update_norm);

    // Check for Convergence
    newton_stop = check_newton_convergence(&n_it);

    if (inparam.print_level > 3) {
      // Solution at each timestep
      sprintf(solout,"output/solution_newt%03d.vtu",n_it.current_step);
      dump_blocksol_vtk(solout,varname,&mesh,&FE,n_it.sol->val);
    }
  }

  clock_t clk_newton_end = clock();
  printf(" --> elapsed CPU time for Newton = %f seconds.\n\n",(REAL)
  (clk_newton_end-clk_newton_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  // Dump all data in one vtu file
  if(inparam.print_level > 3) {
    create_pvd("output/solution.pvd",n_it.current_step+1,"solution_newt","timestep");
  }
  /*******************************************************************************************/

  /********************* Compute Errors if you have exact solution ****************************/
  clock_t clk_error_start = clock();
  REAL* solerrL2 = (REAL *) calloc(dim+1, sizeof(REAL));
  REAL* solerrH1 = (REAL *) calloc(dim+1, sizeof(REAL)); // Note: No H1 error for P0 elements
  if(dim==2){
    L2error_block(solerrL2, n_it.sol->val, exact_sol2D, &FE, &mesh, cq, 0.0);
    if(order_p > 0) HDerror_block(solerrH1, n_it.sol->val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
  } else if(dim==3){
    L2error_block(solerrL2, n_it.sol->val, exact_sol3D, &FE, &mesh, cq, 0.0);
    if(order_p > 0) HDerror_block(solerrH1, n_it.sol->val, exact_sol3D, Dexact_sol3D, &FE, &mesh, cq, 0.0);
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

  // // Plotting
  // get_unknown_component(&v_ux,&sol,&FE,0);
  // get_unknown_component(&v_uy,&sol,&FE,1);
  // if(dim==3) get_unknown_component(&v_uz,&sol,&FE,2);
  // get_unknown_component(&v_p,&sol,&FE,dim);
  //
  // if(inparam.print_level > 3){
  //   // Print in Matlab format to show vector field in nice way.
  //   if(dim==3) print_matlab_vector_field(&v_ux,&v_uy,&v_uz,&FE_ux);
  // }

  /************ Free All the Arrays ***********************************************************/
  // CSR
  // bdcsr_free( &A );
  // dcsr_free( &Mp);
  // for(i=0;i<dim+1;i++)
  //   dcsr_free( &A_diag[i] );
  // if(A_diag) free(A_diag);

  // Vectors
  if(solerrL2) free(solerrL2);
  if(solerrH1) free(solerrH1);
  // dvec_free( &b );
  // dvec_free( &sol );
  // dvec_free( &v_ux );
  // dvec_free( &v_uy );
  // if(dim==3) dvec_free( &v_uz );
  // dvec_free( &v_p );

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

  // Newton
  free_newton(&n_it);

  // Strings
  if(varname) free(varname);

  /*******************************************************************************************/
  clock_t clk_overall_end = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",
         (REAL) (clk_overall_end-clk_overall_start)/CLOCKS_PER_SEC);
  return 0;
}  /* End of Program */
/*********************************************************************************************/
