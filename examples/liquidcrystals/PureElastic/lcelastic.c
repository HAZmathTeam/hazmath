/*! \file lcelastic.c
*
*  Created by James Adler on 9/8/17. Updated by Xiaozhe Hu on 8/10/2023
*  Copyright 2015_HAZMATH__. All rights reserved.
*
* \brief This program solves the pure elastic static LC equations
*
*      min 1/2 K1||div n||^2 + 1/2 K3<Z curl n, curl n>
*               s.t.
*                    n*n = 1
*        including in
*        2D slab geometry -> n(x,y) = (n1(x,y),n2(x,y),n3(x,y))
*        3D geometry      -> n(x,y,z) = (n1(x,y, z), n2(x,y, z), n3(x,y, z))
*
* The minimization leads to a nonlinear weak form for two unknowns, n(x,y), lam(x,y)
*
*       L1 = K1*<div n, div v> + K3<Z(n) curl n, curl v> + (K2-K3)*<n*curl n, v*curl n> + 2*int [(n*v)*lam] = 0
*       L2 = int [(n*n-1)*gam] = 0
*
* One step of Newton then leads to
*       K1<div(dn),div(v)> + K3<Z(nk)*curl(dn),curl(v)> +
*       (K2-K3)*[<nk*curl(nk),dn*curl(v)+v*curl(dn)> +
*                <dn*curl(nk),nk*curl(v) + v*curl(nk)> +
*                <nk*curl(dn),v*curl(nk)>] + int(lamk*dn*v)     +  int(dlam*nk*v)    = -L1(nk,lamk,v)
*
*       int(gam*nk*dn)                                          +  0                 = -L2(nk,lamk,gam)
*
* \note We use nested iteration and adaptive refinement with various marking strategies.
*/

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "hazmath.h"
// Change for different test problems and geometries
#include "../data/2D/lc_basic_twist_data.h"
//#include "../data/2D/lc_harmonic_solution.h"
// #include "../data/2D/lc_nano_pattern.h"
// #include "../data/2D/lc_tilt_twist_data.h"
//#include "../data/3D/lc_harmonic_solution3D.h"
#include "../common/lc_system_combined.h" // assembly
#include "../common/lc_estimator.h"       // error estimator
/*********************************************************************************/

/****** MAIN DRIVER **************************************************************/
int main (int argc, char* argv[])
{
  srand(803087000);
  printf("\n===========================================================================\n");
  printf("Beginning Program to solve the Frank-Oseen Elastic Energy Equation.\n");
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
  printf("\nCreating mesh and FEM spaces:\n");
  clock_t clk_mesh_start = clock();

  // Use HAZMATH built in functions for an initial uniform mesh in 2D or 3D
  // Make sure to save the simplicial complex for adaptive refinement and marking
  scomplex *sc;
  INT dim = inparam.spatial_dim;                 // dimension of computational domain
  INT mesh_ref_levels = inparam.refinement_levels; // refinement levels of initial (coarsest uniform) mesh
  INT mesh_ref_type = inparam.refinement_type;     // refinement type (>10 uniform or <10 other)
  INT set_bndry_codes = inparam.boundary_codes;    // set flags for the boundary DoF (1-16 are Dirichlet) (this isn't used??)
  sc = make_uniform_mesh(dim,mesh_ref_levels,mesh_ref_type,set_bndry_codes);

  // Get Quadrature Nodes for the Mesh
  INT nq1d = inparam.nquad; /* Quadrature points per dimension */
  qcoordinates *cq = get_quadrature(sc,nq1d);

  // Get info for and create FEM spaces
  // Order of elements: 0 - P0; 1 - P1; 2 - P2; 20 - Nedlec; 30 - Raviart-Thomas
  INT order_n = 2;
  INT order_lam = 1;

  // Need Spaces for each component of the director plus lambda
  fespace FE_nx; // Director in x direction
  create_fespace(&FE_nx,sc,order_n);
  fespace FE_ny; // Director in y direction
  create_fespace(&FE_ny,sc,order_n);
  fespace FE_nz; // Director in z direction
  create_fespace(&FE_nz,sc,order_n);
  fespace FE_lam; // Lagrange multiplier, lambda
  create_fespace(&FE_lam,sc,order_lam);

  // Decide whether you want Dirichlet or Periodic boundary conditions.
  // periodic_flag = 1 means it's periodic bc
  // dirichlet_flag = 1 means it's dirichlet bc
  // set in data (*.h) file
  INT bc_flag = characterize_boundary_conditions();
  INT periodic_flag = bc_flag;
  INT dirichlet_flag = 0;
  if (bc_flag == 0){
    dirichlet_flag = 1;
  }
  // Lambda has Neumann boundaries regardless
  set_dirichlet_bdry(&FE_lam,sc,-10,-10);

  // Dirichlet Boundaries (Dirichlet all around for n)
  if (dirichlet_flag == 1){
    set_dirichlet_bdry(&FE_nx,sc,1,6);
    set_dirichlet_bdry(&FE_ny,sc,1,6);
    set_dirichlet_bdry(&FE_nz,sc,1,6);
  }

  // Periodic Boundaries
  if (periodic_flag == 1){
    // Relabel boundary flags for FE space to distinguish
    // relabel_boundary is found in LCSystem.h
    // Boundary codes
    // 1: x=0  2: x=1   3: y=0    4: y=1, 5: z=0, 6: z=1
    relabel_boundary(&FE_nx, dim);
    relabel_boundary(&FE_ny, dim);
    relabel_boundary(&FE_nz, dim);
    relabel_boundary(&FE_lam, dim);

    if(dim==2) {
      // For 2D, the y bounds are Dirichlet while x is periodic.
      set_periodic_bdry(&FE_nx,sc,0.0,1.0,0.0,0.0,0.0,0.0);
      set_periodic_bdry(&FE_ny,sc,0.0,1.0,0.0,0.0,0.0,0.0);
      set_periodic_bdry(&FE_nz,sc,0.0,1.0,0.0,0.0,0.0,0.0);
      set_periodic_bdry(&FE_lam,sc,0.0,1.0,0.0,0.0,0.0,0.0);

      set_dirichlet_bdry(&FE_nx,sc,3,4);
      set_dirichlet_bdry(&FE_ny,sc,3,4);
      set_dirichlet_bdry(&FE_nz,sc,3,4);
    } else if(dim==3) {
      // For 3D, the z bounds are Dirichlet while x & y are periodic
      set_periodic_bdry(&FE_nx,sc,0.0,1.0,0.0,1.0,0.0,0.0);
      set_periodic_bdry(&FE_ny,sc,0.0,1.0,0.0,1.0,0.0,0.0);
      set_periodic_bdry(&FE_nz,sc,0.0,1.0,0.0,1.0,0.0,0.0);
      set_periodic_bdry(&FE_lam,sc,0.0,1.0,0.0,1.0,0.0,0.0);

      // Setting Dirichlet boundary conditions
      set_dirichlet_bdry(&FE_nx,sc,5,6);
      set_dirichlet_bdry(&FE_ny,sc,5,6);
      set_dirichlet_bdry(&FE_nz,sc,5,6);
    } else {
      check_error(ERROR_DIM, __FUNCTION__);
    }
  }

  // ------------------------------------------------------------------

  // Create Block System with ordering (n,lam)
  INT ndof = FE_nx.ndof + FE_ny.ndof + FE_nz.ndof +FE_lam.ndof;
  INT nspaces = 4;
  INT nun = 4;
  // Get Global FE Space
  block_fespace FE;
  initialize_fesystem(&FE,nspaces,nun,ndof,sc->fem->ns_leaf);
  FE.var_spaces[0] = &FE_nx;
  FE.var_spaces[1] = &FE_ny;
  FE.var_spaces[2] = &FE_nz;
  FE.var_spaces[3] = &FE_lam;

  // Set Dirichlet Boundaries
  set_dirichlet_bdry_block(&FE,sc);

  // Set Periodic Boundaries
  block_dCSRmat P_periodic;
  if (periodic_flag == 1){
    generate_periodic_P_blockFE(&FE, &P_periodic);
  }

  clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
  printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
         (REAL) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
  /*******************************************************************************/

  printf("***********************************************************************************\n");
  printf("%lld Dimensional Problem\n",(long long)sc->dim);
  printf("Number of Elements = %lld\tOrder of Quadrature = %d\n",(long long)sc->fem->ns_leaf,2*nq1d-1);
  printf("\n\t--- Element Type ---\n");
  printf("Director Element Type = %d\tLagrange Multiplier Element Type = %d\n",order_n,order_lam);
  printf("\n\t--- Degrees of Freedom ---\n");
  printf("Vertices: %-7lld\tEdges: %-7lld\tFaces: %-7lld",(long long)sc->nv,(long long)sc->fem->nedge,(long long)sc->fem->nface);
  printf("\t--> DOF: %d\n",FE.ndof);
  printf("***********************************************************************************\n\n");

  // Set parameters for linear iterative methods
  linear_itsolver_param linear_itparam;
  param_linear_solver_init(&linear_itparam);
  param_linear_solver_set(&linear_itparam, &inparam);
  INT solver_flag=-20;

  // Initialize Newton Stepping
  printf("Performing Newton Steps:\n");
  clock_t clk_newton_start = clock();
  newton n_it;
  n_it.isblock = 1; // All matrices are blocks
  initialize_newton(&n_it,&inparam,ndof,FE.nspaces);

  // Change Newton weighting if desired
  n_it.step_length = 1.0;

  // Allocate the Initial Guess
  srand(803087000);
  blockFE_Evaluate(n_it.sol->val,initial_guess,&FE,sc,0.0);

  // Error Estimator Stuff
  REAL* errest = (REAL *) calloc(sc->fem->ns_leaf,sizeof(REAL));
  // Plotting stuff
  char estout[40];
  char** estname = malloc(50*sizeof(char *));
  estname[0] = "err_est";
  fespace FE_est; // P0 space to store estimator for plotting
  create_fespace(&FE_est,sc,0);
  // Get estimate of initial guess
  LCerror_estimator(errest,n_it.sol->val,&FE,sc,cq);


  // Dump Initial Guess
  char solout[40];
  char** varname = malloc(50*FE.nspaces*sizeof(char *));
  varname[0] = "n1";
  varname[1] = "n2";
  varname[2] = "n3";
  varname[3] = "lambda";
  if (inparam.print_level > 3) {
    sprintf(solout,"output/solution_newt000.vtu");
    dump_blocksol_vtk(solout,varname,sc,&FE,n_it.sol->val);
    // Dump estimator into vtk
    sprintf(estout,"output/errest_newt000.vtu");
    dump_sol_vtk(estout,estname[0],sc,&FE_est,errest);
  }

  // Perform initial Jacobian assembly
  assemble_global_system(n_it.Jac_block,n_it.rhs,&FE,sc,cq,local_assembly_LCelastic,n_it.sol,NULL,NULL,0.0);
  // Eliminate Dirichlet boundary conditions in matrix and rhs
  eliminate_DirichletBC_blockFE_blockA(bc,&FE,sc,n_it.rhs,n_it.Jac_block,0.0);

  // Deal with Periodic Boundary conditions, but not mess up n_it.Jac_block
  block_dCSRmat Jac_per;
  //dvector* rhs_per = malloc(sizeof(struct dvector));
  dvector rhs_per;
  dvector update_per;
  INT nx_ndof_periodic;
  INT ny_ndof_periodic;
  INT nz_ndof_periodic;
  INT lam_ndof_periodic;
  INT ndof_periodic;

  if (periodic_flag == 1) {
    bdcsr_alloc(4,4,&Jac_per);
    dvec_null(&rhs_per);
    eliminate_PeriodicBC_blockFE_nonoverwrite(&P_periodic, n_it.Jac_block, n_it.rhs, &Jac_per, &rhs_per);

    // Get the new number of degrees of freedom after applying periodic boundary conditions
    nx_ndof_periodic = Jac_per.blocks[0]->row;
    ny_ndof_periodic = Jac_per.blocks[nun+1]->row;
    nz_ndof_periodic = Jac_per.blocks[2*(nun+1)]->row;
    lam_ndof_periodic = Jac_per.blocks[3*(nun+1)]->row;
    ndof_periodic = nx_ndof_periodic+ny_ndof_periodic+nz_ndof_periodic+lam_ndof_periodic;

    // Allocate solution dvector
    update_per = dvec_create(ndof_periodic);
  }

  // -------------------------------------------------------------------------

  // Compute Initial Nonlinear Residual + scale for size
  get_residual_norm(&n_it);
  REAL res_norm_scaled = n_it.res_norm/sqrt(n_it.rhs->row);
  printf("\nInitial Residuals:\nl2-norm of Nonlinear Residual = %25.16e\n",n_it.res_norm);
  printf("               Scaled Version = %25.16e\n",res_norm_scaled);

  // Compute Initial Energy and Length of Director
  REAL* energy = (REAL *) calloc(4,sizeof(REAL));
  compute_LCelastic_energy(energy,n_it.sol->val,&FE,sc,cq);
  REAL unitlength=0.0;
  compute_LCelastic_unitlength(&unitlength,n_it.sol->val,&FE,sc,cq);
  printf("\nInitial Energies & Length of Director:\nTotal Energy: %25.16e\n",energy[0]);
  printf("Splay:        %25.16e\nTwist:        %25.16e\nBend:         %25.16e\n\n",energy[1],energy[2],energy[3]);
  printf("L2 norm of n: %25.16e\n\n",unitlength);

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

    if(periodic_flag==1) {
      if(linear_itparam.linear_itsolver_type == 0) { // Direct Solver
        printf(" --> using UMFPACK's Direct Solver:\n");
        solver_flag = block_directsolve_HAZ(&Jac_per,&rhs_per,&update_per,linear_itparam.linear_print_level);
      } else { // Iterative Solver
        solver_flag = linear_solver_bdcsr_krylov(&Jac_per,&rhs_per,&update_per,&linear_itparam);
      }
      bdcsr_mxv(&P_periodic,update_per.val,n_it.update->val);
    } else {
      if(linear_itparam.linear_itsolver_type == 0) { // Direct Solver
        printf(" --> using UMFPACK's Direct Solver:\n");
        solver_flag = block_directsolve_HAZ(n_it.Jac_block,n_it.rhs,n_it.update,linear_itparam.linear_print_level);
      } else { // Iterative Solver
        solver_flag = linear_solver_bdcsr_krylov(n_it.Jac_block,n_it.rhs,n_it.update,&linear_itparam);
      }
    }

    // Error Check
    if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n", solver_flag);

    clock_t clk_solve_end = clock();
    printf("Elapsed CPU Time for Solve = %f seconds.\n\n",(REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);

    // Update Newton Solution
    update_sol_newton(&n_it);

    // Get norm of update
    get_blockupdate_norm(&n_it,&FE,sc,cq);

    // Update Jacobian and nonlinear residual with new solution
    assemble_global_system(n_it.Jac_block,n_it.rhs,&FE,sc,cq,local_assembly_LCelastic,n_it.sol,NULL,NULL,0.0);

    // Eliminate Dirichlet boundary conditions in matrix and rhs
    eliminate_DirichletBC_blockFE_blockA(bc,&FE,sc,n_it.rhs,n_it.Jac_block,0.0);
    // Eliminate Periodic boundary CONDITIONS
    if (periodic_flag == 1) {
      eliminate_PeriodicBC_blockFE_nonoverwrite(&P_periodic, n_it.Jac_block, n_it.rhs, &Jac_per, &rhs_per);
      // Blow nonlinear residual back up to calculate norm
      bdcsr_mxv(&P_periodic,rhs_per.val,n_it.rhs->val);
    }

    // Compute Nonlinear Residual and scaled version
    get_residual_norm(&n_it);
    res_norm_scaled = n_it.res_norm/sqrt(n_it.rhs->row);

    // Compute Energies
    compute_LCelastic_energy(energy,n_it.sol->val,&FE,sc,cq);

    // Checking Unit Length Constraint
    compute_LCelastic_unitlength(&unitlength,n_it.sol->val,&FE,sc,cq);

    // Print Data
    printf("Residuals:\nl2-norm of Nonlinear Residual = %25.16e\n",n_it.res_norm);
    printf("               Scaled Version = %25.16e\n",res_norm_scaled);
    printf("\nL2-norm of Update             = %25.16e\n",n_it.update_norm);
    printf("\nEnergies & Length of Director:\nTotal Energy: %25.16e\n",energy[0]);
    printf("Splay:        %25.16e\nTwist:        %25.16e\nBend:         %25.16e\n\n",energy[1],energy[2],energy[3]);
    printf("L2-norm of n: %25.16e\n\n",unitlength);

    // Check for Convergence
    newton_stop = check_newton_convergence(&n_it);

    // Compute estimator
    LCerror_estimator(errest,n_it.sol->val,&FE,sc,cq);

    if (inparam.print_level > 3) {
      // Solution at each timestep
      sprintf(solout,"output/solution_newt%03d.vtu",n_it.current_step);
      dump_blocksol_vtk(solout,varname,sc,&FE,n_it.sol->val);
      // Dump estimator into vtk
      sprintf(estout,"output/errest_newt%03d.vtu",n_it.current_step);
      dump_sol_vtk(estout,estname[0],sc,&FE_est,errest);
    }
  }

  clock_t clk_newton_end = clock();
  printf(" --> elapsed CPU time for Newton = %f seconds.\n\n",(REAL)
         (clk_newton_end-clk_newton_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  /********************* Compute Energies and Plots ****************************/
  if(inparam.print_level > 3) {
    create_pvd("output/solution.pvd",n_it.current_step+1,"solution_newt","timestep");
    create_pvd("output/errest.pvd",n_it.current_step+1,"errest_newt","timestep");
  }
  /************ Free All the Arrays ***********************************************************/
  // Arrays
  if(energy) free(energy);

  // FE Spaces
  free_fespace(&FE_nx);
  free_fespace(&FE_ny);
  free_fespace(&FE_nz);
  free_fespace(&FE_lam);
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

  // Free Periodic matrices
  if(periodic_flag) {
    bdcsr_free(&Jac_per);
    bdcsr_free(&P_periodic);
    dvec_free(&rhs_per);
    dvec_free(&update_per);
  }
  /*******************************************************************************************/
  clock_t clk_overall_end = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",
         (REAL) (clk_overall_end-clk_overall_start)/CLOCKS_PER_SEC);
  return 0;
}  /* End of Program */
/*********************************************************************************************/
