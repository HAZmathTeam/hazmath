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
//#include "../data/2D/lc_basic_twist_data.h"
#include "../data/2D/lc_harmonic_solution.h"
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
  // We will also create the nested iteration struct here which will contain the simplicial complex for adaptive refinement and marking
  mesh_struct mesh;
  nested_it ni;
  dvector sol; // Solution on next level mesh
  initialize_ni(&ni);
  INT dim = inparam.spatial_dim;                 // dimension of computational domain
  INT init_ref_levels = inparam.refinement_levels; // refinement levels of initial (coarsest uniform) mesh
  INT mesh_ref_type = inparam.refinement_type;     // refinement type (>10 uniform or <10 other)
  INT set_bndry_codes = inparam.boundary_codes;    // set flags for the boundary DoF (1-16 are Dirichlet) (this isn't used??)
  get_initial_mesh_ni(&ni,dim,init_ref_levels);
  mesh = ni.mesh[0];

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

  // Dump the mesh to file
  char meshout[40];
  if (inparam.print_level > 3) {
    sprintf(meshout,"output/mesh00.vtu");
    dump_mesh_vtk(meshout,&mesh);
  }

  // Get Quadrature Nodes for the initial Mesh
  INT nq1d = inparam.nquad; /* Quadrature points per dimension */
  qcoordinates *cq = get_quadrature(&mesh,nq1d);
  ni.cq = cq;

  // Get info for and create FEM spaces
  // Order of elements: 0 - P0; 1 - P1; 2 - P2; 20 - Nedlec; 30 - Raviart-Thomas
  INT order_n = 2;
  INT order_lam = 1;
  fespace FE_nx; // Director in x direction
  fespace FE_ny; // Director in y direction
  fespace FE_nz; // Director in z direction
  fespace FE_lam; // Lagrange multiplier, lambda
  block_fespace FE; // Global block FE space
  block_dCSRmat P_periodic; // Periodic Boundary Matrix (if needed)
  setup_FEspaces(order_n,order_lam,&FE_nx,&FE_ny,&FE_nz,&FE_lam,&FE,&P_periodic,&mesh,dim);
  INT ndof = FE.ndof; // Number of degrees of freedom
  INT nspaces = FE.nspaces; // Number of spaces in the block FE space
  INT nun = FE.nun; // Number of unknowns in the block FE space
  ni.FE[0] = FE;

  clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
  printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
         (REAL) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
  /*******************************************************************************/

  printf("***********************************************************************************\n");
  printf("%d Dimensional Problem\n",mesh.dim);
  printf("Number of Elements = %d\tOrder of Quadrature = %d\n",mesh.nelm,2*nq1d-1);
  printf("\n\t--- Element Type ---\n");
  printf("Director Element Type = %d\tLagrange Multiplier Element Type = %d\n",order_n,order_lam);
  printf("\n\t--- Degrees of Freedom ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nv,mesh.nedge,mesh.nface);
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
  newton newt;
  newt.isblock = 1; // All matrices are blocks
  initialize_newton(&newt,&inparam,ndof,FE.nspaces);

  // Change Newton weighting if desired
  newt.step_length = 1.0;

  // Allocate the Initial Guess
  srand(803087000);
  blockFE_Evaluate(newt.sol->val,initial_guess,&FE,&mesh,0.0);

  // Error Estimator Stuff
  REAL* errest = (REAL *) calloc(mesh.nelm,sizeof(REAL));
  // Plotting stuff
  char estout[40];
  char** estname = malloc(50*sizeof(char *));
  estname[0] = "err_est";
  fespace FE_est; // P0 space to store estimator for plotting
  create_fespace(&FE_est,&mesh,0);
  // Get estimate of initial guess
  LCerror_estimator(errest,newt.sol->val,&FE,&mesh,cq);


  // Dump Initial Guess
  char solout[40];
  char** varname = malloc(50*FE.nspaces*sizeof(char *));
  varname[0] = "n1";
  varname[1] = "n2";
  varname[2] = "n3";
  varname[3] = "lambda";
  if (inparam.print_level > 3) {
    sprintf(solout,"output/solution_newt000.vtu");
    dump_blocksol_vtk(solout,varname,&mesh,&FE,newt.sol->val);
    // Dump estimator into vtk
    sprintf(estout,"output/errest_newt000.vtu");
    dump_sol_vtk(estout,estname[0],&mesh,&FE_est,errest);
  }

  // Perform initial Jacobian assembly
  assemble_global_Jacobian(newt.Jac_block,newt.rhs,newt.sol,local_assembly_LCelastic,&FE,&mesh,cq,NULL,0.0);
  // Eliminate Dirichlet boundary conditions in matrix and rhs
  eliminate_DirichletBC_blockFE_blockA(bc,&FE,&mesh,newt.rhs,newt.Jac_block,0.0);

  // Deal with Periodic Boundary conditions, but not mess up newt.Jac_block
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
    eliminate_PeriodicBC_blockFE_nonoverwrite(&P_periodic, newt.Jac_block, newt.rhs, &Jac_per, &rhs_per);

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
  get_residual_norm(&newt);
  REAL res_norm_scaled = newt.res_norm/sqrt(newt.rhs->row);
  printf("\nInitial Residuals:\nl2-norm of Nonlinear Residual = %25.16e\n",newt.res_norm);
  printf("               Scaled Version = %25.16e\n",res_norm_scaled);

  // Compute Initial Energy and Length of Director
  REAL* energy = (REAL *) calloc(4,sizeof(REAL));
  compute_LCelastic_energy(energy,newt.sol->val,&FE,&mesh,cq);
  REAL unitlength=0.0;
  compute_LCelastic_unitlength(&unitlength,newt.sol->val,&FE,&mesh,cq);
  printf("\nInitial Energies & Length of Director:\nTotal Energy: %25.16e\n",energy[0]);
  printf("Splay:        %25.16e\nTwist:        %25.16e\nBend:         %25.16e\n\n",energy[1],energy[2],energy[3]);
  printf("L2 norm of n: %25.16e\n\n",unitlength);

  // Check Convergence before starting
  INT newton_stop=0;
  if(res_norm_scaled < newt.tol) {
    printf("The initial nonlinear residual is below the tolerance.  Not doing any stepping.\n\n");
    newton_stop=1;
  }

  // Perform Newton Steps
  while(!newton_stop) {

    // Update Newton Step Data
    update_newtonstep(&newt);
    printf("==================================\n");
    printf("\t\tNewton step %d\t\t\n",newt.current_step);
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
      bdcsr_mxv(&P_periodic,update_per.val,newt.update->val);
    } else {
      if(linear_itparam.linear_itsolver_type == 0) { // Direct Solver
        printf(" --> using UMFPACK's Direct Solver:\n");
        solver_flag = block_directsolve_HAZ(newt.Jac_block,newt.rhs,newt.update,linear_itparam.linear_print_level);
      } else { // Iterative Solver
        solver_flag = linear_solver_bdcsr_krylov(newt.Jac_block,newt.rhs,newt.update,&linear_itparam);
      }
    }

    // Error Check
    if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n", solver_flag);

    clock_t clk_solve_end = clock();
    printf("Elapsed CPU Time for Solve = %f seconds.\n\n",(REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);

    // Update Newton Solution
    update_sol_newton(&newt);

    // Get norm of update
    get_blockupdate_norm(&newt,&FE,&mesh,cq);

    // Update Jacobian and nonlinear residual with new solution
    assemble_global_Jacobian(newt.Jac_block,newt.rhs,newt.sol,local_assembly_LCelastic,&FE,&mesh,cq,NULL,0.0);

    // Eliminate Dirichlet boundary conditions in matrix and rhs
    eliminate_DirichletBC_blockFE_blockA(bc,&FE,&mesh,newt.rhs,newt.Jac_block,0.0);
    // Eliminate Periodic boundary CONDITIONS
    if (periodic_flag == 1) {
      eliminate_PeriodicBC_blockFE_nonoverwrite(&P_periodic, newt.Jac_block, newt.rhs, &Jac_per, &rhs_per);
      // Blow nonlinear residual back up to calculate norm
      bdcsr_mxv(&P_periodic,rhs_per.val,newt.rhs->val);
    }

    // Compute Nonlinear Residual and scaled version
    get_residual_norm(&newt);
    res_norm_scaled = newt.res_norm/sqrt(newt.rhs->row);

    // Compute Energies
    compute_LCelastic_energy(energy,newt.sol->val,&FE,&mesh,cq);

    // Checking Unit Length Constraint
    compute_LCelastic_unitlength(&unitlength,newt.sol->val,&FE,&mesh,cq);

    // Print Data
    printf("Residuals:\nl2-norm of Nonlinear Residual = %25.16e\n",newt.res_norm);
    printf("               Scaled Version = %25.16e\n",res_norm_scaled);
    printf("\nL2-norm of Update             = %25.16e\n",newt.update_norm);
    printf("\nEnergies & Length of Director:\nTotal Energy: %25.16e\n",energy[0]);
    printf("Splay:        %25.16e\nTwist:        %25.16e\nBend:         %25.16e\n\n",energy[1],energy[2],energy[3]);
    printf("L2-norm of n: %25.16e\n\n",unitlength);

    // Check for Convergence
    newton_stop = check_newton_convergence(&newt);

    // Compute estimator
    LCerror_estimator(errest,newt.sol->val,&FE,&mesh,cq);

    if (inparam.print_level > 3) {
      // Solution at each timestep
      sprintf(solout,"output/solution_newt%03d.vtu",newt.current_step);
      dump_blocksol_vtk(solout,varname,&mesh,&FE,newt.sol->val);
      // Dump estimator into vtk
      sprintf(estout,"output/errest_newt%03d.vtu",newt.current_step);
      dump_sol_vtk(estout,estname[0],&mesh,&FE_est,errest);
    }
  }

  clock_t clk_newton_end = clock();
  printf(" --> elapsed CPU time for Newton = %f seconds.\n\n",(REAL)
         (clk_newton_end-clk_newton_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  /********************* Adaptive refinement ****************************/
  //-------------------------------------------------
  // Added by XH for mesh refinement
  // Moved to HAZMATH nested_iteration.c and marking.c
  //-------------------------------------------------

  // First we refine the mesh and the FE spaces, but keep the previous ones to interpolate the solution.
 
  // Refine the mesh to get to next level and update quadrature on the new mesh
  ni.mark_type=1;
  ni.mark_param=0.0;
  ni.err_est=errest;
  ni_refine_mesh(&ni);
  // Dump Mesh for plotting
  if (inparam.print_level > 3) {
    sprintf(meshout,"output/mesh%02d.vtu",1);
    dump_mesh_vtk(meshout,ni.mesh);
  }

  // Update FE Spaces on new mesh
  // This is problem dependent so must be done by user
  setup_FEspaces(order_n,order_lam,&FE_nx,&FE_ny,&FE_nz,&FE_lam,&FE,&P_periodic,&mesh,dim);
  dvec_alloc(FE.ndof,&sol);

  // icsr_print_matlab(stdout,ni.sc_all[0]->parent_v);
  // // old vertices come first ALWAYS

  // Update solution on new mesh
  next_update_sol(&ni,newt.sol,&sol,&FE,&mesh);

  // Update Local Variables to current data
  mesh = ni.mesh[0];
  cq = ni.cq;
  ndof = FE.ndof; // Number of degrees of freedom
  nspaces = FE.nspaces; // Number of spaces in the block FE space
  nun = FE.nun; // Number of unknowns in the block FE space
  ni.FE[0] = FE;


  //-------------------------------------------------

  /********************* Compute Energies and Plots ****************************/
  if(inparam.print_level > 3) {
    create_pvd("output/solution.pvd",newt.current_step+1,"solution_newt","timestep");
    create_pvd("output/errest.pvd",newt.current_step+1,"errest_newt","timestep");
  }
  /************ Free All the Arrays ***********************************************************/
  // Arrays
  if(energy) free(energy);

  // FE Spaces
  free_fespace(&FE_nx);
  free_fespace(&FE_ny);
  free_fespace(&FE_nz);
  free_fespace(&FE_lam);
  //free_blockfespace(&FE);

  // Quadrature
  // if(cq){
  //   free_qcoords(cq);
  //   free(cq);
  //   cq=NULL;
  // }

  // // Mesh
  // free_mesh(&mesh);

  // Strings
  if(varname) free(varname);

  // Nested Iteration
  free_ni(&ni);

  // Newton
  free_newton(&newt);

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
