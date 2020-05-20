/*! \file LCelastic.c
 *
 *  Created by James Adler on 9/8/17.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program solves the pure elastic static LC equations
 *
 *      min 1/2 K1||div n||^2 + 1/2 K3<Z curl n, curl n>
 *               s.t.
 *                    n*n = 1
 *
 *        in 2D slab geometry -> n(x,y) = (n1(x,y),n2(x,y),n3(x,y))
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
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "hazmath.h"
#include "LCData.h"
#include "LCSystem.h"
/*********************************************************************************/

/****** MAIN DRIVER **************************************************************/
int main (int argc, char* argv[])
{

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

  // Open gridfile for reading
  printf("\nCreating mesh and FEM spaces:\n");
  FILE* gfid = HAZ_fopen(inparam.gridfile,"r");

  // Create the mesh (now we assume triangles in 2D or tetrahedra in 3D)
  // File types possible are 0 - old format; 1 - vtk format
  INT mesh_type = 0;
  clock_t clk_mesh_start = clock(); // Time mesh generation FE setup
  mesh_struct mesh;
  printf(" --> loading grid from file: %s\n",inparam.gridfile);
  creategrid_fread(gfid,mesh_type,&mesh);
  fclose(gfid);
  INT dim = mesh.dim;

  // Get Quadrature Nodes for the Mesh
  INT nq1d = inparam.nquad; /* Quadrature points per dimension */
  qcoordinates *cq = get_quadrature(&mesh,nq1d);

  // Get info for and create FEM spaces
  // Order of elements: 0 - P0; 1 - P1; 2 - P2; 20 - Nedlec; 30 - Raviart-Thomas
  INT order_n = 2;
  INT order_lam = 1;

  // Need Spaces for each component of the director plus lambda
  fespace FE_nx; // Director in x direction
  create_fespace(&FE_nx,&mesh,order_n);
  fespace FE_ny; // Director in y direction
  create_fespace(&FE_ny,&mesh,order_n);
  fespace FE_nz; // Director in z direction
  create_fespace(&FE_nz,&mesh,order_n);
  fespace FE_lam; // Lagrange multiplier, lambda
  create_fespace(&FE_lam,&mesh,order_lam);

 // Decide whether you want Dirichlet or Periodic boundary conditions.
   // periodic_flag = 1 means it's periodic bc
   // dirichlet_flag = 1 means it's dirichlet bc
  INT periodic_flag = 1;
  INT dirichlet_flag = 0;

// DIRICHLET BOUNDARY CONDITIONS
  if (dirichlet_flag == 1){
  // Set Dirichlet Boundaries (Dirichlet all around for n, and Neumann for lam)
    set_dirichlet_bdry(&FE_nx,&mesh,1,1);
    set_dirichlet_bdry(&FE_ny,&mesh,1,1);
    set_dirichlet_bdry(&FE_nz,&mesh,1,1);
    set_dirichlet_bdry(&FE_lam,&mesh,-10,-10);
}
  // Set Periodic Boundaries

  if (periodic_flag == 1){
    // Using function from HAZMATH, setting periodic boundaries
    set_periodic_bdry(&FE_nx,&mesh,0.0,1.0,0.0,0.0,0.0,0.0);
    set_periodic_bdry(&FE_ny,&mesh,0.0,1.0,0.0,0.0,0.0,0.0);
    set_periodic_bdry(&FE_nz,&mesh,0.0,1.0,0.0,0.0,0.0,0.0);
    set_periodic_bdry(&FE_lam,&mesh,0.0,1.0,0.0,0.0,0.0,0.0);

    // Relabel boundary flags for FE space to distinguish between x and y
    // relabel_boundary is found in LCSystem.h
    relabel_boundary(&FE_nx);
    relabel_boundary(&FE_ny);
    relabel_boundary(&FE_nz);
    relabel_boundary(&FE_lam);

    // Setting dirichlet  boundary conditions for the y bounds while the
    // x bounds are periodic.
    // The bounds (3,4) say that on the 2D square, the y lines are DBC
    // 1: x=0  2: x=1   3: y=0    4: y=1
    set_dirichlet_bdry(&FE_nx,&mesh,3,4);
    set_dirichlet_bdry(&FE_ny,&mesh,3,4);
    set_dirichlet_bdry(&FE_nz,&mesh,3,4);
    set_dirichlet_bdry(&FE_lam,&mesh,-10,-10);

  }

  // ------------------------------------------------------------------

  // Create Block System with ordering (n,lam)
  INT ndof = FE_nx.ndof + FE_ny.ndof + FE_nz.ndof +FE_lam.ndof;
  // Get Global FE Space
  block_fespace FE;
  FE.nun = 4;
  FE.ndof = ndof;
  FE.nbdof = FE_nx.nbdof + FE_ny.nbdof + FE_nz.nbdof+ FE_lam.nbdof;
  FE.nspaces = 4;
  FE.var_spaces = (fespace **) calloc(FE.nspaces,sizeof(fespace *));

  FE.var_spaces[0] = &FE_nx;
  FE.var_spaces[1] = &FE_ny;
  FE.var_spaces[2] = &FE_nz;
  FE.var_spaces[3] = &FE_lam;

  // Set Dirichlet Boundaries
  set_dirichlet_bdry_block(&FE,&mesh);

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
  printf("Number of Elements = %d\tOrder of Quadrature = %d\n",mesh.nelm,2*nq1d-1);
  printf("\n\t--- Element Type ---\n");
  printf("Director Element Type = %d\tLagrange Multiplier Element Type = %d\n",order_n,order_lam);
  printf("\n\t--- Degrees of Freedom ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nv,mesh.nedge,mesh.nface);
  printf("\t--> DOF: %d\n",FE.ndof);
  printf("\n\t--- Boundaries ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nbv,mesh.nbedge,mesh.nbface);
  printf("\t--> Boundary DOF: %d\n",FE.nbdof);
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
  blockFE_Evaluate(n_it.sol->val,initial_guess,&FE,&mesh,0.0);

  // Dump Initial Guess
  char solout[40];
  char** varname = malloc(50*FE.nspaces*sizeof(char *));
  varname[0] = "n1";
  varname[1] = "n2";
  varname[2] = "n3";
  varname[3] = "lambda";
  if (inparam.print_level > 3) {
    sprintf(solout,"output/solution_newt000.vtu");
    dump_blocksol_vtk(solout,varname,&mesh,&FE,n_it.sol->val);
  }

  // Perform initial Jacobian assembly
  // assemble_global_Jacobian(n_it.Jac_block,n_it.rhs,n_it.sol,local_assembly_LC_One_Constant_elastic,&FE,&mesh,cq,NULL,0.0);
  // The non-one constant one
  assemble_global_Jacobian(n_it.Jac_block,n_it.rhs,n_it.sol,local_assembly_LCelastic,&FE,&mesh,cq,NULL,0.0);
  // Eliminate Dirichlet boundary conditions in matrix and rhs
  eliminate_DirichletBC_blockFE_blockA(bc,&FE,&mesh,n_it.rhs,n_it.Jac_block,0.0);

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
    ny_ndof_periodic = Jac_per.blocks[5]->row;
    nz_ndof_periodic = Jac_per.blocks[10]->row;
    lam_ndof_periodic = Jac_per.blocks[15]->row;
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
  compute_LCelastic_energy(energy,n_it.sol->val,&FE,&mesh,cq);
  REAL unitlength=0.0;
  compute_LCelastic_unitlength(&unitlength,n_it.sol->val,&FE,&mesh,cq);
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
        solver_flag = block_directsolve_UMF(&Jac_per,&rhs_per,&update_per,linear_itparam.linear_print_level);
      } else { // Iterative Solver
        solver_flag = linear_solver_bdcsr_krylov(&Jac_per,&rhs_per,&update_per,&linear_itparam);
      }
      bdcsr_mxv(&P_periodic,update_per.val,n_it.update->val);
    } else {
      if(linear_itparam.linear_itsolver_type == 0) { // Direct Solver
        printf(" --> using UMFPACK's Direct Solver:\n");
        solver_flag = block_directsolve_UMF(n_it.Jac_block,n_it.rhs,n_it.update,linear_itparam.linear_print_level);
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
    get_blockupdate_norm(&n_it,&FE,&mesh,cq);

    // Update Jacobian and nonlinear residual with new solution
    // assemble_global_Jacobian(n_it.Jac_block,n_it.rhs,n_it.sol,local_assembly_LC_One_Constant_elastic,&FE,&mesh,cq,NULL,0.0);
    // The non-one constant
    assemble_global_Jacobian(n_it.Jac_block,n_it.rhs,n_it.sol,local_assembly_LCelastic,&FE,&mesh,cq,NULL,0.0);


    // Eliminate Dirichlet boundary conditions in matrix and rhs
    eliminate_DirichletBC_blockFE_blockA(bc,&FE,&mesh,n_it.rhs,n_it.Jac_block,0.0);
    // Eliminate Periodic boundary CONDITIONS
    if (periodic_flag == 1) {
      eliminate_PeriodicBC_blockFE_nonoverwrite(&P_periodic, n_it.Jac_block, n_it.rhs, &Jac_per, &rhs_per);
    }

    // Compute Nonlinear Residual and scaled version
    get_residual_norm(&n_it);
    res_norm_scaled = n_it.res_norm/sqrt(n_it.rhs->row);

    // Compute Energies
    compute_LCelastic_energy(energy,n_it.sol->val,&FE,&mesh,cq);

    // Checking Unit Length Constraint
    compute_LCelastic_unitlength(&unitlength,n_it.sol->val,&FE,&mesh,cq);

    // Print Data
    printf("Residuals:\nl2-norm of Nonlinear Residual = %25.16e\n",n_it.res_norm);
    printf("               Scaled Version = %25.16e\n",res_norm_scaled);
    printf("\nL2-norm of Update             = %25.16e\n",n_it.update_norm);
    printf("\nEnergies & Length of Director:\nTotal Energy: %25.16e\n",energy[0]);
    printf("Splay:        %25.16e\nTwist:        %25.16e\nBend:         %25.16e\n\n",energy[1],energy[2],energy[3]);
    printf("L2-norm of n: %25.16e\n\n",unitlength);

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

  /********************* Compute Energies and Plots ****************************/
  if(inparam.print_level > 3) {
    create_pvd("output/solution.pvd",n_it.current_step+1,"solution_newt","timestep");
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
  free_mesh(&mesh);

  // Strings
  if(varname) free(varname);

  // Newton
  free_newton(&n_it);

  // Free Periodic matrices
  bdcsr_free(&Jac_per);
  bdcsr_free(&P_periodic);
  dvec_free(&rhs_per);
  dvec_free(&update_per);

  /*******************************************************************************************/
  clock_t clk_overall_end = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",
         (REAL) (clk_overall_end-clk_overall_start)/CLOCKS_PER_SEC);
  return 0;
}  /* End of Program */
/*********************************************************************************************/
