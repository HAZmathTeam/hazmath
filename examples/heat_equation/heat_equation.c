/*! \file examples/heat_equation/heat_equation.c
 *
 *  Created by James Adler and Xiaozhe Hu on 10/16/16.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program solves the following PDE using finite elements
 *
 *        du/dt - div(a(x)grad(u)) = 0
 *
 *        where du/dt is discretized with Crank-Nicolson,
 *        BDF-1 (Backward Euler), or BDF-2 in 2D or 3D.
 *
 *        Along the boundary of the region, Dirichlet conditions are imposed:
 *
 *          u = 0 for P1 or P2 elements
 *
 * \note This example extends the one in basic_elliptic, while also showing
 *       how to implement a variety of time discretizations, including
 *       Backward Euler (BDF1), BDF2, and Crank-Nicolson. The timestepper
 *       struct is introduced here.
 *
 */

/************* HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
#include "heat_equation_data.h"
/***********************************************************************/

/****** MAIN DRIVER **************************************************/
int main (int argc, char* argv[])
{

  printf("\n===========================================================================\n");
  printf("Beginning Program to solve the Heat Equation.\n");
  printf("===========================================================================\n");

  /****** INITIALIZE PARAMETERS **************************************/
  // Flag for errors
  SHORT status;

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
  INT nq1d = inparam.nquad; // Quadrature points per dimension
  qcoordinates *cq = get_quadrature(&mesh,nq1d);

  // Get info for and create FEM spaces
  // Order of Elements: 0 - P0; 1 - P1; 2 - P2
  INT order = inparam.FE_type;
  fespace FE;
  create_fespace(&FE,&mesh,order);
  // Strings for printing
  char elmtype[8];
  sprintf(elmtype,"P%d",order);

  // Set Dirichlet Boundaries
  // Assume physical boundaries are Dirichlet
  // The mesh is set up so that flag values 1-16 are Dirichlet and 17-32 are Neumann
  set_dirichlet_bdry(&FE,&mesh,1,1);

  // Dump some of the data
  if(inparam.print_level > 5 && inparam.output_dir!=NULL) {
    // FE space
    char varu[10];
    char dir[20];
    sprintf(dir,"output");
    sprintf(varu,"u");
    dump_fespace(&FE,varu,dir);

    // Mesh
    char* namevtk = "output/mesh.vtu";
    dump_mesh_vtk(namevtk,&mesh);
  }

  clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
  printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
         (REAL) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
  /*******************************************************************/

  printf("***********************************************************************************\n");
  printf("\t--- %d-dimensional grid ---\n",dim);
  printf("Number of Elements = %d\tElement Type = %s\tOrder of Quadrature = %d\n",mesh.nelm,elmtype,2*nq1d-1);
  printf("\n\t--- Degrees of Freedom ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nv,mesh.nedge,mesh.nface);
  printf("\t--> DOF: %d\n",FE.ndof);
  printf("\n\t--- Boundaries ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nbv,mesh.nbedge,mesh.nbface);
  printf("\t--> Boundary DOF: %d\n",FE.nbdof);
  printf("***********************************************************************************\n\n");

  /*** Assemble the matrix and right hand side ***********************/
  printf("Assembling the matrix and right-hand side:\n");
  clock_t clk_assembly_start = clock();

  // Allocate the right-hand side and declare the csr matrix
  dvector b;
  dCSRmat A;
  dCSRmat M;

  // Assemble the matrix without BC
  // Diffusion block
  assemble_global(&A,&b,assemble_DuDv_local,&FE,&mesh,cq,myrhs,
                  diffusion_coeff,0.0);

  // Time-Derivative block
  assemble_global(&M,NULL,assemble_mass_local,&FE,&mesh,cq,NULL,
                  one_coeff_scal,0.0);

  // Create Time Operator (one with BC and one without)
  // We solve operators of the form du/dt + L(u) = f
  // Thus our discretization would be of the form d/dt Mu + Au = b
  // Note that since this is linear, L(u) = Au, so we set Ldata to A
  timestepper time_stepper;
  // The RHS is not time-dependent, so this makes our life a bit easier
  // The code can skip some steps and not build <f,v> each time step
  INT rhs_timedep = 0;
  initialize_timestepper(&time_stepper,&inparam,rhs_timedep,b.row);
  time_stepper.A = &A;
  time_stepper.Ldata=&A;
  time_stepper.M = &M;
  // Create the time operator.  Note this is the first time we create it,
  // so the 2nd argument is 1, and we need to copy the operator to apply
  // the boundary conditions correctly at each time step, so the 3rd argument
  // is also 1.
  get_timeoperator(&time_stepper,1,1);

  clock_t clk_assembly_end = clock();
  printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL)
         (clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);
  /*******************************************************************/

  /**************** Solve ********************************************/

  // Create Solution Vector
  dvector sol = dvec_create(FE.ndof);
  dvector exact_sol = dvec_create(FE.ndof);
  REAL current_time = 0.0;
  if(dim==2) FE_Evaluate(exact_sol.val,exactsol2D,&FE,&mesh,current_time);
  if(dim==3) FE_Evaluate(exact_sol.val,exactsol3D,&FE,&mesh,current_time);

  // Set parameters for linear iterative methods
  linear_itsolver_param linear_itparam;
  param_linear_solver_init(&linear_itparam);
  param_linear_solver_set(&linear_itparam, &inparam);
  INT solver_flag=-20;

  // For direct solver we can factorize the matrix ahead of time and not each time step
  void* Numeric = NULL;

  // Set parameters for algebriac multigrid methods
  AMG_param amgparam;
  param_amg_init(&amgparam);
  param_amg_set(&amgparam, &inparam);
  param_amg_print(&amgparam);

  // Get Initial Conditions
  if(dim==2) FE_Evaluate(sol.val,initial_conditions2D,&FE,&mesh,current_time);
  if(dim==3) FE_Evaluate(sol.val,initial_conditions3D,&FE,&mesh,current_time);
  time_stepper.sol = &sol;

  // Dump Solution
  char solout[40];
  char exactout[40];
  if (inparam.output_dir!=NULL) {
    sprintf(solout,"output/solution_ts000.vtu");
    dump_sol_vtk(solout,"u",&mesh,&FE,time_stepper.sol->val);
  }

  // Store current RHS
  time_stepper.rhs = &b;

  // Compute initial errors and norms
  REAL* uerr = (REAL *) calloc(time_stepper.tsteps+1,sizeof(REAL));
  if(dim==2) uerr[0] = L2error(time_stepper.sol->val,exactsol2D,&FE,&mesh,cq,time_stepper.time);
  if(dim==3) uerr[0] = L2error(time_stepper.sol->val,exactsol3D,&FE,&mesh,cq,time_stepper.time);
  REAL* unorm = (REAL *) calloc(time_stepper.tsteps+1,sizeof(REAL));
  unorm[0] = L2norm(time_stepper.sol->val,&FE,&mesh,cq);
  REAL* utnorm = (REAL *) calloc(time_stepper.tsteps+1,sizeof(REAL));
  utnorm[0] = L2norm(exact_sol.val,&FE,&mesh,cq);

  printf("Performing %d Time Steps with step size dt = %1.3f\n",time_stepper.tsteps,time_stepper.dt);
  printf("--------------------------------------------------------------\n\n");
  printf("============================\n");
  printf("Time Step %d: Time = %1.8f\n",time_stepper.current_step,time_stepper.time);
  printf("============================\n");
  printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
  printf("L2 Norm of u            = %26.13e\n",unorm[0]);
  printf("L2 Norm of exact u       = %26.13e\n",utnorm[0]);
  printf("L2 Norm of u error      = %26.13e\n",uerr[0]);
  printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n\n");

  clock_t clk_timeloop_start = clock();
  // Begin Timestepping Loop
  INT j; // Time step counter
  for(j=0;j<time_stepper.tsteps;j++) {
    clock_t clk_timestep_start = clock();

    // Update Time Step Data (includes time, counters, solution, and rhs)
    update_timestep(&time_stepper);

    printf("============================\n");
    printf("Time Step %d: Time = %1.3f\n",time_stepper.current_step,time_stepper.time);
    printf("============================\n");

    // Recompute RHS if it's time-dependent
    if(time_stepper.rhs_timedep) {
      assemble_global_RHS(time_stepper.rhs,&FE,&mesh,cq,myrhs,time_stepper.time);
    }

    // Update RHS
    update_time_rhs(&time_stepper);

    // For first time step eliminate boundary conditions in matrix and rhs
    if(j==0) {
      eliminate_DirichletBC(bc,&FE,&mesh,time_stepper.rhs_time,time_stepper.At,time_stepper.time);
      // If Direct Solver used only factorize once
      if(linear_itparam.linear_itsolver_type == 0) { // Direct Solver
#if WITH_SUITESPARSE
        printf(" --> using UMFPACK's Direct Solver: factorization \n");
#else
        printf(" --> using HAZMATH's Direct Solver: factorization \n");
	//        error_extlib(255,__FUNCTION__,"SuiteSparse");
	//        return 0;
#endif
        Numeric = factorize_HAZ(time_stepper.At,linear_itparam.linear_print_level);
      }
    } else {
      eliminate_DirichletBC_RHS(bc,&FE,&mesh,time_stepper.rhs_time,time_stepper.At_noBC,time_stepper.time);
    }

    // Solve
    clock_t clk_solve_start = clock();
    if(linear_itparam.linear_itsolver_type == 0) { // Direct Solver
#if WITH_SUITESPARSE
      printf(" --> using UMFPACK's Direct Solver: solve\n");
#else
      printf(" --> using HAZMATH's Direct Solver: solve\n");
      //      error_extlib(255,__FUNCTION__,"SuiteSparse");
      //      return 0;
#endif
      solver_flag = solve_HAZ(time_stepper.At,time_stepper.rhs_time,time_stepper.sol,Numeric,linear_itparam.linear_print_level);
    } else { // Iterative Solver
      // Use AMG as iterative solver
      if (linear_itparam.linear_itsolver_type == SOLVER_AMG){
        solver_flag = linear_solver_amg(time_stepper.At,time_stepper.rhs_time,time_stepper.sol, &amgparam);
      } else { // Use Krylov Iterative Solver
        // Determine Preconditioner
        switch (linear_itparam.linear_precond_type) {
        case PREC_DIAG:  // Diagonal Preconditioner
          solver_flag = linear_solver_dcsr_krylov_diag(time_stepper.At,time_stepper.rhs_time,time_stepper.sol,&linear_itparam);
          break;
        case PREC_AMG:  // AMG preconditioner
          solver_flag = linear_solver_dcsr_krylov_amg(time_stepper.At,time_stepper.rhs_time,time_stepper.sol, &linear_itparam, &amgparam);
          break;
        default:  // No Preconditioner
          solver_flag = linear_solver_dcsr_krylov(time_stepper.At,time_stepper.rhs_time,time_stepper.sol,&linear_itparam);
          break;
        }
      }
    }

    // Error Check
    if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n", solver_flag);

    clock_t clk_solve_end = clock();
    printf("Elapsed CPU Time for Solve = %f seconds.\n\n",(REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);

    clock_t clk_timestep_end = clock();
    printf("Elapsed CPU Time for Time Step = %f seconds.\n\n",(REAL) (clk_timestep_end-clk_timestep_start)/CLOCKS_PER_SEC);

    /**************** Compute Errors if you have exact solution *******/
    clock_t clk_error_start = clock();

    if(dim==2) uerr[j+1] = L2error(time_stepper.sol->val,exactsol2D,&FE,&mesh,cq,time_stepper.time);
    if(dim==3) uerr[j+1] = L2error(time_stepper.sol->val,exactsol3D,&FE,&mesh,cq,time_stepper.time);
    unorm[j+1] = L2norm(time_stepper.sol->val,&FE,&mesh,cq);
    if(dim==2) FE_Evaluate(exact_sol.val,exactsol2D,&FE,&mesh,time_stepper.time);
    if(dim==3) FE_Evaluate(exact_sol.val,exactsol3D,&FE,&mesh,time_stepper.time);
    utnorm[j+1] = L2norm(exact_sol.val,&FE,&mesh,cq);
    clock_t clk_error_end = clock();
    printf("Elapsed CPU time for getting errors = %lf seconds.\n\n",(REAL)
           (clk_error_end-clk_error_start)/CLOCKS_PER_SEC);
    printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
    printf("L2 Norm of u            = %26.13e\n",unorm[j+1]);
    printf("L2 Norm of exact u       = %26.13e\n",utnorm[j+1]);
    printf("L2 Norm of u error      = %26.13e\n",uerr[j+1]);
    printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
    /*******************************************************************/

    if (inparam.output_dir!=NULL) {
      sprintf(solout,"output/solution_ts%03d.vtu",time_stepper.current_step);
      dump_sol_vtk(solout,"u",&mesh,&FE,time_stepper.sol->val);
      sprintf(exactout,"output/exact_solution_ts%03d.vtu",time_stepper.current_step);
      dump_sol_vtk(exactout,"ut",&mesh,&FE,exact_sol.val);
    }
    printf("\n");
  } // End Timestepping Loop
  printf("----------------------- Timestepping Complete ---------------------------------------\n\n");
  clock_t clk_timeloop_end = clock();
  printf("Elapsed CPU Time ALL Time Steps = %f seconds.\n\n",
         (REAL) (clk_timeloop_end-clk_timeloop_start)/CLOCKS_PER_SEC);
  /*******************************************************************/

  /******** Summary Print ********************************************/
  printf("Summary of Timestepping\n");
  printf("Time Step\tTime\t\t\t||u||\t\t\t\t||u_exact||\t\t\t||error||\n\n");
  for(j=0;j<=time_stepper.tsteps;j++) {
    printf("%02d\t\t%f\t%25.16e\t%25.16e\t%25.16e\n",j,j*time_stepper.dt,unorm[j],utnorm[j],uerr[j]);
  }

  // Combine all timestep vtks in one file
  if (inparam.output_dir!=NULL) {
    create_pvd("output/solution.pvd",time_stepper.tsteps+1,"solution_ts","timestep");
  }

  /******** Free All the Arrays **************************************/
  if(unorm) free(unorm);
  if(utnorm) free(utnorm);
  if(uerr) free(uerr);
  dvec_free(&exact_sol);
  free_timestepper(&time_stepper);
  free_fespace(&FE);
  if(cq) {
    free_qcoords(cq);
    free(cq);
    cq = NULL;
  }
  free_mesh(&mesh);
  //#if WITH_SUITESPARSE
  if (Numeric) hazmath_free_numeric(&Numeric);
  //#endif

  /*******************************************************************/

  clock_t clk_overall_end = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",
         (REAL) (clk_overall_end-clk_overall_start)/CLOCKS_PER_SEC);
  return 0;

}	/* End of Program */
/*******************************************************************/
