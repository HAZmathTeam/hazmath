/*! \file examples/DarcyFlow/DarcyFlow.c
 *
 *  Created by Adler, Hu, Zikatanov on 8/30/16.
 *  Copyright 2016_HAZMATH__. All rights reserved.
 *
 * \brief This program creates the Discretization of the Darcy Flow
 *        System using Mixed Finite Elements
 *
 *          Ss dh/dt - div K(x)grad h = W
 *
 *        ->  Ss dh/dt - div q = W
 *            K^(-1) q - grad h = 0
 *
 *       For now we assume K is a diagonal matrix.  The weak form is:
 *
 *        <K^(-1) q, r> + <h, div r> = <g,r*n>_boundary
 *        -<Ss dh/dt, v> + <div q, v> = -<W,v>
 *
 *       where <*,*> is the L2 inner product and we use Crank-Nicolson
 *       or BDF1 or 2 for the time stepping
 *
 * \note This example shows how to build your own bilinear form for a
 *       system.  The forms are found in Darcy_assemble.h and
 *       coefficients, right hand side etc problem data is found in
 *       Darcy_functions.h
 *
 *
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "global.h"
#include "hazmath.h"
#include "Darcy_functions.h"
#include "Darcy_assemble.h"
/*********************************************************************************/

/****** MAIN DRIVER **************************************************************/
int main (int argc, char* argv[]) 
{

  printf("\n===========================================================================\n");
  printf("\nBeginning Program to solve Darcy Flow eqn by RT0-P0 mixed FE method\n");
  printf("\n===========================================================================\n");

  /****** INITIALIZE PARAMETERS **************************************************/
  // Loop Indices
  INT i,j,ii;

  // Overall Timing
  clock_t clk_overall_start = clock();

  // Set Parameters from Reading in Input File
  input_param inparam;
  param_input_init(&inparam);
  param_input("./input.dat", &inparam);
  
  // Open gridfile for reading
  printf("\nCreating mesh and FEM spaces:\n");
  FILE* gfid = HAZ_fopen(inparam.gridfile,"r");

  // Dimension is needed for all this to work

  // Create the mesh (now we assume triangles in 2D or tetrahedra in 3D)
  // File types possible are 0 - old format; 1 - vtk format
  clock_t clk_mesh_start = clock(); // Time mesh generation FE setup
  INT mesh_type = 0;
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
  // Order of Elements: 0 - P0; 1 - P1; 2 - P2; 20 - Nedelec; 30 - Raviart-Thomas
  INT order_h = 0;
  INT order_q = 30;
  fespace FE_h;
  create_fespace(&FE_h,&mesh,order_h);
  fespace FE_q;
  create_fespace(&FE_q,&mesh,order_q);

  // Time stepping parameters
  block_timestepper time_stepper;

  if(inparam.print_level > 3) {
    // Dump Mesh and FE data
    char varh[10];
    char dir[20];
    sprintf(dir,"output");
    sprintf(varh,"h");
    dump_fespace(&FE_h,varh,dir);
    FILE* fid = fopen("output/coords.dat","w");
    dump_coords(fid,mesh.cv);
    fclose(fid);
    
    char* namevtk = "output/mesh.vtu";
    dump_mesh_vtk(namevtk,&mesh);
  }

  // Set ESSENTIAL (Dirichlet) Boundaries. For mixed formulation
  // essential = Neumann.
  set_dirichlet_bdry(&FE_q,&mesh,(INT )MARKER_NEUMANN,(INT )MARKER_ROBIN);
  /*(ltz) unclear why weneed this here, it should be gone
    for(i=0;i<FE_q.ndof;i++) {
    if(FE_q.dirichlet[i]==1 && (mesh.f_mid[i*dim+2]!=1 && mesh.f_mid[i*dim+2]!=0)) {
    FE_q.dirichlet[i] = 0;
    }
    }
  */
  // Create Block System with ordering (q,h)
  INT ndof = FE_q.ndof + FE_h.ndof;
  // Get Global FE Space
  block_fespace FE;
  FE.ndof = FE_q.ndof + FE_h.ndof;
  FE.nbdof = FE_q.nbdof + FE_h.nbdof;
  FE.nspaces = 2;
  FE.nun = dim+1;
  FE.var_spaces = (fespace **) calloc(2,sizeof(fespace *));
  FE.var_spaces[0] = &FE_q;
  FE.var_spaces[1] = &FE_h;
  set_dirichlet_bdry_block(&FE,&mesh);

  clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
  printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
         (REAL) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
  /*******************************************************************************/

  printf("***********************************************************************************\n");
  printf("Element Type: RT for q and P0 for h.\n");
  printf("Number of Elements = %d\tOrder of Quadrature = %d\n",mesh.nelm,2*nq1d-1);
  printf("\n\t--- Degrees of Freedom ---\n");
  printf("Vertices: %-7d  Edges: %-7d  Faces: %-7d",mesh.nv,mesh.nedge,mesh.nface);
  printf("  --> DOF: %d (%d for q and %d for h).\n",FE.ndof,FE_q.ndof,FE_h.ndof);
  printf("\n\t--- Boundaries ---\n");
  printf("Vertices: %-7d  Edges: %-7d  Faces: %-7d",mesh.nbv,mesh.nbedge,mesh.nbface);
  printf("  --> Boundary DOF: %d (%d for q and %d for h).\n",FE.nbdof,FE_q.nbdof,FE_h.nbdof);
  printf("***********************************************************************************\n\n");

  /*** Assemble the matrix and right hand side *******************************/
  /* Here we assemble the discrete system:
   *  For now we assume K is a diagonal matrix.  The weak form in steady state is:
   *
   *  <K^(-1) q, r> + <h, div r> = <g,r*n>_boundary
   *  -<Ss dh/dt, v> + <div q, v> = -<W,v>
   */

  printf("Assembling the matrix and right-hand side:\n");
  clock_t clk_assembly_start = clock();

  // Assemble the matrices

  // Get RHS ready
  dvector b;
  dvector b_bdry = dvec_create(FE_q.ndof);
  dvec_set(FE_q.ndof,&b_bdry,0.0);

  // Steady-State Part
  // A is a 2 x 2 block system
  // A[0]=qq-block  A[1]=qh-block
  // A[2]=hq-block  A[3]=hh-block
  block_dCSRmat A;
  bdcsr_alloc(2,2,&A);

  // All terms but boundary integral
  // Note that this allocates every block based on FEM stencil
  assemble_global_block(&A,&b,steady_state_Darcy,steady_state_Darcy_RHS,&FE,&mesh,cq,source,0.0);

  // Boundary Integral <g,r*n>_boundary
  // Flag is which boundary you want to compute this
  INT flag = 1; //this flag should be irrelevant here. but it is used 
  assemble_global_RHS_face(&b_bdry,NULL,steady_state_Darcy_bdryRHS,&FE_q,&mesh,cq,myg,0.0,flag);

  // Add RHS vectors together
  for(i=0;i<FE_q.ndof;i++) {
    b.val[i] += b_bdry.val[i];
  }

  // Create Global Mass Matrix for time-stepping
  block_dCSRmat M;
  bdcsr_alloc(2,2,&M);

  // Since M matrix is mostly zeros except for (1,1) block, we'll set the remaining blocks to
  // a zero matrix.  We are not calling a block assembly so this is necessary.
  dcsr_set_zeromatrix(M.blocks[0], A.blocks[0]->row,A.blocks[0]->col,1);
  dcsr_set_zeromatrix(M.blocks[1], A.blocks[1]->row,A.blocks[1]->col,1);
  dcsr_set_zeromatrix(M.blocks[2], A.blocks[2]->row,A.blocks[2]->col,1);

  // Get Mass Matrix for h
  assemble_global(M.blocks[3],NULL,assemble_mass_local,&FE_h,&mesh,cq,NULL,storage_coeff,0.0);

  // Create Time Operator (one with BCOND and one without)
  initialize_blktimestepper(&time_stepper,&inparam,0,FE.ndof,2);
  time_stepper.A = &A;
  time_stepper.M = &M;
  time_stepper.Ldata=&A;
  get_blktimeoperator(&time_stepper,1,1);

  //-----------------------------------------------
  // prepare for preconditioners
  // prepare diagonal blocks
  dCSRmat *A_diag;
  A_diag = (dCSRmat *)calloc(2, sizeof(dCSRmat));
  
  dCSRmat BTB;
  dvector el_vol;
  //-----------------------------------------------
  
  clock_t clk_assembly_end = clock();
  printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL)
         (clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  /**************** Solve ********************************************************************/

  dvector sol = dvec_create(ndof);

  // Set parameters for linear iterative methods
  linear_itsolver_param linear_itparam;
  param_linear_solver_init(&linear_itparam);
  param_linear_solver_set(&linear_itparam, &inparam);

  // Set parameters for AMG
  AMG_param amgparam;
  param_amg_init(&amgparam);
  param_amg_set(&amgparam, &inparam);
  param_amg_print(&amgparam);

  // Solver flag
  INT solver_flag;

  // Get Initial Conditions
  blockFE_Evaluate(sol.val,initial_conditions,&FE,&mesh,0.0);
  time_stepper.sol = &sol;

  // Store current RHS
  time_stepper.rhs = &b;

  // Dump Initial Condition
  char solout[40];
  char** varname = malloc(50*FE.nspaces*sizeof(char *));
  varname[0] = "q";
  varname[1] = "h";
  if (inparam.output_dir!=NULL) {
    sprintf(solout,"output/solution_ts000.vtu");
    dump_blocksol_vtk(solout,varname,&mesh,&FE,sol.val);

    //if(dim==3) print_matlab_vector_field(&v_ux,&v_uy,&v_uz,&FE_ux);
  }

  // Begin Timestepping Loop
  clock_t clk_timeloop_start = clock();
  for(j=0;j<time_stepper.tsteps;j++) {
    clock_t clk_timestep_start = clock();
    // Update Time Step Data (includes time, counters, solution, and rhs)
    update_blktimestep(&time_stepper);

    // Recompute RHS if it's time-dependent
    if(time_stepper.rhs_timedep) {
      // Interior
      assemble_global_RHS_block(time_stepper.rhs,steady_state_Darcy_RHS,&FE,&mesh,cq,source,time_stepper.time);
      // Boundary Integral <g,r*n>_boundary
      assemble_global_RHS_face(&b_bdry,NULL,steady_state_Darcy_bdryRHS,&FE_q,&mesh,cq,myg,0.0,flag);
      // Add RHS vectors together
      for(i=0;i<FE_q.ndof;i++) {
        time_stepper.rhs->val[i] += b_bdry.val[i];
      }
    }

    // Update RHS
    update_blktime_rhs(&time_stepper);

    // For first time step eliminate boundary conditions in matrix and rhs
    if(j==0) {
      eliminate_DirichletBC_blockFE_blockA(bcond,&FE,&mesh,time_stepper.rhs_time,time_stepper.At,time_stepper.time);
    } else {
      eliminate_DirichletBC_RHS_blockFE_blockA(bcond,&FE,&mesh,time_stepper.rhs_time,time_stepper.At_noBC,time_stepper.time);
    }

    // Solve
    clock_t clk_solve_start = clock();

    // Solve the linear system
    // solve in block CSR format
    bdcsr_shift(time_stepper.At, -1);  // shift A

    if (linear_itparam.linear_precond_type >= 10 & linear_itparam.linear_precond_type < 13) {

      // get preconditioner diagonal blocks
      dcsr_alloc(time_stepper.At->blocks[0]->row, time_stepper.At->blocks[0]->col, time_stepper.At->blocks[0]->nnz, &A_diag[0]);
      dcsr_mxm(time_stepper.At->blocks[1],time_stepper.At->blocks[2],&BTB);
      dcsr_add(&BTB, 100.0, time_stepper.At->blocks[0], 1.0, &A_diag[0]);

      dcsr_alloc(FE_h.ndof, FE_h.ndof, FE_h.ndof, &A_diag[1]);
      for (i=0; i<=A_diag[1].row; i++)
        A_diag[1].IA[i] = i;
      for (i=0; i<A_diag[1].row; i++)
        A_diag[1].JA[i] = i;
      for (i=0; i<A_diag[1].row; i++)
        A_diag[1].val[i] = mesh.el_vol[i];

      // solve
      solver_flag = linear_solver_bdcsr_krylov_block_2(time_stepper.At, time_stepper.rhs_time, time_stepper.sol, &linear_itparam, &amgparam, A_diag);

      // free spaces
      dcsr_free(&A_diag[0]);
      dcsr_free(&A_diag[1]);

    }
    // solver diaognal blocks inexactly
    else if ( (linear_itparam.linear_precond_type >= 20 & linear_itparam.linear_precond_type < 23) || (linear_itparam.linear_precond_type >= 30 & linear_itparam.linear_precond_type < 33) ) {

      // element volume (only for RT0-P0)
      el_vol.row = time_stepper.At->blocks[2]->row;
      el_vol.val = mesh.el_vol;

      // solve
      solver_flag = linear_solver_bdcsr_krylov_mixed_darcy(time_stepper.At, time_stepper.rhs_time, time_stepper.sol, &linear_itparam, &amgparam, &el_vol);

    }
    // solver without preconditioner
    else {

      // solve
      solver_flag = linear_solver_bdcsr_krylov(time_stepper.At, time_stepper.rhs_time, time_stepper.sol, &linear_itparam);

    }

    bdcsr_shift(time_stepper.At, 1);   // shift A back

    // Error Check
    if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n", solver_flag);

    clock_t clk_solve_end = clock();
    printf("Elapsed CPU Time for Solve = %f seconds.\n\n",(REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);

    clock_t clk_timestep_end = clock();
    printf("Elapsed CPU Time for Time Step = %f seconds.\n\n",(REAL) (clk_timestep_end-clk_timestep_start)/CLOCKS_PER_SEC);

    // Output Solutions
    if (inparam.output_dir!=NULL) {
      // Solution at each timestep
      sprintf(solout,"output/solution_ts%03d.vtu",time_stepper.current_step);
      dump_blocksol_vtk(solout,varname,&mesh,&FE,time_stepper.sol->val);
    }
  } // End Timestepping Loop
  clock_t clk_timeloop_end = clock();
  printf("Elapsed CPU Time ALL Time Steps = %f seconds.\n\n",(REAL) (clk_timeloop_end-clk_timeloop_start)/CLOCKS_PER_SEC);

  /*******************************************************************************************/

  /******** Compute Errors or Plot *************************************************************/

  // Combine all timestep vtks in one file
  if(inparam.print_level > 3) {
    create_pvd("output/solution.pvd",time_stepper.tsteps+1,"solution_ts","timestep");
  }
  /*******************************************************************************************/
  


  /******** Free All the Arrays *************************************************************/

  // Time Stepper (and all the matrices inside)
  free_blktimestepper(&time_stepper);
  dvec_free(&b);
  dvec_free(&b_bdry);

  // Preconditioners
  if(A_diag) free(A_diag);

  // Extra Vectors and Arrays
  dvec_free(&sol);
  if(varname) free(varname);

  // FE Spaces
  free_fespace(&FE_h);
  free_fespace(&FE_q);
  free_blockfespace(&FE);

  // Quadrature
  if(cq) {
    free_qcoords(cq);
    free(cq);
    cq=NULL;
  }

  // Mesh
  free_mesh(&mesh);
  
  /*****************************************************************************************/

  clock_t clk_overall_end = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",(REAL) (clk_overall_end-clk_overall_start)/CLOCKS_PER_SEC);
  return 0;

}	/* End of Program */
/*******************************************************************************************/
