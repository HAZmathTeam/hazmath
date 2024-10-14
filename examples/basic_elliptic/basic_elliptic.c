/*! \file examples/basic_elliptic/basic_elliptic.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 2015/01/09.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program solves the following PDE using finite elements
 *
 *        D^*(a(x)D(u)) + c(x)u = 0
 *
 *        where D = grad, D^* = -div for P1,P2 elements,
 *              D = curl, D^* = curl for Nedelec elements
 *              D = div,  D^* = -grad for Raviart-Thomas elements
 *
 *        in 1D, 2D, or 3D
 *
 *        Along the boundary of the region, Dirichlet conditions are
 *        imposed:
 *
 *          u = 0 for P1, P2
 *          u*t = 0 for Nedelec
 *          u*n = 0 for Raviart-Thomas
 *
 * \note This example highlights some of the basic features of HAZmath,
 * including how to create a basic uniform mesh, set up finite-element spaces on
 * a given mesh, create linear systems, and solve those systems using a variety
 * of Krylov and/or multigrid solvers (a direct solver can also be implemented).
 * It also illustrates how to set the problem data, such as boundary conditions
 * and right-hand sides, and how to output the solution in VTK formats.
 *
 * \note This is intended to give you different examples in different dimensions
 * and with different types of coefficients.  There are lots of if statements
 * that would be unnecessary in your own program.
 *
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
#include "basic_elliptic_data.h"
/*********************************************************************/

/****** MAIN DRIVER **************************************************/
int main (int argc, char* argv[])
{

  printf("\n===========================================================================\n");
  printf("Beginning Program to solve H(D) problem: <D u, D v> + <u,v> = <f,v>.\n");
  printf("===========================================================================\n");

  /****** INITIALIZE PARAMETERS **************************************/
  // Flag for errors
  SHORT status;

  // Overall CPU Timing
  clock_t clk_overall_start = clock();

  // Set Parameters by Reading in Input File
  input_param inparam;
  param_input_init(&inparam);
  param_input("./input.dat", &inparam);

  // Create the mesh
  // We have two options:
  // 1) Use Hazmath's built-in functions to create a uniform mesh from a
  //    simplicial complex
  // 2) Read in a user-define mesh defined in .haz format  See the example 2D
  //    mesh that is provided for formatting guidelines
  INT read_mesh_from_file=inparam.read_mesh_from_file;

  printf("\nCreating mesh and FEM spaces:\n");

  // Time the mesh generation and FE setup
  clock_t clk_mesh_start = clock();

  // Initialize some mesh parameters
  mesh_struct mesh;
  // Note the following are overwritten by a mesh file if provided and read in
  INT dim = inparam.spatial_dim;                 // dimension of computational domain
  INT mesh_ref_levels=inparam.refinement_levels; // refinement levels
  INT mesh_ref_type=inparam.refinement_type;     // refinement type (>10 uniform or <10 other)
  INT set_bndry_codes=inparam.boundary_codes;    // set flags for the boundary DoF ([1:16] are Dirichlet, [17:32] are Neumann, etc)

  if(read_mesh_from_file) {
    // Reading from file
    FILE* gfid = HAZ_fopen(inparam.gridfile,"r");
    // File types possible are 0 - HAZ format; 1 - VTK format (not implemented)
    INT mesh_type = 0;
    printf(" --> loading grid from file: %s\n",inparam.gridfile);
    creategrid_fread(gfid,mesh_type,&mesh);
    fclose(gfid);
  } else {
    // Use HAZMATH built in functions for a uniform mesh in 2D or 3D
    mesh=make_uniform_mesh(dim,mesh_ref_levels,mesh_ref_type,set_bndry_codes);
  }

  // Get Quadrature Nodes for the Mesh
  INT nq1d = inparam.nquad; // Quadrature points per dimension
  qcoordinates *cq = get_quadrature(&mesh,nq1d);

  // Get info for and create FEM spaces
  // Order of Elements:
  //    0 - P0; 1 - P1; 2 - P2; 20 - Nedelec; 30 - Raviart-Thomas
  INT order = inparam.FE_type;
  fespace FE;
  create_fespace(&FE,&mesh,order);

  // Set Dirichlet Boundaries
  // Assume physical boundaries are Dirichlet
  // The mesh is set up so that flag values 1-16 are Dirichlet and 17-32 are Neumann
  set_dirichlet_bdry(&FE,&mesh,1,16);

  // Strings for printing
  char elmtype[8];
  if(order>=0 && order<10) {
    sprintf(elmtype,"P%lld",(long long )order);
    printf(" --> using P%lld elements => D = grad\n",(long long )order);
  } else if(order==20) {
    sprintf(elmtype,"Ned");
    printf(" --> using Nedelec elements => D = curl\n");
  } else if(order==30) {
    sprintf(elmtype,"RT");
    printf(" --> using Raviart-Thomas elements => D = div\n");
  } else {
    printf("ERROR: Unknown Finite Element Type\n");
    exit(0);
  }

  // Dump some of the data
  if(inparam.print_level > 3) {
    // FE space
    char varu[10];
    char dir[20];
    sprintf(dir,"output");
    sprintf(varu,"u");
    dump_fespace(&FE,varu,dir);

    // Mesh
    char* namevtk = strdup("output/mesh.vtu");
    dump_mesh_vtk(namevtk,&mesh);
    char* namemesh = strdup("output/mesh.haz");
    dump_mesh_haz(namemesh,&mesh);
    free(namemesh);
    free(namevtk);
  }

  clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
  printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
         (REAL) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
  /*******************************************************************/

  printf("***********************************************************************************\n");
  printf("\t--- %lld-dimensional grid ---\n",(long long )dim);
  printf("Number of Elements = %lld\tElement Type = %s\tOrder of Quadrature = %lld\n",(long long )mesh.nelm,elmtype,2*(long long )nq1d-1);
  printf("\n\t--- Degrees of Freedom ---\n");
  printf("Vertices: %-7lld\tEdges: %-7lld\tFaces: %-7lld",(long long )mesh.nv,(long long )mesh.nedge,(long long )mesh.nface);
  printf("\t--> DOF: %lld\n",(long long )FE.ndof);
  printf("\n\t--- Boundaries ---\n");
  printf("Vertices: %-7lld\tEdges: %-7lld\tFaces: %-7lld",(long long )mesh.nbv,(long long )mesh.nbedge,(long long )mesh.nbface);
  printf("\t--> Boundary DOF: %lld\n",(long long )FE.nbdof);
  printf("***********************************************************************************\n\n");

  /*** Assemble the matrix and right hand side ***********************/
  printf("Assembling the matrix and right-hand side:\n");
  clock_t clk_assembly_start = clock();

  // Allocate the right-hand side and declare the csr matrix
  dvector b;
  dCSRmat Diff;
  dCSRmat A;
  dCSRmat Mass;

  // Assemble the matrix without BC
  // Different cases for dimension and FE of test problem

  // Diffusion block
  if(dim==1) {
    if(FE.FEtype>0 && FE.FEtype<10) { // PX
      assemble_global(&Diff,&b,assemble_DuDv_local,&FE,&mesh,cq,
                      rhs_1D_PX,diffusion_coeff,0.0);
    } else {
      status = ERROR_FE_TYPE;
      check_error(status, __FUNCTION__);
    }
  } else if(dim==2) {
    if(FE.FEtype>0 && FE.FEtype<10) { // PX
      assemble_global(&Diff,&b,assemble_DuDv_local,&FE,&mesh,cq,
                      rhs_2D_PX,diffusion_coeff,0.0);
    } else if(FE.FEtype==20) { // Nedelec
      assemble_global(&Diff,&b,assemble_DuDv_local,&FE,&mesh,cq,
                      rhs_2D_Ned,diffusion_coeff,0.0);
    } else if(FE.FEtype==30) { // RT
      assemble_global(&Diff,&b,assemble_DuDv_local,&FE,&mesh,cq,
                      rhs_2D_RT,diffusion_coeff,0.0);
    } else {
      status = ERROR_FE_TYPE;
      check_error(status, __FUNCTION__);
    }
  } else if(dim==3) {
    if(FE.FEtype>0 && FE.FEtype<10) { // PX
      assemble_global(&Diff,&b,assemble_DuDv_local,&FE,&mesh,cq,
                      rhs_3D_PX,diffusion_coeff,0.0);
    } else if(FE.FEtype==20) { // Nedelec
      assemble_global(&Diff,&b,assemble_DuDv_local,&FE,&mesh,cq,
                      rhs_3D_Ned,diffusion_coeff,0.0);
    } else if(FE.FEtype==30) { // RT
      assemble_global(&Diff,&b,assemble_DuDv_local,&FE,&mesh,cq,
                      rhs_3D_RT,diffusion_coeff,0.0);
    } else {
      status = ERROR_FE_TYPE;
      check_error(status, __FUNCTION__);
    }
  } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
  }

  // Reaction block
  assemble_global(&Mass,NULL,assemble_mass_local,&FE,&mesh,cq,NULL,
                  reaction_coeff,0.0);

  // Add M + D
  dcsr_add(&Diff,1.0,&Mass,1.0,&A);
  dcsr_free(&Diff);
  dcsr_free(&Mass);

  // Eliminate Dirichlet BC
  // Different cases for dimension and FE of test problem
  if(dim==1) {
    if(FE.FEtype>0 && FE.FEtype<10) { // PX
      eliminate_DirichletBC(bc_1D_PX,&FE,&mesh,&b,&A,0.0);
    } else {
      status = ERROR_FE_TYPE;
      check_error(status, __FUNCTION__);
    }
  } else if(dim==2) {
    if(FE.FEtype>0 && FE.FEtype<10) { // PX
      eliminate_DirichletBC(bc_2D_PX,&FE,&mesh,&b,&A,0.0);
    } else if(FE.FEtype==20) { // Nedelec
      eliminate_DirichletBC(bc_2D_Ned,&FE,&mesh,&b,&A,0.0);
    } else if(FE.FEtype==30) { // RT
      eliminate_DirichletBC(bc_2D_RT,&FE,&mesh,&b,&A,0.0);
    } else {
      status = ERROR_FE_TYPE;
      check_error(status, __FUNCTION__);
    }
  } else if(dim==3) {
    if(FE.FEtype>0 && FE.FEtype<10) { // PX
      eliminate_DirichletBC(bc_3D_PX,&FE,&mesh,&b,&A,0.0);
    } else if(FE.FEtype==20) { // Nedelec
      eliminate_DirichletBC(bc_3D_Ned,&FE,&mesh,&b,&A,0.0);
    } else if(FE.FEtype==30) { // RT
      eliminate_DirichletBC(bc_3D_RT,&FE,&mesh,&b,&A,0.0);
    } else {
      status = ERROR_FE_TYPE;
      check_error(status, __FUNCTION__);
    }
  } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
  }

  // Dump matrices for testing
  if(inparam.print_level > 3) {
    FILE* matid = HAZ_fopen("output/mat.dat","w");
    csr_print_matlab(matid,&A);
    fclose(matid);
    FILE* rhsid = HAZ_fopen("output/rhs.dat","w");
    dvector_print(rhsid,&b);
    fclose(rhsid);
  }

  clock_t clk_assembly_end = clock();
  printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL)
           (clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);
  /*******************************************************************/

  /**************** Solve ********************************************/
  printf("Solving the System:\n");
  clock_t clk_solve_start = clock();

  // Create Solution Vector
  dvector u = dvec_create(FE.ndof);

  // Set initial guess to be all zero
  dvec_set(u.row, &u, 0.0);

  // Set Solver Parameters
  INT solver_flag=-20;

  // Set parameters for linear iterative methods
  linear_itsolver_param linear_itparam;
  param_linear_solver_set(&linear_itparam, &inparam);

  // Set parameters for algebriac multigrid methods
  AMG_param amgparam;
  param_amg_init(&amgparam);
  param_amg_set(&amgparam, &inparam);
  param_amg_print(&amgparam);

  // choose preconditioner for iterative methods
  if (linear_itparam.linear_itsolver_type != 0)  // if using direct solver, do nothing
  {
    if (FE.FEtype == 20) // Nedelec element
    {
      linear_itparam.linear_precond_type = PREC_HX_CURL_M; // use HX preconditioner for curl
    }
    else if (FE.FEtype == 30) // RT element
    {
      if (dim == 3)
      {
        linear_itparam.linear_precond_type = PREC_HX_DIV_M; // use HX preconditioner for div
      }
      else
      {
        linear_itparam.linear_precond_type = PREC_AMG;
      }
    }
    else // Lagrange element
    {
      linear_itparam.linear_precond_type = PREC_AMG; // use AMG
    }
  }

  //=================================================================//

  // Solve the linear system
  if(linear_itparam.linear_itsolver_type == 0) { // Direct Solver
    //#if WITH_SUITESPARSE
    printf(" --> using Direct Solver (UMFPACK|HAZMATH):\n");
    solver_flag = directsolve_HAZ(&A,&b,&u,linear_itparam.linear_print_level);
    //#else
    //    error_extlib(255,__FUNCTION__,"SuiteSparse");
    //    return 0;
    //#endif
  } else { // Iterative Solver

    // Use AMG as iterative solver
    if (linear_itparam.linear_itsolver_type == SOLVER_AMG){
      printf("Start calling linear solver for A:\n");
      solver_flag = linear_solver_amg(&A, &b, &u, &amgparam);
      printf("Start calling linear solver for A again:\n");
      solver_flag = linear_solver_amg(&A, &b, &u, &amgparam);
      printf("Start calling linear solver for A (3):\n");
      solver_flag = linear_solver_amg(&A, &b, &u, &amgparam);
      printf("Start calling linear solver for A (4):\n");
      solver_flag = linear_solver_amg(&A, &b, &u, &amgparam);
      printf("Start calling linear solver for A (5):\n");
      solver_flag = linear_solver_amg(&A, &b, &u, &amgparam);
      printf("Start calling linear solver for A (6):\n");
      solver_flag = linear_solver_amg(&A, &b, &u, &amgparam);
    } else { // Use Krylov Iterative Solver
      // Determine Preconditioner
      // Diagonal preconditioner
      if (linear_itparam.linear_precond_type == PREC_DIAG) {
          solver_flag = linear_solver_dcsr_krylov_diag(&A, &b, &u, &linear_itparam);
      }
      // AMG preconditioner
      else if (linear_itparam.linear_precond_type == PREC_AMG){
          solver_flag = linear_solver_dcsr_krylov_amg(&A, &b, &u, &linear_itparam, &amgparam);
      }
      // HX preconditioner (for H(curl) problem only)
      else if (linear_itparam.linear_precond_type == PREC_HX_CURL_A || linear_itparam.linear_precond_type == PREC_HX_CURL_M){

          // declare data
          dCSRmat P_curl;
          dCSRmat Grad;

          // get P_curl and Grad
          get_Pigrad_H1toNed(&P_curl,&mesh);
          get_grad_H1toNed(&Grad,&mesh);

          solver_flag = linear_solver_dcsr_krylov_hx_curl(&A, &b, &u, &linear_itparam, &amgparam, &P_curl, &Grad);

          // clean
          dcsr_free(&P_curl);
          dcsr_free(&Grad);

      }
      // HX preconditioner (for H(div) problem only)
      else if (linear_itparam.linear_precond_type == PREC_HX_DIV_A || linear_itparam.linear_precond_type == PREC_HX_DIV_M){

          // declare data
          dCSRmat P_curl;
          dCSRmat Curl;
          dCSRmat P_div;

          // get P_curl and Grad
          get_Pigrad_H1toNed(&P_curl,&mesh);
          get_curl_NedtoRT(&Curl,&mesh); // needs to have a 2D version for this--XH
          get_Pigrad_H1toRT(&P_div,&P_curl,&Curl,&mesh);

          solver_flag = linear_solver_dcsr_krylov_hx_div(&A, &b, &u, &linear_itparam, &amgparam,&P_curl,&P_div,&Curl);

          // clean
          dcsr_free(&P_curl);
          dcsr_free(&Curl);
          dcsr_free(&P_div);

      }
      // No preconditoner
      else{
          solver_flag = linear_solver_dcsr_krylov(&A, &b, &u, &linear_itparam);
      }

    }
  }

  // Error Check
  if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %lld!\n", (long long )solver_flag);

  clock_t clk_solve_end = clock();
  printf("Elapsed CPU Time for Solve = %f seconds.\n\n",
         (REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);
  /*******************************************************************/

  /**************** Compute Errors if you have exact solution *********/
  // Again this depends on dimension and FE type
  printf("Computing Exact Solution and Errors:\n");
  clock_t clk_error_start = clock();
  REAL uerr=0.0;
  REAL graduerr=0.0;
  if(dim==1) {
    if(FE.FEtype>0 && FE.FEtype<10) { // PX
      uerr = L2error(u.val,exactsol_1D_PX,&FE,&mesh,cq,0.0);
      graduerr = HDsemierror(u.val,D_exactsol_1D_PX,&FE,&mesh,cq,0.0);
    } else {
      status = ERROR_FE_TYPE;
      check_error(status, __FUNCTION__);
    }
  } else if(dim==2) {
    if(FE.FEtype>0 && FE.FEtype<10) { // PX
      uerr = L2error(u.val,exactsol_2D_PX,&FE,&mesh,cq,0.0);
      graduerr = HDsemierror(u.val,D_exactsol_2D_PX,&FE,&mesh,cq,0.0);
    } else if(FE.FEtype==20) { // Nedelec
      uerr = L2error(u.val,exactsol_2D_Ned,&FE,&mesh,cq,0.0);
      graduerr = HDsemierror(u.val,D_exactsol_2D_Ned,&FE,&mesh,cq,0.0);
    } else if(FE.FEtype==30) { // RT
      uerr = L2error(u.val,exactsol_2D_RT,&FE,&mesh,cq,0.0);
      graduerr = HDsemierror(u.val,D_exactsol_2D_RT,&FE,&mesh,cq,0.0);
    } else {
      status = ERROR_FE_TYPE;
      check_error(status, __FUNCTION__);
    }
  } else if(dim==3) {
    if(FE.FEtype>0 && FE.FEtype<10) { // PX
      uerr = L2error(u.val,exactsol_3D_PX,&FE,&mesh,cq,0.0);
      graduerr = HDsemierror(u.val,D_exactsol_3D_PX,&FE,&mesh,cq,0.0);
    } else if(FE.FEtype==20) { // Nedelec
      uerr = L2error(u.val,exactsol_3D_Ned,&FE,&mesh,cq,0.0);
      graduerr = HDsemierror(u.val,D_exactsol_3D_Ned,&FE,&mesh,cq,0.0);
    } else if(FE.FEtype==30) { // RT
      uerr = L2error(u.val,exactsol_3D_RT,&FE,&mesh,cq,0.0);
      graduerr = HDsemierror(u.val,D_exactsol_3D_RT,&FE,&mesh,cq,0.0);
    } else {
      status = ERROR_FE_TYPE;
      check_error(status, __FUNCTION__);
    }
  } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
  }
  REAL uH1err = sqrt(uerr*uerr + graduerr*graduerr);

  printf("************************************************************************************\n");
  printf("L2 Norm of u error      = %26.13e\n",uerr);
  printf("H1 Semi-Norm of u error = %26.13e\n",graduerr);
  printf("H1 Norm of u error      = %26.13e\n",uH1err);
  printf("************************************************************************************\n");
  clock_t clk_error_end = clock();
  printf("Elapsed CPU time for getting errors = %lf seconds.\n\n",(REAL)
         (clk_error_end-clk_error_start)/CLOCKS_PER_SEC);
  /*******************************************************************/

  /**************** Print Results or Dump Results ********************/
  char solout[128];
  strncpy(solout,inparam.output_dir,128);
  strcat(solout,"sol.vtu");

  dump_sol_vtk(solout,"u",&mesh,&FE,u.val);

  dvector exact_sol = dvec_create(FE.ndof);
  if(dim==1) {
    if(FE.FEtype>=0 && FE.FEtype<10) { // PX
      FE_Evaluate(exact_sol.val,exactsol_1D_PX,&FE,&mesh,0.0);
    } else {
      status = ERROR_FE_TYPE;
      check_error(status, __FUNCTION__);
    }
  } else if(dim==2) {
      if(FE.FEtype>=0 && FE.FEtype<10) { // PX
        FE_Evaluate(exact_sol.val,exactsol_2D_PX,&FE,&mesh,0.0);
      } else if(FE.FEtype==20) { // Nedelec
        FE_Evaluate(exact_sol.val,exactsol_2D_Ned,&FE,&mesh,0.0);
      } else if(FE.FEtype==30) { // RT
        FE_Evaluate(exact_sol.val,exactsol_2D_RT,&FE,&mesh,0.0);
      } else {
        status = ERROR_FE_TYPE;
        check_error(status, __FUNCTION__);
      }
  } else if(dim==3) {
    if(FE.FEtype>=0 && FE.FEtype<10) { // PX
      FE_Evaluate(exact_sol.val,exactsol_3D_PX,&FE,&mesh,0.0);
    } else if(FE.FEtype==20) { // Nedelec
      FE_Evaluate(exact_sol.val,exactsol_3D_Ned,&FE,&mesh,0.0);
    } else if(FE.FEtype==30) { // RT
      FE_Evaluate(exact_sol.val,exactsol_3D_RT,&FE,&mesh,0.0);
    } else {
      status = ERROR_FE_TYPE;
      check_error(status, __FUNCTION__);
    }
  } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
  }

  char exactout[128];
  strncpy(exactout,inparam.output_dir,128);
  strcat(exactout,"exact.vtu");
  dump_sol_vtk(exactout,"ut",&mesh,&FE,exact_sol.val);

  dvec_free(&exact_sol);

  /*******************************************************************/

  /******** Free All the Arrays **************************************/
  dcsr_free(&A);
  if(b.val) free(b.val);
  if(u.val) free(u.val);
  free_fespace(&FE);
  if(cq) {
    free_qcoords(cq);
    free(cq);
    cq = NULL;
  }
  free_mesh(&mesh);
  /*******************************************************************/

  clock_t clk_overall_end = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",
         (REAL) (clk_overall_end-clk_overall_start)/CLOCKS_PER_SEC);
  return 0;

}	/* End of Program */
/*******************************************************************/
