/*! \file trabajo/periodic_test/periodic_test.c
 *
 *  Created by James Adler and Xiaozhe Hu 2019/06/14.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program solves the following PDE using finite elements
 *
 *        -div(a(x)*grad(u)) + c(x)*u = 0
 *
 *        with periodic boundary conditions.
 *
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
/*********************************************************************/

/******** Data Input *************************************************/
// PDE Coefficients
void diffusion_coeff(REAL *val,REAL* x,REAL time,void *param) {
  // a(x)
  *val = 1.0;
}
void reaction_coeff(REAL *val,REAL* x,REAL time,void *param) {
  // c(x)
  *val = 0.0;
}

// Exact Solution (if you have one)
void exactsol(REAL *val,REAL* x,REAL time,void *param) {
  // 2D
  *val = sin(M_PI*x[0])*sin(M_PI*x[1]);
}

// Derivative of Exact Solution (if you have one)
void D_exactsol(REAL *val,REAL* x,REAL time,void *param) {
  // 2D
  val[0] = M_PI*cos(M_PI*x[0])*sin(M_PI*x[1]);
  val[1] = M_PI*sin(M_PI*x[0])*cos(M_PI*x[1]);
}

// Right-hand Side
void rhs(REAL *val,REAL* x,REAL time,void *param) {
  // 2D
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu=-666.6;
  reaction_coeff(&myc,x,time,param);
  diffusion_coeff(&mya,x,time,param);
  exactsol(&myu,x,time,param);
  *val = (mya*2*M_PI*M_PI + myc)*myu;
}

// Boundary Conditions
void bc(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu;
  exactsol(&myu,x,time,param);
  *val= myu;
}
/*********************************************************************/

/****** MAIN DRIVER **************************************************/
int main (int argc, char* argv[])
{

  printf("\n===========================================================================\n");
  printf("Beginning Program to solve periodic test problem: <grad u, grad v> + <u,v> = <f,v>.\n");
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

  // Open gridfile for reading
  printf("\nCreating mesh and FEM spaces:\n");
  FILE* gfid = HAZ_fopen(inparam.gridfile,"r");

  // Create the mesh
  // File types possible are 0 - HAZ format; 1 - VTK format
  INT mesh_type = 0;
  clock_t clk_mesh_start = clock(); // Time mesh generation FE setup
  mesh_struct mesh;
  printf(" --> loading grid from file: %s\n",inparam.gridfile);
  creategrid_fread(gfid,mesh_type,&mesh);
  fclose(gfid);

  // Dimension is needed for all this to work
  INT dim = mesh.dim;

  // Get Quadrature Nodes for the Mesh
  INT nq1d = inparam.nquad; // Quadrature points per dimension
  qcoordinates *cq = get_quadrature(&mesh,nq1d);

  // Get info for and create FEM spaces
  // Order of Elements:
  //    0 - P0; 1 - P1; 2 - P2; -1 - Nedelec; -2 - Raviart-Thomas
  INT order = inparam.FE_type;
  fespace FE;
  create_fespace(&FE,&mesh,order);

  // Set Dirichlet Boundaries
  // Assume physical boundaries (flag of 1 in mesh file) are Dirichlet
  set_dirichlet_bdry(&FE,&mesh,1,1);

  // Set periodic Boundaries
  // TODO
  set_periodic_bdry(&FE,&mesh,0.0,1.0,0.0,1.0,0.0,1.0);
  for(INT i=0;i<FE.ndof;i++)
    printf("periodic[%d]=%d\n",i,FE.periodic[i]);
  exit(0);

  // Strings for printing
  char elmtype[8];
  if(order>=0 && order<10) {
    sprintf(elmtype,"P%d",order);
    printf(" --> using P%d elements => D = grad\n",order);
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
    char* namevtk = "output/mesh.vtu";
    dump_mesh_vtk(namevtk,&mesh);
  }

  clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
  printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
         (REAL) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
  /*******************************************************************/

  printf("***********************************************************************************\n");
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
  dCSRmat Diff;
  dCSRmat A;
  dCSRmat Mass;

  // Assemble the matrix without BC
  // Different cases for dimension and FE of test problem

  // Diffusion block
  assemble_global(&Diff,&b,assemble_DuDv_local,&FE,&mesh,cq,rhs,diffusion_coeff,0.0);

  // Reaction block
  assemble_global(&Mass,NULL,assemble_mass_local,&FE,&mesh,cq,NULL,reaction_coeff,0.0);

  // Add M + D
  dcsr_add(&Diff,1.0,&Mass,1.0,&A);
  dcsr_free(&Diff);
  dcsr_free(&Mass);

  // Eliminate Dirichlet BC
  eliminate_DirichletBC(bc,&FE,&mesh,&b,&A,0.0);

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
  //param_amg_print(&amgparam);
  //=================================================================//

  // Solve the linear system
  if(linear_itparam.linear_itsolver_type == 0) { // Direct Solver
#if WITH_SUITESPARSE
    printf(" --> using UMFPACK's Direct Solver:\n");
    solver_flag = directsolve_UMF(&A,&b,&u,linear_itparam.linear_print_level);
#else
    error_extlib(255,__FUNCTION__,"SuiteSparse");
    return 0;
#endif
  } else { // Iterative Solver

    // Use AMG as iterative solver
    if (linear_itparam.linear_itsolver_type == SOLVER_AMG){
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
          get_curl_NedtoRT(&Curl,&mesh);
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
  if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n", solver_flag);

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
  uerr = L2error(u.val,exactsol,&FE,&mesh,cq,0.0);
  graduerr = HDsemierror(u.val,D_exactsol,&FE,&mesh,cq,0.0);
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
  if (inparam.output_dir != NULL) {
    char solout[128];
    strncpy(solout,inparam.output_dir,128);
    strcat(solout,"sol.vtu");

    dump_sol_vtk(solout,"u",&mesh,&FE,u.val);

    dvector exact_sol = dvec_create(FE.ndof);
    FE_Evaluate(exact_sol.val,exactsol,&FE,&mesh,0.0);

    char exactout[128];
    strncpy(exactout,inparam.output_dir,128);
    strcat(exactout,"exact.vtu");
    dump_sol_vtk(exactout,"ut",&mesh,&FE,exact_sol.val);

    dvec_free(&exact_sol);
  }
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
