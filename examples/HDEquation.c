/*
 *  HDEquation.c
 *
 *  Created by James Adler and Xiaozhe Hu on 1/9/15.
 *  Copyright 2015_HAZMAT__. All rights reserved.
 *
 *  Discussion:
 *
 *    This program solves the following PDE using finite elements
 *
 *      D^*(a(x)D(u)) + c(x)u = 0
 *
 *    where D = grad, D^* = -div for P1,P2 elements,
 *          D = curl, D^* = curl for Nedelec elements
 *          D = div, D^* = -grad for Raviart-Thomas elements
 *
 *    in 2D or 3D
 *
 *   Along the boundary of the region, Dirichlet conditions are imposed:
 *
 *      u = 0 for P1, P2
 *    u*t = 0 for Nedelec
 *    u*n = 0 for Raviart-Thomas
 */

/*********** HAZMAT FUNCTIONS and INCLUDES **********************************************/
#include "hazmat.h"
/****************************************************************************************/

/******** Data Input ********************************************************************/
// Boundary Conditions
void bc(REAL *val,REAL* x,REAL time) {
  *val = 0.0;
}

// PDE Coefficients
void diffusion_coeff(REAL *val,REAL* x,REAL time) {
  *val = 1.0;
}

void reaction_coeff(REAL *val,REAL* x,REAL time) {
  *val = 1.0;
}

// True Solution (if you have one)
void truesol(REAL *val,REAL* x,REAL time) {
  // 2D - grad grad
  //*val = sin(M_PI*x[0])*sin(M_PI*x[1]);
  // 2D - curl curl
  val[0] = x[1]*(1-x[1]);
  val[1] = x[0]*(1-x[0]);
  // 3D
  //*val = sin(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
}

// Right-hand Side
void myrhs(REAL *val,REAL* x,REAL time) {
  REAL myc=-666.6;
  REAL mya=-666.6;
  //REAL myu=-666.6;
  REAL myu[2];
  reaction_coeff(&myc,x,time);
  diffusion_coeff(&mya,x,time);
  truesol(myu,x,time);
  // 2D - grad grad
  //*val = mya*2*M_PI*M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]) + myc*myu;
  // 2D - curl curl
  val[0] = mya*2.0 + myc*myu[0];
  val[2] = mya*2.0 + myc*myu[1];
  // 3D
  //*val = mya*3*M_PI*M_PI*sin(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]) + myc*myu;
}


/****** MAIN DRIVER *********************************************************************/
int main (int argc, char* argv[])
{
  printf("\n===========================================================================\n");
  printf("Beginning Program to solve H(D) problem: <D u, D v> + <u,v> = <f,v>.\n");
  printf("===========================================================================\n");
  /****** INITIALIZE PARAMETERS ********************************************************/
  // Timing Parameters
  clock_t clk_start,clk_end,clk1,clk2;
  clk_start = clock();
    
  //------------------------//
  // Step 0. Set parameters //
  //------------------------//
  input_param     inparam; // parameters from input files
  param_input("./input.dat", &inparam); // read in
    
  // Open gridfile for reading
  printf("\nCreating mesh and FEM spaces:\n");
  FILE* gfid = fopen(inparam.gridfile,"r");
  if( gfid == NULL ) {
    printf("\nError opening Grid File!!!\n");
    printf("File (%s) probably doesn't exist!\n\n",inparam.gridfile);
    return 0;
  }
    
  // Dimension is needed for all this to work
  INT dim = inparam.dim;
    
  // Create the mesh (now we assume triangles in 2D or tetrahedra in 3D)
  clk1 = clock();
  trimesh mesh;
  printf(" --> loading grid from file: %s\n",inparam.gridfile);
  initialize_mesh(&mesh);
  creategrid(gfid,dim,0,&mesh);
  fclose(gfid);
    
  // Get Quadrature Nodes for the Mesh
  INT nq1d = inparam.nquad; /* Quadrature points per dimension */
  qcoordinates *cq = get_quadrature(&mesh,nq1d);
    
  // Get info for and create FEM spaces
  // Order of Elements: 0 - P0; 1 - P1; 2 - P2; -1 - Nedelec; -2 - Raviart-Thomas
  INT order = inparam.FE_type;
  fespace FE;
  initialize_fespace(&FE);
  create_fespace(&FE,&mesh,order);
  char elmtype[8];
  if(order>=0) {
    sprintf(elmtype,"P%d",order);
    printf(" --> using P%d elements => D = grad\n",order);
  } else if(order==-1) {
    sprintf(elmtype,"Ned");
    printf(" --> using Nedelec elements => D = curl\n");
  } else if(order==-2) {
    sprintf(elmtype,"RT");
    printf(" --> using Raviart-Thomas elements => D = div\n");
  } else {
    printf("ERROR: Unknown Finite Element Type\n");
    exit(0);
  }
    
  if(inparam.print_level > 3) {
    dump_fespace(&FE);
  }
    
  clk2 = clock();
  printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
	 (REAL) (clk2 - clk1)/CLOCKS_PER_SEC);
  /**************************************************************************************/
    
  printf("***********************************************************************************\n");
  printf("Number of Elements = %d\tElement Type = %s\tOrder of Quadrature = %d\n",mesh.nelm,elmtype,2*nq1d-1);
  printf("\n\t--- Degrees of Freedom ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nv,mesh.nedge,mesh.nface);
  printf("\t--> DOF: %d\n",FE.ndof);
  printf("\n\t--- Boundaries ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nbv,mesh.nbedge,mesh.nbface);
  printf("\t--> Boundary DOF: %d\n",FE.nbdof);
  printf("***********************************************************************************\n\n");
    
  /*** Assemble the matrix and right hand side *******************************/
  printf("Assembling the matrix and right-hand side:\n");
  clk1 = clock();
    
  // Allocate the right-hand side and declare the csr matrix
  dvector b;
  dvector bnull;
  dCSRmat Diff;
  dCSRmat A;
  dCSRmat Mass;
    
  // Assemble the matrix without BC
  // Diffusion block
  assemble_global(&Diff,&b,assemble_DuDv_local,&FE,&mesh,cq,myrhs,diffusion_coeff,0.0);
  // Reaction block
  assemble_global(&Mass,&bnull,assemble_mass_local,&FE,&mesh,cq,myrhs,reaction_coeff,0.0);

  // Add the M + D
  dcsr_add_1(&Diff,1.0,&Mass,1.0,&A);
  if(bnull.val) free(bnull.val);
  dcsr_free(&Diff);
  dcsr_free(&Mass);
  
  // Eliminate Dirichlet BC
  eliminate_DirichletBC(bc,&FE,&mesh,&b,&A,0.0);

  FILE* matid = fopen("mat.dat","w");
  csr_print_matlab(matid,&A);
  fclose(matid);
  FILE* rhsid = fopen("rhs.dat","w");
  dvector_print(rhsid,&b);
  fclose(rhsid);
    
  clk2 = clock();
  printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL) (clk2-clk1)/CLOCKS_PER_SEC);
  /*******************************************************************************************/
    
    
  /**************** Solve ********************************************************************/
  printf("Solving the System:\n");
  clk1 = clock();
  // parameters
  INT solver_flag=-20;
  INT solver_type = inparam.linear_itsolver_type;
  REAL tol = inparam.linear_itsolver_tol;
  INT MaxIt = inparam.linear_itsolver_maxit;
  SHORT restart = 5;
  SHORT stop_type = STOP_REL_RES;
  SHORT print_level = PRINT_MORE;
    
  // Allocate the solution
  dvector u = dvec_create(b.row);
    
  // set initial guess to be all zero
  dvec_set(u.row, &u, 0.0);

  // shift A
  dcsr_shift(&A, -1);
    
  // solve the linear system
  if(solver_type==1) {
    printf(" --> using Conjugate Gradient Method:\n");
    solver_flag = dcsr_pcg(&A, &b, &u, NULL, tol, MaxIt, stop_type, print_level);
  } else if(solver_type==2) {
    printf(" --> using MINRES:\n");
  } else if(solver_type==3) {
    printf(" --> using GMRES:\n");
    solver_flag = dcsr_pvgmres(&A, &b, &u, NULL, tol, MaxIt, restart, stop_type, print_level);
  } else {
    printf("Unknown Solver Type\n");
    exit(0);
  }
  
  // Error Check
  if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n", solver_flag);
    
  // shift A back
  dcsr_shift(&A, 1);

  clk2 = clock();
  printf(" --> elapsed CPU time for solve = %f seconds.\n\n",(REAL) (clk2-clk1)/CLOCKS_PER_SEC);
  /*******************************************************************************************/    
       
  /**************** Compute Errors if you have true solution *********************************/
  printf("Computing True Solution and Errors:\n");
  clk2 = clock();
  REAL uerr = L2error(u.val,truesol,&FE,&mesh,cq,0.0);
  REAL graduerr = HDsemierror(u.val,truesol,&FE,&mesh,cq,0.0);
  REAL uH1err = HDerror(u.val,truesol,&FE,&mesh,cq,0.0);
  
  printf("************************************************************************************\n"); 
  printf("L2 Norm of u error      = %25.17g\n",uerr);
  printf("H1 Semi-Norm of u error = %25.17g\n",graduerr);
  printf("H1 Norm of u error      = %25.17g\n",uH1err);
  printf("************************************************************************************\n");
  printf(" --> elapsed CPU time for getting errors = %f seconds.\n\n",(REAL) (clk2-clk1)/CLOCKS_PER_SEC);
  /*******************************************************************************************/
    
  /**************** Print Results or Dump Results *******************************************/
  if (inparam.output_type==2) {
    FILE* uid = fopen("sol.dat","w");
    dvector_print(uid,&u);
    fclose(uid);
  }
  /*******************************************************************************************/
    
  /******** Free All the Arrays ***********************************************************/
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
  /****************************************************************************************/
    
  clk_end = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",(REAL) (clk_end-clk_start)/CLOCKS_PER_SEC);
  return 0;
    
}	/* End of Program */
/******************************************************************************************/


