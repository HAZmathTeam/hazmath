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
// PDE Coefficients
void diffusion_coeff(REAL *val,REAL* x,REAL time) {
  *val = 1.0;
}

void reaction_coeff(REAL *val,REAL* x,REAL time) {
  *val = 1.0;
}

// True Solution (if you have one)
// Pick one of these and rename it truesol
void truesol_2D_PX(REAL *val,REAL* x,REAL time) {
  //void truesol(REAL *val,REAL* x,REAL time) {
  // 2D - grad grad
  *val = sin(M_PI*x[0])*sin(M_PI*x[1]);
}
void truesol_3D_PX(REAL *val,REAL* x,REAL time) {
  //void truesol(REAL *val,REAL* x,REAL time) {
  // 3D - grad grad
  *val = sin(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
}
void truesol_2D_Ned(REAL *val,REAL* x,REAL time) {
  //void truesol(REAL *val,REAL* x,REAL time) {
  // 2D - curl curl
  val[0] = cos(M_PI*x[0])*sin(M_PI*x[1]);
  val[1] = -sin(M_PI*x[0])*cos(M_PI*x[1]);
}
//void truesol_3D_Ned(REAL *val,REAL* x,REAL time) {
void truesol(REAL *val,REAL* x,REAL time) {
  // 3D - curl curl
  val[0] = cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
  val[1] = sin(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]);
  val[2] = -sin(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*x[2]);
}
void truesol_2D_RT(REAL *val,REAL* x,REAL time) {
  //void truesol(REAL *val,REAL* x,REAL time) {
  // 2D - grad div
  val[0] = sin(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = cos(M_PI*x[0])*sin(M_PI*x[1]);
}
void truesol_3D_RT(REAL *val,REAL* x,REAL time) {
  //void truesol(REAL *val,REAL* x,REAL time) {
  // 3D - grad div
  val[0] = sin(M_PI*x[0])*cos(M_PI*x[1])*cos(M_PI*x[2]);
  val[1] = cos(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*x[2]);
  val[2] = cos(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]);
}

// Right-hand Side
// Pick one of these and rename it myrhs
void rhs_2D_PX(REAL *val,REAL* x,REAL time) {
  //void myrhs(REAL *val,REAL* x,REAL time) {
  // 2D - grad grad
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu=-666.6;
  reaction_coeff(&myc,x,time);
  diffusion_coeff(&mya,x,time);
  truesol(&myu,x,time);
  *val = (mya*2*M_PI*M_PI + myc)*myu;
}
void rhs_3D_PX(REAL *val,REAL* x,REAL time) {
  //void myrhs(REAL *val,REAL* x,REAL time) {
  // 3D - grad grad
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu=-666.6;
  reaction_coeff(&myc,x,time);
  diffusion_coeff(&mya,x,time);
  truesol(&myu,x,time);
  *val = (mya*3*M_PI*M_PI + myc)*myu;
}
void rhs_2D_Ned(REAL *val,REAL* x,REAL time) {
  //void myrhs(REAL *val,REAL* x,REAL time) {
  // 2D - curl curl
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu[2];
  reaction_coeff(&myc,x,time);
  diffusion_coeff(&mya,x,time);
  truesol(myu,x,time);
  val[0] = (mya*2.0*M_PI*M_PI + myc)*myu[0];
  val[1] = (mya*2.0*M_PI*M_PI + myc)*myu[1];
}
//void rhs_3D_Ned(REAL *val,REAL* x,REAL time) {
void myrhs(REAL *val,REAL* x,REAL time) {
  // 3D - curl curl
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu[3];
  reaction_coeff(&myc,x,time);
  diffusion_coeff(&mya,x,time);
  truesol(myu,x,time);
  val[0] = (2*mya*M_PI*M_PI + myc)*myu[0];
  val[1] = (2*mya*M_PI*M_PI + myc)*myu[1];
  val[2] = (4*mya*M_PI*M_PI + myc)*myu[2];
}
void rhs_2D_RT(REAL *val,REAL* x,REAL time) {
  //void myrhs(REAL *val,REAL* x,REAL time) {
  // 2D - grad div
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu[2];
  reaction_coeff(&myc,x,time);
  diffusion_coeff(&mya,x,time);
  truesol(myu,x,time);
  val[0] = (mya*2.0*M_PI*M_PI + myc)*myu[0];
  val[1] = (mya*2.0*M_PI*M_PI + myc)*myu[1];
}
void rhs_3D_RT(REAL *val,REAL* x,REAL time) {
  //void myrhs(REAL *val,REAL* x,REAL time) {
  // 3D - grad div
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu[3];
  reaction_coeff(&myc,x,time);
  diffusion_coeff(&mya,x,time);
  truesol(myu,x,time);
  val[0] = (mya*3.0*M_PI*M_PI + myc)*myu[0];
  val[1] = (mya*3.0*M_PI*M_PI + myc)*myu[1];
  val[2] = (mya*3.0*M_PI*M_PI + myc)*myu[2];
}

// Boundary Conditions
// Switch one to bc
void bc_PX(REAL *val,REAL* x,REAL time) {
  //void bc(REAL *val,REAL* x,REAL time) {
  REAL myu;
  truesol(&myu,x,time);
  *val= myu;
}
void bc_2Dvec(REAL *val,REAL* x,REAL time) {
  //void bc(REAL *val,REAL* x,REAL time) {
  REAL myu[2];
  truesol(myu,x,time);
  val[0] = myu[0];
  val[1] = myu[1];
}
//void bc_3Dvec(REAL *val,REAL* x,REAL time) {
  void bc(REAL *val,REAL* x,REAL time) {
  REAL myu[3];
  truesol(myu,x,time);
  val[0] = myu[0];
  val[1] = myu[1];
  val[2] = myu[2];
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

  /* INT i,a1,a2,j; */
  /* for(i=0;i<mesh.nface;i++) { */
  /*   a1 = mesh.f_v->IA[i]-1; */
  /*   a2 = mesh.f_v->IA[i+1]-1; */
  /*   for(j=a1;j<a2;j++) { */
  /*     printf("(%f,%f) ",mesh.cv->x[mesh.f_v->JA[j]-1],mesh.cv->y[mesh.f_v->JA[j]-1]); */
  /*   } */
  /*   printf("A = %f, n = (%f,%f), m = (%f,%f)\n",mesh.f_area[i],mesh.f_norm[i*dim],mesh.f_norm[i*dim+1],mesh.f_mid[i*dim],mesh.f_mid[i*dim+1]); */
  /* } */
    
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
    
  // solve the linear system
  if(solver_type==0) {
    printf(" --> using UMFPACK's Direct Solver:\n");
    solver_flag = directsolve_UMF_symmetric(&A,&b,u.val,print_level);
  } else if(solver_type==1) {
    printf(" --> using Conjugate Gradient Method:\n");
    //dcsr_shift(&A, -1);  // shift A
    solver_flag = dcsr_pcg(&A, &b, &u, NULL, tol, MaxIt, stop_type, print_level);
    //dcsr_shift(&A, 1);   // shift A back
  } else if(solver_type==2) {
    printf(" --> using MINRES:\n");
    //dcsr_shift(&A, -1);  // shift A
    printf(" NOTHING IMPLEMENTED FOR MINRES\n");
    //dcsr_shift(&A, 1);   // shift A back
  } else if(solver_type==3) {
    printf(" --> using GMRES:\n");
    //dcsr_shift(&A, -1);  // shift A
    solver_flag = dcsr_pvgmres(&A, &b, &u, NULL, tol, MaxIt, restart, stop_type, print_level);
    //dcsr_shift(&A, 1);   // shift A back
  } else {
    printf("Unknown Solver Type\n");
    exit(0);
  }
  
  // Error Check
  if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n", solver_flag);  

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


