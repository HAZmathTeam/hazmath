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

/*********** HAZMAT FUNCTIONS and INCLUDES ****************************************/
#include "hazmat.h"
/**********************************************************************************/

/******** Data Input **************************************************************/
// PDE Coefficients
void diffusion_coeff(REAL *val,REAL* x,REAL time) {
  *val = 1.0;
}

void reaction_coeff(REAL *val,REAL* x,REAL time) {
  *val = 1.0;
}

// True Solution (if you have one)
void truesol_2D_PX(REAL *val,REAL* x,REAL time) {
  // 2D - grad grad
  *val = sin(M_PI*x[0])*sin(M_PI*x[1]);
}
void truesol_3D_PX(REAL *val,REAL* x,REAL time) {
  // 3D - grad grad
  *val = sin(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
}
void truesol_2D_Ned(REAL *val,REAL* x,REAL time) {
  // 2D - curl curl
  val[0] = cos(M_PI*x[0])*sin(M_PI*x[1]);
  val[1] = -sin(M_PI*x[0])*cos(M_PI*x[1]);
}
void truesol_3D_Ned(REAL *val,REAL* x,REAL time) {
  // 3D - curl curl
  val[0] = cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
  val[1] = sin(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]);
  val[2] = -sin(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*x[2]);
}
void truesol_2D_RT(REAL *val,REAL* x,REAL time) {
  // 2D - grad div
  val[0] = sin(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = cos(M_PI*x[0])*sin(M_PI*x[1]);
}
void truesol_3D_RT(REAL *val,REAL* x,REAL time) {
  // 3D - grad div
  val[0] = sin(M_PI*x[0])*cos(M_PI*x[1])*cos(M_PI*x[2]);
  val[1] = cos(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*x[2]);
  val[2] = cos(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]);
}

// Derivative of True Solution (if you have one)
void D_truesol_2D_PX(REAL *val,REAL* x,REAL time) {
  // 2D - grad grad
  val[0] = M_PI*cos(M_PI*x[0])*sin(M_PI*x[1]);
  val[1] = M_PI*sin(M_PI*x[0])*cos(M_PI*x[1]);
}
void D_truesol_3D_PX(REAL *val,REAL* x,REAL time) {
  // 3D - grad grad
  val[0] = M_PI*cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
  val[1] = M_PI*sin(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]);
  val[2] = M_PI*sin(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*x[2]);
}
void D_truesol_2D_Ned(REAL *val,REAL* x,REAL time) {
  // 2D - curl curl
  *val = -2*M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
}
void D_truesol_3D_Ned(REAL *val,REAL* x,REAL time) {
  // 3D - curl curl
  val[0] = -2*M_PI*sin(M_PI*x[0])*cos(M_PI*x[1])*cos(M_PI*x[2]);
  val[1] = 2*M_PI*cos(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*x[2]);
  val[2] = 0;
}
void D_truesol_2D_RT(REAL *val,REAL* x,REAL time) {
  // 2D - grad div
  *val = 2*M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
}
void D_truesol_3D_RT(REAL *val,REAL* x,REAL time) {
  // 3D - grad div
  *val = 3*M_PI*cos(M_PI*x[0])*cos(M_PI*x[1])*cos(M_PI*x[2]);
}

// Right-hand Side
void rhs_2D_PX(REAL *val,REAL* x,REAL time) {
  // 2D - grad grad
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu=-666.6;
  reaction_coeff(&myc,x,time);
  diffusion_coeff(&mya,x,time);
  truesol_2D_PX(&myu,x,time);
  *val = (mya*2*M_PI*M_PI + myc)*myu;
}
void rhs_3D_PX(REAL *val,REAL* x,REAL time) {
  // 3D - grad grad
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu=-666.6;
  reaction_coeff(&myc,x,time);
  diffusion_coeff(&mya,x,time);
  truesol_3D_PX(&myu,x,time);
  *val = (mya*3*M_PI*M_PI + myc)*myu;
}
void rhs_2D_Ned(REAL *val,REAL* x,REAL time) {
  // 2D - curl curl
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu[2];
  reaction_coeff(&myc,x,time);
  diffusion_coeff(&mya,x,time);
  truesol_2D_Ned(myu,x,time);
  val[0] = (mya*2.0*M_PI*M_PI + myc)*myu[0];
  val[1] = (mya*2.0*M_PI*M_PI + myc)*myu[1];
}
void rhs_3D_Ned(REAL *val,REAL* x,REAL time) {
  // 3D - curl curl
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu[3];
  reaction_coeff(&myc,x,time);
  diffusion_coeff(&mya,x,time);
  truesol_3D_Ned(myu,x,time);
  val[0] = (2*mya*M_PI*M_PI + myc)*myu[0];
  val[1] = (2*mya*M_PI*M_PI + myc)*myu[1];
  val[2] = (4*mya*M_PI*M_PI + myc)*myu[2];
}
void rhs_2D_RT(REAL *val,REAL* x,REAL time) {
  // 2D - grad div
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu[2];
  reaction_coeff(&myc,x,time);
  diffusion_coeff(&mya,x,time);
  truesol_2D_RT(myu,x,time);
  val[0] = (mya*2.0*M_PI*M_PI + myc)*myu[0];
  val[1] = (mya*2.0*M_PI*M_PI + myc)*myu[1];
}
void rhs_3D_RT(REAL *val,REAL* x,REAL time) {
  // 3D - grad div
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu[3];
  reaction_coeff(&myc,x,time);
  diffusion_coeff(&mya,x,time);
  truesol_3D_RT(myu,x,time);
  val[0] = (mya*3.0*M_PI*M_PI + myc)*myu[0];
  val[1] = (mya*3.0*M_PI*M_PI + myc)*myu[1];
  val[2] = (mya*3.0*M_PI*M_PI + myc)*myu[2];
}

// Boundary Conditions
void bc_2D_PX(REAL *val,REAL* x,REAL time) {
  REAL myu;
  truesol_2D_PX(&myu,x,time);
  *val= myu;
}
void bc_3D_PX(REAL *val,REAL* x,REAL time) {
  REAL myu;
  truesol_3D_PX(&myu,x,time);
  *val= myu;
}
void bc_2D_Ned(REAL *val,REAL* x,REAL time) {
  REAL myu[2];
  truesol_2D_Ned(myu,x,time);
  val[0] = myu[0];
  val[1] = myu[1];
}
void bc_3D_Ned(REAL *val,REAL* x,REAL time) {
  REAL myu[3];
  truesol_3D_Ned(myu,x,time);
  val[0] = myu[0];
  val[1] = myu[1];
  val[2] = myu[2];
}
void bc_2D_RT(REAL *val,REAL* x,REAL time) {
  REAL myu[2];
  truesol_2D_RT(myu,x,time);
  val[0] = myu[0];
  val[1] = myu[1];
}
void bc_3D_RT(REAL *val,REAL* x,REAL time) {
  REAL myu[3];
  truesol_3D_RT(myu,x,time);
  val[0] = myu[0];
  val[1] = myu[1];
  val[2] = myu[2];
}
/**********************************************************************************/

/****** MAIN DRIVER ***************************************************************/
int main (int argc, char* argv[])
{
  
  printf("\n===========================================================================\n");
  printf("Beginning Program to solve H(D) problem: <D u, D v> + <u,v> = <f,v>.\n");
  printf("===========================================================================\n");
  
  /****** INITIALIZE PARAMETERS **************************************************/
  // flag for errors
  SHORT status;

  // Timing Parameters
  clock_t clk_start,clk_end,clk1,clk2;
  clk_start = clock();
    
  // Set Parameters from Reading in Input File
  input_param     inparam;
  param_input_init(&inparam);
  param_input("./input.dat", &inparam); 
    
  // Open gridfile for reading
  printf("\nCreating mesh and FEM spaces:\n");
  FILE* gfid = HAZ_fopen(inparam.gridfile,"r");
    
  // Dimension is needed for all this to work
  INT dim = inparam.dim;
    
  // Create the mesh (now we assume triangles in 2D or tetrahedra in 3D)
  // File types possible are 0 - old format; 1 - vtk format (doesn't work yet)
  INT mesh_type = 0;
  clk1 = clock();
  trimesh mesh;
  printf(" --> loading grid from file: %s\n",inparam.gridfile);
  creategrid_fread(gfid,mesh_type,&mesh);
  fclose(gfid);

  // Dump mesh for testing
  char* namevtk = "output/mesh.vtu";
  dump_mesh_vtk(namevtk,&mesh);
    
  // Get Quadrature Nodes for the Mesh
  INT nq1d = inparam.nquad; /* Quadrature points per dimension */
  qcoordinates *cq = get_quadrature(&mesh,nq1d);
    
  // Get info for and create FEM spaces
  // Order of Elements: 0 - P0; 1 - P1; 2 - P2; -1 - Nedelec; -2 - Raviart-Thomas
  INT order = inparam.FE_type;
  fespace FE;
  create_fespace(&FE,&mesh,order);
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

if(inparam.print_level > 3) {
    char varu[10];
    char dir[20];

    sprintf(dir,"output");
    sprintf(varu,"u");

    dump_fespace(&FE,varu,dir);
  }
    
  clk2 = clock();
  printf(" --> elapsed CPU time for mesh and FEM space construction = %lf seconds.\n\n", \
	 (REAL) (clk2 - clk1)/CLOCKS_PER_SEC);
  /*******************************************************************************/
    
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
  dCSRmat Diff;
  dCSRmat A;
  dCSRmat Mass;
    
  // Assemble the matrix without BC
  // Diffusion block
  if(dim==2) {
    if(FE.FEtype>0 && FE.FEtype<10) { // PX
      assemble_global(&Diff,&b,assemble_DuDv_local,&FE,&mesh,cq,rhs_2D_PX,diffusion_coeff,0.0);
    } else if(FE.FEtype==20) { // Nedelec
      assemble_global(&Diff,&b,assemble_DuDv_local,&FE,&mesh,cq,rhs_2D_Ned,diffusion_coeff,0.0);
    } else if(FE.FEtype==30) { // RT
      assemble_global(&Diff,&b,assemble_DuDv_local,&FE,&mesh,cq,rhs_2D_RT,diffusion_coeff,0.0);
    } else {
      printf("Unsure of what elements you are using\n");
      exit(0);
    }
  } else if(dim==3) {
    if(FE.FEtype>0 && FE.FEtype<10) { // PX
      assemble_global(&Diff,&b,assemble_DuDv_local,&FE,&mesh,cq,rhs_3D_PX,diffusion_coeff,0.0);
    } else if(FE.FEtype==20) { // Nedelec
      assemble_global(&Diff,&b,assemble_DuDv_local,&FE,&mesh,cq,rhs_3D_Ned,diffusion_coeff,0.0);
    } else if(FE.FEtype==30) { // RT
      assemble_global(&Diff,&b,assemble_DuDv_local,&FE,&mesh,cq,rhs_3D_RT,diffusion_coeff,0.0);
    } else {
      printf("Unsure of what elements you are using\n");
      exit(0);
    }
  } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
  }

  // Reaction block
  assemble_global(&Mass,NULL,assemble_mass_local,&FE,&mesh,cq,NULL,reaction_coeff,0.0);

  // Add the M + D
  dcsr_add_1(&Diff,1.0,&Mass,1.0,&A);
  dcsr_free(&Diff);
  dcsr_free(&Mass);
  
  // Eliminate Dirichlet BC
  if(dim==2) {
    if(FE.FEtype>0 && FE.FEtype<10) { // PX
      eliminate_DirichletBC(bc_2D_PX,&FE,&mesh,&b,&A,0.0);
    } else if(FE.FEtype==20) { // Nedelec
      eliminate_DirichletBC(bc_2D_Ned,&FE,&mesh,&b,&A,0.0);
    } else if(FE.FEtype==30) { // RT
      eliminate_DirichletBC(bc_2D_RT,&FE,&mesh,&b,&A,0.0);
    } else {
      printf("Unsure of what elements you are using\n");
      exit(0);
    }
  } else if(dim==3) {
    if(FE.FEtype>0 && FE.FEtype<10) { // PX
      eliminate_DirichletBC(bc_3D_PX,&FE,&mesh,&b,&A,0.0);
    } else if(FE.FEtype==20) { // Nedelec
      eliminate_DirichletBC(bc_3D_Ned,&FE,&mesh,&b,&A,0.0);
    } else if(FE.FEtype==30) { // RT
      eliminate_DirichletBC(bc_3D_RT,&FE,&mesh,&b,&A,0.0);
    } else {
      printf("Unsure of what elements you are using\n");
      exit(0);
    }
  } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
  }

  if(inparam.print_level > 3) {
    FILE* matid = HAZ_fopen("output/mat.dat","w");
    csr_print_matlab(matid,&A);
    fclose(matid);
    FILE* rhsid = HAZ_fopen("output/rhs.dat","w");
    dvector_print(rhsid,&b);
    fclose(rhsid);
  }
    
  clk2 = clock();
  printf(" --> elapsed CPU time for assembly = %lf seconds.\n\n",(REAL) (clk2-clk1)/CLOCKS_PER_SEC);
  /*******************************************************************************/
    
  /**************** Solve ********************************************************/
  printf("Solving the System:\n");
  clk1=clock();
  
  // Parameters
  INT solver_flag=-20;
  
  //============= SHOULD THIS BE IF STATEMENTS? =====================//
  // Set parameters for linear iterative methods
  linear_itsolver_param linear_itparam;
  param_linear_solver_set(&linear_itparam, &inparam);
  param_linear_solver_print(&linear_itparam);
    
  // Set parameters for algebriac multigrid methods
  AMG_param amgparam;
  param_amg_init(&amgparam);
  param_amg_set(&amgparam, &inparam);
  param_amg_print(&amgparam);
    
  // Set parameters for ILU methods
  ILU_param iluparam;
  param_ilu_init(&iluparam);
  param_ilu_set(&iluparam, &inparam);
  param_ilu_print(&iluparam);

  // Data for HX preconditioner
  dCSRmat P_curl;
  dCSRmat Grad;
  //=================================================================//
    
  // Allocate the solution
  dvector u = dvec_create(b.row);
    
  // Set initial guess to be all zero
  dvec_set(u.row, &u, 0.0);
    
  // Solve the linear system
  if(linear_itparam.linear_itsolver_type == 0) { // Direct Solver
    printf(" --> using UMFPACK's Direct Solver:\n");
    //solver_flag = directsolve_UMF_symmetric(&A,&b,u.val,linear_itparam.linear_print_level);
  } else { // Iterative Solver
    dcsr_shift(&A, -1);  // shift A
      
    // Use AMG as iterative solver
    if (linear_itparam.linear_itsolver_type == SOLVER_AMG){
      solver_flag = linear_solver_amg(&A, &b, &u, &amgparam);
    } else { // Use Krylov Iterative Solver
      // Determine Preconditioner
      switch (linear_itparam.linear_precond_type) {      
      case PREC_DIAG:  // Diagonal Preconditioner
	solver_flag = linear_solver_dcsr_krylov_diag(&A, &b, &u, &linear_itparam);
	break;        
      case PREC_AMG:  // AMG preconditioner
	solver_flag = linear_solver_dcsr_krylov_amg(&A, &b, &u, &linear_itparam, &amgparam);
	break;            
      case PREC_ILU:  // ILU preconditioner
	solver_flag = linear_solver_dcsr_krylov_ilu(&A, &b, &u, &linear_itparam, &iluparam);
	break;            
      case PREC_HX_CURL_A: // HX precondtioner
	get_Pigrad_H1toNed(&P_curl,&mesh);
	get_grad_H1toNed(&Grad,&mesh);
	dcsr_shift(&P_curl, -1);  // shift A
	dcsr_shift(&Grad, -1);  // shift A      
	solver_flag = linear_solver_dcsr_krylov_hx_curl(&A, &b, &u, &linear_itparam, &amgparam, &P_curl, &Grad);
	dcsr_shift(&P_curl, 1);   // shift A back
	dcsr_shift(&Grad, 1);   // shift A back
	dcsr_free(&P_curl);
	dcsr_free(&Grad);     
	break;       
      case PREC_HX_CURL_M: // HX precondtioner              
	get_Pigrad_H1toNed(&P_curl,&mesh);
	get_grad_H1toNed(&Grad,&mesh);
	dcsr_shift(&P_curl, -1);  // shift A
	dcsr_shift(&Grad, -1);  // shift A      
	solver_flag = linear_solver_dcsr_krylov_hx_curl(&A, &b, &u, &linear_itparam, &amgparam, &P_curl, &Grad);
	dcsr_shift(&P_curl, 1);   // shift A back
	dcsr_shift(&Grad, 1);   // shift A back            
	dcsr_free(&P_curl);
	dcsr_free(&Grad);     
	break;       
      default:  // No Preconditioner
	solver_flag = linear_solver_dcsr_krylov(&A, &b, &u, &linear_itparam);
	break;
      }
    }    
    dcsr_shift(&A, 1);   // shift A back
  }

  // Error Check
  if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n", solver_flag);
    
  clk2=clock();
  printf(" --> elapsed CPU time for solve = %lf seconds.\n\n",(REAL) (clk2-clk1)/CLOCKS_PER_SEC);  
  /******************************************************************************/
    
  /**************** Compute Errors if you have true solution ********************/
  printf("Computing True Solution and Errors:\n");
  clk2 = clock();
  REAL uerr=0.0;
  REAL graduerr=0.0;
  if(dim==2) {
    if(FE.FEtype>0 && FE.FEtype<10) { // PX
      uerr = L2error(u.val,truesol_2D_PX,&FE,&mesh,cq,0.0);
      graduerr = HDsemierror(u.val,D_truesol_2D_PX,&FE,&mesh,cq,0.0);
    } else if(FE.FEtype==20) { // Nedelec
      uerr = L2error(u.val,truesol_2D_Ned,&FE,&mesh,cq,0.0);
      graduerr = HDsemierror(u.val,D_truesol_2D_Ned,&FE,&mesh,cq,0.0);
    } else if(FE.FEtype==30) { // RT
      uerr = L2error(u.val,truesol_2D_RT,&FE,&mesh,cq,0.0);
      graduerr = HDsemierror(u.val,D_truesol_2D_RT,&FE,&mesh,cq,0.0);
    } else {
      printf("Unsure of what elements you are using\n");
      exit(0);
    }
  } else if(dim==3) {
    if(FE.FEtype>0 && FE.FEtype<10) { // PX
      uerr = L2error(u.val,truesol_3D_PX,&FE,&mesh,cq,0.0);
      graduerr = HDsemierror(u.val,D_truesol_3D_PX,&FE,&mesh,cq,0.0);
    } else if(FE.FEtype==20) { // Nedelec
      uerr = L2error(u.val,truesol_3D_Ned,&FE,&mesh,cq,0.0);
      graduerr = HDsemierror(u.val,D_truesol_3D_Ned,&FE,&mesh,cq,0.0);
    } else if(FE.FEtype==30) { // RT
      uerr = L2error(u.val,truesol_3D_RT,&FE,&mesh,cq,0.0);
      graduerr = HDsemierror(u.val,D_truesol_3D_RT,&FE,&mesh,cq,0.0);
    } else {
      printf("Unsure of what elements you are using\n");
      exit(0);
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
  printf(" --> elapsed CPU time for getting errors = %lf seconds.\n\n",(REAL) (clk2-clk1)/CLOCKS_PER_SEC);
  /******************************************************************************/
    
  /**************** Print Results or Dump Results *******************************/
  if (inparam.output_type==2) {
    FILE* uid = HAZ_fopen("output/sol.dat","w");
    dvector_print(uid,&u);
    fclose(uid);
  }
  /******************************************************************************/
    
  /******** Free All the Arrays *************************************************/
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
  /******************************************************************************/
    
  clk_end = clock();
  printf("\nEnd of Program: Total CPU Time = %lf seconds.\n\n",(REAL) (clk_end-clk_start)/CLOCKS_PER_SEC);
  return 0;
    
}	/* End of Program */
/**********************************************************************************/


