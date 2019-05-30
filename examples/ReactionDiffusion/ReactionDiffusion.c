/*! \file examples/ReactionDiffusion/ReactionDiffusion.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 20150109.
 *  Edited by Casey Cavanaugh 20181028.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program solves the reaction-diffusion equation
 *
 *        -div(A(x)grad(u)) + c(x)u = f(x)
 *
 *        in 1D, 2D, or 3D using Lagrange finite elements.
 *
 *        Along the boundary of the region, Dirichlet conditions are
 *        imposed.
 *
 */
 
/*********** HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
#include "ReactionDiffusionData.h"
#include "ReactionDiffusionSystem.h"
/*********************************************************************/

/****** MAIN DRIVER **************************************************/
int main (int argc, char* argv[])
{
  
  printf("\n===========================================================================\n");
  printf("Beginning Program to solve  <a D u, D v> + <c u,v> = <f,v>.\n");
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
  //  1 - P1; 2 - P2; 
  INT order = inparam.FE_type;
  fespace FE;
  create_fespace(&FE,&mesh,order);

  // Set Dirichlet Boundaries
  // Assume physical boundaries (flag of 1 in mesh file) are Dirichlet
  set_dirichlet_bdry(&FE,&mesh,1,1);
  
  // Strings for printing
  char elmtype[8];
  if(order>0 && order<10) {
    sprintf(elmtype,"P%d",order);
    printf(" --> using P%d elements \n",order);
  } else {
    printf("ERROR: Unsupported Finite Element Type for this Problem\n");
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
   dCSRmat A;

    
   // Assemble the entire LHS without BC
   //Last argument of assemble global is unused extra parameter (usually time)
   if(FE.FEtype>0 && FE.FEtype<10) { // PX
    if(dim==1) {
		assemble_global(&A,&b,assemble_local_diffusion,&FE,&mesh,cq,
							rhs_1D,pde_coeff,0.0);
    } else if(dim==2) {
        assemble_global(&A,&b,assemble_local_diffusion,&FE,&mesh,cq,
							  rhs_2D,pde_coeff,0.0);
	} else if(dim==3) {
        assemble_global(&A,&b,assemble_local_diffusion,&FE,&mesh,cq,
							  rhs_3D,pde_coeff,0.0);
    } else {
	  status = ERROR_DIM;
      check_error(status, __FUNCTION__);
	}
  } else {
      status = ERROR_FE_TYPE;
      printf("ERROR: Finite element type not valid \n");
  }
 //**************************************************************

	
  // Eliminate Dirichlet BC
  // Different cases for dimension 
  if(dim==1) {
		eliminate_DirichletBC(bc_1D,&FE,&mesh,&b,&A,0.0);
  } else if(dim==2) {
		eliminate_DirichletBC(bc_2D,&FE,&mesh,&b,&A,0.0);
  } else if(dim==3) {
		eliminate_DirichletBC(bc_3D,&FE,&mesh,&b,&A,0.0);
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
  AMG_param amgparam;  
    
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
		// Set parameters for algebriac multigrid methods
		param_amg_init(&amgparam);
		param_amg_set(&amgparam, &inparam);
        solver_flag = linear_solver_amg(&A, &b, &u, &amgparam);
    } else { // Use Krylov Iterative Solver
      // Determine Preconditioner
      // Diagonal preconditioner
      if (linear_itparam.linear_precond_type == PREC_DIAG) {
          solver_flag = linear_solver_dcsr_krylov_diag(&A, &b, &u, &linear_itparam);
      }
      // AMG preconditioner
      else if (linear_itparam.linear_precond_type == PREC_AMG){
		  // Set parameters for algebriac multigrid methods
		  param_amg_init(&amgparam);
		  param_amg_set(&amgparam, &inparam);
		  solver_flag = linear_solver_dcsr_krylov_amg(&A, &b, &u, &linear_itparam, &amgparam);
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
  // Again this depends on dimension
  printf("Computing Exact Solution and Errors:\n");
  clock_t clk_error_start = clock();
  REAL uL2err=0.0;
  REAL uH1err=0.0;
  REAL energyerr=0.0;
  
   if(dim==1) {
		uL2err = L2error(u.val,exactsol_1D,&FE,&mesh,cq,0.0);
		uH1err = HDerror(u.val,exactsol_1D,D_exactsol_1D,&FE,&mesh,cq,0.0);
		energyerr = energyerror(u.val, exactsol_1D, D_exactsol_1D, pde_coeff, &FE,&mesh,cq,0.0);

  } else if(dim==2) {
		uL2err = L2error(u.val,exactsol_2D,&FE,&mesh,cq,0.0);
		uH1err = HDerror(u.val,exactsol_2D,D_exactsol_2D,&FE,&mesh,cq,0.0);
		energyerr = energyerror(u.val, exactsol_2D, D_exactsol_2D, pde_coeff, &FE,&mesh,cq,0.0);

  } else if(dim==3) {
		uL2err = L2error(u.val,exactsol_3D,&FE,&mesh,cq,0.0);
		uH1err = HDerror(u.val,exactsol_3D,D_exactsol_3D,&FE,&mesh,cq,0.0);
		energyerr = energyerror(u.val, exactsol_3D, D_exactsol_3D, pde_coeff, &FE,&mesh,cq,0.0);

  } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
  }

  
  printf("************************************************************************************\n");
  printf("L2 Norm of u error      = %26.13e\n",uL2err);
  printf("H1 Norm of u error      = %26.13e\n",uH1err);
  printf("Energy Norm of u error  = %26.13e\n",energyerr);
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
	
  if(dim==1) {
	    FE_Evaluate(exact_sol.val,exactsol_1D,&FE,&mesh,0.0);

  } else if(dim==2) {
	    FE_Evaluate(exact_sol.val,exactsol_2D,&FE,&mesh,0.0);

  } else if(dim==3) {
		    FE_Evaluate(exact_sol.val,exactsol_3D,&FE,&mesh,0.0);

  } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
  }

   
 
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

    
	