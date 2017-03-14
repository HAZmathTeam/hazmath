/*! \file examples/ConvectionDiffusion/ConvectionDiffusion.c
 *
 *  Created by Xiaozhe Hu, James Adler, and Ludmil Zikatanov 2017-03-09
 *  (originally 1994-02-02)  
 *
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program solves the following PDE using EAFE FE discretization: 
 *
 *        -div(a(x)*grad(u) - u*b) = f
 *
 *        b is an advection vector; a(x) is a diffusion matrix. 
 *        u = 0 on the boundary.
 *        
 *        Details about the EAFE discretization are found in: 
 *        Jinchao Xu and Ludmil Zikatanov: A monotone finite element
 *        scheme for convection-diffusion equations. Math. Comp. 68
 *        (1999), no. 228, 1429–1446.
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
/* local include */
#include "ConvectionDiffusion.h"
/*********************************************************************/
int main (int argc, char* argv[])
{
  
  printf("\n===========================================================================\n");
  printf("Beginning Program to solve a Convection Diffusion Equation.\n");
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
  clock_t clk_mesh_start = clock(); // Time mesh generation FE setup
  INT mesh_type = 0;
  trimesh mesh;
  printf(" --> loading grid from file: %s\n",inparam.gridfile);
  creategrid_fread(gfid,mesh_type,&mesh);
  fclose(gfid);

  // Dimension is needed for all this to work
  INT dim = mesh.dim;

  // Get Quadrature Nodes for the Mesh
  INT nq1d = inparam.nquad; // Quadrature points per dimension
  qcoordinates *cq = get_quadrature(&mesh,nq1d);

  // Get info for and create FEM spaces
  // Order of Elements: P1 for this. 
  INT order = 1; //inparam.FE_type;
  fespace FE;
  create_fespace(&FE,&mesh,order);
  // Strings for printing
  char elmtype[8];
  sprintf(elmtype,"P%d",order);

  // Set Dirichlet Boundaries
  // Assume the physical boundaries (flag of 1 in mesh file) are Dirichlet
  set_dirichlet_bdry(&FE,&mesh,1);

  // Dump some of the data
  if(inparam.print_level > 3) { // && inparam.output_dir != NULL) {
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
    
  // Assemble the matrix with natural BC. 
  // Diffusion block
  /*  
      assemble_global(&A,&b,assemble_DuDv_local,&FE,&mesh,cq,f_rhs,
                  poisson_coeff,0.0);
  */
  eafe(&A,&b,assemble_DuDv_local,			\
       mesh,FE,cq,					\
       diffusion_coeff,f_rhs,advection,bc_any,0.0);
  clock_t clk_assembly_end = clock();
  printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL)
         (clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);
  /*******************************************************************/

  /**************** Solve ********************************************/

  // Create Solution Vector
  dvector sol = dvec_create(FE.ndof);
  dvector exact_sol = dvec_create(FE.ndof);
  FE_Evaluate(exact_sol.val,exactsol,&FE,&mesh,0.0);

  // Set parameters for linear iterative methods
  linear_itsolver_param linear_itparam;
  param_linear_solver_init(&linear_itparam);
  param_linear_solver_set(&linear_itparam, &inparam);
  
  // Set parameters for AMG methods
  AMG_param amgparam;
  param_amg_init(&amgparam);
  param_amg_set(&amgparam, &inparam);
  
  INT solver_flag=-20;
  

  eliminate_DirichletBC(&bc_any,&FE,&mesh,&b,&A,0.0);
  // Solve
  clock_t clk_solve_start = clock();
  dcsr_shift(&A, -1);  // shift A
  switch (linear_itparam.linear_itsolver_type) {      
  case 3:
      
      switch (linear_itparam.linear_precond_type) {
       
        case PREC_AMG: // GMRES + AMG precondtioner
          solver_flag = linear_solver_dcsr_krylov_amg(&A, &b, &sol, &linear_itparam, &amgparam);
          break;
          
        default: // GMRES+ No Preconditioner
          solver_flag = linear_solver_dcsr_krylov(&A,&b,&sol,&linear_itparam);
          break;
      }
      
    dcsr_shift(&A, 1);   // shift A back
    // Error Check
    if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n", solver_flag);
    break;        
  case 4:  ; //empty first statement
    //#ifdef MGRAPH
    // Multigraph preconditioner
    fprintf(stdout, "\nSolve using multigraph\n");
    INT *ka = NULL;
    INT *jareb=NULL;
    REAL *areb=NULL;
    INT idoilu = 1; 
    mgraph_wrap(idoilu, A.row,A.IA,A.JA,A.val,b.val,sol.val,jareb,areb,ka);
    if(ka) free(ka);
    if(jareb) free(jareb);
    if(areb) free(areb);
    //#endif
    break;                        
  default:  //GMRES+ No Preconditioner
    solver_flag = linear_solver_dcsr_krylov(&A,&b,&sol,&linear_itparam);
    // Error Check
    if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n", solver_flag);
    break;
  }
  dcsr_shift(&A, 1);   // shift A back
   


  clock_t clk_solve_end = clock();
  printf("Elapsed CPU Time for Solve = %f seconds.\n\n",(REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);
  // Dump Solution
  // Compute norms of errors
  REAL uerr = L2error(sol.val,exactsol,&FE,&mesh,cq,0.0);
  REAL unorm = L2norm(sol.val,&FE,&mesh,cq);
  printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
  printf("L2 Norm of u            = %26.13e\n",unorm);
  printf("L2 Norm of u error      = %26.13e\n",uerr);
  printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n\n");
  char solout[40];
  char exactout[40];
  //  if (inparam.output_dir!=NULL) { always true. 
    sprintf(solout,"output/solution_ts%03d.vtu",0);
    dump_sol_onV_vtk(solout,&mesh,sol.val,1);
    sprintf(exactout,"output/exact_solution_ts%03d.vtu",0);
    dump_sol_onV_vtk(exactout,&mesh,exact_sol.val,1);
    //  }
  printf("\n");
  /*******************************************************************/
  dvec_free(&exact_sol);
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
