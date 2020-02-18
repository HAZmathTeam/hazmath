/*! \file Maxwell.c
 *
 *  Copyright 2016_HAZMAT__. All rights reserved.
 *
 *  \brief This program solves the Maxwell's equation system for disappearing solutions
 *    using Raviart-Thomas, Nedelec, and P1 elements:
 *
 *      eps* dE/dt - curl(1/mu B) = -j,
 *      dB/dt + curl(E) = 0,
 *      div(eps*E) = 0,
 *      div(B) = 0.
 *
 *    This is solved in a 3D rectangular region in the plane.
 *
 *    Along the boundary of the region perfect conductor conditions:
 *
 *      Exterior: n * B = 0, n x E = f(x,t),
 *      *
 * \note This example illustrates how to combine finite-element spaces into a block system
 *       for a mixed-FEM problem. 
 *
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "hazmath.h"
#include "Maxwell_data.h"
/*********************************************************************************/

/****** MAIN DRIVER **************************************************************/
int main (int argc, char* argv[]) 
{

  printf("\n===========================================================================\n");
  printf("\nBeginning Program to Solve Maxwell's Equations with Disappearing Solutions\n");
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

  // Create the mesh (now we assume triangles in 2D or tetrahedra in 3D)
  // File types possible are 0 - old format; 1 - vtk format
  clock_t clk_mesh_start = clock(); // Time mesh generation FE setup
  INT mesh_type = 0;
  mesh_struct mesh;
  printf(" --> loading grid from file: %s\n",inparam.gridfile);
  creategrid_fread(gfid,mesh_type,&mesh);
  fclose(gfid);
  INT dim = mesh.dim;

  
  // Get Quadrature Nodes for the Mesh
  INT nq1d = inparam.nquad; /* Quadrature points per dimension */
  qcoordinates *cq = get_quadrature(&mesh,nq1d);

  // Time stepping parameters
  block_timestepper time_stepper;
  
  // Get info for and create FEM spaces
  // Order of Elements: 0 - P0; 1 - P1; 2 - P2; 20 - Nedelec; 30 - Raviart-Thomas
  INT order_E = 20;
  INT order_B = 30;
  INT order_p = 1;
  fespace FE_E, FE_B, FE_p;
  create_fespace(&FE_E,&mesh,order_E);
  create_fespace(&FE_B,&mesh,order_B);
  create_fespace(&FE_p,&mesh,order_p);

  if(inparam.print_level > 3) {
    // Dump Mesh
    char* namevtk = "output/mesh.vtu";
    dump_mesh_vtk(namevtk,&mesh);
  }
  
  // Set boundary conditions
  set_dirichlet_bdry(&FE_E,&mesh,1,1);
  set_dirichlet_bdry(&FE_B,&mesh,1,1);
  set_dirichlet_bdry(&FE_p,&mesh,1,1);

  // Create Block System with ordering (E,B,p)
  INT ndof = FE_E.ndof + FE_B.ndof + FE_p.ndof;
  // Get Global FE Space
  block_fespace FE;
  FE.ndof = ndof;
  FE.nbdof = FE_E.nbdof + FE_B.nbdof + FE_p.nbdof;
  FE.nspaces = 3;
  FE.nun = 2*dim+1;
  FE.var_spaces = (fespace **) calloc(FE.nspaces,sizeof(fespace *));
  FE.var_spaces[0] = &FE_E;
  FE.var_spaces[1] = &FE_B;
  FE.var_spaces[2] = &FE_p;
  set_dirichlet_bdry_block(&FE,&mesh);

  clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
  printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
         (REAL) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
  /*******************************************************************************/

  printf("***********************************************************************************\n");
  printf("Element Type: Nedelec for E, Raviart-Thomas for B, and Linears for p.\n");
  printf("Number of Elements = %d\tOrder of Quadrature = %d\n",mesh.nelm,2*nq1d-1);
  printf("\n\t--- Degrees of Freedom ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nv,mesh.nedge,mesh.nface);
  printf("\t--> DOF: %d\n",FE_E.ndof+FE_B.ndof+FE_p.ndof);
  printf("\n\t--- Boundaries ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nbv,mesh.nbedge,mesh.nbface);
  printf("\t--> Boundary DOF: %d\n",FE_E.nbdof+FE_B.nbdof+FE_p.ndof);
  printf("***********************************************************************************\n\n");

  /*** Assemble the matrix and right hand side *******************************/
  // Here we assemble the discrete system (notice a slack variable, p, is used for p in H^1_0:
  //		<eps dE/dt,F> + <eps grad p,F> - <(1/mu)B,curl F> = -<j,F>
  //            <dp/dt,q> - <eps E, grad q> = 0
  //            <(1/mu) dB/dt,C> + <(1/mu) curl E, C> = 0
  //
  //            p are H1 elements with 0 Dirichlet boundaries
  //            B is RT elements with div B = 0, and B*n = 0 on boundary 
  //            E is Nedelec elements with n x E = f(x,t) on boundary 

  // We do this using only mass and incidence matrices
  //  Me = <eps E,F>   Mf = <1/mu B,C>   Mn = <p,q>
  //  De = diag(|e|)    Df = diag(|f|)
  //  G = De^(-1) ed_v    K = Df^(-1) face_ed De
  //
  
  printf("Assembling the matrix and right-hand side:\n");
  clock_t clk_assembly_start = clock();

  // Allocate the right-hand side and declare the csr matrices
  dvector b_E;  /* RHS vector for E-equation -> <-j,F> */
  dvector b_B;  /* RHS vector for B-equation -> 0 */
  dvector b_p;  /* RHS vector for p-equation -> 0 */
  dCSRmat K;                 /* Incidence matrix of face to edge map (includes signs) |E|*K*|F|^(-1) = Curl Operator */
  dCSRmat G;                 /* Incidence matrix of edge to node map (includes signs) |E|^(-1)*Ggrad = grad operator */
  dCSRmat Me;                /* Mass Matrix for Edge DOF */
  dCSRmat Mf;                /* Mass Matrix for Face DOF */
  dCSRmat Mv;                /* Mass Matrix for Vertex DOF */
  dCSRmat MG;                /* Me*G */
  dCSRmat MGt;               /* G'*Me */
  dCSRmat MK;                /* Mf*K */
  dCSRmat MKt;               /* K'*Mf */

  // Assemble the matrices without BC first

  // Mass matrices
  assemble_global(&Me,&b_E,assemble_mass_local,&FE_E,&mesh,cq,current_density,permitivity,0.0);
  assemble_global(&Mf,&b_B,assemble_mass_local,&FE_B,&mesh,cq,zero_coeff_vec3D,oneovermu,0.0);
  assemble_global(&Mv,&b_p,assemble_mass_local,&FE_p,&mesh,cq,zero_coeff_scal,one_coeff_scal,0.0);
  
  
  coordinates* cv_vor;
  cv_vor =  allocatecoords(mesh.nelm,mesh.dim);
  dvector vor_edge_length = dvec_create(mesh.nface);
  
  dvector vor_face_area = dvec_create(mesh.nedge);
  
  REAL* pt_on_face=(REAL *)calloc(3*mesh.nedge,sizeof(REAL));
  
  
 dvector vor_el_vol = dvec_create(mesh.nv);

 compute_Voronoi_nodes(&mesh, cv_vor);

 compute_Voronoi_edges(&mesh, cv_vor, &vor_edge_length);

 compute_Voronoi_faces(&mesh, cv_vor, pt_on_face, &vor_face_area);

 compute_Voronoi_volumes(&mesh, cv_vor, &vor_face_area, pt_on_face, &vor_el_vol);

 
 //Matrices with mesh info
 dCSRmat Vv;
 dCSRmat Va;
 dCSRmat Ve;
 dCSRmat Dv;
 dCSRmat Da;
 dCSRmat De;
 
 //lumped mass matrices
 dCSRmat MeL;
 dCSRmat MbL;
 dCSRmat MpL;
 
Vv = dcsr_create_diagonal_matrix(&vor_el_vol);
Va = dcsr_create_diagonal_matrix(&vor_face_area);
Ve = dcsr_create_diagonal_matrix(&vor_edge_length);

dvector del_el_vol;
del_el_vol.row = mesh.nelm;
del_el_vol.val = mesh.el_vol;

dvector del_face_area;
del_face_area.row = mesh.nface;
del_face_area.val = mesh.f_area;

dvector del_edge_length;
del_edge_length.row = mesh.nedge;
del_edge_length.val = mesh.ed_len;


Dv = dcsr_create_diagonal_matrix(&del_el_vol);
Da = dcsr_create_diagonal_matrix(&del_face_area);
De = dcsr_create_diagonal_matrix(&del_edge_length);
 
  
  // Grad and Curl matrices
  get_grad_H1toNed(&G,&mesh);
  get_curl_NedtoRT(&K,&mesh);
  

  // Start Combining
  //  A =   0        -K^T Mf          Med G
  //      Mf K         0                0
  //    -G^T Med       0                0
  //
  //  M ->    <eps E,F>
  //         <(1/mu) B,C>
  //            <pt,q>
  //
  //    f ->     <-j,F>
  //               0
  //               0
  // M du/dt + Au = f
  
  
  //create lumped matrices
  dcsr_mxm(&Va,&De,&MeL); //MeL = vor area * del edge
  dcsr_mxm(&Da,&Ve,&MbL); //MbL = del face * vor edge
  MpL = Vv;				//MpL = vor volumes
  
  dcsr_mxm(&MeL,&G,&MG); // Me*G
  dcsr_trans(&MG,&MGt); // G'*Me
  dcsr_axm(&MGt,-1); // -G'*Me
  dcsr_mxm(&MbL,&K,&MK); // Mf*K
  dcsr_trans(&MK,&MKt); // K'*Mf
  dcsr_axm(&MKt,-1); // -K'*Mf
  
/*   dcsr_mxm(&Me,&G,&MG); // Me*G
  dcsr_trans(&MG,&MGt); // G'*Me
  dcsr_axm(&MGt,-1); // -G'*Me
  dcsr_mxm(&Mf,&K,&MK); // Mf*K
  dcsr_trans(&MK,&MKt); // K'*Mf
  dcsr_axm(&MKt,-1); // -K'*Mf */

/*   // Block Matrix M;
  block_dCSRmat Mb;
  bdcsr_alloc_minimal(3,3,&Mb);
  Mb.blocks[0] = &Me;
  Mb.blocks[1] = NULL;
  Mb.blocks[2] = NULL;
  Mb.blocks[3] = NULL;
  Mb.blocks[4] = &Mf;
  Mb.blocks[5] = NULL;
  Mb.blocks[6] = NULL;
  Mb.blocks[7] = NULL;
  Mb.blocks[8] = &Mv; */
  
  
/*     // Block Matrix M;
  block_dCSRmat Mb;
  bdcsr_alloc_minimal(3,3,&Mb);
  Mb.blocks[0] = &Me;
  Mb.blocks[1] = NULL;
  Mb.blocks[2] = NULL;
  Mb.blocks[3] = NULL;
  Mb.blocks[4] = &Mf;
  Mb.blocks[5] = NULL;
  Mb.blocks[6] = NULL;
  Mb.blocks[7] = NULL;
  Mb.blocks[8] = &Mv;

  // Block Matrix A(shifts needed)
  block_dCSRmat Ab;
  bdcsr_alloc_minimal(3,3,&Ab);
  Ab.blocks[0] = NULL;
  Ab.blocks[1] = &MKt;
  Ab.blocks[2] = &MG;
  Ab.blocks[3] = &MK;
  Ab.blocks[4] = NULL;
  Ab.blocks[5] = NULL;
  Ab.blocks[6] = &MGt;
  Ab.blocks[7] = NULL;
  Ab.blocks[8] = NULL; */
  
  
      // Block Matrix M;
  block_dCSRmat Mb;
  bdcsr_alloc_minimal(3,3,&Mb);
  Mb.blocks[0] = &MeL;
  Mb.blocks[1] = NULL;
  Mb.blocks[2] = NULL;
  Mb.blocks[3] = NULL;
  Mb.blocks[4] = &MbL;
  Mb.blocks[5] = NULL;
  Mb.blocks[6] = NULL;
  Mb.blocks[7] = NULL;
  Mb.blocks[8] = &MpL;

  // Block Matrix A(shifts needed)
  block_dCSRmat Ab;
  bdcsr_alloc_minimal(3,3,&Ab);
  Ab.blocks[0] = NULL;
  Ab.blocks[1] = &MKt;
  Ab.blocks[2] = &MG;
  Ab.blocks[3] = &MK;
  Ab.blocks[4] = NULL;
  Ab.blocks[5] = NULL;
  Ab.blocks[6] = &MGt;
  Ab.blocks[7] = NULL;
  Ab.blocks[8] = NULL;

  // Block RHS (need the current RHS, the updated one for time stepping, and the function evaluted at the new time step (if time-dependent))
  dvector b = dvec_create(ndof);
  for(i=0;i<mesh.nedge;i++) b.val[i] = b_E.val[i];
  for(i=0;i<mesh.nface;i++) b.val[i+mesh.nedge] = b_B.val[i];
  for(i=0;i<mesh.nv;i++) b.val[i+mesh.nedge+mesh.nv] = b_p.val[i];
  
  
  /*
  
  Scale before timestepping
  
  */

  // Create Time Operator
  initialize_blktimestepper(&time_stepper,&inparam,1,FE.ndof,FE.nspaces);
  time_stepper.A = &Ab;
  time_stepper.M = &Mb;
  time_stepper.Ldata = &Ab;
  get_blktimeoperator(&time_stepper,1,1);

  clock_t clk_assembly_end = clock();
  printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL)
         (clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/


  /**************** Perform Timestepping ********************************************************************/
  printf("===========================================================================\n");
  printf("The system is time-dependent.  Performing %d Time Step(s) of %s\n",time_stepper.tsteps,time_stepper.time_scheme_str);
  printf("===========================================================================\n\n");
  clock_t clk1 = clock();

  // Set up Linear Solver
  INT solver_flag=-20;

  // Set parameters for linear iterative methods
  linear_itsolver_param linear_itparam;
  param_linear_solver_init(&linear_itparam);
  param_linear_solver_set(&linear_itparam, &inparam);
  

  
  // Set parameters for algebriac multigrid methods
  AMG_param amgparam;
  param_amg_init(&amgparam);
  param_amg_set(&amgparam, &inparam);

  // Allocate the solution: u = (E, B, p)^T
  dvector uE = dvec_create(mesh.nedge);
  dvector uB = dvec_create(mesh.nface);
  dvector up = dvec_create(mesh.nv);
  dvector u = dvec_create(ndof);
  dvector tsol = dvec_create(ndof);
  
  // Initial Condition and RHS
  blockFE_Evaluate(u.val,truesol,&FE,&mesh,0.0);
  
  char** varname = malloc(50*FE.nspaces*sizeof(char*));
  varname[0] = "E";
  varname[1] = "B";
  varname[2] = "p";
  char trueout[40]; 
  
  // Get components for error computing
  get_unknown_component(&uE,&u,&FE,0);
/*   for (i=0;i<uE.row;i++) {
	  printf("uE[%d] = %f, mid_x=%f, mid_y=%f, mid_z=%f, tan_x = %f, tan_y = %f, tan_z = %f\n", i, uE.val[i], mesh.ed_mid[i*dim+0], mesh.ed_mid[i*dim+1], mesh.ed_mid[i*dim+2], mesh.ed_tau[i*dim+0], mesh.ed_tau[i*dim+1], mesh.ed_tau[i*dim+2]);
	  } 
   */
  
  get_unknown_component(&uB,&u,&FE,1);
  get_unknown_component(&up,&u,&FE,2);
  

  // Set initial condition and RHS in timestepper
  time_stepper.sol = &u;
  time_stepper.rhs = &b;

  // Dump Initial Condition
  char solout[40];


  if(inparam.output_dir!=NULL) {
    sprintf(solout,"output/solution_ts000.vtu");
    dump_blocksol_vtk(solout,varname,&mesh,&FE,u.val);
	sprintf(trueout,"output/truesol_ts000.vtu");
    dump_blocksol_vtk(trueout,varname,&mesh,&FE,u.val);
  }
  

  printf("*******************************************\n");
  printf("Initial Condition:\t Actual Time = 0.0\n");
  printf("*******************************************\n\n");
  // Print L2 norms and Errors
  REAL Eerr = L2error(uE.val,Etrue,&FE_E,&mesh,cq,0.0);
  REAL Berr = L2error(uB.val,Btrue,&FE_B,&mesh,cq,0.0);
  REAL perr = L2error(up.val,ptrue,&FE_p,&mesh,cq,0.0);
  REAL EL2 = L2norm(uE.val,&FE_E,&mesh,cq);
  REAL BL2 = L2norm(uB.val,&FE_B,&mesh,cq);
  REAL pL2 = L2norm(up.val,&FE_p,&mesh,cq);
  printf("\n------------- Norms and Errors -----------------------------\n");
  printf("||E||_0 = %25.16e\n",EL2);
  printf("||B||_0 = %25.16e\n",BL2);
  printf("||p||_0 = %25.16e\n",pL2);
  printf("||E-Etrue||_0 = %25.16e\n",Eerr);
  printf("||B-Btrue||_0 = %25.16e\n",Berr);
  printf("||p-ptrue||_0 = %25.16e\n",perr);
  printf("------------------------------------------------------------\n");
  


  clock_t clk_timeloop_start = clock();
  // Begin Time Stepping
  for (i=0; i<time_stepper.tsteps; i++) {
    clock_t clk_timestep_start = clock();
    // Update Time Step Data (includes time, counters, solution, and rhs)
    update_blktimestep(&time_stepper);

    printf("\n*******************************************\n");
    printf("TIME STEP: %d\t Actual Time = %f\n",i+1,time_stepper.time);
    printf("*******************************************\n\n");

    // Update RHS
    if(time_stepper.rhs_timedep) { // If RHS is time-dependent, get new version
      assemble_global_RHS(&b_E,&FE_E,&mesh,cq,current_density,time_stepper.time);	  
      assemble_global_RHS(&b_B,&FE_B,&mesh,cq,zero_coeff_vec3D,time_stepper.time);
      assemble_global_RHS(&b_p,&FE_p,&mesh,cq,zero_coeff_scal,time_stepper.time);
      for(j=0;j<mesh.nedge;j++) time_stepper.rhs->val[j] = b_E.val[j];
      for(j=0;j<mesh.nface;j++) time_stepper.rhs->val[j+mesh.nedge] = b_B.val[j];
      for(j=0;j<mesh.nv;j++) time_stepper.rhs->val[j+mesh.nedge+mesh.nface] = b_p.val[j];
    }

    // Update RHS
    update_blktime_rhs(&time_stepper);

    // For first time step eliminate boundary conditions in matrix and rhs
    if(i==0) {
      eliminate_DirichletBC_blockFE_blockA(bc,&FE,&mesh,time_stepper.rhs_time,time_stepper.At,time_stepper.time);
    } else {
      eliminate_DirichletBC_RHS_blockFE_blockA(bc,&FE,&mesh,time_stepper.rhs_time,time_stepper.At_noBC,time_stepper.time);
    }

    // Solve
    clock_t clk_solve_start = clock();
    //bdcsr_shift(time_stepper.At,-1);
	
	//solver_flag = block_directsolve_UMF(time_stepper.At, time_stepper.rhs_time,time_stepper.sol, 3);
	
    solver_flag = linear_solver_bdcsr_krylov(time_stepper.At, time_stepper.rhs_time, time_stepper.sol, &linear_itparam);
    //bdcsr_shift(time_stepper.At,1);
    clock_t clk_solve_end = clock();
    printf("Elapsed CPU Time for Solve = %f seconds.\n\n",(REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);

//    // Update time steps
    get_unknown_component(&uE,time_stepper.sol,&FE,0);
    get_unknown_component(&uB,time_stepper.sol,&FE,1);
    get_unknown_component(&up,time_stepper.sol,&FE,2);

	REAL plinf;
    // Compute Errors
    printf("\n------------- Norms and Errors -----------------------------\n");
    Eerr = L2error(uE.val,Etrue,&FE_E,&mesh,cq,time_stepper.time);
    Berr = L2error(uB.val,Btrue,&FE_B,&mesh,cq,time_stepper.time);
    perr = L2error(up.val,ptrue,&FE_p,&mesh,cq,time_stepper.time);
    EL2 = L2norm(uE.val,&FE_E,&mesh,cq);
    BL2 = L2norm(uB.val,&FE_B,&mesh,cq);
    pL2 = L2norm(up.val,&FE_p,&mesh,cq);
	plinf = dvec_norminf(&up);
    printf("||E||_0 = %25.16e\n",EL2);
    printf("||B||_0 = %25.16e\n",BL2);
    printf("||p||_0 = %25.16e\n",pL2);
    printf("||E-Etrue||_0 = %25.16e\n",Eerr);
    printf("||B-Btrue||_0 = %25.16e\n",Berr);
    printf("||p-ptrue||_0 = %25.16e\n",perr);
	printf("||p||_inf = %25.16e\n",plinf);
    printf("------------------------------------------------------------\n");
    //fprintf(nid,"%d&%1.3f&%1.3f&%1.3f\\\\ \n",i+1,EL2,BL2,pL2);

	if (inparam.output_dir!=NULL) {
      sprintf(solout,"output/solution_ts%03d.vtu",time_stepper.current_step);
      dump_blocksol_vtk(solout,varname,&mesh,&FE,time_stepper.sol->val);
      blockFE_Evaluate(tsol.val,truesol,&FE,&mesh,time_stepper.time);
      sprintf(trueout,"output/truesol_ts%03d.vtu",time_stepper.current_step);
      dump_blocksol_vtk(trueout,varname,&mesh,&FE,tsol.val);
    }

	
    clock_t clk_timestep_end = clock();
    printf("Elapsed CPU Time for Time Step = %f seconds.\n\n",(REAL) (clk_timestep_end-clk_timestep_start)/CLOCKS_PER_SEC);
  } // End of Timestepping Loop
  
  clock_t clk_timeloop_end = clock();
  printf("Elapsed CPU Time ALL Time Steps = %f seconds.\n\n",(REAL) (clk_timeloop_end-clk_timeloop_start)/CLOCKS_PER_SEC);

 if (inparam.output_dir!=NULL) {
    create_pvd("output/solution.pvd",time_stepper.tsteps+1,"solution_ts","timestep");
    create_pvd("output/truesol.pvd",time_stepper.tsteps+1,"truesol_ts","timestep");
  }
  /*******************************************************************************************/


  /******** Free All the Arrays *************************************************************/

 // Time Stepper
  free_blktimestepper(&time_stepper, 0);

  // CSR Matrices
  dcsr_free(&K);
  dcsr_free(&G);
  dcsr_free(&MG);
  dcsr_free(&MGt);
  dcsr_free(&MK);
  dcsr_free(&MKt);
  
  // Vectors
  dvec_free(&b_E);
  dvec_free(&b_B);
  dvec_free(&b_p);
  dvec_free(&uE);
  dvec_free(&uB);
  dvec_free(&up);
  
  // FE Spaces
  free_fespace(&FE_E);
  free_fespace(&FE_B);
  free_fespace(&FE_p);
  free_blockfespace(&FE);
  
  // Quadrature
  if(cq) {
    free_qcoords(cq);
    free(cq);
    cq = NULL;
  }
  
  // Mesh
  free_mesh(&mesh);

  // Arrays
  if(varname) free(varname);
  
  
  /*****************************************************************************************/

  clock_t clk_overall_end = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",(REAL) (clk_overall_end-clk_overall_start)/CLOCKS_PER_SEC);
  return 0;

}	/* End of Program */
/*******************************************************************************************/
