/*!
 *
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program solves Biot's PDE for poroelasticity using finite elements
 * Locking-Free enrichecd Galerkin is used for the mechanics
 * Locally-conservative enriched Galerkin in used for the pressure
 *
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "hazmath.h"

#include "movetohazmath.h"
#include "eg_stokes_params.h"
#include "eg_stokes_error.h"
#include "eg_stokes_system.h"
#include "eg_stokes_precond.h"

/****************************************************************/

/****** MAIN DRIVER **************************************************************/
INT main(int argc, char* argv[])
{

  printf("\n===========================================================================\n");
  printf("Beginning Program to solve Elasticity Equation.\n");
  printf("===========================================================================\n");
  // Aug.3.2020 SLEE
  // Define variables forthe error convergence test
  int total_num_cycle = TOTAL_NUM_CYCLES_GLOBAL ;
  // SLEE initialize the vectors to save the errors for each cycle
  double L2_error_per_cycle[total_num_cycle];
  double H1_error_per_cycle[total_num_cycle];
  double H1_stress_error_per_cycle[total_num_cycle];

  double L2_EG_error_per_cycle[total_num_cycle];
  double H1_EG_error_per_cycle[total_num_cycle];
  double H1_stress_EG_error_per_cycle[total_num_cycle];
  double H1_energy_EG_error_per_cycle[total_num_cycle];


  // For Pressure
  double L2_error_p_per_cycle[total_num_cycle];
  double H1_error_p_per_cycle[total_num_cycle];

  //double L2_p_EG_error_per_cycle[total_num_cycle];
  //double H1_p_EG_error_per_cycle[total_num_cycle];


  // SLEE vector to save the DOF
  int dof_per_cycle_CG[total_num_cycle];
  int dof_per_cycle_EG[total_num_cycle];

  // For Pressure
  int dof_per_cycle_CG_p[total_num_cycle];
  int dof_per_cycle_EG_p[total_num_cycle];

  double mesh_size_per_cycle[total_num_cycle];

  // SLEE vector to save the convergence rate
  double L2_conv_rate_per_cycle[total_num_cycle];
  double H1_conv_rate_per_cycle[total_num_cycle];
  double H1_stress_conv_rate_per_cycle[total_num_cycle];


  double L2_EG_conv_rate_per_cycle[total_num_cycle];
  double H1_EG_conv_rate_per_cycle[total_num_cycle];
  double H1_stress_EG_conv_rate_per_cycle[total_num_cycle];
  double H1_energy_EG_conv_rate_per_cycle[total_num_cycle];

  // For Pressure

  double L2_p_conv_rate_per_cycle[total_num_cycle];
  double H1_p_conv_rate_per_cycle[total_num_cycle];

  //double L2_p_EG_conv_rate_per_cycle[total_num_cycle];
  //double H1_p_EG_conv_rate_per_cycle[total_num_cycle];

  int global_dim_space = 0;

  // ALL TIME STEPPING ALGORITHMS /////
  REAL time = 0.;
  REAL timestep = 0.01;
  INT timestep_number = 0;
  INT total_timestep = 10;

  for(int cycle=0; cycle<total_num_cycle; ++cycle) {
    //Aug.3.2020 SLEE initilize
    L2_error_per_cycle[cycle] = 0.;
    H1_error_per_cycle[cycle] = 0.;
    H1_stress_error_per_cycle[cycle] = 0.;

    L2_EG_error_per_cycle[cycle] = 0.;
    H1_EG_error_per_cycle[cycle] = 0.;
    H1_stress_EG_error_per_cycle[cycle] = 0.;
    H1_energy_EG_error_per_cycle[cycle] = 0.;

    // For Pressure
    L2_error_p_per_cycle[cycle] = 0.;
    H1_error_p_per_cycle[cycle] = 0.;

    dof_per_cycle_CG[cycle]=0;
    dof_per_cycle_EG[cycle]=0;

    // For Pressure
    dof_per_cycle_CG_p[cycle]=0;
    dof_per_cycle_EG_p[cycle]=0;

    mesh_size_per_cycle[cycle]=0.;

    L2_conv_rate_per_cycle[cycle]=0.;
    H1_conv_rate_per_cycle[cycle]=0.;
    H1_stress_conv_rate_per_cycle[cycle]=0.;

    L2_EG_conv_rate_per_cycle[cycle]=0.;
    H1_EG_conv_rate_per_cycle[cycle]=0.;
    H1_stress_EG_conv_rate_per_cycle[cycle]=0.;
    H1_energy_EG_conv_rate_per_cycle[cycle]=0.;

    // For Pressure
    L2_p_conv_rate_per_cycle[cycle]=0.;
    H1_p_conv_rate_per_cycle[cycle]=0.;


    printf("************ CYCLE   %d  /   %d  ************** \n", cycle, total_num_cycle);

    /****** INITIALIZE PARAMETERS **************************************************/
    // loops
    INT i;
    // Overall CPU Timing
    clock_t clk_overall_start = clock();

    // Set Parameters from Reading in Input File
    input_param inparam;
    param_input_init(&inparam);
    param_input("./input.dat", &inparam);

    // Open gridfile for reading
    printf("\nCreating mesh and FEM spaces:\n");
    //FILE* gfid = HAZ_fopen(inparam.gridfile,"r");
    //SLEE
    FILE* gfid;

    //Jul.10.2020 SLEE: setup the code to read the different mesh files for each cycle
    //gfid = HAZ_fopen(inparam.gridfile,"r");
    char filename_per_cycle[512]={'\0'};
    //sprintf(filename_per_cycle, "%s%d.haz", inparam.gridfile,cycle);

    //DEBUG SIMPLE MESH
    //    sprintf(filename_per_cycle, "%s%d.haz", inparam.gridfile,cycle+1);
    INT ncy=1+(1<<(cycle+1)); // number of nodes in one direction
    sprintf(filename_per_cycle, "%s%d.haz", inparam.gridfile,ncy);
    fprintf(stdout,"\n%s\n",filename_per_cycle);
    fflush(stdout);
    /* exit(11); */
    //continue;
    gfid = HAZ_fopen(filename_per_cycle,"r");
    //  gfid = HAZ_fopen(inparam.gridfile,"r");
    if(gfid == NULL){
      perror("Could not find and open the file !!!! ");
    }

    // Create the mesh (now we assume triangles in 2D or tetrahedra in 3D)
    // File types possible are 0 - old format; 1 - vtk format
    INT mesh_type = 0;
    clock_t clk_mesh_start = clock(); // Time mesh generation FE setup
    mesh_struct mesh;

    //printf(" --> loading grid from file: %s\n",inparam.gridfile);
    //Jul.10. 2020 SLEE
    printf(" --> loading grid from file: %s\n",filename_per_cycle);

    initialize_mesh(&mesh);   // Q1. Why only here?
    creategrid_fread(gfid,mesh_type,&mesh);
    fclose(gfid);

    INT dim = mesh.dim;
    // Jul.10.2020 SLEE : for further use in the convergence test
    global_dim_space = dim;

    // Get Quadrature Nodes for the Mesh
    INT nq1d = inparam.nquad; /* Quadrature points per dimension */
    qcoordinates *cq = get_quadrature(&mesh,nq1d);

    // Get info for and create FEM spaces
    // Order of elements: 0 - P0; 1 - P1; 2 - P2; 20 - Nedlec; 30 - Raviart-Thomas
    INT order_u = 1;
    INT order_u_eg = 0;

    // Need Spaces for each component of the Mechanics plus pressure
    fespace FE_ux; // Mechanics in x direction
    create_fespace(&FE_ux,&mesh,order_u);
    fespace FE_uy; // Mechanics in y direction
    create_fespace(&FE_uy,&mesh,order_u);
    fespace FE_uz; // Mechanics in z direction
    if(dim==3) create_fespace(&FE_uz,&mesh,order_u);
    fespace FE_u_eg; // Mechanics EG
    create_fespace(&FE_u_eg,&mesh,order_u_eg);

    INT order_p = 0;
    //INT order_p_eg = 0;
    fespace FE_p; // Pressuer
    create_fespace(&FE_p,&mesh,order_p);
    //fespace FE_p_eg; // Pressuer
    //create_fespace(&FE_p_eg,&mesh,order_p_eg);
    // Set Dirichlet Boundaries

    // TODO: I think this does the same trick
    if(BOOL_WEAKLY_IMPOSED_BC) {
      set_dirichlet_bdry(&FE_ux,&mesh,-1,-1);
      set_dirichlet_bdry(&FE_uy,&mesh,-1,-1);
      if(dim==3) set_dirichlet_bdry(&FE_uz,&mesh,-1,-1);
      set_dirichlet_bdry(&FE_u_eg,&mesh,-1,-1);
      set_dirichlet_bdry(&FE_p,&mesh,-1,-1);
    } else {
      set_dirichlet_bdry(&FE_ux,&mesh,1,1);
      set_dirichlet_bdry(&FE_uy,&mesh,1,1);
      //if(dim==3) set_dirichlet_bdry(&FE_uz,&mesh,1,1);
      set_dirichlet_bdry(&FE_u_eg,&mesh,1,1);
      set_dirichlet_bdry(&FE_p,&mesh,1,1);
    }
    //
    //
    // //set_dirichlet_bdry(&FE_p_eg,&mesh,1,1);
    //
    // if(BOOL_WEAKLY_IMPOSED_BC){
    //   for(i=0;i<FE_u_eg.ndof;i++) {
    //     FE_u_eg.dirichlet[i] = 0;
    //   }
    //   for(i=0;i<FE_ux.ndof;i++) {
    //     FE_ux.dirichlet[i] = 0;
    //   }
    //   for(i=0;i<FE_uy.ndof;i++) {
    //     FE_uy.dirichlet[i] = 0;
    //   }
    //
    //   for(i=0;i<FE_p.ndof;i++) {
    //     FE_p.dirichlet[i] = 0;
    //   }
    // }

    // Create Block System with ordering (u,p)
    INT u_ndof = FE_ux.ndof + FE_uy.ndof;
    if(dim==3) u_ndof += FE_uz.ndof;
    INT eg_ndof = FE_u_eg.ndof;
    INT p_ndof = FE_p.ndof;

    //p_debug
    INT ndof = u_ndof + eg_ndof + p_ndof;

    // Get Global FE Space
    // TODO: Is there an issue here with how many unkonwns?  Is EG really an unkonwn?
    INT nun = dim+2;
    INT nspaces = dim+2;
    block_fespace FE;
    FE.nun = nun;  // p_debug
    FE.ndof = ndof;
    FE.nbdof = FE_ux.nbdof + FE_uy.nbdof + FE_u_eg.nbdof +FE_p.nbdof; //+FE_p_eg.nbdof;
    //if(dim==3) FE.nbdof += FE_uz.nbdof;
    FE.nspaces = nspaces; // SLEE?
    FE.var_spaces = (fespace **) calloc(nspaces,sizeof(fespace *));

    /*
      FE.nun = dim+1;  // p_debug
      FE.ndof = ndof;
      FE.nbdof = FE_ux.nbdof + FE_uy.nbdof + FE_u_eg.nbdof;// +FE_p.nbdof+FE_p_eg.nbdof;
      FE.nspaces = dim+1;
      FE.var_spaces = (fespace **) calloc(dim+1,sizeof(fespace *));
    */

    FE.var_spaces[0] = &FE_ux;
    FE.var_spaces[1] = &FE_uy;
    if(dim==3) FE.var_spaces[2] = &FE_uz;
    FE.var_spaces[dim] = &FE_u_eg;
    //p_debug
    FE.var_spaces[dim+1]   = &FE_p;
    //FE.var_spaces[dim+1+1] = &FE_p_eg;

    // Set Dirichlet Boundaries
    // TODO: I don't think the if statement is necessary, if indiviual spaces BC were set correctly this will just grab those
    //if(!BOOL_WEAKLY_IMPOSED_BC)
    set_dirichlet_bdry_block(&FE,&mesh);

    clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
    printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
	   (REAL) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
    /*******************************************************************************/

    printf("***********************************************************************************\n");
    printf("Number of Elements = %d\tOrder of Quadrature = %d\n",mesh.nelm,2*nq1d-1);
    printf("\n\t--- Element Type ---\n");
    printf("Mechanics Element Type = %d\t Meachnics EG  Type = %d\n",order_u,order_u_eg);
    printf("\n\t--- Degrees of Freedom ---\n");
    printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nv,mesh.nedge,mesh.nface);
    printf("\t--> DOF: %d\n",FE.ndof);
    printf("\n\t--- Boundaries ---\n");
    printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nbv,mesh.nbedge,mesh.nbface);
    printf("\t--> Boundary DOF: %d\n",FE.nbdof);
    printf("***********************************************************************************\n\n");

    //Jul.10.2020 SLEE: insert the total number of DOF for convergnece computation
    dof_per_cycle_CG[cycle]  = u_ndof;////FE.ndof + FE.nbdof;
    dof_per_cycle_EG[cycle] =  u_ndof + eg_ndof;////FE.ndof + FE.nbdof;

    dof_per_cycle_CG_p[cycle]  = FE_p.ndof;
    dof_per_cycle_EG_p[cycle] =  FE_p.ndof;// + FE_p_eg.ndof;


    printf("FE.ndof = %d | u_ndof = %d | FE_ux.ndof = %d | FE_uy.ndof = %d |  FE_u_eg.ndof = %d \n",
	   FE.ndof , u_ndof, FE_ux.ndof  , FE_uy.ndof , FE_u_eg.ndof );
    printf("FE.ndof = %d | p_ndof = %d | FE_p.ndof = %d  \n",
	   FE.ndof , p_ndof, FE_p.ndof);


    printf("FE.nbdof = %d \n", FE.nbdof);

    printf("###########\n");
    printf("CG DOF for Mechanics = %d\n",   dof_per_cycle_CG[cycle] );
    printf("EG DOF for Mechanics = %d\n",   dof_per_cycle_EG[cycle] );
    printf("###########\n");

    printf("###########\n");
    printf("CG DOF for Pressure = %d\n",   dof_per_cycle_CG_p[cycle] );
    printf("EG DOF for Pressure = %d\n",   dof_per_cycle_EG_p[cycle] );
    printf("###########\n");


    // Get the minimum mesh size
    int zz=0;
    double min_mesh_size = 10000000.;
    double tmp_size =0;
    mesh_struct *mesh_2 = &mesh;

    for (zz=0; zz<mesh_2->nface; zz++) {

      tmp_size=mesh_2->f_area[zz];

      if(tmp_size < min_mesh_size)
	min_mesh_size = tmp_size;
    }

    mesh_size_per_cycle[cycle] = min_mesh_size;
    if(time == 0)
      {
	//Set Initial Condition
	//Inside of the assemble for p0 ...
      }

    dvector sol = dvec_create(ndof);
    dvector old_timestep_sol = dvec_create(ndof);

    //for(time = timestep; timestep_number < total_timestep; time += timestep){
    dvec_cp(&sol, &old_timestep_sol);
    timestep_number = timestep_number+1;
    printf(" ---  CYCLE = %d  ----------- TIME STEPPING TIME = %f  (timestep # = %d | %d) ------------------- \n", \
	   cycle,time,timestep_number,total_timestep);
    //printf(" [Time Step = %f]  \n", timestep);
    //fflush(stdout);

    clock_t clk_assembly_start = clock();

    // Allocate the right-hand and declare the csr matrix
    dvector b;
    // Put into Block Form
    block_dCSRmat A;
    bdcsr_alloc(nspaces,nspaces,&A);
    /*** Assemble the matrix and right hand side *******************************/
    printf("Assembling the matrix and right-hand side:\n");fflush(stdout);


    if(dim==2) assemble_global_block_neighbor(&A,&b,old_timestep_sol.val,local_assembly_Elasticity_FACE,local_assembly_Elasticity,FEM_Block_RHS_Local_Elasticity,&FE,&mesh,cq,source2D,exact_sol2D,Dexact_sol2D,exact_sol2D_dt, time,timestep);
    if(dim==3) assemble_global_block_neighbor(&A,&b,old_timestep_sol.val,local_assembly_Elasticity_FACE,local_assembly_Elasticity,FEM_Block_RHS_Local_Elasticity,&FE,&mesh,cq,source3D,exact_sol3D,Dexact_sol3D,exact_sol3D_dt, time,timestep);
    printf("\n------ Assemble done: \n");fflush(stdout);
    //printf("cycle = %d -- total = %d \n", cycle, total_num_cycle);
    /*
      if((cycle) == total_num_cycle-1)
      {

      FILE* fid;
      fid = fopen("matrix.txt","w");

      dCSRmat Amerge = bdcsr_2_dcsr(&A);
      csr_print_matlab(fid,&Amerge);
      dcsr_free(&Amerge);
      exit(0);
      }
    */
    // Eliminate boundary conditions in matrix and rhs
    if(!BOOL_WEAKLY_IMPOSED_BC) eliminate_DirichletBC_blockFE_blockA(bc2D,&FE,&mesh,&b,&A,time);
    /**************************************************/
    //  Apply Pressure "BCs" (removes singularity)
    /*
      REAL pressureval =0.;
      INT pressureloc = 0;

      clock_t clk_assembly_end = clock();
      printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL)
      (clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);

      // Prepare diagonal blocks
      dCSRmat *A_diag;
      A_diag = (dCSRmat *)calloc(dim+1, sizeof(dCSRmat));

      for(i=0;i<dim;i++){ // copy block diagonal to A_diag
      dcsr_alloc(A.blocks[i*(dim+2)]->row, A.blocks[i*(dim+2)]->col, A.blocks[i*(dim+2)]->nnz, &A_diag[i]);
      dcsr_cp(A.blocks[i*(dim+2)], &A_diag[i]);
      }

      // Get Mass Matrix for p
      dCSRmat Mp;
      assemble_global(&Mp,NULL,assemble_mass_local,&FE_p,&mesh,cq,NULL,one_coeff_scal,0.0);
      dcsr_alloc(Mp.row, Mp.col, Mp.nnz, &A_diag[dim]);
      dcsr_cp(&Mp, &A_diag[dim]);

      printf("Solving the System:\n");
      clock_t clk_solve_start = clock();

      INT solver_flag = -20;

      // Allocate solution
      dvector sol = dvec_create(ndof);
      dvector v_ux = dvec_create(FE_ux.ndof);
      dvector v_uy = dvec_create(FE_uy.ndof);
      dvector v_uz;
      if(dim==3) v_uz = dvec_create(FE_uz.ndof);
      dvector v_p = dvec_create(FE_p.ndof);

      // Set initial guess to be all zero
      dvec_set(sol.row, &sol, 0.0);
      // Set initial guess for pressure to match the known "boundary condition" for pressure
      if(dim==2) sol.val[FE_ux.ndof + FE_uy.ndof + pressureloc]  = pressureval;
      if(dim==3) sol.val[FE_ux.ndof + FE_uy.ndof + FE_uz.ndof + pressureloc]  = pressureval;

      // Set parameters for linear iterative methods
      linear_itsolver_param linear_itparam;
      param_linear_solver_init(&linear_itparam);
      param_linear_solver_set(&linear_itparam,&inparam);
      INT solver_type = linear_itparam.linear_itsolver_type;
      INT solver_printlevel = linear_itparam.linear_print_level;

      // Solve
      if(solver_type==0) { // Direct Solver
      solver_flag = block_directsolve_UMF(&A,&b,&sol,solver_printlevel);
      } else { // Iterative Solver
      if (linear_itparam.linear_precond_type == PREC_NULL) {
      solver_flag = linear_solver_bdcsr_krylov(&A, &b, &sol, &linear_itparam);
      } else {
      if(dim==2) solver_flag = linear_solver_bdcsr_krylov_block_3(&A, &b, &sol, &linear_itparam, NULL, A_diag);
      if(dim==3) solver_flag = linear_solver_bdcsr_krylov_block_4(&A, &b, &sol, &linear_itparam, NULL, A_diag);
      }
      }

      // Error Check
      if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n",solver_flag);

      clock_t clk_solve_end = clock();
      printf("Elapsed CPU Time for Solve = %f seconds.\n\n",
      (REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);

    */

    /**************************************************/
    //  Apply Pressure "BCs" (removes singularity)

    //	  REAL pressureval =0.;
    //	  INT pressureloc = 0;

    clock_t clk_assembly_end = clock();
    printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL)
	   (clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);

    printf("Solving the System:\n");fflush(stdout);
    clock_t clk_solve_start = clock();

    //*****************************************
    //    SOLVER-SOLVER AND average P	=0
    //*****************************************
    INT jj=-10,solver_flag=-20;
    REAL pmin=1e20,pmax=-1e20,sum=-1e20;
    void *numeric=NULL;
    sol.row++;
    sol.val=realloc(sol.val,sol.row*sizeof(REAL));
    dvec_set(sol.row, &sol, 0.0);
    b.row++;
    b.val=realloc(b.val,b.row*sizeof(REAL));
    b.val[b.row-1]=0e0;
    //sol.val[FE_ux.ndof + FE_uy.ndof + pressureloc]  = pressureval;
    //
    //	  add row and column for the pressure block
    // which is the pressure block?Ans: dim+1;
    // extend

    //-----------------------
    // set paramters for linear solver
    //-----------------------
    linear_itsolver_param linear_itparam;
    param_linear_solver_init(&linear_itparam);
    param_linear_solver_set(&linear_itparam,&inparam);
    INT solver_type = linear_itparam.linear_itsolver_type;
    INT solver_printlevel = linear_itparam.linear_print_level;

    // Set parameters for AMG
    AMG_param amgparam;
    param_amg_init(&amgparam);
    param_amg_set(&amgparam, &inparam);
    param_amg_print(&amgparam);

    // Prepare diagonal blocks for block preconditioner
    dCSRmat *A_diag;

    // get the blocks
    // for ux, uy, and u_eg, grab the blocks directly
    // TODO THis is done twice so I'm leaving the latter one
    // for(i=0;i<dim+1;i++){ // copy block diagonal to A_diag
    //   dcsr_alloc(A.blocks[i*(dim+3)]->row, A.blocks[i*(dim+3)]->col, A.blocks[i*(dim+3)]->nnz, &A_diag[i]);
    //   dcsr_cp(A.blocks[i*(dim+3)], &A_diag[i]);
    // }
    // // for pressure, use the mass matrix
    // dcsr_alloc(p_ndof, p_ndof, p_ndof, &A_diag[dim+1]);
    // for (i=0; i<=A_diag[dim+1].row; i++) A_diag[dim+1].IA[i] = i;
    // for (i=0; i<A_diag[dim+1].row; i++) A_diag[dim+1].JA[i] = i;
    // for (i=0; i<A_diag[dim+1].row; i++) A_diag[dim+1].val[i] = mesh.el_vol[i];

    // Linear Solver
    if(solver_type==0) { // Direct Solver
      bdcsr_extend(&A,mesh.el_vol,mesh.el_vol,(dim+1),1.0,1.0);
      printf("nblocks = %d\n", A.brow);
      getchar();
      numeric=block_factorize_UMF(&A,0);//inparam.print_level);
      // solve
      solver_flag=(INT )block_solve_UMF(&A,
					&b, // rhs.
					&sol,  // solution.
					numeric,
					0);//     inparam.print_level);
      free(numeric);
    } else { // Iterative Solver
      if (linear_itparam.linear_precond_type == PREC_NULL) {
	solver_flag = linear_solver_bdcsr_krylov(&A, &b, &sol, &linear_itparam);
      }
      else if (linear_itparam.linear_precond_type >=10 && linear_itparam.linear_precond_type <100 )
	{ // do not merge velocity unknowns, directly apply solvers
	  // prepare the preconditioner
	  A_diag = (dCSRmat *)calloc(dim+2, sizeof(dCSRmat)); // number of blocks = dim+2

	  // get the blocks
	  // for ux, uy, and u_eg, grab the blocks directly
	  for(i=0;i<dim+1;i++){ // copy block diagonal to A_diag
	    dcsr_alloc(A.blocks[i*(dim+3)]->row, A.blocks[i*(dim+3)]->col, A.blocks[i*(dim+3)]->nnz, &A_diag[i]);
	    dcsr_cp(A.blocks[i*(dim+3)], &A_diag[i]);
	  }
	  // for pressure, use the mass matrix
	  dcsr_alloc(p_ndof, p_ndof, p_ndof, &A_diag[dim+1]);
	  for (i=0; i<=A_diag[dim+1].row; i++) A_diag[dim+1].IA[i] = i;
	  for (i=0; i<A_diag[dim+1].row; i++) A_diag[dim+1].JA[i] = i;
	  for (i=0; i<A_diag[dim+1].row; i++) A_diag[dim+1].val[i] = mesh.el_vol[i];

	  // now ready to solve
	  //if(dim==2) solver_flag = linear_solver_bdcsr_krylov_block_3(&A, &b, &sol, &linear_itparam, NULL, A_diag);
	  solver_flag = linear_solver_bdcsr_krylov_block_4(&A, &b, &sol, &linear_itparam, &amgparam, A_diag);

	}
      else   // use specific solver for Stokes eg dicretization
	{

	  printf("hello\n");

	  // get preconditioner data
	  precond_block_data *precdata = get_precond_block_data_eg_stokes(&A, p_ndof, mesh.el_vol, &linear_itparam, &amgparam);

	  // setup the preconditioner
	  precond prec;
	  prec.data = precdata;
	  //prec.fct = precond_block_diag_eg_stokes_additive;
	  prec.fct = precond_block_diag_eg_stokes_multiplicative;

	  // solve
	  solver_flag = solver_bdcsr_linear_itsolver(&A, &b, &sol, &prec, &linear_itparam);

	  //TODO: clean data

	}
    }



    //// printing for debug
    /* FILE *fptmp=NULL; */
    /* fptmp=fopen("debug/u.data","w"); dvector_print(fptmp,&sol); fclose(fptmp); */
    /* fptmp=fopen("debug/b.data","w"); dvector_print(fptmp,&b); fclose(fptmp); */
    //// end printing for debug
    jj=sol.row - mesh.nelm - 1; // begin index of pressure block
    // after extension
    pmin=sol.val[jj];
    pmax=sol.val[jj];
    sum=0e0;
    for(i=0;i<mesh.nelm;i++){
      //	    fprintf(stdout,"\nnel:%7d, dof=%7d: sol=%e",	\
      //		    i,jj,sol.val[jj]);
      sum+=mesh.el_vol[i]*sol.val[jj];
      if(sol.val[jj]<pmin) pmin=sol.val[jj];
      if(sol.val[jj]>pmax) pmax=sol.val[jj];
      jj++;
    }
    fprintf(stdout,"\nINTEGRAL(p)=%13.6e, min(p)=%11.4e, max(p)=%11.4e\n",sum,pmin,pmax);
    // Error Check
    if (solver_flag < 0) fprintf(stdout,"### ERROR: Solver does not converge with error code = %d!\n", solver_flag);
    b.row--;
    b.val=realloc(b.val,b.row*sizeof(REAL));
    sol.row--;
    sol.val=realloc(sol.val,sol.row*sizeof(REAL));
    //*****************************************
    //    END-SOLVER-SOLVER AND average P=0
    //*****************************************


    /*
    // Prepare diagonal blocks
    dCSRmat *A_diag;
    A_diag = (dCSRmat *)calloc(dim+1, sizeof(dCSRmat));

    for(i=0;i<dim;i++){ // copy block diagonal to A_diag
    dcsr_alloc(A.blocks[i*(dim+2)]->row, A.blocks[i*(dim+2)]->col, A.blocks[i*(dim+2)]->nnz, &A_diag[i]);
    dcsr_cp(A.blocks[i*(dim+2)], &A_diag[i]);
    }
    // Get Mass Matrix for p
    dCSRmat Mp;
    assemble_global(&Mp,NULL,assemble_mass_local,&FE_p,&mesh,cq,NULL,one_coeff_scal,0.0);
    dcsr_alloc(Mp.row, Mp.col, Mp.nnz, &A_diag[dim]);
    dcsr_cp(&Mp, &A_diag[dim]);

    // Set parameters for linear iterative methods
    linear_itsolver_param linear_itparam;
    param_linear_solver_init(&linear_itparam);
    param_linear_solver_set(&linear_itparam,&inparam);
    INT solver_type = linear_itparam.linear_itsolver_type;
    INT solver_printlevel = linear_itparam.linear_print_level;

    // Solve
    if(solver_type==0) { // Direct Solver
    solver_flag = block_directsolve_UMF(&A,&b,&sol,solver_printlevel);
    } else { // Iterative Solver
    if (linear_itparam.linear_precond_type == PREC_NULL) {
    solver_flag = linear_solver_bdcsr_krylov(&A, &b, &sol, &linear_itparam);
    } else {
    if(dim==2) solver_flag = linear_solver_bdcsr_krylov_block_3(&A, &b, &sol, &linear_itparam, NULL, A_diag);
    if(dim==3) solver_flag = linear_solver_bdcsr_krylov_block_4(&A, &b, &sol, &linear_itparam, NULL, A_diag);
    }
    }
    */
    /////////////////////////////////////////////////////////////////
    // Error Check
    if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n",solver_flag);



    clock_t clk_solve_end = clock();
    printf("Elapsed CPU Time for Solve = %f seconds.\n\n",
	   (REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);


    //if(timestep_number == total_timestep){
    {
      printf("Compute Error at Time = %f \n", time);

      //////////////////////////////////////////////////////////////////////////////////////////////
      /********************* Compute Errors if you have exact solution ****************************/
      clock_t clk_error_start = clock();

      REAL* solerrL2 = (REAL *) calloc(nspaces, sizeof(REAL));
      REAL* solerrH1 = (REAL *) calloc(nspaces, sizeof(REAL)); // Note: No H1 error for P0 elements
      REAL* solerr_stress = (REAL *) calloc(nspaces, sizeof(REAL)); // Note: No H1 error for P0 elements

      if(dim==2) {
	L2error_block(solerrL2, sol.val, exact_sol2D, &FE, &mesh, cq, time);
	//HDerror_block(solerrH1, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
	HDsemierror_block(solerrH1, sol.val, Dexact_sol2D, &FE, &mesh, cq, time);
	//NEW SLEE Aug 17 2020
	HDsemierror_block_Stress(solerr_stress, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, time);
      }
      if(dim==3) {
	L2error_block(solerrL2, sol.val, exact_sol3D, &FE, &mesh, cq, time);
	//HDerror_block(solerrH1, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
	HDsemierror_block(solerrH1, sol.val, Dexact_sol3D, &FE, &mesh, cq, time);
	//NEW SLEE Aug 17 2020
	HDsemierror_block_Stress(solerr_stress, sol.val, exact_sol3D, Dexact_sol3D, &FE, &mesh, cq, time);
      }

      REAL uerrL2 = 0;
      REAL uerrH1 = 0;
      REAL uerr_stressH1 = 0;

      //For Pressure
      REAL uerrL2_p = 0;
      REAL uerrH1_p = 0;

      for(i=0;i<dim;i++) uerrL2 += solerrL2[i]*solerrL2[i];
      for(i=0;i<dim;i++) uerrH1 += solerrH1[i]*solerrH1[i];
      for(i=0;i<dim;i++) uerr_stressH1 += solerr_stress[i]*solerr_stress[i];

      uerrL2 = sqrt(uerrL2);
      uerrH1 = sqrt(uerrH1);
      uerr_stressH1 = sqrt(uerr_stressH1);

      //For Pressure
      //uerrL2_p += solerrL2[3]*solerrL2[3];
      //uerrH1_p += solerrL2[4]*solerrL2[4];
      //uerrL2_p = sqrt(uerrL2_p);
      //uerrH1_p = sqrt(uerrH1_p);

      uerrL2_p = solerrL2[dim+1];
      uerrH1_p = solerrH1[dim+1];


      //REAL perrL2 = solerrL2[dim];
      //REAL perrH1 = solerrH1[dim];

      printf("************* MECHANCIS   *****************************\n");
      printf("[CG] L2 Norm of u error    = %26.13e\n",uerrL2);
      printf("[CG] H1 Norm of u error    = %26.13e\n",uerrH1);
      printf("[CG] Stress Norm of u error    = %26.13e\n",uerrH1);
      printf("*******************************************************\n\n");

      //Jul. 10. 2020 SLEE save the errors for convergence computation
      L2_error_per_cycle[cycle] = uerrL2;
      H1_error_per_cycle[cycle] = uerrH1;
      H1_stress_error_per_cycle[cycle] = uerr_stressH1;

      //For Pressure
      L2_error_p_per_cycle[cycle] = uerrL2_p;
      H1_error_p_per_cycle[cycle] = uerrH1_p;

      printf("************* PRESSURE    *****************************\n");
      printf("[CG] L2 Norm of u error    = %26.13e\n",uerrL2_p);
      printf("[CG] H1 Norm of u error    = %26.13e\n",uerrH1_p);
      //printf("[CG] Stress Norm of u error    = %26.13e\n",uerrH1);
      printf("*******************************************************\n\n");


      //L2_error_p_per_cycle[cycle] = perrL2;
      //NEW SLEE Aug 23 2020
      //NEW ERROR FOR EG

      REAL* solerrL2_EG = (REAL *) calloc(nspaces, sizeof(REAL));
      REAL* solerrH1_EG = (REAL *) calloc(nspaces, sizeof(REAL)); // Note: No H1 error for P0 elements
      REAL* solerr_stress_EG = (REAL *) calloc(nspaces, sizeof(REAL)); // Note: No H1 error for P0 elements
      REAL* solerr_energy_EG = (REAL *) calloc(nspaces, sizeof(REAL)); // Note: No H1 error for P0 elements

      if(dim==2) {
	L2error_block_EG(solerrL2_EG, sol.val, exact_sol2D, &FE, &mesh, cq, time);
	HDerror_block_EG(solerrH1_EG, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, time);
	HDsemierror_block_Stress_EG(solerr_stress_EG, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, time);
	//HDsemierror_block_EnergyNorm_EG(solerr_energy_EG, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
	HDsemierror_block_EnergyNorm_EG_FaceLoop(solerr_energy_EG, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, time);
      }
      if(dim==3) {
	L2error_block_EG(solerrL2_EG, sol.val, exact_sol3D, &FE, &mesh, cq, time);
	HDerror_block_EG(solerrH1_EG, sol.val, exact_sol3D, Dexact_sol3D, &FE, &mesh, cq, time);
	HDsemierror_block_Stress_EG(solerr_stress_EG, sol.val, exact_sol3D, Dexact_sol3D, &FE, &mesh, cq, time);
	//HDsemierror_block_EnergyNorm_EG(solerr_energy_EG, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
	HDsemierror_block_EnergyNorm_EG_FaceLoop(solerr_energy_EG, sol.val, exact_sol3D, Dexact_sol3D, &FE, &mesh, cq, time);
      }

      //NEW SLEE Aug 17 2020
      //NEW ERROR FOR STRESS, \mu \epsilon(u) + \lambda \div u
      REAL uerrL2_EG = 0;
      REAL uerrH1_EG = 0;
      REAL uerr_stress_EG = 0;
      REAL uerr_energy_EG = 0;

      // For Pressure
      REAL uerrL2_EG_p = 0;
      REAL uerrH1_EG_p = 0;


      for(i=0;i<dim;i++) uerrL2_EG += solerrL2_EG[i]*solerrL2_EG[i];
      for(i=0;i<dim;i++) uerrH1_EG += solerrH1_EG[i]*solerrH1_EG[i];
      for(i=0;i<dim;i++) uerr_stress_EG += solerr_stress_EG[i]*solerr_stress_EG[i];
      //for(i=0;i<dim;i++)
      //uerr_energy_EG += solerr_energy_EG[i]*solerr_energy_EG[i];

      // TODO: Not sure I understand this calculation
      uerr_energy_EG = solerr_energy_EG[0]+solerr_energy_EG[1];
      if(dim==3) uerr_energy_EG += solerr_energy_EG[2];
      //DEBUG
      uerr_energy_EG += uerrH1_EG;

      // For Pressure
      uerrL2_EG_p = solerrL2_EG[dim+1];//*solerrL2_EG[3];
      uerrH1_EG_p = solerrH1_EG[dim+1];//*solerrL2_EG[4];

      // DEUBG L2 - H1
      uerrL2_EG = sqrt(uerrL2_EG);
      uerrH1_EG = sqrt(uerrH1_EG);
      uerr_stress_EG = sqrt(uerr_stress_EG);
      uerr_energy_EG = sqrt(uerr_energy_EG);

      // For Pressure
      uerrL2_EG_p = sqrt(uerrL2_EG_p);
      uerrH1_EG_p = sqrt(uerrH1_EG_p);

      L2_EG_error_per_cycle[cycle] = uerrL2_EG;
      H1_EG_error_per_cycle[cycle] = uerrH1_EG;
      H1_stress_EG_error_per_cycle[cycle] = uerr_stress_EG;
      H1_energy_EG_error_per_cycle[cycle] = uerr_energy_EG;

      // For Pressure
      //L2_p_EG_error_per_cycle[cycle] = uerrL2_EG_p;
      //H1_p_EG_error_per_cycle[cycle] = uerrH1_EG_p;

      printf("# of elements = %d \n", mesh.nelm);
      printf("*************     MECHANCIS   **************************\n");
      printf("L2 Norm of u (EG) error    = %26.13e\n",uerrL2_EG);
      printf("H1 Norm of u (EG) error    = %26.13e\n",uerrH1_EG);
      printf("Stress Norm of u (EG) error    = %26.13e\n",uerr_stress_EG);
      printf("Energy Norm of u (EG) error    = %26.13e\n",uerr_energy_EG);
      printf("*******************************************************\n");
      printf("*************     PRESSURE   **************************\n");
      printf("L2 Norm of u (EG) error    = %26.13e\n",uerrL2_EG_p);
      printf("H1 Norm of u (EG) error    = %26.13e\n",uerrH1_EG_p);


      /////////////////////////////////////
      /*fprintf(stdout,"\n%d,%d\n\n",A.brow,A.bcol);
	INT j;
	for(i=0;i<A.brow;i++){
	for(j=0;j<A.brow;j++){
	fprintf(stdout,"\n(%d,%d):::%d,%d,%d\n\n",i,j,
	A.blocks[j*A.bcol+i]->row,
	A.blocks[j*A.bcol+i]->col,
	A.blocks[j*A.bcol+i]->nnz);
	}
	}

	fflush(stdout);
      */
      clock_t clk_error_end = clock();
      printf("Elapsed CPU time for getting errors = %lf seconds.\n\n",(REAL)
	     (clk_error_end-clk_error_start)/CLOCKS_PER_SEC);

      /*******************************************************************************************/
      /// Plotting
      // Allocate solution
      dvector v_ux = dvec_create(FE_ux.ndof);
      dvector v_uy = dvec_create(FE_uy.ndof);
      dvector v_uz;
      if(dim==3) v_uz = dvec_create(FE_uz.ndof);
      dvector v_u_eg = dvec_create(FE_u_eg.ndof);

      dvector v_p = dvec_create(FE_p.ndof);
      //dvector v_p_eg = dvec_create(FE_p_eg.ndof);


      get_unknown_component(&v_ux,&sol,&FE,0);
      get_unknown_component(&v_uy,&sol,&FE,1);
      if(dim==3) get_unknown_component(&v_uz,&sol,&FE,2);
      get_unknown_component(&v_u_eg,&sol,&FE,dim);

      get_unknown_component(&v_p,&sol,&FE,dim+1);
      //get_unknown_component(&v_p_eg,&sol,&FE,4);




      //printf("OUTPUT?\n");


      if(inparam.print_level > 3){
	char** varname;
	char output_filename_per_cycle[512]={'\0'};
	sprintf( output_filename_per_cycle, "output/solution_%d_%d.vtu", cycle,timestep_number);
	char* soldump = output_filename_per_cycle;//"output/solution.vtu";
	//char* soldump = "output/solution.vtu";

	// TODO: Is this the right size?
	varname = malloc(5*FE.nspaces*sizeof(char *));
	varname[0] = "ux";
	varname[1] = "uy";
	if(dim==3) varname[2] = "uz";
	varname[dim] = "u_eg ";

	varname[dim+1] = "p";
	varname[dim+2] = "p_eg";

	dump_blocksol_vtk(soldump,varname,&mesh,&FE,sol.val);
	if(varname) free(varname);


	// Print in Matlab format to show vector field in nice way.
	//if(dim==3)
	//print_matlab_vector_field(&v_ux,&v_uy,&v_uz,&FE_ux);
      }

      if(solerrL2) free(solerrL2);
      if(solerrL2_EG) free(solerrL2_EG);
      if( solerr_stress ) free( solerr_stress);
      if(solerrH1) free(solerrH1);
      if(solerrH1_EG) free(solerrH1_EG);
      if(solerr_stress_EG) free(solerr_stress_EG);
      if(solerr_energy_EG) free(solerr_energy_EG);
      dvec_free( &v_ux );
      dvec_free( &v_uy );
      if(dim==3) dvec_free( &v_uz );
      dvec_free( &v_u_eg );
      dvec_free( &v_p );
      //dvec_free( &v_p_eg );


    }//if timestep == 10

    //
    bdcsr_free( &A );
    dvec_free( &b );
    // Quadrature


    //	}//Time Loop

    //RESET TIME
    time = 0.;
    timestep_number = 0;
    //timestep = timestep/2.;
    //total_timestep = total_timestep *2;
    /************ Free All the Arrays ***********************************************************/
    // CSR
    /*
      bdcsr_free( &A );
      if(solerrL2) free(solerrL2);
      if(solerrL2_EG) free(solerrL2_EG);
      if( solerr_stress ) free( solerr_stress);
      if(solerrH1) free(solerrH1);

      dvec_free( &b );
      dvec_free( &sol );
      dvec_free( &v_ux );
      dvec_free( &v_uy );
      if(dim==3) dvec_free( &v_uz );
      dvec_free( &v_u_eg );
      dvec_free( &v_p );
      dvec_free( &v_p_eg );
    */

    dvec_free( &sol );
    dvec_free( &old_timestep_sol );
    // FE Spaces
    free_fespace(&FE_ux);
    free_fespace(&FE_uy);
    if(dim==3) free_fespace(&FE_uz);
    free_fespace(&FE_u_eg);
    free_fespace(&FE_p);
    //free_fespace(&FE_p_eg);

    free_blockfespace(&FE);

    // Quadrature
    if(cq){
      free_qcoords(cq);
      free(cq);
      cq=NULL;
    }

    // Mesh
    free_mesh(&mesh);
    //*/
    // Strings
    /*******************************************************************************************/
    clock_t clk_overall_end = clock();
    printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",
	   (REAL) (clk_overall_end-clk_overall_start)/CLOCKS_PER_SEC);


  }//SLEE :: cycle loop end


  // Jul.10.2020 SLEE: Error Computation
  //SLEE compute the convergence rate and print
  INT tmp;
  for(tmp=0; tmp<total_num_cycle; ++tmp)
    {
      if(tmp == 0){
	L2_conv_rate_per_cycle[tmp] = 0;
	H1_conv_rate_per_cycle[tmp] = 0;
	H1_stress_conv_rate_per_cycle[tmp] = 0;

	L2_EG_conv_rate_per_cycle[tmp] = 0;
	H1_EG_conv_rate_per_cycle[tmp] = 0;
	H1_stress_EG_conv_rate_per_cycle[tmp] = 0;
	H1_energy_EG_conv_rate_per_cycle[tmp] = 0;

	//For Pressure
	L2_p_conv_rate_per_cycle[tmp] = 0;
	H1_p_conv_rate_per_cycle[tmp] = 0;

	//L2_p_EG_conv_rate_per_cycle[tmp] = 0;
	//H1_p_EG_conv_rate_per_cycle[tmp] = 0;

      }
      else{
	// multiplied dim since we use DOFs not h here.
	L2_conv_rate_per_cycle[tmp] = global_dim_space * (log(L2_error_per_cycle[tmp]) -log(L2_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_CG[tmp-1]) -log(dof_per_cycle_CG[tmp]) );
	H1_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_error_per_cycle[tmp]) -log(H1_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_CG[tmp-1]) -log(dof_per_cycle_CG[tmp]) );
	H1_stress_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_stress_error_per_cycle[tmp]) -log(H1_stress_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_CG[tmp-1]) -log(dof_per_cycle_CG[tmp]) );

	L2_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(L2_EG_error_per_cycle[tmp]) -log(L2_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG[tmp-1]) -log(dof_per_cycle_EG[tmp]) );
	H1_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_EG_error_per_cycle[tmp]) -log(H1_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG[tmp-1]) -log(dof_per_cycle_EG[tmp]) );

	H1_stress_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_stress_EG_error_per_cycle[tmp]) -log(H1_stress_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG[tmp-1]) -log(dof_per_cycle_EG[tmp]) );
	H1_energy_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_energy_EG_error_per_cycle[tmp]) -log(H1_energy_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG[tmp-1]) -log(dof_per_cycle_EG[tmp]) );

	//For Pressure
	L2_p_conv_rate_per_cycle[tmp] = global_dim_space * (log(L2_error_p_per_cycle[tmp]) -log(L2_error_p_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_CG_p[tmp-1]) -log(dof_per_cycle_CG_p[tmp]) );
	H1_p_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_error_p_per_cycle[tmp]) -log(H1_error_p_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_CG_p[tmp-1]) -log(dof_per_cycle_CG_p[tmp]) );

	//L2_p_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(L2_p_EG_error_per_cycle[tmp]) -log(L2_p_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG_p[tmp-1]) -log(dof_per_cycle_EG_p[tmp]) );
	//H1_p_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_p_EG_error_per_cycle[tmp]) -log(H1_p_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG_p[tmp-1]) -log(dof_per_cycle_EG_p[tmp]) );


      }


      printf("****** CYCLE = %d ****** \n", tmp);
      printf("----- MECHANCICS   ----------------------------------\n");

      printf("L2 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, L2_error_per_cycle[tmp], dof_per_cycle_CG[tmp],L2_conv_rate_per_cycle[tmp]);
      printf("H1 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_error_per_cycle[tmp], dof_per_cycle_CG[tmp],H1_conv_rate_per_cycle[tmp]);
      printf("H1 Stress Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_stress_error_per_cycle[tmp], dof_per_cycle_CG[tmp],H1_stress_conv_rate_per_cycle[tmp]);

      printf("EG - L2 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, L2_EG_error_per_cycle[tmp], dof_per_cycle_EG[tmp],
	     L2_EG_conv_rate_per_cycle[tmp]);
      printf("EG - H1 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_EG_error_per_cycle[tmp], dof_per_cycle_EG[tmp],
	     H1_EG_conv_rate_per_cycle[tmp]);

      printf("EG - Stress_Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_stress_EG_error_per_cycle[tmp], dof_per_cycle_EG[tmp],
	     H1_stress_EG_conv_rate_per_cycle[tmp]);

      printf("EG - EnergyNorm_Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_energy_EG_error_per_cycle[tmp], dof_per_cycle_EG[tmp],
	     H1_energy_EG_conv_rate_per_cycle[tmp]);

      printf("----------------------------------------------------\n");


      printf("-----  PRESSURE   ----------------------------------\n");

      printf("L2 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, L2_error_p_per_cycle[tmp], dof_per_cycle_CG_p[tmp],L2_p_conv_rate_per_cycle[tmp]);
      printf("H1 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_error_p_per_cycle[tmp], dof_per_cycle_CG_p[tmp],H1_p_conv_rate_per_cycle[tmp]);

      //printf("EG - L2 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, L2_p_EG_error_per_cycle[tmp], dof_per_cycle_EG_p[tmp],L2_p_EG_conv_rate_per_cycle[tmp]);
      //printf("EG - H1 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_p_EG_error_per_cycle[tmp], dof_per_cycle_EG_p[tmp],H1_p_EG_conv_rate_per_cycle[tmp]);

      printf("----------------------------------------------------\n");

    }

  /*
  //for Latex Print Table
  printf("** LATEX TABLE CG MECHANICS ** \n");
  for(int tmp=0; tmp<total_num_cycle; ++tmp)
  {
  printf("%d & %f & %f &  %f &  %f   &  %f &  %f \\\\ \\hline \n", dof_per_cycle_CG[tmp], L2_error_per_cycle[tmp], L2_conv_rate_per_cycle[tmp],H1_error_per_cycle[tmp], H1_conv_rate_per_cycle[tmp],H1_stress_error_per_cycle[tmp], H1_stress_conv_rate_per_cycle[tmp] );
  }

  printf("** LATEX TABLE CG PRESSURE ** \n");
  for(int tmp=0; tmp<total_num_cycle; ++tmp)
  {
  printf("%d & %f & %f &  %f &  %f   \\\\ \\hline \n", dof_per_cycle_CG[tmp], L2_error_p_per_cycle[tmp], L2_p_conv_rate_per_cycle[tmp],H1_error_p_per_cycle[tmp], H1_p_conv_rate_per_cycle[tmp]);
  }

  printf("** LATEX TABLE ALL CG PRESSURE ** \n");
  for(int tmp=0; tmp<total_num_cycle; ++tmp)
  {
  printf("%d & %f & %f &  %f &  %f   &  %f &  %f \\\\ \\hline \n", dof_per_cycle_CG[tmp], L2_error_per_cycle[tmp], L2_conv_rate_per_cycle[tmp],H1_error_per_cycle[tmp], H1_conv_rate_per_cycle[tmp],H1_error_p_per_cycle[tmp], H1_p_conv_rate_per_cycle[tmp] );
  }


  for(int tmp=0; tmp<total_num_cycle; ++tmp)
  {
  printf("%f \n",  mesh_size_per_cycle[tmp]);
  }
  */

  /*
    printf("** LATEX TABLE EG MECHANICS ** \n");
    for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
    printf("%d & %f & %f &  %f &  %f   &  %f &  %f  & %f & %f \\\\ \\hline \n",  dof_per_cycle_EG[tmp],L2_EG_error_per_cycle[tmp], L2_EG_conv_rate_per_cycle[tmp],H1_EG_error_per_cycle[tmp], H1_EG_conv_rate_per_cycle[tmp],
    H1_stress_EG_error_per_cycle[tmp], H1_stress_EG_conv_rate_per_cycle[tmp], H1_energy_EG_error_per_cycle[tmp], H1_energy_EG_conv_rate_per_cycle[tmp] );
    }

    printf("** LATEX TABLE EG PRESSURE ** \n");
    for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
    printf("%d & %f & %f &  %f &  %f  \\\\ \\hline \n", dof_per_cycle_EG[tmp], L2_p_EG_error_per_cycle[tmp], L2_p_EG_conv_rate_per_cycle[tmp],H1_p_EG_error_per_cycle[tmp], H1_p_EG_conv_rate_per_cycle[tmp]);
    }
  */
  /*
    printf("** LATEX TABLE ALL CG MECHANICS && EG PRESSURE ** \n");
    for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
    printf("%f & %d & %f & %f &  %f &  %f  & %d  & %f  & %f  \\\\ \\hline \n",mesh_size_per_cycle[tmp],  dof_per_cycle_CG[tmp],L2_error_per_cycle[tmp], L2_conv_rate_per_cycle[tmp],H1_error_per_cycle[tmp], H1_conv_rate_per_cycle[tmp],
    dof_per_cycle_EG[tmp], H1_p_EG_error_per_cycle[tmp], H1_p_EG_conv_rate_per_cycle[tmp]);
    }


    printf("** LATEX TABLE ALL EG MECHANICS && EG PRESSURE ** \n");
    for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
    printf("%f & %d & %f & %f &  %f &  %f  & %f  & %f  \\\\ \\hline \n", mesh_size_per_cycle[tmp], dof_per_cycle_EG[tmp],L2_EG_error_per_cycle[tmp], L2_EG_conv_rate_per_cycle[tmp],H1_energy_EG_error_per_cycle[tmp], H1_energy_EG_conv_rate_per_cycle[tmp],
    H1_p_EG_error_per_cycle[tmp], H1_p_EG_conv_rate_per_cycle[tmp]);
    }

    printf("** NO L2 CG ** \n");
    for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
    printf("%f & %d & %f &  %f  & %d  & %f  & %f  \\\\ \\hline \n", mesh_size_per_cycle[tmp], dof_per_cycle_CG[tmp],H1_error_per_cycle[tmp], H1_conv_rate_per_cycle[tmp],
    dof_per_cycle_EG[tmp], H1_p_EG_error_per_cycle[tmp], H1_p_EG_conv_rate_per_cycle[tmp]);
    }

    printf("** NO L2 EG ** \n");
    for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
    printf("%f & %d &  %f &  %f & %d  & %f  & %f  \\\\ \\hline \n", mesh_size_per_cycle[tmp], dof_per_cycle_EG[tmp],H1_energy_EG_error_per_cycle[tmp], H1_energy_EG_conv_rate_per_cycle[tmp],
    dof_per_cycle_EG[tmp], H1_p_EG_error_per_cycle[tmp], H1_p_EG_conv_rate_per_cycle[tmp]);
    }
  */

  return 0;


}  /* End of Program */
/*********************************************************************************************/
