/*! \file Maxwell_XJ.c
 *
 *  Created by James Adler and Xiaozhe Hu on 2/1/16.
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
 *    This is solved in a 3D rectangular region in the plane, with an interior
 *    sphere (or cube) removed.
 *
 *    Along the boundary of the region, impedence conditions are
 *    are imposed in the interior boundary and perfect conductor conditions on the
 *    exterior boundary:
 *
 *      Exterior: n * B = 0, n x E = 0,
 *      Interior: -n x (n x E) = -n x (1/mu B).
 *
 * \note This example illustrates how to combine finte-element spaces into a block system
 *       for a mixed-FEM problem.  Additionally, it shows how to deal with unconventional
 *       meshes, and weak boundary conditions.
 *
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "hazmath.h"
#include "Maxwell_data.h"
#include "Maxwell_assemble.h"
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
  trimesh mesh;
  printf(" --> loading grid from file: %s\n",inparam.gridfile);
  initialize_mesh(&mesh);
  creategrid_fread(gfid,mesh_type,&mesh);
  fclose(gfid);
  INT dim = mesh.dim;

  // Get Quadrature Nodes for the Mesh
  INT nq1d = inparam.nquad; /* Quadrature points per dimension */
  qcoordinates *cq = get_quadrature(&mesh,nq1d);

  // Time stepping parameters
  INT tsteps = inparam.time_steps;	    /* Max number of Time stepping steps */
  REAL dt = inparam.time_step_size;	    /* Size of time step */
  REAL time = 0.0;
  INT time_scheme = inparam.time_step_type; /* Indicates type of timestepping 0-CN k-BDFk */
  char time_scheme_str[30];
  if(time_scheme==0) {
    sprintf(time_scheme_str,"Crank-Nicolson");
  } else {
    sprintf(time_scheme_str,"BDF-%d",time_scheme);
  }
  INT rhs_timedep = 0;                       /* Indicates if RHS function is time-dependent */
  
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
    // Dump Mesh and FE data
    char varE[10];
    char varB[10];
    char varp[10];
    char dir[20];
    sprintf(dir,"output");
    sprintf(varE,"E");
    sprintf(varB,"B");
    sprintf(varp,"p");
    dump_fespace(&FE_E,varE,dir);
    dump_fespace(&FE_B,varB,dir);
    dump_fespace(&FE_p,varp,dir);
    FILE* fid = fopen("output/coords.dat","w");
    dump_coords(fid,mesh.cv);
    fclose(fid);
  }
  
  // Fix boundary conditions for interior obstacle and far-away boundary
  // Far Boundary (marked as 1 in mesh): E, p, B all Dirichlet
  // Obstacle Boundary (marked as -1 in mesh): p, Dirichlet, E & B Impedance
  set_dirichlet_bdry(&FE_E,&mesh,1,1);
  set_dirichlet_bdry(&FE_B,&mesh,1,1);
  set_dirichlet_bdry(&FE_p,&mesh,1,1);
  set_dirichlet_bdry(&FE_p,&mesh,-1,-1);

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
  //		<eps dE/dt,F> + <eps grad p,F> - <(1/mu)B,curl F> + <(n x E),(n x B)>_bdry = -<j,F>
  //            <dp/dt,q> - <eps E, grad q> = 0
  //            <(1/mu) dB/dt,C> + <(1/mu) curl E, C> = 0
  //
  //            p are H1 elements with 0 Dirichlet boundaries
  //            B is RT elements with div B = 0, and B*n = 0 on outer boundary and impedence condition on interior
  //            E is Nedelec elements with n x E = 0 on outer boundary and impedence condition on interior

  // We do this using only mass and incidence matrices
  //  Me = <eps E,F>   Mf = <1/mu B,C>   Mn = <p,q>
  //  De = diag(|e|)    Df = diag(|f|)
  //  G = De^(-1) ed_v    K = Df^(-1) face_ed De
  //
  
  printf("Assembling the matrix and right-hand side:\n");
  clock_t clk_assembly_start = clock();

  // Allocate the right-hand side and declare the csr matrices
  dvector b_E; b_E.val=NULL; /* RHS vector for E-equation -> <-j,F> */
  dvector b_B; b_B.val=NULL; /* RHS vector for B-equation -> 0 */
  dvector b_p; b_p.val=NULL; /* RHS vector for p-equation -> 0 */
  dvector b_E_old;           /* RHS vector for E-equation -> <-j,F> at prev time step*/
  dvector b_B_old;           /* RHS vector for B-equation -> 0 at prev time step*/
  dvector b_p_old;           /* RHS vector for p-equation -> 0 at prev time step*/
  dCSRmat Z;                 /* Impedence Boundary Matrix */
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
  assemble_global(&Me,&b_E_old,assemble_mass_local,&FE_E,&mesh,cq,current_density,permitivity,0.0);
  assemble_global(&Mf,&b_B_old,assemble_mass_local,&FE_B,&mesh,cq,zero_vec,oneovermu,0.0);
  assemble_global(&Mv,&b_p_old,assemble_mass_local,&FE_p,&mesh,cq,zero_scal,one_scal,0.0);

  // Grad and Curl matrices
  get_grad_H1toNed(&G,&mesh);
  get_curl_NedtoRT(&K,&mesh);
  
  // Build Z (actually (1+gamma)*Z for gamma NOT time-dependent)
  assemble_global_face(&Z,NULL,NULL,impedancebdry_local,NULL,&FE_E,&mesh,cq,oneplusgamm,NULL,0.0,-1,-1);

  // Start Combining
  //  A =   0        -K^T Mf          Med G
  //      Mf K         0                0
  //    -G^T Med       0                0
  //
  //  M ->    <eps E,F>
  //         <(1/mu) B,C>
  //            <pt,q>
  //
  //    Z ->  <n x E, n x F>_bdry_i
  //               0
  //               0
  //
  //    f ->     -<j,F>
  //               0
  //               0
  // M du/dt + (A+Z)u = f
  dcsr_mxm_1(&Me,&G,&MG); // Me*G
  dcsr_trans_1(&MG,&MGt); // G'*Me
  dcsr_axm(&MGt,-1); // -G'*Me
  dcsr_mxm_1(&Mf,&K,&MK); // Mf*K
  dcsr_trans_1(&MK,&MKt); // K'*Mf
  dcsr_axm(&MKt,-1); // -K'*Mf

  // Shift Matrices to put in block form
  dcsr_shift(&Z,-1);
  dcsr_shift(&MKt,-1);
  dcsr_shift(&MG,-1);
  dcsr_shift(&MK,-1);
  dcsr_shift(&MGt,-1);
  dcsr_shift(&Me,-1);
  dcsr_shift(&Mf,-1);
  dcsr_shift(&Mv,-1);

  // Block Matrix M;
  block_dCSRmat Mb;
  Mb.brow = 3; Mb.bcol = 3;
  Mb.blocks = (dCSRmat **) calloc(9,sizeof(dCSRmat *));
  Mb.blocks[0] = &Me;
  Mb.blocks[1] = NULL;
  Mb.blocks[2] = NULL;
  Mb.blocks[3] = NULL;
  Mb.blocks[4] = &Mf;
  Mb.blocks[5] = NULL;
  Mb.blocks[6] = NULL;
  Mb.blocks[7] = NULL;
  Mb.blocks[8] = &Mv;
  // Convert back to CSR and shift back
  dCSRmat M = bdcsr_2_dcsr(&Mb);
  dcsr_shift(&M,1);

  // Block Matrix AZ = A + Z (shifts needed)
  block_dCSRmat AZb;
  AZb.brow = 3; AZb.bcol = 3;
  AZb.blocks = (dCSRmat **) calloc(9,sizeof(dCSRmat *));
  AZb.blocks[0] = &Z;
  AZb.blocks[1] = &MKt;
  AZb.blocks[2] = &MG;
  AZb.blocks[3] = &MK;
  AZb.blocks[4] = NULL;
  AZb.blocks[5] = NULL;
  AZb.blocks[6] = &MGt;
  AZb.blocks[7] = NULL;
  AZb.blocks[8] = NULL;
  // Convert back to CSR and shift back
  dCSRmat AZ = bdcsr_2_dcsr(&AZb);
  dcsr_shift(&AZ,1);

  // Block RHS (need the current RHS, the updated one for time stepping, and the function evaluted at the new time step (if time-dependent))
  dvector b = dvec_create(ndof);
  dvector b_update = dvec_create(ndof);
  dvec_set(ndof,&b_update,0.0);
  for(i=0;i<mesh.nedge;i++) b.val[i] = b_E_old.val[i];
  for(i=0;i<mesh.nface;i++) b.val[i+mesh.nedge] = b_B_old.val[i];
  for(i=0;i<mesh.nv;i++) b.val[i+mesh.nedge+mesh.nv] = b_p_old.val[i];
  dvector bnew = dvec_create(ndof);

  // Global Matrix (fixed for all time (unless Z is time-dependent))
  dCSRmat Atime;
  get_timeoperator_old(&M,&AZ,time_scheme,dt,&Atime);
  dCSRmat Atime_noBC = dcsr_create(Atime.row,Atime.col,Atime.nnz);
  dcsr_cp(&Atime,&Atime_noBC);


  clock_t clk_assembly_end = clock();
  printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL)
         (clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  /**************** Perform Timestepping ********************************************************************/
  printf("===========================================================================\n");
  printf("The system is time-dependent.  Performing %d Time Step(s) of %s\n",tsteps,time_scheme_str);
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
  dvector uprev = dvec_create(ndof);
  
  blockFE_Evaluate(u.val,truesol,&FE,&mesh,0.0);
  dvec_cp(&u,&uprev);
  
  // Get components for error computing
  get_unknown_component(&uE,&u,&FE,0);
  get_unknown_component(&uB,&u,&FE,1);
  get_unknown_component(&up,&u,&FE,2);
  
  // Clean up divergence of E
  ProjectOut_Grad(&uE,&FE_p,&FE_E,&mesh,cq,&G);
  // Update B accordingly B = 1/r curl E
  dcsr_mxv_1(&K,uE.val,uB.val);
  REAL gam;
  gamm(&gam,NULL,time,NULL);
  REAL myr = 0.5*(1-sqrt(1+4/gam));
  array_ax(mesh.nface,-1.0/myr,uB.val);
  set_unknown_component(&uE,&u,&FE,0);
  set_unknown_component(&uB,&u,&FE,1);
  set_unknown_component(&uE,&uprev,&FE,0);
  set_unknown_component(&uB,&uprev,&FE,1);
  

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
  
  FILE* nid = fopen("output/norms.tex","w");
  fprintf(nid,"\\begin{table}\n");
  fprintf(nid,"\\begin{tabular}{|c|ccc|}\n");
  fprintf(nid,"\\hline\n");
  fprintf(nid,"Step&$||\\vec{E}||_{L_2(\\Omega)}$&$||\\vec{B}||_{L_2(\\Omega)}$&$||p||_{L_2(\\Omega)}$\\\\ \n");
  fprintf(nid,"\\hline\n");
  fprintf(nid,"0&%1.3f&%1.3f&%1.3f\\\\ \n",EL2,BL2,pL2);

  // ----------------------------------
  // eliminate the boundary of matrix
  // ----------------------------------
  // In this test, the matrix does not depends on time, we can do this outside the timestepping
  eliminate_DirichletBC_blockFE(NULL,&FE,&mesh,NULL,&Atime,0.0);

  // ----------------------------------
  // data for solvers
  // ----------------------------------
  // data for block csr format solver
  block_dCSRmat Atime_bcsr;
  dvector b_update_bcsr = dvec_create(ndof);
  dvector u_bcsr = dvec_create(ndof);

  // indecis for different unknows
  ivector E_idx = ivec_create(mesh.nedge);
  ivector B_idx = ivec_create(mesh.nface);
  ivector p_idx = ivec_create(mesh.nv);

  for(i=0;i<mesh.nedge;i++) E_idx.val[i] = i;
  for(i=0;i<mesh.nface;i++) B_idx.val[i] = mesh.nedge + i;
  for(i=0;i<mesh.nv;i++)    p_idx.val[i] = mesh.nedge + mesh.nface + i;

  // change to block CSR format
  // shift
  dcsr_shift(&Atime,-1);
  // allocate
  Atime_bcsr.brow = 3; Atime_bcsr.bcol = 3;
  Atime_bcsr.blocks = (dCSRmat **) calloc(9,sizeof(dCSRmat *));
  for (ii=0; ii<9 ;ii++) {
    Atime_bcsr.blocks[ii] = (dCSRmat *)calloc(1, sizeof(dCSRmat));
  }

  // assign each block
  // A_BB
  dcsr_getblk(&Atime, B_idx.val, B_idx.val, mesh.nface, mesh.nface, Atime_bcsr.blocks[0]);
  // A_BE
  dcsr_getblk(&Atime, B_idx.val, E_idx.val, mesh.nface, mesh.nedge, Atime_bcsr.blocks[1]);
  // A_Bp
  dcsr_getblk(&Atime, B_idx.val, p_idx.val, mesh.nface, mesh.nv,    Atime_bcsr.blocks[2]);
  // A_EB
  dcsr_getblk(&Atime, E_idx.val, B_idx.val, mesh.nedge, mesh.nface, Atime_bcsr.blocks[3]);
  // A_EE
  dcsr_getblk(&Atime, E_idx.val, E_idx.val, mesh.nedge, mesh.nedge, Atime_bcsr.blocks[4]);
  // A_Ep
  dcsr_getblk(&Atime, E_idx.val, p_idx.val, mesh.nedge, mesh.nv,    Atime_bcsr.blocks[5]);
  // A_pB
  dcsr_getblk(&Atime, p_idx.val, B_idx.val, mesh.nv,    mesh.nface, Atime_bcsr.blocks[6]);
  // A_pE
  dcsr_getblk(&Atime, p_idx.val, E_idx.val, mesh.nv,    mesh.nedge, Atime_bcsr.blocks[7]);
  // A_pp
  dcsr_getblk(&Atime, p_idx.val, p_idx.val, mesh.nv,    mesh.nv,    Atime_bcsr.blocks[8]);

  // shift
  dcsr_shift(&Atime,1);

  //-------------------------------
  // prepare block preconditioner
  //-------------------------------
  //-------------------------------------
  // lower block triangular for LU solve
  //-------------------------------------
  block_dCSRmat Lb;
  Lb.brow = 3; Lb.bcol = 3;
  Lb.blocks = (dCSRmat **) calloc(9,sizeof(dCSRmat *));

  // get Gt and Kt
  dCSRmat Gt;
  dcsr_trans_1(&G,&Gt);
  dCSRmat Kt;
  dcsr_trans_1(&K,&Kt);

  dCSRmat IB = dcsr_create(mesh.nface, mesh.nface, mesh.nface);
  dCSRmat IE = dcsr_create(mesh.nedge, mesh.nedge, mesh.nedge);
  dCSRmat Ip = dcsr_create(mesh.nv, mesh.nv, mesh.nv);

  for (ii=0; ii<mesh.nface; ii++) {
    IB.IA[ii] = ii;
    IB.JA[ii] = ii;
    IB.val[ii] = 1.0;
  }
  IB.IA[mesh.nface] = mesh.nface;

  for (ii=0; ii<mesh.nedge; ii++) {
    IE.IA[ii] = ii;
    IE.JA[ii] = ii;
    IE.val[ii] = 1.0;
  }
  IE.IA[mesh.nedge] = mesh.nedge;

  for (ii=0; ii<mesh.nv; ii++) {
    Ip.IA[ii] = ii;
    Ip.JA[ii] = ii;
    Ip.val[ii] = 1.0;
  }
  Ip.IA[mesh.nv] = mesh.nv;

  dcsr_shift(&Gt,-1);
  dcsr_shift(&Kt,-1);

  Lb.blocks[0] = &IE;
  Lb.blocks[1] = &Kt;
  Lb.blocks[2] = NULL;
  Lb.blocks[3] = NULL;
  Lb.blocks[4] = &IB;
  Lb.blocks[5] = NULL;
  Lb.blocks[6] = &Gt;
  Lb.blocks[7] = NULL;
  Lb.blocks[8] = &Ip;

  // Convert back to CSR and shift back
  dCSRmat L = bdcsr_2_dcsr(&Lb);
  dcsr_shift(&L,1);

  // eliminate the boundary of matrix
  eliminate_DirichletBC_blockFE(NULL,&FE,&mesh,NULL,&L,0.0);
  dcsr_shift(&L,-1);

  // get G and K without boundary
  dCSRmat Gtb;
  dCSRmat Ktb;

  dcsr_getblk(&L, E_idx.val, B_idx.val, mesh.nedge, mesh.nface, &Ktb);
  dcsr_getblk(&L, p_idx.val, E_idx.val, mesh.nv,    mesh.nedge, &Gtb);

  dCSRmat Kb;
  dcsr_trans(&Ktb,&Kb);
  dCSRmat Gb;
  dcsr_trans(&Gtb, &Gb);

  // scale offdiagonal blocks
  dcsr_axm(&Kb,  dt/2.0);
  dcsr_axm(&Gb,  dt/2.0);
  dcsr_axm(&Ktb, dt/2.0);
  dcsr_axm(&Gtb, dt/2.0);

  // clean up
  bdcsr_free(&Lb);
  dcsr_free(&L);
  dcsr_free(&Gt);
  dcsr_free(&Kt);
  dcsr_free(&IB);
  dcsr_free(&IE);
  dcsr_free(&Ip);
  //-------------------------------------

  // prepare diagonal blocks
  dCSRmat *A_diag;
  A_diag = (dCSRmat *)calloc(3, sizeof(dCSRmat));

  // first diagonal block: 2/dt * M_f
  dcsr_alloc(Mf.row, Mf.col, Mf.nnz, &A_diag[0]);
  dcsr_cp(&Mf, &A_diag[0]);
  dcsr_axm(&A_diag[0], 2.0/dt);

  dcsr_shift(&A_diag[0], 1);

  eliminate_DirichletBC(NULL,&FE_B,&mesh,NULL,&A_diag[0],0.0);

  dcsr_shift(&A_diag[0], -1);

  // second diagonal block: dt/2 * K^tMfK + 2/dt*Me + Z
  dCSRmat KT;
  dCSRmat KTMfK;
  dcsr_shift(&K,-1);
  dcsr_trans(&K, &KT);
  dcsr_rap(&KT, &Mf, &K, &KTMfK);
  dcsr_add(&KTMfK, dt/2.0, &Me, 2.0/dt, &A_diag[1]);
  dcsr_free(&KTMfK);
  dcsr_alloc(A_diag[1].row, A_diag[1].col, A_diag[1].nnz, &KTMfK);
  dcsr_cp(&A_diag[1],&KTMfK);
  dcsr_free(&A_diag[1]);
  dcsr_add(&KTMfK, 1.0, &Z, 1.0, &A_diag[1]);
  dcsr_free(&KTMfK);
  dcsr_free(&KT);
  dcsr_shift(&K,1);

  dcsr_shift(&A_diag[1], 1);

  eliminate_DirichletBC(NULL,&FE_E,&mesh,NULL,&A_diag[1],0.0);

  dcsr_shift(&A_diag[1], -1);

  // data for HX preconditioner
  dCSRmat P_curl;
  dCSRmat Grad;

  get_Pigrad_H1toNed(&P_curl,&mesh);
  get_grad_H1toNed(&Grad,&mesh);

  dcsr_shift(&P_curl, -1);  // shift
  dcsr_shift(&Grad, -1);  // shift

  // third diagonal block: dt/2 G^tMeG + 2/dt Mv
  dCSRmat GT;
  dCSRmat GTMeG;
  dcsr_shift(&G,-1);
  dcsr_trans(&G, &GT);
  dcsr_rap(&GT, &Me, &G, &GTMeG);
  dcsr_add(&GTMeG, dt/2.0, &Mv, 2.0/dt, &A_diag[2]);
  dcsr_free(&GTMeG);
  dcsr_free(&GT);
  dcsr_shift(&G,1);

  dcsr_shift(&A_diag[2], 1);

  eliminate_DirichletBC(NULL,&FE_p,&mesh,NULL,&A_diag[2],0.0);

  dcsr_shift(&A_diag[2], -1);

  clock_t clka, clkb, clk_linear_total, clk_linear_start, clk_linear_stop;
  // Begin Time Stepping
  for (i=0; i<tsteps; i++) {
    clka = clock();
    time = (i+1)*dt;
    printf("\n*******************************************\n");
    printf("TIME STEP: %d\t Actual Time = %f\n",i+1,time);
    printf("*******************************************\n\n");


    // Update RHS
    if(rhs_timedep) { // If RHS is time-dependent, get new version
      assemble_global_RHS(&b_E,&FE_E,&mesh,cq,current_density,time);
      assemble_global_RHS(&b_B,&FE_B,&mesh,cq,zero_vec,time);
      assemble_global_RHS(&b_p,&FE_p,&mesh,cq,zero_scal,time);
      for(j=0;j<mesh.nedge;j++) bnew.val[j] = b_E.val[j];
      for(j=0;j<mesh.nface;j++) bnew.val[j+mesh.nedge] = b_B.val[j];
      for(j=0;j<mesh.nv;j++) bnew.val[j+mesh.nedge+mesh.nface] = b_p.val[j];
    } else {
      dvec_cp(&b,&bnew);
    }
    fixrhs_time(&bnew, &b, &M, &AZ, &uprev, time_scheme, dt, &b_update);

    eliminate_DirichletBC_RHS_blockFE(bc,&FE,&mesh,&b_update,&Atime_noBC,time);

    
    //Let's dump matrix for comparison
    if(i==0) {
      FILE* matid = fopen("output/mat.dat","w");
      csr_print_matlab(matid,&Atime);
      fclose(matid);
    }

    // Solve
  
    //---------------------
    // solve in bcsr format
    //---------------------

    // assign right hand side
    for(ii=0;ii<mesh.nface;ii++) b_update_bcsr.val[ii] = b_update.val[ii+mesh.nedge]; // B
    for(ii=0;ii<mesh.nedge;ii++) b_update_bcsr.val[ii+mesh.nface] = b_update.val[ii]; // E
    for(ii=0;ii<mesh.nv;ii++)    b_update_bcsr.val[ii+mesh.nface+mesh.nedge] = b_update.val[ii+mesh.nface+mesh.nedge]; // p

    // assign initial guess
    for(ii=0;ii<mesh.nface;ii++) u_bcsr.val[ii] = u.val[ii+mesh.nedge]; // B
    for(ii=0;ii<mesh.nedge;ii++) u_bcsr.val[ii+mesh.nface] = u.val[ii]; // E
    for(ii=0;ii<mesh.nv;ii++)    u_bcsr.val[ii+mesh.nface+mesh.nedge] = u.val[ii+mesh.nface+mesh.nedge]; // p

    clk_linear_start = clock();

    // call solver
    solver_flag = linear_solver_bdcsr_krylov_maxwell(&Atime_bcsr, &b_update_bcsr, &u_bcsr, &linear_itparam, &amgparam, A_diag, &P_curl, &Grad, &Gb, &Kb, &Gtb, &Ktb);

    clk_linear_stop = clock();
    printf("Elapsed CPU Time for Linear solver at Time Step %d = %f seconds.\n\n\n",i+1,(REAL) (clk_linear_stop-clk_linear_start)/CLOCKS_PER_SEC);

    clk_linear_total = clk_linear_total + (clk_linear_stop-clk_linear_start);

    // assign solution back
    for(ii=0;ii<mesh.nedge;ii++) u.val[ii] = u_bcsr.val[mesh.nface + ii];
    for(ii=0;ii<mesh.nface;ii++) u.val[ii+mesh.nedge] = u_bcsr.val[ii];
    for(ii=0;ii<mesh.nv;ii++)    u.val[ii+mesh.nface+mesh.nedge] = u_bcsr.val[ii+mesh.nface+mesh.nedge];

    
    // Update time steps
    dvec_cp(&u,&uprev);
    dvec_cp(&bnew,&b);
    get_unknown_component(&uE,&u,&FE,0);
    get_unknown_component(&uB,&u,&FE,1);
    get_unknown_component(&up,&u,&FE,2);

    // Compute Errors
    printf("\n------------- Norms and Errors -----------------------------\n");
    Eerr = L2error(uE.val,Etrue,&FE_E,&mesh,cq,time);
    Berr = L2error(uB.val,Btrue,&FE_B,&mesh,cq,time);
    perr = L2error(up.val,ptrue,&FE_p,&mesh,cq,time);
    EL2 = L2norm(uE.val,&FE_E,&mesh,cq);
    BL2 = L2norm(uB.val,&FE_B,&mesh,cq);
    pL2 = L2norm(up.val,&FE_p,&mesh,cq);
    printf("||E||_0 = %25.16e\n",EL2);
    printf("||B||_0 = %25.16e\n",BL2);
    printf("||p||_0 = %25.16e\n",pL2);
    printf("||E-Etrue||_0 = %25.16e\n",Eerr);
    printf("||B-Btrue||_0 = %25.16e\n",Berr);
    printf("||p-ptrue||_0 = %25.16e\n",perr);
    printf("------------------------------------------------------------\n");
    fprintf(nid,"%d&%1.3f&%1.3f&%1.3f\\\\ \n",i+1,EL2,BL2,pL2);

    clkb = clock();
    printf("Elapsed CPU Time for Time Step %d = %f seconds.\n\n\n",i+1,(REAL) (clkb-clka)/CLOCKS_PER_SEC);
  } // End of Timestepping Loop
  fprintf(nid,"\\hline\n");
  fprintf(nid,"\\end{tabular}\n");
  fprintf(nid,"\\caption{Sphere Obstacle. $\\gamma = %f$.  $h = 1/8$.}\n",0.05);
  fprintf(nid,"\\end{table}");
  fclose(nid);
  
  clock_t clk2 = clock();

  printf("Elapsed CPU Time for Assembling = %f seconds.\n\n",(REAL) (clk2-clk1 - clk_linear_total)/CLOCKS_PER_SEC);

  printf("Elapsed CPU Time for Linear solvers  = %f seconds.\n\n",(REAL) clk_linear_total/CLOCKS_PER_SEC);

  printf("Elapsed CPU Time for Solving All Time Steps = %f seconds.\n\n",(REAL) (clk2-clk1)/CLOCKS_PER_SEC);
  /*******************************************************************************************/


  /******** Free All the Arrays *************************************************************/
  // CSR Matrices
  dcsr_free(&Atime);
  dcsr_free(&Atime_noBC);
  dcsr_free(&Z);
  dcsr_free(&K);
  dcsr_free(&G);
  dcsr_free(&Me);
  dcsr_free(&Mf);
  dcsr_free(&Mv);
  dcsr_free(&MG);
  dcsr_free(&MGt);
  dcsr_free(&MK);
  dcsr_free(&MKt);
  dcsr_free(&M);
  dcsr_free(&AZ);
  
  dcsr_free(&A_diag[0]);
  dcsr_free(&A_diag[1]);
  dcsr_free(&A_diag[2]);
  
  if (A_diag) free(A_diag);
  
  dcsr_free(&P_curl);
  dcsr_free(&Grad);

  dcsr_free(&Gtb);
  dcsr_free(&Ktb);
  dcsr_free(&Kb);
  dcsr_free(&Gb);
  
  // block CSR Matrices
  bdcsr_free(&Mb);
  bdcsr_free(&AZb);
  
  bdcsr_free(&Atime_bcsr);
  
  
  // Vectors
  if(b_E.val) free(b_E.val);
  if(b_B.val) free(b_B.val);
  if(b_p.val) free(b_p.val);
  if(b_E_old.val) free(b_E_old.val);
  if(b_B_old.val) free(b_B_old.val);
  if(b_p_old.val) free(b_p_old.val);
  if(b.val) free(b.val);
  if(b_update.val) free(b_update.val);
  if(bnew.val) free(bnew.val);
  if(u.val) free(u.val);
  if(uprev.val) free(uprev.val);
  if(uE.val) free(uE.val);
  if(uB.val) free(uB.val);
  if(up.val) free(up.val);

  if (B_idx.val) free(B_idx.val);
  if (E_idx.val) free(E_idx.val);
  if (p_idx.val) free(p_idx.val);

  if (u_bcsr.val) free(u_bcsr.val);
  if (b_update_bcsr.val) free(b_update_bcsr.val);
  
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
  
  
  /*****************************************************************************************/

  clock_t clk_overall_end = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",(REAL) (clk_overall_end-clk_overall_start)/CLOCKS_PER_SEC);
  return 0;

}	/* End of Program */
/*******************************************************************************************/
