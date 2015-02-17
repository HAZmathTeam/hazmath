/*
 *  ReactionAdvectionDiffusion.c
 *  
 *  Created by James Adler and Xiaozhe Hu on 1/9/15.
 *  Copyright 2015_HAZMAT__. All rights reserved.
 *
 *  Discussion:
 *
 *    This program solves the Advection-Reaction-Diffusion Equation using finite elements
 *
 *      -Div(a(x)grad(u)) + b(x)*grad(u) + c(x)u = 0
 *
 *    in a 2D or 3D
 *
 *   Along the boundary of the region, Dirichlet conditions are imposed:
 *
 *      u = 0
 */

/*********** EXTERNAL FUNCTIONS *****************/
// Standard Includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <getopt.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
// Internal Includes
#include "sparse.h"
#include "grid.h"
#include "vec.h"
#include "quad.h"
#include "fem.h"
/************************************************/

/****** MAIN DRIVER *****************************/
int main (int argc, char* argv[]) 
{	
  printf("\nBeginning Program to create Finite Element Matrices for ReactionAdvectionDiffusion.c\n");	
  /****** INITIALIZE PARAMETERS ****************/
  // Loop Indices
  INT i,j;

  // Timing Parameters
  clock_t clk_start,clk_end,clk1,clk2;
  clk_start = clock();
	
  // Grab parameters from input.dat file
  INT ipar[55];
  REAL fpar[20];	
  char gridfile[50];
  getinput(gridfile,fpar,ipar);
  // Open gridfile for reading
  gfid = fopen(gridfile,"r");
  if( gfid == NULL ) { 
    printf("\nError opening Grid File!!!\nFile (%s) probably doesn't exist...\nAborting Program.\n\n",gridfile);
    return 0;
  }
	
  // Dimension is needed for all this to work
  INT mydim = ipar[0];

  // Create the mesh (now we assume triangles in 2D or tetrahedra in 3D)
  clk1 = clock();
  trimesh mesh;
  printf("\nLoading grid from file and creating mesh and all its properties: %s->\n",gridfile);
  initialize_mesh(&mesh);
  creategrid(gfid,mydim,mesh);
  fclose(gfid);

  // Get Quadrature Nodes for the Mesh
  INT nq1d = ipar[1];	/* Quadrature points per dimension */
  qcoordinates cq = get_quadrature(&mesh,nq1d);
  	
  // Get info for and create FEM spaces
  INT poly = ipar[2];	/* Order of Elements: 0 - P0; 1 - P1; 2 - P2; -1 - Nedelec; -2 - Raviart-Thomas */
  fespace FE;
  create_fespace(&FE,&mesh,poly);

  if(dumpmesh==1) {
    dumpmeshdata(el_n.IA,el_n.JA,el_v.IA,el_v.JA,cn.x,cn.y,cn.z,nelm,n,nvert,element_order,mydim+1,mydim,40,n_bdry,n);
  }
	
  clk1 = clock();
  printf("Elapsed CPU Time for Conversion = %f seconds.\n\n",(REAL) (clk1 - clk0)/CLOCKS_PER_SEC);
  /*******************************************************************************************/
	
  printf("***********************************************************************************************\n");
  printf("Number of Elements = %d\n",nelm);
  printf("Number of Nodes = %d\t\tNumber of Edges = %d\n",n,nedge);
  printf("Number of Boundary Nodes = %d\tNumber of Boundary Edges = %d\n",nbvert,nbedge);
  printf("Order of Elements = %d\t\tOrder of Quadrature = %d\n",poly,2*nq1d-1);
  printf("************************************************************************************************\n\n");
	
  /*** Assemble the matrix and right hand side *******************************/
  printf("Assembling the Matrix and Right-Hand Side:\n");
  clk0 = clock();
	
  // Number of degrees of freedom equals number of nodes plus edges for P2 elements for each component of u plus the number of elements for P0 elements for p
  //ndof = mydim*n+nelm;
  ndof = n;
  A.row = ndof;
  A.col = ndof;
  f = calloc(ndof,sizeof(REAL));
	
  // First get P2 Laplacian Matrix for each component of u
  //iA = calloc(ndof+1,sizeof(INT));
  A.IA = calloc(ndof+1,sizeof(INT));

  // Get non-zeros of A (ignores cancellations, so maybe more than necessary)
  // stiffG_nnz(iA,&nnzA,el_n.IA,el_n.JA,n,nelm,n_bdry);
  stiffG_nnz(A.IA,&A.nnz,el_n.IA,el_n.JA,n,nelm,n_bdry);

	
  // Build Stiffness Matrix and right-hand side vector
  A.val = calloc(A.nnz,sizeof(REAL));
  for (i=0; i<A.nnz; i++) {
    A.val[i] = 0;
  }
  A.JA = calloc(A.nnz,sizeof(INT));

  H1_assemble_Lap(A.IA,A.JA,A.val,f,cn.x,cn.y,cn.z,cq.x,cq.y,cq.z,cq.w,el_n.IA,el_n.JA,element_order,nq1d,coef,mydim,nelm,n,n_bdry,compRHS,compB,0,NULL,NULL,NULL,0.0);
  /* INT* iM = calloc(ndof+1,sizeof(INT)); */
  /* stiffG_nnz(iM,&nnzA,el_n.IA,el_n.JA,n,nelm,n_bdry); */
  /* INT* jM = calloc(nnzA,sizeof(INT)); */
  /* REAL* M = calloc(nnzA,sizeof(REAL)); */
  /* for(i=0;i<nnzA;i++) { M[i] = 0.0; } */
  /* H1_assemble_Mass(iM,jM,M,cn.x,cn.y,cn.z,cq.x,cq.y,cq.z,cq.w,el_n.IA,el_n.JA,element_order,nq1d,mydim,nelm,n,n_bdry,compB,1.0); */

  clk1 = clock();
  printf("Elapsed CPU Time for Assembly = %f seconds.\n\n",(REAL) (clk1-clk0)/CLOCKS_PER_SEC);
  /*******************************************************************************************/
	
  /**************** Solve ********************************************************************/
  printf("Solving the System: Using Preconditioned Minimal Residual Method (CG for all inversions)\n");
  clk0 = clock();
  u = calloc(ndof,sizeof(REAL));
	
  // Initialize solution to be rhs / diag(A)
  for (i=0; i<A.row; i++) {
    cola = A.IA[i]-1;
    colb = A.IA[i+1]-1;
    for (jaa=cola; jaa<colb; jaa++) {
      j = A.JA[jaa]-1;
      if (j==i) {
	diag = A.val[jaa];
      }
    }
    u[i] = f[i]/diag;
  }

  REAL* res = (REAL *) calloc(ndof,sizeof(REAL));
  if (solvertype==0) {
    cg(A.IA,A.JA,A.val,f,ndof,u,pcgtol,maxit,outprint);
  } else if (solvertype==1) {
    minres(A.IA,A.JA,A.val,f,ndof,u,pcgtol,maxit,outprint);
  } else if (solvertype==2) {
    gmres(A.IA,A.JA,A.val,f,ndof,u,pcgtol,maxit,outprint,ndof);
  } else if (solvertype==6) {
    // Compute Residuals
    bminax(f,A.IA,A.JA,A.val,u,&ndof,res);

    // Loop over subdomains, restrict residuals, assemble restricted matrices A and B and solve constrained system
    //INT macroelms1d = sqrt(nelm/2);
    INT nnew,k;
    // Choose size of subdomains (few big elements subsize = 2 or lots of small elements subsize = log_2{macroelms1d})
    dCSRmat Asub;
    Asub.IA = (INT *) calloc(ndof+1,sizeof(INT));
    Asub.JA = (INT *) calloc(A.nnz,sizeof(INT));
    Asub.val = (REAL *) calloc(A.nnz,sizeof(REAL));
    Asub.row = ndof;
    Asub.col = ndof;
    Asub.nnz = A.nnz;
    INT subsize = 2;
    INT subdomains = ((INT) pow(2,subsize))-1;
    //subdomains = 2;
    REAL overlap = 1.0/(subdomains+1);
    // overlap = 0.5;
    REAL domainsize = 2.0/(subdomains+1);
    //  domainsize=0.5;
    REAL* fAsub = (REAL *) calloc(ndof,sizeof(REAL));
      
    REAL* unew;
    REAL* resnew;
    INT* newnodes = (INT *) calloc(ndof,sizeof(INT));
    

    // Find the centers of the elements if needed
    resnew = (REAL *) calloc(ndof,sizeof(REAL));
    unew = (REAL *) calloc(ndof,sizeof(REAL));
    for(i=0;i<ndof;i++) {
      resnew[i] = res[i];
      unew[i] = 0.0;
    }

    INT s;
    for(s=0;s<128;s++) {
      for(i=0;i<subdomains;i++) {
	for(j=0;j<subdomains;j++) {
	  // Restrict Atmp,B,ftmp,fB
	  restrictmat2(&Asub,fAsub,A,res,cn,n,subdomains,overlap,domainsize,j,i,&nnew,newnodes);

	  for(k=0;k<nnew;k++) {
	    resnew[k] = fAsub[k];
	    unew[k] = 0.0;
	  }
	 
	  // Solve constrained system via schur_complement GMRES
	  cg(Asub.IA,Asub.JA,Asub.val,fAsub,nnew,unew,pcgtol,maxit,1);
	  /* maxa = (INT *) calloc(ndofA+1,sizeof(INT)); */
	  /* upperi_(Asub.IA,Asub.JA,maxa,&ndofAnew,&nwk); */
	  /* U = (REAL *) calloc(nwk,sizeof(REAL)); */
	  /* upper_(Asub.IA,Asub.JA,Asub.val,maxa,U,&ndofAnew,&nwk); */
	  /* dolu_(U,maxa,&ndofAnew); */
	  /* gausselim(Asub.IA,Asub.JA,Asub.val,U,maxa,fAsub,ndofAnew,unew); */
	  /* free(U); */
	  /* free(maxa); */
	  // minres(Asub.IA,Asub.JA,Asub.val,fAsub,nnew,unew,pcgtol,maxit,0);
	  // schur_solve(Asub.IA,Asub.JA,Asub.val,Bsub.IA,Bsub.JA,Bsub.val,resnew,nun*nnew,elmnew,unew,solvetol,maxit,0,2);
	  
	  // Update utmp
	  for(k=0;k<nnew;k++) {
	    u[newnodes[k]-1]+=unew[k];
	  }

	  /* if(i*subdomains+j==0) exit(0); */
	  // Recompute Residuals
	  bminax(f,A.IA,A.JA,A.val,u,&ndof,res);
	} 
      }
    }
  }
  clk1 = clock();
  printf("Elapsed CPU Time for Solve = %f seconds.\n\n",(REAL) (clk1-clk0)/CLOCKS_PER_SEC);
  /*******************************************************************************************/
	
  /**************** Get True Solution and Error***********************************************/
  utrue = calloc(ndof,sizeof(REAL));
  myerr = calloc(ndof,sizeof(REAL));
  if(havetrue) {
    printf("Computing True Solution and Errors (if available):\n");
    clk0 = clock();
		
    REAL udott = -6.66;
    if(mydim==3) {
      for (i=0; i<n; i++) {
	getknownfunction(&udott,cn.x[i],cn.y[i],cn.z[i],0.0,mydim,1,0.0,ut);
	utrue[i]=udott;
	myerr[i] = fabs(u[i]-utrue[i]);
	myerr[i] = u[i]-utrue[i];
      }
    } else if (mydim==2) {
      for (i=0; i<n; i++) {
	getknownfunction(&udott,cn.x[i],cn.y[i],0.0,0.0,mydim,1,coef,ut);
	utrue[i]=udott;
	myerr[i] = fabs(u[i]-utrue[i]);
	myerr[i] = u[i]-utrue[i];
      }
    } else {
      printf("bad dimension\n");
      exit(0);
    }
		
    xdy(&ndof,myerr,myerr,&errnorml2);
    L2normMassassemble(&errnorm,myerr,el_n.IA,el_n.JA,el_n.IA,el_n.JA,cn.x,cn.y,cn.z,cq.x,cq.y,cq.z,cq.w,nelm,nq1d,mydim,element_order,element_order,NULL,NULL,NULL,NULL,0);
    clk1 = clock();
    printf("Elapsed CPU Time for Getting Errors = %f seconds.\n\n",(REAL) (clk1-clk0)/CLOCKS_PER_SEC);
  }
  /*******************************************************************************************/
	
  /**************** Print Results or Dump Results *******************************************/
  if (dumpsol>=2) { 
    uid = fopen("output/sol.dat","w"); 
    if(havetrue) { truid = fopen("output/true.dat","w"); }
  }
  dumpsolution(u,utrue,uid,truid,ndof,n,nelm,dumpsol,1,mydim,havetrue);
  if(dumpmat>=2) {
    matid = fopen("output/mat.dat","w");
    rhsid = fopen("output/rhs.dat","w");
  }
  dumpmatrices(A.IA,A.JA,A.val,f,matid,rhsid,A.row,dumpmat);
  if (havetrue) {
    printf("***********************************************************************************************\n");
    printf("L2 Norm of u error = %25.17g\n",errnorm);
    printf("l2 Norm of error =   %25.17g\n",sqrt(errnorml2));
    printf("***********************************************************************************************\n");
    if(dumpsol>=2) {fclose(truid);}
  }
  if (dumpsol>=2) {
    fclose(uid);
  }
  if (dumpmat>=2) {
    fclose(rhsid);
    fclose(matid);
  }
  /*******************************************************************************************/
	
  /******** Free All the Arrays ***************************************************************/
  free_mesh(mesh);
  free_qcoords(cq);
  free_fespace(FE);
  /*******************************************************************************************/
	
  clkb = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",(REAL) (clkb-clka)/CLOCKS_PER_SEC);
  return 0;	
	
}	/* End of Program */
/*******************************************************************************************/


