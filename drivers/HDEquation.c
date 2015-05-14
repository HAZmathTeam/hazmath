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

/*********** EXTERNAL FUNCTIONS *********************************************************/
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
#include "macro.h"
#include "grid.h"
#include "sparse.h"
#include "vec.h"
#include "functs.h"
#include "fem.h"
/****************************************************************************************/

/******** Data Input ********************************************************************/
// Right-hand Side
void myrhs(REAL *val,REAL* x,REAL time) {
    *val = 2*M_PI*M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);
}

// Boundary Conditions
void bc(REAL *val,REAL* x,REAL time) {
    *val = 0.0;
}

// PDE Coefficients
void diffcoeff(REAL *val,REAL* x,REAL time) {
    *val = 1.0;
}

// True Solution (if you have one)
void truesol(REAL *val,REAL* x,REAL time) {
    *val = sin(M_PI*x[0])*sin(M_PI*x[1]);
}


/****** MAIN DRIVER *********************************************************************/
int main (int argc, char* argv[])
{
    printf("\nBeginning Program to solve Reaction-Advection-Diffusion Problem.\n");
    /****** INITIALIZE PARAMETERS ********************************************************/
    // Timing Parameters
    clock_t clk_start,clk_end,clk1,clk2;
    clk_start = clock();
    
    // Grab parameters from input.dat file
    INT ipar[55];
    REAL fpar[20];
    char gridfile[50];
    getinput(gridfile,fpar,ipar);
    // Open gridfile for reading
    FILE* gfid = fopen(gridfile,"r");
    if( gfid == NULL ) {
        printf("\nError opening Grid File!!!\n");
        printf("File (%s) probably doesn't exist!\n\n",gridfile);
        return 0;
    }
    
    // Dimension is needed for all this to work
    INT dim = ipar[0];
    
    // Create the mesh (now we assume triangles in 2D or tetrahedra in 3D)
    clk1 = clock();
    trimesh mesh;
    printf("\nLoading grid from file: %s\n",gridfile);
    printf("--> creating mesh and all its properties.\n");
    initialize_mesh(&mesh);
    creategrid(gfid,dim,0,&mesh);
    fclose(gfid);
    
    // Get Quadrature Nodes for the Mesh
    INT nq1d = ipar[1];	/* Quadrature points per dimension */
    qcoordinates *cq = get_quadrature(&mesh,nq1d);
    
    // Get info for and create FEM spaces
    // Order of Elements: 0 - P0; 1 - P1; 2 - P2; -1 - Nedelec; -2 - Raviart-Thomas
    INT order = ipar[2];
    fespace FE;
    initialize_fespace(&FE);
    create_fespace(&FE,&mesh,order);
    
    // Dump some data if needed
    INT dumpmesh=ipar[37];
    if(dumpmesh==1) {
        dump_fespace(&FE);
    }
    
    clk2 = clock();
    printf("Elapsed CPU Time for Mesh and FEM Space Construction = %f seconds.\n\n",
           (REAL) (clk2 - clk1)/CLOCKS_PER_SEC);
    /**************************************************************************************/
    
    printf("***************************************************************************\n");
    printf("Number of Elements = %d\tOrder of Elements = %d\tOrder of Quadrature = %d\n",mesh.nelm,FE.FEtype,2*nq1d-1);
    printf("\n--- Degrees of Freedom ---\n");
    printf("Vertices: %d\tEdges: %d\tFaces: %d",mesh.nv,mesh.nedge,mesh.nface);
    printf("\t--> Total DOF: %d\n",FE.ndof);
    printf("\n--- Boundaries ---\n");
    printf("Vertices: %d\tEdges: %d\tFaces: %d",mesh.nbv,mesh.nbedge,mesh.nbface);
    printf("\t--> Total Boundary DOF: %d\n",FE.nbdof);
    printf("***************************************************************************\n\n");
    
    /*** Assemble the matrix and right hand side *******************************/
    printf("Assembling the Matrix and Right-Hand Side:\n");
    clk1 = clock();
    
    // Allocate the right-hand side and declare the csr matrix
    dvector b;
    dCSRmat A;
    
    // Assemble the matrix
    assemble_global(&A,&b,assemble_mass_local,&FE,&mesh,cq,myrhs,bc,diffcoeff,0.0);
    FILE* matid = fopen("mat.dat","w");
    csr_print_matlab(matid,&A);
    fclose(matid);
    FILE* rhsid = fopen("rhs.dat","w");
    dvector_print(rhsid,&b);
    fclose(rhsid);
    
    clk2 = clock();
    printf("Elapsed CPU Time for Assembly = %f seconds.\n\n",(REAL) (clk2-clk1)/CLOCKS_PER_SEC);
    /*******************************************************************************************/
    
    
    /**************** Solve ********************************************************************/
    printf("Solving the System: Using Krylov Subspace Methos \n");
    
    // parameters
    INT solver_flag;
    REAL tol = 1e-6;
    INT MaxIt = 100;
    SHORT restart = 50;
    SHORT stop_type = STOP_REL_RES;
    SHORT print_level = PRINT_MORE;
    
    // Allocate the solution
    dvector u = dvec_create(b.row);
    
    // set initial guess to be all zero
    dvec_set(u.row, &u, 0.0);

    // shift A
    dcsr_shift(&A, -1);
    
    // solve the linear system
    
    printf("\n");
    
    printf("Conjugate gradient method:\n");
    solver_flag = dcsr_pcg(&A, &b, &u, NULL, tol, MaxIt, stop_type, print_level);
    
    printf("\n");
    
    /* printf("GMRes method:\n"); */
    /* solver_flag = dcsr_pvgmres(&A, &b, &u, NULL, tol, MaxIt, restart, stop_type, print_level); */
    
    if (solver_flag < 0) printf("### ERROR: CG does not converge with error code = %d!\n", solver_flag);
    
    // shift A back
    dcsr_shift(&A, 1);
    
    /*******************************************************************************************/
    
    
    /* /\**************** Solve ********************************************************************\/ */
    /* printf("Solving the System: Using Preconditioned Minimal Residual Method (CG for all inversions)\n"); */
    /* clk0 = clock(); */
    /* u = calloc(ndof,sizeof(REAL)); */
    
    /* // Initialize solution to be rhs / diag(A) */
    /* for (i=0; i<A.row; i++) { */
    /*   cola = A.IA[i]-1; */
    /*   colb = A.IA[i+1]-1; */
    /*   for (jaa=cola; jaa<colb; jaa++) { */
    /*     j = A.JA[jaa]-1; */
    /*     if (j==i) { */
    /* 	diag = A.val[jaa]; */
    /*     } */
    /*   } */
    /*   u[i] = f[i]/diag; */
    /* } */
    
    /* REAL* res = (REAL *) calloc(ndof,sizeof(REAL)); */
    /* if (solvertype==0) { */
    /*   cg(A.IA,A.JA,A.val,f,ndof,u,pcgtol,maxit,outprint); */
    /* } else if (solvertype==1) { */
    /*   minres(A.IA,A.JA,A.val,f,ndof,u,pcgtol,maxit,outprint); */
    /* } else if (solvertype==2) { */
    /*   gmres(A.IA,A.JA,A.val,f,ndof,u,pcgtol,maxit,outprint,ndof); */
    /* } else if (solvertype==6) { */
    /*   // Compute Residuals */
    /*   bminax(f,A.IA,A.JA,A.val,u,&ndof,res); */
    
    /*   // Loop over subdomains, restrict residuals, assemble restricted matrices A and B and solve constrained system */
    /*   //INT macroelms1d = sqrt(nelm/2); */
    /*   INT nnew,k; */
    /*   // Choose size of subdomains (few big elements subsize = 2 or lots of small elements subsize = log_2{macroelms1d}) */
    /*   dCSRmat Asub; */
    /*   Asub.IA = (INT *) calloc(ndof+1,sizeof(INT)); */
    /*   Asub.JA = (INT *) calloc(A.nnz,sizeof(INT)); */
    /*   Asub.val = (REAL *) calloc(A.nnz,sizeof(REAL)); */
    /*   Asub.row = ndof; */
    /*   Asub.col = ndof; */
    /*   Asub.nnz = A.nnz; */
    /*   INT subsize = 2; */
    /*   INT subdomains = ((INT) pow(2,subsize))-1; */
    /*   //subdomains = 2; */
    /*   REAL overlap = 1.0/(subdomains+1); */
    /*   // overlap = 0.5; */
    /*   REAL domainsize = 2.0/(subdomains+1); */
    /*   //  domainsize=0.5; */
    /*   REAL* fAsub = (REAL *) calloc(ndof,sizeof(REAL)); */
    
    /*   REAL* unew; */
    /*   REAL* resnew; */
    /*   INT* newnodes = (INT *) calloc(ndof,sizeof(INT)); */
    
    
    /*   // Find the centers of the elements if needed */
    /*   resnew = (REAL *) calloc(ndof,sizeof(REAL)); */
    /*   unew = (REAL *) calloc(ndof,sizeof(REAL)); */
    /*   for(i=0;i<ndof;i++) { */
    /*     resnew[i] = res[i]; */
    /*     unew[i] = 0.0; */
    /*   } */
    
    /*   INT s; */
    /*   for(s=0;s<128;s++) { */
    /*     for(i=0;i<subdomains;i++) { */
    /* 	for(j=0;j<subdomains;j++) { */
    /* 	  // Restrict Atmp,B,ftmp,fB */
    /* 	  restrictmat2(&Asub,fAsub,A,res,cn,n,subdomains,overlap,domainsize,j,i,&nnew,newnodes); */
    
    /* 	  for(k=0;k<nnew;k++) { */
    /* 	    resnew[k] = fAsub[k]; */
    /* 	    unew[k] = 0.0; */
    /* 	  } */
    
    /* 	  // Solve constrained system via schur_complement GMRES */
    /* 	  cg(Asub.IA,Asub.JA,Asub.val,fAsub,nnew,unew,pcgtol,maxit,1); */
    /* 	  /\* maxa = (INT *) calloc(ndofA+1,sizeof(INT)); *\/ */
    /* 	  /\* upperi_(Asub.IA,Asub.JA,maxa,&ndofAnew,&nwk); *\/ */
    /* 	  /\* U = (REAL *) calloc(nwk,sizeof(REAL)); *\/ */
    /* 	  /\* upper_(Asub.IA,Asub.JA,Asub.val,maxa,U,&ndofAnew,&nwk); *\/ */
    /* 	  /\* dolu_(U,maxa,&ndofAnew); *\/ */
    /* 	  /\* gausselim(Asub.IA,Asub.JA,Asub.val,U,maxa,fAsub,ndofAnew,unew); *\/ */
    /* 	  /\* free(U); *\/ */
    /* 	  /\* free(maxa); *\/ */
    /* 	  // minres(Asub.IA,Asub.JA,Asub.val,fAsub,nnew,unew,pcgtol,maxit,0); */
    /* 	  // schur_solve(Asub.IA,Asub.JA,Asub.val,Bsub.IA,Bsub.JA,Bsub.val,resnew,nun*nnew,elmnew,unew,solvetol,maxit,0,2); */
    
    /* 	  // Update utmp */
    /* 	  for(k=0;k<nnew;k++) { */
    /* 	    u[newnodes[k]-1]+=unew[k]; */
    /* 	  } */
    
    /* 	  /\* if(i*subdomains+j==0) exit(0); *\/ */
    /* 	  // Recompute Residuals */
    /* 	  bminax(f,A.IA,A.JA,A.val,u,&ndof,res); */
    /* 	}  */
    /*     } */
    /*   } */
    /* } */
    /* clk1 = clock(); */
    /* printf("Elapsed CPU Time for Solve = %f seconds.\n\n",(REAL) (clk1-clk0)/CLOCKS_PER_SEC); */
    /* /\*******************************************************************************************\/ */
    
    /* /\**************** Get True Solution and Error***********************************************\/ */
    /* utrue = calloc(ndof,sizeof(REAL)); */
    /* myerr = calloc(ndof,sizeof(REAL)); */
    /* if(havetrue) { */
    /*   printf("Computing True Solution and Errors (if available):\n"); */
    /*   clk0 = clock(); */
    
    /*   REAL udott = -6.66; */
    /*   if(mydim==3) { */
    /*     for (i=0; i<n; i++) { */
    /* 	getknownfunction(&udott,cn.x[i],cn.y[i],cn.z[i],0.0,mydim,1,0.0,ut); */
    /* 	utrue[i]=udott; */
    /* 	myerr[i] = fabs(u[i]-utrue[i]); */
    /* 	myerr[i] = u[i]-utrue[i]; */
    /*     } */
    /*   } else if (mydim==2) { */
    /*     for (i=0; i<n; i++) { */
    /* 	getknownfunction(&udott,cn.x[i],cn.y[i],0.0,0.0,mydim,1,coef,ut); */
    /* 	utrue[i]=udott; */
    /* 	myerr[i] = fabs(u[i]-utrue[i]); */
    /* 	myerr[i] = u[i]-utrue[i]; */
    /*     } */
    /*   } else { */
    /*     printf("bad dimension\n"); */
    /*     exit(0); */
    /*   } */
    
    /*   xdy(&ndof,myerr,myerr,&errnorml2); */
    /*   L2normMassassemble(&errnorm,myerr,el_n.IA,el_n.JA,el_n.IA,el_n.JA,cn.x,cn.y,cn.z,cq.x,cq.y,cq.z,cq.w,nelm,nq1d,mydim,element_order,element_order,NULL,NULL,NULL,NULL,0); */
    /*   clk1 = clock(); */
    /*   printf("Elapsed CPU Time for Getting Errors = %f seconds.\n\n",(REAL) (clk1-clk0)/CLOCKS_PER_SEC); */
    /* } */
    /* /\*******************************************************************************************\/ */
    
    /* /\**************** Print Results or Dump Results *******************************************\/ */
    /* if (dumpsol>=2) {  */
    /*   uid = fopen("output/sol.dat","w");  */
    /*   if(havetrue) { truid = fopen("output/true.dat","w"); } */
    /* } */
    /* dumpsolution(u,utrue,uid,truid,ndof,n,nelm,dumpsol,1,mydim,havetrue); */
    /* if(dumpmat>=2) { */
    /*   matid = fopen("output/mat.dat","w"); */
    /*   rhsid = fopen("output/rhs.dat","w"); */
    /* } */
    /* dumpmatrices(A.IA,A.JA,A.val,f,matid,rhsid,A.row,dumpmat); */
    /* if (havetrue) { */
    /*   printf("***********************************************************************************************\n"); */
    /*   printf("L2 Norm of u error = %25.17g\n",errnorm); */
    /*   printf("l2 Norm of error =   %25.17g\n",sqrt(errnorml2)); */
    /*   printf("***********************************************************************************************\n"); */
    /*   if(dumpsol>=2) {fclose(truid);} */
    /* } */
    /* if (dumpsol>=2) { */
    /*   fclose(uid); */
    /* } */
    /* if (dumpmat>=2) { */
    /*   fclose(rhsid); */
    /*   fclose(matid); */
    /* } */
    /* /\*******************************************************************************************\/ */
    
    /******** Free All the Arrays ***********************************************************/
    dcsr_free(&A);
    if(b.val) free(b.val);
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


