/*
 *  stokesmain.c
 *  
 *
 *  Created by James Adler on 6/8/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <getopt.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>


/*
 *  Discussion:
 *
 *    This program solves the Laplace Equation using P1 or P2 finite elements
 *
 *      -Laplace u = 0
 *
 *    in a 2D or 3D rectangular region in the plane.
 *
 *    Along the boundary of the region, homogeneous Dirichlet conditions
 *    are imposed:
 *
 *      u = 0
 *
 *
 *
 *    INPUT:
 *			gridfile = File containing finite element mesh
 *          (nodes,elements,boundaries)
 *
 *    OUTPUT: 
 *         Creates matrices and solves using Minres preconditioned algorithm
 */

/*********** EXTERNAL FUNCTIONS ****************************************************************/
#include "fem_mhd_functs.h"
/***********************************************************************************************/

/****** MAIN DRIVER ****************************************************************************/
int main (int argc, char* argv[]) 
{
	
  printf("\nBeginning Program to create Finite Element Matrices for Laplace\n");
	
  /****** INITIALIZE PARAMETERS **************************************************************/
  INT i,j; /* Loop indices */
  clock_t clk0,clk1,clka,clkb; /* timing parameters */
  clka = clock();
	
  // Grid data
  FILE* gfid;
  char gridfile[50] = "Data/grids/2D/hp52d.dat";	/* File name for grid data file */	
  INT input_grid = 1;								/* Set to 0 for default grid to be used */
	
  // Grab parameters from input.dat file
  INT ipar[55];
  REAL fpar[20];
  getinput(gridfile,fpar,ipar);
	
  // Dimension is needed for all this to work
  INT mydim = ipar[0];
	
  // How many nodes and edges per element?
  // This is default for standard triangles and tetrahedra
  INT poly = ipar[2];				       	/* Order of Elements */
  INT nve = mydim+1;					/* Vertices per element */
  INT edge_order = 3*(mydim-1);		                /* Edges per element */
  INT element_order = nve+(poly-1)*edge_order;	        /* Total Nodes per Element (dim+1 for linears) */
  INT nq1d = ipar[1];					/* Quadrature points per dimension */
  INT nq = pow(nq1d,mydim);			        /* Quadrature points per element */
	
  // Properties of Mesh
  INT nelm = 0;				        	/* Number of Elements */
  INT n = 0;				       		/* Number of Nodes */
  INT nvert = 0;		       			/* Number of vertices total */
  INT nedge = nelm+n-(mydim-1);	                	/* Number of Edges */
  INT nbedge = 0;			       		/* Number of edges on boundary */
  INT nbvert = 0;		       			/* Number of nodes on boundary */
  INT ndof;		       				/* Generic degrees of freedom */
	
  // Discrete System Properties
  INT maxit=ipar[41];	       				/* Max Number of PCG iterations */
  REAL diag=-666.66;                                  /* Diagonal of Matrix */
  INT outprint=ipar[39];		       		/* Display CG Residuals */
  REAL pcgtol=fpar[6];                                /* Tolerance for PCG solving Ax=b */
  INT solvertype = ipar[40];                            /* Indicates type of solver, for now: 1->MINRES 2->GMRES */
	
  // Domain Boundaries (if no domain file given.  Assumes rectangle for now)
  REAL xL = 0.0;
  REAL xR = 1.0;
  REAL yB = 0.0;
  REAL yT = 1.0;
  REAL zD = 0.0;
  REAL zU = 1.0;
	
  // Parameters
  REAL coef = fpar[3];
  INT compRHS = ipar[20];			       	/* Right-hand size = mydim*pi^2*sin(pi*x)*sin(pi*y)*sin(pi*z) */
  INT compB = ipar[13];			       		/* Zero Dirichlet Boundary Condition (2D or 3D) */
  INT ut = ipar[5];
		
  // Solution Stuff
  INT havetrue = ipar[4];	       			/* Indicates whether there exists an analytic solution to compute error. */
  INT cola,colb,jaa;	       				/* Loop indices for prINTing matrix */
  INT dumpmat = ipar[38];      				/* Dump Matrix? 0 -> No  1 -> Just to Screen  2 -> Just to File  3 -> Screen and File */
  INT dumpsol = ipar[36];      				/* Dump Solution? same as above */
  INT dumpmesh = ipar[37];     				/* Dump Mesh data? 1 yes 0 no */

  FILE* matid;	       					/* File to dump matrix A */
  FILE* rhsid;	       					/* File to dump rhs b */
  FILE* truid;	       					/* File to dump true solution */
  FILE* uid;         				       	/* File to dump numerical solution */
	
  // Temporary Grid stuff for my fake grid
  INT nx;		                                /* Number of nodes in each direction */
  INT ny;
  INT nex = 2;                                          /* Number of macroelements (squares) in each direciton */
  INT ney = 2;
  /*******************************************************************************************/
	
  /*** Get Input arguments if any ************************************************************/
  for (i=1; i < argc; i++) {
    if (strcmp("-nex", argv[i])==0) { nex = atoi(argv[++i]);}
    if (strcmp("-ney", argv[i])==0) { ney = atoi(argv[++i]);}
    if (strcmp("-dim", argv[i])==0) { mydim = atoi(argv[++i]);}
    if (strcmp("-xL", argv[i])==0) { xL = atof(argv[++i]);}
    if (strcmp("-xR", argv[i])==0) { xR = atof(argv[++i]);}
    if (strcmp("-yB", argv[i])==0) { yB = atof(argv[++i]);}
    if (strcmp("-yT", argv[i])==0) { yT = atof(argv[++i]);}
    if (strcmp("-poly", argv[i])==0) { poly = atoi(argv[++i]);}
    if (strcmp("-nq", argv[i])==0) { nq1d = atoi(argv[++i]);}
    if (strcmp("-grid", argv[i])==0) { strcpy(gridfile,argv[++i]); input_grid=1;}
    if (strcmp("-true", argv[i])==0) { havetrue = 1; }
    if (strcmp("-dump", argv[i])==0) { dumpmat = atoi(argv[++i]); }
  }
  /*******************************************************************************************/
	
  /** INITIALIZE ANY ARRAYS NEEDED ***********************************************************/
  INT* element_node=NULL;	/* Element to Node map Given as Input (NOT CSR FORMAT) */
  CSRinc el_v;                  /* Element to Vertex Map (CSR FORMAT) */
  coordinates cv;               /* Coordinates of vertices */
  CSRinc el_n;                  /* Element to Node Map */
  coordinates cn;               /* Coordinates of nodes */
  CSRinc ed_n;                  /* Edge to Vertex Map */
  CSRinc el_ed;	        	/*  Element to Edge Map (CSR Format) */
  INT* ed_bdry=NULL;	       	/* Indicates whether an edge is a boundary */
  INT* v_bdry=NULL;	       	/* Indicates whether a vertex is a boundary */
  INT* n_bdry=NULL; 	      	/* Indicates whether a node is a boundary */
  INT* bdry_v=NULL; 		/* Indicates nodes of an edge on boundary */
  qcoordinates cq;              /* Quadrature nodes and weights */
  dCSRmat A;                    /* Global stiffness matrix */
  REAL* f=NULL;	       	/* Global right-hand side vector */
  REAL* u=NULL;	       	/* Solution vector (represents u*tau) */
  REAL* utrue=NULL;       	/* True Solution vector */
  REAL* myerr=NULL;           /* Error Vector */
  REAL errnorm=0.0;     	/* L2 Norm of u1 Error */
  REAL errnorml2=0.0;    	/* Little l2 error */
  /*******************************************************************************************/
	
  /******* Input Grid ********************************************************/
  // Either input from file or use the default grid:
  clk0 = clock();
  if(input_grid==0) {
    printf("\nWhat!? You don't have your OWN grid...Creating Default Grid:->\n");
		
    //Get number of nodes in each direction and total
    nx = nex/2 + 1;
    ny = ney/2 + 1;
    nvert = nx*ny;
		
    // Get total number of elements
    nelm = 2*(nex/2)*(ney/2);
		
    // Get total number of edges
    nedge = nelm+nvert-(mydim-1);	/* Number of Edges */
		
    // Get total number of boundaries
    nbedge = 2*(nx-1)+2*(ny-1);
    nbvert = nbedge;
		
    // Allocate arrays needed
    bdry_v = calloc(nbedge*2,sizeof(INT));
    allocatecoords(&cv,nvert,mydim);
    element_node = (INT *) calloc((nelm) * (nvert),sizeof(INT));
		
    // Create Default Grid
    std_tri_grid(element_node,cv.x,cv.y,nex,ney,xL,xR,yB,yT,bdry_v);
		
  } else {
    printf("\nLoading Grid From File: %s->\n",gridfile);
		
    // Open file for reading
    gfid = fopen(gridfile,"r");
    if( gfid == NULL ) { 
      printf("\nError opening Grid File!!!\nFile (%s) probably doesn't exist...\nAborting Program.\n\n",gridfile);
      return 0;
    }
		
    if(mydim==2) {
      // Get Number of elements, nodes and boundary edges first
      INT* line1 = calloc(4,sizeof(INT));
      INT lenhead = 4;
      rveci_(gfid,line1,&lenhead);
      nelm = line1[0];
      nvert = line1[1];
      nbedge = line1[2];
      nedge = nelm+nvert-(mydim-1); /* Number of Edges */
      free(line1);
			
      // Allocate arrays
      element_node = calloc(nelm*nve,sizeof(INT));
      allocatecoords(&cv,nvert,mydim);
      bdry_v = calloc(nbedge*2,sizeof(INT));
			
      // Read in data
      readgrid2D(gfid,element_node,cv.x,cv.y,bdry_v,nelm,nvert,nbedge,nve);
    } else {
			
      // Get Number of elements, nodes and boundary edges first
      INT* line1 = calloc(4,sizeof(INT));
      INT lenhead = 4;
      rveci_(gfid,line1,&lenhead);
      nelm = line1[0];
      nvert = line1[1];
      free(line1);
			
      // Allocate arrays
      element_node = calloc(nelm*nve,sizeof(INT));
      allocatecoords(&cv,nvert,mydim);
      v_bdry = calloc(nvert,sizeof(INT));
      // Read in data
      readgrid3D(gfid,element_node,cv.x,cv.y,cv.z,v_bdry,nelm,nvert,nve,&nbvert);
    }
    fclose(gfid);
		
  }
  /*******************************************************************************************/
  /*** Convert Map matrices to Sparse Matrix Format ******************************************/
  printf("\nConverting Grid Maps to CSR and Computing Data Structures:\n ");
	
  /* How many Nodes per element? */	
  if (poly==1) {
    element_order = nve;
  } else if (poly==2) {
    element_order = nve+edge_order;
  }
	
  /* Element Vertex Map */
  allocateCSRinc(&el_v,nelm,nvert,nelm*nve);
  convert_elmnode(el_v.IA,el_v.JA,element_node,nelm,nvert,nve);
  if(element_node) free(element_node);
	
  /* Edge to Node Map */
  if (mydim==3) { get_nedge(&nedge,el_v.IA,el_v.JA,nvert,nelm,nve); }
  allocateCSRinc(&ed_n,nedge,nvert,2*nedge);
  get_edge_n(ed_n.IA,ed_n.JA,nedge,el_v.IA,el_v.JA,nvert,nelm,nve);
	
  /* Get Boundary Edges and Nodes */
  if (mydim==2) {
    ed_bdry = calloc(nedge,sizeof(INT));
    v_bdry = calloc(nvert,sizeof(INT));
    isboundary_ed(ed_n.IA,ed_n.JA,nedge,nbedge,bdry_v,ed_bdry);
    isboundary_n(nvert,bdry_v,v_bdry,nbedge,&nbvert);
    if(bdry_v) free(bdry_v);
  } else {
    ed_bdry = calloc(nedge,sizeof(INT));
    isboundary_ed3D(ed_n.IA,ed_n.JA,nedge,cv.x,cv.y,cv.z,xL,xR,yB,yT,zD,zU,&nbedge,v_bdry,ed_bdry);
  }
	
  /* Element to Edge Map */
  INT nnz_eled=0;
  INT* in_ed = calloc(nvert+1,sizeof(INT));
  INT* jn_ed = calloc(2*nedge,sizeof(INT));
  atransps(ed_n.IA,ed_n.JA,nedge,nvert,in_ed,jn_ed);
  get_nnz(el_v.IA,el_v.JA,nelm,nvert,nedge,in_ed,jn_ed,&nnz_eled);
  allocateCSRinc(&el_ed,nelm,nedge,nnz_eled);
  abybs_mult(el_v.IA,el_v.JA,nelm,nvert,nedge,el_ed.IA,el_ed.JA,in_ed,jn_ed,2);
  if(in_ed) free(in_ed);
  if(jn_ed) free(jn_ed);	
	
  /* Get Quadrature Nodes */
  nq = pow(nq1d,mydim);
  allocateqcoords(&cq,nq1d,nelm,mydim);	
  if (mydim==2) {
    get_quadrature(cq.x,cq.y,NULL,cq.w,cv.x,cv.y,NULL,el_v.IA,el_v.JA,nelm,nve,nq1d,mydim);
  } else {
    get_quadrature(cq.x,cq.y,cq.z,cq.w,cv.x,cv.y,cv.z,el_v.IA,el_v.JA,nelm,nve,nq1d,mydim);
  }
	
  // Get Higher Order Grids if necessary
  if(poly==1) {
    n = nvert;
  } else if(poly==2) {
    n = nvert+nedge;
    nbvert = nbvert+nbedge;
  }
  allocatecoords(&cn,n,mydim);
  allocateCSRinc(&el_n,nelm,n,nelm*element_order);
  n_bdry = calloc(n,sizeof(INT));
	
  if (mydim==3) {
    if (poly==1) {
      for(i=0;i<n;i++) {
	cn.x[i] = cv.x[i];
	cn.y[i] = cv.y[i];
	cn.z[i] = cv.z[i];
	n_bdry[i] = v_bdry[i];
      }
      for (i=0; i<nelm+1; i++) {
	el_n.IA[i] = el_v.IA[i];
      }
      for (i=0; i<nelm*element_order; i++) {
	el_n.JA[i] = el_v.JA[i];
      }
    } else if (poly==2) {
      get_P2(cn.x,cn.y,cn.z,el_n.IA,el_n.JA,cv.x,cv.y,cv.z,el_v.IA,el_v.JA,nve,element_order,el_ed.IA,el_ed.JA,ed_n.IA,ed_n.JA,nvert,mydim,nelm,nedge,n_bdry,v_bdry,ed_bdry);
    }	
  } else { // 2D		
    if (poly==1) {
      for	(i=0;i<n;i++) {
	cn.x[i] = cv.x[i];
	cn.y[i] = cv.y[i];
	n_bdry[i] = v_bdry[i];
      }
      for (i=0; i<nelm+1; i++) {
	el_n.IA[i] = el_v.IA[i];
      }
      for (i=0; i<nelm*element_order; i++) {
	el_n.JA[i] = el_v.JA[i];
      }
    } else if (poly==2) {
      get_P2(cn.x,cn.y,cn.z,el_n.IA,el_n.JA,cv.x,cv.y,cv.z,el_v.IA,el_v.JA,nve,element_order,el_ed.IA,el_ed.JA,ed_n.IA,ed_n.JA,nvert,mydim,nelm,nedge,n_bdry,v_bdry,ed_bdry);
    }
  }
  if(dumpmesh==1) {
    dumpmeshdata(el_n.IA,el_n.JA,el_v.IA,el_v.JA,cn.x,cn.y,cn.z,nelm,n,nvert,element_order,mydim+1,mydim,40,n_bdry,n);
  }
	
  freecoords(cv);
  freeCSRinc(el_v);
  if(v_bdry) free(v_bdry);
	
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
  freecoords(cn);
  freeCSRinc(el_n);
  freeCSRinc(ed_n);
  freeCSRinc(el_ed);
  freeqcoords(cq);
  if(f) free(f);
  if(ed_bdry) free(ed_bdry);
  if(n_bdry) free(n_bdry);
  if(u) free(u);
  if(utrue) free(utrue);
  if(myerr) free(myerr);
  freedCSRmat(A);
  /*******************************************************************************************/
	
  clkb = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",(REAL) (clkb-clka)/CLOCKS_PER_SEC);
  return 0;	
	
}	/* End of Program */
/*******************************************************************************************/


