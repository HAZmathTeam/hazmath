/*
 *  Maxwell.c
 *  
 *
 *  Created by Xiaozhe Hu and James Adler on 10/17/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "hazmat.h"

/*
 *  Discussion:
 *
 *    This program solves the Maxwell's equation system for disappearing solutions
 *    using Raviart-Thomas and Nedelec elements (this is modified/simplified by Xiaozhe and I).
 *
 *      eps* dE/dt - curl(1/mu B) = -j
 *      dB/dt + curl(E) = 0
 *      div(eps*E) = 0
 *      div(B) = 0
 *		
 *    in a 3D rectangular region in the plane. with an interior sphere (or cube) removed
 *
 *    Along the boundary of the region, impedence conditions are
 *    are imposed in the interior boundary and perfect conductor conditions on the exterior boundary
 *
 *      Exterior: n * B = 0, n x E = 0 
 *      Interior: -n x (n x E) = -n x (1/mu B)
 *
 *    INPUT:
 *		gridfile = File containing finite element mesh (nodes,elements,boundaries)
 *
 *    OUTPUT: 
 *              Creates matrices and solves
 */

/*********** EXTERNAL FUNCTIONS ****************************************************************/
#include "fem_mhd_functs.h"
/***********************************************************************************************/

/****** MAIN DRIVER ****************************************************************************/
int main (int argc, char* argv[]) 
{
	
  printf("\nBeginning Program to Solve Maxwell's Equations with Disappearing Solutions\n");
	
  /****** INITIALIZE PARAMETERS **************************************************************/
  INT i,j; /* Loop indices */
	
  // Timing Parameters
  clock_t clk0,clk1,clk2,clk3,clka,clkb;
  clka = clock();
	
  // Grid data
  FILE* gfid;
  char gridfile[50] = "Data/grids/2D/hp52d.dat";	/* File name for grid data file */	
  INT input_grid = 1;					/* Set to 0 for default grid to be used */
	
  // Grab parameters from input.dat file
  INT ipar[45];
  REAL fpar[20];
  getinput(gridfile,fpar,ipar);
	
  // Dimension is needed for all this to work
  INT mydim = ipar[0];
	
  // How many nodes and edges per element?
  // This is default for standard triangles and tetrahedra
  INT poly = ipar[2];		       		/* Order of Elements */
  INT nve = mydim+1;	       			/* Vertices per element */
  INT element_order = nve;	       		/* Total Nodes per Element (dim+1 for linears) */
  INT edge_order = 3*(mydim-1);	         	/* Edges per element */
  INT face_order = nve;                         /* Faces per element (in 2D these are edges) */
  INT ndpf = mydim;                             /* Number of nodes per face = dimension */
  INT edpf;                                     /* Number of edges per face */
 
  INT nq1d = ipar[1];	       			/* Quadrature points per dimension */
  INT nq = pow(nq1d,mydim);    		        /* Quadrature points per element */
	
  // Properties of Mesh
  INT nelm;	       			 	/* Number of Elements */
  INT n;       		   	      		/* Number of Nodes */
  INT nvert;          		         	/* Number of vertices total */
  INT nedge;		              	       	/* Number of Edges */
  INT nbedge;		      	       	       	/* Number of edges on boundary */
  INT nbvert;		              		/* Number of vertices on boundary */
  INT ndof;                                     /* Generic degrees of freedom */
  INT nface;                                    /* Number of Faces */
  INT nbface;                                   /* Number of Faces on Boundary */
	
  // Discrete System Properties
  INT maxit=ipar[41];		       		/* Max Number of PCG iterations */
  REAL diag=-666.66;                          /* Diagonal of Matrix */
  INT outprint=ipar[39];		   	/* Display CG Residuals */
  REAL solvetol=fpar[6];                      /* Tolerance for PCG solving Ax=b */
	
  // Solution Stuff
  INT havetrue = ipar[4];		       	/* Indicates whether there exists an analytic solution to compute error. */
  INT E1t = ipar[5];		              	/* Indicates truesolution of E1 */
  INT E2t = ipar[6];	       	       		/* Indicates truesolution of E2 */
  INT E3t = ipar[7];   	       			/* Indicates truesolution of E3 */
  INT B1t = ipar[9];  				/* Indicates truesolution of B1 */
  INT B2t = ipar[10];	      		        /* Indicates truesolution of B2 */
  INT B3t = ipar[11];	      		        /* Indicates truesolution of B3 */
  INT dumpmat = ipar[38];      			/* Dump Matrix? 0 -> No  1 -> Just to Screen  2 -> Just to File  3 -> Screen and File */
  INT dumpsol = ipar[36];      			/* Dump Solution? 0 -> No  1 -> Just to Screen  2 -> Just to File  3 -> Screen and File */
  INT dumpmesh = ipar[37];    		        /* Dump Mesh? 0 -> No 1 -> Yes */
  INT solvertype = ipar[40];   	       		/* Indicates type of solver, for now: 0->CG 1->MINRES 2->GMRES 3->Precond MINRES 4-> Precond GMRES */
  INT restart = ipar[42];                       /* Restart value if GMRES used */
  
  // Time stepping parameters
  INT nsteps = ipar[34];	       	      	/* Max number of Time stepping steps */
  REAL dt = fpar[4];			       	/* Size of time step */
  REAL time = 0.0;
  INT E1init = ipar[26];	       		/* Initial condition for E1 */
  INT E2init = ipar[27];       			/* Iniital condition for E2 */
  INT E3init = ipar[28];       		       	/* Initial condition for E3 */
  INT B1init = ipar[30];       	       		/* Initial condition for B1 */
  INT B2init = ipar[31];       	       		/* Initial condition for B2 */
  INT B3init = ipar[32];       	       		/* Initial condition for B3 */
  INT istimdep = 1;	       	      	       	/* Indicates if problem is time dependent 1 yes 0 no */
  INT timscheme = ipar[33];                     /* Indicates type of timestepping 0-CN 1-BDF1 */
  INT btimdep = ipar[12];                       /* Indicates if boundaries are time-dependent 1-yes 0-no */
	
  // Kappa and other coefficients for Test Problem
  REAL mu = fpar[0];			      	/* Permeability coefficient  */
  REAL oneovermu = 1.0/mu;
  REAL myeps = fpar[2];	      	        	/* Permitivity coefficient */
  REAL alpha = fpar[3];	       	                /* Indicates time-dependence */
  REAL myb = fpar[1];                           /* Parameter in true solution and initial condition */
  REAL gam = fpar[7];                           /* Function for gamma(t) in Z marix */
  INT gamistimedep = 0;                         /* Indicates if gamma is time-dependent */
  if(gam==-666.66) {                            /* That means gamma is time-dependent and we use a function of t */
    gamistimedep = 1;
  }
  REAL* myB=NULL;                               /* B = 1/r curl E */
  REAL* myBt=NULL;
  REAL myr = 0.5*(1-sqrt(1+4/gam));             /* r = 1/2(1 - sqrt(1 + 4/gam)) */

  if (alpha==-666.66) {
    /* -666.66 is the fail safe or test */
    istimdep = 1;
    if (timscheme==1) { // BDF-1 Backward Euler u_n/dt + L(u_n) = f + u_0/dt
      alpha=1.0/dt;
    } else if(timscheme==0) { // Crank-Nicholson 2*u_n/dt + L(u_n) = f + 2*u_0/dt - L(u_0)
      alpha = 2.0/dt;
    }
  }
  INT E1bdry= ipar[13];			/* Dirichlet Boundary Condition for u variables (1 = Homogeneous) */
  INT E2bdry = ipar[14]; 
  INT E3bdry = ipar[15];       	       	// Indicates what type of boundary condition for u1
  INT B1bdry = ipar[17];
  INT B2bdry = ipar[18];
  INT B3bdry = ipar[19];
  INT j1rhs = ipar[20];		       		/* RHS function for Ampere's law = j_1 */
  INT j2rhs = ipar[21];			       	/* RHS function for Ampere's law = j_2 */
  INT j3rhs = ipar[22];		       	       	/* RHS function for Ampere's law = j_3 */
  REAL uval;
	
  // Temporary Grid stuff for my fake grid
  INT nx;		          			       	/* Number of nodes in each direction */
  INT ny;
  INT nex = 2;			       	       			/* Number of macroelements (squares) in each direciton */
  INT ney = 2;
  /***********************************************************************************************/
	
  /** INITIALIZE ANY ARRAYS NEEDED ***********************************************************/
  // Incidence Maps and Coordinates
  INT* element_node=NULL;    /* Element to Node map Given as Input (NOT CSR FORMAT) */
  CSRinc el_v;               /* Element to Vertex Map */
  coordinates cv;            /* Coordinates of vertices */
  CSRinc el_n;               /* ELement to Node map */
  coordinates cn;            /* Coordinates of nodes */
  CSRinc ed_n;               /* Edge to Vertex Map */
  CSRinc el_ed;	             /* Element to Edge Map (CSR Format) */
  REAL* ed_len=NULL;         /* Length of each Edge */
  REAL* ed_tau=NULL;         /* Tangent Vector on each edge nedge x mydim */
  REAL* ed_mid=NULL;         /* MidpoINT of each edge */
  INT* ed_bdry=NULL;	     /* Indicates whether an edge is a boundary */
  INT* v_bdry=NULL;          /* Indicates whether a vertex is a boundary */
  INT* n_bdry=NULL; 	     /* Indicates whether a node is a boundary */
  INT* dof_bdry=NULL;        /* Indicates whether a node is a boundary for each dof */
  INT* bdry_v=NULL; 	     /* Indicates nodes of an edge on boundary */
  CSRinc face_n;             /* Face to Node Map */
  iCSRmat el_face;           /* Element to Face Map */
  INT* fel_order=NULL;       /* Indicates ordering of faces on given element using node ordering */
  REAL* f_area=NULL;         /* Area of each Face (Length in 2D) */
  REAL* f_norm=NULL;         /* Normal vector for each face nface x mydim (based on opposite node from low element to high element or out on bdry) */
  REAL* f_mid=NULL;          /* Midpoint of each face */
  INT* face_bdry=NULL;       /* Indicates whether a face is a boundary */
  CSRinc f_ed;               /* Face to Edge Map */
  qcoordinates cq;           /* Quadrature nodes and weights */				
  dCSRmat Z;                 /* Impedence Boundary Matrix */
  REAL* Ztmp;                /* Holds values of (1+gamma)*Z */
  dCSRmat K;                 /* Incidence matrix of face to edge map (includes signs) |E|*K*|F|^(-1) = Curl Operator */
  dCSRmat G;                 /* Incidence matrix of edge to node map (includes signs) |E|^(-1)*Ggrad = grad operator */
  dCSRmat Me;                /* Mass Matrix for Edge DOF */
  dCSRmat Mf;                /* Mass Matrix for Face DOF */
  dCSRmat Mv;                /* Mass Matrix for Vertex DOF */
  dCSRmat MG;                /* Me*G */
  dCSRmat MGt;               /* G'*Me */
  dCSRmat MK;                /* Mf*K */
  dCSRmat MKt;               /* K'*Mf */
  dCSRmat GMG;               /* G'MeG */
  dCSRmat KMK;               /* K'MfK */
  REAL* f=NULL;	         	/* Global right-hand side vector */
  REAL* u=NULL;           	/* Solution vector */
  REAL* utrue=NULL;        /* True Solution vector */
  REAL* uprev=NULL;       	/* Solution at previous time step */
  REAL* uzero=NULL;        /* Stores the zero vector if needed */
  /*******************************************************************************************/
  
  /******* Input Grid ********************************************************/
  // Either input from file or use the default grid:
  clk0 = clock();
  if(input_grid==0) {
    printf("\nWhat!? You don't have your OWN grid...DYING NOW:->\n");
  } else {
    printf("\nLoading Grid From File:->%s\n",gridfile);
		
    // Open file for reading
    gfid = fopen(gridfile,"r");
    if( gfid == NULL ) { 
      printf("\nError opening Grid File!!!\nFile probably doesn't exist...\nAborting Program.\n\n");
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
      cv.z = NULL;
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
  fflush(stdout);	
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
  }

  /* Element to Edge Map */
  INT nnz_eled=0;
  INT* in_ed = calloc(nvert+1,sizeof(INT));
  INT* jn_ed = calloc(2*nedge,sizeof(INT));
  atransps(ed_n.IA,ed_n.JA,nedge,nvert,in_ed,jn_ed);
  get_nnz(el_v.IA,el_v.JA,nelm,nvert,nedge,in_ed,jn_ed,&nnz_eled);
  allocateCSRinc(&el_ed,nelm,nedge,nnz_eled);
  abybs_mult(el_v.IA,el_v.JA,nelm,nvert,nedge,el_ed.IA,el_ed.JA,in_ed,jn_ed,2);
  ed_len = calloc(nedge,sizeof(REAL));
  ed_tau = calloc(nedge*mydim,sizeof(REAL));
  ed_mid = calloc(nedge*mydim,sizeof(REAL));
  edge_stats_all(ed_len,ed_tau,ed_mid,cv.x,cv.y,cv.z,ed_n.IA,ed_n.JA,mydim,nedge);

  // Face to Node/Element etc stats
  INT euler,nholes;
  if(mydim==2) {
    nface = nedge;
    euler = nvert - nedge + nelm;
  } else if (mydim==3) {
    nface = 1 + nedge-nvert+nelm;
    nholes=0;
    nface = nface + nholes; // add number of holes!
    euler = nvert - nedge + nface - nelm;
  } else {
    baddimension();
  }
  if(euler!=1+nholes) {
    printf("Your simplices are all messed up.  Euler Characteristic doesn't equal 1!");
    exit(0);
  }

  allocateiCSRmat(&el_face,nelm,nface,face_order*nelm);
  face_bdry = calloc(nface,sizeof(INT));
  allocateCSRinc(&face_n,nface,nvert,ndpf*nface);
  fel_order = calloc(face_order*ndpf,sizeof(INT));
  f_area = calloc(nface,sizeof(REAL));
  f_norm = calloc(nface*mydim,sizeof(REAL));
  f_mid = calloc(nface*mydim,sizeof(REAL));
  get_face_maps(el_v.IA,el_v.JA,nelm,nvert,nve,nface,mydim,face_order,el_face.IA,el_face.JA,el_face.val,face_bdry,&nbface,face_n.IA,face_n.JA,fel_order);
  fix_facebdry(face_bdry,v_bdry,nface,nvert,face_n.IA,face_n.JA);
  face_stats(f_area,f_norm,nface,cv.x,cv.y,cv.z,el_face.IA,el_face.JA,el_face.val,face_n.IA,face_n.JA,mydim,face_order,fel_order,el_v.IA,el_v.JA,nve,nelm,f_mid);

  /* Face to Edge Map */
  INT nnz_fed=0;
  if(mydim==2) {
    edpf = 1;
  } else if(mydim==3) {
    edpf = 3;
  } else {
    baddimension();
  }
  get_nnz(face_n.IA,face_n.JA,nface,nvert,nedge,in_ed,jn_ed,&nnz_fed);
  CSRinc face_ed;
  allocateCSRinc(&face_ed,nface,nedge,nnz_fed);
  abybs_mult(face_n.IA,face_n.JA,nface,nvert,nedge,face_ed.IA,face_ed.JA,in_ed,jn_ed,2);
  allocateCSRinc(&f_ed,nface,nedge,nface*edpf);
  for(i=0;i<nface+1;i++)
    f_ed.IA[i] = face_ed.IA[i];
  for(i=0;i<nface*edpf;i++)
    f_ed.JA[i] = face_ed.JA[i];
  freeCSRinc(face_ed);
  ed_bdry = calloc(nedge,sizeof(INT));
  isboundary_ed3D_new(f_ed.IA,f_ed.JA,nface,face_bdry,nedge,&nbedge,ed_bdry);
  if(in_ed) free(in_ed);
  if(jn_ed) free(jn_ed);
 

  /* Get Quadrature Nodes */
  nq = nq1d;  for (j=1; j<mydim; j++) nq=nq*nq1d;
  allocateqcoords(&cq,nq1d,nelm,mydim);
  if(mydim==3) {
    get_quadrature(cq.x,cq.y,cq.z,cq.w,cv.x,cv.y,cv.z,el_v.IA,el_v.JA,nelm,nve,nq1d,mydim);
  } else {
    get_quadrature(cq.x,cq.y,NULL,cq.w,cv.x,cv.y,NULL,el_v.IA,el_v.JA,nelm,nve,nq1d,mydim);
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
      get_P2(cn.x,cn.y,cn.z,el_n.IA,el_n.JA,cv.x,cv.y,cv.z,el_v.IA,el_v.JA,nve,element_order,el_ed.IA,el_ed.JA,ed_n.IA,ed_n.JA,nvert,mydim,nelm,nedge, \
  	     n_bdry,v_bdry,ed_bdry);
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
      get_P2(cn.x,cn.y,cn.z,el_n.IA,el_n.JA,cv.x,cv.y,cv.z,el_v.IA,el_v.JA,nve,element_order,el_ed.IA,el_ed.JA,ed_n.IA,ed_n.JA,nvert,mydim,nelm,nedge, \
  	     n_bdry,v_bdry,ed_bdry);
    }
  }

  // Reorder Face Node map for orientaion
  sync_facenode(face_n.IA,face_n.JA,ndpf,mydim,nface,f_norm,cn.x,cn.y,cn.z);

  // For each unknown get boundary condition
  ndof = n + nedge + nface;
  INT* edzero = calloc(nedge,sizeof(INT));
  INT* fzero = calloc(nface,sizeof(INT));
  INT* nzero = calloc(n,sizeof(INT));
  INT* isdirichletE = (INT *) calloc(nedge,sizeof(INT));
  INT* isdirichletF = (INT *) calloc(nface,sizeof(INT));
  INT* isdirichletV = (INT *) calloc(n,sizeof(INT));
  INT* isdirichlet = (INT *) calloc(ndof,sizeof(INT));
  dof_bdry = (INT *) calloc(n+nedge+nface,sizeof(INT));
  for(j=0;j<nedge;j++) {
    dof_bdry[j] = ed_bdry[j];
    if(ed_bdry[j]<=0) {
      isdirichletE[j]=0;
      isdirichlet[j]=0;
    } else {
      isdirichletE[j]=1;
      isdirichlet[j]=1;
    }
    edzero[j] = 0;
  }
  for(i=0;i<nface;i++) {
    dof_bdry[i+nedge] = face_bdry[i];
    if(face_bdry[i]<=0) {
      isdirichletF[i]=0;
      isdirichlet[i+nedge]=0;
    } else {
      isdirichletF[i]=1;
      isdirichlet[i+nedge]=1;
    }
    fzero[i]=0;
  }
  for(j=0;j<n;j++) {
    dof_bdry[j+nedge+nface] = n_bdry[j];
    if(n_bdry[j]==0) {
      isdirichletV[j]=0;
      isdirichlet[j+nedge+nface]=0;
    } else {
      isdirichletV[j]=1;
      isdirichlet[j+nedge+nface]=1;
    }
    nzero[j] = 0;
  }
  
  if(dumpmesh==1) {
    dumpmeshdata(el_n.IA,el_n.JA,el_v.IA,el_v.JA,cn.x,cn.y,cn.z,nelm,n,nvert,element_order,nve,mydim,40,dof_bdry,n+nedge+nface);
  }
  
  clk1 = clock();
  printf("Elapsed CPU Time for Conversion = %f seconds.\n\n",(REAL) (clk1 - clk0)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  printf("***********************************************************************************************\n");
  printf("Element Type: Nedelec for E, Raviart-Thomas for B\n");
  printf("Dimension = %d\t\t\tNumber of Elements = %d\n",mydim,nelm);
  printf("Number of Nodes = %d\t\tNumber of Edges = %d\t\tNumber of Faces = %d\n",n,nedge,nface);
  printf("Number of Boundary Nodes = %d\tNumber of Boundary Edges = %d\tNumber of Boundary Faces = %d\n",nbvert,nbedge,nbface);
  printf("Order of Quadrature = %d\n",2*nq1d-1);
  printf("************************************************************************************************\n\n");
  fflush(stdout);

  /*** Assemble the Linear matrix and right hand side *******************************/
	
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
  //  A =   0        -K^T Mf          Med G
  //      Mf K         0                0
  //    -G^T Med       0                0
  //
  //  M ->    <eps dE/dt,F>
  //         <(1/mu) dB/dt,C>
  //            <dp/dt,q>
  //
  //    Z ->  <n x E, n x F>_bdry_i
  //               0
  //               0
  //
  //    f ->     -<j,F>
  //               0
  //               0
	
  printf("Assembling Linearized System:\n");
  fflush(stdout);
  clk0 = clock();
	
  // Number of degrees of freedom equals number of nodes for p plus edges for E plus faces for B
  ndof = n + nedge + nface;
  
  // Build individual Mass matrices
  // TODO: Have confirmed that Me is the same as in the giant assembly set up, but Mv and Mf are off only at a few DOF.  Might be counting thing.
  Me.row = nedge;
  Me.col = nedge;
  Mv.row = n;
  Mv.col = n;
  Mf.row = nface;
  Mf.col = nface;
  Me.IA = (INT *) calloc(nedge+1,sizeof(INT));
  Mf.IA = (INT *) calloc(nface+1,sizeof(INT));
  Mv.IA = (INT *) calloc(n+1,sizeof(INT));
  stiffG_nnz(Me.IA,&Me.nnz,el_ed.IA,el_ed.JA,nedge,nelm,isdirichletE);
  stiffG_nnz(Mf.IA,&Mf.nnz,el_face.IA,el_face.JA,nface,nelm,isdirichletF);
  stiffG_nnz(Mv.IA,&Mv.nnz,el_n.IA,el_n.JA,n,nelm,isdirichletV);
  Me.JA = (INT *) calloc(Me.nnz,sizeof(INT));
  Mf.JA = (INT *) calloc(Mf.nnz,sizeof(INT));
  Mv.JA = (INT *) calloc(Mv.nnz,sizeof(INT));
  Me.val = (REAL *) calloc(Me.nnz,sizeof(REAL));
  Mf.val = (REAL *) calloc(Mf.nnz,sizeof(REAL));
  Mv.val = (REAL *) calloc(Mv.nnz,sizeof(REAL));
  for(i=0;i<Me.nnz;i++) { Me.val[i] = 0.0; }
  for(i=0;i<Mf.nnz;i++) { Mf.val[i] = 0.0; }
  for(i=0;i<Mv.nnz;i++) { Mv.val[i] = 0.0; }
  stiffG_cols(Me.JA,el_ed.IA,el_ed.JA,nedge,nelm,isdirichletE);
  stiffG_cols(Mf.JA,el_face.IA,el_face.JA,nface,nelm,isdirichletF);
  stiffG_cols(Mv.JA,el_n.IA,el_n.JA,n,nelm,isdirichletV);

  mass_assemble(Me.IA,Me.JA,Me.val,cn,cq,el_n,el_ed,el_face,ed_n,edge_order,element_order,nq1d,myeps,mydim,nelm,nedge,isdirichletE,1,NULL,NULL);
  mass_assemble(Mf.IA,Mf.JA,Mf.val,cn,cq,el_n,el_ed,el_face,face_n,face_order,element_order,nq1d,oneovermu,mydim,nelm,nface,isdirichletF,2,f_area,f_norm);
  mass_assemble(Mv.IA,Mv.JA,Mv.val,cn,cq,el_n,el_n,el_face,ed_n,element_order,element_order,nq1d,1.0,mydim,nelm,n,isdirichletV,0,NULL,NULL);
    
  // Get gradient + curl: G seems to be correct since G'MeG gives laplacian.  K'MfK doesn't give curl curl so, something is off...
  get_grad_H1toNed(&G,ed_len,ed_n,nedge);
  get_curl_NedtoRT(&K,f_area,f_norm,ed_len,cn,mydim,face_n,f_ed,ed_n,nface,ndpf);
  // Because p is in H^1_0 we annihilate the gradient at the boundary node
  for(i=0;i<G.nnz;i++) {
    if(n_bdry[G.JA[i]-1]!=0) {
      G.val[i] = 0.0;
    }
  }
  // Because E is in H_0(curl) on outer boundary, we zero those edges out.
  for(i=0;i<K.nnz;i++) {
    if(ed_bdry[K.JA[i]-1]!=0) {
      K.val[i] = 0.0;
    }
  }

  // Build Z - this matches other code
  Z.row = nedge;
  Z.col = nedge;
  Z.IA = (INT *) calloc(nedge+1,sizeof(INT));
  stiffGsubset_nnz(Z.IA,&Z.nnz,el_ed.IA,el_ed.JA,nedge,nelm,dof_bdry,-1);
  Z.JA = (INT *) calloc(Z.nnz,sizeof(INT));
  Z.val = (REAL *) calloc(Z.nnz,sizeof(REAL));
  stiffGsubset_cols(Z.JA,el_ed.IA,el_ed.JA,nedge,nelm,dof_bdry,-1);
  get_impedenceZ(Z.IA,Z.JA,Z.val,cn,cq,nq1d,el_n,element_order,el_ed,ed_n,edge_order,f_ed,face_order,el_face,edpf,ndpf,f_norm,f_area,ed_mid,mydim, \
  		 nelm,n,nedge,nface,dof_bdry);
  // Multiply Z by (1+gamma(t))
  if(gamistimedep) {
    getgamma(&gam,0.0);
  }
  Ztmp = (REAL *) calloc(Z.nnz,sizeof(REAL));
  for(i=0;i<Z.nnz;i++) {
    Ztmp[i] = Z.val[i];
    Z.val[i] = (1.0+gam)*Z.val[i];
  }

  if (dumpmat>=2) {
    
    FILE* meid = fopen("output/me.dat","w");
    FILE* mfid = fopen("output/mf.dat","w");
    FILE* zid = fopen("output/z.dat","w");
    FILE* mvid = fopen("output/mv.dat","w");
    FILE* gid = fopen("output/grad.dat","w");
    FILE* kid = fopen("output/curl.dat","w");
    dumpmatrices(Me.IA,Me.JA,Me.val,NULL,meid,NULL,nedge,dumpmat);
    dumpmatrices(Mf.IA,Mf.JA,Mf.val,NULL,mfid,NULL,nface,dumpmat);
    dumpmatrices(Mv.IA,Mv.JA,Mv.val,NULL,mvid,NULL,n,dumpmat);
    dumpmatrices(K.IA,K.JA,K.val,NULL,kid,NULL,nface,dumpmat);
    dumpmatrices(G.IA,G.JA,G.val,NULL,gid,NULL,nedge,dumpmat);
    dumpmatrices(Z.IA,Z.JA,Z.val,NULL,zid,NULL,nedge,dumpmat);
    fclose(meid);
    fclose(mvid);
    fclose(mfid);
    fclose(kid);
    fclose(gid);
    fclose(zid);
  }

  // Get operators for fast algorithm:
  // (dt/2 G'MeG + 2/dt Mv) p_{n+1} = 2 G'Me E_n + (-dt/2 G'MeG + 2/dt Mv) p_n
  // (dt/2 K'MfK + 2/dt Me + Z) E_{n+1} = 2 K'Mf B_n + (2/dt Me - dt/2 K'MfK - Z)E_n - MeG(p_n + p_{n+1})
  // 2/dt Mf B_{n+1} = 2/dt Mf B_n - MfK(E_n + E_{n+1})

  // MeG and G'Me = MGt
  abybCSR(Me,G,&MG);
  atranspCSR(MG,&MGt);
  // MfK and K'Mf = MKt
  abybCSR(Mf,K,&MK);
  atranspCSR(MK,&MKt);
  // G'MeG and K'MfK
  abybCSR(MGt,G,&GMG); /* This is correct, except 0 on the boundaries instead of identity */
  abybCSR(MKt,K,&KMK); /* XXX Not matching curl-curl code yet...*/
  
  //Test for Comparison build laplacian and curl curl
  dCSRmat Lap;
  dCSRmat CurlCurl;
  Lap.row = n;
  Lap.col = n;
  CurlCurl.row = nedge;
  CurlCurl.col = nedge;
  f = calloc(n,sizeof(REAL));
  REAL *fC = calloc(nedge,sizeof(REAL));
  Lap.IA = calloc(n+1,sizeof(INT));
  stiffG_nnz(Lap.IA,&Lap.nnz,el_n.IA,el_n.JA,n,nelm,isdirichletV);
  Lap.val = calloc(Lap.nnz,sizeof(REAL));
  for (i=0; i<Lap.nnz; i++) {
    Lap.val[i] = 0;
  }
  Lap.JA = calloc(Lap.nnz,sizeof(INT));

  H1_assemble_Lap(Lap.IA,Lap.JA,Lap.val,f,cn.x,cn.y,cn.z,cq.x,cq.y,cq.z,cq.w,el_n.IA,el_n.JA,element_order,nq1d,1.0,mydim,nelm,n,n_bdry,0,1,0,NULL,NULL,NULL,0.0);

  CurlCurl.IA = calloc(nedge+1,sizeof(INT));
  stiffG_nnz(CurlCurl.IA,&CurlCurl.nnz,el_ed.IA,el_ed.JA,nedge,nelm,isdirichletE);
  CurlCurl.val = calloc(CurlCurl.nnz,sizeof(REAL));
  for (i=0; i<CurlCurl.nnz; i++) {
    CurlCurl.val[i] = 0;
  }
  CurlCurl.JA = calloc(CurlCurl.nnz,sizeof(INT));
   
  ned_assembleG(CurlCurl.IA,CurlCurl.JA,CurlCurl.val,fC,cn.x,cn.y,cn.z,cq.x,cq.y,cq.z,cq.w,el_n.IA,el_n.JA,el_ed.IA,el_ed.JA,ed_n.IA,ed_n.JA,edge_order, \
  		element_order,nq1d,1.0,0.0,mydim,nelm,nedge,isdirichletE,0,0);

  printf("hello\n\n\n");
  FILE* mylap;
  FILE* gmgid;
  FILE* mycurl;
  FILE* kmkid;
  mylap = fopen("output/lap.dat","w");
  gmgid = fopen("output/gmg.dat","w");
  mycurl = fopen("output/curlcurl.dat","w");
  kmkid = fopen("output/kmk.dat","w");
  dumpmatrices(Lap.IA,Lap.JA,Lap.val,NULL,mylap,NULL,n,dumpmat);
  dumpmatrices(CurlCurl.IA,CurlCurl.JA,CurlCurl.val,NULL,mycurl,NULL,nedge,dumpmat);
  dumpmatrices(GMG.IA,GMG.JA,GMG.val,NULL,gmgid,NULL,n,dumpmat);
  dumpmatrices(KMK.IA,KMK.JA,KMK.val,NULL,kmkid,NULL,nedge,dumpmat);
  fclose(mylap);
  fclose(gmgid);
  fclose(kmkid);
  fclose(mycurl);
  freedCSRmat(Lap);
  freedCSRmat(CurlCurl);
  
  /* // Now combine to get operators: */
  /* // (dt/2 G'MeG + 2/dt Mv) p_{n+1} = 2 G'Me E_n + (-dt/2 G'MeG + 2/dt Mv) p_n */
  /* dCSRmat pop; */
  /* dCSRmat popR; */
  /* REAL* pr1=NULL; */
  /* REAL* pr2=NULL; */
  /* REAL* pavg=NULL; */
  /* REAL* Er1=NULL; */
  /* REAL* Er2=NULL; */
  /* REAL* Er3=NULL; */
  /* REAL* utmp=NULL; */
  /* REAL* uptmp=NULL; */
  /* REAL* ftmp=NULL; */
  /* aplusbCSR(GMG,Mv,&pop,0.5*dt,alpha); */
  /* aplusbCSR(GMG,Mv,&popR,-0.5*dt,alpha); */
  /* // (dt/2 K'MfK + 2/dt Me + Z) E_{n+1} = 2 K'Mf B_n + (2/dt Me - dt/2 K'MfK - Z)E_n - MeG(p_n + p_{n+1}) */
  /* dCSRmat Eop1; */
  /* dCSRmat EopR1; */
  /* dCSRmat Eop; */
  /* dCSRmat EopR; */
  /* aplusbCSR(KMK,Me,&Eop1,0.5*dt,alpha); */
  /* aplusbCSR(Eop1,Z,&Eop,1.0,1.0); */
  /* /\* aplusbCSR(Lap,Z,&Eop,1.0,1.0); *\/ */
  /* aplusbCSR(KMK,Me,&EopR1,-0.5*dt,alpha); */
  /* aplusbCSR(EopR1,Z,&EopR,1.0,-1.0); */
  /* /\* aplusbCSR(LapR,Z,&EopR,1.0,-1.0); *\/ */
  /* // 2/dt Mf B_{n+1} = 2/dt Mf B_n - MfK(E_n + E_{n+1}) */
  /* REAL* Br1=NULL; */
  /* REAL* Br2=NULL; */
  /* REAL* Eavg=NULL; */
  
  /* clk1 = clock(); */
  /* printf("Elapsed CPU Time for Assembly of Linear Systems = %f seconds.\n\n",(REAL) (clk1-clk0)/CLOCKS_PER_SEC); */

  /* /\*******************************************************************************************\/ */

  /* /\**************** Initialize Timestepping if Necessary ********************************************************************\/ */
	
  /* // Get Initial Conditions First */
  /* u = calloc(ndof,sizeof(REAL)); */
  /* uprev = calloc(ndof,sizeof(REAL)); */
  /* clk0 = clock(); */

  /* // Open files for dumping if necessary */
  /* if (havetrue) { utrue = calloc(ndof,sizeof(REAL)); } */

  /* if (dumpsol>=2) { */
  /*   uid = fopen("output/sol.dat","w"); */
  /*   if(havetrue) { truid = fopen("output/true.dat","w"); } */
  /* } */

  /* if (istimdep) {  // If Time stepping is activated initialize solutions */
		
  /*   printf("The system is time-dependent.  Perform %d Time Step(s)\n\n",nsteps); */
		
  /*   // Initialize solution to be initial condition */
  /*   for (i=0; i<nedge; i++) { // tau*E */
  /*     if(ed_bdry[i]==1) { */
  /* 	u[i] = 0.0; */
  /* 	uprev[i] = 0.0; */
  /*     } else { */
  /* 	get_taudot(&uval,ed_tau,E1init,E2init,E3init,ed_mid,ed_len[i],cn,ed_n,mydim,0.0,gam,i+1,-1); */
  /* 	u[i]=uval; */
  /* 	uprev[i] = uval; */
  /*     } */
  /*     if (havetrue) { */
  /* 	if(ed_bdry[i]==1) { */
  /* 	  utrue[i] = 0.0; */
  /* 	} else { */
  /* 	  get_taudot(&uval,ed_tau,E1t,E2t,E3t,ed_mid,ed_len[i],cn,ed_n,mydim,0.0,gam,i+1,-1); */
  /* 	  utrue[i] = uval; */
  /* 	} */
  /*     } */
  /*   } */

  /*   // Clean up divergence of E */
  /*   cleandivE(u,G,n,nedge,nelm,el_n,el_ed,ed_n,cn,cq,nq1d,element_order,edge_order,mydim,n_bdry,ed_len,ed_bdry,0); */

  /*   // Clean up harmonic component of E */
  /*   //cleanharmonic(u,Ggrad,Gulap,&Gulapnorm,n,nedge,nelm,el_n,el_ed,ed_n,cn,cq,nq1d,element_order,edge_order,mydim,n_bdry,ed_len,ed_bdry,1); */

  /*   for(i=0;i<nedge;i++) { uprev[i] = u[i]; } */

  /*   myr = fabs(1.0/myr); */
  /*   myB = (REAL *) calloc(nface,sizeof(REAL)); */
  /*   if(havetrue) myBt = (REAL *) calloc(nface,sizeof(REAL)); */
  /*   ax0(&nface,K.IA,K.JA,K.val,u,myB); */
  /*   if(havetrue) ax0(&nface,K.IA,K.JA,K.val,utrue,myBt); */
  /*   for(i=0;i<nface;i++) { // n*B */
  /*     j = i+nedge; */
  /*     u[j] = myr*myB[i]; */
  /*     uprev[j] = u[j]; */
  /*     if(havetrue) { */
  /* 	utrue[j] = myr*myBt[i]; */
  /*     } */
  /*   } */
  /*   if(myB) free(myB); */
  /*   if(myBt) free(myBt); */

  /*   for (i=0; i<n; i++) { // p */
  /*     j = i+nedge+nface; */
  /*     u[j] = 0.0; */
  /*     uprev[j] = 0.0; */
  /*     if (havetrue) { */
  /* 	utrue[j] = 0.0; */
  /*     } */
  /*   } */
  /* } else { // else no time stepping */
  /*   nsteps = 1; */
  /*   printf("How can you not have any timestepping????\n"); */
  /*   exit(0); */
  /* } */

  /* // Remove Dirichlet DOFS */
  /* clk0 = clock(); */
  /* INT removedof = 0; */
  /* if(removedof) { */
  /*   fprintf(stdout,"\nRemoving Dirichlet dofs:\n "); */
  /*   shrinkDOFv(u,&ndof,&ndofnew,isdirichlet); */
  /*   shrinkDOFv(uprev,&ndof,&ndofnew,isdirichlet); */
  /*   shrinkDOFv(Ztmp,&nedge,&nedgenew,isdirichletE); */
  /*   shrinkDOFm(&n,&n,&nnew,&nnew,pop,&pop,isdirichletV,isdirichletV); */
  /*   shrinkDOFm(&n,&n,&nnew,&nnew,popR,&popR,isdirichletV,isdirichletV); */
  /*   shrinkDOFm(&n,&nedge,&nnew,&nedgenew,MGt,&MGt,isdirichletV,isdirichletE); */
  /*   shrinkDOFm(&nedge,&nedge,&nedgenew,&nedgenew,Eop,&Eop,isdirichletE,isdirichletE); */
  /*   shrinkDOFm(&nedge,&nedge,&nedgenew,&nedgenew,Eop1,&Eop1,isdirichletE,isdirichletE); */
  /*   shrinkDOFm(&nedge,&nedge,&nedgenew,&nedgenew,EopR,&EopR,isdirichletE,isdirichletE); */
  /*   shrinkDOFm(&nedge,&nedge,&nedgenew,&nedgenew,EopR1,&EopR1,isdirichletE,isdirichletE); */
  /*   shrinkDOFm(&nedge,&nedge,&nedgenew,&nedgenew,Z,&Z,isdirichletE,isdirichletE); */
  /*   shrinkDOFm(&nedge,&nface,&nedgenew,&nfacenew,MKt,&MKt,isdirichletE,isdirichletF); */
  /*   shrinkDOFm(&nedge,&n,&nedgenew,&nnew,MG,&MG,isdirichletE,isdirichletV); */
  /*   shrinkDOFm(&nface,&nface,&nfacenew,&nfacenew,Mf,&Mf,isdirichletF,isdirichletF); */
  /*   shrinkDOFm(&nface,&nedge,&nfacenew,&nedgenew,MK,&MK,isdirichletF,isdirichletE); */

  /*   dofbdrynew = calloc(ndofnew,sizeof(INT)); */
  /*   for(i=0;i<ndofnew;i++) { */
  /*     dofbdrynew[i] = 0; */
  /*   } */
  /* } else { */
  /*   nedgenew = nedge; */
  /*   nfacenew = nface; */
  /*   nnew = n; */
  /*   ndofnew = ndof; */
  /* } */
  /* printf("Old DOF: Tot: %d\tEdge: %d\tFace: %d\tVertex: %d\n New DOF: Tot: %d\tEdge: %d\tFace: %d\tVertex: %d\n",ndof,nedge,nface,n,ndofnew,nedgenew,nfacenew,nnew); */
  /* clk1 = clock(); */
  /* printf("Elapsed CPU Time for removing the Dirichlet dofs = %f seconds.\n\n",(REAL) (clk1-clk0)/CLOCKS_PER_SEC); */

  /* // Allocate array for updated right-hand side */
  /* f = (REAL *) calloc(ndofnew,sizeof(REAL)); */
  /* uzero = (REAL *) calloc(ndofnew,sizeof(REAL)); */
  /* pr1 = (REAL *) calloc(nnew,sizeof(REAL)); */
  /* pr2 = (REAL *) calloc(nnew,sizeof(REAL)); */
  /* pavg = (REAL *) calloc(nnew,sizeof(REAL)); */
  /* Er1 = (REAL *) calloc(nedgenew,sizeof(REAL)); */
  /* Er2 = (REAL *) calloc(nedgenew,sizeof(REAL)); */
  /* Er3 = (REAL *) calloc(nedgenew,sizeof(REAL)); */
  /* Eavg = (REAL *) calloc(nedgenew,sizeof(REAL)); */
  /* Br1 = (REAL *) calloc(nfacenew,sizeof(REAL)); */
  /* Br2 = (REAL *) calloc(nfacenew,sizeof(REAL)); */
  /* for (i=0; i<ndofnew; i++) { */
  /*   f[i] = 0.0; */
  /*   uzero[i] = 0.0; */
  /* } */
  /* /\*******************************************************************************************\/ */

  /* /\***************** Call Arpack ***********************************************************\/ */
  /* /\* if(dumpsol>=2 && dumpmat >=2) return 2; *\/ */

  /* /\* fprintf(stdout,"\nStarting ARPACK: "); *\/ */
  /* /\* fflush(stdout);	 *\/ */

  /* /\* clk2 = clock(); *\/ */
  /* /\* callarpack(ndofnew, iAold, jAold, Aold,	\ *\/ */
  /* /\* 	     iMold,jMold, Mold,			\ *\/ */
  /* /\* 	     iZold,jZold, Zold);  *\/ */
  /* /\* clk3 = clock(); *\/ */
  /* /\* printf("Elapsed CPU Time for ARPACK = %f seconds.\n\n",(REAL) (clk3-clk2)/CLOCKS_PER_SEC); *\/ */
  /* /\* fflush(stdout);	 *\/ */
  /* /\* return 3; *\/ */
  /* /\*******************************************************************************************\/ */

  /* /\**************** Begin Timestepping if Necessary *************************************\/ */
  /* char mystring[100]; */
  /* // Begin Time Stepping Loop (make sure if turned off nsteps = 1!!) */
  /* if (solvertype==0) { */
  /*   strcpy(mystring,"Solving the System: Using Conjugate Gradient Method\n\n"); */
  /* } else if (solvertype==1) { */
  /*   strcpy(mystring,"Solving the System: Using the Minimum Residual Method (MINRES)\n\n"); */
  /* } else if (solvertype==2) { */
  /*   strcpy(mystring,"Solving the System: Using the Generalized Minimum Residual Method (GMRES)\n\n"); */
  /* } else if (solvertype==3) { */
  /*   strcpy(mystring,"Solving the System: Using Preconditioned MINRES (M = diag((Lap)^(-1) (Mass_P0)^(-1))\n\n"); */
  /* } else if (solvertype==4) { */
  /*   strcpy(mystring,"Solving the System: UMFPACK\n\n"); */
  /* } else if (solvertype==5) { */
  /*   strcpy(mystring,"Solving the System in 3 parts using Fast Algorithm\n\n"); */
  /* } else { */
  /*   strcpy(mystring,"Not sure what kind of solver you're using.  Good luck!\n\n"); */
  /* } */

  /* if(restart==-1) { restart = ndofnew; } */

  /* if(dumpsol==2) { */
  /*   texout = fopen("output/texout.dat","w"); */
  /*   nmout = fopen("output/maxdisp_norms.out","w"); */
  /* } */

  /* if (istimdep) { // If time dependent */
  /*   clk0 = clock(); */
  /*   // Initial Condition */
  /*   printf("\n********************************************************************************************************\n"); */
  /*   printf("Initial Condition:\t Actual Time = 0.0\n"); */
  /*   printf("********************************************************************************************************\n\n"); */
  /*   if (dumpsol==2) { */
  /*     if (havetrue) { */
  /* 	fflush(stdout); */
  /* 	REAL* unew = (REAL *) calloc(ndof,sizeof(REAL)); */
  /* 	for(j=0;j<ndof;j++) { unew[j] = u[j]; } */
  /* 	shrink0(&ndof,&ndofnew,&nedge,&nface,unew,dof_bdry); */
  /* 	gettrue_max(unew,texout,nmout,uid,truid,ndof,ndofnew,n,nedge,nface,nelm,mydim,cn,cq.x,cq.y,cq.z,cq.w,nq1d, \ */
  /* 		    el_n.IA,el_n.JA,element_order,nve,face_order,	\ */
  /* 		    el_face.IA,el_face.JA,face_n.IA,face_n.JA,edge_order,el_ed.IA,el_ed.JA,ed_n,0.0,0,gam,E1t,E2t,E3t,B1t,B2t,B3t,f_norm, \ */
  /* 		    f_area,f_mid,ed_tau,ed_len,ed_mid,dumpsol,dof_bdry); */
  /* 	fflush(stdout); */
  /* 		if(unew) free(unew); */
  /*     } else { */
  /* 	dumpsolution(u,utrue,uid,truid,ndofnew,n,nelm,2,2,mydim,havetrue); */
  /*     } */
  /*   } */
  /*   printf("\n"); */
  /*   for (i=0; i<nsteps; i++) { */
  /*     clk2 = clock(); */
  /*     time = (i+1)*dt; */
  /*     printf("********************************************************************************************************\n"); */
  /*     printf("TIME STEP: %d\t Actual Time = %f\n",i+1,time); */
  /*     printf("********************************************************************************************************\n\n"); */

  /*     // Clean Divergence of E */
  /*     cleandivE(u,G,n,nedge,nelm,el_n,el_ed,ed_n,cn,cq,nq1d,element_order,edge_order,mydim,n_bdry,ed_len,ed_bdry,removedof); */
  /*     for(j=0;j<nedgenew;j++) { */
  /* 	uprev[j] = u[j]; */
  /*     } */
      
  /*     // Fix (1+gamma(t))*Z and fix operators */
  /*     if(gamistimedep) { */
  /* 	getgamma(&gam,time); */
  /* 	for(i=0;i<Z.IA[nedgenew]-1;i++) { */
  /* 	  Z.val[i] = (1.0+gam)*Ztmp[i]; */
  /* 	} */
  /* 	aplusbCSR(Eop1,Z,&Eop,1.0,1.0); */
  /* 	aplusbCSR(EopR1,Z,&EopR,1.0,-1.0); */
  /*     } */
      
  /*     // Solve */
  /*     // p:  */
  /*     // (dt/2 G'MeG + 2/dt Mv) p_{n+1} = 2 G'Me E_n + (-dt/2 G'MeG + 2/dt Mv) p_n */
  /*     printf("Solving for p\n"); */
  /*     ax0(&nnew,MGt.IA,MGt.JA,MGt.val,uprev,pr1); */
  /*     uptmp = uprev + nedgenew + nfacenew; */
  /*     ax0(&nnew,popR.IA,popR.JA,popR.val,uptmp,pr2); */
  /*     for(j=0;j<nnew;j++) { */
  /* 	f[j+nedgenew+nfacenew] = 2*pr1[j] + pr2[j]; */
  /*     } */
     
  /*     // For now use CG: */
  /*     utmp = u + nedgenew + nfacenew; */
  /*     ftmp = f + nedgenew + nfacenew; */
  /*     cg(pop.IA,pop.JA,pop.val,ftmp,nnew,utmp,solvetol,maxit,outprint); */

  /*     // E: */
  /*     // (dt/2 K'MfK + 2/dt Me + Z) E_{n+1} = 2 K'Mf B_n + (2/dt Me - dt/2 K'MfK - Z)E_n - MeG(p_n + p_{n+1}) */
  /*     printf("Solving for E\n"); */
  /*     utmp = uprev + nedgenew; */
  /*     ax0(&nedgenew,MKt.IA,MKt.JA,MKt.val,utmp,Er1); */
  /*     ax0(&nedgenew,EopR.IA,EopR.JA,EopR.val,uprev,Er2); */
  /*     for(j=0;j<nnew;j++) { */
  /* 	pavg[j] = uprev[nedgenew+nfacenew+j] + u[nedgenew+nfacenew+j]; */
  /*     } */
  /*     ax0(&nedgenew,MG.IA,MG.JA,MG.val,pavg,Er3); */
  /*     for(j=0;j<nedgenew;j++) { */
  /* 	f[j] = 2*Er1[j] + Er2[j] - Er3[j]; */
  /*     } */
  /*     // For now use CG: */
  /*     cg(Eop.IA,Eop.JA,Eop.val,f,nedgenew,u,solvetol,maxit,outprint); */

  /*     // B: */
  /*     // 2/dt Mf B_{n+1} = 2/dt Mf B_n - MfK(E_n + E_{n+1}) (Divide by alpha for ease of implementation) */
  /*     printf("Solving for B\n"); */
  /*     uptmp = uprev + nedgenew; */
  /*     ax0(&nfacenew,Mf.IA,Mf.JA,Mf.val,uptmp,Br1); */
  /*     for(j=0;j<nedgenew;j++) { */
  /* 	Eavg[j] = uprev[j] + u[j]; */
  /*     } */
  /*     ax0(&nfacenew,MK.IA,MK.JA,MK.val,Eavg,Br2); */
  /*     for(j=0;j<nfacenew;j++) { */
  /* 	f[j+nedgenew] = Br1[j] - 0.5*dt*Br2[j]; */
  /*     } */
  /*     // For now use CG: */
  /*     ftmp = f + nedgenew; */
  /*     utmp = u + nedgenew; */
  /*     cg(Mf.IA,Mf.JA,Mf.val,ftmp,nedgenew,utmp,solvetol,maxit,outprint); */
      
  /*     // Update Time Step */
  /*     for (j=0; j<ndofnew; j++) { */
  /* 	uprev[j] = u[j]; */
  /*     } */
  /*     // Dump Solution and get errors if true solution available */
  /*     if (dumpsol==2) { */
  /* 	if (havetrue) { */
  /* 	  fflush(stdout); */
  /* 	  REAL* unew = (REAL *) calloc(ndof,sizeof(REAL)); */
  /* 	  for(j=0;j<ndof;j++) { unew[j] = u[j]; } */
  /* 	  shrink0(&ndof,&ndofnew,&nedge,&nface,unew,dof_bdry); */
  /* 	  gettrue_max(unew,texout,nmout,uid,truid,ndof,ndofnew,n,nedge,nface,nelm,mydim,cn,cq.x,cq.y,cq.z,cq.w,nq1d,el_n.IA,el_n.JA,element_order,nve,face_order, \ */
  /* 		      el_face.IA,el_face.JA,face_n.IA,face_n.JA,edge_order,el_ed.IA,el_ed.JA,ed_n,time,i+1,gam,E1t,E2t,E3t,B1t,B2t,B3t,f_norm, \ */
  /* 		      f_area,f_mid,ed_tau,ed_len,ed_mid,dumpsol,dof_bdry); */
  /* 	  fflush(stdout); */
  /* 	  if(unew) free(unew); */
  /* 	} else { */
  /* 	  dumpsolution(u,utrue,uid,truid,ndofnew,n,nelm,2,2,mydim,havetrue); */
  /* 	} */
  /*     } */
  /*     clk3 = clock(); */
  /*     printf("Elapsed CPU Time for Time Step %d = %f seconds.\n\n\n",i+1,(REAL) (clk3-clk2)/CLOCKS_PER_SEC); */
  /*   } // End of Timestepping Loop */
  /* } else { // if not time dependent */
  /*   printf("How can you not have any timestepping????\n"); */
  /*   exit(0); */
  /* } */
  /* fclose(texout); */
  /* fclose(nmout); */
  /* clk1 = clock(); */
  /* printf("Elapsed CPU Time for Solve = %f seconds.\n\n",(REAL) (clk1-clk0)/CLOCKS_PER_SEC); */
  /* /\*******************************************************************************************\/ */
  /* /\******** Free All the Arrays ***************************************************************\/ */
  
  /* freedCSRmat(Z); */
  /* freedCSRmat(K); */
  /* freedCSRmat(G); */
  /* freedCSRmat(Me); */
  /* freedCSRmat(Mf); */
  /* freedCSRmat(Mv); */
  /* freedCSRmat(MG); */
  /* freedCSRmat(MK); */
  /* freedCSRmat(MGt); */
  /* freedCSRmat(MKt); */
  /* freedCSRmat(GMG); */
  /* freedCSRmat(KMK); */
  /* freedCSRmat(pop); */
  /* freedCSRmat(popR); */
  /* freedCSRmat(Eop); */
  /* freedCSRmat(EopR); */
  /* if(f) free(f); */
  /* if(u) free(u); */
  /* if(pr1) free(pr1); */
  /* if(pr2) free(pr2); */
  /* if(pavg) free(pavg); */
  /* if(Er1) free(Er1); */
  /* if(Er2) free(Er2); */
  /* if(Er3) free(Er3); */
  /* if(Eavg) free(Eavg); */
  /* if(Br1) free(Br1); */
  /* if(Br2) free(Br2); */
  /* if(uprev) free(uprev); */
  /* if(uzero) free(uzero); */
  /* freecoords(cn); */
  /* freecoords(cv); */
  /* freeCSRinc(el_n); */
  /* freeCSRinc(ed_n); */
  /* freeCSRinc(el_ed); */
  /* freeCSRinc(el_v); */
  /* freeiCSRmat(el_face); */
  /* freeCSRinc(face_n); */
  /* freeCSRinc(f_ed); */
  /* if(ed_len) free(ed_len); */
  /* if(ed_tau) free(ed_tau); */
  /* if(ed_mid) free(ed_mid); */
  /* if(fel_order) free(fel_order); */
  /* if(f_area) free(f_area); */
  /* if(f_norm) free(f_norm); */
  /* if(f_mid) free(f_mid); */
  /* freeqcoords(cq); */
  /* if(ed_bdry) free(ed_bdry); */
  /* if(n_bdry) free(n_bdry); */
  /* if(v_bdry) free(v_bdry); */
  /* if(dof_bdry) free(dof_bdry); */
  /* if(face_bdry) free(face_bdry); */
  /* if(dofbdrynew) free(dofbdrynew); */
  /* if(havetrue) { if(utrue) free(utrue); } */
  /* if (dumpsol>=2) { */
  /*   fclose(uid); */
  /*   if(havetrue) { fclose(truid); } */
  /* } */
  /*******************************************************************************************/
	
  clkb = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",(REAL) (clkb-clka)/CLOCKS_PER_SEC);
  return 0;	
	
}	/* End of Program */
/*******************************************************************************************/
