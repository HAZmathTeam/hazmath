/*
 *  assemble_local.c   
 *  
 *  Created by James Adler and Xiaozhe Hu on 4/22/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

/* This code will build local stiffness matrices for various PDE systems */

// Standard Includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// Our Includes
#include "macro.h"
#include "grid.h"
#include "sparse.h"
#include "vec.h"
#include "functs.h"
#include "fem.h"

/******************************************************************************************************/
void assemble_DuDv_local(REAL* ALoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time) 
{
	
  /* Computes the local stiffness matrix for coeff*<Du,Dv> = <f,v> bilinear form using various element types
   * (eg. P1, P2 -> (grad u, grad v), Nedelec <curl u, curl v>, and Raviart-Thomas <div u, div v>).
   * 
   * For this problem we compute:
   *
   * coeff*D*D u = f  ---->   coeff*<D u, D v> = <f,v>
   *
   * which gives Ax = b,
   *
   * A_ij = coeff*<D phi_j,D phi_i>
   * b_i  = <f,phi_i>
   * 
   *    INPUT:
   *            fespace		    Finite-Element Space Struct
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights
   *            dof_on_elm          Specific DOF on element
   *            elm                 Current element
   *            coeff               Function that gives coefficient (for now assume constant)
   *            time                Time if time dependent
   *
   *    OUTPUT:
   *	       	ALoc                Local Stiffness Matrix (Full Matrix)
   */
	
  // Mesh and FE data
  INT dof_per_elm = FE->dof_per_elm;
  INT dim = mesh->dim;
  
  // Loop Indices
  INT quad,test,trial;
  // Quadrature Weights and Nodes
  REAL w; 
  REAL* qx = (REAL *) calloc(3,sizeof(REAL));
  // Stiffness Matrix Entry
  REAL kij;
  // Coefficient Value at Quadrature Nodes
  REAL coeff_val=0.0;
	
  // Basis Functions and its derivatives if necessary
  REAL* phi=NULL;
  REAL* dphix=NULL;
  REAL* dphiy=NULL;
  REAL* dphiz=NULL;

  if(FE->FEtype>0) { // PX elements
    phi = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    dphix = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    dphiy = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    if(dim==3) dphiz = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    
    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {        
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      (*coeff)(&coeff_val,qx,time);
      
      //  Get the Basis Functions at each quadrature node
      PX_H1_basis(phi,dphix,dphiy,dphiz,qx[0],qx[1],qx[2],dof_on_elm,FE->FEtype,mesh);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->dof_per_elm; trial++) {
	  kij = coeff_val*(dphix[test]*dphix[trial] + dphiy[test]*dphiy[trial]);
	  if(dim==3) kij+=dphiz[test]*dphiz[trial];
	  ALoc[test*FE->dof_per_elm+trial] = ALoc[test*FE->dof_per_elm+trial] + w*kij;
	}
      }
    }
  } else if(FE->FEtype==-1) { // Nedelec elements
    phi = (REAL *) calloc(dof_per_elm*dim,sizeof(REAL));
    if(dim==2) {
      dphix = (REAL *) calloc(dof_per_elm,sizeof(REAL)); // Curl of basis function
    } else if (dim==3) {
      dphix = (REAL *) calloc(dof_per_elm*dim,sizeof(REAL)); // Curl of basis function
    } else {
      baddimension();
    }

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {        
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      (*coeff)(&coeff_val,qx,time);

      //  Get the Basis Functions at each quadrature node
      ned_basis(phi,dphix,qx[0],qx[1],qx[2],v_on_elm,dof_on_elm,mesh);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->dof_per_elm; trial++) {
	  if(dim==2) {
	    kij = coeff_val*(dphix[test]*dphix[trial]);
	  } else if (dim==3) {
	    kij = coeff_val*(dphix[test*dim]*dphix[trial*dim] + dphix[test*dim+1]*dphix[trial*dim+1] + dphix[test*dim+2]*dphix[trial*dim+2]);
	  } else {
	    baddimension();
	  }
	  ALoc[test*FE->dof_per_elm+trial] = ALoc[test*FE->dof_per_elm+trial] + w*kij;
	}
      }
    }
  } else if(FE->FEtype==-2) { // Raviart-Thomas elements
    phi = (REAL *) calloc(dof_per_elm*dim,sizeof(REAL));
    dphix = (REAL *) calloc(dof_per_elm,sizeof(REAL)); // Divergence of element 

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {        
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      (*coeff)(&coeff_val,qx,time);

      //  Get the Basis Functions at each quadrature node
      rt_basis(phi,dphix,qx[0],qx[1],qx[2],v_on_elm,dof_on_elm,mesh);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->dof_per_elm; trial++) {
	  kij = coeff_val*(dphix[test]*dphix[trial]);
	  ALoc[test*FE->dof_per_elm+trial] = ALoc[test*FE->dof_per_elm+trial] + w*kij;
	}
      }
    }
  } else {
    printf("Trying to implement non-existent elements!\n\n");
    exit(0);
  }
  
  if (phi) free(phi);
  if(dphix) free(dphix);
  if(dphiy) free(dphiy);
  if(dphiz) free(dphiz);
  if(qx) free(qx);
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void assemble_mass_local(REAL* MLoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time) 
{
	
  /* Computes the global mass matrix for a <u,v> = <f,v> bilinear form using various element types
   * (eg. P1, P2, Nedelec, and Raviart-Thomas).
   * 
   * For this problem we compute:
   *
   * coef * u = f  ---->   coef*< u, v> = <f,v>
   *
   * which gives Mx = b,
   *
   * M_ij = coef *<phi_j,phi_i>
   * b_i  = <f,phi_i>
   * 
   *    INPUT:
   *            fespace		    Finite-Element Space Struct
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights
   *            dof_on_elm          Specific DOF on element
   *            elm                 Current element
   *            coeff               Function that gives coefficient (for now assume constant)
   *            time                Time if time dependent
   *
   *    OUTPUT:
   *	       	MLoc                Local Mass Matrix (Full Matrix)
   */
	
  // Mesh and FE data
  INT dof_per_elm = FE->dof_per_elm;
  INT dim = mesh->dim;
  
  // Loop Indices
  INT quad,test,trial;
  // Quadrature Weights and Nodes
  REAL w; 
  REAL* qx = (REAL *) calloc(3,sizeof(REAL));
  // Stiffness Matrix Entry
  REAL kij;
	
  // Basis Functions and its derivatives if necessary
  REAL* phi=NULL;
  REAL* dphix=NULL;
  REAL* dphiy=NULL;
  REAL* dphiz=NULL;
  
  // Coefficient Value at Quadrature Nodes
  REAL coeff_val=0.0;


  if(FE->FEtype>0) { // PX elements
    phi = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    dphix = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    dphiy = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    if(dim==3) dphiz = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    
    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {        
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      (*coeff)(&coeff_val,qx,time);
      
      //  Get the Basis Functions at each quadrature node
      PX_H1_basis(phi,dphix,dphiy,dphiz,qx[0],qx[1],qx[2],dof_on_elm,FE->FEtype,mesh);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->dof_per_elm; trial++) {
	  kij = coeff_val*(phi[test]*phi[trial]);
	  MLoc[test*FE->dof_per_elm+trial] = MLoc[test*FE->dof_per_elm+trial] + w*kij;
	}
      }
    }
  } else if(FE->FEtype==-1) { // Nedelec elements
    phi = (REAL *) calloc(dof_per_elm*dim,sizeof(REAL));
    if(dim==2) {
      dphix = (REAL *) calloc(dof_per_elm,sizeof(REAL)); // Curl of basis function
    } else if (dim==3) {
      dphix = (REAL *) calloc(dof_per_elm*dim,sizeof(REAL)); // Curl of basis function
    } else {
      baddimension();
    }

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {        
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      (*coeff)(&coeff_val,qx,time);

      //  Get the Basis Functions at each quadrature node
      ned_basis(phi,dphix,qx[0],qx[1],qx[2],v_on_elm,dof_on_elm,mesh);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->dof_per_elm; trial++) {
	  kij = coeff_val*(phi[test*dim]*phi[trial*dim] + phi[test*dim+1]*phi[trial*dim+1]);
	  if(dim==3) kij = kij + coeff_val*phi[test*dim+2]*phi[trial*dim+2];
	  MLoc[test*FE->dof_per_elm+trial] = MLoc[test*FE->dof_per_elm+trial] + w*kij;
	}
      }
    }
  } else if(FE->FEtype==-2) { // Raviart-Thomas elements
    phi = (REAL *) calloc(dof_per_elm*dim,sizeof(REAL));
    dphix = (REAL *) calloc(dof_per_elm,sizeof(REAL)); // Divergence of element 

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {        
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      (*coeff)(&coeff_val,qx,time);

      //  Get the Basis Functions at each quadrature node
      rt_basis(phi,dphix,qx[0],qx[1],qx[2],v_on_elm,dof_on_elm,mesh);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->dof_per_elm; trial++) {
	  kij = coeff_val*(phi[test*dim]*phi[trial*dim] + phi[test*dim+1]*phi[trial*dim+1]);
	  if(dim==3) kij = kij + coeff_val*phi[test*dim+2]*phi[trial*dim+2];
	  MLoc[test*FE->dof_per_elm+trial] = MLoc[test*FE->dof_per_elm+trial] + w*kij;
	}
      }
    }
  } else {
    printf("Trying to implement non-existent elements!\n\n");
    exit(0);
  }
  
  if (phi) free(phi);
  if(dphix) free(dphix);
  if(dphiy) free(dphiy);
  if(dphiz) free(dphiz);
  if(qx) free(qx);
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void FEM_RHS_Local(REAL* bLoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT elm,void (*rhs)(REAL *,REAL *,REAL),REAL time) 
{
	
  /* Computes the local Right hand side vector for Galerkin Finite Elements
   * b_i  = <f,phi_i>
   * 
   *
   *    INPUT:
   *            fespace		    Finite-Element Space Struct
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights
   *            dof_on_elm          Specific DOF on element
   *            elm                 Current element
   *            rhs                 Function that gives coefficient (for now assume constant)
   *            time                Time if time dependent
   *
   *    OUTPUT:
   *	       	ALoc                Local Stiffness Matrix (Full Matrix)
   */

  // Loop Indices
  INT quad,test;
  // Quadrature Weights and Nodes
  REAL w; 
  REAL* qx = (REAL *) calloc(3,sizeof(REAL));
	
  // Basis Functions
  REAL* p = (REAL *) calloc(FE->dof_per_elm,sizeof(REAL));
  REAL* dpx = (REAL *) calloc(FE->dof_per_elm,sizeof(REAL));
  REAL* dpy = (REAL *) calloc(FE->dof_per_elm,sizeof(REAL));
  REAL* dpz = (REAL *) calloc(FE->dof_per_elm,sizeof(REAL));
  
  // Right-hand side function at Quadrature Nodes
  REAL rhs_val=0.0;
  
  //  Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {        
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];
    (*rhs)(&rhs_val,qx,time);
    
    //  Get the Basis Functions at each quadrature node
    if(FE->FEtype>0) { // PX elements
      PX_H1_basis(p,dpx,dpy,dpz,qx[0],qx[1],qx[2],dof_on_elm,FE->FEtype,mesh);
    } else if(FE->FEtype==-1) { // Nedelec elements
      printf("Nedelec Not Implemented Yet\n!");
    } else if(FE->FEtype==-2) { // Raviart-Thomas elements
      printf("RT Not Implemented Yet!\n");
    } else {
      printf("Trying to implement non-existent elements!\n\n");
      exit(0);
    }	

    // Loop over test functions and integrate rhs
    for (test=0; test<FE->dof_per_elm;test++) {
      bLoc[test] = bLoc[test] + w*rhs_val*p[test];
    }
  }
		
  if(p) free(p);
  if(dpx) free(dpx);
  if(dpy) free(dpy);
  if(dpz) free(dpz);
  if(qx) free(qx);
	
  return;
}
/******************************************************************************************************/

