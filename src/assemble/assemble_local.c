/*
 *  assemble_local.c   
 *  
 *  Created by James Adler and Xiaozhe Hu on 4/22/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

/* This code will build local stiffness matrices for various PDE systems */

#include "hazmat.h"

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
      if(coeff!=NULL) {
	(*coeff)(&coeff_val,qx,time);
      } else {
	coeff_val = 1.0;
      }
      
      //  Get the Basis Functions at each quadrature node
      PX_H1_basis(phi,dphix,dphiy,dphiz,qx,dof_on_elm,FE->FEtype,mesh);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->dof_per_elm; trial++) {
	  kij = coeff_val*(dphix[test]*dphix[trial] + dphiy[test]*dphiy[trial]);
	  if(dim==3) kij+=coeff_val*dphiz[test]*dphiz[trial];
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
      if(coeff!=NULL) {
	(*coeff)(&coeff_val,qx,time);
      } else {
	coeff_val = 1.0;
      }

      //  Get the Basis Functions at each quadrature node
      ned_basis(phi,dphix,qx,v_on_elm,dof_on_elm,mesh);

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
      if(coeff!=NULL) {
	(*coeff)(&coeff_val,qx,time);
      } else {
	coeff_val = 1.0;
      }

      //  Get the Basis Functions at each quadrature node
      rt_basis(phi,dphix,qx,v_on_elm,dof_on_elm,mesh);

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
	
  /* Computes the local mass matrix for a <u,v> = <f,v> bilinear form using various element types
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
      if(coeff!=NULL) {
	(*coeff)(&coeff_val,qx,time);
      } else {
	coeff_val = 1.0;
      }
      
      //  Get the Basis Functions at each quadrature node
      PX_H1_basis(phi,dphix,dphiy,dphiz,qx,dof_on_elm,FE->FEtype,mesh);

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
      if(coeff!=NULL) {
	(*coeff)(&coeff_val,qx,time);
      } else {
	coeff_val = 1.0;
      }

      //  Get the Basis Functions at each quadrature node
      ned_basis(phi,dphix,qx,v_on_elm,dof_on_elm,mesh);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->dof_per_elm; trial++) {
	  kij = coeff_val*(phi[test*dim]*phi[trial*dim] + phi[test*dim+1]*phi[trial*dim+1]);
	  if(dim==3) kij += coeff_val*phi[test*dim+2]*phi[trial*dim+2];
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
      if(coeff!=NULL) {
	(*coeff)(&coeff_val,qx,time);
      } else {
	coeff_val = 1.0;
      }

      //  Get the Basis Functions at each quadrature node
      rt_basis(phi,dphix,qx,v_on_elm,dof_on_elm,mesh);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->dof_per_elm; trial++) {
	  kij = coeff_val*(phi[test*dim]*phi[trial*dim] + phi[test*dim+1]*phi[trial*dim+1]);
	  if(dim==3) kij += coeff_val*phi[test*dim+2]*phi[trial*dim+2];
	  MLoc[test*FE->dof_per_elm+trial] = MLoc[test*FE->dof_per_elm+trial] + w*kij;
	}
      }
    }
  } else {
    printf("Trying to implement non-existent elements!\n\n");
    exit(0);
  }

  /* for (test=0; test<FE->dof_per_elm;test++) { */
  /*   // Loop over Trial Functions (Columns) */
  /*   for (trial=0; trial<FE->dof_per_elm; trial++) { */
  /*     fprintf(stdout,"%23.16g ", MLoc[test*FE->dof_per_elm+trial]); */
  /*   } */
  /*   //fprintf(stdout,"\n"); */
  /* } */
  /* fprintf(stdout,"xxxxx\n"); */
  
  
  if (phi) free(phi);
  if(dphix) free(dphix);
  if(dphiy) free(dphiy);
  if(dphiz) free(dphiz);
  if(qx) free(qx);
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void FEM_RHS_Local(REAL* bLoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*rhs)(REAL *,REAL *,REAL),REAL time) 
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

  // Mesh and FE data
  INT dof_per_elm = FE->dof_per_elm;
  INT dim = mesh->dim;
  
  // Loop Indices
  INT quad,test;
  
  // Quadrature Weights and Nodes
  REAL w; 
  REAL* qx = (REAL *) calloc(3,sizeof(REAL));
	
  // Basis Functions and its derivatives if necessary
  REAL* phi=NULL;
  REAL* dphix=NULL;
  REAL* dphiy=NULL;
  REAL* dphiz=NULL;
  
  // Right-hand side function at Quadrature Nodes
  REAL* rhs_val=NULL;
  
  if(FE->FEtype>0) { // PX elements
    phi = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    dphix = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    dphiy = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    if(dim==3) dphiz = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    rhs_val = (REAL *) calloc(1,sizeof(REAL));
    
    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {        
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      (*rhs)(rhs_val,qx,time);
    
      //  Get the Basis Functions at each quadrature node
      PX_H1_basis(phi,dphix,dphiy,dphiz,qx,dof_on_elm,FE->FEtype,mesh);

      // Loop over test functions and integrate rhs
      for (test=0; test<FE->dof_per_elm;test++) {
	bLoc[test] = bLoc[test] + w*rhs_val[0]*phi[test];
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
    rhs_val = (REAL *) calloc(dim,sizeof(REAL));

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {        
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      (*rhs)(rhs_val,qx,time);

      //  Get the Basis Functions at each quadrature node
      ned_basis(phi,dphix,qx,v_on_elm,dof_on_elm,mesh);

      // Loop over test functions and integrate rhs
      for (test=0; test<FE->dof_per_elm;test++) {
	bLoc[test] = bLoc[test] + w*(rhs_val[0]*phi[test*dim] + rhs_val[1]*phi[test*dim+1]);
	if(dim==3) bLoc[test] = bLoc[test] + w*rhs_val[2]*phi[test*dim+2];
      }
    }
  } else if(FE->FEtype==-2) { // Raviart-Thomas elements
    phi = (REAL *) calloc(dof_per_elm*dim,sizeof(REAL));
    dphix = (REAL *) calloc(dof_per_elm,sizeof(REAL)); // Divergence of element
    rhs_val = (REAL *) calloc(dim,sizeof(REAL));

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {        
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      (*rhs)(rhs_val,qx,time);

      //  Get the Basis Functions at each quadrature node
      rt_basis(phi,dphix,qx,v_on_elm,dof_on_elm,mesh);

      // Loop over test functions and integrate rhs
      for (test=0; test<FE->dof_per_elm;test++) {
	bLoc[test] = bLoc[test] + w*(rhs_val[0]*phi[test*dim] + rhs_val[1]*phi[test*dim+1]);
	if(dim==3) bLoc[test] = bLoc[test] + w*rhs_val[2]*phi[test*dim+2];
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
  if(rhs_val) free(rhs_val);
	
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void assemble_DuDvplusmass_local(REAL* ALoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time) 
{
	
  /* Computes the local stiffness matrix for coeff1*<Du,Dv> + coeff2*<u,v> = <f,v> bilinear form using various element types
   * (eg. P1, P2 -> (grad u, grad v) + (u,v), Nedelec <curl u, curl v> + <u,v>, and Raviart-Thomas <div u, div v> + <u,v>).
   * 
   * For this problem we compute:
   *
   * coeff1*D*D u + coeff2*u = f  ---->   coeff1*<D u, D v> + coeff2*<u,v> = <f,v>
   *
   * which gives Ax = b,
   *
   * A_ij = coeff[0]*<D phi_j,D phi_i> + coeff[1]*<phi_j,phi_i>
   * b_i  = <f,phi_i>
   * 
   *    INPUT:
   *            fespace		    Finite-Element Space Struct
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights
   *            dof_on_elm          Specific DOF on element
   *            elm                 Current element
   *            coeff               Function that gives coefficient (should be vector of 2: coeff[0] for <Du,Dv> term and coeff[1] for <u,v> term)
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
  REAL* coeff_val = (REAL *) calloc(2,sizeof(REAL));
	
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
      if(coeff!=NULL) {
	(*coeff)(coeff_val,qx,time);
      } else {
	coeff_val[0] = 1.0;
	coeff_val[1] = 1.0;
      }
      
      //  Get the Basis Functions at each quadrature node
      PX_H1_basis(phi,dphix,dphiy,dphiz,qx,dof_on_elm,FE->FEtype,mesh);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->dof_per_elm; trial++) {
	  kij = coeff_val[0]*(dphix[test]*dphix[trial] + dphiy[test]*dphiy[trial]) + coeff_val[1]*(phi[test]*phi[trial]);
	  if(dim==3) kij+=coeff_val[0]*dphiz[test]*dphiz[trial];
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
      if(coeff!=NULL) {
	(*coeff)(coeff_val,qx,time);
      } else {
	coeff_val[0] = 1.0;
	coeff_val[1] = 1.0;
      }

      //  Get the Basis Functions at each quadrature node
      ned_basis(phi,dphix,qx,v_on_elm,dof_on_elm,mesh);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->dof_per_elm; trial++) {
	  if(dim==2) {
	    kij = coeff_val[0]*(dphix[test]*dphix[trial]) + coeff_val[1]*(phi[test*dim]*phi[trial*dim] + phi[test*dim+1]*phi[trial*dim+1]);
	  } else if (dim==3) {
	    kij = coeff_val[0]*(dphix[test*dim]*dphix[trial*dim] + dphix[test*dim+1]*dphix[trial*dim+1] + dphix[test*dim+2]*dphix[trial*dim+2]) + 
	      coeff_val[1]*(phi[test*dim]*phi[trial*dim] + phi[test*dim+1]*phi[trial*dim+1] + phi[test*dim+2]*phi[trial*dim+2]);
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
      if(coeff!=NULL) {
	(*coeff)(coeff_val,qx,time);
      } else {
	coeff_val[0] = 1.0;
	coeff_val[1] = 1.0;
      }

      //  Get the Basis Functions at each quadrature node
      rt_basis(phi,dphix,qx,v_on_elm,dof_on_elm,mesh);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->dof_per_elm; trial++) {
	  kij = coeff_val[0]*(dphix[test]*dphix[trial]) + 
	    coeff_val[1]*(phi[test*dim]*phi[trial*dim] + phi[test*dim+1]*phi[trial*dim+1]);
	  if(dim==3) kij += coeff_val[1]*phi[test*dim+2]*phi[trial*dim+2];
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
  if(coeff_val) free(coeff_val);
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void impedancebdry_local(REAL* ZLoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *ed_on_f, \
			 INT *ed_on_elm,INT *v_on_elm,INT face,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time) 
{

  /* Computes the local weak formulation of the Impedance boundary condition
   *
   * Uses midpoint rule to integrate on edges of boundary face
   * ASSUMING 3D ONLY
   * For this problem we compute the left-hand side of:
   *
   *    <n x E,n x F>_bdryobstacle    for all F in H_imp(curl) (Nedelec)
   * 
   *
   *    INPUT:
   *            FE		    Finite-Element Space Struct for E
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights
   *            ed_on_f             Specific Edge on boundary face
   *            v_on_elm            Vertices on current element associated with current face
   *            elm                 Current element associated with current face
   *            face                Current face on boundary
   *            coeff               Function that gives coefficient (for now assume constant)
   *            time                Time if time dependent
   *
   *    OUTPUT:
   *	       	ZLoc                Local Boundary Matrix (Full Matrix)
   */
	
  // Mesh and FE data
  INT ed_per_elm = FE->dof_per_elm;
  INT dim = mesh->dim;
  
  // Loop Indices
  INT j,quad,test,trial,ed,edt,edb;

  // Quadrature Weights and Nodes
  INT nq = 2*dim-3; // = ed_per_face
  REAL* qx = (REAL *) calloc(nq,sizeof(REAL));
  // 3D: Using triangle midpoint rule, so qx is midpoint of edges and w is |F|/3
  REAL w = mesh->f_area[face]/3.0; 

  // Get normal vector components on face
  REAL nx = mesh->f_norm[face*dim];
  REAL ny = mesh->f_norm[face*dim+1];
  REAL nz = mesh->f_norm[face*dim+2];

  // Stiffness Matrix Entry
  REAL kij,kij1,kij2,kij3,kij4,kij5,kij6;
	
  // Basis Functions and its curl
  REAL* phi= (REAL *) calloc(ed_per_elm*dim,sizeof(REAL));
  REAL* cphi = (REAL *) calloc(ed_per_elm*dim,sizeof(REAL));
  
  // Coefficient Value at Quadrature Nodes
  REAL coeff_val=0.0;

  //  Sum over midpoints of edges
  for (quad=0;quad<nq;quad++) { 
    ed = ed_on_f[quad]-1;
    qx[0] = mesh->ed_mid[ed*dim];
    qx[1] = mesh->ed_mid[ed*dim+1];
    qx[2] = mesh->ed_mid[ed*dim+2];

    if(coeff!=NULL) {
	(*coeff)(&coeff_val,qx,time);
      } else {
	coeff_val = 1.0;
      }

    //  Get the Basis Functions at each quadrature node
    ned_basis(phi,cphi,qx,v_on_elm,ed_on_elm,mesh);

    // Loop over Test Functions (Rows - edges)
    for (test=0; test<nq;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<nq; trial++) {
	// Make sure ordering for global matrix is right
	for(j=0;j<mesh->ed_per_elm;j++) {
	  if(ed_on_f[test]==ed_on_elm[j]) {
	    edt = j;
	  }
	  if(ed_on_f[trial]==ed_on_elm[j]) {
	    edb = j;
	  }
	}
	kij1 = phi[edb*dim+1]*nz - phi[edb*dim+2]*ny;
	kij2 = phi[edb*dim+2]*nx - phi[edb*dim]*nz;
	kij3 = phi[edb*dim]*ny - phi[edb*dim+1]*nx;
	kij4 = phi[edt*dim+1]*nz - phi[edt*dim+2]*ny;
	kij5 = phi[edt*dim+2]*nx - phi[edt*dim]*nz;
	kij6 = phi[edt*dim]*ny - phi[edt*dim+1]*nx;
	kij = coeff_val*(kij1*kij4+kij2*kij5+kij3*kij6);
	ZLoc[test*nq+trial]+=w*kij;
      }
    }
  }
  
  if (phi) free(phi);
  if(cphi) free(cphi);
  if(qx) free(qx);
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void boundary_mass_local(REAL* ALoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_f, \
			 INT *dof_on_elm,INT *v_on_elm,INT face,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time)
{
  /* Computes the local weak formulation of the mass matrix on a boundary face (3D -> 2D surface; 2D -> 1D curve)
   *
   * For this problem we compute the left-hand side of:
   *
   *    <u,v>_bdryobstacle    for all v
   *
   *
   *    INPUT:
   *            FE		    Finite-Element Space Struct for E
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights on surface/curve
   *            dof_on_f            Specific DOF on boundary face
   *            v_on_elm            Vertices on current element associated with current face
   *            elm                 Current element associated with current face
   *            face                Current face on boundary
   *            coeff               Function that gives coefficient (for now assume constant)
   *            time                Time if time dependent
   *
   *    OUTPUT:
   *	       	ALoc                Local Boundary Matrix (Full Matrix)
   */
	
  // Mesh and FE data
  INT dof_per_elm = FE->dof_per_elm;
  INT dim = mesh->dim;
  INT dof_per_f = 0;
  
  // Loop Indices
  INT j,quad,test,trial,doft,dofb;

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

    // Get DOF Per Face
    if(dim==2) {
      dof_per_f = FE->FEtype+1;
    } else if(dim==3) {
      dof_per_f = 3*FE->FEtype;
    } else {
      baddimension();
    }

    //  Sum over midpoints of edges
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[face*cq->nq_per_elm+quad];
      qx[1] = cq->y[face*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[face*cq->nq_per_elm+quad];
      w = cq->w[face*cq->nq_per_elm+quad];

      if(coeff!=NULL) {
	(*coeff)(&coeff_val,qx,time);
      } else {
	coeff_val = 1.0;
      }
    
      //  Get the Basis Functions at each quadrature node
      PX_H1_basis(phi,dphix,dphiy,dphiz,qx,dof_on_elm,FE->FEtype,mesh);

      // Loop over Test Functions (Rows - vertices)
      for (test=0; test<dof_per_f;test++) {
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<dof_per_f; trial++) {
	  // Make sure ordering for global matrix is right
	  for(j=0;j<FE->dof_per_elm;j++) {
	    if(dof_on_f[test]==dof_on_elm[j]) {
	      doft = j;
	    }
	    if(dof_on_f[trial]==dof_on_elm[j]) {
	      dofb = j;
	    }
	  }
	  kij = coeff_val*(phi[dofb]*phi[doft]);
	  ALoc[test*dof_per_f+trial]+=w*kij;
	}
      }
    }
  } else if(FE->FEtype==-1) {

    // Basis Functions and its curl
    phi = (REAL *) calloc(dof_per_elm*dim,sizeof(REAL));
    if(dim==2) {
      dphix = (REAL *) calloc(dof_per_elm,sizeof(REAL)); // Curl of basis function
    } else if (dim==3) {
      dphix = (REAL *) calloc(dof_per_elm*dim,sizeof(REAL)); // Curl of basis function
    } else {
      baddimension();
    }

    // Get DOF Per Face
    if(dim==2) {
      dof_per_f = 1;
    } else if(dim==3) {
      dof_per_f = 3;
    } else {
      baddimension();
    }

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {        
      qx[0] = cq->x[face*cq->nq_per_elm+quad];
      qx[1] = cq->y[face*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[face*cq->nq_per_elm+quad];
      w = cq->w[face*cq->nq_per_elm+quad];

      if(coeff!=NULL) {
	(*coeff)(&coeff_val,qx,time);
      } else {
	coeff_val = 1.0;
      }

      //  Get the Basis Functions at each quadrature node
      ned_basis(phi,dphix,qx,v_on_elm,dof_on_elm,mesh);

      // Loop over Test Functions (Rows - edges)
      for (test=0; test<dof_per_f;test++) {
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<dof_per_f; trial++) {
	  // Make sure ordering for global matrix is right
	  for(j=0;j<FE->dof_per_elm;j++) {
	    if(dof_on_f[test]==dof_on_elm[j]) {
	      doft = j;
	    }
	    if(dof_on_f[trial]==dof_on_elm[j]) {
	      dofb = j;
	    }
	  }
	  kij = coeff_val*(phi[dofb*dim]*phi[doft*dim] + phi[dofb*dim+1]*phi[doft*dim+1]);
	  if(dim==3) kij += coeff_val*(phi[dofb*dim+2]*phi[doft*dim+2]);
	  ALoc[test*dof_per_f+trial]+=w*kij;
	}
      }
    }
  } else if(FE->FEtype==-2) { // Raviart-Thomas elements
    phi = (REAL *) calloc(FE->dof_per_elm*dim,sizeof(REAL));
    dphix = (REAL *) calloc(FE->dof_per_elm,sizeof(REAL)); // Divergence of element

    // Get DOF Per Face
    if(dim==2) {
      dof_per_f = 1;
    } else if(dim==3) {
      dof_per_f = 1;
    } else {
      baddimension();
    }

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {        
      qx[0] = cq->x[face*cq->nq_per_elm+quad];
      qx[1] = cq->y[face*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[face*cq->nq_per_elm+quad];
      w = cq->w[face*cq->nq_per_elm+quad];

      if(coeff!=NULL) {
	(*coeff)(&coeff_val,qx,time);
      } else {
	coeff_val = 1.0;
      }

      //  Get the Basis Functions at each quadrature node
      rt_basis(phi,dphix,qx,v_on_elm,dof_on_elm,mesh);

      /// Loop over Test Functions (Rows - edges)
      for (test=0; test<dof_per_f;test++) {
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<dof_per_f; trial++) {
	  // Make sure ordering for global matrix is right
	  for(j=0;j<FE->dof_per_elm;j++) {
	    if(dof_on_f[test]==dof_on_elm[j]) {
	      doft = j;
	    }
	    if(dof_on_f[trial]==dof_on_elm[j]) {
	      dofb = j;
	    }
	  }
	  kij = coeff_val*(phi[dofb*dim]*phi[doft*dim] + phi[dofb*dim+1]*phi[doft*dim+1]);
	  kij += coeff_val*(phi[dofb*dim+2]*phi[doft*dim+2]);
	  ALoc[test*dof_per_f+trial]+=w*kij;
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
void Ned_GradH1_RHS_local(REAL* bLoc,fespace *FE_H1,fespace *FE_Ned,trimesh *mesh,qcoordinates *cq,INT *ed_on_elm,INT *v_on_elm,INT elm,dvector* u)  
{

  /* Computes the local weak formulation of <E,grad(q)> where E is a given
   * Nedelec approximation and q in H_0^1 (linears)
   *
   *    INPUT:
   *            FE_H1		    Finite-Element Space Struct for H1 elements
   *            FE_Ned              Finite-Element Space Struct for Nedelec elements
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights
   *            ed_on_elm           Specific Edges on element
   *            v_on_elm            Vertices on current element 
   *            elm                 Current element 
   *            u                   FEM Function that gives coefficient
   *
   *    OUTPUT:
   *	       	bLoc                Local rhs
   */
	
   // Mesh and FE data
  INT v_per_elm = FE_H1->dof_per_elm;
  INT dim = mesh->dim;
  
  // Loop Indices
  INT quad,test;
  
  // Quadrature Weights and Nodes
  REAL w; 
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));
	
  // Basis Functions and its derivatives if necessary
  REAL* phi=NULL;
  REAL* dphix=NULL;
  REAL* dphiy=NULL;
  REAL* dphiz=NULL;
  phi = (REAL *) calloc(v_per_elm,sizeof(REAL));
  dphix = (REAL *) calloc(v_per_elm,sizeof(REAL));
  dphiy = (REAL *) calloc(v_per_elm,sizeof(REAL));
  if(dim==3) dphiz = (REAL *) calloc(v_per_elm,sizeof(REAL));
  
  // Right-hand side function at Quadrature Nodes
  REAL* ucoeff = (REAL *) calloc(dim,sizeof(REAL));
  
  //  Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {        
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];
    
    // Get FEM function at quadrature nodes
    FE_Interpolation(ucoeff,u->val,qx,ed_on_elm,v_on_elm,FE_Ned,mesh,mesh->nedge,1);
    
    //  Get the Basis Functions at each quadrature node
    PX_H1_basis(phi,dphix,dphiy,dphiz,qx,v_on_elm,FE_H1->FEtype,mesh);

    // Loop over test functions and integrate rhs
    for (test=0; test<FE_H1->dof_per_elm;test++) {
      bLoc[test] += w*(ucoeff[0]*dphix[test]+ucoeff[1]*dphiy[test]);
      if(dim==3) bLoc[test] += w*ucoeff[2]*dphiz[test];
    }
  }
  
  if(phi) free(phi);
  if(dphix) free(dphix);
  if(dphiy) free(dphiy);
  if(dphiz) free(dphiz);
  if(qx) free(qx);
  if(ucoeff) free(ucoeff);
	
  return;
}
/******************************************************************************************************/


/////////////////////// TEMPORARY ////////////////////////////////////////////////////////////////////

/********************** Local Assembly **************************************************/
void myJacobian_MHD(REAL* ALoc,dvector* old_sol,block_fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,REAL time) 
{
	
  /* Computes the local stiffness matrix for the MHD system 
   * 
   * For this problem we compute first component of:
   *
   * <2/Rs ep(u),ep(v)> + <un*grad u + u*grad un, v> = <f,v>
   *
   * 2/Rs(<u1x, v1x> + 0.5<u2x+u1y,v2x+v1y> + <u2y,v2y>) + <un*grad u1 + u1*un1x + u2*un1y,v1> + <un*grad u2 + u1*un2x + u2*un2y,v2> = <f1,v1> + <f2,v2>
   *
   * u1-u1 block:
   *
   * 2/Rs(<u1x, v1x> + 0.5<u1y,v1y>) + <un*grad u1 + u1*un1x,v1> = <f1,v1>
   *
   * 
   *    INPUT:
   *            uprev               Solution at previous Newton step
   *            FE		    Finite-Element Space Struct in block form
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights
   *            dof_on_elm          Specific DOF on element
   *            v_on_elm            Vertices on element
   *            elm                 Current element
   *            time                Time if time dependent
   *
   *    OUTPUT:
   *	       	ALoc                Local Stiffness Matrix (Full Matrix)
   */
	
  INT i,j;

  // Mesh and FE data
  INT dof_per_elm = 0;
  for(i=0;i<FE->nspaces;i++) 
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  INT* local_dof_on_elm = NULL;
  REAL* local_uprev = NULL;

  INT dim = mesh->dim;
  
  // Loop Indices
  INT quad,test,trial;
  // Quadrature Weights and Nodes
  REAL w; 
  REAL* qx = (REAL *) calloc(3,sizeof(REAL));
  // Stiffness Matrix Entry
  REAL kij = 0.0;

  // Data
  REAL Rs=0.0;
  REAL Rm=0.0;
  
	
  // Basis Functions and its derivatives if necessary
  // u
  REAL* u1_phi = (REAL *) calloc(FE->var_spaces[0]->dof_per_elm,sizeof(REAL));
  REAL* u1x_phi = (REAL *) calloc(FE->var_spaces[0]->dof_per_elm,sizeof(REAL));
  REAL* u1y_phi = (REAL *) calloc(FE->var_spaces[0]->dof_per_elm,sizeof(REAL));
  REAL* u2_phi = (REAL *) calloc(FE->var_spaces[1]->dof_per_elm,sizeof(REAL));
  REAL* u2x_phi = (REAL *) calloc(FE->var_spaces[1]->dof_per_elm,sizeof(REAL));
  REAL* u2y_phi = (REAL *) calloc(FE->var_spaces[1]->dof_per_elm,sizeof(REAL));
  REAL u1_prev;
  REAL u2_prev;
  REAL du1_prev[dim];
  REAL du2_prev[dim];
  // p
  REAL* p_phi = (REAL *) calloc(FE->var_spaces[2]->dof_per_elm,sizeof(REAL));
  REAL* px_phi = (REAL *) calloc(FE->var_spaces[2]->dof_per_elm,sizeof(REAL));
  REAL* py_phi = (REAL *) calloc(FE->var_spaces[2]->dof_per_elm,sizeof(REAL));
  REAL p_prev;
  // B
  REAL* B_phi = (REAL *) calloc(FE->var_spaces[3]->dof_per_elm*dim,sizeof(REAL));
  REAL* curlB_phi = (REAL *) calloc(FE->var_spaces[3]->dof_per_elm,sizeof(REAL));
  REAL B_prev[dim];
  REAL curlB_prev;
  // r
  REAL* r_phi = (REAL *) calloc(FE->var_spaces[4]->dof_per_elm,sizeof(REAL));
  REAL* rx_phi = (REAL *) calloc(FE->var_spaces[4]->dof_per_elm,sizeof(REAL));
  REAL* ry_phi = (REAL *) calloc(FE->var_spaces[4]->dof_per_elm,sizeof(REAL));
  REAL dr_prev[dim];

  INT local_row_index, local_col_index;
      
  //  Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {        
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];
    
    //  Get the Basis Functions and old solutions at each quadrature node
    // u1
    PX_H1_basis(u1_phi,u1x_phi,u1y_phi,NULL,qx,dof_on_elm,FE->var_spaces[0]->FEtype,mesh);
    FE_Interpolation(&u1_prev,old_sol->val,qx,dof_on_elm,v_on_elm,FE->var_spaces[0],mesh,FE->var_spaces[0]->ndof,1);
    FE_DerivativeInterpolation(du1_prev,old_sol->val,qx,dof_on_elm,v_on_elm,FE->var_spaces[0],mesh,FE->var_spaces[0]->ndof,1);
    local_uprev = old_sol->val + FE->var_spaces[0]->ndof;
    local_dof_on_elm = dof_on_elm + FE->var_spaces[0]->dof_per_elm;
    // u2
    PX_H1_basis(u2_phi,u2x_phi,u2y_phi,NULL,qx,local_dof_on_elm,FE->var_spaces[1]->FEtype,mesh);
    FE_Interpolation(&u2_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh,FE->var_spaces[1]->ndof,1);
    FE_DerivativeInterpolation(du1_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh,FE->var_spaces[1]->ndof,1);
    local_uprev = local_uprev + FE->var_spaces[1]->ndof;
    local_dof_on_elm = local_dof_on_elm + FE->var_spaces[1]->dof_per_elm;
    // p
    PX_H1_basis(p_phi,px_phi,py_phi,NULL,qx,local_dof_on_elm,FE->var_spaces[2]->FEtype,mesh);
    FE_Interpolation(&p_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[2],mesh,FE->var_spaces[2]->ndof,1);
    local_uprev = local_uprev + FE->var_spaces[2]->ndof;
    local_dof_on_elm = local_dof_on_elm + FE->var_spaces[2]->dof_per_elm;
    // B
    ned_basis(B_phi,curlB_phi,qx,v_on_elm,local_dof_on_elm,mesh);
    FE_Interpolation(B_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[3],mesh,FE->var_spaces[3]->ndof,1);
    FE_DerivativeInterpolation(&curlB_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[3],mesh,FE->var_spaces[3]->ndof,1);
    local_uprev = local_uprev + FE->var_spaces[3]->ndof;
    local_dof_on_elm = local_dof_on_elm + FE->var_spaces[3]->dof_per_elm;
    // r
    PX_H1_basis(r_phi,rx_phi,ry_phi,NULL,qx,local_dof_on_elm,FE->var_spaces[4]->FEtype,mesh);
    FE_DerivativeInterpolation(dr_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[4],mesh,FE->var_spaces[4]->ndof,1);
    // Data
    Rs = 1;
    Rm = 1;

    // u1-u1 block
    local_row_index = 0;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm; trial++) {
	kij = (2.0/Rs)*(u1x_phi[trial]*u1x_phi[test] + 0.5*u1y_phi[trial]*u1y_phi[test]) 
	  + (u1_prev*u1x_phi[trial] + u2_prev*u1y_phi[trial] + u1_phi[trial]*du1_prev[0])*u1_phi[test];
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }
    
    // u1-u2 block
    local_col_index += FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm; trial++) {
	kij = (2.0/Rs)*(0.5*u2x_phi[trial]*u1y_phi[test]) 
	  + (u2_phi[trial]*du1_prev[1])*u1_phi[test];
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }

    // u1-p block
    local_col_index += FE->var_spaces[1]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[2]->dof_per_elm; trial++) {
	kij = -p_phi[trial]*(u1x_phi[test]);
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }

    // u1-B block
    local_col_index += FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[3]->dof_per_elm; trial++) {
	kij = -(-curlB_prev*B_phi[trial*dim+1] - curlB_phi[trial]*B_prev[1])*(u1_phi[test]);
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }

    // u2-u1 block
    local_row_index += FE->var_spaces[0]->dof_per_elm;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm; trial++) {
	kij = (2.0/Rs)*(0.5*u1y_phi[trial]*u2x_phi[test]) 
	  + (u1_phi[trial]*du2_prev[0])*u1_phi[test];
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }
    
    // u2-u2 block
    local_col_index += FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm; trial++) {
	kij = (2.0/Rs)*(0.5*u2x_phi[trial]*u2x_phi[test] + u2y_phi[trial]*u2y_phi[test]) 
	  + (u1_prev*u2x_phi[trial]+u2_prev*u2y_phi[trial] + u2_phi[trial]*du2_prev[1])*u2_phi[test];
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }

    // u2-p block
    local_col_index += FE->var_spaces[1]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[2]->dof_per_elm; trial++) {
	kij = -p_phi[trial]*(u2y_phi[test]);
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }

    // u2-B block
    local_col_index += FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[3]->dof_per_elm; trial++) {
	kij = -(curlB_prev*B_phi[trial*dim] + curlB_phi[trial]*B_prev[0])*(u2_phi[test]);
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }

    // p-u1 block
    local_row_index += FE->var_spaces[1]->dof_per_elm;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm; trial++) {
	kij = p_phi[test]*u1x_phi[trial];
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }
    
    // p-u2 block
    local_col_index += FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm; trial++) {
	kij = p_phi[test]*u2y_phi[trial];
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }
    
    // B-u1 block
    local_row_index += FE->var_spaces[2]->dof_per_elm;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm; trial++) {
	kij = -u1_phi[trial]*B_prev[1]*curlB_phi[test];
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }
    
    // B-u2 block
    local_col_index += FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm; trial++) {
	kij = u2_phi[trial]*B_prev[0]*curlB_phi[test];
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }

    // B-p block is empty
    local_col_index += FE->var_spaces[1]->dof_per_elm;

    // B-B
    local_col_index += FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[3]->dof_per_elm; trial++) {
	kij = (1.0/Rm)*curlB_phi[trial]*curlB_phi[test] 
	  - (u1_prev*B_phi[trial*dim+1] - u2_prev*B_phi[trial*dim])*curlB_phi[test];
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }

    // B-r block
    local_col_index += FE->var_spaces[3]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[4]->dof_per_elm; trial++) {
	kij = -rx_phi[trial]*B_phi[test*dim] - ry_phi[trial]*B_phi[test*dim+1];
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }
    
    // r-u1 block
    local_row_index += FE->var_spaces[3]->dof_per_elm;
    local_col_index = 0;
  
    // r-u2 block
    local_col_index += FE->var_spaces[0]->dof_per_elm;
   
    // r-p block is empty
    local_col_index += FE->var_spaces[1]->dof_per_elm;

    // r-B
    local_col_index += FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[4]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[3]->dof_per_elm; trial++) {
	kij = B_phi[trial*dim]*rx_phi[test] + B_phi[trial*dim+1]*ry_phi[test];
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }
  }
    
  
  
  if (u1_phi) free(u1_phi);
  if(u1x_phi) free(u1x_phi);
  if(u1y_phi) free(u1y_phi);
  if (u2_phi) free(u2_phi);
  if(u2x_phi) free(u2x_phi);
  if(u2y_phi) free(u2y_phi);
  if (p_phi) free(p_phi);
  if(px_phi) free(px_phi);
  if(py_phi) free(py_phi);
  if (r_phi) free(r_phi);
  if(rx_phi) free(rx_phi);
  if(ry_phi) free(ry_phi);
  if (B_phi) free(B_phi);
  if(curlB_phi) free(curlB_phi);

  if(qx) free(qx);
  return;
}

void nonlin_res_MHD(REAL* bLoc,dvector* old_sol,block_fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*rhs)(REAL *,REAL *,REAL),REAL time) 
{
	
  /* Computes the local rhs for the linearized MHD system 
   * 
   *
   *
   * 
   *    INPUT:
   *            uprev               Solution at previous Newton step
   *            FE		    Finite-Element Space Struct in block form
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights
   *            dof_on_elm          Specific DOF on element
   *            v_on_elm            Vertices on element
   *            elm                 Current element
   *            time                Time if time dependent
   *
   *    OUTPUT:
   *	       	bLoc                Local RHS
   */
	
  INT i,j;

  // Mesh and FE data
  INT dof_per_elm = 0;
  for(i=0;i<FE->nspaces;i++) 
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  INT* local_dof_on_elm = NULL;
  REAL* local_uprev = NULL;

  INT dim = mesh->dim;
  
  // Loop Indices
  INT quad,test;
  // Quadrature Weights and Nodes
  REAL w; 
  REAL* qx = (REAL *) calloc(3,sizeof(REAL));
  // Stiffness Matrix Entry
  REAL kij = 0.0;

  // Data
  REAL Rs=0.0;
  REAL Rm=0.0;
  
	
  // Basis Functions and its derivatives if necessary
  // u
  REAL* u1_phi = (REAL *) calloc(FE->var_spaces[0]->dof_per_elm,sizeof(REAL));
  REAL* u1x_phi = (REAL *) calloc(FE->var_spaces[0]->dof_per_elm,sizeof(REAL));
  REAL* u1y_phi = (REAL *) calloc(FE->var_spaces[0]->dof_per_elm,sizeof(REAL));
  REAL* u2_phi = (REAL *) calloc(FE->var_spaces[1]->dof_per_elm,sizeof(REAL));
  REAL* u2x_phi = (REAL *) calloc(FE->var_spaces[1]->dof_per_elm,sizeof(REAL));
  REAL* u2y_phi = (REAL *) calloc(FE->var_spaces[1]->dof_per_elm,sizeof(REAL));
  REAL u1_prev;
  REAL u2_prev;
  REAL du1_prev[dim];
  REAL du2_prev[dim];
  // p
  REAL* p_phi = (REAL *) calloc(FE->var_spaces[2]->dof_per_elm,sizeof(REAL));
  REAL* px_phi = (REAL *) calloc(FE->var_spaces[2]->dof_per_elm,sizeof(REAL));
  REAL* py_phi = (REAL *) calloc(FE->var_spaces[2]->dof_per_elm,sizeof(REAL));
  REAL p_prev;
  // B
  REAL* B_phi = (REAL *) calloc(FE->var_spaces[3]->dof_per_elm*dim,sizeof(REAL));
  REAL* curlB_phi = (REAL *) calloc(FE->var_spaces[3]->dof_per_elm,sizeof(REAL));
  REAL B_prev[dim];
  REAL curlB_prev;
  // r
  REAL* r_phi = (REAL *) calloc(FE->var_spaces[4]->dof_per_elm,sizeof(REAL));
  REAL* rx_phi = (REAL *) calloc(FE->var_spaces[4]->dof_per_elm,sizeof(REAL));
  REAL* ry_phi = (REAL *) calloc(FE->var_spaces[4]->dof_per_elm,sizeof(REAL));
  REAL dr_prev[dim];

  REAL rhs_val[6];

  INT local_row_index;
      
  //  Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {        
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];
    
    //  Get the Basis Functions and old solutions at each quadrature node
    // u1
    PX_H1_basis(u1_phi,u1x_phi,u1y_phi,NULL,qx,dof_on_elm,FE->var_spaces[0]->FEtype,mesh);
    FE_Interpolation(&u1_prev,old_sol->val,qx,dof_on_elm,v_on_elm,FE->var_spaces[0],mesh,FE->var_spaces[0]->ndof,1);
    FE_DerivativeInterpolation(du1_prev,old_sol->val,qx,dof_on_elm,v_on_elm,FE->var_spaces[0],mesh,FE->var_spaces[0]->ndof,1);
    local_uprev = old_sol->val + FE->var_spaces[0]->ndof;
    local_dof_on_elm = dof_on_elm + FE->var_spaces[0]->dof_per_elm;
    // u2
    PX_H1_basis(u2_phi,u2x_phi,u2y_phi,NULL,qx,local_dof_on_elm,FE->var_spaces[1]->FEtype,mesh);
    FE_Interpolation(&u2_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh,FE->var_spaces[1]->ndof,1);
    FE_DerivativeInterpolation(du1_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh,FE->var_spaces[1]->ndof,1);
    local_uprev = local_uprev + FE->var_spaces[1]->ndof;
    local_dof_on_elm = local_dof_on_elm + FE->var_spaces[1]->dof_per_elm;
    // p
    PX_H1_basis(p_phi,px_phi,py_phi,NULL,qx,local_dof_on_elm,FE->var_spaces[2]->FEtype,mesh);
    FE_Interpolation(&p_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[2],mesh,FE->var_spaces[2]->ndof,1);
    local_uprev = local_uprev + FE->var_spaces[2]->ndof;
    local_dof_on_elm = local_dof_on_elm + FE->var_spaces[2]->dof_per_elm;
    // B
    ned_basis(B_phi,curlB_phi,qx,v_on_elm,local_dof_on_elm,mesh);
    FE_Interpolation(B_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[3],mesh,FE->var_spaces[3]->ndof,1);
    FE_DerivativeInterpolation(&curlB_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[3],mesh,FE->var_spaces[3]->ndof,1);
    local_uprev = local_uprev + FE->var_spaces[3]->ndof;
    local_dof_on_elm = local_dof_on_elm + FE->var_spaces[3]->dof_per_elm;
    // r
    PX_H1_basis(r_phi,rx_phi,ry_phi,NULL,qx,local_dof_on_elm,FE->var_spaces[4]->FEtype,mesh);
    FE_DerivativeInterpolation(dr_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[4],mesh,FE->var_spaces[4]->ndof,1);
    // Data
    Rs = 1;
    Rm = 1;

    (*rhs)(rhs_val,qx,0.0);

    // u1 block
    local_row_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++) {
      kij = rhs_val[0]*u1_phi[test] - (1.0/Rs)*(du1_prev[0]*u1x_phi[test] + 0.5*(du2_prev[0]*u1y_phi[test] + du1_prev[1]*u1y_phi[test]))
	+ p_prev*u1x_phi[test] - (u1_prev*du1_prev[0] + u2_prev*du1_prev[1] 
				  -curlB_prev*B_prev[1])*u1_phi[test];
      bLoc[(local_row_index+test)] += w*kij;
    }

    // u2 block
    local_row_index += FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++) {
      kij = rhs_val[1]*u2_phi[test] - (1.0/Rs)*(du2_prev[1]*u2y_phi[test] + 0.5*(du2_prev[0]*u2x_phi[test] + du1_prev[1]*u2x_phi[test]))
	+ p_prev*u2y_phi[test] - (u1_prev*du2_prev[0] + u2_prev*du2_prev[1] 
				  + curlB_prev*B_prev[0])*u2_phi[test];
      bLoc[(local_row_index+test)] += w*kij;
    }

    // p block
    local_row_index += FE->var_spaces[1]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++) {
      kij = rhs_val[2]*p_phi[test] - (du1_prev[0]+du2_prev[1])*p_phi[test];
      bLoc[(local_row_index+test)] += w*kij;
    }

    // B block
    local_row_index += FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++) {
      kij = rhs_val[3]*B_phi[test*dim] + rhs_val[4]*B_phi[test*dim+1] 
	- (1.0/Rm)*curlB_prev*curlB_phi[test] + (u1_prev*B_prev[1] - u2_prev*B_prev[0])*curlB_phi[test]
	+ dr_prev[0]*B_phi[test*dim] + dr_prev[1]*B_phi[test*dim+1];
      bLoc[(local_row_index+test)] += w*kij;
    }
    
    // r block
    local_row_index += FE->var_spaces[3]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[4]->dof_per_elm;test++) {
      kij = rhs_val[5]*r_phi[test] - B_prev[0]*rx_phi[test] - B_prev[1]*ry_phi[test];
      bLoc[(local_row_index+test)] += w*kij;
    }
      
  }
    
  
  
  if (u1_phi) free(u1_phi);
  if(u1x_phi) free(u1x_phi);
  if(u1y_phi) free(u1y_phi);
  if (u2_phi) free(u2_phi);
  if(u2x_phi) free(u2x_phi);
  if(u2y_phi) free(u2y_phi);
  if (p_phi) free(p_phi);
  if(px_phi) free(px_phi);
  if(py_phi) free(py_phi);
  if (r_phi) free(r_phi);
  if(rx_phi) free(rx_phi);
  if(ry_phi) free(ry_phi);
  if (B_phi) free(B_phi);
  if(curlB_phi) free(curlB_phi);

  if(qx) free(qx);
  return;
}
/******************************************************************************************************/


