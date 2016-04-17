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

