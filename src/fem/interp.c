/*
 *  interp.c
 *  
 *  Created by James Adler and Xiaozhe Hu on 2/1/15.
 *  Copyright 2015_JXLcode__. All rights reserved.
 *
 */

#include "hazmat.h"

/****************************************************************************************************************************/
// INTERPOLATION AND EVALUATION ROUTINES for Finite Element Basis Functions
/****************************************************************************************************************************/

/****************************************************************************************************************************/
void FE_Interpolation(REAL* val,REAL *u,REAL x,REAL y,REAL z,INT *dof_on_elm,INT *v_on_elm,fespace *FE,trimesh *mesh,INT ndof,INT nun)
{

/* Interpolate a finite-element approximation to any other point in the given element using the given type of elements 
   *    INPUT:
   *		       u       Approximation to interpolate
   *               x,y,z       Coordinates where to compute value
   *          dof_on_elm       DOF belonging to particular element
   *          v_on_elm         Vertices belonging to particular element
   *                  FE       FE Space struct
   *                mesh       Mesh struct
   *		    ndof       Total DOF for each unknown
   *		     nun       Number of unknowns in u (u1,u2,etc?) 1 is a scalar...
   *    OUTPUT:
   *          val	       Value of approximation at given values
   *
   */
	
	
  INT i,nd,j;

  // Get FE and Mesh data
  INT dof_per_elm = FE->dof_per_elm;
  INT FEtype = FE->FEtype;
  INT dim = mesh->dim;

  // Basis Functions and its derivatives if necessary
  REAL* phi=NULL;
  REAL* dphix=NULL;
  REAL* dphiy=NULL;
  REAL* dphiz=NULL;

  if(FEtype>0) { // Lagrange Elements
    phi = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    dphix = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    dphiy = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    if(dim==3) dphiz = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    PX_H1_basis(phi,dphix,dphiy,dphiz,x,y,z,dof_on_elm,FEtype,mesh);
    REAL coef;

    for (i=0; i<nun; i++) {
      coef = 0.0;
      
      for (j=0; j<dof_per_elm; j++) {
	nd = i*ndof + dof_on_elm[j] - 1;
	coef = coef + u[nd]*phi[j];
      }
      val[i] = coef;
    }
  } else if (FEtype==-1) { // Nedelec
    REAL coef1,coef2,coef3;
    INT edge;
    phi = (REAL *) calloc(dof_per_elm*dim,sizeof(REAL));
    if(dim==2) {
      dphix = (REAL *) calloc(dof_per_elm,sizeof(REAL)); // Curl of basis function
    } else if (dim==3) {
      dphix = (REAL *) calloc(dof_per_elm*dim,sizeof(REAL)); // Curl of basis function
    } else {
      baddimension();
    }
    
    ned_basis(phi,dphix,x,y,z,v_on_elm,dof_on_elm,mesh);
    	
    coef1 = 0.0;
    coef2 = 0.0;
    coef3 = 0.0;
		
    if (dim==2) {
      for (j=0; j<dof_per_elm; j++) {
	edge = dof_on_elm[j]-1;
	coef1 = coef1 + u[edge]*phi[j*dim+0];
	coef2 = coef2 + u[edge]*phi[j*dim+1];
      }
      val[0] = coef1;
      val[1] = coef2;
    } else if (dim==3) {
      for (j=0; j<dof_per_elm; j++) {
	edge = dof_on_elm[j]-1;
	coef1 = coef1 + u[edge]*phi[j*dim+0];
	coef2 = coef2 + u[edge]*phi[j*dim+1];
	coef3 = coef3 + u[edge]*phi[j*dim+2];
      }
      val[0] = coef1;
      val[1] = coef2;
      val[2] = coef3;
    }
  } else if (FEtype==-2) { // Raviart-Thomas
    REAL coef1,coef2,coef3;
    INT face;
    phi = (REAL *) calloc(dof_per_elm*dim,sizeof(REAL));
    dphix = (REAL *) calloc(dof_per_elm,sizeof(REAL)); // Divergence of element
    
    rt_basis(phi,dphix,x,y,z,v_on_elm,dof_on_elm,mesh);
    	
    coef1 = 0.0;
    coef2 = 0.0;
    coef3 = 0.0;
		
    if (dim==2) {
      for (j=0; j<dof_per_elm; j++) {
	face = dof_on_elm[j]-1;
	coef1 = coef1 + u[face]*phi[j*dim+0];
	coef2 = coef2 + u[face]*phi[j*dim+1];
      }
      val[0] = coef1;
      val[1] = coef2;
    } else if (dim==3) {
      for (j=0; j<dof_per_elm; j++) {
	face = dof_on_elm[j]-1;
	coef1 = coef1 + u[face]*phi[j*dim+0];
	coef2 = coef2 + u[face]*phi[j*dim+1];
	coef3 = coef3 + u[face]*phi[j*dim+2];
      }
      val[0] = coef1;
      val[1] = coef2;
      val[2] = coef3;
    }
  }
    
    if (phi) free(phi);
    if(dphix) free(dphix);
    if(dphiy) free(dphiy);
    if(dphiz) free(dphiz);
    return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
void FE_DerivativeInterpolation(REAL* val,REAL *u,REAL x,REAL y,REAL z,INT *dof_on_elm,INT *v_on_elm,fespace *FE,trimesh *mesh,INT ndof,INT nun)
{

/* Interpolate the "derivative" of a finite-element approximation to any other point in the given element using the given type of elements 
 * Note that for Lagrange Elements this means the Gradient, grad u, for Nedelec it means the Curl, curl u, and for RT it is the Divergence, div u.
   *    INPUT:
   *		       u       Approximation to interpolate
   *               x,y,z       Coordinates where to compute value
   *          dof_on_elm       DOF belonging to particular element
   *          v_on_elm         Vertices belonging to particular element
   *                  FE       FE Space struct
   *                mesh       Mesh struct
   *		    ndof       Total DOF for each unknown
   *		     nun       Number of unknowns in u (u1,u2,etc?) 1 is a scalar...
   *    OUTPUT:
   *          val	       Value of derivative of approximation at given values
   *
   */
	
	
  INT i,nd,j;

  // Get FE and Mesh data
  INT dof_per_elm = FE->dof_per_elm;
  INT FEtype = FE->FEtype;
  INT dim = mesh->dim;

  // Basis Functions and its derivatives if necessary
  REAL* phi=NULL;
  REAL* dphix=NULL;
  REAL* dphiy=NULL;
  REAL* dphiz=NULL;
  REAL* coef=NULL;

  if(FEtype>0) { // Lagrange Elements
    phi = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    dphix = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    dphiy = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    if(dim==3) dphiz = (REAL *) calloc(dof_per_elm,sizeof(REAL));
    PX_H1_basis(phi,dphix,dphiy,dphiz,x,y,z,dof_on_elm,FEtype,mesh);
    coef = (REAL *) calloc(dim,sizeof(REAL));

    for (i=0; i<nun; i++) {
      for(j=0;j<dim;j++) coef[0] = 0.0;
      
      for (j=0; j<dof_per_elm; j++) {
	nd = i*ndof + dof_on_elm[j] - 1;
	coef[0] = coef[0] + u[nd]*dphix[j];
	coef[1] = coef[1] + u[nd]*dphiy[j];
	if(dim==3) coef[2] = coef[2] + u[nd]*dphiz[j];
      }
      val[i] = coef[0];
      val[nun+i] = coef[1];
      if(dim==3) val[2*nun+i] = coef[2];
    }
  } else if (FEtype==-1) { // Nedelec
    INT edge;
    phi = (REAL *) calloc(dof_per_elm*dim,sizeof(REAL));
    if(dim==2) {
      dphix = (REAL *) calloc(dof_per_elm,sizeof(REAL)); // Curl of basis function
    } else if (dim==3) {
      dphix = (REAL *) calloc(dof_per_elm*dim,sizeof(REAL)); // Curl of basis function
    } else {
      baddimension();
    }
    ned_basis(phi,dphix,x,y,z,v_on_elm,dof_on_elm,mesh);
    	
    coef = (REAL *) calloc(dim,sizeof(REAL));
    for(j=0;j<dim;j++) coef[0] = 0.0;
		
    if (dim==2) {
      for (j=0; j<dof_per_elm; j++) {
	edge = dof_on_elm[j]-1;
	coef[0] = coef[0] + u[edge]*dphix[j];
      }
      val[0] = coef[0];
    } else if (dim==3) {
      for (j=0; j<dof_per_elm; j++) {
	edge = dof_on_elm[j]-1;
	coef[0] = coef[0] + u[edge]*dphix[j*dim+0];
	coef[1] = coef[1] + u[edge]*dphix[j*dim+1];
	coef[2] = coef[2] + u[edge]*dphix[j*dim+2];
      }
      val[0] = coef[0];
      val[1] = coef[1];
      val[2] = coef[2];
    }
  } else if (FEtype==-2) { // Raviart-Thomas
    INT face;
    phi = (REAL *) calloc(dof_per_elm*dim,sizeof(REAL));
    dphix = (REAL *) calloc(dof_per_elm,sizeof(REAL)); // Divergence of element
    
    rt_basis(phi,dphix,x,y,z,v_on_elm,dof_on_elm,mesh);
    	
    coef = (REAL *) calloc(1,sizeof(REAL));
    coef[0] = 0.0;
		
    for (j=0; j<dof_per_elm; j++) {
      face = dof_on_elm[j]-1;
      coef[0] = coef[0] + u[face]*dphix[j];
    }
    val[0] = coef[0];
  }
    
    if (phi) free(phi);
    if(dphix) free(dphix);
    if(dphiy) free(dphiy);
    if(dphiz) free(dphiz);
    if(coef) free(coef);
    return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
void FE_Evaluate(REAL* val,void (*expr)(REAL *,REAL *,REAL),fespace *FE,trimesh *mesh,REAL time)
{

/* Evaluate a given analytic function on the finite-element space given
   *    INPUT:
   *		    expr       Function call to analytic expression expr(FE approx, x values, time)
   *                  FE       FE Space struct
   *                mesh       Mesh struct
   *                time       Time to evaluate function if time-dependent
   *    OUTPUT:
   *          val	       FE approximation of function on fespace
   *
   */

  int i,j;
  REAL* x = (REAL *) calloc(mesh->dim,sizeof(REAL));
  REAL* valx = NULL;
  INT dim = mesh->dim;
  INT FEtype = FE->FEtype;
  
  if(FEtype>0) { // Lagrange Elements u[dof] = u[x_i}
    for(i=0;i<FE->ndof;i++) {
      valx = (REAL *) calloc(1,sizeof(REAL));
      x[0] = FE->cdof->x[i];
      x[1] = FE->cdof->y[i];
      if(dim==3) x[2] = FE->cdof->z[i];
      (*expr)(valx,x,time);
      val[i] = valx[0];
    }
  } else if (FEtype==-1) { // Nedelec u[dof] = (1/elen) \int_edge u*t_edge
    for(i=0;i<FE->ndof;i++) {
      valx = (REAL *) calloc(dim,sizeof(REAL));
      x[0] = mesh->ed_mid[i*dim];
      x[1] = mesh->ed_mid[i*dim+1];
      if(dim==3) x[2] = mesh->ed_mid[i*dim+1];
      (*expr)(valx,x,time);
      val[i] = 0.0;
      for(j=0;j<dim;j++) val[i]+=mesh->ed_tau[i*dim+j]*valx[j];
    }
  } else if (FEtype==-2) { // Raviart-Thomas u[dof] = 1/farea \int_face u*n_face
    for(i=0;i<FE->ndof;i++) {
      valx = (REAL *) calloc(dim,sizeof(REAL));
      x[0] = mesh->f_mid[i*dim];
      x[1] = mesh->f_mid[i*dim+1];
      if(dim==3) x[2] = mesh->f_mid[i*dim+1];
      (*expr)(valx,x,time);
      val[i] = 0.0;
      for(j=0;j<dim;j++) val[i]+=mesh->f_norm[i*dim+j]*valx[j];
    }
  }
  
  if (x) free(x);
  if(valx) free(valx);
  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
REAL FE_Evaluate_DOF(void (*expr)(REAL *,REAL *,REAL),fespace *FE,trimesh *mesh,REAL time,INT DOF)
{

/* Evaluate a given analytic function on the specific degree of freedom of the finite-element space given
   *    INPUT:
   *		    expr       Function call to analytic expression expr(FE approx, x values, time)
   *                  FE       FE Space struct
   *                mesh       Mesh struct
   *                time       Time to evaluate function if time-dependent
   *                 DOF       DOF index to evaluate (start at 0)
   *    OUTPUT:
   *          val	       FE approximation of function on fespace at DOF
   *
   */

  INT j;
  REAL* x = (REAL *) calloc(mesh->dim,sizeof(REAL));
  REAL* valx = NULL;
  INT dim = mesh->dim;
  INT FEtype = FE->FEtype;
  REAL val=-666e+00;
  
  if(FEtype>0) { // Lagrange Elements u[dof] = u[x_i}
    valx = (REAL *) calloc(1,sizeof(REAL));
    x[0] = FE->cdof->x[DOF];
    x[1] = FE->cdof->y[DOF];
    if(dim==3) x[2] = FE->cdof->z[DOF];
    (*expr)(valx,x,time);
    val = valx[0];
  } else if (FEtype==-1) { // Nedelec u[dof] = (1/elen) \int_edge u*t_edge
    valx = (REAL *) calloc(dim,sizeof(REAL));
    x[0] = mesh->ed_mid[DOF*dim];
    x[1] = mesh->ed_mid[DOF*dim+1];
    if(dim==3) x[2] = mesh->ed_mid[DOF*dim+2];
    (*expr)(valx,x,time);
    val = 0.0;
    for(j=0;j<dim;j++) val+=mesh->ed_tau[DOF*dim+j]*valx[j];
  } else if (FEtype==-2) { // Raviart-Thomas u[dof] = 1/farea \int_face u*n_face
    valx = (REAL *) calloc(dim,sizeof(REAL));
    x[0] = mesh->f_mid[DOF*dim];
    x[1] = mesh->f_mid[DOF*dim+1];
    if(dim==3) x[2] = mesh->f_mid[DOF*dim+2];
    (*expr)(valx,x,time);
    val = 0.0;
    for(j=0;j<dim;j++) val+=mesh->f_norm[DOF*dim+j]*valx[j];
  }
  
  if (x) free(x);
  if(valx) free(valx);

  return val;
}
/****************************************************************************************************************************/
