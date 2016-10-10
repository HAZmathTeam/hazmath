/*! \file src/fem/error.c
 *  
 *  Created by James Adler and Xiaozhe Hu on 5/13/15.
 *  Copyright 2015_JXLcode__. All rights reserved.
 *
 */

#include "hazmat.h"

/****************************************************************************************************************************/
// Error and Norm ROUTINES
/****************************************************************************************************************************/

/***************************************************************************/
REAL L2norm(REAL *u,fespace *FE,trimesh *mesh,qcoordinates *cq)
{

  /* Computes the L2 Norm of a FE approximation using the mass matrix assembly any type of element.
   *
   * Input:		u	Numerical Solution at DOF
   *                    FE      fespace struct
   *                    mesh    mesh struct
   *                    cq      Quadrature Nodes
   *
   * Output:		norm	L2 Norm
   *
   */
		
  INT i,j,k,rowa,rowb,jcntr;
  REAL sum = 0.0;

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* MLoc = calloc(local_size,sizeof(REAL));
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  /* Loop over all Elements */
  for (i=0; i<FE->nelm; i++) {

    // Zero out local matrices
    for (j=0; j<local_size; j++) MLoc[j] = 0.0;
    
    // Find DOF for given Element
    rowa = FE->el_dof->IA[i]-1;
    rowb = FE->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm[jcntr] = FE->el_dof->JA[j];
      jcntr++;
    }

    //Find Nodes for given Element if not H1 elements
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr=0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }
    		
    // Compute Local Stiffness Matrix for given Element
    assemble_mass_local(MLoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,constant_coeff_scal,1.0);

    for(j=0;j<dof_per_elm;j++) {
      for(k=0;k<dof_per_elm;k++) {
	sum+=u[dof_on_elm[j]-1]*MLoc[j*dof_per_elm+k]*u[dof_on_elm[k]-1];
      }
    }
  }

  if(MLoc) free(MLoc);
  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);

  return sqrt(sum);
}
/*******************************************************************************************************************************************************/

/***************************************************************************/
REAL L2_InnerProduct(REAL *u,REAL *v,fespace *FE,trimesh *mesh,qcoordinates *cq)
{

  /* Computes the L2 inner product of two FE approximations using the mass matrix assembly any type of element.
   *
   * Input:		u,v	Numerical Solutions at DOF
   *                    FE      fespace struct
   *                    mesh    mesh struct
   *                    cq      Quadrature Nodes
   *
   * Output:		product	L2 Inner Product of u and v, <u,v>
   *
   */
		
  INT i,j,k,rowa,rowb,jcntr;
  REAL sum = 0.0;

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* MLoc = calloc(local_size,sizeof(REAL));
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  /* Loop over all Elements */
  for (i=0; i<FE->nelm; i++) {

    // Zero out local matrices
    for (j=0; j<local_size; j++) MLoc[j] = 0.0;
    
    // Find Nodes for given Element
    rowa = FE->el_dof->IA[i]-1;
    rowb = FE->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm[jcntr] = FE->el_dof->JA[j];
      jcntr++;
    }

    //Find Nodes for given Element if not H1 elements
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr=0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }
    		
    // Compute Local Stiffness Matrix for given Element
    assemble_mass_local(MLoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,constant_coeff_scal,1.0);

    for(j=0;j<dof_per_elm;j++) {
      for(k=0;k<dof_per_elm;k++) {
	sum+=v[dof_on_elm[j]-1]*MLoc[j*dof_per_elm+k]*u[dof_on_elm[k]-1];
      }
    }
  }

  if(MLoc) free(MLoc);
  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);

  return sum;
}
/*******************************************************************************************************************************************************/

/***************************************************************************/
REAL L2error(REAL *u,void (*truesol)(REAL *,REAL *,REAL),fespace *FE,trimesh *mesh,qcoordinates *cq,REAL time)
{

  /* Computes the L2 Norm of the error of a FE approximation and a true solution given by a function using quadrature for any type of element.
   *
   * Input:		u	Numerical Solution at DOF
   *                  truesol   Function to get true solution at any given point
   *                    FE      fespace struct
   *                    mesh    mesh struct
   *                    cq      Quadrature Nodes
   *                  time      Time to evaluate solution
   *
   * Output:		error	L2 Error
   *
   */

  INT dim = mesh->dim;

  // Quadrature Weights and Nodes
  REAL w; 
  REAL* qx = (REAL *) calloc(mesh->dim,sizeof(REAL));
  
  // FE Stuff
  INT FEtype = FE->FEtype;
  INT elm,quad,j,rowa,rowb,jcntr;
  REAL sum = 0.0;
  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  // True Solution and FE Solution at Quadrature Nodes
  INT ncomp = 0;
  if(FEtype>=0) { // Lagrange Elements
    ncomp = 1;
  } else { // Vector Elements
    ncomp = dim;
  }
  REAL* val_true = (REAL *) calloc(ncomp,sizeof(REAL));
  REAL* val_sol = (REAL *) calloc(ncomp,sizeof(REAL));

  /* Loop over all Elements */
  for (elm=0; elm<FE->nelm; elm++) {
    // Loop over quadrature nodes on element
    for (quad=0;quad<cq->nq_per_elm;quad++) {        
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];

      // Get True Solution at Quadrature Nodes
      (*truesol)(val_true,qx,time);
    
      // Find DOF for given Element
      rowa = FE->el_dof->IA[elm]-1;
      rowb = FE->el_dof->IA[elm+1]-1;
      jcntr = 0;
      for (j=rowa; j<rowb; j++) {
	dof_on_elm[jcntr] = FE->el_dof->JA[j];
	jcntr++;
      }

      //Find Vertices for given Element if not H1 elements
      rowa = mesh->el_v->IA[elm]-1;
      rowb = mesh->el_v->IA[elm+1]-1;
      jcntr=0;
      for (j=rowa; j<rowb; j++) {
	v_on_elm[jcntr] = mesh->el_v->JA[j];
	jcntr++;
      }

      // Interpolate FE solution to quadrature point
      FE_Interpolation(val_sol,u,qx,dof_on_elm,v_on_elm,FE,mesh,FE->ndof,1);

      // Compute Error on Element
      for(j=0;j<ncomp;j++) {
	sum+=w*(ABS(val_sol[j] - val_true[j]))*(ABS(val_sol[j] - val_true[j]));
      }
    }
  }


  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(qx) free(qx);
  if(val_true) free(val_true);
  if(val_sol) free(val_sol);

  return sqrt(sum);
}
/*******************************************************************************************************************************************************/

/***************************************************************************/
REAL L2error_mass(REAL *u,void (*truesol)(REAL *,REAL *,REAL),fespace *FE,trimesh *mesh,qcoordinates *cq,REAL time)
{

  /* Computes the L2 Norm of the error of a FE approximation and a true solution given by a function using the mass matrix assembly any type of element.
   *
   * Input:		u	Numerical Solution at DOF
   *                  truesol   Function to get true solution at any given point
   *                    FE      fespace struct
   *                    mesh    mesh struct
   *                    cq      Quadrature Nodes
   *                  time      Time to evaluate solution
   *
   * Output:		error	L2 Error
   *
   */
		
  INT i,j,k,rowa,rowb,jcntr;
  REAL sum = 0.0;
  REAL utk,utj,erk,erj;

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* MLoc = calloc(local_size,sizeof(REAL));
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  /* Loop over all Elements */
  for (i=0; i<FE->nelm; i++) {

    // Zero out local matrices
    for (j=0; j<local_size; j++) MLoc[j] = 0.0;
    
    // Find Nodes for given Element
    rowa = FE->el_dof->IA[i]-1;
    rowb = FE->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm[jcntr] = FE->el_dof->JA[j];
      jcntr++;
    }

    //Find Nodes for given Element if not H1 elements
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr=0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }
    		
    // Compute Local Stiffness Matrix for given Element
    assemble_mass_local(MLoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,constant_coeff_scal,1.0);

    for(j=0;j<dof_per_elm;j++) {
      for(k=0;k<dof_per_elm;k++) {
	utj = FE_Evaluate_DOF(truesol,FE,mesh,time,dof_on_elm[j]-1);
	utk = FE_Evaluate_DOF(truesol,FE,mesh,time,dof_on_elm[k]-1);
	erj = ABS(utj - u[dof_on_elm[j]-1]);
	erk = ABS(utk - u[dof_on_elm[k]-1]);
	sum+=erj*MLoc[j*dof_per_elm+k]*erk;
      }
    }
  }

  if(MLoc) free(MLoc);
  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);

  return sqrt(sum);
}
/*******************************************************************************************************************************************************/

/***************************************************************************/
REAL HDseminorm(REAL *u,fespace *FE,trimesh *mesh,qcoordinates *cq)
{
  /* Computes the H(D) semi-Norm of a FE approximation using the <Du,Dv> matrix assembly for any type of element.
   *          Nodal - <grad u, grad u> - |u|_1
   *          RT - <div u, div u>      - |u|_(H(div))
   *          Nedelec - <curl u, curl u>  - |u|_(H(curl))
   *
   * Input:		u	Numerical Solution at DOF
   *                    FE      fespace struct
   *                    mesh    mesh struct
   *                    cq      Quadrature Nodes
   *
   * Output:		norm	Semi-Norm
   *
   */
		
  INT i,j,k,rowa,rowb,jcntr;
  REAL sum = 0.0;
  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* ALoc = calloc(local_size,sizeof(REAL));
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  /* Loop over all Elements */
  for (i=0; i<FE->nelm; i++) {
		
    // Zero out local matrices
    for (j=0; j<local_size; j++) ALoc[j] = 0.0;
		
    // Find DOF for given Element
    rowa = FE->el_dof->IA[i]-1;
    rowb = FE->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm[jcntr] = FE->el_dof->JA[j];
      jcntr++;
    }

    //Find Nodes for given Element if not H1 elements
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr=0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }
    		
    // Compute Local Stiffness Matrix for given Element
    assemble_DuDv_local(ALoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,constant_coeff_scal,1.0);

    for(j=0;j<dof_per_elm;j++) {
      for(k=0;k<dof_per_elm;k++) {
	sum+=u[dof_on_elm[j]-1]*ALoc[j*dof_per_elm+k]*u[dof_on_elm[k]-1];
      }
    }
  }

  if(ALoc) free(ALoc);
  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);

  if(sum<0.0) {
    printf("Your H1 Semi Norm Squared is negative!  Outputting the square itself\n");
    return sum;
  } else {
    return sqrt(sum);
  }

}
/*******************************************************************************************************************************************************/

/***************************************************************************/
REAL HDsemierror(REAL *u,void (*D_truesol)(REAL *,REAL *,REAL),fespace *FE,trimesh *mesh,qcoordinates *cq,REAL time)
{

  /* Computes the H(D) semi-Norm of the error of a FE approximation and a true solution given by a function using quadrature for any type of element.
   *          Nodal - <grad u, grad u> - |u|_1
   *          RT - <div u, div u>      - |u|_(H(div))
   *          Nedelec - <curl u, curl u>  - |u|_(H(curl))
   *
   * Input:		u	Numerical Solution at DOF
   *                  truesol   Function to get derivative of true solution at any given point
   *                    FE      fespace struct
   *                    mesh    mesh struct
   *                    cq      Quadrature Nodes
   *                  time      Time to evaluate solution
   *
   * Output:		error	Semi-Norm Error
   *
   */

  INT dim = mesh->dim;

  // Quadrature Weights and Nodes
  REAL w; 
  REAL* qx = (REAL *) calloc(mesh->dim,sizeof(REAL));
  
  // FE Stuff
  INT FEtype = FE->FEtype;
  INT elm,quad,j,rowa,rowb,jcntr;
  REAL sum = 0.0;
  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  // Derivative of True Solution and FE Solution at Quadrature Nodes
  INT ncomp = 0;
  if(FEtype>=0) { // Lagrange Elements -> Grad
    ncomp = dim;
  } else if(FEtype==-1 && dim==3) { // Nedelec Elements in 3D -> Curl is a Vector
    ncomp = dim;
  } else { // 2D Nedelec -> Curl is a scalar or RT -> Div is a scalar
    ncomp = 1;
  }
  REAL* val_true = (REAL *) calloc(ncomp,sizeof(REAL));
  REAL* val_sol = (REAL *) calloc(ncomp,sizeof(REAL));

  /* Loop over all Elements */
  for (elm=0; elm<FE->nelm; elm++) {
    // Loop over quadrature nodes on element
    for (quad=0;quad<cq->nq_per_elm;quad++) {        
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];

      // Get True Solution at Quadrature Nodes
      (*D_truesol)(val_true,qx,time);
    
      // Find DOF for given Element
      rowa = FE->el_dof->IA[elm]-1;
      rowb = FE->el_dof->IA[elm+1]-1;
      jcntr = 0;
      for (j=rowa; j<rowb; j++) {
	dof_on_elm[jcntr] = FE->el_dof->JA[j];
	jcntr++;
      }

      //Find Vertices for given Element if not H1 elements
      rowa = mesh->el_v->IA[elm]-1;
      rowb = mesh->el_v->IA[elm+1]-1;
      jcntr=0;
      for (j=rowa; j<rowb; j++) {
	v_on_elm[jcntr] = mesh->el_v->JA[j];
	jcntr++;
      }

      // Interpolate FE solution to quadrature point
      FE_DerivativeInterpolation(val_sol,u,qx,dof_on_elm,v_on_elm,FE,mesh,FE->ndof,1);

      // Compute Error on Element
      for(j=0;j<ncomp;j++) {
	sum+=w*(ABS(val_sol[j] - val_true[j]))*(ABS(val_sol[j] - val_true[j]));
      }
    }
  }


  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(qx) free(qx);
  if(val_true) free(val_true);
  if(val_sol) free(val_sol);

  return sqrt(sum);
}
/*******************************************************************************************************************************************************/

/***************************************************************************/
REAL HDsemierror_stiff(REAL *u,void (*truesol)(REAL *,REAL *,REAL),fespace *FE,trimesh *mesh,qcoordinates *cq,REAL time)
{
  /* Computes the H(D) semi-Norm of the error of a FE approximation and a true solution given by a function using the <Du,Dv> matrix assembly for any type of element.
   *          Nodal - <grad u, grad u> - |u|_1
   *          RT - <div u, div u>      - |u|_(H(div))
   *          Nedelec - <curl u, curl u>  - |u|_(H(curl))
   *
   * Input:		u	Numerical Solution at DOF
   *                  truesol   Function to get true solution at any given point
   *                    FE      fespace struct
   *                    mesh    mesh struct
   *                    cq      Quadrature Nodes
   *                  time      Time to evaluate solution
   *
   * Output:		error	Semi-Norm Error
   *
   */
		
  INT i,j,k,rowa,rowb,jcntr;
  REAL sum = 0.0;
  REAL utk,utj,erk,erj;

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* ALoc = calloc(local_size,sizeof(REAL));
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  /* Loop over all Elements */
  for (i=0; i<FE->nelm; i++) {
		
    // Zero out local matrices
    for (j=0; j<local_size; j++) ALoc[j] = 0.0;
		
    // Find DOF for given Element
    rowa = FE->el_dof->IA[i]-1;
    rowb = FE->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm[jcntr] = FE->el_dof->JA[j];
      jcntr++;
    }

    //Find Nodes for given Element if not H1 elements
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr=0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }
    		
    // Compute Local Stiffness Matrix for given Element
    assemble_DuDv_local(ALoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,constant_coeff_scal,1.0);

    for(j=0;j<dof_per_elm;j++) {
      for(k=0;k<dof_per_elm;k++) {
	utj = FE_Evaluate_DOF(truesol,FE,mesh,time,dof_on_elm[j]-1);
	utk = FE_Evaluate_DOF(truesol,FE,mesh,time,dof_on_elm[k]-1);
	erj = ABS(utj - u[dof_on_elm[j]-1]);
	erk = ABS(utk - u[dof_on_elm[k]-1]);
	sum+=erj*ALoc[j*dof_per_elm+k]*erk;
      }
    }
  }

  if(ALoc) free(ALoc);
  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);

  if(sum<0.0) {
    printf("Your H1 Semi Norm Error Squared is negative!  Outputting the square itself\n");
    return sum;
  } else {
    return sqrt(sum);
  }

}
/*******************************************************************************************************************************************************/

/***************************************************************************/
REAL HDnorm(REAL *u,fespace *FE,trimesh *mesh,qcoordinates *cq)
{
  /* Computes the H(D) norm of a FE approximation using the <u,v> and <Du,Dv> matrix assembly for any type of element.
   *          Nodal - <u,u> + <grad u, grad u> - ||u||_1
   *          RT - <u,u> + <div u, div u>      - ||u||_(H(div))
   *          Nedelec - <u,u> + <curl u, curl u>  - ||u||_(H(curl))
   *
   * Input:		u	Numerical Solution at DOF
   *                    FE      fespace struct
   *                    mesh    mesh struct
   *                    cq      Quadrature Nodes
   *
   * Output:		norm	Norm
   *
   */
		

  REAL sumL2 = L2norm(u,FE,mesh,cq);
  REAL sumSemi = HDseminorm(u,FE,mesh,cq);

  return sqrt(sumL2*sumL2 + sumSemi*sumSemi);

}
/*******************************************************************************************************************************************************/

/***************************************************************************/
REAL HDerror(REAL *u,void (*truesol)(REAL *,REAL *,REAL),void (*D_truesol)(REAL *,REAL *,REAL),fespace *FE,trimesh *mesh,qcoordinates *cq,REAL time)
{
  /* Computes the H(D) Norm of the error of a FE approximation and a true solution given by a function using the <u,v> and <Du,Dv> matrix assembly for any type of element.
   *          Nodal - <u,u> + <grad u, grad u> - ||u||_1
   *          RT - <u,u> + <div u, div u>      - ||u||_(H(div))
   *          Nedelec - <u,u> + <curl u, curl u>  - ||u||_(H(curl))
   *
   * Input:		u	Numerical Solution at DOF
   *                  truesol   Function to get true solution at any given point
   *                  D_truesol Derivative of Function to get true solution at any given point
   *                    FE      fespace struct
   *                    mesh    mesh struct
   *                    cq      Quadrature Nodes
   *                  time      Time to evaluate solution
   *
   * Output:		error	Norm Error
   *
   */
		
  REAL sumL2 = L2error(u,truesol,FE,mesh,cq,time);
  REAL sumSemi = HDsemierror(u,D_truesol,FE,mesh,cq,time);

  return sqrt(sumL2*sumL2 + sumSemi*sumSemi);

}
/*******************************************************************************************************************************************************/
