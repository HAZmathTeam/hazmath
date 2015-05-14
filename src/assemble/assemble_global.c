/*
 *  assemble_global.c   
 *  
 *  Created by James Adler and Xiaozhe Hu on 4/22/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

/* This code will build global stiffness matrices for various PDE systems */

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
void assemble_global(dCSRmat* A,dvector *b,void (*local_assembly)(REAL *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),void (*bc)(REAL *,REAL *,REAL),void (*coeff)(REAL *,REAL *,REAL),REAL time) 
{
	
  /* Computes the global stiffness matrix and rhs for any a(u,v) = <f,v> bilinear form using various element types
   * (eg. P1, P2, Nedelec, and Raviart-Thomas).
   * Also takes care of Dirichlet boundary conditions.  If the node is a boundary the row will be
   * zeroed out except for the diagonal entry being 1.  The corresponding column will also be 0 and 
   * the right-hand side adjusted.
   *
   * For this problem we compute:
   *
   * Lu = f  ---->   a(u,v) = <f,v>
   *
   * which gives Ax = b,
   *
   * A_ij = a( phi_j, phi_i)
   * b_i  = <f,phi_i>
   * 
   *    INPUT:
   *            local_assembly      Function to local assembly routine
   *            fespace		    Finite-Element Space Struct
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights
   *            rhs                 Function that gives right-hand side at given coordinates rhs(&val,x,time)
   *            bc                  Function that gives right-hand side at given coordinates bc(&val,x,time)
   *            coeff               Function that gives coefficient (for now assume constant)
   *            time                Time if time dependent
   *
   *    OUTPUT:
   *	       	A                   Global Stiffness Matrix (CSR Format)
   *            b	       	    Global Right hand side vector
   */

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT i,j;

  // Allocate Row Array
  A->row = FE->ndof;
  A->col = FE->ndof;
  A->IA = (INT *) calloc(FE->ndof+1,sizeof(INT));
  b->row = FE->ndof;
  b->val = (REAL *) calloc(b->row,sizeof(REAL));

  // Get Sparsity Structure First
  // Non-zeros of A and IA (ignores cancellations, so maybe more than necessary)
  stiffG_nnz(A,FE);
	
  // Columns of A -> JA
  A->JA = (INT *) calloc(A->nnz,sizeof(INT));
  stiffG_cols(A,FE);
  
  // Set values
  A->val = (REAL *) calloc(A->nnz,sizeof(REAL));
  for (i=0; i<A->nnz; i++) {
    A->val[i] = 0;
  }
  	
  // Now Build Global Matrix entries
  // First deal with boundary rows
  for (i=0; i<FE->ndof; i++) {
    if (FE->dof_bdry[i]==1) { /* This is a boundary row.  Just make identity and fix right hand side */
      A->val[A->IA[i]-1] = 1.0;
      b->val[i] = FE_Evaluate_DOF(bc,FE,mesh,time,i);
    }
  }

  // Now adjust other rows
  /* Loop over all Elements and build local matrix and rhs */
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* ALoc = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* bLoc = (REAL *) calloc(dof_per_elm,sizeof(REAL));
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT rowa,rowb,jcntr;
  for (i=0; i<FE->nelm; i++) {	
    // Zero out local matrices
    for (j=0; j<local_size; j++) {
      ALoc[j]=0;
    }
    for (j=0; j<dof_per_elm; j++) {
      bLoc[j]=0;
    }
		
    // Find DOF for given Element
    rowa = FE->el_dof->IA[i]-1;
    rowb = FE->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm[jcntr] = FE->el_dof->JA[j];
      jcntr++;
    }

    // Find vertices for given Element
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }
		
    // Compute Local Stiffness Matrix for given Element
    (*local_assembly)(ALoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,coeff,time);
    FEM_RHS_Local(bLoc,FE,mesh,cq,dof_on_elm,i,rhs,time);

    // Loop over DOF and place in appropriate slot globally
    LocaltoGlobal(dof_on_elm,FE,b,A,ALoc,bLoc);
  }	

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(ALoc) free(ALoc);
  if(bLoc) free(bLoc);
		
  return;
}
/******************************************************************************************************/

