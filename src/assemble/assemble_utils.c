/*
 *  assemble_utils.c   
 *  
 *  Created by James Adler and Xiaozhe Hu on 4/22/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

/* This code will contain all the tools needed to build stiffness matrices */

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
void stiffG_nnz(dCSRmat *A, fespace *FE) 
{
  /* Computes the "possible" number of nonzeros for the global stiffness matrix
   * Also takes into account Dirichlet boundary conditions.  If the DOF of the corresponding row
   * is a boundary, only include 1 nonzero for that row (it will be identity!)
   *
   * Entry will be nonzero if degree of freedom belongs to the element (i.e. row)
   * 
   *
   *    INPUT:
   *          fespace		    Finite-Element Space Struct
   *
   *    OUTPUT:
   *          A->IA(ndof)	    Rows of Global A
   *	      A->nnzA		    Nonzeroes of A (ignoring cancellations)
   */
	
  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;
		
  // We will need the DOF to element map first
  iCSRmat dof_el;
  icsr_trans_1(FE->el_dof,&dof_el);

  INT ndof = FE->ndof;
	
  INT* ix = (INT *) calloc(ndof,sizeof(INT));
  for (i=0; i<ndof; i++) {
    ix[i] = 0;	
  }
	
  // Loop over all DOF and count possible nonzeros in A
  // Also build A->IA, while you're at it...
  icp=1;
  for (i=0; i<ndof; i++) {
    A->IA[i] = icp;
    // Check if this is a boundary edge
    if (FE->dof_bdry[i]!=0) {	/* Only 1 nonzero this row */
      icp++;
    } else {			/* It's not a boundary and compute as normal */
      // Loop over all Elements connected to particular edge
      j_a = dof_el.IA[i];
      j_b = dof_el.IA[i+1]-1;
      for (j=j_a; j<=j_b; j++) {
	if1 = dof_el.JA[j-1];
	k_a = FE->el_dof->IA[if1-1];
	k_b = FE->el_dof->IA[if1]-1;
	for (k=k_a; k<=k_b; k++) {
	  mydof = FE->el_dof->JA[k-1];
	  if (ix[mydof-1]!=i+1 && FE->dof_bdry[mydof-1]==0) { /* We haven't been here AND it's not a boundary */
	    icp++;
	    ix[mydof-1] = i+1;
	  }
	}
				
      }
    }
  }
  A->IA[ndof] = icp;
  A->nnz = icp-1;
		
  if(ix) free(ix);
  icsr_free(&dof_el);
	
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void stiffG_cols(dCSRmat *A, fespace *FE) 
{
  /* Finds the column sparsity structure of the  Global Stiffness Matrix
   * Also takes into account Dirichlet boundary conditions.  If the DOF of the corresponding row
   * is a boundary, only include 1 nonzero for that row (it will be identity!)
   *
   * Entry will be nonzero if degree of freedom belongs to the element (i.e. row)
   * 
   *
   *    INPUT:
   *          fespace		    Finite-Element Space Struct
   *
   *    OUTPUT:
   *          A->IA(ndof)	    Rows of Global A
   *	      A->nnzA		    Nonzeroes of A (ignoring cancellations)
   */
	
  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;
	
  // We will need the DOF to element map first
  iCSRmat dof_el;
  icsr_trans_1(FE->el_dof,&dof_el);

  INT ndof = FE->ndof;
	
  INT* ix = (INT *) calloc(ndof,sizeof(INT));
  for (i=0; i<ndof; i++) {
    ix[i] = 0;	
  }
	
  // Loop over all DOF and build A->JA
  icp=1;
  for (i=0; i<ndof; i++) {
    // Check if this is a boundary edge
    if (FE->dof_bdry[i]==1) {	/* Only 1 nonzero this row */
      A->JA[icp-1]=i+1;
      icp++;
    } else {			/* It's not a boundary and compute as normal */
      // Loop over all Elements connected to particular edge
      j_a = dof_el.IA[i];
      j_b = dof_el.IA[i+1]-1;
      for (j=j_a; j<=j_b; j++) {
	if1 = dof_el.JA[j-1];
	k_a = FE->el_dof->IA[if1-1];
	k_b = FE->el_dof->IA[if1]-1;
	for (k=k_a; k<=k_b; k++) {
	  mydof = FE->el_dof->JA[k-1];
	  if (ix[mydof-1]!=i+1 && FE->dof_bdry[mydof-1]==0) { /* We haven't been here AND it's not a boundary */
	    A->JA[icp-1] = mydof;
	    icp++;
	    ix[mydof-1] = i+1;
	  }
	}
				
      }
    }
  }
	
	
  if(ix) free(ix);
  icsr_free(&dof_el);
	
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void LocaltoGlobal(INT *dof_on_elm,fespace *FE,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc) 
{
  /********* Maps the Local Matrix to Global Matrix considering boundaries *********************
   *
   *	Input:
   *		dof_on_elm	Array of DOF entries on the given element
   *		FE              finite-element space struct
   *		b               Global right hand side
   *		A		Global Stiffness Matrix
   *		ALoc		Local Stiffness Matrix
   *		bLoc		Local right hand side
   *
   *	Output:		
   *            A		Adjusted Global Matrix
   *		b		Adjusted Global right-hand side
   *
   */
	
  INT i,j,k,row,col,col_a,col_b,acol;
	
  for (i=0; i<FE->dof_per_elm; i++) { /* Rows of Local Stiffness */
    row = dof_on_elm[i]-1;
    if (FE->dof_bdry[row]==0) { /* Only if not on a boundary */		
      // Adjust Right-hand side globally
      b->val[row] = b->val[row] + bLoc[i];
			
      for (j=0; j<FE->dof_per_elm; j++) { /* Columns of Local Stiffness */
	col = dof_on_elm[j]-1;
	if (FE->dof_bdry[col]==0) { /* Only do stuff if hit a non-boundary edge */
	  col_a = A->IA[row]-1;
	  col_b = A->IA[row+1]-1;
	  for (k=col_a; k<col_b; k++) { /* Columns of A */
	    acol = A->JA[k]-1;				
	    if (acol==col) {	/* If they match, put it in the global matrix */
	      A->val[k] = A->val[k] + ALoc[i*FE->dof_per_elm+j];
	    }
	  }
	} else { /* If boundary adjust Right hand side */
	  b->val[row] = b->val[row]-(ALoc[i*FE->dof_per_elm+j])*b->val[col];
	}
      }
    }
  }
  return;
}
/******************************************************************************************************/
