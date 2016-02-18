/*
 *  assemble_utils.c   
 *  
 *  Created by James Adler and Xiaozhe Hu on 4/22/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *  Modified 2015-08-08 (ltz)
 */

/* This code will contain all the tools needed to build stiffness matrices */

#include "hazmat.h"

/******************************************************************************************************/
void stiffG_nnz(dCSRmat *A, fespace *FE) 
{
  /* Computes the "possible" number of nonzeros for the global stiffness matrix
   * Ignores Dirichlet boundary conditions, which can be eliminated later.  
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
    // Loop over all Elements connected to particular DOF
    j_a = dof_el.IA[i];
    j_b = dof_el.IA[i+1]-1;
    for (j=j_a; j<=j_b; j++) {
      if1 = dof_el.JA[j-1];
      k_a = FE->el_dof->IA[if1-1];
      k_b = FE->el_dof->IA[if1]-1;
      for (k=k_a; k<=k_b; k++) {
	mydof = FE->el_dof->JA[k-1];
	if (ix[mydof-1]!=i+1) { /* We haven't been here  */
	  icp++;
	  ix[mydof-1] = i+1;
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
   * Ignores Dirichlet boundary conditions to be eliminated later.  
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
    // Loop over all Elements connected to particular edge
    j_a = dof_el.IA[i];
    j_b = dof_el.IA[i+1]-1;
    for (j=j_a; j<=j_b; j++) {
      if1 = dof_el.JA[j-1];
      k_a = FE->el_dof->IA[if1-1];
      k_b = FE->el_dof->IA[if1]-1;
      for (k=k_a; k<=k_b; k++) {
	mydof = FE->el_dof->JA[k-1];
	if (ix[mydof-1]!=i+1) { /* We haven't been here */
	  A->JA[icp-1] = mydof;
	  icp++;
	  ix[mydof-1] = i+1;
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
  /********* Maps the Local Matrix to Global Matrix NOT considering boundaries *********************
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
    // Adjust Right-hand side globally
    if(bLoc!=NULL)
      b->val[row] = b->val[row] + bLoc[i];
			
    for (j=0; j<FE->dof_per_elm; j++) { /* Columns of Local Stiffness */
      col = dof_on_elm[j]-1;
    
      col_a = A->IA[row]-1;
      col_b = A->IA[row+1]-1;
      for (k=col_a; k<col_b; k++) { /* Columns of A */
	acol = A->JA[k]-1;				
	if (acol==col) {	/* If they match, put it in the global matrix */
	  A->val[k] = A->val[k] + ALoc[i*FE->dof_per_elm+j];
	}
      }
     
    }
  }
  return;
}
/******************************************************************************************************/

// The following are similar routines above, but eliminates Dirichlet boundary conditions during the
// construction of the matrices.

/******************************************************************************************************/
void stiffG_nnzBC(dCSRmat *A, fespace *FE) 
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
void stiffG_colsBC(dCSRmat *A, fespace *FE) 
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
void LocaltoGlobalBC(INT *dof_on_elm,fespace *FE,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc) 
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
      if(bLoc!=NULL)
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
	  if(bLoc!=NULL)
	    b->val[row] = b->val[row]-(ALoc[i*FE->dof_per_elm+j])*b->val[col];
	}
      }
    }
  }
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void stiffG_nnz_subset(dCSRmat *A, fespace *FE,INT flag) 
{
  /* Computes the "possible" number of nonzeros for the global stiffness matrix
   * Also takes into account special boundary conditions.  Here if the boundary is equal to "flag" (whatever that corresponds to)
   * Then we add something to the matrix
   * 
   *
   *    INPUT:
   *          fespace		    Finite-Element Space Struct
   *             flag               Indicates boundary term we want
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
    // Check if this is a boundary dof we want
    if (FE->dof_bdry[i]==flag) {
      // Loop over all Elements connected to particular edge
      j_a = dof_el.IA[i];
      j_b = dof_el.IA[i+1]-1;
      for (j=j_a; j<=j_b; j++) {
	if1 = dof_el.JA[j-1];
	k_a = FE->el_dof->IA[if1-1];
	k_b = FE->el_dof->IA[if1]-1;
	for (k=k_a; k<=k_b; k++) {
	  mydof = FE->el_dof->JA[k-1];
	  if (ix[mydof-1]!=i+1 && FE->dof_bdry[mydof-1]==flag) { /* We haven't been here AND it is a boundary */
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
void stiffG_cols_subset(dCSRmat *A, fespace *FE,INT flag) 
{
  /* Finds the column sparsity structure of the  Global Stiffness Matrix
   * Also takes into account special boundary conditions.  Here if the boundary is equal to "flag" (whatever that corresponds to)
   * Then we add something to the matrix
   *
   *    INPUT:
   *          fespace		    Finite-Element Space Struct
   *             flag               Indicates boundary term we want
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
    if (FE->dof_bdry[i]==flag) {
      // Loop over all Elements connected to particular edge
      j_a = dof_el.IA[i];
      j_b = dof_el.IA[i+1]-1;
      for (j=j_a; j<=j_b; j++) {
	if1 = dof_el.JA[j-1];
	k_a = FE->el_dof->IA[if1-1];
	k_b = FE->el_dof->IA[if1]-1;
	for (k=k_a; k<=k_b; k++) {
	  mydof = FE->el_dof->JA[k-1];
	  if (ix[mydof-1]!=i+1 && FE->dof_bdry[mydof-1]==flag) { /* We haven't been here AND it is a boundary */
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
void LocaltoGlobal_face(INT *dof_on_f,INT dof_per_f,fespace *FE,dCSRmat *A,REAL *ALoc,INT flag) 
{
  /********* Maps the Local Matrix to Global Matrix considering "special" boundaries *****************
   *
   *	Input:
   *		dof_on_elm	Array of DOF entries on the given element
   *		FE              finite-element space struct
   *		A		Global Stiffness Matrix
   *		ALoc		Local Stiffness Matrix
   *		flag		Indicates which boundary DOF to grab
   *
   *	Output:		
   *            A		Adjusted Global Matrix
   *
   */
	
  INT i,j,k,row,col,col_a,col_b,acol;
	
  for (i=0; i<dof_per_f; i++) { /* Rows of Local Stiffness */
    row = dof_on_f[i]-1;
    if (FE->dof_bdry[row]==flag) { /* Only if on special boundary */		
			
      for (j=0; j<dof_per_f; j++) { /* Columns of Local Stiffness */
	col = dof_on_f[j]-1;
	if (FE->dof_bdry[col]==flag) { /* Only do stuff if hit a special boundary */
	  col_a = A->IA[row]-1;
	  col_b = A->IA[row+1]-1;
	  for (k=col_a; k<col_b; k++) { /* Columns of A */
	    acol = A->JA[k]-1;				
	    if (acol==col) {	/* If they match, put it in the global matrix */
	      A->val[k] = A->val[k] + ALoc[i*dof_per_f+j];
	    }
	  }
	}
      }
    }
  }
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void eliminate_DirichletBC_old(void (*bc)(REAL *,REAL *,REAL),fespace *FE,trimesh *mesh,dvector *b,dCSRmat *A,REAL time) 
{
  /********* Eliminates the Dirichlet boundaries from the global matrix *********************
   * 
   *  for each row in A that corresponds to a Dirichlet boundary, make the diagonal 1 and 
   *  off-diagonals zero.  Then, making the corresponding column entries 0.
   *  For the rhs, if it's a boundary DOF, set to the actual boundary value
   *  If it's DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
   *
   *	Input:
   *            bc              Function that gives boundary condition at given coordinates bc(&val,x,time)
   *		FE              finite-element space struct
   *	      	mesh            Mesh Struct
   *		b               Global right hand side
   *		A		Global Stiffness Matrix
   *            time            Time if time dependent
   *
   *	Output:		
   *            A		Adjusted Global Matrix
   *		b		Adjusted Global right-hand side
   *
   */

  INT i,j,cola,colb;
  INT ndof = FE->ndof;

  // Error Check
  if(A->row!=ndof || A->col!=ndof) {
    printf("\nYou're matrix size doesn't match your number of DOF in Boundary Elimination\n");
    exit(0);
  }

  
  // Loop over rows of A
  for (i=0; i<ndof; i++) {
    cola = A->IA[i]-1;
    colb = A->IA[i+1]-1;
    if(FE->dof_bdry[i]==1) { // Boundary Row
      // Set rhs to BC
      if(b!=NULL)
	b->val[i] = FE_Evaluate_DOF(bc,FE,mesh,time,i); 
      // Loop over columns and 0 out row except for diagonal
      for(j=cola; j<colb; j++) { 
	if(A->JA[j]==i+1)
	  A->val[j] = 1.0;
	else
	  A->val[j] = 0.0;
      }
    } else { // Non-boundary-row
      // Loop over columns and 0 out if column is boundary
      for(j=cola; j<colb; j++) { 
	if(FE->dof_bdry[A->JA[j]-1]==1) {
	  // Adjust RHS accordingly as well
	  if(b!=NULL)
	    b->val[i] = b->val[i] - A->val[j]*FE_Evaluate_DOF(bc,FE,mesh,time,A->JA[j]-1);
	  // Zero out matrix entry
	  A->val[j] = 0.0;
	}
      }
    }
  }
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void eliminate_DirichletBC(void (*bc)(REAL *,REAL *,REAL),fespace *FE,trimesh *mesh,dvector *b,dCSRmat *A,REAL time) 
{
  /********* Eliminates the Dirichlet boundaries from the global matrix *********************
   * 
   *  for each row in A that corresponds to a Dirichlet boundary, make the diagonal 1 and 
   *  off-diagonals zero.  Then, making the corresponding column entries 0.
   *  For the rhs, if it's a boundary DOF, set to the actual boundary value
   *  If it's DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
   *
   *	Input:
   *            bc              Function that gives boundary condition at given coordinates bc(&val,x,time)
   *		FE              finite-element space struct
   *	      	mesh            Mesh Struct
   *		b               Global right hand side
   *		A		Global Stiffness Matrix
   *            time            Time if time dependent
   *
   *	Output:		
   *            A		Adjusted Global Matrix
   *		b		Adjusted Global right-hand side
   *
   */

  INT i,j,cola,colb;
  INT ndof = FE->ndof;

  // Error Check
  if(A->row!=ndof || A->col!=ndof) {
    printf("\nYour matrix size doesn't match your number of DOF in Boundary Elimination\n");
    exit(0);
  }

  // Fix RHS First b_interior = (b - A(0;u_bdry)^T)_interior
  //               b_bdry = u_bdry
  if(b!=NULL) 
    eliminate_DirichletBC_RHS(bc,FE,mesh,b,A,time);

  // Now fix A 
  // Loop over rows of A
  for (i=0; i<ndof; i++) {
    cola = A->IA[i]-1;
    colb = A->IA[i+1]-1;
    if(FE->dof_bdry[i]==1) { // Boundary Row 
      // Loop over columns and 0 out row except for diagonal
      for(j=cola; j<colb; j++) { 
	if(A->JA[j]==i+1)
	  A->val[j] = 1.0;
	else
	  A->val[j] = 0.0;
      }
    } else { // Non-boundary-row
      // Loop over columns and 0 out if column is boundary
      for(j=cola; j<colb; j++) { 
	if(FE->dof_bdry[A->JA[j]-1]==1) {
	  // Zero out matrix entry
	  A->val[j] = 0.0;
	}
      }
    }
  }
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void eliminate_DirichletBC_RHS(void (*bc)(REAL *,REAL *,REAL),fespace *FE,trimesh *mesh,dvector *b,dCSRmat *A,REAL time) 
{
  /********* Eliminates the Dirichlet boundaries from the global RHS vector *********************
   * 
   *  For the rhs, if it's a boundary DOF, set to the actual boundary value
   *  If it's DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
   *
   *	Input:
   *            bc              Function that gives boundary condition at given coordinates bc(&val,x,time)
   *		FE              finite-element space struct
   *	      	mesh            Mesh Struct
   *		b               Global right hand side
   *		A		Global Stiffness Matrix (unmodified)
   *            time            Time if time dependent
   *
   *	Output:		
   *		b		Adjusted Global right-hand side
   *
   */

  INT i,j,cola,colb;
  INT ndof = FE->ndof;
  REAL* ub = (REAL *) calloc(ndof,sizeof(REAL));

  // Error Check
  if(A->row!=ndof || A->col!=ndof) {
    printf("\nYou're matrix size doesn't match your number of DOF in Boundary Elimination\n");
    exit(0);
  }

  // Get solution vector that's 0 on interior and boundary value on boundary
  for(i=0; i<ndof; i++) {
    if(FE->dof_bdry[i]==1) {
      ub[i] = FE_Evaluate_DOF(bc,FE,mesh,time,i);
    } else {
      ub[i] = 0.0;
    }
  }

  // b = b - Aub
  dcsr_aAxpy_1(-1.0,A,ub,b->val);

  // Fix boundary values
  for(i=0;i<ndof;i++) {
    if(FE->dof_bdry[i]==1) {
      b->val[i] = ub[i];
    }
  }

  if(ub) free(ub);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void eliminate_DirichletBC_blockFE(void (*bc)(REAL *,REAL *,REAL),block_fespace *FE,trimesh *mesh,dvector *b,dCSRmat *A,REAL time) 
{
  /********* Eliminates the Dirichlet boundaries from the global matrix *********************
   *     Assumes block structure of FE space
   *  for each row in A that corresponds to a Dirichlet boundary, make the diagonal 1 and 
   *  off-diagonals zero.  Then, making the corresponding column entries 0.
   *  For the rhs, if it's a boundary DOF, set to the actual boundary value
   *  If it's DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
   *
   *	Input:
   *            bc              Function that gives boundary condition at given coordinates bc(&val,x,time)
   *		FE              finite-element space struct (block format)
   *	      	mesh            Mesh Struct
   *		b               Global right hand side
   *		A		Global Stiffness Matrix
   *            time            Time if time dependent
   *
   *	Output:		
   *            A		Adjusted Global Matrix
   *		b		Adjusted Global right-hand side
   *
   */

  INT i,j,cola,colb;
  INT ndof = 0;
  for(i=0;i<FE->nspaces;i++) ndof+=FE->var_spaces[i]->ndof;

  // Error Check
  if(A->row!=ndof || A->col!=ndof) {
    printf("\nYour matrix size doesn't match your number of DOF in Boundary Elimination\n");
    exit(0);
  }

  // Fix RHS First b_interior = (b - A(0;u_bdry)^T)_interior
  //               b_bdry = u_bdry
  if(b!=NULL) 
    eliminate_DirichletBC_RHS_blockFE(bc,FE,mesh,b,A,time);

  // Now fix A 
  // Loop over rows of A
  for (i=0; i<ndof; i++) {
    cola = A->IA[i]-1;
    colb = A->IA[i+1]-1;
    if(FE->dof_bdry[i]==1) { // Boundary Row 
      // Loop over columns and 0 out row except for diagonal
      for(j=cola; j<colb; j++) { 
	if(A->JA[j]==i+1)
	  A->val[j] = 1.0;
	else
	  A->val[j] = 0.0;
      }
    } else { // Non-boundary-row
      // Loop over columns and 0 out if column is boundary
      for(j=cola; j<colb; j++) { 
	if(FE->dof_bdry[A->JA[j]-1]==1) {
	  // Zero out matrix entry
	  A->val[j] = 0.0;
	}
      }
    }
  }
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void eliminate_DirichletBC_RHS_blockFE(void (*bc)(REAL *,REAL *,REAL),block_fespace *FE,trimesh *mesh,dvector *b,dCSRmat *A,REAL time) 
{
  /********* Eliminates the Dirichlet boundaries from the global RHS vector *********************
   *    Assumes block structure of FE space
   *  For the rhs, if it's a boundary DOF, set to the actual boundary value
   *  If it's DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
   *
   *	Input:
   *            bc              Function that gives boundary condition at given coordinates bc(&val,x,time)
   *		FE              finite-element space struct (block structure)
   *	      	mesh            Mesh Struct
   *		b               Global right hand side
   *		A		Global Stiffness Matrix (unmodified)
   *            time            Time if time dependent
   *
   *	Output:		
   *		b		Adjusted Global right-hand side
   *
   */

  INT i,j,cola,colb;
  INT ndof = FE->ndof;
  REAL* ub = (REAL *) calloc(ndof,sizeof(REAL));

  // Error Check
  if(A->row!=ndof || A->col!=ndof) {
    printf("\nYou're matrix size doesn't match your number of DOF in Boundary Elimination\n");
    exit(0);
  }

  // Get solution vector that's 0 on interior and boundary value on boundary
  for(i=0; i<ndof; i++) {
    if(FE->dof_bdry[i]==1) {
      ub[i] = blockFE_Evaluate_DOF(bc,FE,mesh,time,i);
    } else {
      ub[i] = 0.0;
    }
  }

  // b = b - Aub
  dcsr_aAxpy_1(-1.0,A,ub,b->val);

  // Fix boundary values
  for(i=0;i<ndof;i++) {
    if(FE->dof_bdry[i]==1) {
      b->val[i] = ub[i];
    }
  }

  if(ub) free(ub);

  return;
}
/******************************************************************************************************/
