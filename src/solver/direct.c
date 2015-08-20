/*
 *  direct.c
 *
 *  Created by James Adler and Xiaozhe Hu on 8/19/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

#include "hazmat.h"
#include "umfpack.h"


/*! \brief Direct Solver methods -- for now just using UMFPACK Solvers
 *
 */

/***************************************************************************************************************************/
INT directsolve_UMF_symmetric(dCSRmat *A,dvector *f,REAL *x) 
{
	
  /* Performs Gaussian Elmination on Ax = f, using UMFPACK  
   * Assumes A is symmetric (I think) and in CSR format
   * Also assumes counting of arrays starts at 1.
   *
   * Input:		
   *        	A:	       	Matrix A to be solved
   *	       	f:	       	Right hand side vector
   *
   * Output:	
   *	        x:	       	Solution
   *
   */

  // UMFPACK pointers
  void *Symbolic=NULL;
  void *Numeric=NULL;

  // Size of A
  INT ndof = f->row;
  if(ndof!=A->row) {
    printf("RHS and Matrix don't match in Direct Solver");
    exit(0);
  }

  /* Symbolic analysis */
  umfpack_di_symbolic(ndof,ndof,A->IA,A->JA,A->val,Symbolic,NULL,NULL);
  printf("SYMBOLIC =%d\n",Symbolic);

  /* LU factorization */
  umfpack_di_numeric(A->IA,A->JA,A->val,Symbolic,Numeric,NULL,NULL);
  umfpack_di_free_symbolic(Symbolic);

  /* solve system */
  umfpack_di_solve(UMFPACK_A,A->IA,A->JA,A->val,x,f->val,Numeric,NULL,NULL);
  printf("SYMBOLIC =%d\n",Numeric);
  
  umfpack_di_free_numeric(Numeric);

  // Need to add the error checks
  return 1;
}
/***************************************************************************************************************************/
