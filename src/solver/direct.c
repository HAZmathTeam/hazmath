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
INT directsolve_UMF_symmetric(dCSRmat *A,dvector *f,REAL *x,INT print_level) 
{
	
  /* Performs Gaussian Elmination on Ax = f, using UMFPACK  
   * Assumes A is symmetric (I think) and in CSR format
   *
   * Input:		
   *        	A:	       	Matrix A to be solved
   *	       	f:	       	Right hand side vector
   *  print_level:              Flag to print messages to screen
   *
   * Output:	
   *	        x:	       	Solution
   *
   */
  
  INT i;
  INT shift_flag = 0;

  // UMFPACK requires transpose of A
  dCSRmat At;
  // Check if counting from 1 or 0 for CSR arrays
  if(A->IA[0]==1) {
    dcsr_shift(A, -1);  // shift A
    shift_flag = 1;
  }
  // Transpose A
  dcsr_trans(A,&At);

  // Fix A back to correct counting
  if(shift_flag==1) {
    dcsr_shift(A, 1);  // shift A back
  }

  // Reformat to longs
  long mydof2 = (long) At.row;
  long* iA = (long *) calloc(At.row+1,sizeof(long));
  long* jA = (long *) calloc(At.nnz,sizeof(long));
  for(i=0;i<At.row+1;i++) iA[i] = (long) At.IA[i];
  for(i=0;i<At.nnz;i++) jA[i] = (long) At.JA[i];

  // Get UMFPACK pointers
  void *Symbolic=NULL, *Numeric=NULL;
  double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO];

  if(print_level>=PRINT_SOME) {
    printf("\n=======================================================================");
    printf("\nUsing UMFPACK LU Factorization for Direct Solve of the System.\n");
    printf("=======================================================================\n\n");
  }

  /* symbolic analysis */
  umfpack_dl_symbolic(mydof2,mydof2,iA,jA,At.val,&Symbolic,Control,Info);
  if(print_level>=PRINT_SOME) {
    printf("\tSymbolic factorization done  -> ");
    if(Info[0]==0.0) {
      printf("No Errors\n");
    }
  }
  if(Info[0]!=0.0) {
    printf("Error in UMFPACK Symbolic Factorization: Info[0] = %f\n",Info[0]);
    return Info[0];
  }
  
  /* LU factorization */
  umfpack_dl_numeric(iA,jA,At.val,Symbolic,&Numeric,Control,Info);
  umfpack_dl_free_symbolic(&Symbolic);

  if(print_level>=PRINT_SOME) {
    printf("\tNumeric factorization done   -> ");
    if(Info[0]==0.0) {
      printf("No Errors\n");
    }
  }
  if(Info[0]!=0.0) {
    printf("Error in UMFPACK Numeric Factorization: Info[0] = %f\n",Info[0]);
    return Info[0];
  }

  /* solve system */
  umfpack_dl_solve(UMFPACK_A,iA,jA,At.val,x,f->val,Numeric,Control,Info);
  umfpack_dl_free_numeric(&Numeric);
  
  // Clean up
  dcsr_free(&At);
  if(iA) free(iA);
  if(jA) free(jA);

  if(print_level>=PRINT_SOME) {
    printf("\n=======================================================================");
    printf("\n UMFPACK Solve Complete.\n");
    printf("=======================================================================\n\n");
  }
 return 0;
}

/***************************************************************************************************************************/
