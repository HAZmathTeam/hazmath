/*! \file src/solver/direct.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 8/19/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 */

#include "hazmath.h"
#if WITH_SUITESPARSE
#include "umfpack.h"
#endif


/*! \brief Direct Solver methods -- for now just using UMFPACK Solvers
 *
 */
#if WITH_SUITESPARSE
/***************************************************************************************************************************/
/**
 * \fn directsolve_UMF(dCSRmat *A,dvector *f,REAL *x,INT print_level)
 *
 * \brief Performs Gaussian Elmination on Ax = f, using UMFPACK (Assumes A is in CSR format)
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
INT directsolve_UMF(dCSRmat *A,
                    dvector *f,
                    REAL *x,
                    INT print_level)
{ 
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
/**
 * \fn void* umfpack_factorize (dCSRmat *ptrA, const SHORT prtlvl)
 * \brief factorize A by UMFpack
 *
 * \param ptrA      Pointer to stiffness matrix of levelNum levels
 * \param Numeric   Pointer to the numerical factorization
 *
 * \author Xiaozhe Hu
 * \date   10/20/2010
 *
 * Modified by Xiaozhe Hu on 05/02/2014
 */
void* umfpack_factorize (dCSRmat *ptrA,
                         const SHORT prtlvl)
{   
    const INT n = ptrA->col;
    
    INT *Ap = ptrA->IA;
    INT *Ai = ptrA->JA;
    double *Ax = ptrA->val;
    void *Symbolic;
    void *Numeric;
    INT status = SUCCESS;
    
    clock_t start_time = clock();
    
    status = umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);
    if(status<0) {
	fprintf(stderr,"UMFPACK ERROR in Symbolic, status = %d\n\n",status);
    }
    status = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
    if(status<0) {
	fprintf(stderr,"UMFPACK ERROR in Numeric, status = %d\n\n",status);
    }
    umfpack_di_free_symbolic (&Symbolic);
        
    if ( prtlvl > PRINT_MIN ) {
        clock_t end_time = clock();
        double fac_time = (double)(end_time - start_time)/(double)(CLOCKS_PER_SEC);
        printf("UMFPACK factorize costs %f seconds.\n", fac_time);
    }
    
    return Numeric;
}

/***************************************************************************************************************************/
/**
 * \fn INT umfpack_solve (dCSRmat *ptrA, dvector *b, dvector *u,
 *                             void *Numeric, const SHORT prtlvl)
 * \brief Solve Au=b by UMFpack, numerical factorization is given
 *
 * \param ptrA      Pointer to stiffness matrix of levelNum levels
 * \param b         Pointer to the dvector of right hand side term
 * \param u         Pointer to the dvector of dofs
 * \param Numeric   Pointer to the numerical factorization
 * \param prtlvl    Output level
 *
 * \author Xiaozhe Hu
 * \date   10/20/2010
 *
 * Modified by Xiaozhe on 05/10/2014
 */
INT umfpack_solve (dCSRmat *ptrA,
                        dvector *b,
                        dvector *u,
                        void *Numeric,
                        const SHORT prtlvl)
{   

    INT *Ap = ptrA->IA;
    INT *Ai = ptrA->JA;
    double *Ax = ptrA->val;
    INT status = SUCCESS;
    
    clock_t start_time = clock();
    
    status = umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, u->val, b->val, Numeric, NULL, NULL);
    
    if ( prtlvl > PRINT_NONE ) {
        clock_t end_time = clock();
        double solve_time = (double)(end_time - start_time)/(double)(CLOCKS_PER_SEC);
        printf("UMFPACK costs %f seconds.\n", solve_time);
    }
    
    return status;
}

/***************************************************************************************************************************/
/**
 * \fn INT umfpack_free_numeric (void *Numeric)
 * \brief Solve Au=b by UMFpack
 *
 * \param Numeric   Pointer to the numerical factorization
 *
 * \author Xiaozhe Hu
 * \date   10/20/2010
 */
INT umfpack_free_numeric (void *Numeric)
{   
    INT status = SUCCESS;
    
    umfpack_di_free_numeric (&Numeric);
    
    return status;
}
#else
INT directsolve_UMF(dCSRmat *A, dvector *f,  REAL *x, INT print_level)
{
  fprintf(stderr,"\n\n*** FATAL ERROR: direct solver (umfpack) requested,\n");
  fprintf(stderr,"*** but suitesparse support is not compiled in the hazmath library\n");
  fprintf(stderr, "*** Either change the solver type ... OR ...\n\n"); 
  fprintf(stderr,"*** IF you have UMFPACK installed, THEN\n");
  fprintf(stderr,"*** try RECOMPILING the hazmath library with suitesparse support.\n\n");
  //  fprintf(stderr,"***by issuing the commands: make config suitesparse=yes; make install\n\n");
  exit(-1); 
}
#endif

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
