/*! \file src/solver/direct.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 8/19/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 * \brief Direct Solver methods -- for now just using UMFPACK Solvers
 *
 * \note   Done cleanup for releasing -- Xiaozhe Hu 08/28/2021
 *
 */

#include "hazmath.h"

#if WITH_SUITESPARSE
#include "umfpack.h"
#endif

/***********************************************************************************************/
/**
 * \fn directsolve_UMF(dCSRmat *A,dvector *f,dvector *x,INT print_level)
 *
 * \brief Performs Gaussian Elmination on Ax = f, using UMFPACK (Assumes A is in CSR format)
 * \note This routine does everything, factorization and solve in one shot.  So it is not
 *       efficient if you want to solve the same matrix system many times, since it will not
 *       save the factorization.
 *
 * \param A	       	        Matrix A to be solved
 * \param f	       	        Right hand side vector
 * \param print_level       Flag to print messages to screen
 *
 * \return x	         	Solution
 *
 */
INT directsolve_UMF(dCSRmat *A,
                    dvector *f,
                    dvector *x,
                    INT print_level)
{

#if WITH_SUITESPARSE
  INT err_flag;
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

  // Data for Numerical factorization
  void* Numeric;

  // Factorize
  Numeric = umfpack_factorize(&At,print_level);

  // Solve
  err_flag = umfpack_solve(&At,f,x,Numeric,print_level);

  // Error Check
  if(err_flag<0) {
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- UMFPACK SOLVE ERROR!!!\n\n",__FUNCTION__);
    exit(err_flag);
  }

  // Clean up
  dcsr_free(&At);
  err_flag = umfpack_free_numeric(Numeric);
  if(err_flag<0) {
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- UMFPACK FREE ERROR!!!\n\n",__FUNCTION__);
    exit(err_flag);
  }

  return err_flag;
#else
  error_extlib(252, __FUNCTION__, "SuiteSparse");
  return 0;
#endif
}

/***********************************************************************************************/
/**
 * \fn void* factorize_UMF(dCSRmat *A,INT print_level)
 *
 * \brief Performs factorization of A using UMFPACK (Assumes A is in CSR format)
 * \note This routine does factorization only.
 *
 * \param A	       	        Matrix A to be solved (no need of transpose)
 * \param print_level       Flag to print messages to screen
 *
 * \return Numeric	        Stores LU decomposition of A
 *
 * \author    Xiaozhe Hu
 *
 */
void* factorize_UMF(dCSRmat *A,
                    INT print_level)
{

#if WITH_SUITESPARSE
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

  // Data for Numerical factorization
  void* Numeric;

  // Factorize
  Numeric = umfpack_factorize(&At,print_level);

  // Clean up
  dcsr_free(&At);

  return Numeric;
#else
  error_extlib(252, __FUNCTION__, "SuiteSparse");
  return 0;
#endif
}

/***********************************************************************************************/
/**
 * \fn INT solve_UMF(dCSRmat *A,dvector *f,dvector *x, void* Numeric, INT print_level)
 *
 * \brief Performs Gaussian Elmination on Ax = f, using UMFPACK
 *        (Assumes A is in CSR format and factorization has been done)
 *
 * \note This routine does solve only
 *
 * \param A	       	        Matrix A to be solved
 * \param f	       	        Right hand side vector
 * \param Numeric           LU decomposition of A
 * \param print_level       Flag to print messages to screen
 *
 * \return x	         	    Solution
 *
 * \author Xiaozhe Hu
 *
 */
INT solve_UMF(dCSRmat *A,
              dvector *f,
              dvector *x,
              void *Numeric,
              INT print_level)
{

#if WITH_SUITESPARSE
  INT err_flag;
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

  // Solve
  err_flag = umfpack_solve(&At,f,x,Numeric,print_level);

  // Error Check
  if(err_flag<0) {
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- UMFPACK SOLVE ERROR!!!\n\n",__FUNCTION__);
    exit(err_flag);
  }

  // Clean up
  dcsr_free(&At);

  return err_flag;
#else
  error_extlib(252, __FUNCTION__, "SuiteSparse");
  return 0;
#endif
}

/***********************************************************************************************/
/**
 * \fn block_directsolve_UMF(block_dCSRmat *bA,dvector *f,dvector *x,INT print_level)
 *
 * \brief Performs Gaussian Elmination on Ax = f, using UMFPACK (Assumes A is in block_dCSR format)
 * \note This routine does everything, factorization and solve in one shot.  So it is not
 *       efficient if you want to solve the same matrix system many times, since it will not
 *       save the factorization.
 *
 * \param bA	       	      Matrix A to be solved
 * \param f	       	        Right hand side vector
 * \param print_level       Flag to print messages to screen
 *
 * \return x	         	Solution
 *
 */
INT block_directsolve_UMF(block_dCSRmat *bA,
                    dvector *f,
                    dvector *x,
                    INT print_level)
{
#if WITH_SUITESPARSE
  INT err_flag;

  // Convert block matrix to regular matrix
  dCSRmat A = bdcsr_2_dcsr(bA);

  dcsr_compress_inplace(&A, 1.0e-14);

  // Call regular solve
  err_flag = directsolve_UMF(&A,f,x,print_level);

  // Clean up
  dcsr_free(&A);

  return err_flag;
#else
  error_extlib(252, __FUNCTION__, "SuiteSparse");
  return 0;
#endif
}

/***********************************************************************************************/
/**
 * \fn void* block_factorize_UMF(block_dCSRmat *bA,INT print_level)
 *
 * \brief Performs factorization of A using UMFPACK (Assumes A is in block_dCSR format)
 * \note This routine does factorization only.
 *
 * \param bA	       	        Matrix A to be solved (no need of transpose)
 *
 * \return Numeric	        Stores LU decomposition of A
 *
 */
void* block_factorize_UMF(block_dCSRmat *bA,
                    INT print_level)
{

#if WITH_SUITESPARSE

  // Data for Numerical factorization
  void* Numeric = NULL;

  // Convert block matrix to regular matrix
  dCSRmat A = bdcsr_2_dcsr(bA);

  // Factorize
  Numeric = factorize_UMF(&A,print_level);

  // Cleanup
  dcsr_free(&A);

  return Numeric;
#else
  error_extlib(252, __FUNCTION__, "SuiteSparse");
  return 0;
#endif
}

/***********************************************************************************************/
/**
 * \fn INT block_solve_UMF(block_dCSRmat *bA,dvector *f,dvector *x, void* Numeric, INT print_level)
 *
 * \brief Performs Gaussian Elmination on Ax = f, using UMFPACK
 *        (Assumes A is in block_dCSR format and factorization has been done)
 *
 * \note This routine does solve only
 *
 * \param A	       	        Matrix A to be solved
 * \param f	       	        Right hand side vector
 * \param Numeric           LU decomposition of A
 * \param print_level       Flag to print messages to screen
 *
 * \return x	         	Solution
 *
 */
INT block_solve_UMF(block_dCSRmat *bA,
              dvector *f,
              dvector *x,
              void *Numeric,
              INT print_level)
{

#if WITH_SUITESPARSE
  INT err_flag;

  // Convert block matrix to regular matrix
  dCSRmat A = bdcsr_2_dcsr(bA);

  // Call regular solve
  err_flag = solve_UMF(&A,f,x,Numeric,print_level);

  // Clean up
  dcsr_free(&A);

  return err_flag;
#else
  error_extlib(252, __FUNCTION__, "SuiteSparse");
  return 0;
#endif
}

// Suitesparse routines - assumes conversions have been done
/***********************************************************************************************/
/**
 * \fn void* umfpack_factorize (dCSRmat *ptrA, const SHORT prtlvl)
 * \brief factorize A by UMFpack
 *
 * \param ptrA      Pointer to dCSRmat matrix (transpose has been done!)
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

#if WITH_SUITESPARSE
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
#else
  error_extlib(253, __FUNCTION__, "SuiteSparse");
  return NULL;
#endif

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

#if WITH_SUITESPARSE
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
#else
  error_extlib(254, __FUNCTION__, "SuiteSparse");
  return 0;
#endif
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
#if WITH_SUITESPARSE
  INT status = SUCCESS;

  umfpack_di_free_numeric (&Numeric);

  return status;
#else
  error_extlib(255, __FUNCTION__, "SuiteSparse");
  return 0;
#endif
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
