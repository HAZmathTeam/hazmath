/*! \file src/solver/direct.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 8/19/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 * \brief Direct Solver methods
 *
 * \note   Done cleanup for releasing -- Xiaozhe Hu 08/28/2021
 * \note   modified -- Ludmil 20220802
 *
 */

#include "hazmath.h"

#if WITH_SUITESPARSE
#include "umfpack.h"
#endif

/************************************************************/
/**
 * \fn directsolve_HAZ(dCSRmat *A,dvector *f,dvector *x,INT print_level)
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
 * \note modified by Ludmil on 20220802 to include hazmath direct
 *       solve option.
 *
 */
INT directsolve_HAZ(dCSRmat *A,
                    dvector *f,
                    dvector *x,
                    INT print_level)
{
  // Data for Numerical factorization
  void* Numeric;
  INT err_flag_s,err_flag_f;

#if WITH_SUITESPARSE
  INT shift_flag = 0;
  // Check if counting from 1 or 0 for CSR arrays
  if(A->IA[0]==1) {
    dcsr_shift(A, -1);  // shift A
    shift_flag = 1;
  }
  // UMFPACK requires transpose of A
  dCSRmat At;
  // Transpose A
  dcsr_trans(A,&At);
  
  // Fix A back to correct counting
  if(shift_flag==1) {
    dcsr_shift(A, 1);  // shift A back
  }
  // Factorize
  //  fprintf(stdout,"\n%s: USING UMFPACK\n\n",__FUNCTION__);fflush(stdout);
  Numeric = hazmath_factorize(&At,print_level);
  // Solve
  err_flag_s = hazmath_solve(&At,f,x,Numeric,print_level);
  // Clean up
  dcsr_free(&At);
  err_flag_f = hazmath_free_numeric(&Numeric);
#else
  // HAZ Factorize
  //  fprintf(stdout,"\n%s: USING HAZMATH\n\n",__FUNCTION__);fflush(stdout);
  Numeric = run_hazmath_factorize(A,print_level);
  // HAZ Solve
  err_flag_s = run_hazmath_solve(A,f,x,Numeric,print_level);
  // clean up.
  err_flag_f = hazmath_free_numeric(&Numeric);
  //  error_extlib(252, __FUNCTION__, "SuiteSparse not present");
  //  return 0;
#endif
  // Error Check
  if(err_flag_s<0) {
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- DIRECT SOLVE ERROR!!!\n\n",__FUNCTION__);
    exit(err_flag_s);
  }
  if(err_flag_f<0) {
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- DIRECT SOLVE FREE ERROR!!!\n\n",__FUNCTION__);
    exit(err_flag_f);
  }
  return err_flag_s;
}

/***********************************************************************************************/
/**
 * \fn void* factorize_HAZ(dCSRmat *A,INT print_level)
 *
 * \brief Performs factorization of A using UMFPACK (Assumes A is in CSR format)
 * \note This routine does factorization only.
 *
 * \param A	       	        Matrix A to be solved (no need of transpose)
 * \param print_level       Flag to print messages to screen
 *
 * \return Numeric	        Stores LU decomposition of A
 *
 * \note modified by Ludmil on 20220802 to include hazmath direct
 *       solve option.
 *
 * \author    Xiaozhe Hu
 *
 * \note modified by Ludmil on 20220802 to include hazmath direct
 *       solve option.
 *
 */
void* factorize_HAZ(dCSRmat *A,
                    INT print_level)
{
  // Data for Numerical factorization
  void* Numeric;

#if WITH_SUITESPARSE

  INT shift_flag = 0;  
  // Check if counting from 1 or 0 for CSR arrays
  if(A->IA[0]==1) {
    dcsr_shift(A, -1);  // shift A
    shift_flag = 1;
  }

  // UMFPACK requires transpose of A
  dCSRmat At;
  // Transpose A
  dcsr_trans(A,&At);

  // Fix A back to correct counting
  if(shift_flag==1) {
    dcsr_shift(A, 1);  // shift A back
  }

  // Factorize
  Numeric = hazmath_factorize(&At,print_level);
  // Clean up
  dcsr_free(&At);
#else
  // Factorize
  Numeric = hazmath_factorize(A,print_level);
  //  error_extlib(252, __FUNCTION__, "SuiteSparse");
  //  return 0;
#endif
  return Numeric;
}

/***********************************************************************************************/
/**
 * \fn INT solve_HAZ(dCSRmat *A,dvector *f,dvector *x, void* Numeric, INT print_level)
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
 * \note modified by Ludmil on 20220802 to include hazmath direct
 *       solve option.
 *
 */
INT solve_HAZ(dCSRmat *A,
              dvector *f,
              dvector *x,
              void *Numeric,
              INT print_level)
{

  INT err_flag;
  // UMFPACK requires transpose of A
#if WITH_SUITESPARSE
  dCSRmat At;
  // Check if counting from 1 or 0 for CSR arrays
  INT shift_flag = 0;
  if(A->IA[0]==1) {
    dcsr_shift(A, -1);  // shift A
    shift_flag = 1;
  }
  // Transpose A
  dcsr_trans(A,&At);
  // Solve
  err_flag = hazmath_solve(&At,f,x,Numeric,print_level);
  // Fix A back to correct counting
  if(shift_flag==1) {
    dcsr_shift(A, 1);  // shift A back
  }
  // Clean up
  dcsr_free(&At);
#else
  err_flag = run_hazmath_solve(A,f,x,Numeric,print_level);
  //  error_extlib(252, __FUNCTION__, "SuiteSparse");
  //  return 0;
#endif
  // Error Check
  if(err_flag<0) {
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- DIRECT SOLVE ERROR!!!\n\n",__FUNCTION__);
    exit(err_flag);
  }
  return err_flag;
}
/***********************************************************************************************/
/**
 * \fn block_directsolve_HAZ(block_dCSRmat *bA,dvector *f,dvector *x,INT print_level)
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
 * \note modified by Ludmil on 20220802 to include hazmath direct
 *       solve option.
 *
 */
INT block_directsolve_HAZ(block_dCSRmat *bA,
                    dvector *f,
                    dvector *x,
                    INT print_level)
{
  INT err_flag;

  // Convert block matrix to regular matrix
  dCSRmat A = bdcsr_2_dcsr(bA);

  dcsr_compress_inplace(&A, 1.0e-14);

  // Call regular solve (either umpfpack or hazmath
  err_flag = directsolve_HAZ(&A,f,x,print_level);

  // Clean up
  dcsr_free(&A);

  return err_flag;
  //#if WITH_SUITESPARSE
  //#else
  //  error_extlib(252, __FUNCTION__, "SuiteSparse");
  //  return 0;
  //#endif
}

/***********************************************************************************************/
/**
 * \fn void* block_factorize_HAZ(block_dCSRmat *bA,INT print_level)
 *
 * \brief Performs factorization of A using UMFPACK (if suitesparse libraries are included) or HAZMath (Assumes A is in block_dCSR format)
 * \note This routine does factorization only.
 *
 * \param bA	       	        Matrix A to be solved (no need of transpose)
 *
 * \return Numeric	        Stores LU decomposition of A
 *
 * \note modified by Ludmil on 20220802 to include hazmath direct
 *       solve option.
 *
 */
void* block_factorize_HAZ(block_dCSRmat *bA,
                    INT print_level)
{

  // Data for Numerical factorization
  void* Numeric = NULL;

  // Convert block matrix to regular matrix
  dCSRmat A = bdcsr_2_dcsr(bA);

  // Factorize
  Numeric = factorize_HAZ(&A,print_level);

  // Cleanup
  dcsr_free(&A);

  return Numeric;
  //#if WITH_SUITESPARSE
  //#else
  //  error_extlib(252, __FUNCTION__, "SuiteSparse");
  //  return 0;
  //#endif
}

/***********************************************************************************************/
/**
 * \fn INT block_solve_HAZ(block_dCSRmat *bA,dvector *f,dvector *x, void* Numeric, INT print_level)
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
 * \note modified by Ludmil on 20220802 to include hazmath direct
 *       solve option.
 *
 */
INT block_solve_HAZ(block_dCSRmat *bA,
              dvector *f,
              dvector *x,
              void *Numeric,
              INT print_level)
{
  INT err_flag;

  // Convert block matrix to regular matrix
  dCSRmat A = bdcsr_2_dcsr(bA);

  // Call regular solve
  err_flag = solve_HAZ(&A,f,x,Numeric,print_level);

  // Clean up
  dcsr_free(&A);

  return err_flag;

}
// Suitesparse routines - assumes conversions have been done
// NOT ANYMORE: hazmath calls suitesparse only if compiled.
/*******************************************************************/
/**
 * \fn void* hazmath_factorize (dCSRmat *ptrA, const SHORT prtlvl)
 * \brief factorize A using either HAZmath or UMFpack
 *
 * \param ptrA Pointer to dCSRmat matrix (for UMFPACK transpose must
 *             have been done!)
 *
 * \param Numeric   Pointer to the numerical factorization
 *
 * \author Xiaozhe Hu
 * \date   10/20/2010
 * Modified by Xiaozhe Hu on 05/02/2014
 *
 * \note Modified by Ludmil (20220802)
 *
 */
void* hazmath_factorize (dCSRmat *ptrA,
                         const SHORT prtlvl)
{
  void *Numeric;
  clock_t start_time = clock();
#if WITH_SUITESPARSE
  const INT n = ptrA->col;

  INT *Ap = ptrA->IA;
  INT *Ai = ptrA->JA;
  double *Ax = ptrA->val;
  void *Symbolic;
  INT status = SUCCESS;
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
    fprintf(stdout,"\nUMFPACK: ");
  }
#else
  Numeric=run_hazmath_factorize(ptrA,(INT )prtlvl);
  //  error_extlib(253, __FUNCTION__, "SuiteSparse");
  //  return NULL;
  if ( prtlvl > PRINT_MIN ) {
    fprintf(stdout,"\nHAZMATH: ");
  }
#endif
  if ( prtlvl > PRINT_MIN ) {
    clock_t end_time = clock();
    double fac_time = (double)(end_time - start_time)/(double)(CLOCKS_PER_SEC);
    fprintf(stdout,"factorization costs %f seconds.\n", fac_time);
  }
  return Numeric;// in all cases. 
}

/***************************************************************************************************************************/
/**
 * \fn INT hazmath_solve (dCSRmat *ptrA, dvector *b, dvector *u,
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
 *
 * \modified Ludmil (20220802)
 *
 */
INT hazmath_solve (dCSRmat *ptrA,
                   dvector *b,
                   dvector *u,
                   void *Numeric,
                   const SHORT prtlvl)
{

  clock_t start_time = clock();
  INT status = SUCCESS;
#if WITH_SUITESPARSE
  INT *Ap = ptrA->IA;
  INT *Ai = ptrA->JA;
  double *Ax = ptrA->val;

  status = umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, u->val, b->val, Numeric, NULL, NULL);
  if ( prtlvl > PRINT_NONE ) {
    fprintf(stdout,"UMFPACK: ");fflush(stdout);
  }
#else
  //  error_extlib(254, __FUNCTION__, "SuiteSparse");
  //  return 0;
  if ( prtlvl > PRINT_NONE ) {
    fprintf(stdout,"HAZMATH: ");fflush(stdout);
  }
  status = run_hazmath_solve(ptrA, b, u, Numeric, (INT )prtlvl);
#endif
  if ( prtlvl > PRINT_NONE ) {
    clock_t end_time = clock();
    double solve_time = (double)(end_time - start_time)/(double)(CLOCKS_PER_SEC);
    printf("direct solve costs %f seconds.\n", solve_time);
  }

  return status;
}
/*****************************************************************/
/**
 * \fn INT hazmath_free_numeric (void **Numeric)
 * \brief Solve Au=b by UMFpack
 *
 * \param Numeric   double pointer to the numerical factorization
 *
 * \author Xiaozhe Hu
 * \date   10/20/2010
 * \note Modified by Ludmil (20220805)
 */
//NO: INT umfpack_free_numeric (void *Numeric)
INT hazmath_free_numeric (void **Numeric)
{
  INT status = SUCCESS;
#if WITH_SUITESPARSE
  
  umfpack_di_free_numeric (&Numeric[0]);
  
#else
  dCSRmat *U; dvector *dinv; SHORT *extra;
  hazmath_get_numeric(Numeric[0], &U, &dinv,&extra);
  //
  dcsr_free(U);
  dvec_free(dinv);
  free(Numeric[0]);
  Numeric[0]=NULL;
  //  error_extlib(255, __FUNCTION__, "SuiteSparse");
  //  return 0;
#endif
  return (INT )status; 
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
