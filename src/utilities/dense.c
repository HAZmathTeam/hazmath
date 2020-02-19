/*! \file src/utilities/dense.c
 *
 *  Created by Xiaozhe Hu on 01/04/20.
 *  Copyright 2020__HAZMATH__. All rights reserved.
 *
 *
 */

#include "hazmath.h"

/***********************************************************************************************/
/*!
 * \fn dDENSEmat ddense_create (const INT n, const INT m)
 *
 * \brief Create a dDENSEmat dense matrix
 *
 * \param n    Number of rows
 * \param m    Number of columns
 *
 * \return A   the dDENSEmat matrix
 *
 * \author Xiaozhe Hu
 * \data   01/04/2020
 */
dDENSEmat ddense_create(const INT n,
                        const INT m)
{
  // local variable
  dDENSEmat A;

  // set size
  A.row = n; A.col = m;

  // allocate
  if ( n*m > 0 ) {
    A.val = (REAL *)calloc(n*m, sizeof(REAL));
  }
  else {
    A.val = NULL;
  }

  // return
  return A;

}

/***********************************************************************************************/
/*!
 * \fn void ddense_alloc(const INT n, const INT m, dDENSEmat *A)
 *
 * \brief Allocate dDENSEmat dense matrix memory space
 *
 * \param n      Number of rows
 * \param m      Number of columns
 * \param A      Pointer to the dDENSEmat matrix
 *
 */
void ddense_alloc(const INT n,
                  const INT m,
                  dDENSEmat *A)
{

  // set size
  A->row=n; A->col=m;

  // allocate
  if ( n*m > 0 ) {
    A->val=(REAL*)calloc(n*m,sizeof(REAL));
  }
  else {
    A->val = NULL;
  }

  // return
  return;
}


/***********************************************************************************************/
/*!
 * \fn void ddense_set(dDENSEmat *A, REAL val)
 *
 * \brief Initialize dDENSEmat A and set all the entries to be a given value
 *
 * \param A      Pointer to dDENSEmat
 * \param val    Given REAL value
 *
 * \note Memory space has to be allocated before calling this function -- Xiaozhe Hu
 *
 */
void ddense_set(dDENSEmat *A,
                REAL val)
{
    // local variables
    INT i;
    const INT n = A->row;
    const INT m = A->col;
    const INT nnz = n*m;

    if (val == 0.0) {
        memset(A->val, 0x0, sizeof(REAL)*nnz);
    }
    else {
        for (i=0; i<nnz; ++i) A->val[i]=val;
    }
}

/***********************************************************************************************/
/*!
 * \fn void ddense_free(dDENSEmat *A)
 *
 * \brief Free dDENSEmat A
 *
 * \param A      Pointer to dDENSEmat
 *
 */
void ddense_free(dDENSEmat *A)
{

  if (A->val==NULL) return;

  free(A->val);
  A->row = 0; A->col = 0; A->val = NULL;

}


/***********************************************************************************************/
/*!
 * \fn void ddense_cp(dDENSEmat *A, dDENSEmat *B)
 *
 * \brief copy dDENSEmat A and to dDENSEmat B (B=A)
 *
 * \param A  Pointer to dDENSEmat
 * \param B  Pointer to dDENSEmat
 *
 * \note Memory space has to be allocated before calling this function -- Xiaozhe Hu
 *
 */
void ddense_cp(dDENSEmat *A,
               dDENSEmat *B)
{

    // copy size
    B->row = A->row;
    B->col = A->col;

    // copy
    memcpy(B->val, A->val, (A->row)*(A->col)*sizeof(REAL));

}

/***********************************************************************************************/
/*!
 * \fn dDENSEmat ddense_random_JL(const INT k, const INT d)
 *
 * \brief generate a dense random matrix that satisifes Johnson-Lindenstrauss Lemma
 *
 * \param k   low dimension
 * \param d   high dimension
 *
 */
dDENSEmat ddense_random_JL(const INT k,
                           const INT d)
{
    // local variables
    INT i;
    const INT nnz = k*d;
    REAL random_number;

    dDENSEmat Q = ddense_create(k, d);

    // uncomment this if you really need randomness
    //srand ( time(NULL) );

    // main loop
    for (i=0; i<nnz; i++){

      // generate random number between 0 and 1
      random_number = ((double) rand() / (RAND_MAX));

      if (random_number > 0.5){
        Q.val[i] = 1.0/sqrt(k);
      }
      else {
        Q.val[i] = -1.0/sqrt(k);
      }

    }

    // return
    return Q;

}

/***********************************************************************************************/
/*!
* \fn void find_det_4( REAL* A, REAL deta)
*
* \brief find det of 4x4 matrix using expansion by minors once, 
			then direct computation of det of 3x3 matrices
*
* \param A            vectorized matrix (row-wise)
*
* \return deta        determinant
*
*/

void find_det_4( REAL* A, REAL* deta)
{

  REAL s0p, s0m, s1p, s1m, s2p, s2m, s3p, s3m;

  s0p = A[5]*A[10]*A[15] + A[6]*A[11]*A[13] + A[7]*A[9]*A[14];
  s0m = A[13]*A[10]*A[7] + A[14]*A[11]*A[5] + A[15]*A[9]*A[6];

  s1p = A[4]*A[10]*A[15] + A[6]*A[11]*A[12] + A[7]*A[8]*A[14];
  s1m = A[12]*A[10]*A[7] + A[14]*A[11]*A[4] + A[15]*A[8]*A[6];

  s2p = A[4]*A[9]*A[15] + A[5]*A[11]*A[12] + A[7]*A[8]*A[13];
  s2m = A[12]*A[9]*A[7] + A[13]*A[11]*A[4] + A[15]*A[8]*A[5];

  s3p = A[4]*A[9]*A[14] + A[5]*A[10]*A[12] + A[6]*A[8]*A[13];
  s3m = A[12]*A[9]*A[6] + A[13]*A[10]*A[4] + A[14]*A[8]*A[5];

  *deta = A[0]*(s0p-s0m) - A[1]*(s1p-s1m) + A[2]*(s2p-s2m) - A[3]*(s3p-s3m);

  return;

}





/************************************** END ***************************************************/
