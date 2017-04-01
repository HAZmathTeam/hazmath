/*! \file src/utilities/sparse.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 3/6/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/30/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/31/2016
 *
 */

#include "hazmath.h"

/***********************************************************************************************/
/*!
 * \fn dCSRmat dcsr_create (const INT m, const INT n, const INT nnz)
 *
 * \brief Create a dCSRmat sparse matrix
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 *
 * \return A   the dCSRmat matrix
 *
 */
dCSRmat dcsr_create (const INT m,
                     const INT n,
                     const INT nnz)
{
  dCSRmat A;

  if ( m > 0 ) {
    A.IA = (INT *)calloc(m+1, sizeof(INT));
  }
  else {
    A.IA = NULL;
  }

  if ( n > 0 ) {
    A.JA = (INT *)calloc(nnz, sizeof(INT));
  }
  else {
    A.JA = NULL;
  }

  if ( nnz > 0 ) {
    A.val = (REAL *)calloc(nnz, sizeof(REAL));
  }
  else {
    A.val = NULL;
  }

  A.row = m; A.col = n; A.nnz = nnz;

  return A;
}

/***********************************************************************************************/
/*!
 * \fn dCSRmat dcsr_create_zeromatrix (const INT m, const INT n, const INT nnz)
 *
 * \brief Create a CSR sparse matrix that is all zeros
 *
 * \param  m             Number of rows
 * \param  n             Number of columns
 * \param  index_start   Number from which memory is indexed (1 for fem/assembly, 0 for solver)
 *
 * \return A             the new dCSRmat matrix
 *
 */
dCSRmat dcsr_create_zeromatrix(const INT m,
                               const INT n,
                               const INT index_start)
{
  dCSRmat A;

  A.IA = (INT *)calloc(m+1, sizeof(INT));
  A.JA = (INT *)calloc(1, sizeof(INT));
  A.val = (REAL *)calloc(1, sizeof(REAL));

  A.row = m; A.col = n; A.nnz = 1;

  INT i;
  A.IA[0]=index_start;
  for(i=1;i<m+1;i++) A.IA[i]=index_start+1;
  A.JA[0]=index_start;
  A.val[0] = 0.0;

  return A;
}

/***********************************************************************************************/
/**
 * \fn dCSRmat dcsr_create_single_nnz_matrix (const INT m, const INT n, const INT row,
 *                           const INT col, const REAL val, const INT index_start)
 *
 * \brief Create a dCSRmat sparse matrix that is all zeros except for one non zero element
 *
 * \param m             Number of rows
 * \param n             Number of columns
 * \param row           Row index of nonzero value
 * \param col           Column index of nonzero value
 * \param val           Value of nonzero value
 * \param index_start   Number from which memory is indexed (1 for fem/assembly, 0 for solver)
 *
 * \return A   the new dCSRmat matrix
 *
 * \note row and col are indexed consistently with index_start.
 *
 */
dCSRmat dcsr_create_single_nnz_matrix(const INT m,
                                      const INT n,
                                      const INT row,
                                      const INT col,
                                      const REAL val,
                                      const INT index_start)
{
  dCSRmat A;

  A.IA = (INT *)calloc(m+1, sizeof(INT));
  A.JA = (INT *)calloc(1, sizeof(INT));
  A.val = (REAL *)calloc(1, sizeof(REAL));

  A.row = m; A.col = n; A.nnz = 1;

  INT i;
  for(i=0;i<row+1-index_start;i++) A.IA[i]=index_start;
  for(i=row+1-index_start;i<m+1;i++) A.IA[i]=index_start+1;
  A.JA[0]=col;
  A.val[0] = val;

  return A;
}

/***********************************************************************************************/
/*!
 * \fn dCSRmat dcsr_create_identity_matrix (const INT m, const INT n, const INT index_start)
 *
 * \brief Create a dCSRmat sparse matrix that is the identity matrix
 *
 * \param m             Number of rows
 * \param index_start   Number from which memory is indexed (1 for fem/assembly, 0 for solver)
 *
 * \return A            the new dCSRmat matrix
 *
 */
dCSRmat dcsr_create_identity_matrix(const INT m,
                                    const INT index_start)
{
  dCSRmat A;

  A.IA = (INT *)calloc(m+1, sizeof(INT));
  A.JA = (INT *)calloc(m, sizeof(INT));
  A.val = (REAL *)calloc(m, sizeof(REAL));

  A.row = m; A.col = m; A.nnz = m;

  INT i;
  for(i=0;i<m;i++){
    A.IA[i]=i + index_start;
    A.JA[i]=i + index_start;
    A.val[i] = 1.0;
  }
  A.IA[m] = m + index_start;

  return A;
}

/***********************************************************************************************/
/*!
 * \fn void dcsr_alloc (const INT m, const INT n, const INT nnz, dCSRmat *A)
 *
 * \brief Allocate dCSRmat sparse matrix memory space
 *
 * \param m      Number of rows
 * \param n      Number of columns
 * \param nnz    Number of nonzeros
 * \param A      Pointer to the dCSRmat matrix
 *
 */
void dcsr_alloc(const INT m,
                const INT n,
                const INT nnz,
                dCSRmat *A)
{
  if ( m > 0 ) {
    A->IA=(INT*)calloc(m+1,sizeof(INT));
  }
  else {
    A->IA = NULL;
  }

  if ( n > 0 ) {
    A->JA=(INT*)calloc(nnz,sizeof(INT));
  }
  else {
    A->JA = NULL;
  }

  if ( nnz > 0 ) {
    A->val=(REAL*)calloc(nnz,sizeof(REAL));
  }
  else {
    A->val = NULL;
  }

  A->row=m; A->col=n; A->nnz=nnz;

  return;
}

/***********************************************************************************************/
/**
 * \fn iCSRmat icsr_create (const INT m, const INT n, const INT nnz)
 *
 * \brief Create iCSRmat sparse matrix
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 *
 * \return A   the new iCSRmat matrix
 *
 */
iCSRmat icsr_create(const INT m,
                    const INT n,
                    const INT nnz)
{
  iCSRmat A;

  if ( m > 0 ) {
    A.IA = (INT *)calloc(m+1, sizeof(INT));
  }
  else {
    A.IA = NULL;
  }

  if ( n > 0 ) {
    A.JA = (INT *)calloc(nnz, sizeof(INT));
  }
  else {
    A.JA = NULL;
  }

  if ( nnz > 0 ) {
    A.val = (INT *)calloc(nnz, sizeof(INT));
  }
  else {
    A.val = NULL;
  }

  A.row = m; A.col = n; A.nnz = nnz;

  return A;
}

/***********************************************************************************************/
/**
 * \fn iCSRmat icsr_create_identity (const INT m, const INT index_start)
 *
 * \brief Create a iCSRmat sparse matrix that is the identity matrix
 *
 * \param m             Number of rows
 * \param index_start   Number from which memory is indexed (1 for fem/assembly, 0 for solver)
 *
 * \return A            the new dCSRmat matrix
 *
 */
iCSRmat icsr_create_identity(const INT m,
                             const INT index_start)
{
     INT i;
  iCSRmat A;

  if ( m > 0 ) {
    A.IA = (INT *)calloc(m+1, sizeof(INT));
    A.JA = (INT *)calloc(m, sizeof(INT));
    A.val = NULL;
  }
  else {
    A.IA = NULL;
    A.JA = NULL;
    A.val = NULL;
  }

  A.row = m; A.col = m; A.nnz = m;

  for(i=0;i<m;i++) {
    A.IA[i] = i+index_start;
    A.JA[i] = i+index_start;
  }
  A.IA[m] = m+index_start;

  return A;
}

/***********************************************************************************************/
/*!
 * \fn void dcsr_free (dCSRmat *A)
 *
 * \brief Free dCSRmat sparse matrix
 *
 * \param A   Pointer to the dCSRmat matrix
 *
 */
void dcsr_free(dCSRmat *A)
{
  if ( A == NULL ) return;

  if (A->IA) {
    free(A->IA);
    A->IA  = NULL;
  }

  if (A->JA) {
    free(A->JA);
    A->JA  = NULL;
  }

  if (A->val) {
    free(A->val);
    A->val = NULL;
  }

}

/***********************************************************************************************/
/*!
 * \fn void icsr_free (iCSRmat *A)
 *
 * \brief Free iCSRmat sparse matrix
 *
 * \param A   Pointer to the iCSRmat matrix
 *
 */
void icsr_free(iCSRmat *A)
{
  if ( A == NULL ) return;

  if (A->IA) {
    free(A->IA);
    A->IA  = NULL;
  }


  if (A->JA) {
    free(A->JA);
    A->JA  = NULL;
  }

  if (A->val) {
    free(A->val);
    A->val = NULL;
  }
}

/***********************************************************************************************/
/*!
 * \fn void dcsr_null (dCSRmat *A)
 *
 * \brief Initialize dCSRmat sparse matrix (set arrays to NULL)
 *
 * \param A   Pointer to the dCSRmat matrix
 *
 */
void dcsr_null(dCSRmat *A)
{
  A->row = A->col = A->nnz = 0;
  A->IA  = A->JA  = NULL;
  A->val = NULL;
}

/***********************************************************************************************/
/*!
 * \fn void icsr_null (iCSRmat *A)
 *
 * \brief Initialize iCSRmat sparse matrix (set arrays to NULL)
 *
 * \param A   Pointer to the iCSRmat matrix
 *
 */
void icsr_null(iCSRmat *A)
{
  A->row = A->col = A->nnz = 0;
  A->IA  = A->JA  = NULL;
  A->val = NULL;
}

/***********************************************************************************************/
/*!
 * \fn dCSRmat dcsr_perm (dCSRmat *A, INT *P)
 *
 * \brief Apply permutation to a dCSRmat matrix A, i.e. Aperm=PAP', for a given ordering P
 *
 * \param A  Pointer to the original dCSRmat matrix
 * \param P  Pointer to orders
 *
 * \return   The new ordered dCSRmat matrix if succeed, NULL if fail
 *
 * \note   P[i] = k means k-th row and column become i-th row and column!
 *
 */
dCSRmat dcsr_perm(dCSRmat *A,
                  INT *P)
{
  const INT n=A->row,nnz=A->nnz;
  const INT *ia=A->IA, *ja=A->JA;
  const REAL *Aval=A->val;
  INT i,j,k,jaj,i1,i2,start;

  dCSRmat Aperm = dcsr_create(n,n,nnz);

  // form the transpose of P
  INT *Pt = (INT*)calloc(n,sizeof(INT));

  for (i=0; i<n; ++i) Pt[P[i]] = i;

  // compute IA of P*A (row permutation)
  Aperm.IA[0] = 0;
  for (i=0; i<n; ++i) {
    k = P[i];
    Aperm.IA[i+1] = Aperm.IA[i]+(ia[k+1]-ia[k]);
  }

  // perform actual P*A
  for (i=0; i<n; ++i) {
    i1 = Aperm.IA[i]; i2 = Aperm.IA[i+1]-1;
    k = P[i];
    start = ia[k];
    for (j=i1; j<=i2; ++j) {
      jaj = start+j-i1;
      Aperm.JA[j] = ja[jaj];
      Aperm.val[j] = Aval[jaj];
    }
  }

  // perform P*A*P' (column permutation)
  for (k=0; k<nnz; ++k) {
    j = Aperm.JA[k];
    Aperm.JA[k] = Pt[j];
  }

  free(Pt);

  return(Aperm);
}

/***********************************************************************************************/
/*!
 * \fn void icsr_cp (iCSRmat *A, iCSRmat *B)
 *
 * \brief Copy a iCSRmat to a new one B=A
 *
 * \param A   Pointer to the iCSRmat matrix
 * \param B   Pointer to the iCSRmat matrix
 *
 */
void icsr_cp(iCSRmat *A,
             iCSRmat *B)
{
  B->row=A->row;
  B->col=A->col;
  B->nnz=A->nnz;

  iarray_cp (A->row+1, A->IA, B->IA);
  iarray_cp (A->nnz, A->JA, B->JA);
  iarray_cp (A->nnz, A->val, B->val);
}

/***********************************************************************************************/
/*!
 * \fn void dcsr_cp (dCSRmat *A, dCSRmat *B)
 *
 * \brief copy a dCSRmat to a new one B=A
 *
 * \param A   Pointer to the dCSRmat matrix
 * \param B   Pointer to the dCSRmat matrix
 *
 */
void dcsr_cp(dCSRmat *A,
             dCSRmat *B)
{
  B->row=A->row;
  B->col=A->col;
  B->nnz=A->nnz;

  iarray_cp(A->row+1, A->IA, B->IA);
  iarray_cp(A->nnz, A->JA, B->JA);
  array_cp(A->nnz, A->val, B->val);
}

/***********************************************************************************************/
/*!
 * \fn void dcsr_trans (dCSRmat *A, dCSRmat *AT)
 *
 * \brief Transpose a dCSRmat matrix A
 *
 * \param A   Pointer to the dCSRmat matrix
 * \param AT  Pointer to the transpose of dCSRmat matrix A (output)
 *
 * \note Index of the arrays start at 0!!! -- Xiaozhe Hu
 *
 */
INT dcsr_trans(dCSRmat *A,
               dCSRmat *AT)
{
  const INT n=A->row, m=A->col, nnz=A->nnz;

  // Local variables
  INT i,j,k,p;

  AT->row=m;
  AT->col=n;
  AT->nnz=nnz;

  AT->IA=(INT*)calloc(m+1,sizeof(INT));

  AT->JA=(INT*)calloc(nnz,sizeof(INT));

  if (A->val) {
    AT->val=(REAL*)calloc(nnz,sizeof(REAL));

  }
  else { AT->val=NULL; }

  // first pass: find the Number of nonzeros in the first m-1 columns of A
  // Note: these Numbers are stored in the array AT.IA from 1 to m-1

  memset(AT->IA, 0, sizeof(INT)*(m+1));

  for (j=0;j<nnz;++j) {
    i=A->JA[j]; // column Number of A = row Number of A'
    if (i<m-1) AT->IA[i+2]++;
  }

  for (i=2;i<=m;++i) AT->IA[i]+=AT->IA[i-1];

  // second pass: form A'
  if (A->val) {
    for (i=0;i<n;++i) {
      INT ibegin=A->IA[i], iend=A->IA[i+1];
      for (p=ibegin;p<iend;p++) {
        j=A->JA[p]+1;
        k=AT->IA[j];
        AT->JA[k]=i;
        AT->val[k]=A->val[p];
        AT->IA[j]=k+1;
      } // end for p
    } // end for i
  }
  else {
    for (i=0;i<n;++i) {
      INT ibegin=A->IA[i], iend1=A->IA[i+1];
      for (p=ibegin;p<iend1;p++) {
        j=A->JA[p]+1;
        k=AT->IA[j];
        AT->JA[k]=i;
        AT->IA[j]=k+1;
      } // end for p
    } // end of i
  } // end if

  return 0;
}

/***********************************************************************************************/
/*!
 * \fn void dcsr_trans_1 (dCSRmat *A, dCSRmat *AT)
 *
 * \brief Transpose a dCSRmat matrix A
 *
 * \param A   Pointer to the dCSRmat matrix
 * \param AT  Pointer to the transpose of dCSRmat matrix A (output)
 *
 * \note  Index of the arrays start at 1!!! -- Xiaozhe Hu
 *
 */
void dcsr_trans_1(dCSRmat *A,
                  dCSRmat *AT)
{
  // shift A
  dcsr_shift(A, -1);

  // transpose
  dcsr_trans(A, AT);

  // shift A and AT back
  dcsr_shift(A, 1);
  dcsr_shift(AT, 1);

}

/***********************************************************************************************/
/*!
 * \fn void icsr_trans (iCSRmat *A, iCSRmat *AT)
 *
 * \brief Transpose an iCSRmat matrix A
 *
 * \param A   Pointer to the iCSRmat matrix
 * \param AT  Pointer to the transpose of iCSRmat matrix A (output)
 *
 * \note Index of the arrays start at 0!!! -- Xiaozhe Hu
 *
 */
void icsr_trans(iCSRmat *A,
                iCSRmat *AT)
{
  const INT n=A->row, m=A->col, nnz=A->nnz, m1=m-1;

  // Local variables
  INT i,j,k,p;
  INT ibegin, iend;

  AT->row=m; AT->col=n; AT->nnz=nnz;

  AT->IA=(INT*)calloc(m+1,sizeof(INT));

  AT->JA=(INT*)calloc(nnz,sizeof(INT));

  if (A->val) {
    AT->val=(INT*)calloc(nnz,sizeof(INT));
  }
  else {
    AT->val=NULL;
  }

  // first pass: find the Number of nonzeros in the first m-1 columns of A
  // Note: these Numbers are stored in the array AT.IA from 1 to m-1
  memset(AT->IA, 0, sizeof(INT)*(m+1));

  for (j=0;j<nnz;++j) {
    i=A->JA[j]; // column Number of A = row Number of A'
    if (i<m1) AT->IA[i+2]++;
  }

  for (i=2;i<=m;++i) AT->IA[i]+=AT->IA[i-1];

  // second pass: form A'
  if (A->val != NULL) {
    for (i=0;i<n;++i) {
      ibegin=A->IA[i], iend=A->IA[i+1];
      for (p=ibegin;p<iend;p++) {
        j=A->JA[p]+1;
        k=AT->IA[j];
        AT->JA[k]=i;
        AT->val[k]=A->val[p];
        AT->IA[j]=k+1;
      } // end for p
    } // end for i
  }
  else {
    for (i=0;i<n;++i) {
      ibegin=A->IA[i], iend=A->IA[i+1];
      for (p=ibegin;p<iend;p++) {
        j=A->JA[p]+1;
        k=AT->IA[j];
        AT->JA[k]=i;
        AT->IA[j]=k+1;
      } // end for p
    } // end for i
  } // end if
}

/***********************************************************************************************/
/*!
 * \fn void icsr_trans_1 (iCSRmat *A, iCSRmat *AT)
 *
 * \brief Transpose an iCSRmat matrix A
 *
 * \param A   Pointer to the iCSRmat matrix
 * \param AT  Pointer to the transpose of iCSRmat matrix A (output)
 *
 * \note  Index of the arrays start at 1!!! -- Xiaozhe Hu
 *
 */
void icsr_trans_1(iCSRmat *A,
                  iCSRmat *AT)
{
  // shift A and B
  icsr_shift(A, -1);

  icsr_trans(A, AT);

  // shift A and AT back
  icsr_shift(A, 1);
  icsr_shift(AT, 1);


}

/***********************************************************************************************/
/**
 * \fn void dcsr_compress (dCSRmat *A, dCSRmat *B, REAL dtol)
 *
 * \brief Compress a CSR matrix A and store in CSR matrix B by
 *        dropping small entries such that abs(aij)<=dtol
 *
 * \param A     Pointer to dCSRmat CSR matrix
 * \param B     Pointer to dCSRmat CSR matrix (OUTPUT)
 * \param dtol  Drop tolerance
 *
 */
void dcsr_compress(dCSRmat *A,
                   dCSRmat *B,
                   REAL dtol)
{
    // local variables
    INT i, j, k;
    INT ibegin,iend1;

    // allocate
    INT *index=(INT*)calloc(A->nnz,sizeof(INT));

    B->row=A->row; B->col=A->col;
    B->IA=(INT*)calloc(A->row+1,sizeof(INT));
    B->IA[0]=A->IA[0];

    // first step: determine the size of B
    k=0;
    for (i=0;i<A->row;++i) {
        ibegin=A->IA[i]; iend1=A->IA[i+1];
        for (j=ibegin;j<iend1;++j)
            if (ABS(A->val[j])>dtol) {
                index[k]=j;
                ++k;
            } /* end of j */
        B->IA[i+1]=k;
    } /* end of i */
    B->nnz=k;

    // allocate
    B->JA=(INT*)calloc(B->nnz,sizeof(INT));
    B->val=(REAL*)calloc(B->nnz,sizeof(REAL));

    // second step: generate the index and element to B
    for (i=0;i<B->nnz;++i) {
        B->JA[i]=A->JA[index[i]];
        B->val[i]=A->val[index[i]];
    }

    free(index);
}

/***********************************************************************************************/
/**
 * \fn SHORT dcsr_compress_inplace (dCSRmat *A, REAL dtol)
 *
 * \brief Compress a CSR matrix A by dropping small entries such abs(aij)<=dtol (still stores in A)
 *
 * \param A     Pointer to dCSRmat CSR matrix (OUTPUT)
 * \param dtol  Drop tolerance
 *
 */
SHORT dcsr_compress_inplace(dCSRmat *A,
                            REAL dtol)
{
    // local variables
    const INT row=A->row;
    const INT nnz=A->nnz;

    INT i, j, k;
    INT ibegin, iend = A->IA[0];
    SHORT status = SUCCESS;

    k = 0;
    for ( i=0; i<row; ++i ) {
        ibegin = iend; iend = A->IA[i+1];
        for ( j=ibegin; j<iend; ++j )
            if ( ABS(A->val[j]) > dtol ) {
                A->JA[k]  = A->JA[j];
                A->val[k] = A->val[j];
                ++k;
            } /* end of j */
        A->IA[i+1] = k;
    } /* end of i */

    if ( k <= nnz ) {
        A->nnz=k;
        A->JA  = (INT  *)realloc(A->JA,  k*sizeof(INT));
        A->val = (REAL *)realloc(A->val, k*sizeof(REAL));
    }
    else {
        printf("### HAZMATH ERROR: Size of compressed matrix is larger than the original!\n");
        status = ERROR_UNKNOWN;
    }

    return (status);
}

/***********************************************************************************************/
/*!
 * \fn void dcsr_shift (dCSRmat *A, INT offset)
 *
 * \brief Re-index a dCSRmat matrix to make the index starting from 0 or 1
 *
 * \param A         Pointer to dCSRmat matrix
 * \param offset    Size of offset (1 or -1)
 *
 */
void dcsr_shift(dCSRmat *A,
                INT offset)
{
  if (A == NULL) return;

  const INT nnz=A->nnz;
  const INT n=A->row+1;
  INT i, *ai=A->IA, *aj=A->JA;

  for (i=0; i<n; ++i) ai[i]+=offset;

  for (i=0; i<nnz; ++i) aj[i]+=offset;

  return;

}

/***********************************************************************************************/
/*!
 * \fn void icsr_shift (iCSRmat *A, INT offset)
 *
 * \brief Re-index an iCSRmat matrix to make the index starting from 0 or 1
 *
 * \param A         Pointer to iCSRmat matrix
 * \param offset    Size of offset (1 or -1)
 *
 */
void icsr_shift(iCSRmat *A,
                INT offset)
{
  const INT nnz=A->nnz;
  const INT n=A->row+1;
  INT i, *ai=A->IA, *aj=A->JA;

  for (i=0; i<n; ++i) ai[i]+=offset;

  for (i=0; i<nnz; ++i) aj[i]+=offset;

}

/***********************************************************************************************/
/*!
 * \fn void dcsr_add (dCSRmat *A, const REAL alpha, dCSRmat *B,
 *                              const REAL beta, dCSRmat *C)
 *
 * \brief compute C = alpha*A + beta*B in dCSRmat format
 *
 * \param A      Pointer to dCSRmat matrix
 * \param alpha  REAL factor alpha
 * \param B      Pointer to dCSRmat matrix
 * \param beta   REAL factor beta
 * \param C      Pointer to dCSRmat matrix (OUTPUT)
 *
 * \return Flag of whether the adding is succesful or not (SUCCUESS: 0; FAIL: <0)
 *
 */
INT dcsr_add(dCSRmat *A,
             const REAL alpha,
             dCSRmat *B,
             const REAL beta,
             dCSRmat *C)
{
  INT i,j,k,l;
  INT count=0, added, countrow;
  INT status = SUCCESS;

  // both matrices A and B are NULL
  if (A == NULL && B == NULL) {
    C->row=0; C->col=0; C->nnz=0;
    status=SUCCESS;
    goto FINISHED;
  }

  // matrix A is NULL but B is not
  if (A == NULL) {

     // matrix B is empty
     if (B->nnz == 0) {
         C->row=B->row; C->col=B->col; C->nnz=0;
         status=SUCCESS;
         goto FINISHED;
     }
     // matrix B is not empty
     else {
         dcsr_alloc(B->row,B->col,B->nnz,C);
         memcpy(C->IA,B->IA,(B->row+1)*sizeof(INT));
         memcpy(C->JA,B->JA,(B->nnz)*sizeof(INT));

         for (i=0;i<B->nnz;++i) C->val[i]=B->val[i]*beta;

         status = SUCCESS;
         goto FINISHED;
     }

  }

  // matrix B is NULL but A is not
  if (B == NULL) {

      // matrix A is empty
     if (A->nnz == 0) {
         C->row=A->row; C->col=A->col; C->nnz=0;
         status=SUCCESS;
         goto FINISHED;
     }
     // matrix A is not empty
     else {
         dcsr_alloc(A->row,A->col,A->nnz,C);
         memcpy(C->IA,A->IA,(A->row+1)*sizeof(INT));
         memcpy(C->JA,A->JA,(A->nnz)*sizeof(INT));

         for (i=0;i<A->nnz;++i) C->val[i]=A->val[i]*alpha;

         status = SUCCESS;
         goto FINISHED;
     }

  }

  // neither matrix A or B is NULL
  // size does not match!
  if (A->row != B->row || A->col != B->col) {
    printf("### ERROR HAZMATH DANGER: Dimensions of matrices do not match!!! %s\n", __FUNCTION__);
    status = ERROR_MAT_SIZE;
    goto FINISHED;
  }

  // both matrices A and B are empty
  if (A->nnz == 0 && B->nnz == 0) {
    C->row=A->row; C->col=A->col; C->nnz=A->nnz;
    status = SUCCESS;
    goto FINISHED;
  }

  // empty matrix A
  if (A->nnz == 0) {
    dcsr_alloc(B->row,B->col,B->nnz,C);
    memcpy(C->IA,B->IA,(B->row+1)*sizeof(INT));
    memcpy(C->JA,B->JA,(B->nnz)*sizeof(INT));

    for (i=0;i<B->nnz;++i) C->val[i]=B->val[i]*beta;

    status = SUCCESS;
    goto FINISHED;
  }

  // empty matrix B
  if (B->nnz == 0) {
    dcsr_alloc(A->row,A->col,A->nnz,C);
    memcpy(C->IA,A->IA,(A->row+1)*sizeof(INT));
    memcpy(C->JA,A->JA,(A->nnz)*sizeof(INT));

    for (i=0;i<A->nnz;++i) C->val[i]=A->val[i]*alpha;

    status = SUCCESS;
    goto FINISHED;
  }

  // Both matrices A and B are neither NULL or empty
  C->row=A->row; C->col=A->col;

  C->IA=(INT*)calloc(C->row+1,sizeof(INT));

  // allocate work space for C->JA and C->val
  C->JA=(INT *)calloc(A->nnz+B->nnz,sizeof(INT));

  C->val=(REAL *)calloc(A->nnz+B->nnz,sizeof(REAL));

  // initial C->IA
  memset(C->IA, 0, sizeof(INT)*(C->row+1));
  memset(C->JA, -1, sizeof(INT)*(A->nnz+B->nnz));

  for (i=0; i<A->row; ++i) {
    countrow = 0;
    for (j=A->IA[i]; j<A->IA[i+1]; ++j) {
      C->val[count] = alpha * A->val[j];
      C->JA[count] = A->JA[j];
      C->IA[i+1]++;
      count++;
      countrow++;
    } // end for js

    for (k=B->IA[i]; k<B->IA[i+1]; ++k) {
      added = 0;

      for (l=C->IA[i]; l<C->IA[i]+countrow+1; l++) {
        if (B->JA[k] == C->JA[l]) {
          C->val[l] = C->val[l] + beta * B->val[k];
          added = 1;
          break;
        }
      } // end for l

      if (added == 0) {
        C->val[count] = beta * B->val[k];
        C->JA[count] = B->JA[k];
        C->IA[i+1]++;
        count++;
      }

    } // end for k

    C->IA[i+1] += C->IA[i];

  }

  C->nnz = count;
  C->JA  = (INT *)realloc(C->JA, (count)*sizeof(INT));
  C->val = (REAL *)realloc(C->val, (count)*sizeof(REAL));

FINISHED:
  return status;
}

/***********************************************************************************************/
/*!
 * \fn void dcsr_add_1 (dCSRmat *A, const REAL alpha, dCSRmat *B,
 *                              const REAL beta, dCSRmat *C)
 *
 * \brief Compute C = alpha*A + beta*B in dCSRmat format (index starts at 1!!!)
 *
 * \param A      Pointer to dCSRmat matrix
 * \param alpha  REAL factor alpha
 * \param B      Pointer to dCSRmat matrix
 * \param beta   REAL factor beta
 * \param C      Pointer to dCSRmat matrix
 *
 * \return Flag of whether the adding is succesful or not (SUCCUESS: 0; FAIL: <0)
 *
 * \note Index starts at 1!!! -- Xiaozhe Hu
 *
 */
INT dcsr_add_1(dCSRmat *A,
               const REAL alpha,
               dCSRmat *B,
               const REAL beta,
               dCSRmat *C)
{
  INT status = SUCCESS;

  // shift A
  dcsr_shift(A, -1);
  // shift B
  dcsr_shift(B, -1);

  // add
  status = dcsr_add(A,alpha,B,beta,C);

  // shift A back
  dcsr_shift(A, 1);
  // shift B back
  dcsr_shift(B, 1);
  // shift C back
  dcsr_shift(C, 1);

  return status;

}

/***********************************************************************************************/
/*!
 * \fn void dcsr_axm (dCSRmat *A, const REAL alpha)
 *
 * \brief Multiply a dCSRmat format sparse matrix A by a scalar number alpha.
 *
 * \param A      Pointer to dCSRmat matrix A
 * \param alpha  Scalar REAL number
 *
 * \note This works for both index cases (starts at 0 or 1) -- Xiaozhe Hu
 *
 */
void dcsr_axm(dCSRmat *A,
              const REAL alpha)
{
  const INT nnz=A->nnz;

  // A direct calculation can be written as:
  array_ax(nnz, alpha, A->val);

}


/***********************************************************************************************/
/*!
 * \fn void dcsr_mxv (dCSRmat *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x (index starts at 0!!)
 *
 * \param A   Pointer to dCSRmat matrix A
 * \param x   Pointer to array x
 * \param y   Pointer to array y
 *
 * \note Index starts at 0!!! -- Xiaozhe Hu
 *
 */
void dcsr_mxv(dCSRmat *A,
              REAL *x,
              REAL *y)
{
  const INT m=A->row;
  const INT *ia=A->IA, *ja=A->JA;
  const REAL *aj=A->val;
  INT i, k, begin_row, end_row, nnz_num_row;
  register REAL temp;

  for (i=0;i<m;++i) {
    temp=0.0;
    begin_row=ia[i];
    end_row=ia[i+1];
    nnz_num_row = end_row - begin_row;
    switch(nnz_num_row) {
    case 1:
      k=begin_row;
      temp+=aj[k]*x[ja[k]];
      break;
    case 2:
      k=begin_row;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      break;
    case 3:
      k=begin_row;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      break;
    case 4:
      k=begin_row;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      break;
    case 5:
      k=begin_row;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      break;
    case 6:
      k=begin_row;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      break;
    case 7:
      k=begin_row;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      break;
    case 8:
      k=begin_row;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      break;
    case 9:
      k=begin_row;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      k ++;
      temp+=aj[k]*x[ja[k]];
      break;
    default:
      for (k=begin_row; k<end_row; ++k) {
        temp+=aj[k]*x[ja[k]];
      }
      break;
    }

    y[i]=temp;

  }
}

/***********************************************************************************************/
/*!
 * \fn void dcsr_mxv_1 (dCSRmat *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x (index starts with 1!!)
 *
 * \param A   Pointer to dCSRmat matrix A
 * \param x   Pointer to array x
 * \param y   Pointer to array y
 *
 * \note Index starts at 1!!! -- Xiaozhe Hu
 *
 */
void dcsr_mxv_1(dCSRmat *A,
                REAL *x,
                REAL *y)
{
  // Shift A
  dcsr_shift(A,-1);

  // Perform Matrix Vector Multiply
  dcsr_mxv(A,x,y);

  // Shift A back
  dcsr_shift(A,1);
}

/***********************************************************************************************/
/*!
 * \fn dcsr_mxv_1_forts (void *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x (index starts with 1!!)
 *        Used for time-stepping routines
 *
 * \param A   Pointer to dCSRmat matrix A
 * \param x   Pointer to array x
 * \param y   Pointer to array y
 *
 * \note Index starts at 1!!! -- Xiaozhe Hu
 *
 */
void dcsr_mxv_1_forts(void *At,
                      REAL *x,
                      REAL *y)
{
  // Declare A as dCSRmat
  dCSRmat *A = (dCSRmat *) At;

  // Shift A
  dcsr_shift(A,-1);

  // Perform Matrix Vector Multiply
  dcsr_mxv(A,x,y);

  // Shift A back
  dcsr_shift(A,1);
}


/***********************************************************************************************/
/*!
 * \fn void dcsr_mxv_agg (dCSRmat *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x, where the entries of A are all ones.
 *
 * \param A   Pointer to dCSRmat matrix A
 * \param x   Pointer to array x
 * \param y   Pointer to array y
 *
 * \note This subroutine is used only for unsmoothed aggregation AMG!!! -- Xiaozhe Hu
 *
 */
void dcsr_mxv_agg(dCSRmat *A,
                  REAL *x,
                  REAL *y)
{
  const INT  m  = A->row;
  const INT *ia = A->IA, *ja = A->JA;
  INT i, k, begin_row, end_row;
  register REAL temp;

  for (i=0;i<m;++i) {
    temp=0.0;
    begin_row=ia[i]; end_row=ia[i+1];
    for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
    y[i]=temp;
  }

}

/***********************************************************************************************/
/*!
 * \fn void dcsr_aAxpy (const REAL alpha, dCSRmat *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y
 *
 * \param alpha  REAL factor alpha
 * \param A      Pointer to dCSRmat matrix A
 * \param x      Pointer to array x
 * \param y      Pointer to array y
 *
 * \note Index starts at 0!!! -- Xiaozhe Hu
 *
 */
void dcsr_aAxpy(const REAL alpha,
                dCSRmat *A,
                REAL *x,
                REAL *y)
{
  const INT  m  = A->row;
  const INT *ia = A->IA, *ja = A->JA;
  const REAL *aj = A->val;
  INT i, k, begin_row, end_row;
  register REAL temp;

  if ( alpha == 1.0 ) {
    for (i=0;i<m;++i) {
      temp=0.0;
      begin_row=ia[i]; end_row=ia[i+1];
      for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
      y[i]+=temp;
    }
  }

  else if ( alpha == -1.0 ) {
    for (i=0;i<m;++i) {
      temp=0.0;
      begin_row=ia[i]; end_row=ia[i+1];
      for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
      y[i]-=temp;
    }
  }

  else {
    for (i=0;i<m;++i) {
      temp=0.0;
      begin_row=ia[i]; end_row=ia[i+1];
      for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
      y[i]+=temp*alpha;
    }
  }
}

/***********************************************************************************************/

/***********************************************************************************************/
/*!
 * \fn void dcsr_aAxpy (const REAL alpha, dCSRmat *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y
 *
 * \param alpha  REAL factor alpha
 * \param A      Pointer to dCSRmat matrix A
 * \param x      Pointer to array x
 * \param y      Pointer to array y
 *
 * \note Index starts at 1!!! -- Xiaozhe Hu
 *
 */
void dcsr_aAxpy_1(const REAL alpha,
                  dCSRmat *A,
                  REAL *x,
                  REAL *y)
{
  // SHIFT A First
  dcsr_shift(A, -1);

  // multiply
  dcsr_aAxpy(alpha, A, x, y);

  // shift A back
  dcsr_shift(A, 1);
}

/***********************************************************************************************/
/*!
 * \fn void dcsr_aAxpy_agg (const REAL alpha, dCSRmat *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y (the entries of A are all ones)
 *
 * \param alpha  REAL factor alpha
 * \param A      Pointer to dCSRmat matrix A
 * \param x      Pointer to array x
 * \param y      Pointer to array y
 *
 * \note This subroutine is used only for unsmoothed aggregation AMG!!! -- Xiaozhe Hu
 *
 */
void dcsr_aAxpy_agg(const REAL alpha,
                    dCSRmat *A,
                    REAL *x,
                    REAL *y)
{
  const INT  m  = A->row;
  const INT *ia = A->IA, *ja = A->JA;

  INT i, k, begin_row, end_row;
  register REAL temp;

  if ( alpha == 1.0 ) {
    for (i=0;i<m;++i) {
      temp=0.0;
      begin_row=ia[i]; end_row=ia[i+1];
      for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
      y[i]+=temp;
    }
  }
  else if ( alpha == -1.0 ) {
    for (i=0;i<m;++i) {
      temp=0.0;
      begin_row=ia[i]; end_row=ia[i+1];
      for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
      y[i]-=temp;
    }
  }

  else {
    for (i=0;i<m;++i) {
      temp=0.0;
      begin_row=ia[i]; end_row=ia[i+1];
      for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
      y[i]+=temp*alpha;
    }
  }

}

/***********************************************************************************************/
/*!
 * \fn REAL dcsr_vmv (dCSRmat *A, REAL *x, REAL *y)
 *
 * \brief vector-Matrix-vector multiplication alpha = y'*A*x
 *
 * \param A   Pointer to dCSRmat matrix A
 * \param x   Pointer to array x
 * \param y   Pointer to array y
 *
 * \note Index starts at 0!!! -- Xiaozhe Hu
 *
 */
REAL dcsr_vmv(dCSRmat *A,
              REAL *x,
              REAL *y)
{
  register REAL value=0.0;
  const INT m=A->row;
  const INT *ia=A->IA, *ja=A->JA;
  const REAL *aj=A->val;
  INT i, k, begin_row, end_row;
  register REAL temp;

  for (i=0;i<m;++i) {
    temp=0.0;
    begin_row=ia[i]; end_row=ia[i+1];
    for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
    value+=y[i]*temp;
  }

  return value;
}

/***********************************************************************************************/
/*!
 * \fn REAL dcsr_vmv_1 (dCSRmat *A, REAL *x, REAL *y)
 *
 * \brief vector-Matrix-vector multiplication alpha = y'*A*x
 *
 * \param A   Pointer to dCSRmat matrix A
 * \param x   Pointer to array x
 * \param y   Pointer to array y
 *
 * \note Index starts at 1!!! -- Xiaozhe Hu
 *
 */
REAL dcsr_vmv_1(dCSRmat *A,
                REAL *x,
                REAL *y)
{
  REAL value=0.0;

  // SHIFT A First
  dcsr_shift(A, -1);

  // multiply
  value = dcsr_vmv(A, x, y);

  // shift A back
  dcsr_shift(A, 1);

  return value;
}


/***********************************************************************************************/
/*!
 * \fn void dcsr_mxm (dCSRmat *A, dCSRmat *B, dCSRmat *C)
 *
 * \brief Sparse matrix multiplication C=A*B (index starts with 0!!)
 *
 * \param A   Pointer to the dCSRmat matrix A
 * \param B   Pointer to the dCSRmat matrix B
 * \param C   Pointer to dCSRmat matrix equal to A*B
 *
 * \note Index starts at 0!!! -- Xiaozhe Hu
 *
 */
void dcsr_mxm(dCSRmat *A,
              dCSRmat *B,
              dCSRmat *C)
{
  INT i,j,k,l,count;

  INT *JD = (INT *)calloc(B->col,sizeof(INT));

  C->row=A->row;
  C->col=B->col;
  C->val = NULL;
  C->JA  = NULL;
  C->IA  = (INT*)calloc(C->row+1,sizeof(INT));

  for (i=0;i<B->col;++i) JD[i]=-1;

  // step 1: Find first the structure IA of C
  for(i=0;i<C->row;++i) {
    count=0;

    for (k=A->IA[i];k<A->IA[i+1];++k) {
      for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
        for (l=0;l<count;l++) {
          if (JD[l]==B->JA[j]) break;
        }

        if (l==count) {
          JD[count]=B->JA[j];
          count++;
        }
      }
    }
    C->IA[i+1]=count;
    for (j=0;j<count;++j) {
      JD[j]=-1;
    }
  }

  for (i=0;i<C->row;++i) C->IA[i+1]+=C->IA[i];

  // step 2: Find the structure JA of C
  INT countJD;

  C->JA=(INT*)calloc(C->IA[C->row],sizeof(INT));

  for (i=0;i<C->row;++i) {
    countJD=0;
    count=C->IA[i];
    for (k=A->IA[i];k<A->IA[i+1];++k) {
      for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
        for (l=0;l<countJD;l++) {
          if (JD[l]==B->JA[j]) break;
        }

        if (l==countJD) {
          C->JA[count]=B->JA[j];
          JD[countJD]=B->JA[j];
          count++;
          countJD++;
        }
      }
    }

    //for (j=0;j<countJD;++j) JD[j]=-1;
    memset(JD, -1, sizeof(INT)*countJD);
  }

  free(JD);

  // step 3: Find the structure A of C
  C->val=(REAL*)calloc(C->IA[C->row],sizeof(REAL));

  for (i=0;i<C->row;++i) {
    for (j=C->IA[i];j<C->IA[i+1];++j) {
      C->val[j]=0;
      for (k=A->IA[i];k<A->IA[i+1];++k) {
        for (l=B->IA[A->JA[k]];l<B->IA[A->JA[k]+1];l++) {
          if (B->JA[l]==C->JA[j]) {
            C->val[j]+=A->val[k]*B->val[l];
          } // end if
        } // end for l
      } // end for k
    } // end for j
  }    // end for i

  C->nnz = C->IA[C->row]-C->IA[0];

}

/***********************************************************************************************/
/*!
   * \fn void dcsr_mxm_1(dCSRmat *A, dCSRmat *B, dCSRmat *C)
   *
   * \brief Sparse matrix multiplication C=A*B (index starts with 1!!)
   *
   * \param A   Pointer to the dCSRmat matrix A
   * \param B   Pointer to the dCSRmat matrix B
   * \param C   Pointer to dCSRmat matrix equal to A*B
   *
   * \note Index starts at 1!!! -- Xiaozhe Hu
   *
   */
void dcsr_mxm_1(dCSRmat *A,
                dCSRmat *B,
                dCSRmat *C)
{
  // shift A and B
  dcsr_shift(A, -1);
  dcsr_shift(B, -1);

  // multiply A and B
  dcsr_mxm(A, B, C);

  // shift A, B, and C back
  dcsr_shift(A, 1);
  dcsr_shift(B, 1);
  dcsr_shift(C, 1);

}

/***********************************************************************************************/
/*!
   * \fn void icsr_mxm (iCSRmat *A, iCSRmat *B, iCSRmat *C)
   *
   * \brief Sparse matrix multiplication C=A*B (index starts with 0!!)
   *
   * \param A   Pointer to the iCSRmat matrix A
   * \param B   Pointer to the iCSRmat matrix B
   * \param C   Pointer to iCSRmat matrix equal to A*B
   *
   * \note Index starts at 0!!! -- Xiaozhe Hu
   *
   */
void icsr_mxm(iCSRmat *A,
              iCSRmat *B,
              iCSRmat *C)
{
  INT i,j,k,l,count;

  INT *JD = (INT *)calloc(B->col,sizeof(INT));

  C->row = A->row;
  C->col = B->col;
  C->val = NULL;
  C->JA  = NULL;
  C->IA  = (INT*)calloc(C->row+1,sizeof(INT));

  for (i=0;i<B->col;++i) JD[i]=-1;

  // step 1: Find first the structure IA of C
  for(i=0;i<C->row;++i) {
    count=0;

    for (k=A->IA[i];k<A->IA[i+1];++k) {
      for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
        for (l=0;l<count;l++) {
          if (JD[l]==B->JA[j]) break;
        }

        if (l==count) {
          JD[count]=B->JA[j];
          count++;
        }
      }
    }
    C->IA[i+1]=count;
    for (j=0;j<count;++j) {
      JD[j]=-1;
    }
  }

  for (i=0;i<C->row;++i) C->IA[i+1]+=C->IA[i];

  // step 2: Find the structure JA of C
  INT countJD;

  C->JA=(INT*)calloc(C->IA[C->row],sizeof(INT));

  for (i=0;i<C->row;++i) {
    countJD=0;
    count=C->IA[i];
    for (k=A->IA[i];k<A->IA[i+1];++k) {
      for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
        for (l=0;l<countJD;l++) {
          if (JD[l]==B->JA[j]) break;
        }

        if (l==countJD) {
          C->JA[count]=B->JA[j];
          JD[countJD]=B->JA[j];
          count++;
          countJD++;
        }
      }
    }

    //for (j=0;j<countJD;++j) JD[j]=-1;
    memset(JD, -1, sizeof(INT)*countJD);
  }

  free(JD);

  // step 3: Find the structure A of C
  C->val=(INT*)calloc(C->IA[C->row],sizeof(INT));

  for (i=0;i<C->row;++i) {
    for (j=C->IA[i];j<C->IA[i+1];++j) {
      C->val[j]=0;
      for (k=A->IA[i];k<A->IA[i+1];++k) {
        for (l=B->IA[A->JA[k]];l<B->IA[A->JA[k]+1];l++) {
          if (B->JA[l]==C->JA[j]) {
            C->val[j]+=A->val[k]*B->val[l];
          } // end if
        } // end for l
      } // end for k
    } // end for j
  }    // end for i

  C->nnz = C->IA[C->row]-C->IA[0];


}

/***********************************************************************************************/
/*!
   * \fn void icsr_mxm_1 (iCSRmat *A, iCSRmat *B, iCSRmat *C)
   *
   * \brief Sparse matrix multiplication C=A*B (index starts with 1!!)
   *
   * \param A   Pointer to the iCSRmat matrix A
   * \param B   Pointer to the iCSRmat matrix B
   * \param C   Pointer to iCSRmat matrix equal to A*B
   *
   * \note Index starts at 1!!! -- Xiaozhe Hu
   *
   */
void icsr_mxm_1(iCSRmat *A,
                iCSRmat *B,
                iCSRmat *C)
{
  // shift A and B
  icsr_shift(A, -1);
  icsr_shift(B, -1);

  icsr_mxm(A, B, C);

  // shift A, B, and C back
  icsr_shift(A, 1);
  icsr_shift(B, 1);
  icsr_shift(C, 1);

}

/***********************************************************************************************/
/*!
   * \fn void icsr_mxm_symb (iCSRmat *A, iCSRmat *B, iCSRmat *C)
   *
   * \brief Symbolic sparse matrix multiplication C=A*B (index starts with 0!!)
   *
   * \param A   Pointer to the iCSRmat matrix A
   * \param B   Pointer to the iCSRmat matrix B
   * \param C   Pointer to iCSRmat matrix equal to A*B
   *
   * \note Index starts at 0!!! -- Xiaozhe Hu
   *
   */
void icsr_mxm_symb(iCSRmat *A,
                   iCSRmat *B,
                   iCSRmat *C)
{
  INT i,j,k,l,count;

  INT *JD = (INT *)calloc(B->col,sizeof(INT));

  C->row = A->row;
  C->col = B->col;
  C->val = NULL;
  C->JA  = NULL;
  C->IA  = (INT*)calloc(C->row+1,sizeof(INT));

  for (i=0;i<B->col;++i) JD[i]=-1;

  // step 1: Find first the structure IA of C
  for(i=0;i<C->row;++i) {
    count=0;

    for (k=A->IA[i];k<A->IA[i+1];++k) {
      for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
        for (l=0;l<count;l++) {
          if (JD[l]==B->JA[j]) break;
        }

        if (l==count) {
          JD[count]=B->JA[j];
          count++;
        }
      }
    }
    C->IA[i+1]=count;
    for (j=0;j<count;++j) {
      JD[j]=-1;
    }
  }

  for (i=0;i<C->row;++i) C->IA[i+1]+=C->IA[i];

  // step 2: Find the structure JA of C
  INT countJD;

  C->JA=(INT*)calloc(C->IA[C->row],sizeof(INT));

  for (i=0;i<C->row;++i) {
    countJD=0;
    count=C->IA[i];
    for (k=A->IA[i];k<A->IA[i+1];++k) {
      for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
        for (l=0;l<countJD;l++) {
          if (JD[l]==B->JA[j]) break;
        }

        if (l==countJD) {
          C->JA[count]=B->JA[j];
          JD[countJD]=B->JA[j];
          count++;
          countJD++;
        }
      }
    }

    //for (j=0;j<countJD;++j) JD[j]=-1;
    memset(JD, -1, sizeof(INT)*countJD);
  }

  free(JD);

  C->nnz = C->IA[C->row]-C->IA[0];

}

/***********************************************************************************************/
/*!
   * \fn void icsr_mxm_symb_1 (iCSRmat *A, iCSRmat *B, iCSRmat *C)
   *
   * \brief Symbolic sparse matrix multiplication C=A*B (index starts with 1!!)
   *
   * \param A   Pointer to the iCSRmat matrix A
   * \param B   Pointer to the iCSRmat matrix B
   * \param C   Pointer to iCSRmat matrix equal to A*B
   *
   * \note Index starts at 0!!! -- Xiaozhe Hu
   *
   */
void icsr_mxm_symb_1(iCSRmat *A,
                     iCSRmat *B,
                     iCSRmat *C)
{
  // shift A and B
  icsr_shift(A, -1);
  icsr_shift(B, -1);

  // symbolic multiplication
  icsr_mxm_symb(A, B, C);

  // shift back
  icsr_shift(A, 1);
  icsr_shift(B, 1);
  icsr_shift(C, 1);

}

/***********************************************************************************************/
/*!
   * \fn void icsr_mxm_symb_max (iCSRmat *A, iCSRmat *B, iCSRmat *C, INT multmax)
   *
   * \brief symbolic sparse matrix multiplication C=A*B (index starts with 0!!)
   *
   * \param A         Pointer to the iCSRmat matrix A
   * \param B         Pointer to the iCSRmat matrix B
   * \param C         Pointer to iCSRmat matrix equal to A*B
   * \param multimax  value allowed in the iCSRmat matrix C, any entry that is not equal to multimax will be deleted
   *
   * \note Index starts at 0!!! -- Xiaozhe Hu
   *
   */
void icsr_mxm_symb_max(iCSRmat *A,
                       iCSRmat *B,
                       iCSRmat *C,
                       INT multmax)
{
  INT i,j,k,l,count;

  INT *JD = (INT *)calloc(B->col,sizeof(INT));
  INT *entry_count = (INT *)calloc(B->col,sizeof(INT));

  C->row = A->row;
  C->col = B->col;
  C->val = NULL;
  C->JA  = NULL;
  C->IA  = (INT*)calloc(C->row+1,sizeof(INT));

  for (i=0;i<B->col;++i) {
    JD[i] = -1;
    entry_count[i] = 0;
  }

  // step 1: Find first the structure IA of C
  for(i=0;i<C->row;++i) {
    count=0;

    for (k=A->IA[i];k<A->IA[i+1];++k) {
      for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
        for (l=0;l<count;l++) {
          if (JD[l]==B->JA[j]) {
            entry_count[l] = entry_count[l]+1;
            break;
          }
        }

        if (l==count) {
          JD[count]=B->JA[j];
          entry_count[count] = 1;
          count++;
        }
      }


    }


    C->IA[i+1]=count;

    for (j=0;j<count;++j) {

      JD[j]=-1;

      if (entry_count[j] != multmax) C->IA[i+1] = C->IA[i+1]-1;

      entry_count[j] = 0;

    }


  }

  for (i=0;i<C->row;++i) C->IA[i+1]+=C->IA[i];


  // step 2: Find the structure JA of C
  INT countJD;

  C->JA=(INT*)calloc(C->IA[C->row],sizeof(INT));

  for (i=0;i<C->row;++i) {
    countJD=0;
    count=C->IA[i];

    for (k=A->IA[i];k<A->IA[i+1];++k) {
      for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
        for (l=0;l<countJD;l++) {
          if (JD[l]==B->JA[j]) {
            entry_count[l] = entry_count[l]+1;
            break;
          }
        }

        if (l==countJD) {
          //C->JA[count]=B->JA[j];
          JD[countJD]=B->JA[j];
          entry_count[countJD] = 1;
          //count++;
          countJD++;
        }
      }

    }


    for (j=0;j<countJD;++j) {


      if (entry_count[j] == multmax) {
        C->JA[count]=JD[j];
        count++;
      }

      JD[j]=-1;
      entry_count[j] = 0;

    }



  }

  // free
  free(JD);
  free(entry_count);

  C->nnz = C->IA[C->row]-C->IA[0];
}

/***********************************************************************************************/
/*!
   * \fn void icsr_mxm_symb_max_1 (iCSRmat *A, iCSRmat *B, iCSRmat *C INT multmax)
   *
   * \brief symbolic sparse matrix multiplication C=A*B (index starts with 1!!)
   *
   * \param A         Pointer to the iCSRmat matrix A
   * \param B         Pointer to the iCSRmat matrix B
   * \param C         Pointer to iCSRmat matrix equal to A*B
   * \param multimax  value allowed in the iCSRmat matrix C, any entry that is not equal to multimax will be deleted
   *
   * \note Index starts at 1!!! -- Xiaozhe Hu
   *
   */
void icsr_mxm_symb_max_1(iCSRmat *A,
                         iCSRmat *B,
                         iCSRmat *C,
                         INT multmax)
{
  // shift A and B
  icsr_shift(A, -1);
  icsr_shift(B, -1);

  // multiply
  icsr_mxm_symb_max(A, B, C, multmax);

  // shift back
  icsr_shift(A, 1);
  icsr_shift(B, 1);
  icsr_shift(C, 1);

}

/***********************************************************************************************/
/*!
   * \fn void dcsr_getdiag (INT n, dCSRmat *A, dvector *diag)
   *
   * \brief Get first n diagonal entries of a dCSRmat sparse matrix A
   *
   * \param n     Number of diagonal entries to get (if n=0, then get all diagonal entries)
   * \param A     Pointer to dCSRmat CSR matrix
   * \param diag  Pointer to the diagonal as a dvector
   *
   * \note Index starts at 0!!! -- Xiaozhe Hu
   *
   */
void dcsr_getdiag(INT n,
                  dCSRmat *A,
                  dvector *diag)
{
  INT i,k,j,ibegin,iend;

  if ( n==0 || n>A->row || n>A->col ) n = MIN(A->row,A->col);

  dvec_alloc(n, diag);

  for (i=0;i<n;++i) {
    ibegin=A->IA[i]; iend=A->IA[i+1];
    for (k=ibegin;k<iend;++k) {
      j=A->JA[k];
      if ((j-i)==0) {
        diag->val[i] = A->val[k];
        break;
      } // end if
    } // end for k
  } // end for i

}

/***********************************************************************************************/
/*!
   * \fn void dcsr_diagpref ( dCSRmat *A )
   *
   * \brief Re-order the column and data arrays of a dCSRmat matrix row by row,
   *        so that the first entry in each row is the diagonal entry of that row
   *
   * \param A   Pointer to the matrix to be re-ordered
   *
   * \note Reordering is done in place.
   * \note Index starts at 0!!! -- Xiaozhe Hu
   *
   */
void dcsr_diagpref(dCSRmat *A)
{
  const INT    m    = A->row;
  INT         *IA   = A->IA;
  INT         *JA   = A->JA;
  REAL        *valA = A->val;


  // Local variable
  INT    i, j;
  INT    tempi, row_size;
  REAL   tempd;

  for (i=0; i<m; i ++) {

    row_size = IA[i+1] - IA[i];

    // check whether the first entry is already diagonal
    if (JA[0] != i) {
      for (j=1; j<row_size; j ++) {
        if (JA[j] == i) {
          //swap JA
          tempi  = JA[0];
          JA[0] = JA[j];
          JA[j] = tempi;

          // swap valA
          tempd     = valA[0];
          valA[0] = valA[j];
          valA[j] = tempd;

          break;
        }
      }
      if (j==row_size) {
        check_error(ERROR_DATA_ZERODIAG, __FUNCTION__);
      }
    }

    JA    += row_size;
    valA  += row_size;
  }

}

/***********************************************************************************************/
/*!
   * \fn void dcsr_rap(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *RAP)
   *
   * \brief Triple sparse matrix multiplication B=R*A*P
   *
   * \param R   Pointer to the dCSRmat matrix R
   * \param A   Pointer to the dCSRmat matrix A
   * \param P   Pointer to the dCSRmat matrix P
   * \param RAP Pointer to dCSRmat matrix equal to R*A*P
   *
   * \note Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package.
   *       Advances in Computational Mathematics, 1 (1993), pp. 127-137.
   * \note Index starts at 0!!! -- Xiaozhe Hu
   *
   */
void dcsr_rap(dCSRmat *R,
              dCSRmat *A,
              dCSRmat *P,
              dCSRmat *RAP)
{
  INT n_coarse = R->row;
  INT *R_i = R->IA;
  INT *R_j = R->JA;
  REAL *R_data = R->val;

  INT n_fine = A->row;
  INT *A_i = A->IA;
  INT *A_j = A->JA;
  REAL *A_data = A->val;

  INT *P_i = P->IA;
  INT *P_j = P->JA;
  REAL *P_data = P->val;

  INT  RAP_size;
  INT *RAP_i = NULL;
  INT *RAP_j = NULL;
  REAL *RAP_data = NULL;

  INT *Ps_marker = NULL;
  INT *As_marker = NULL;

  INT ic, i1, i2, i3, jj1, jj2, jj3;
  INT jj_counter, jj_row_begining;
  REAL r_entry, r_a_product, r_a_p_product;

  INT nthreads = 1;

  INT coarse_mul_nthreads = n_coarse * nthreads;
  INT fine_mul_nthreads = n_fine * nthreads;
  INT coarse_add_nthreads = n_coarse + nthreads;
  INT minus_one_length = coarse_mul_nthreads + fine_mul_nthreads;
  INT total_calloc = minus_one_length + coarse_add_nthreads + nthreads;

  Ps_marker = (INT *)calloc(total_calloc, sizeof(INT));
  As_marker = Ps_marker + coarse_mul_nthreads;

  /*------------------------------------------------------*
     *  First Pass: Determine size of RAP and set up RAP_i  *
     *------------------------------------------------------*/
  RAP_i = (INT *)calloc(n_coarse+1, sizeof(INT));

  iarray_set(minus_one_length, Ps_marker, -1);

  jj_counter = 0;
  for (ic = 0; ic < n_coarse; ic ++) {
    Ps_marker[ic] = jj_counter;
    jj_row_begining = jj_counter;
    jj_counter ++;

    for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++) {
      i1 = R_j[jj1];

      for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++) {
        i2 = A_j[jj2];
        if (As_marker[i2] != ic) {
          As_marker[i2] = ic;
          for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++) {
            i3 = P_j[jj3];
            if (Ps_marker[i3] < jj_row_begining) {
              Ps_marker[i3] = jj_counter;
              jj_counter ++;
            }
          }
        }
      }
    }

    RAP_i[ic] = jj_row_begining;
  }

  RAP_i[n_coarse] = jj_counter;
  RAP_size = jj_counter;


  RAP_j = (INT *)calloc(RAP_size, sizeof(INT));
  RAP_data = (REAL *)calloc(RAP_size, sizeof(REAL));

  iarray_set(minus_one_length, Ps_marker, -1);

  jj_counter = 0;
  for (ic = 0; ic < n_coarse; ic ++) {
    Ps_marker[ic] = jj_counter;
    jj_row_begining = jj_counter;
    RAP_j[jj_counter] = ic;
    RAP_data[jj_counter] = 0.0;
    jj_counter ++;

    for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++) {
      r_entry = R_data[jj1];

      i1 = R_j[jj1];
      for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++) {
        r_a_product = r_entry * A_data[jj2];

        i2 = A_j[jj2];
        if (As_marker[i2] != ic) {
          As_marker[i2] = ic;
          for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++) {
            r_a_p_product = r_a_product * P_data[jj3];

            i3 = P_j[jj3];
            if (Ps_marker[i3] < jj_row_begining) {
              Ps_marker[i3] = jj_counter;
              RAP_data[jj_counter] = r_a_p_product;
              RAP_j[jj_counter] = i3;
              jj_counter ++;
            }
            else {
              RAP_data[Ps_marker[i3]] += r_a_p_product;
            }
          }
        }
        else {
          for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++) {
            i3 = P_j[jj3];
            r_a_p_product = r_a_product * P_data[jj3];
            RAP_data[Ps_marker[i3]] += r_a_p_product;
          }
        }
      }
    }
  }

  RAP->row = n_coarse;
  RAP->col = n_coarse;
  RAP->nnz = RAP_size;
  RAP->IA = RAP_i;
  RAP->JA = RAP_j;
  RAP->val = RAP_data;

  free(Ps_marker);
}

/***********************************************************************************************/
/*!
   * \fn void dcsr_rap_agg (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *RAP)
   *
   * \brief Triple sparse matrix multiplication B=R*A*P (all the entries in R and P are ones)
   *
   * \param R   Pointer to the dCSRmat matrix R
   * \param A   Pointer to the dCSRmat matrix A
   * \param P   Pointer to the dCSRmat matrix P
   * \param RAP Pointer to dCSRmat matrix equal to R*A*P
   *
   * \note Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package.
   *       Advances in Computational Mathematics, 1 (1993), pp. 127-137.
   * \note Index starts at 0!!! -- Xiaozhe Hu
   * \note Only used in unsmoothed aggregation AMG!!! -- Xiaozhe Hu
   *
   */
void dcsr_rap_agg(dCSRmat *R,
                  dCSRmat *A,
                  dCSRmat *P,
                  dCSRmat *RAP)
{
  INT n_coarse = R->row;
  INT *R_i = R->IA;
  INT *R_j = R->JA;

  INT n_fine = A->row;
  INT *A_i = A->IA;
  INT *A_j = A->JA;
  REAL *A_data = A->val;

  INT *P_i = P->IA;
  INT *P_j = P->JA;

  INT RAP_size;
  INT *RAP_i = NULL;
  INT *RAP_j = NULL;
  REAL *RAP_data = NULL;

  INT *Ps_marker = NULL;
  INT *As_marker = NULL;

  INT ic, i1, i2, i3, jj1, jj2, jj3;
  INT jj_counter, jj_row_begining;

  INT nthreads = 1;

  INT coarse_mul_nthreads = n_coarse * nthreads;
  INT fine_mul_nthreads = n_fine * nthreads;
  INT coarse_add_nthreads = n_coarse + nthreads;
  INT minus_one_length = coarse_mul_nthreads + fine_mul_nthreads;
  INT total_calloc = minus_one_length + coarse_add_nthreads + nthreads;

  Ps_marker = (INT *)calloc(total_calloc, sizeof(INT));
  As_marker = Ps_marker + coarse_mul_nthreads;

  /*------------------------------------------------------*
     *  First Pass: Determine size of RAP and set up RAP_i  *
     *------------------------------------------------------*/
  RAP_i = (INT *)calloc(n_coarse+1, sizeof(INT));

  iarray_set(minus_one_length, Ps_marker, -1);

  jj_counter = 0;
  for (ic = 0; ic < n_coarse; ic ++) {
    Ps_marker[ic] = jj_counter;
    jj_row_begining = jj_counter;
    jj_counter ++;

    for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++) {
      i1 = R_j[jj1];

      for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++) {
        i2 = A_j[jj2];
        if (As_marker[i2] != ic) {
          As_marker[i2] = ic;
          for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++) {
            i3 = P_j[jj3];
            if (Ps_marker[i3] < jj_row_begining) {
              Ps_marker[i3] = jj_counter;
              jj_counter ++;
            }
          }
        }
      }
    }

    RAP_i[ic] = jj_row_begining;
  }

  RAP_i[n_coarse] = jj_counter;
  RAP_size = jj_counter;

  RAP_j = (INT *)calloc(RAP_size, sizeof(INT));
  RAP_data = (REAL *)calloc(RAP_size, sizeof(REAL));

  iarray_set(minus_one_length, Ps_marker, -1);

  jj_counter = 0;
  for (ic = 0; ic < n_coarse; ic ++) {
    Ps_marker[ic] = jj_counter;
    jj_row_begining = jj_counter;
    RAP_j[jj_counter] = ic;
    RAP_data[jj_counter] = 0.0;
    jj_counter ++;

    for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++) {
      i1 = R_j[jj1];
      for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++) {
        i2 = A_j[jj2];
        if (As_marker[i2] != ic) {
          As_marker[i2] = ic;
          for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++) {
            i3 = P_j[jj3];
            if (Ps_marker[i3] < jj_row_begining) {
              Ps_marker[i3] = jj_counter;
              RAP_data[jj_counter] = A_data[jj2];
              RAP_j[jj_counter] = i3;
              jj_counter ++;
            }
            else {
              RAP_data[Ps_marker[i3]] += A_data[jj2];
            }
          }
        }
        else {
          for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++) {
            i3 = P_j[jj3];
            RAP_data[Ps_marker[i3]] += A_data[jj2];
          }
        }
      }
    }
  }

  RAP->row = n_coarse;
  RAP->col = n_coarse;
  RAP->nnz = RAP_size;
  RAP->IA = RAP_i;
  RAP->JA = RAP_j;
  RAP->val = RAP_data;

  free(Ps_marker);
}

/***********************************************************************************************/
/*!
   * \fn SHORT dcsr_getblk (dCSRmat *A, INT *Is, INT *Js, const INT m,
   *                             const INT n, dCSRmat *B)
   *
   * \brief Get a sub dCSRmat sparse matrix from a dCSRmat sparse matrix A with specified rows and columns
   *
   * \param A     Pointer to dCSRmat matrix
   * \param B     Pointer to dCSRmat matrix
   * \param Is    Pointer to selected rows
   * \param Js    Pointer to selected columns
   * \param m     Number of selected rows
   * \param n     Number of selected columns
   *
   * \return      SUCCESS if succeeded, otherwise return error information.
   *
   * \note Index starts at 0!!! -- Xiaozhe Hu
   *
   */
SHORT dcsr_getblk(dCSRmat *A,
                  INT *Is,
                  INT *Js,
                  const INT m,
                  const INT n,
                  dCSRmat *B)
{
  INT status = SUCCESS;

  INT i,j,k,nnz=0;
  INT *col_flag;

  // create column flags
  col_flag = (INT*)calloc(A->col,sizeof(INT));

  B->row = m; B->col = n;

  B->IA  = (INT*)calloc(m+1,sizeof(INT));
  B->JA  = (INT*)calloc(A->nnz,sizeof(INT));
  B->val = (REAL*)calloc(A->nnz,sizeof(REAL));

  for (i=0;i<n;++i) col_flag[Js[i]]=i+1;

  // Count nonzeros for sub matrix and fill in
  B->IA[0]=0;
  for (i=0;i<m;++i) {
    for (k=A->IA[Is[i]];k<A->IA[Is[i]+1];++k) {
      j=A->JA[k];
      if (col_flag[j]>0) {
        B->JA[nnz]=col_flag[j]-1;
        B->val[nnz]=A->val[k];
        nnz++;
      }
    } /* end for k */
    B->IA[i+1]=nnz;
  } /* end for i */
  B->nnz=nnz;


  // re-allocate memory space
  B->JA=(INT*)realloc(B->JA, sizeof(INT)*nnz);
  B->val=(REAL*)realloc(B->val, sizeof(REAL)*nnz);

  free(col_flag);

  return(status);
}

/***********************************************************************************************/
/*!
   * \fn dcsr_bandwith (dCSRmat *A, INT *bndwith)
   *
   * \brief Get bandwith of a dCSRmat format sparse matrix
   *
   * \param A       pointer to the dCSRmat matrix
   * \param bndwith pointer to the bandwith
   *
   */
void dcsr_bandwith(dCSRmat *A,
                   INT *bndwith)
{
  INT row = A->row;
  INT *IA = A->IA;

  INT i, max;
  max = 0;

  for (i=0; i<row; ++i) max = MAX(max, IA[i+1]-IA[i]);

  *bndwith = max;
}


/***********************************************************************************************/
// Block_dCSRmat subroutines starts here!
/***********************************************************************************************/
/*!
   * \fn void bdcsr_alloc (const INT brow, const INT bcol, block_dCSRmat *A)
   *
   * \brief Allocate block dCSRmat sparse matrix memory space
   *
   * \param brow      Number of block rows
   * \param bcol      Number of block columns
   * \param A         Pointer to the dCSRmat matrix
   *
   * \note: this only allocates A->blocks, but does not allocate each block
   *
   */
void bdcsr_alloc(const INT brow,
                 const INT bcol,
                 block_dCSRmat *A)
{
    SHORT i;

    A->brow = brow;
    A->bcol = bcol;

    if (brow == 0 || bcol == 0){
        A->blocks = NULL;
    }
    else {
        A->blocks = (dCSRmat **) calloc(brow*bcol,sizeof(dCSRmat *));
        for (i=0; i<brow*bcol; i++) A->blocks[i] = (dCSRmat *)calloc(1,sizeof(dCSRmat));
    }

  return;
}

/***********************************************************************************************/
/*!
   * \fn void bdcsr_free (block_dCSRmat *A)
   *
   * \brief Free block dCSR sparse matrix data
   *
   * \param A   Pointer to the block_dCSRmat matrix
   *
   */
void bdcsr_free(block_dCSRmat *A)
{
  if (A == NULL) return; // Nothing need to be freed!

  INT i;
  INT num_blocks = (A->brow)*(A->bcol);

  for ( i=0; i<num_blocks; i++ ) {
    dcsr_free(A->blocks[i]);
    A->blocks[i] = NULL;
  }

  free(A->blocks);
  A->blocks = NULL;

  return;
}

/***********************************************************************************************/
/*!
   * \fn void bdcsr_dealloc (block_dCSRmat *A)
   *
   * \brief Free block dCSR sparse matrix data, which is created from bdcsr_alloc
   * \note This involves allocating space for each block, since size isn't known ahead of time
   *
   * \param A   Pointer to the block_dCSRmat matrix
   *
   */
void bdcsr_dealloc(block_dCSRmat *A)
{
  if (A == NULL) return; // Nothing need to be freed!

  INT i;
  INT num_blocks = (A->brow)*(A->bcol);

  for ( i=0; i<num_blocks; i++ ) {
    dcsr_free(A->blocks[i]);
    if(A->blocks[i]) free(A->blocks[i]);
    A->blocks[i] = NULL;
  }

  free(A->blocks);
  A->blocks = NULL;

  return;
}
/***********************************************************************************************/
/*!
   * \fn void bdcsr_cp (block_dCSRmat *A, block_dCSRmat *B)
   *
   * \brief copy a block_dCSRmat to a new one B=A
   *
   * \param A   Pointer to the block_dCSRmat matrix
   * \param B   Pointer to the block_dCSRmat matrix
   *
   * \note: B->blocks has been allocated outsie but each block will be allocated here
   *
   */
void bdcsr_cp(block_dCSRmat *A,
              block_dCSRmat *B)
{
    SHORT i;

    for (i=0; i<(A->brow*A->bcol); i++){

        if (A->blocks[i] == NULL){
            B->blocks[i] = NULL;
        } else {
            dcsr_alloc(A->blocks[i]->row, A->blocks[i]->col, A->blocks[i]->nnz, B->blocks[i]);
            dcsr_cp(A->blocks[i], B->blocks[i]);
        }

    }

}


/***********************************************************************************************/
/*!
   * \fn void bdcsr_add (block_dCSRmat *A, const REAL alpha, block_dCSRmat *B,
   *                              const REAL beta, block_dCSRmat *C)
   *
   * \brief compute C = alpha*A + beta*B in block_dCSRmat format
   *
   * \param A      Pointer to block dCSRmat matrix
   * \param alpha  REAL factor alpha
   * \param B      Pointer to block_dCSRmat matrix
   * \param beta   REAL factor beta
   * \param C      Pointer to block_dCSRmat matrix
   *
   * \return Flag of whether the adding is succesful or not (SUCCUESS: 0; FAIL: <0)
   *
   */
INT bdcsr_add(block_dCSRmat *A,
              const REAL alpha,
              block_dCSRmat *B,
              const REAL beta,
              block_dCSRmat *C)
{
    INT i,j;
    INT status = SUCCESS;

    // both matrices A and B are NULL
    if (A == NULL && B == NULL) {
      C->brow=0; C->bcol=0; C->blocks=NULL;
      status=SUCCESS;
      goto FINISHED;
    }

    // only matrices A is NULL
    if (A == NULL) {
        for (i=0; i<B->brow; i++){
            for (j=0; j<B->bcol; j++){
                status = dcsr_add_1(NULL, alpha, B->blocks[i*A->brow+j], beta, C->blocks[i*A->brow+j]);
                if (status < 0) {goto FINISHED;}
            }
        }
    }

    // only matrices B is NULL
    if (B == NULL) {
        for (i=0; i<A->brow; i++){
            for (j=0; j<A->bcol; j++){
                status = dcsr_add_1(A->blocks[i*A->brow+j], alpha, NULL, beta, C->blocks[i*A->brow+j]);
                if (status < 0) {goto FINISHED;}
            }
        }
    }

    if (A->brow != B->brow || A->bcol != B->bcol) {
      printf("### ERROR HAZMATH DANGER: Dimensions of block matrices do not match!!! %s\n", __FUNCTION__);
      status = ERROR_MAT_SIZE;
      goto FINISHED;
    }

    // blockwise addition
    for (i=0; i<A->brow; i++){
        for (j=0; j<A->bcol; j++){
            status = dcsr_add(A->blocks[i*A->brow+j], alpha, B->blocks[i*A->brow+j], beta, C->blocks[i*A->brow+j]);
            if (status < 0) {goto FINISHED;}
        }
    }

FINISHED:
  return status;

}

/***********************************************************************************************/
/*!
   * \fn void bdcsr_add_1 (block_dCSRmat *A, const REAL alpha, block_dCSRmat *B,
   *                              const REAL beta, block_dCSRmat *C)
   *
   * \brief compute C = alpha*A + beta*B in block_dCSRmat format
   *
   * \param A      Pointer to block dCSRmat matrix
   * \param alpha  REAL factor alpha
   * \param B      Pointer to block_dCSRmat matrix
   * \param beta   REAL factor beta
   * \param C      Pointer to block_dCSRmat matrix
   *
   * \return Flag of whether the adding is succesful or not (SUCCUESS: 0; FAIL: <0)
   *
   */
INT bdcsr_add_1(block_dCSRmat *A,
                const REAL alpha,
                block_dCSRmat *B,
                const REAL beta,
                block_dCSRmat *C)
{
    INT i,j;
    INT status = SUCCESS;

    // both matrices A and B are NULL
    if (A == NULL && B == NULL) {
      C->brow=0; C->bcol=0; C->blocks=NULL;
      status=SUCCESS;
      goto FINISHED;
    }

    // only matrices A is NULL
    if (A == NULL) {
        for (i=0; i<B->brow; i++){
            for (j=0; j<B->bcol; j++){
                status = dcsr_add_1(NULL, alpha, B->blocks[i*B->brow+j], beta, C->blocks[i*B->brow+j]);
                if (status < 0) {goto FINISHED;}
            }
        }
    }

    // only matrices B is NULL
    if (B == NULL) {
        for (i=0; i<A->brow; i++){
            for (j=0; j<A->bcol; j++){
                status = dcsr_add_1(A->blocks[i*A->brow+j], alpha, NULL, beta, C->blocks[i*A->brow+j]);
                if (status < 0) {goto FINISHED;}
            }
        }
    }

    // none of the matrices are NULL
    if (A->brow != B->brow || A->bcol != B->bcol) {
      printf("### ERROR HAZMATH DANGER: Dimensions of block matrices do not match!!! %s\n", __FUNCTION__);
      status = ERROR_MAT_SIZE;
      goto FINISHED;
    }

    // blockwise addition
    for (i=0; i<A->brow; i++){
        for (j=0; j<A->bcol; j++){
            status = dcsr_add_1(A->blocks[i*A->brow+j], alpha, B->blocks[i*A->brow+j], beta, C->blocks[i*A->brow+j]);
            if (status < 0) {goto FINISHED;}
        }
    }

FINISHED:
  return status;

}

/***********************************************************************************************/
/*!
   * \fn void bdcsr_aAxpy (const REAL alpha, block_dCSRmat *A,
   *                                 REAL *x, REAL *y)
   *
   * \brief Matrix-vector multiplication y = alpha*A*x + y in block_dCSRmat format
   *
   * \param alpha  REAL factor a
   * \param A      Pointer to block_dCSRmat matrix A
   * \param x      Pointer to array x
   * \param y      Pointer to array y
   *
   */
void bdcsr_aAxpy(const REAL alpha,
                 block_dCSRmat *A,
                 REAL *x,
                 REAL *y)
{
  // information of A
  INT brow = A->brow;

  // local variables
  register dCSRmat *A11, *A12, *A21, *A22;
  register dCSRmat *A13, *A23, *A31, *A32, *A33;

  unsigned INT row1, col1;
  unsigned INT row2, col2;

  register REAL *x1, *x2, *y1, *y2;
  register REAL *x3, *y3;

  INT i,j;
  INT start_row;
  INT start_col;

  switch (brow) {

  case 2:
    A11 = A->blocks[0];
    A12 = A->blocks[1];
    A21 = A->blocks[2];
    A22 = A->blocks[3];

    row1 = A11->row;
    col1 = A11->col;

    x1 = x;
    x2 = &(x[col1]);
    y1 = y;
    y2 = &(y[row1]);

    // y1 = alpha*A11*x1 + alpha*A12*x2 + y1
    if (A11) dcsr_aAxpy(alpha, A11, x1, y1);
    if (A12) dcsr_aAxpy(alpha, A12, x2, y1);

    // y2 = alpha*A21*x1 + alpha*A22*x2 + y2
    if (A21) dcsr_aAxpy(alpha, A21, x1, y2);
    if (A22) dcsr_aAxpy(alpha, A22, x2, y2);

    break;

  case 3:
    A11 = A->blocks[0];
    A12 = A->blocks[1];
    A13 = A->blocks[2];
    A21 = A->blocks[3];
    A22 = A->blocks[4];
    A23 = A->blocks[5];
    A31 = A->blocks[6];
    A32 = A->blocks[7];
    A33 = A->blocks[8];

    row1 = A11->row;
    col1 = A11->col;
    row2 = A22->row;
    col2 = A22->col;

    x1 = x;
    x2 = &(x[col1]);
    x3 = &(x[col1+col2]);
    y1 = y;
    y2 = &(y[row1]);
    y3 = &(y[row1+row2]);

    // y1 = alpha*A11*x1 + alpha*A12*x2 + alpha*A13*x3 + y1
    if (A11) dcsr_aAxpy(alpha, A11, x1, y1);
    if (A12) dcsr_aAxpy(alpha, A12, x2, y1);
    if (A13) dcsr_aAxpy(alpha, A13, x3, y1);

    // y2 = alpha*A21*x1 + alpha*A22*x2 + alpha*A23*x3 + y2
    if (A21) dcsr_aAxpy(alpha, A21, x1, y2);
    if (A22) dcsr_aAxpy(alpha, A22, x2, y2);
    if (A23) dcsr_aAxpy(alpha, A23, x3, y2);

    // y3 = alpha*A31*x1 + alpha*A32*x2 + alpha*A33*x3 + y2
    if (A31) dcsr_aAxpy(alpha, A31, x1, y3);
    if (A32) dcsr_aAxpy(alpha, A32, x2, y3);
    if (A33) dcsr_aAxpy(alpha, A33, x3, y3);

    break;

  default:

    start_row = 0;
    start_col = 0;

    for (i=0; i<brow; i++) {

      for (j=0; j<brow; j++){

        if (A->blocks[i*brow+j]){
          dcsr_aAxpy(alpha, A->blocks[i*brow+j], &(x[start_col]), &(y[start_row]));
        }
        start_col = start_col + A->blocks[j*brow+j]->col;
      }

      start_row = start_row + A->blocks[i*brow+i]->row;
      start_col = 0;
    }

    break;

  } // end of switch

}

/***********************************************************************************************/
/*!
   * \fn void bdcsr_mxv (block_dCSRmat *A, REAL *x, REAL *y)
   *
   * \brief Matrix-vector multiplication y = A*x in block_dCSRmat format
   *
   * \param A      Pointer to block_dCSRmat matrix A
   * \param x      Pointer to array x
   * \param y      Pointer to array y
   *
   */
void bdcsr_mxv(block_dCSRmat *A,
               REAL *x,
               REAL *y)
{
  // information of A
  INT brow = A->brow;

  // local variables
  register dCSRmat *A11, *A12, *A21, *A22;
  register dCSRmat *A13, *A23, *A31, *A32, *A33;

  unsigned INT row1, col1;
  unsigned INT row2, col2;

  register REAL *x1, *x2, *y1, *y2;
  register REAL *x3, *y3;

  INT i,j,k;
  INT start_row;
  INT start_col;

  switch (brow) {

  case 2:
    A11 = A->blocks[0];
    A12 = A->blocks[1];
    A21 = A->blocks[2];
    A22 = A->blocks[3];

    row1 = A11->row;
    col1 = A11->col;

    x1 = x;
    x2 = &(x[col1]);
    y1 = y;
    y2 = &(y[row1]);

    // y1 = A11*x1 + A12*x2
    if (A11) dcsr_mxv(A11, x1, y1);
    if (A12) dcsr_aAxpy(1.0, A12, x2, y1);

    // y2 = A21*x1 + A22*x2
    if (A21) dcsr_mxv(A21, x1, y2);
    if (A22) dcsr_aAxpy(1.0, A22, x2, y2);

    break;

  case 3:
    A11 = A->blocks[0];
    A12 = A->blocks[1];
    A13 = A->blocks[2];
    A21 = A->blocks[3];
    A22 = A->blocks[4];
    A23 = A->blocks[5];
    A31 = A->blocks[6];
    A32 = A->blocks[7];
    A33 = A->blocks[8];

    row1 = A11->row;
    col1 = A11->col;
    row2 = A22->row;
    col2 = A22->col;

    x1 = x;
    x2 = &(x[col1]);
    x3 = &(x[col1+col2]);
    y1 = y;
    y2 = &(y[row1]);
    y3 = &(y[row1+row2]);

    // y1 = A11*x1 + A12*x2 + A13*x3 + y1
    if (A11) dcsr_mxv(A11, x1, y1);
    if (A12) dcsr_aAxpy(1.0, A12, x2, y1);
    if (A13) dcsr_aAxpy(1.0, A13, x3, y1);

    // y2 = A21*x1 + A22*x2 + A23*x3 + y2
    if (A21) dcsr_mxv(A21, x1, y2);
    if (A22) dcsr_aAxpy(1.0, A22, x2, y2);
    if (A23) dcsr_aAxpy(1.0, A23, x3, y2);

    // y3 = A31*x1 + A32*x2 + A33*x3 + y2
    if (A31) dcsr_mxv(A31, x1, y3);
    if (A32) dcsr_aAxpy(1.0, A32, x2, y3);
    if (A33) dcsr_aAxpy(1.0, A33, x3, y3);

    break;

  default:

    start_row = 0;
    start_col = 0;

    for (i=0; i<brow; i++) {

      for (j=0; j<brow; j++){

        if (j==0) {
          if (A->blocks[i*brow+j]){
            dcsr_mxv(A->blocks[i*brow+j], &(x[start_col]), &(y[start_row]));
          }
        }
        else {
          if (A->blocks[i*brow+j]){
            dcsr_aAxpy(1.0, A->blocks[i*brow+j], &(x[start_col]), &(y[start_row]));
          }
        }

        for (k=0; k<brow; k++){
           if(A->blocks[k*brow+j])
           {
               start_col = start_col + A->blocks[k*brow+j]->col;
               break;
           }
        }
      }

      for (k=0; k<brow; k++){
          if(A->blocks[i*brow+k])
          {
            start_row = start_row + A->blocks[i*brow+k]->row;
            break;
          }
      }

      start_col = 0;

    }

    break;

  } // end of switch

}

/***********************************************************************************************/

/***********************************************************************************************/
/*!
   * \fn void bdcsr_mxv_1 (block_dCSRmat *A, REAL *x, REAL *y)
   *
   * \brief Matrix-vector multiplication y = A*x in block_dCSRmat format
   * \note Assume each dscr blocks indexes at 1.
   *
   * \param A      Pointer to block_dCSRmat matrix A
   * \param x      Pointer to array x
   * \param y      Pointer to array y
   *
   */
void bdcsr_mxv_1(block_dCSRmat *A,
               REAL *x,
               REAL *y)
{
  // information of A
  INT brow = A->brow;

  INT i,j,k;
  INT start_row = 0;
  INT start_col = 0;


  for (i=0; i<brow; i++) {

    for (j=0; j<brow; j++){

      if (j==0) {
        if (A->blocks[i*brow+j]){
          dcsr_mxv_1(A->blocks[i*brow+j], &(x[start_col]), &(y[start_row]));
        }
      }
      else {
        if (A->blocks[i*brow+j]){
          dcsr_aAxpy_1(1.0, A->blocks[i*brow+j], &(x[start_col]), &(y[start_row]));
        }
      }

      for (k=0; k<brow; k++){
         if(A->blocks[k*brow+j])
         {
             start_col = start_col + A->blocks[k*brow+j]->col;
             break;
         }
      }

    }

    for (k=0; k<brow; k++){
        if(A->blocks[i*brow+k])
        {
          start_row = start_row + A->blocks[i*brow+k]->row;
          break;
        }
    }

    start_col = 0;
  }

}

/***********************************************************************************************/

/***********************************************************************************************/
/*!
   * \fn void bdcsr_shift (block_dCSRmat *A, INT shift)
   *
   * \brief Shift indexing in block_dCSRmat format
   *
   * \param A      Pointer to block_dCSRmat matrix A
   *
   */
void bdcsr_shift(block_dCSRmat *A,INT shift)
{
  // information of A
  INT brow = A->brow;

  INT i,j;


  for (i=0; i<brow; i++) {

    for (j=0; j<brow; j++){

      if (A->blocks[i*brow+j]) {
        dcsr_shift(A->blocks[i*brow+j],shift);
      }

    }

  }
}

/***********************************************************************************************/

/*!
   * \fn void bdcsr_mxv_forts (block_dCSRmat *A, REAL *x, REAL *y)
   *
   * \brief Matrix-vector multiplication y = A*x in block_dCSRmat format
   *
   * \param A      Pointer to block_dCSRmat matrix A
   * \param x      Pointer to array x
   * \param y      Pointer to array y
   *
   */
void bdcsr_mxv_forts(void *At,
                     REAL *x,
                     REAL *y)
{
  // Declare A as block_dCSRmat
  block_dCSRmat *A = (block_dCSRmat *) At;

  // information of A
  INT brow = A->brow;

  INT i,j,k;
  INT start_row = 0;
  INT start_col = 0;

  for (i=0; i<brow; i++) {

    for (j=0; j<brow; j++){

      if (j==0) {
        if (A->blocks[i*brow+j]){
          dcsr_mxv_1(A->blocks[i*brow+j], &(x[start_col]), &(y[start_row]));
        }
      }
      else {
        if (A->blocks[i*brow+j]){
          dcsr_aAxpy_1(1.0, A->blocks[i*brow+j], &(x[start_col]), &(y[start_row]));
        }
      }

      for (k=0; k<brow; k++){
         if(A->blocks[k*brow+j])
         {
             start_col = start_col + A->blocks[k*brow+j]->col;
             break;
         }
      }

    }

    for (k=0; k<brow; k++){
        if(A->blocks[i*brow+k])
        {
          start_row = start_row + A->blocks[i*brow+k]->row;
          break;
        }
    }

    start_col = 0;
  }


}

/************************************** END ***************************************************/

