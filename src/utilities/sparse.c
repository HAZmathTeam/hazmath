/*! \file src/utilities/sparse.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 3/6/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/30/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/31/2016 & 08/28/2021
 *  \note: modified by James Adler on 02/22/2019 for 0-1 fix
 *  \note: modified by ludmil zikatanov on 20200412
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
 * \param  index_start   Number from which memory is indexed (1 for fortran, 0 for C)
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
/*!
 * \fn dCSRmat dcsr_create_fullmatrix (const INT m, const INT n)
 *
 * \brief Create a CSR sparse matrix that is actually full
 *
 * \param  m             Number of rows
 * \param  n             Number of columns
 *
 * \return A             the new dCSRmat matrix with zero entries
 *
 */
dCSRmat dcsr_create_fullmatrix(const INT m,
                               const INT n)
{
  dCSRmat A;
  INT nnz = m*n;

  A.IA = (INT *)calloc(m+1, sizeof(INT));
  A.JA = (INT *)calloc(nnz, sizeof(INT));
  A.val = (REAL *)calloc(nnz, sizeof(REAL));

  A.row = m; A.col = n; A.nnz = nnz;

  INT i,j;
  for(i=0;i<m+1;i++) A.IA[i]=i*n;
  for(i=0;i<m;i++) {
    for(j=0;j<n;j++) {
      A.JA[i*n+j] = j;
    }
  }

  return A;
}

/***********************************************************************************************/
/*!
 * \fn void dcsr_set_zeromatrix(dCSRmat *A,const INT m,const INT n,const INT index_start)
{
 *
 * \brief Create a CSR sparse matrix that is all zeros
 *
 * \param  A             pointer to the dCSRmat
 * \param  m             Number of rows
 * \param  n             Number of columns
 * \param  index_start   Number from which memory is indexed (1 for fortran, 0 for C)
 *
 *
 */
void dcsr_set_zeromatrix(dCSRmat *A,
                         const INT m,
                         const INT n,
                         const INT index_start)
{

  A->IA = (INT *)calloc(m+1, sizeof(INT));
  A->JA = (INT *)calloc(1, sizeof(INT));
  A->val = (REAL *)calloc(1, sizeof(REAL));

  A->row = m; A->col = n; A->nnz = 1;

  INT i;
  A->IA[0]=index_start;
  for(i=1;i<m+1;i++) A->IA[i]=index_start+1;
  A->JA[0]=index_start;
  A->val[0] = 0.0;

  return;
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
 * \param index_start   Number from which memory is indexed (1 for fortran, 0 for C)
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
 * \param index_start   Number from which memory is indexed (1 for fortran, 0 for C)
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

/***********************************************************************/
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
/***********************************************************************/
/*!
 * \fn void dcsr_realloc (const INT m, const INT n, const INT nnz, dCSRmat *A)
 *
 * \brief RE-allocate dCSRmat sparse matrix. Just use realloc to
 *        expand, shrink the arrays in A
 *
 * \param m      Number of rows
 * \param n      Number of columns
 * \param nnz    Number of nonzeros
 * \param A      Pointer to the dCSRmat matrix
 *
 */
void dcsr_realloc(const INT m,
		  const INT n,
		  const INT nnz,
		  dCSRmat *A)
{
  if ( m > 0 ) {
    A->IA=(INT*)realloc(A->IA,(m+1)*sizeof(INT));
  } else {
    if(A->IA) {
      free(A->IA);
      A->IA = NULL;
    }
  }
  //
  if ( n > 0 ) {
    A->JA=(INT*)realloc(A->JA,nnz*sizeof(INT));
  } else {
    if(A->JA){
      free(A->JA);
      A->JA = NULL;
    }
  }
  //
  if ( nnz > 0 ) {
    A->val=(REAL*)realloc(A->val,nnz*sizeof(REAL));
  } else {
    if(A->val){
      free(A->val);
      A->val = NULL;
    }
  }

  A->row=m; A->col=n; A->nnz=nnz;

  return;
}
/***********************************************************************/
/*!
 * \fn void icsr_realloc (const INT m, const INT n, const INT nnz, iCSRmat *A)
 *
 * \brief RE-allocate iCSRmat sparse matrix. Using realloc to expand,
 *        shrink the arrays in A. If upon entry the arrays are NULL it
 *        will allocate them (by the definition of realloc).
 *
 * \param m      Number of rows
 * \param n      Number of columns
 * \param nnz    Number of nonzeros
 * \param A      Pointer to the iCSRmat matrix
 *
 */
void icsr_realloc(const INT m,
		  const INT n,
		  const INT nnz,
		  iCSRmat *A)
{
  if ( m > 0 ) {
    A->IA=(INT*)realloc(A->IA,(m+1)*sizeof(INT));
  } else {
    if(A->IA) {
      free(A->IA);
      A->IA = NULL;
    }
  }
  //
  if ( n > 0 ) {
    A->JA=(INT*)realloc(A->JA,nnz*sizeof(INT));
  } else {
    if(A->JA){
      free(A->JA);
      A->JA = NULL;
    }
  }
  //
  if ( nnz > 0 ) {
    A->val=(INT *)realloc(A->val,nnz*sizeof(INT));
  } else {
    if(A->val){
      free(A->val);
      A->val = NULL;
    }
  }

  A->row=m; A->col=n; A->nnz=nnz;

  return;
}
/***********************************************************************************************/
/**
 * \fn dCOOmat dcoo_create (INT m, INT n, INT nnz)
 *
 * \brief Create IJ sparse matrix data memory space
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 *
 * \return A   The new dCOOmat matrix
 *
 */
dCOOmat dcoo_create (INT m,
                     INT n,
                     INT nnz)
{
    dCOOmat A;

    A.rowind = (INT *)calloc(nnz, sizeof(INT));
    A.colind = (INT *)calloc(nnz, sizeof(INT));
    A.val    = (REAL *)calloc(nnz, sizeof(REAL));

    A.row = m; A.col = n; A.nnz = nnz;

    return A;
}
/***********************************************************************************************/
/**
 * \fn void dcoo_alloc (const INT m, const INT n, const INT nnz, dCOOmat *A)
 *
 * \brief Allocate COO sparse matrix memory space
 *
 * \param m      Number of rows
 * \param n      Number of columns
 * \param nnz    Number of nonzeros
 * \param A      Pointer to the dCOOmat matrix
 *
 */
void dcoo_alloc (const INT m,
                 const INT n,
                 const INT nnz,
                 dCOOmat *A)
{

    if ( nnz > 0 ) {
        A->rowind = (INT *)calloc(nnz, sizeof(INT));
        A->colind = (INT *)calloc(nnz, sizeof(INT));
        A->val    = (REAL*)calloc(nnz,sizeof(REAL));
    }
    else {
        A->rowind = NULL;
        A->colind = NULL;
        A->val    = NULL;
    }

    A->row = m; A->col = n; A->nnz = nnz;

    return;
}

/***********************************************************************************************/
/**
 * \fn void dcoo_free (dCOOmat *A)
 *
 * \brief Free IJ sparse matrix data memory space
 *
 * \param A  Pointer to the dCOOmat matrix
 *
 */
void dcoo_free (dCOOmat *A)
{
    if (A==NULL) return;

    free(A->rowind); A->rowind= NULL;
    free(A->colind); A->colind = NULL;
    free(A->val);    A->val = NULL;
}

/***********************************************************************************************/
/**
 * \fn void icoo_free (dCOOmat *A)
 *
 * \brief Free IJ sparse matrix data memory space
 *
 * \param A  Pointer to the dCOOmat matrix
 *
 */
void icoo_free (iCOOmat *A)
{
    if (A==NULL) return;

    free(A->rowind); A->rowind= NULL;
    free(A->colind); A->colind = NULL;
    free(A->val);    A->val = NULL;
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
 * \param index_start   Number from which memory is indexed (1 for fortran, 0 for C)
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
  const unsigned INT n=A->row,nnz=A->nnz;
  const INT *ia=A->IA, *ja=A->JA;
  const REAL *Aval=A->val;
  INT i,j,k,jaj,i1,i2,start;

  dCSRmat Aperm = dcsr_create(n,n,nnz);

  // form the transpose of P
  INT *Pt = NULL;
  Pt =(INT *) calloc(n,sizeof(INT));

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
 * \fn void dcsr_transz (dCSRmat *A,  INT *p, dCSRmat *AT)
 *
 * \brief Generalized transpose of A: (n x m) matrix given in dCSRmat format
 *
 * \param A     Pointer to matrix in dCSRmat for transpose, INPUT
 * \param p     Permutation, INPUT
 * \param AT    Pointer to matrix AT = transpose(A)  if p = NULL, OR
 *                                AT = transpose(A)p if p is not NULL
 *
 * \note The storage for all pointers in AT should already be allocated,
 *       i.e. AT->IA, AT->JA and AT->val should be allocated before calling
 *       this function. If A.val=NULL, then AT->val[] is not changed.
 *
 * \note performs AT=transpose(A)p, where p is a permutation. If
 *       p=NULL then p=I is assumed. Applying twice this procedure one
 *       gets At=transpose(transpose(A)p)p = transpose(p)Ap, which is the
 *       same A with rows and columns permutted according to p.
 *
 * \note If A=NULL, then only transposes/permutes the structure of A.
 *
 * \note For p=NULL, applying this two times A-->AT-->A orders all the
 *       row indices in A in increasing order.
 *
 * Reference: Fred G. Gustavson. Two fast algorithms for sparse
 *            matrices: multiplication and permuted transposition.
 *            ACM Trans. Math. Software, 4(3):250–269, 1978.
 *
 * \author Ludmil Zikatanov
 * \date   19951219 (Fortran), 20150912 (C)
 */
void dcsr_transz (dCSRmat *A,
                  INT *p,
                  dCSRmat *AT)
{
    /* tested for permutation and transposition */
    /* transpose or permute; if A.val is null ===> transpose the
       structure only */
    const INT   n=A->row, m=A->col, nnz=A->nnz;
    const INT *ia=NULL,*ja=NULL;
    const REAL *a=NULL;
    INT m1=m+1;
    ia=A->IA; ja=A->JA; a=A->val;
    /* introducing few extra pointers hould not hurt too much the speed */
    INT *iat=NULL, *jat=NULL;
    REAL *at=NULL;

    /* loop variables */
    INT i,j,jp,pi,iabeg,iaend,k;

    /* initialize */
    AT->row=m; AT->col=n; AT->nnz=nnz;

    /* all these should be allocated or change this to allocate them here */
    iat=AT->IA; jat=AT->JA; at=AT->val;
    for (i = 0; i < m1; ++i) iat[i] = 0;
    iaend=ia[n];
    for (i = 0; i < iaend; ++i) {
        j = ja[i] + 2;
        if (j < m1) iat[j]++;
    }
    iat[0] = 0;
    iat[1] = 0;
    if (m != 1) {
        for (i= 2; i < m1; ++i) {
            iat[i] += iat[i-1];
        }
    }

    if (p && a) {
        /* so we permute and also use matrix entries */
        for (i = 0; i < n; ++i) {
            pi=p[i];
            iabeg = ia[pi];
            iaend = ia[pi+1];
            if (iaend > iabeg){
                for (jp = iabeg; jp < iaend; ++jp) {
                    j = ja[jp]+1;
                    k = iat[j];
                    jat[k] = i;
                    at[k] = a[jp];
                    iat[j] = k+1;
                }
            }
        }
    } else if (a && !p) {
        /* transpose values, no permutation */
        for (i = 0; i < n; ++i) {
            iabeg = ia[i];
            iaend = ia[i+1];
            if (iaend > iabeg){
                for (jp = iabeg; jp < iaend; ++jp) {
                    j = ja[jp]+1;
                    k = iat[j];
                    jat[k] = i;
                    at[k] = a[jp];
                    iat[j] = k+1;
                }
            }
        }
    } else if (!a && p) {
        /* Only integers and permutation (only a is null) */
        for (i = 0; i < n; ++i) {
            pi=p[i];
            iabeg = ia[pi];
            iaend = ia[pi+1];
            if (iaend > iabeg){
                for (jp = iabeg; jp < iaend; ++jp) {
                    j = ja[jp]+1;
                    k = iat[j];
                    jat[k] = i;
                    iat[j] = k+1;
                }
            }
        }
    } else {
        /* Only integers and no permutation (both a and p are null */
        for (i = 0; i < n; ++i) {
            iabeg = ia[i];
            iaend = ia[i+1];
            if (iaend > iabeg){
                for (jp = iabeg; jp < iaend; ++jp) {
                    j = ja[jp]+1;
                    k = iat[j];
                    jat[k] = i;
                    iat[j] = k+1;
                }
            }
        }
    }

    return;
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
 * \fn void icsr_concat(iCSRmat* A, iCSRmat* B, iCSRmat* C)
 *
 * \brief Concatenate two iCSRmat matrices of the same size.
 *
 * \param A Pointer to the first iCSRmat matrix
 * \param B Pointer to the second iCSRmat matrix
 * \param C Pointer to the concatenated matrix (Memory allocated in function)
 *
 * \note Concatenates columns, so the number of rows does not change.
 * \note Index of the arrays start at 0 -- Peter Ohm
 *
 */
void icsr_concat(iCSRmat* A,
                 iCSRmat* B,
                 iCSRmat* C)
{
  INT i,j;
  INT Arl, Brl;

  // Allocate for concatenation
  C->IA = (INT *)calloc(A->row + 1, sizeof(INT));//Concat so number of rows does not change
  C->JA = (INT *)calloc(A->nnz + B->nnz, sizeof(INT));
  C->val = (INT *)calloc(A->nnz + B->nnz, sizeof(INT));//allocation might be unneeded

  INT Astart = 0;
  INT Bstart = 0;
  INT cnt = 0;
  for(i=0; i<A->row; i++){
    // Get row length
    Arl = A->IA[i+1] - A->IA[i];
    Brl = B->IA[i+1] - B->IA[i];

    // Fill the row index pointer
    C->IA[i] = A->IA[i] + B->IA[i];

    // Fill the column index and val from A
    for(j=Astart; j<Astart+Arl; j++){
      C->JA[cnt] = A->JA[j];
      //C->val[cnt] = A->val[j];
      cnt++;
    }
    Astart = Astart+Arl;

    // Fill the column index and val from B
    for(j=Bstart; j<Bstart+Brl; j++){
      C->JA[cnt] = B->JA[j] + (A->col);
      //C->val[cnt] = B->val[j];
      cnt++;
    }
    Bstart = Bstart+Brl;
  }

  C->row = A->row;
  C->col = A->col + B->col;
  C->nnz = A->nnz + B->nnz;
  C->IA[C->row] = C->nnz;
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
 * modified ltz (20200811): replaced memset (-1) with an explicit loop)
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
    printf("### ERROR HAZMATH DANGER: both matrices are null, not sure if dimensions match!!! %s\n", __FUNCTION__);
    status = ERROR_MAT_SIZE;
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

  // initialize C->IA
  memset(C->IA, 0, sizeof(INT)*(C->row+1));
  memset(C->JA, -1, sizeof(INT)*(A->nnz+B->nnz));
  for (i=0; i<(A->nnz+B->nnz); ++i) {
    C->JA[i]=-1;
  }
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
 * \fn dcsr_mxv_forts (void *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *        Used for time-stepping routines
 *
 * \param A   Pointer to dCSRmat matrix A
 * \param x   Pointer to array x
 * \param y   Pointer to array y
 *
 *
 */
void dcsr_mxv_forts(void *At,
                      REAL *x,
                      REAL *y)
{
  // Declare A as dCSRmat
  dCSRmat *A = (dCSRmat *) At;

  // Perform Matrix Vector Multiply
  dcsr_mxv(A,x,y);
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
 * \fn void dcsr_mxm (dCSRmat *A, dCSRmat *B, dCSRmat *C)
 *
 * \brief Sparse matrix multiplication C=A*B (index starts with 0!!)
 *
 * \param A   Pointer to the dCSRmat matrix A
 * \param B   Pointer to the dCSRmat matrix B
 * \param C   Pointer to dCSRmat matrix equal to A*B
 *
 * modified ltz (20200811). replaces memset to (JD,-1...) with
 *                          explicit loop as memset may be cannot set
 *                          to -1.
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

    for (j=0;j<countJD;++j) JD[j]=-1;
    //    memset(JD, -1, sizeof(INT)*countJD);
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
   * \fn void icsr_mxm_symb (iCSRmat *A, iCSRmat *B, iCSRmat *C)
   *
   * \brief Symbolic sparse matrix multiplication C=A*B (index starts with 0!!)
   *
   * \param A   Pointer to the iCSRmat matrix A
   * \param B   Pointer to the iCSRmat matrix B
   * \param C   Pointer to iCSRmat matrix equal to A*B
   *
   */
void icsr_mxm_symb(iCSRmat *A,iCSRmat *B,iCSRmat *C)
{
  INT np=A->row,nr=B->col; // INT nq=A->col;
  INT i,j,k,nnz,jp,kp;
  C->val=NULL;C->JA=NULL;
  INT *ia=A->IA, *ja=A->JA,*ib=B->IA, *jb=B->JA;
  INT *ic, *jc;
  INT *ix=(INT *)calloc(nr,sizeof(INT));
  nnz = 0;
  for(i=0;i<nr;i++)ix[i]=-1;
  for(i=0;i<np;i++){
    for(jp=ia[i];jp<ia[i+1];jp++){
      j = ja[jp];
      for(kp = ib[j];kp<ib[j+1];kp++){
	k = jb[kp];
	if(ix[k]!=i) {
	  nnz++;
	  ix[k] = i;
	}
      }
    }
  }
  ic=(INT *)calloc((np+1),sizeof(INT));
  jc=(INT *)calloc(nnz,sizeof(INT));
  nnz=0;
  for(i=0;i<nr;++i) ix[i]=-1;
  for(i=0;i<np;i++){
    ic[i] = nnz;
    for(jp=ia[i];jp<ia[i+1];jp++){
      j = ja[jp];
      for(kp = ib[j];kp<ib[j+1];kp++){
	k = jb[kp];
	if(ix[k]!=i) {
	  jc[nnz] = k;
	  nnz++;
	  ix[k] = i;
	}
      }
    }
  }
  ic[np] = nnz;
  free(ix);
  C->nnz=nnz;
  C->row=np; C->col=nr;
  /**/
  C->IA=ic;
  C->JA=jc;
  C->val=NULL;
  return;
}
/**************************************************************************/
/*!
   * \fn void icsr_mxm (iCSRmat *A, iCSRmat *B, iCSRmat *C, INT multmax)
   *
   * \brief symbolic sparse matrix multiplication C=A*B (index starts with 0!!)
   *
   * \param A         Pointer to the iCSRmat matrix A
   * \param B         Pointer to the iCSRmat matrix B
   * \param C         Pointer to iCSRmat matrix equal to A*B
   *
   */
void icsr_mxm (iCSRmat *A,iCSRmat *B,iCSRmat *C)
{
  /* C-------------------------------------------------------------------- */
  /* call the symbolic muultiplication */
  icsr_mxm_symb(A,B,C);
  /**/
  INT np=C->row,nr=C->col,nnz=C->nnz; // INT nq=A->col;
  INT i,j,k,jp,kp,aij;
  INT *ia=A->IA, *ja=A->JA,*ib=B->IA, *jb=B->JA;
  INT *ic=C->IA, *jc=C->JA;
  INT *a=A->val, *b=B->val;
  /**/
  INT *c=(INT *)calloc(nnz,sizeof(INT));
  INT *ix=(INT *)calloc(nr,sizeof(INT));
  /*do the multipplication*/
  /* C===================================================================== */
/*       subroutine abyb(ia,ja,ib,jb,np,ic,jc,an,bn,cn,x,nq) */
/* C===================================================================== */
/*       implicit real*8 (a-h,o-z) */
/*       dimension ia(1),ja(1),ib(1),jb(1),ic(1),jc(1) */
/*       dimension an(1),bn(1),cn(1),x(nq) */
/* C-------------------------------------------------------------------- */
/* C...  MULTIPLICATION OF TWO GENERAL SPARSE MATRICIES: C = A*B */
/* C-------------------------------------------------------------------- */
//  ix=realloc(ix,nr*sizeof(INT));
  for(i=0;i<np;i++){
    for(j = ic[i];j<ic[i+1];j++) ix[jc[j]] = 0;
    for(jp=ia[i];jp<ia[i+1];jp++){
      j = ja[jp];
      aij = a[jp];
      for(kp=ib[j];kp<ib[j+1];kp++){
	k = jb[kp];
        ix[k]+=aij*b[kp];
      }
    }
    for(j=ic[i];j<ic[i+1];j++) c[j] = ix[jc[j]];
  }
  ic[np]=nnz;
  free(ix);
  C->val=c;
  C->nnz=nnz;
  return;
}
/**************************************************************************/
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
   */
void icsr_mxm_symb_max(iCSRmat *A,iCSRmat *B,iCSRmat *C,	\
			INT multmax)
{
  /*No need of C->val here*/
  /* call the symbolic muultiplication. this is same as A*B, but we do
     not call A*B bc there is no multiplication here aij=bij=1 */
  icsr_mxm_symb(A,B,C);
  /**/
  INT np=C->row,nr=C->col,nnz=C->nnz; // INT nq=A->col;
  INT izz,i,j,k,jp,kp;
  INT *ia=A->IA, *ja=A->JA,*ib=B->IA, *jb=B->JA;
  INT *ic=C->IA, *jc=C->JA;
  /**/
  INT *ix=(INT *)calloc(nr,sizeof(INT));
  // skip everything that is not equal to multmax.
  nnz=0;
  for(i=0;i<np;i++){
    izz=ic[i];
    ic[i]=nnz;
    for(j = izz;j<ic[i+1];j++){
      ix[jc[j]] = 0;
    }
    for(jp=ia[i];jp<ia[i+1];jp++){
      j = ja[jp];
      //      aij = a(jp)
      for(kp=ib[j];kp<ib[j+1];kp++){
	k = jb[kp];
        ix[k]++; // +=aij*b(kp)
      }
    }
    for(j=izz;j<ic[i+1];j++){
      if(ix[jc[j]]!=multmax) continue;
      jc[nnz]=jc[j];
      nnz++;
    }
  }
  ic[np]=nnz;
  free(ix);
  /**/
  C->JA=realloc(jc,nnz*sizeof(INT));
  C->nnz=nnz;
  return;
}
/***********************************************************************************************/
/*!
   * \fn dCSRmat dcsr_create_diagonal_matrix(dvector *diag)
   *
   * \brief create diagonal matrix using diag as diagonal entres
   *
   * \param diag  Pointer to the diagonal as a dvector
   *
   * \return D     Pointer to dCSRmat CSR diagonal matrix
   *
   */
dCSRmat dcsr_create_diagonal_matrix(dvector *diag)
{
  //local variable
  INT n = diag->row;
  INT i;

  // form the diaongal matrix
  dCSRmat D = dcsr_create(n,n,n);
  for (i=0;i<n;i++)
  {
    D.IA[i] = i;
    D.JA[i] = i;
    D.val[i]   = diag->val[i];
  }
  D.IA[n] = n;

  // return
  return D;
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
   * \fn void dcsr_getdiag_pow (INT n, REAL p, dCSRmat *A, dvector *diag)
   *
   * \brief Get first n diagonal entries of a dCSRmat sparse matrix A
   *        to the power p
   *
   * \param n     Number of diagonal entries to get (if n=0, then get all diagonal entries)
   * \param p     Exponential power
   * \param A     Pointer to dCSRmat CSR matrix
   * \param diag  Pointer to the diagonal as a dvector
   *
   * \note Added by Ana Budisa 2020-05-20; for fractional AMG smoothers only
   */
void dcsr_getdiag_pow(INT n,
                      REAL p,
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
        diag->val[i] = pow(A->val[k], p);
        break;
      } // end if
    } // end for k
  } // end for i

}

/***********************************************************************************************/
/*!
   * \fn void dcsr_shiftdiag (INT n, dCSRmat *A, dvector *shifts)
   *
   * \brief Shift first n diagonal entries of a dCSRmat sparse matrix A by some real values
   *
   * \param n      Number of diagonal entries to get (if n=0, then get all diagonal entries)
   * \param A      Pointer to dCSRmat CSR matrix
   * \param shifts dvector of length n with REAL value shifts
   *
   */
void dcsr_shiftdiag(INT n,
                    dCSRmat *A,
                    dvector *shifts)
{
  INT i,k,j,ibegin,iend;

  if ( n==0 || n>A->row || n>A->col ) n = MIN(A->row,A->col);

  for (i=0;i<n;++i) {
    ibegin=A->IA[i]; iend=A->IA[i+1];
    for (k=ibegin;k<iend;++k) {
      j=A->JA[k];
      if ((j-i)==0) {
        A->val[k] += shifts->val[i];
        break;
      } // end if
    } // end for k
  } // end for i

}
/***********************************************************************************************/
/*!
   * \fn void dcsr_row_scale(dCSRmat *A, dvector *row_scale)
   *
   * \brief row scaling of matrix A
   *
   * \param A           Pointer to dCSRmat CSR matrix
   * \param diag_scale  Pointer to the diagonal scaling
   *
   */
void dcsr_row_scale(dCSRmat *A,
                    dvector *row_scale)
{
  // local variables
  INT i,k,ibegin,iend;
  INT n= A->row;

  // TODO: need to check the size of A and row_scale -- Xiaozhe Hu

  for (i=0;i<n;++i) {
    ibegin=A->IA[i]; iend=A->IA[i+1];
    for (k=ibegin;k<iend;++k) {
      A->val[k] = row_scale->val[i]*A->val[k];
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
   * \fn SHORT dcsr_delete_rowcol(dCSRmat *A,
                                  INT *delete_row,
                                  INT *delete_col,
                                  dCSRmat *B)
   *
   * \brief Delete selected specified rows and columns from a dCSRmat sparse matrix A
   *
   * \param A             Pointer to dCSRmat matrix
   * \param B             Pointer to dCSRmat matrix
   * \param delete_row    Pointer to rows that need to be deleted
   * \param delete_col    Pointer to cols that need to be deleted
   *
   * \return      SUCCESS if succeeded, otherwise return error information.
   *
   * \note Index starts at 0!!! -- Xiaozhe Hu
   *
   */
SHORT dcsr_delete_rowcol(dCSRmat *A,
                         INT *delete_row,
                         INT *delete_col,
                         dCSRmat *B)
{

  // return variable
  INT status = SUCCESS;

  // local variables
  INT i;
  INT row_count, col_count;
  INT *row_stay, *col_stay;

  // allocate for row and columns that stays
  row_stay = (INT*)calloc(A->row,sizeof(INT));
  col_stay = (INT*)calloc(A->col,sizeof(INT));

  // generate row_stay
  row_count = 0;
  for (i=0; i<A->row; i++) {

    // get rows that stay
    if (delete_row[i] == 0) {
      row_stay[row_count] = i;
      row_count++;
    }

  }
  // reallocate
  row_stay = (INT*)realloc(row_stay, sizeof(INT)*row_count);

  // generate col_stay
  col_count = 0;
  for (i=0; i<A->col; i++) {

    // get rows that stay
    if (delete_col[i] == 0) {
      col_stay[col_count] = i;
      col_count++;
    }

  }
  // reallocate
  col_stay = (INT*)realloc(col_stay, sizeof(INT)*col_count);

  // get submatrix from original matrix
  status = dcsr_getblk(A, row_stay, col_stay, row_count, col_count, B);

  // free memory
  free(row_stay);
  free(col_stay);

  // return
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
/**
 * \fn dCSRmat dcsr_sympat(dCSRmat *A)
 * \brief Get symmetric part of a dCSRmat matrix
 *
 * \param *A      pointer to the dCSRmat matrix
 *
 * \return symmetrized the dCSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date 03/21/2011
 */
dCSRmat dcsr_sympat (dCSRmat *A)
{
    //local variable
    dCSRmat AT;

    //return variable
    dCSRmat SA;

    // get the transpose of A
    dcsr_trans(A,  &AT);

    // get symmetrized A
    dcsr_add(A, 0.5, &AT, 0.5, &SA);

    // clean
    dcsr_free(&AT);

    // return
    return SA;
}

/***********************************************************************************************/
/**
 * \fn dCSRmat dcsr_reorder(dCSRmat *A, INT *order)
 * \brief Reorder a dCSRmat matrix
 *
 * \param *A      pointer to the dCSRmat matrix
 * \param *order  pointer to the new order
 *
 * \return reordered dCSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date 05/16/2011
 */
dCSRmat dcsr_reorder(dCSRmat *A,
                      INT *order)
 {
    // local variables
    dCOOmat A_coo;
    dCSRmat B;
    const INT nnz = A->nnz;

    INT i;

    // convert to dCOOmat
    dcsr_2_dcoo(A, &A_coo);

    // map
    INT *map = (INT *)calloc(A->row, sizeof(INT));
    for (i=0; i< A->row; i++){
        map[order[i]] = i;
    }

    // reorder dCOOmat
    for (i=0; i<nnz; i++){
        A_coo.rowind[i] = map[A_coo.rowind[i]];
        A_coo.colind[i] = map[A_coo.colind[i]];
    }

    // convert back to dCSRmat
    dcoo_2_dcsr(&A_coo, &B);

    // clean
    if(map) free(map);
    dcoo_free(&A_coo);

    // return
    return B;

 }

/**********************************************************************/
/*removing diagonal or extracting upper/lower triangle of a sparse
  matrix*/
/**********************************************************************/
/***********************************************************************************************/
/**
 * \fn icsr_nodiag(iCSRmat *a)
 *
 * \brief removing the diagonal of an icsr matrix and cleans the off
 * diagonal zeros as well (done in-inplace, a IS OVERWRITTEN)
 *
 * \param *a      pointer to the iCSRmat matrix
 *
 * \return void
 *
 * \author Ludmil
 * \date 20210101
 */
void icsr_nodiag(iCSRmat *a)
{
  INT k,j,kj,kj0,kj1;
  a->nnz=a->IA[0];
  for(k=0;k<a->row;k++){
    kj0=a->IA[k];
    kj1=a->IA[k+1];
    a->IA[k]=a->nnz;
    for(kj=kj0;kj<kj1;kj++){
      j=a->JA[kj];
      if((k!=j) && (a->val[kj]!=0)){
	a->JA[a->nnz]=j;
	a->val[a->nnz]=a->val[kj];
	a->nnz++;
      }
    }
  }
  a->IA[a->row]=a->nnz;
  if(a->nnz>0){
    a->JA=realloc(a->JA,a->nnz*sizeof(INT));
    a->val=realloc(a->val,a->nnz*sizeof(INT));
  }else{
    free(a->JA);
    free(a->val);
    a->JA=NULL;
    a->val=NULL;
  }
  return;
}
/***********************************************************************************************/
/**
 * \fn void icsr_tri(iCSRmat *a, const char loup)
 *
 * \brief extracting the lower/upper triangle of an icsr matrix (done
 * in-place, a is overwritten also removes all zeros).
 *
 * \param *a      pointer to the iCSRmat matrix
 *
 * \param character corresponding to a different action: loup='u' or
 *        'U': upper triangle; loup='l' or 'L': extract lower
 *        triangle; if loup is anything else: do nothing; the diagonal
 *        * is included and not removed (only zero elements are
 *        removed).
 *
 *
 * \return void; a is modified.
 *
 * \author Ludmil
 * \date 20210101
 */
/**********************************************************************/
void icsr_tri(iCSRmat *a,const char loup)
{
  INT lu;
  if(loup=='u' || loup=='U')
    lu=1;
  else if(loup=='l' || loup=='L')
    lu=-1;
  else
    return;
  INT k,j,kj,kj0,kj1;
  a->nnz=a->IA[0];
  for(k=0;k<a->row;k++){
    kj0=a->IA[k];
    kj1=a->IA[k+1];
    a->IA[k]=a->nnz;
    for(kj=kj0;kj<kj1;kj++){
      j=a->JA[kj];
      //      if(k<=j) for upper; (k>=j) for lower;
      if( ((k-j)*lu>0) || (a->val[kj]==0) ) continue; // lower is rowind>=colind
      a->JA[a->nnz]=j;
      a->val[a->nnz]=a->val[kj];
      a->nnz++;
    }
  }
  a->IA[a->row]=a->nnz;
  if(a->nnz>0){
    a->JA=realloc(a->JA,a->nnz*sizeof(INT));
    a->val=realloc(a->val,a->nnz*sizeof(INT));
  }else{
    free(a->JA);
    free(a->val);
    a->JA=NULL;
    a->val=NULL;
  }
  return;
}
/***********************************************************************************************/
/**
 * \fn void icsr_rm_value(iCSRmat *a,const INT value)
 *
 * \brief  REMOVES all entries in an iscr matrix EQUAL to a given value (done
     in-inplace, a is overwritten)
 *
 * \param *a      pointer to the iCSRmat matrix
 *
 * \param value ; removes all a(i,j)=value;
 *
 *
 * \return void; a is modified.
 *
 * \author Ludmil
 * \date 20210101
 */
/**********************************************************************/
void icsr_rm_value(iCSRmat *a,const INT value)
{
  INT k,j,kj,kj0,kj1;
  a->nnz=a->IA[0];
  for(k=0;k<a->row;k++){
    kj0=a->IA[k];
    kj1=a->IA[k+1];
    a->IA[k]=a->nnz;
    for(kj=kj0;kj<kj1;kj++){
      j=a->JA[kj];
      if(a->val[kj]!=value){
	a->JA[a->nnz]=j;
	a->val[a->nnz]=a->val[kj];
	a->nnz++;
      }
    }
  }
  a->IA[a->row]=a->nnz;
  if(a->nnz>0){
    a->JA=realloc(a->JA,a->nnz*sizeof(INT));
    a->val=realloc(a->val,a->nnz*sizeof(INT));
  }else{
    a->JA=NULL;
    a->val=NULL;
  }
  return;
}
/***********************************************************************************************/
/**
 * \fn void icsr_keep_value(iCSRmat *a,const INT value)
 *
 * \brief  KEEPS ONLY entries in an iscr matrix EQUAL to a given value (done
     in-inplace, a is overwritten)
 *
 * \param *a      pointer to the iCSRmat matrix
 *
 * \param value ; KEEPS all a(i,j)=value. Removes the rest;
 *
 * \return void; a is modified.
 *
 * \author Ludmil
 * \date 20210101
 */
/**********************************************************************/
void icsr_keep_value(iCSRmat *a,const INT value)
{
  /* KEEPS only entries in an iscr matrix equal to value (done
     in-inplace, a is overwritten) */
  INT k,j,kj,kj0,kj1;
  a->nnz=a->IA[0];
  for(k=0;k<a->row;k++){
    kj0=a->IA[k];
    kj1=a->IA[k+1];
    a->IA[k]=a->nnz;
    for(kj=kj0;kj<kj1;kj++){
      j=a->JA[kj];
      if(a->val[kj]!=value)
	continue;
      a->JA[a->nnz]=j;
      a->val[a->nnz]=a->val[kj];
      a->nnz++;
    }
  }
  a->IA[a->row]=a->nnz;
  if(a->nnz>0){
    a->JA=realloc(a->JA,a->nnz*sizeof(INT));
    a->val=realloc(a->val,a->nnz*sizeof(INT));
  }else{
    a->JA=NULL;
    a->val=NULL;
  }
  return;
}
/***********************************************************************************************/
// Block_dCSRmat subroutines starts here!
/***********************************************************************************************/
/*!
   * \fn void bdcsr_alloc_minimal(const INT brow, const INT bcol, block_dCSRmat *A)
   *
   * \brief Allocate block dCSRmat sparse matrix memory space (minimal spaces that are needed)
   *
   * \param brow      Number of block rows
   * \param bcol      Number of block columns
   * \param A         Pointer to the dCSRmat matrix
   *
   * \note: this only allocates A->blocks, but does not allocate either each A->block[i] nor each dcsr matrix
   *
   */
void bdcsr_alloc_minimal(const INT brow,
                         const INT bcol,
                         block_dCSRmat *A)
{
    //SHORT i;

    A->brow = brow;
    A->bcol = bcol;

    if (brow == 0 || bcol == 0){
        A->blocks = NULL;
    }
    else {
        A->blocks = (dCSRmat **) calloc(brow*bcol,sizeof(dCSRmat *));
    }

  return;
}

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
   * \note: this allocates A->blocks and each A->block[i], but does not allocate each dcsr matrix
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
   * \fn void bdcsr_free_minimal (block_dCSRmat *A)
   *
   * \brief Free block dCSR sparse matrix data (which is allocated by bdcsr_alloc_minimal)
   *
   * \param A   Pointer to the block_dCSRmat matrix
   *
   */
void bdcsr_free_minimal(block_dCSRmat *A)
{
  if (A == NULL) return; // Nothing need to be freed!

  INT i;
  INT num_blocks = (A->brow)*(A->bcol);

  for ( i=0; i<num_blocks; i++ ) {
    dcsr_free(A->blocks[i]);
  }

  if(A->blocks) {
      free(A->blocks);
      A->blocks = NULL;
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
    if(A->blocks[i]) {
        free(A->blocks[i]);
        A->blocks[i] = NULL;
    }
  }

  if(A->blocks) {
      free(A->blocks);
      A->blocks = NULL;
  }
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
            if (B->blocks[i]) {
                free(B->blocks[i]);
                B->blocks[i] = NULL;
            }
        } else {
            dcsr_alloc(A->blocks[i]->row, A->blocks[i]->col, A->blocks[i]->nnz, B->blocks[i]);
            dcsr_cp(A->blocks[i], B->blocks[i]);
        }
    }
}

/***********************************************************************************************/
/*!
   * \fn void bdcsr_trans (block_dCSRmat *A, block_dCSRmat *AT)
   *
   * \brief Transpose a block_dCSRmat matrix A
   *
   * \param A   Pointer to the block_dCSRmat matrix
   * \param AT  Pointer to the transpose of block_dCSRmat matrix A (output)
   *
   */
void bdcsr_trans(block_dCSRmat *A,
                 block_dCSRmat *AT)
{

  // local variables
  INT i,j;

  // allocate AT
  bdcsr_alloc_minimal(A->bcol, A->brow, AT);

  // transpose
  for (i=0; i<AT->brow; i++)
  {

    for (j=0; j<AT->bcol; j++)
    {

      if (A->blocks[j*A->bcol+i] == NULL)
      {
        AT->blocks[i*AT->bcol+j] = NULL;
      }
      else
      {
        AT->blocks[i*AT->bcol+j] = (dCSRmat *)calloc(1,sizeof(dCSRmat));
        dcsr_trans(A->blocks[j*A->bcol+i], AT->blocks[i*AT->bcol+j]);
      }

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
   * \return Flag of whether the adding is succesful or not (SUCCESS: 0; FAIL: <0)
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
                status = dcsr_add(NULL, alpha, B->blocks[i*A->brow+j], beta, C->blocks[i*A->brow+j]);
                if (status < 0) {goto FINISHED;}
            }
        }
    }

    // only matrices B is NULL
    if (B == NULL) {
        for (i=0; i<A->brow; i++){
            for (j=0; j<A->bcol; j++){
                status = dcsr_add(A->blocks[i*A->brow+j], alpha, NULL, beta, C->blocks[i*A->brow+j]);
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
            if (A->blocks[i*A->brow+j]==NULL && B->blocks[i*A->brow+j] == NULL) {
                free(C->blocks[i*A->brow+j]);
                C->blocks[i*A->brow+j]=NULL;
                status = SUCCESS;
            }
            else {
                status = dcsr_add(A->blocks[i*A->brow+j], alpha, B->blocks[i*A->brow+j], beta, C->blocks[i*A->brow+j]);
            }
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

  k=0;
  for(i=0; i<brow; i++){
    for(j=0; j<A->blocks[i*brow+i]->row; j++){
      y[k] = 0.0;
      k++;
    }
  }


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
}

/***********************************************************************************************/
/*!
   * \fn SHORT bdcsr_delete_rowcol(block_dCSRmat *A,
                                   INT *delete_row,
                                   INT *delete_col,
                                   block_dCSRmat *B)
   *
   * \brief Delete selected specified rows and columns from a block_dCSRmat sparse matrix A
   *
   * \param A             Pointer to block_dCSRmat matrix
   * \param B             Pointer to block_dCSRmat matrix
   * \param delete_row    Pointer to rows that need to be deleted (global index)
   * \param delete_col    Pointer to cols that need to be deleted (global index)
   *
   * \return      SUCCESS if succeeded, otherwise return error information.
   *
   * \note Index starts at 0!!! -- Xiaozhe Hu
   *
   */
SHORT bdcsr_delete_rowcol(block_dCSRmat *A,
                         INT *delete_row,
                         INT *delete_col,
                         block_dCSRmat *B)
{
  // return variable
  INT status = SUCCESS;

  // local variables
  INT brow = A->brow;
  INT bcol = A->bcol;

  INT i,j;
  INT *row_start, *col_start;

  // allocate row_start and col_start
  row_start = (INT *)calloc(brow,sizeof(INT));
  col_start = (INT *)calloc(bcol,sizeof(INT));

  // allocate B
  bdcsr_alloc(brow, bcol, B);

  // get row_start
  row_start[0] = 0;
  for (i=0; i<brow-1; i++){
    for (j=0; j<bcol; j++){
      if (A->blocks[i*brow+j]) {
        row_start[i+1] = row_start[i] + A->blocks[i*brow+j]->row;
        break;
      }
    }
  }

  // get col_start
  col_start[0] = 0;
  for (j=0; j<bcol-1; j++){
    for (i=0; i<brow; i++){
      if (A->blocks[i*brow+j]){
        col_start[j+1] = col_start[j] + A->blocks[i*brow+j]->col;
        break;
      }
    }
  }

  // main loop
  for (i=0; i<brow; i++){
    for (j=0; j<bcol; j++){
      if (A->blocks[i*brow+j]) {
        status = dcsr_delete_rowcol(A->blocks[i*brow+j], delete_row+row_start[i],delete_col+col_start[j], B->blocks[i*brow+j]);
      }
      else {
        B->blocks[i*brow+j] = NULL;
      }
    }
  }

  // free memory
  free(row_start);
  free(col_start);

  // return
  return(status);

}

/***********************************************************************************************/
/*!
 * \fn void bdcsr_mxm (block_dCSRmat *A, block_dCSRmat *B, block_dCSRmat *C)
 *
 * \brief Sparse matrix multiplication C=A*B (block_dCSRmat format)
 *
 * \param A   Pointer to the block_dCSRmat matrix A
 * \param B   Pointer to the block_dCSRmat matrix B
 * \param C   Pointer to the block_dCSRmat matrix equal to A*B
 *
 */
// void bdcsr_mxm(block_dCSRmat *A,
//                block_dCSRmat *B,
//                block_dCSRmat *C)
// {
//     // get size
//     INT A_brow = A->brow;
//     INT A_bcol = A->bcol;
//     INT B_brow = B->brow;
//     INT B_bcol = B->bcol;
//
//     // check A->bcol = B->brow
//     if (A_bcol != B_brow){
//         printf("HAZMATH DANGER!!! Matrice sizes do not match!!\n");
//         return;
//     }
//
//     // local variables
//     INT i,j,k;
//     dCSRmat temp_mat;
//
//     // allocate C
//     bdcsr_alloc(A->brow, B->bcol, C);
//
//     // main loop
//     for(i=0; i<C->brow; i++){
//         for(j=0; j<C->bcol; j++){
//             // C[i][j] = \sum_k A[i][k]*B[k][j]
//             for (k=0; k<A->bcol; k++){
//                 //A[i][k]*B[k][j]
//                 dcsr_mxm(A->blocks[i*bcol+k], B->blocks[k*bcol+j], &temp_mat);
//                 //
//
//             }
//         }
//     }
//
// }

/***********************************************************************************************/


/* sparse matrix functions returning pointers. */
/*******************************************************************/
/*!
 * \fn dCSRmat *dcsr_create_p(const INT m, const INT n, const INT nnz)
 *
 * \brief Create a dCSRmat sparse matrix. Uses void array for the
 * whole matrix. the void array contains in first position the struct.
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 *
 * \return A a pointer to a dCSRmat matrix. All of this matrix can be
 *         freed by free((void *)A) or even free(A).
 *
 *  \note: modified by ludmil zikatanov on 20200412
 */
dCSRmat *dcsr_create_p (const INT m,		\
			const INT n,		\
			const INT nnz)
{
  dCSRmat *A=NULL;
  size_t structby=sizeof(dCSRmat);// size of the struct
  size_t realby=sizeof(REAL),intby=sizeof(INT);// size of ints and reals
  size_t total=1*structby+3*intby; //at least space for structure.
  if ( m > 0 )
    total+=(m+1)*intby;
  if ( n > 0 )
    total+=nnz*intby;
  if ( nnz > 0 )
    total+=nnz*realby;
  void *w=(void *)calloc(total/sizeof(char),sizeof(char));
  A=(dCSRmat *)w;
  w+=1*structby;
  A->IA = NULL;
  A->JA = NULL;
  A->val = NULL;
  INT *mn_nnz=(INT *)w;
  A->row=mn_nnz[0]=m; A->col=mn_nnz[1]=n; A->nnz=mn_nnz[2]=nnz;
  w+=3*intby;
  if ( m > 0 ) {
    A->IA = (INT *)w;
    w+=(m+1)*intby;
  }
  if ( n > 0 ) {
    A->JA = (INT *)w;
    w+=nnz*intby;
  }
  if ( nnz > 0 ) {
    A->val = (REAL *)w;
    w+=nnz*realby;// end of it.
  }
  return A;
}
/***********************************************************************************************/
/**
 * \fn dCOOmat *dcoo_create_p(INT m, INT n, INT nnz)
 *
 * \brief Create IJ sparse matrix data memory space using one
 * contguous void array for all data including the structure itself
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 *
 * \return A   A pointer to a dCOOmat matrix
 *
 */
dCOOmat *dcoo_create_p(INT m,			\
		       INT n,			\
		       INT nnz)
{
  size_t structby=sizeof(dCOOmat),realby=sizeof(REAL),intby=sizeof(INT);
  size_t total=(1*structby+(3+2*nnz)*intby+nnz*realby)/sizeof(char);
  void *w=(void *)calloc(total, sizeof(char));
  //sturture
  dCOOmat *A=(dCOOmat *)w;
  w+=1*structby;
  INT *mn_nnz=(INT *)w;
  A->row=mn_nnz[0]=m; A->col=mn_nnz[1]=n; A->nnz=mn_nnz[2]=nnz;
  w+=3*intby;
  // arrays;
  A->rowind = (INT *)w;
  w+=nnz*intby;
  A->colind = (INT *)w;
  w+=nnz*intby;
  A->val    = (REAL *)w;
  w+=nnz*realby; // end of it....
  return A;
}

/******************************************************************/
/*!
 * \fn dvector dvec_create_p(const INT m)
 *
 * \brief Create a dvector of given length: use void array to pack the
 * structure in.
 *
 * \param m    length of the dvector
 *
 * \return pointer to u   The new dvector
 *
 */
dvector *dvec_create_p(const INT m)
{
  dvector *u=NULL;
  size_t structby=sizeof(dvector);// size of the struct
  size_t realby=sizeof(REAL),intby=sizeof(INT);// size of ints and reals
  size_t total=1*structby+1*intby; //space for structure and size.
  if (m > 0 )
    total+=m*realby;
  void *w=(void *)calloc(total/sizeof(char),sizeof(char));
  u=(dvector *)w;
  w+=1*structby;
  INT *mm=(INT *)w;
  u->row = mm[0]=m;
  w+=1*intby;
  u->val = NULL;
  if ( m > 0 ) {
    u->val = (REAL *)w;
    w+=m*realby;//end
  }
  return u;
}
/************************************************************************/
/**
 * \fn iCSRmat icsr_create_p (const INT m, const INT n, const INT nnz)
 *
 * \brief Create iCSRmat sparse matrix using a void array; the
 *        structure itself is part of the void array.
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 *
 * \return A   a pointer to an iCSRmat matrix
 *
 */
iCSRmat *icsr_create_p(const INT m,		\
		       const INT n,		\
		       const INT nnz)
{
  iCSRmat *A=NULL;
  size_t structby=sizeof(iCSRmat);// size of the struct
  size_t intby=sizeof(INT);// size of ints
  size_t total=1*structby+3*intby; //space for the structure.
  if ( m > 0 )
    total+=(m+1)*intby;
  if ( n > 0 )
    total+=nnz*intby;
  if ( nnz > 0 )
    total+=nnz*intby;
  void *w=(void *)calloc(total/sizeof(char),sizeof(char));
  A=(iCSRmat *)w;
  w+=1*structby;
  A->IA = NULL;
  A->JA = NULL;
  A->val = NULL;
  INT *mn_nnz=(INT *)w;
  A->row=mn_nnz[0]=m; A->col=mn_nnz[1]=n; A->nnz=mn_nnz[2]=nnz;
  w+=3*intby;
  if ( m > 0 ) {
    A->IA = (INT *)w;
    w+=(m+1)*intby;
  }
  if ( n > 0 ) {
    A->JA = (INT *)w;
    w+=nnz*intby;
  }
  if ( nnz > 0 ) {
    A->val = (INT *)w;
    w+=nnz*intby; // end of it
  }
  return A;
}
/************************************************************/
/*!
 * \fn ivector *ivec_create_p(const INT m)
 *
 * \brief Create an ivector of given length
 *
 * \param m   length of the ivector
 *
 * \return u  The new ivector
 *
 */
ivector *ivec_create_p(const INT m)
{
  ivector *u=NULL;
  size_t structby=sizeof(ivector);// size of the struct
  size_t intby=sizeof(INT);// size of ints and reals
  size_t total=1*structby+1*intby; //space for structure.
  if (m > 0 )
    total+=m*intby;
  void *w=(void *)calloc(total/sizeof(char),sizeof(char));
  u=(ivector *)w;
  w+=1*structby;
  INT *mm=(INT *)w;
  u->row = mm[0]=m;
  w+=1*intby;
  u->val = NULL;
  if ( m > 0 ) {
    u->val = (INT *)w;
    w+=m*intby;// end of struct.
  }
  return u;
}
/**********************************************************************/
/*!
 * \fn void dcsr_alloc_p(const INT m, const INT n, const INT nnz,
 * dCSRmat **A)
 *
 * \brief Allocate a dCSRmat sparse matrix and put the result in the
 *        pointer pointed by A. Same as dCSRmat_create_p, but uses
 *        realloc() to reallocate A.
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 *
 * \return A a pointer to a an array of dCSRmat matrices with one
 *         element. All of this matrix can be freed by free((void *)A)
 *         or even just free(A).
 *
 *  \note: call as dCSRmat *X; dcsr_alloc_p(... &X) to reallocate *X
 *         to desired length.
 *         modified by ludmil zikatanov on 20200412
 */
void dcsr_alloc_p (const INT m,			\
		   const INT n,			\
		   const INT nnz,		\
		   dCSRmat **A)
{
  *A=dcsr_create_p(m,n,nnz);
  return;
}
/******************************************************************/
/*!
 * \fn dvector dvec_alloc_p(const INT m, dvector **u)
 *
 * \brief allocate an "one element array" of dvectors **u.
 *
 * \param m    length of the dvector
 * \param u    pointer to a dvector that is to be reallocated.
 *
 * \note: call as: dvector *u; dvec_alloc_p(m,&u);
 *
 */
void dvec_alloc_p(const INT m, dvector **u)
{
  *u=dvec_create_p(m);
  return;
}
/************************************************************/
/*!
 * \fn void ivec_alloc_p(const INT m, ivector **u)
 *
 * \brief Create an ivector of given length
 *
 * \param m   length of the ivector
 *
 *
 */
void ivec_alloc_p(const INT m,ivector **u)
{
  *u=ivec_create_p(m);
  return;
}

/*******************************************************************/
/*!
 * \fn dCSRmat *dcsr_create_plus(const INT m, const INT n, const INT nnz
                                 void *ia, void *ja, void *aij)
 *
 * \brief Create a dCSRmat sparse matrix. Uses void array for the arrays of the matrix
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 * \param ia   Pointer to the row pointer array of a dCSRmat matrix
 * \param ja   Pointer to the column index array of a dCSRmat matrix
 * \param aij  Pointer to the nonzeros array of a dCSRmat matrix
 *
 * \return A a pointer to a dCSRmat matrix. All of this matrix can be
 *         freed by free((void *)A) or even free(A).
 *
 */
dCSRmat *dcsr_create_plus(INT m,		\
			  INT n,		\
			  INT nnz,		\
			  void *ia,		\
			  void *ja,		\
			  void *aij)
{
  dCSRmat *a=malloc(1*sizeof(dCSRmat));// size of the struct
  a->row=m;  a->col=n;  a->nnz=nnz;
  if(ia)
    a->IA=ia;
  else
    a->IA=(INT *)realloc(NULL,(n+1)*sizeof(INT));
  if(ja)
    a->JA=ja;
  else
    a->JA=(INT *)realloc(NULL,nnz*sizeof(INT));
  if(aij)
    a->val=aij;
  else
    a->val=(REAL *)realloc(NULL,nnz*sizeof(REAL));
  return a;
}

/**********************************************************************/
/*!
 * \fn void dcsr_alloc_plus(const INT m, const INT n, const INT nnz,
 *                          void *ia, void *ja, void *aij, dCSRmat **A)
 *
 * \brief Allocate a dCSRmat sparse matrix and put the result in the
 *        pointer pointed by A. Same as dCSRmat_create_p;us, but uses
 *        realloc() to reallocate A.
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 * \param ia   Pointer to the row pointer array of a dCSRmat matrix
 * \param ja   Pointer to the column index array of a dCSRmat matrix
 * \param aij  Pointer to the nonzeros array of a dCSRmat matrix
 *
 * \return A a pointer to a an array of dCSRmat matrices with one
 *         element. All of this matrix can be freed by free((void *)A)
 *         or even just free(A).
 *
 *  \note: call as dCSRmat *X; dcsr_alloc_p(... &X) to reallocate *X
 *         to desired length.
 *         modified by ludmil zikatanov on 20200412
 */
void dcsr_alloc_plus(INT m, INT n, INT nnz, void *ia, void *ja, void *aij, dCSRmat *A)
{
    A = malloc(1 * sizeof(dCSRmat)); // size of the struct
    A->row = m;  A->col = n;  A->nnz = nnz;
    if(ia)
        A->IA = ia;
    else
        A->IA = (INT *)realloc(NULL, (n + 1) * sizeof(INT));
    if(ja)
        A->JA = ja;
    else
        A->JA = (INT *)realloc(NULL, nnz * sizeof(INT));
    if(aij)
        A->val = aij;
    else
        A->val = (REAL *)realloc(NULL, nnz * sizeof(REAL));
}

/******************************************************************/
/*!
 * \fn dvector dvec_create_plus(const INT m, void *vi)
 *
 * \brief Create a dvector of given length: use void array to pack the
 * structure in.
 *
 * \param m    length of the dvector
 * \param vi   pointer to the dvector
 *
 * \return pointer to u   The new dvector
 *
 */
dvector *dvec_create_plus(INT n,		\
			 void *vi)
{
  dvector *v=malloc(1*sizeof(dvector));// size of the struct
  v->row=n;
  if(vi)
    v->val=vi;
  else
    v->val=(REAL *)realloc(NULL,n*sizeof(REAL));
  return v;
 }

 /******************************************************************/
/**
 * \fn dcsr2full(dCSRmat *A, REAL *Afull)
 *
 * \brief converts a dcsr matrix to a full matrix
 * \note
 *
 * \param A	   Matrix A to be converted
 * \param Afull	   Matrix Afull: full(A). it must be allocated
 *                               before entering here.
 * \return
 *
 */
void dcsr2full(dCSRmat *A,REAL *Afull)
{
  // Afull must have enopugh space for (A->row*A->col) doubkles.
  if(Afull){
    INT n=A->row,m=A->col,i,j,jk,im;
    INT *ia=A->IA,*ja=A->JA;
    REAL *a=A->val;
    memset(Afull,0,n*m*sizeof(REAL));
    for(i=0;i<n;i++)  {
      im=i*m;
      for(jk=ia[i];jk<ia[i+1];jk++)  {
	j=ja[jk];
	Afull[im+j]=a[jk];
      }
    }
    return;
  } else {
    fprintf(stderr,"ERROR in %s: Afull not allocated",__FUNCTION__);
    exit(16);
  }
}

/***********************************************************************************************/
/*!
   * \fn void bdcsr_get_total_size(block_dCSRmat *A, INT total_row, INT total_col, INT total_nnz)
   *
   * \brief get total number of rows and columns of a block_dCSRmat matrix
   *
   * \param A          Pointer to the block_dCSRmat matrix
   * \param total_row  total number of rows
   * \param total_col  total number of cols
   * \param total_nnz  total number of nonzeros
   *
   * \author  Xiaozhe Hu
   *
   *
   */
void bdcsr_get_total_size(block_dCSRmat *A, INT *total_row, INT *total_col, INT *total_nnz)
{
    INT i, j;
    *total_row = 0;
    *total_col = 0;
    *total_nnz = 0;

    INT n = A->brow;
    INT m = A->bcol;

    // get total row
    for (i = 0; i < n; ++i){
        for (j = 0; j < m; ++j) {
            if (A->blocks[i*m + j] != NULL) {
                // get the number row for i-th row
                *total_row = *total_row + A->blocks[i*m + j]->row;
            }
            break;
        }
    }

    // get total row
    for (j = 0; j < m; ++j){
        for (i = 0; i < n; ++i) {
            if (A->blocks[i*m + j] != NULL) {
                // get the number columns for j-th column
                *total_col = *total_col + A->blocks[i*m + j]->col;
            }
            break;
        }
    }

    // get total nnz
    for (i = 0; i < n; ++i){
        for (j = 0; j < m; ++j) {
            if (A->blocks[i*m + j] != NULL) {
                // get the number rows of the ij-th block
                *total_nnz = *total_nnz + A->blocks[i*m + j]->nnz;
            }
        }
    }

}

/***********************************************************************************************/
/*!
   * \fn void bdcsr_getdiagblk (block_dCSRmat *A, dCSRmat *A_diag)
   *
   * \brief get diagonal blocks of a block_dCSRmat matrix and stores them in a dCSRmat array
   *
   * \param A       Pointer to the block_dCSRmat matrix
   * \param A_diag  Array of the dCSRmat matrices (allocated!)
   *
   * \author  Xiaozhe Hu
   *
   * \note Assume the dCSRmat *A_diag has been allocated!!
   *
   */
void bdcsr_getdiagblk(block_dCSRmat *A, dCSRmat *A_diag)
{
    INT i;
    INT n = A->brow;

    // allocate array of dCSRmat
    // dCSRmat *A_diag = (dCSRmat*)calloc(A->brow, sizeof(dCSRmat));

    // copy diag blocks
    for (i = 0; i < n; ++i){
        dcsr_alloc(A->blocks[i*n + i]->row, A->blocks[i*n + i]->col, A->blocks[i*n + i]->nnz, &(A_diag[i]));
        dcsr_cp(A->blocks[i*n + i], &(A_diag[i]));
    }

    // return A_diag;
}

/*********************************************************************************/
/*!
 * \fn dvector *bdcsr_getdiag(block_dCSRmat *Ab, const INT n1, const INT n2)
 *
 * \brief   extracts the diagonal entries of Ab(n1:n2,n1:n2) and stored them in a dvector
 *
 * \param Ab    Point to a block_dCSRmat matrix
 * \param n1    staring index for the blocks
 * \param n2    ending index for the blocks
 *
 * \note index starts with 0
 *
 * \author Ludmil Zikatanov & Xiaozhe Hu
 *
 */
dvector *bdcsr_getdiag(block_dCSRmat *Ab,
                              const INT n1,
                              const INT n2)
{

  // local variables
  INT i;
  INT total_size = 0;
  INT nb = Ab->brow;

  // loop 1: get size
  for (i=n1; i<n2; i++)
  {
    total_size = total_size + Ab->blocks[i*(nb+1)]->row;
  }

  // loop 2: get diagonals
  dvector *diag_A = dvec_create_p(total_size); // allocate
  dvector temp_vec;
  INT temp_n;

  // reset total_size
  total_size = 0;

  for (i=n1; i<n2; i++)
  {
     //printf("i=%d\n",i);
     // get size for current block
     temp_n = Ab->blocks[i*(nb+1)]->row;

     // get the diagonal entries of the current block
     dcsr_getdiag(temp_n, Ab->blocks[i*(nb+1)], &temp_vec);

     // copy diagonal entry to the correct place
     array_cp(temp_n, temp_vec.val, diag_A->val+total_size);

     // update total size
     total_size = total_size + temp_n;

     // free temp_vec
     dvec_free(&temp_vec);

  }

  return diag_A;

}

/*********************************************************************************/
/*!
 * \fn dCSRmat *bdcsr_getdiagblk_dcsr(block_dCSRmat *Ab, const INT n10, const INT n20)
 *
 * \brief   get the diagonal blocks Ab(n1:n2,n1:n2), merge them, and store them in a dCSR matrix;
 *
 * \param Ab    Point to a block_dCSRmat matrix
 * \param n10    staring index for the blocks
 * \param n20    ending index for the blocks
 *
 * \note    Memory space for the dCSRmat matrix is allocated inside this function! -- Xiaozhe Hu
 * \note    modeled on bdcsr_2_dcsr from utilities/format.c -- Ludmil
 *
 * \author Ludmil Zikatanov
 */
INT bdcsr_getdiagblk_dcsr(block_dCSRmat *Ab,
                                 const INT n10,
                                 const INT n20,
                                 dCSRmat *A)
{
  // local variables
  INT m=0,n=0,nnz=0;
  const INT mb=Ab->brow, nb=Ab->bcol;//, n_blocks=mb*nb;
  dCSRmat **blockptr=Ab->blocks, *blockptrij;
  INT i,j,ij,ir,i1,length,ilength,start,irmrow,irmrowp1;
  INT *row, *col;
  INT n1=n10,n2=n20;
  if(n10<0) n1 = 0;
  if(n20>mb) n2=mb;
  if(n2<n1) {j=n2;n2=n1;n1=j;}
  // flag for errors
  SHORT status = SUCCESS;
  row = (INT *)calloc(mb+1,sizeof(INT));
  col = (INT *)calloc(nb+1,sizeof(INT));
  // get the size of A
  row[0]=0; col[0]=0;

  // count number of rows
  for (i=n1;i<n2;++i) {
    status = ERROR_BLKMAT_ZERO;
    for (j=n1; j<n2; ++j){
      if (blockptr[i*nb+j]) {
	m+=blockptr[i*nb+j]->row;
	row[i+1]=m;
	status = SUCCESS;
	break;
      }
    }
    // check error
    if (status < SUCCESS) check_error(ERROR_BLKMAT_ZERO, __FUNCTION__);
  }

  // count number of columns
  for (i=n1;i<n2;++i) {
    status = ERROR_BLKMAT_ZERO;
    for (j=n1;j<n2;++j){
      if (blockptr[j*mb+i]) {
	n+=blockptr[j*mb+i]->col;
	col[i+1]=n;
	status = SUCCESS;
	break;
      }
    }
    // check error
    if (status < SUCCESS) check_error(ERROR_BLKMAT_ZERO, __FUNCTION__);
  }
  // count number of nonzeros
  for (i=n1;i<n2;++i) {
    for (j=n1;j<n2;++j){
      if (blockptr[i*mb+j]) {
	nnz+=blockptr[i*mb+j]->nnz;
      }
    }
  }
  // memory space allocation
  //A = dcsr_create_p(m,n,nnz);
  dcsr_alloc(m,n,nnz,A);
  // set dCSRmat for A
  A->IA[0]=0;
  for (i=n1;i<n2;++i) {
    for (ir=row[i];ir<row[i+1];ir++) {
      for (length=j=n1;j<n2;++j) {
	ij=i*nb+j;
	blockptrij=blockptr[ij];
	if (blockptrij && blockptrij->nnz>0) {
	  start=A->IA[ir]+length;
	  irmrow=ir-row[i];irmrowp1=irmrow+1;
	  ilength=blockptrij->IA[irmrowp1]-blockptrij->IA[irmrow];
	  if (ilength>0) {
	    memcpy((A->val+start),(blockptrij->val+blockptrij->IA[irmrow]),ilength*sizeof(REAL));
	    memcpy((A->JA+start),(blockptrij->JA+blockptrij->IA[irmrow]), ilength*sizeof(INT));
	    // shift column index
	    for (i1=0;i1<ilength;i1++) A->JA[start+i1]+=col[j];
	    length+=ilength;
	  }
	}
      } // end for j
      A->IA[ir+1]=A->IA[ir]+length;
    } // end for ir
  } // end for i
  A->nnz=A->IA[row[n2]];
  /* INT nblk=n2-n1+1; // number of blocks */
  /* for(i=n1;i<=n2;i++){ */
  /*   fprintf(stdout,"\nblk=%d,row=%d",i,row[i]); */
  /* }   */
  /* for(i=n1;i<=n2;i++){ */
  /*   fprintf(stdout,"\nblk=%d,row=%d",i,col[i]); */
  /* } */
  /* fprintf(stdout,"\n*** IAend=%d\n",A->IA[row[n2]]); */
  /* fprintf(stdout,"\nA11 data:(%d,%d,%d):rows:(%d,%d)\n",A->row,A->col,A->nnz,row[n2-1],row[n2]); */
  free(row);
  free(col);
  return 0;
}

/***********************************************************************************************/
/*!
 * \fn SHORT dcsr_sparse (dCSRmat *A, ivector *ii, ivector *jj, dvector *kk, INT m, INT n)
 *
 * \brief Form a dCSRmat A with m rows, n columns. Input vectors ii, jj, kk must be of equal lengths.
 *        For each index l, the (ii(l),jj(l))-entry of A is kk(l). In addition, if (ii(l_1),jj(l_1)),
 *        (ii(l_2),jj(l_2)),...,(ii(l_m),jj(l_m)) are identical, then the (ii(l_1),jj(l_1))-entry
 *        of A is the sum of kk(l_1),kk(l_2),...,kk(l_m).
 *        This function emulates the function "sparse" in Octave/Matlab.
 *
 * \param A   Pointer to dCSRmat matrix
 *
 * \return SUCCESS if successful;
 *
 * \note Ludmil, Yuwen 20210606.
 */
SHORT dcsr_sparse (dCSRmat *A, ivector *ii, ivector *jj, dvector *kk, INT m, INT n)
{
      // get size
    INT nnz=kk->row;
    if ((ii->row!=jj->row)||(ii->row!=kk->row)){
      printf("The input vectors are not of equal length!");
      return 0;
    }
    dCOOmat A0;
    A0.row = m;A0.col = n;A0.nnz = nnz;
    A0.rowind = ii->val;A0.colind = jj->val;A0.val = kk->val;
    dcoo_2_dcsr(&A0,A);

    // copy pointers for easier reference
    INT *ia = A->IA;
    INT *ja = A->JA;
    REAL *a = A->val;
    INT i, ij,j,ih,iai;
    SHORT norepeat;
    INT maxdeg=ia[1]-ia[0];

    INT *ind = (INT *) calloc(n,sizeof(INT));
    for (i=0; i<n; ++i) ind[i]=-1;
    // clean up because there might be some repeated indices
    // compute max degree of all vertices (for memory allocation):
    for (i=1;i<m;++i){
      ih=ia[i+1]-ia[i];
      if(maxdeg<ih) maxdeg=ih;
    }
    REAL *atmp=calloc(maxdeg,sizeof(REAL));
    INT *jatmp=calloc(maxdeg,sizeof(INT));
    nnz=0;
    for (i=0;i<m;++i){
      // loop over each row. first find the length of the row:
      ih=ia[i+1]-ia[i];
      // copy the indices in tmp arrays.
      memcpy(jatmp,(ja+ia[i]),ih*sizeof(INT));
      memcpy(atmp,(a+ia[i]),ih*sizeof(REAL));
      norepeat=1;
      for(ij=0;ij<ih;++ij){
	j=jatmp[ij];
	if(ind[j]<0){
	  ind[j]=ij;
	} else {
	  norepeat=0; // we have a repeated index.
	  atmp[ind[j]]+=atmp[ij];
	  jatmp[ij]=-abs(jatmp[ij]+1);
	}
      }
      for(ij=0;ij<ih;++ij){
	j=jatmp[ij];
	if(j<0) continue;// if j is negative, this has repeated somewhere. do nothing
	ind[j]=-1;// make, for all all visited j, ind[j]=-1;
      }
      if(norepeat) continue; // do nothing if no indices repeat.
      // put everything back, but now we have negative column indices
      // on the repeated column indices and we have accumulated the
      // values in the first position of j on every row.
      memcpy(&ja[ia[i]],jatmp,ih*sizeof(INT));
      memcpy(&a[ia[i]],atmp,ih*sizeof(REAL));
    }
    if (ind) free(ind);
    if(atmp) free(atmp);
    if(jatmp) free(jatmp);
    // run over the matrix and remove all negative column indices.
    iai=ia[0];
    nnz=0;
    for (i=0;i<m;++i){
      for(ij=iai;ij<ia[i+1];++ij){
	j=ja[ij];
	if(j<0) continue;
	ja[nnz]=ja[ij];
	a[nnz]=a[ij];
	++nnz;
      }
      iai=ia[i+1];
      ia[i+1]=nnz;
    }
    A->nnz=nnz;
    A->val=realloc(A->val,A->nnz*sizeof(REAL));
    A->JA=realloc(A->JA,A->nnz*sizeof(INT));
    return SUCCESS;
}

/*****************************************************************************/
/*!
 * \fn void icsr_uniqueij(iCSRmat *U, ivector *ii, ivector *jj)
 *
 * \brief Form an iCSRmat U. Input column vectors ii, jj must be of
 *        equal lengths and be index vectors ranging from 0 to
 *        max(ii,jj).  It removes the duplicated (i,j)-pairs from the
 *        matrix in (ii,jj), rearrange the condensed (ii,jj) in
 *        ascending lexicographically order, and returns it in the
 *        icsr format U. This function is also used to obtain edge
 *        structures in a simplicial complax.
 *
 * \param U   Pointer to iCSRmat matrix
 *
 * \note  Ludmil, Yuwen 20210606.
 */
void icsr_uniqueij(iCSRmat *U, ivector *ii, ivector *jj)
{
    // transform to CSR format
    INT i, j, nr = ii->row, nv=-1;
    for (i=0;i<nr;i++) {
      if (nv<ii->val[i]) {nv=ii->val[i];}
      if (nv<jj->val[i]) {nv=jj->val[i];}
    }
    nv = nv + 1;
    INT nnz = ii->row;
    U->IA=realloc(U->IA,(nv+1)*sizeof(INT));
    U->JA=realloc(U->JA,nnz*sizeof(INT));
    if(U->val) {
      free(U->val);U->val=NULL;
    }
    //
    //    INT *ia=U->IA; INT *ja=U->JA;//alias
    INT *rowind = ii->val;
    INT *colind = jj->val;
    INT iind, jind;

    INT *ind0 = (INT *) calloc(nv+1,sizeof(INT));

    for (i=0; i<nnz; ++i) ind0[rowind[i]+1]++;

    U->IA[0] = 0;
    for (i=1; i<=nv; ++i) {
        U->IA[i] = U->IA[i-1]+ind0[i];
        ind0[i] = U->IA[i];
    }

    for (i=0; i<nnz; ++i) {
        iind = rowind[i];
        jind = ind0[iind];
        U->JA[jind] = colind[i];
        ind0[iind] = ++jind;
    }

    if (ind0) free(ind0);

    // delete duplicated entries
    INT ij,ih,iai;
    SHORT norepeat;
    INT maxdeg=U->IA[1]-U->IA[0];

    INT *ind = (INT *) calloc(nv,sizeof(INT));
    for (i=0; i<nv; ++i) ind[i]=-1;
    // clean up because there might be some repeated indices
    // compute max degree of all vertices (for memory allocation):
    for (i=1;i<nv;++i){
      ih=U->IA[i+1]-U->IA[i];
      if(maxdeg<ih) maxdeg=ih;
    }
    INT *jatmp=calloc(maxdeg,sizeof(INT));
    nnz=0;
    for (i=0;i<nv;++i){
      // loop over each row. first find the length of the row:
      ih=U->IA[i+1]-U->IA[i];
      // copy the indices in tmp arrays.
      memcpy(jatmp,(U->JA+U->IA[i]),ih*sizeof(INT));
      norepeat=1;
      for(ij=0;ij<ih;++ij){
      	j=jatmp[ij];
      	if(ind[j]<0){
      	  ind[j]=ij;
      	} else {
    	   norepeat=0; // we have a repeated index.
	       jatmp[ij]=-abs(jatmp[ij]+1);
        	}
      }
      for(ij=0;ij<ih;++ij){
      	j=jatmp[ij];
	      if(j<0) continue;// if j is negative, this has repeated somewhere. do nothing
	      ind[j]=-1;// make, for all all visited j, ind[j]=-1;
      }
      if(norepeat) continue; // do nothing if no indices repeat.
      // put everything back, but now we have negative column indices
      // on the repeated column indices and we have accumulated the
      // values in the first position of j on every row.
      memcpy(&U->JA[U->IA[i]],jatmp,ih*sizeof(INT));
    }
    if (ind) free(ind);
    if(jatmp) free(jatmp);
    // run over the matrix and remove all negative column indices.
    iai=U->IA[0];
    nnz=0;
    for (i=0;i<nv;++i){
      for(ij=iai;ij<U->IA[i+1];++ij){
      	j=U->JA[ij];
      	if(j<0) continue;
      	U->JA[nnz]=U->JA[ij];
	++nnz;
      }
      iai=U->IA[i+1];
      U->IA[i+1]=nnz;
    }
    U->JA=realloc(U->JA,nnz*sizeof(INT));
    U->row = nv; U->col = nv; U->nnz = nnz;
    //
    /* sorting the i-th row of U to get U->JA[U->IA[i]]:U->JA[U->IA[i+1]-1] in
       ascending lexicographic order */
    iCSRmat UT;
    icsr_trans(U,&UT);
    // free is needed because U is already alloc'd, so valgrind complains rightfully that we leak here.
    icsr_free(U);
    icsr_trans(&UT,U);
    icsr_free(&UT);
    /**END sorting**/
    return;
}

/**
 * \fn dBSRmat dbsr_create (const INT ROW, const INT COL, const INT NNZ,
 *                          const INT nb, const INT storage_manner)
 *
 * \brief Create a BSR sparse matrix (allocate memory)
 *
 * \param ROW             Number of rows of block
 * \param COL             Number of columns of block
 * \param NNZ             Number of nonzero blocks
 * \param nb              Dimension of each block
 * \param storage_manner  Storage manner for each sub-block
 *
 * \return A              The new dBSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date   10/26/2010
 */
dBSRmat dbsr_create (const INT  ROW,
                     const INT  COL,
                     const INT  NNZ,
                     const INT  nb,
                     const INT  storage_manner)
{
    dBSRmat A;

    if ( ROW > 0 ) {
        A.IA = (INT*)calloc(ROW+1, sizeof(INT));
    }
    else {
        A.IA = NULL;
    }

    if ( NNZ > 0 ) {
        A.JA = (INT*)calloc(NNZ ,sizeof(INT));
    }
    else {
        A.JA = NULL;
    }

    if ( nb > 0 && NNZ > 0) {
        A.val = (REAL*)calloc(NNZ*nb*nb, sizeof(REAL));
    }
    else {
        A.val = NULL;
    }

    A.storage_manner = storage_manner;
    A.ROW = ROW;
    A.COL = COL;
    A.NNZ = NNZ;
    A.nb  = nb;

    return A;
}

/**
 * \fn void dbsr_alloc (const INT ROW, const INT COL, const INT NNZ,
 *                      const INT nb, const INT storage_manner, dBSRmat *A)
 *
 * \brief Allocate memory space for a BSR format sparse matrix
 *
 * \param ROW             Number of rows of block
 * \param COL             Number of columns of block
 * \param NNZ             Number of nonzero blocks
 * \param nb              Dimension of each block
 * \param storage_manner  Storage manner for each sub-block
 * \param A               Pointer to new dBSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date   10/26/2010
 */
void dbsr_alloc (const INT  ROW,
                 const INT  COL,
                 const INT  NNZ,
                 const INT  nb,
                 const INT  storage_manner,
                 dBSRmat   *A)
{
    if ( ROW > 0 ) {
        A->IA = (INT*)calloc(ROW+1, sizeof(INT));
    }
    else {
        A->IA = NULL;
    }

    if ( NNZ > 0 ) {
        A->JA = (INT*)calloc(NNZ, sizeof(INT));
    }
    else {
        A->JA = NULL;
    }

    if ( nb > 0 ) {
        A->val = (REAL*)calloc(NNZ*nb*nb, sizeof(REAL));
    }
    else {
        A->val = NULL;
    }

    A->storage_manner = storage_manner;
    A->ROW = ROW;
    A->COL = COL;
    A->NNZ = NNZ;
    A->nb  = nb;

    return;
}


/**
 * \fn void dbsr_free (dBSRmat *A)
 *
 * \brief Free memory space for a BSR format sparse matrix
 *
 * \param A   Pointer to the dBSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date   10/26/2010
 */
void dbsr_free (dBSRmat *A)
{
    if (A==NULL) return;

    if (A->IA){
      free(A->IA);
      A->IA  = NULL;
    }
    if (A->JA){
      free(A->JA);
      A->JA  = NULL;
    }
    if (A->val){
      free(A->val);
      A->val = NULL;
    }

    A->ROW = 0;
    A->COL = 0;
    A->NNZ = 0;
    A->nb  = 0;
    A->storage_manner = 0;
}

/**
 * \fn void dbsr_cp (const dBSRmat *A, dBSRmat *B)
 *
 * \brief copy a dCSRmat to a new one B=A
 *
 * \param A   Pointer to the dBSRmat matrix
 * \param B   Pointer to the dBSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date   08/07/2011
 */
void dbsr_cp (const dBSRmat *A,
              dBSRmat       *B)
{
    B->ROW = A->ROW;
    B->COL = A->COL;
    B->NNZ = A->NNZ;
    B->nb  = A->nb;
    B->storage_manner = A->storage_manner;

    memcpy(B->IA,A->IA,(A->ROW+1)*sizeof(INT));
    memcpy(B->JA,A->JA,(A->NNZ)*sizeof(INT));
    memcpy(B->val,A->val,(A->NNZ)*(A->nb)*(A->nb)*sizeof(REAL));
}

/**
 * \fn INT dbsr_trans (const dBSRmat *A, dBSRmat *AT)
 *
 * \brief Find A^T from given dBSRmat matrix A
 *
 * \param A   Pointer to the dBSRmat matrix
 * \param AT  Pointer to the transpose of dBSRmat matrix A
 *
 * \note Xiaozhe Hu (08/06/2011)
 */
INT dbsr_trans (const dBSRmat *A,
                dBSRmat       *AT)
{
    const INT n = A->ROW, m = A->COL, nnz = A->NNZ, nb = A->nb;

    INT status = SUCCESS;
    INT i,j,k,p,inb,jnb,nb2;

    AT->ROW = m;
    AT->COL = n;
    AT->NNZ = nnz;
    AT->nb  = nb;
    AT->storage_manner = A->storage_manner;

    AT->IA  = (INT*)calloc(m+1,sizeof(INT));
    AT->JA  = (INT*)calloc(nnz,sizeof(INT));
    nb2     = nb*nb;

    if (A->val) {
        AT->val = (REAL*)calloc(nnz*nb2,sizeof(REAL));
    }
    else {
        AT->val = NULL;
    }

    // first pass: find the number of nonzeros in the first m-1 columns of A
    // Note: these numbers are stored in the array AT.IA from 1 to m-1
    iarray_set(m+1, AT->IA, 0);

    for ( j=0; j<nnz; ++j ) {
        i=A->JA[j]; // column number of A = row number of A'
        if (i<m-1) AT->IA[i+2]++;
    }

    for ( i=2; i<=m; ++i ) AT->IA[i]+=AT->IA[i-1];

    // second pass: form A'
    if ( A->val ) {
        for ( i=0; i<n; ++i ) {
            INT ibegin=A->IA[i], iend1=A->IA[i+1];
            for ( p=ibegin; p<iend1; p++ ) {
                j=A->JA[p]+1;
                k=AT->IA[j];
                AT->JA[k]=i;
                for ( inb=0; inb<nb; inb++ )
                    for ( jnb=0; jnb<nb; jnb++ )
                        AT->val[nb2*k + inb*nb + jnb] = A->val[nb2*p + jnb*nb + inb];
                AT->IA[j]=k+1;
            } // end for p
        } // end for i

    }
    else {
        for ( i=0; i<n; ++i ) {
            INT ibegin=A->IA[i], iend1=A->IA[i+1];
            for ( p=ibegin; p<iend1; p++ ) {
                j=A->JA[p]+1;
                k=AT->IA[j];
                AT->JA[k]=i;
                AT->IA[j]=k+1;
            } // end for p
        } // end of i

    } // end if

    return status;
}

/**
 * \fn void dbsr_axm (dBSRmat *A, const REAL alpha)
 *
 * \brief Multiply a sparse matrix A in BSR format by a scalar alpha.
 *
 * \param A      Pointer to dBSRmat matrix A
 * \param alpha  REAL factor alpha
 *
 * \author Xiaozhe Hu
 * \date   05/26/2014
 */
void dbsr_axm (dBSRmat     *A,
               const REAL   alpha)
{
    const INT nnz = A->NNZ;
    const INT nb  = A->nb;

    // A direct calculation can be written as:
    array_ax(nnz*nb*nb, alpha, A->val);
}

/*!
 * \fn void dbsr_aAxpby (const REAL alpha, dBSRmat *A,
 *                       REAL *x, const REAL beta, REAL *y)
 *
 * \brief Compute y := alpha*A*x + beta*y
 *
 * \param alpha  REAL factor alpha
 * \param A      Pointer to the dBSRmat matrix
 * \param x      Pointer to the array x
 * \param beta   REAL factor beta
 * \param y      Pointer to the array y
 *
 * \note Works for general nb (Xiaozhe)
 */
void dbsr_aAxpby (const REAL   alpha,
                  dBSRmat     *A,
                  REAL        *x,
                  const REAL   beta,
                  REAL        *y )
{
    /* members of A */
    INT  ROW  = A->ROW;
    INT  nb   = A->nb;
    INT  *IA  = A->IA;
    INT  *JA  = A->JA;
    REAL *val = A->val;

    /* local variables */
    INT     size = ROW*nb;
    INT     jump = nb*nb;
    INT     i,j,k,iend;
    REAL    temp;
    REAL   *pA  = NULL;
    REAL   *px0 = NULL;
    REAL   *py0 = NULL;
    REAL   *py  = NULL;

    //----------------------------------------------
    //   Treat (alpha == 0.0) computation
    //----------------------------------------------
    if (alpha == 0.0) {
        array_ax(size, beta, y);
        return;
    }

    //-------------------------------------------------
    //   y = (beta/alpha)*y
    //-------------------------------------------------
    temp = beta / alpha;
    if (temp != 1.0) {
        if (temp == 0.0) {
            memset(y, 0X0, size*sizeof(REAL));
        }
        else {
            //for (i = size; i--; ) y[i] *= temp; // modified by Xiaozhe, 03/11/2011
            array_ax(size, temp, y);
        }
    }

    //-----------------------------------------------------------------
    //   y += A*x (Core Computation)
    //   each non-zero block elements are stored in row-major order
    //-----------------------------------------------------------------
    for (i = 0; i < ROW; ++i) {
        py0 = &y[i*nb];
        iend = IA[i+1];
        for (k = IA[i]; k < iend; ++k) {
            j = JA[k];
            pA = val+k*jump; // &val[k*jump];
            px0 = x+j*nb; // &x[j*nb];
            py = py0;
            ddense_ypAx( pA, px0, py, nb );
        }
    }

    //------------------------------------------
    //   y = alpha*y
    //------------------------------------------

    if (alpha != 1.0) {
        array_ax(size, alpha, y);
    }
}


/*!
 * \fn void dbsr_aAxpy (const REAL alpha, const dBSRmat *A,
 *                      const REAL *x, REAL *y)
 *
 * \brief Compute y := alpha*A*x + y
 *
 * \param alpha  REAL factor alpha
 * \param A      Pointer to the dBSRmat matrix
 * \param x      Pointer to the array x
 * \param y      Pointer to the array y
 *
 * \note Works for general nb (Xiaozhe)
 */
void dbsr_aAxpy (const REAL      alpha,
                 const dBSRmat  *A,
                 const REAL     *x,
                 REAL           *y)
{
    /* members of A */
    const INT    ROW = A->ROW;
    const INT    nb  = A->nb;
    const INT   *IA  = A->IA;
    const INT   *JA  = A->JA;
    const REAL  *val = A->val;

    /* local variables */
    const REAL *pA   = NULL;
    const REAL *px0  = NULL;
    REAL *py0        = NULL;
    REAL *py         = NULL;

    REAL  temp = 0.0;
    INT   size = ROW*nb;
    INT   jump = nb*nb;
    INT   i, j, k, iend;

    //----------------------------------------------
    //   Treat (alpha == 0.0) computation
    //----------------------------------------------
    if (alpha == 0.0){
        return; // Nothing to compute
    }

    //-------------------------------------------------
    //   y = (1.0/alpha)*y
    //-------------------------------------------------
    if (alpha != 1.0){
        temp = 1.0 / alpha;
        array_ax(size, temp, y);
    }

    //-----------------------------------------------------------------
    //   y += A*x (Core Computation)
    //   each non-zero block elements are stored in row-major order
    //-----------------------------------------------------------------
    for (i = 0; i < ROW; ++i) {
        py0 = &y[i*nb];
        iend = IA[i+1];
        for (k = IA[i]; k < iend; ++k) {
            j = JA[k];
            pA = val+k*jump; // &val[k*jump];
            px0 = x+j*nb; // &x[j*nb];
            py = py0;
            ddense_ypAx( pA, px0, py, nb );
        }
    }

    //------------------------------------------
    //   y = alpha*y
    //------------------------------------------
    if (alpha != 1.0){
        array_ax(size, alpha, y);
    }
    return;
}


/*!
 * \fn void dbsr_aAxpy_agg (const REAL alpha, const dBSRmat *A,
 *                          const REAL *x, REAL *y)
 *
 * \brief Compute y := alpha*A*x + y where each small block matrix is an identity matrix
 *
 * \param alpha  REAL factor alpha
 * \param A      Pointer to the dBSRmat matrix
 * \param x      Pointer to the array x
 * \param y      Pointer to the array y
 *
 * \author Xiaozhe Hu
 * \date   01/02/2014
 *
 * \note Works for general nb (Xiaozhe)
 */
void dbsr_aAxpy_agg (const REAL      alpha,
                     const dBSRmat  *A,
                     REAL     *x,
                     REAL           *y)
{
    /* members of A */
    const INT   ROW = A->ROW;
    const INT   nb  = A->nb;
    const INT  *IA  = A->IA;
    const INT  *JA  = A->JA;

    /* local variables */
    REAL       *px0 = NULL;
    REAL       *py0 = NULL, *py = NULL;
    SHORT       nthreads = 1, use_openmp = FALSE;

    INT         size = ROW*nb;
    INT         i, j, k, iend;
    REAL        temp = 0.0;

    //----------------------------------------------
    //   Treat (alpha == 0.0) computation
    //----------------------------------------------
    if (alpha == 0.0){
        return; // Nothing to compute
    }

    //-------------------------------------------------
    //   y = (1.0/alpha)*y
    //-------------------------------------------------
    if (alpha != 1.0){
        temp = 1.0 / alpha;
        array_ax(size, temp, y);
    }

    //-----------------------------------------------------------------
    //   y += A*x (Core Computation)
    //   each non-zero block elements are stored in row-major order
    //-----------------------------------------------------------------
    for (i = 0; i < ROW; ++i) {
        py0 = &y[i*nb];
        iend = IA[i+1];
        for (k = IA[i]; k < iend; ++k) {
            j = JA[k];
            px0 = x+j*nb; // &x[j*nb];
            py = py0;
            array_axpy(nb, 1.0, px0, py);
        }

    }

    //------------------------------------------
    //   y = alpha*y
    //------------------------------------------
    if ( alpha != 1.0 ) array_ax(size, alpha, y);

    return;
}

/*!
 * \fn void dbsr_mxv (const dBSRmat *A, const REAL *x, REAL *y)
 *
 * \brief Compute y := A*x
 *
 * \param A      Pointer to the dBSRmat matrix
 * \param x      Pointer to the array x
 * \param y      Pointer to the array y
 *
 *
 * \note Works for general nb (Xiaozhe)
 *
 */
void dbsr_mxv (const dBSRmat  *A,
               const REAL     *x,
               REAL           *y)
{
    /* members of A */
    const INT   ROW = A->ROW;
    const INT   nb  = A->nb;
    const INT  *IA  = A->IA;
    const INT  *JA  = A->JA;
    const REAL *val = A->val;

    /* local variables */
    INT     size = ROW*nb;
    INT     jump = nb*nb;
    INT     i,j,k, num_nnz_row;

    const REAL *pA  = NULL;
    const REAL *px0 = NULL;
    REAL       *py0 = NULL;
    REAL        *py  = NULL;


    //-----------------------------------------------------------------
    //  zero out 'y'
    //-----------------------------------------------------------------
    array_set(size, y, 0.0);

    //-----------------------------------------------------------------
    //   y = A*x (Core Computation)
    //   each non-zero block elements are stored in row-major order
    //-----------------------------------------------------------------
    for (i = 0; i < ROW; ++i)
    {
        py0 = &y[i*nb];
        for (k = IA[i]; k < IA[i+1]; ++k)
        {
            j = JA[k];
            pA = val+k*jump; // &val[k*jump];
            px0 = x+j*nb; // &x[j*nb];
            py = py0;
            ddense_ypAx( pA, px0, py, nb );
        }
    }
}


/*!
 * \fn void dbsr_mxv_agg (const dBSRmat *A, REAL *x, REAL *y)
 *
 * \brief Compute y := A*x, where each small block matrices of A is an identity
 *
 * \param A      Pointer to the dBSRmat matrix
 * \param x      Pointer to the array x
 * \param y      Pointer to the array y
 *
 * \author Xiaozhe Hu
 * \date   01/02/2014
 *
 * \note Works for general nb (Xiaozhe)
 */
void dbsr_mxv_agg (const dBSRmat  *A,
                   REAL     *x,
                   REAL           *y)
{
    /* members of A */
    const INT  ROW  = A->ROW;
    const INT  nb   = A->nb;
    const INT  size = ROW*nb;
    const INT *IA   = A->IA;
    const INT *JA   = A->JA;

    /* local variables */
    REAL  *px0 = NULL;
    REAL        *py0 = NULL, *py = NULL;
    INT          i,j,k, num_nnz_row;
    SHORT        use_openmp = FALSE;

    //-----------------------------------------------------------------
    //  zero out 'y'
    //-----------------------------------------------------------------
    array_set(size, y, 0.0);

    //-----------------------------------------------------------------
    //   y = A*x (Core Computation)
    //   each non-zero block elements are stored in row-major order
    //-----------------------------------------------------------------
    for (i = 0; i < ROW; ++i) {
        py0 = &y[i*nb];
        for (k = IA[i]; k < IA[i+1]; ++k) {
            j = JA[k];
            px0 = x+j*nb; // &x[j*nb];
            py = py0;
            array_axpy (nb, 1.0, px0, py);
        }
    }
}


/**
 * \fn void dbsr_mxm (const dBSRmat *A, const dBSRmat *B, dBSRmat *C)
 *
 * \brief Sparse matrix multiplication C=A*B
 *
 * \param A   Pointer to the dBSRmat matrix A
 * \param B   Pointer to the dBSRmat matrix B
 * \param C   Pointer to dBSRmat matrix equal to A*B
 *
 * \author Xiaozhe Hu
 * \date   05/26/2014
 *
 * \note This fct will be replaced! -- Xiaozhe
 */
void dbsr_mxm (const dBSRmat  *A,
               const dBSRmat  *B,
               dBSRmat        *C)
{

    INT i,j,k,l,count;
    INT *JD = (INT *)calloc(B->COL,sizeof(INT));

    const INT nb  = A->nb;
    const INT nb2 = nb*nb;

    // check A and B see if there are compatible for multiplication
    if ( (A->COL != B->ROW) && (A->nb != B->nb ) ) {
        printf("### HAZMATH ERROR: Matrix sizes do not match!\n");
        check_error(ERROR_MAT_SIZE, __FUNCTION__);
    }

    C->ROW = A->ROW;
    C->COL = B->COL;
    C->nb  = A->nb;
    C->storage_manner = A->storage_manner;

    C->val = NULL;
    C->JA  = NULL;
    C->IA  = (INT*)calloc(C->ROW+1,sizeof(INT));

    REAL *temp = (REAL *)calloc(nb2, sizeof(REAL));

    for (i=0;i<B->COL;++i) JD[i]=-1;

    // step 1: Find first the structure IA of C
    for (i=0;i<C->ROW;++i) {
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

    for (i=0;i<C->ROW;++i) C->IA[i+1]+=C->IA[i];

    // step 2: Find the structure JA of C
    INT countJD;

    C->JA=(INT*)calloc(C->IA[C->ROW],sizeof(INT));

    for (i=0;i<C->ROW;++i) {
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
        iarray_set(countJD, JD, -1);
    }

    free(JD); JD = NULL;

    // step 3: Find the structure A of C
    C->val=(REAL*)calloc((C->IA[C->ROW])*nb2,sizeof(REAL));

    for (i=0;i<C->ROW;++i) {
        for (j=C->IA[i];j<C->IA[i+1];++j) {

            array_set(nb2, C->val+(j*nb2), 0x0);

            for (k=A->IA[i];k<A->IA[i+1];++k) {
                for (l=B->IA[A->JA[k]];l<B->IA[A->JA[k]+1];l++) {
                    if (B->JA[l]==C->JA[j]) {
                        ddense_mul (A->val+(k*nb2), B->val+(l*nb2), temp, nb);
                        array_axpy (nb2, 1.0, temp, C->val+(j*nb2));
                    } // end if
                } // end for l
            } // end for k
        } // end for j
    }    // end for i

    C->NNZ = C->IA[C->ROW]-C->IA[0];

    free(temp); temp = NULL;

}


/**
 * \fn void dbsr_rap (const dBSRmat *R, const dBSRmat *A,
 *                              const dBSRmat *P, dBSRmat *B)
 *
 * \brief dBSRmat sparse matrix multiplication B=R*A*P
 *
 * \param R   Pointer to the dBSRmat matrix
 * \param A   Pointer to the dBSRmat matrix
 * \param P   Pointer to the dBSRmat matrix
 * \param B   Pointer to dBSRmat matrix equal to R*A*P (output)
 *
 * \author Xiaozhe Hu
 * \date   10/24/2012
 *
 * \note Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package.
 *            Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 */
void dbsr_rap (const dBSRmat  *R,
               const dBSRmat  *A,
               const dBSRmat  *P,
               dBSRmat        *B)
{
    const INT row=R->ROW, col=P->COL, nb=A->nb, nb2=A->nb*A->nb;

    const REAL *rj=R->val, *aj=A->val, *pj=P->val;
    const INT  *ir=R->IA,  *ia=A->IA,  *ip=P->IA;
    const INT  *jr=R->JA,  *ja=A->JA,  *jp=P->JA;

    REAL       *acj;
    INT        *iac, *jac;

    INT *Ps_marker = NULL;
    INT *As_marker = NULL;

    INT i, i1, i2, i3, jj1, jj2, jj3;
    INT counter, jj_row_begining;

    INT nthreads = 1;

    INT n_coarse = row;
    INT n_fine   = A->ROW;
    INT coarse_mul_nthreads = n_coarse * nthreads;
    INT fine_mul_nthreads = n_fine * nthreads;
    INT coarse_add_nthreads = n_coarse + nthreads;
    INT minus_one_length = coarse_mul_nthreads + fine_mul_nthreads;
    INT total_calloc = minus_one_length + coarse_add_nthreads + nthreads;

    Ps_marker = (INT *)calloc(total_calloc, sizeof(INT));
    As_marker = Ps_marker + coarse_mul_nthreads;

    /*------------------------------------------------------*
     *  First Pass: Determine size of B and set up B_i  *
     *------------------------------------------------------*/
    iac = (INT *)calloc(n_coarse+1, sizeof(INT));

    iarray_set(minus_one_length, Ps_marker, -1);

    REAL *tmp=(REAL *)calloc(2*nthreads*nb2, sizeof(REAL));

    counter = 0;
    for (i = 0; i < row; ++ i) {
        Ps_marker[i] = counter;
        jj_row_begining = counter;
        counter ++;

        for (jj1 = ir[i]; jj1 < ir[i+1]; ++jj1) {
            i1 = jr[jj1];
            for (jj2 = ia[i1]; jj2 < ia[i1+1]; ++jj2) {
                i2 = ja[jj2];
                if (As_marker[i2] != i) {
                    As_marker[i2] = i;
                    for (jj3 = ip[i2]; jj3 < ip[i2+1]; ++jj3) {
                        i3 = jp[jj3];
                        if (Ps_marker[i3] < jj_row_begining) {
                            Ps_marker[i3] = counter;
                            counter ++;
                        }
                    }
                }
            }
        }
        iac[i] = jj_row_begining;
    }

    iac[row] = counter;

    jac=(INT*)calloc(iac[row], sizeof(INT));

    acj=(REAL*)calloc(iac[row]*nb2, sizeof(REAL));

    iarray_set(minus_one_length, Ps_marker, -1);

    /*------------------------------------------------------*
     *  Second Pass: compute entries of B=R*A*P             *
     *------------------------------------------------------*/
     counter = 0;
     for (i = 0; i < row; ++i) {
         Ps_marker[i] = counter;
         jj_row_begining = counter;
         jac[counter] = i;
         array_set(nb2, &acj[counter*nb2], 0x0);
         counter ++;

         for (jj1 = ir[i]; jj1 < ir[i+1]; ++jj1) {
             i1 = jr[jj1];
             for (jj2 = ia[i1]; jj2 < ia[i1+1]; ++jj2) {
                 ddense_mul(&rj[jj1*nb2],&aj[jj2*nb2], tmp, nb);
                 i2 = ja[jj2];
                 if (As_marker[i2] != i) {
                     As_marker[i2] = i;
                     for (jj3 = ip[i2]; jj3 < ip[i2+1]; ++jj3) {
                         i3 = jp[jj3];
                         ddense_mul(tmp, &pj[jj3*nb2], tmp+nb2, nb);
                         if (Ps_marker[i3] < jj_row_begining) {
                             Ps_marker[i3] = counter;
                             array_cp(nb2, tmp+nb2, &acj[counter*nb2]);
                             jac[counter] = i3;
                             counter ++;
                         }
                         else {
                             array_axpy(nb2, 1.0, tmp+nb2, &acj[Ps_marker[i3]*nb2]);
                             }
                         }
                     }
                     else {
                         for (jj3 = ip[i2]; jj3 < ip[i2+1]; jj3 ++) {
                             i3 = jp[jj3];
                             ddense_mul(tmp, &pj[jj3*nb2], tmp+nb2, nb);
                             array_axpy(nb2, 1.0, tmp+nb2, &acj[Ps_marker[i3]*nb2]);
                         }
                     }
                 }
             }
         }

    // setup coarse matrix B
    B->ROW=row; B->COL=col;
    B->IA=iac; B->JA=jac; B->val=acj;
    B->NNZ=B->IA[B->ROW]-B->IA[0];
    B->nb=A->nb;
    B->storage_manner = A->storage_manner;

    free(Ps_marker); Ps_marker = NULL;
    free(tmp);       tmp       = NULL;
}


/**
 * \fn void dbsr_rap_agg (const dBSRmat *R, const dBSRmat *A,
 *                                  const dBSRmat *P, dBSRmat *B)
 *
 * \brief dBSRmat sparse matrix multiplication B=R*A*P, where small block matrices in
 *        P and R are identity matrices!
 *
 * \param R   Pointer to the dBSRmat matrix
 * \param A   Pointer to the dBSRmat matrix
 * \param P   Pointer to the dBSRmat matrix
 * \param B   Pointer to dBSRmat matrix equal to R*A*P (output)
 *
 * \author Xiaozhe Hu
 * \date   10/24/2012
 */
void dbsr_rap_agg (const dBSRmat  *R,
                   const dBSRmat  *A,
                   const dBSRmat  *P,
                   dBSRmat        *B)
{
    const INT row=R->ROW, col=P->COL, nb2=A->nb*A->nb;

    REAL *aj=A->val;
    INT *ir=R->IA, *ia=A->IA, *ip=P->IA;
    INT *jr=R->JA, *ja=A->JA, *jp=P->JA;

    INT  *iac, *jac;
    REAL *acj;
    INT  *Ps_marker = NULL;
    INT  *As_marker = NULL;

    INT i, i1, i2, i3, jj1, jj2, jj3;
    INT counter, jj_row_begining;

    INT nthreads = 1;

    INT n_coarse = row;
    INT n_fine   = A->ROW;
    INT coarse_mul_nthreads = n_coarse * nthreads;
    INT fine_mul_nthreads = n_fine * nthreads;
    INT coarse_add_nthreads = n_coarse + nthreads;
    INT minus_one_length = coarse_mul_nthreads + fine_mul_nthreads;
    INT total_calloc = minus_one_length + coarse_add_nthreads + nthreads;

    Ps_marker = (INT *)calloc(total_calloc, sizeof(INT));
    As_marker = Ps_marker + coarse_mul_nthreads;

    /*------------------------------------------------------*
     *  First Pass: Determine size of B and set up B_i  *
     *------------------------------------------------------*/
    iac = (INT *)calloc(n_coarse+1, sizeof(INT));

    iarray_set(minus_one_length, Ps_marker, -1);

    counter = 0;
    for (i = 0; i < row; ++ i) {
        Ps_marker[i] = counter;
        jj_row_begining = counter;
        counter ++;

        for (jj1 = ir[i]; jj1 < ir[i+1]; ++jj1) {
            i1 = jr[jj1];
            for (jj2 = ia[i1]; jj2 < ia[i1+1]; ++jj2) {
                i2 = ja[jj2];
                if (As_marker[i2] != i) {
                    As_marker[i2] = i;
                    for (jj3 = ip[i2]; jj3 < ip[i2+1]; ++jj3) {
                        i3 = jp[jj3];
                        if (Ps_marker[i3] < jj_row_begining) {
                            Ps_marker[i3] = counter;
                            counter ++;
                        }
                    }
                }
            }
        }
        iac[i] = jj_row_begining;
    }

    iac[row] = counter;

    jac=(INT*)calloc(iac[row], sizeof(INT));

    acj=(REAL*)calloc(iac[row]*nb2, sizeof(REAL));

    iarray_set(minus_one_length, Ps_marker, -1);

    /*------------------------------------------------------*
     *  Second Pass: compute entries of B=R*A*P             *
     *------------------------------------------------------*/
     counter = 0;
     for (i = 0; i < row; ++i) {
         Ps_marker[i] = counter;
         jj_row_begining = counter;
         jac[counter] = i;
         array_set(nb2, &acj[counter*nb2], 0x0);
         counter ++;

         for (jj1 = ir[i]; jj1 < ir[i+1]; ++jj1) {
             i1 = jr[jj1];
             for (jj2 = ia[i1]; jj2 < ia[i1+1]; ++jj2) {

                 i2 = ja[jj2];
                 if (As_marker[i2] != i) {
                     As_marker[i2] = i;
                     for (jj3 = ip[i2]; jj3 < ip[i2+1]; ++jj3) {
                         i3 = jp[jj3];
                         if (Ps_marker[i3] < jj_row_begining) {
                             Ps_marker[i3] = counter;
                             array_cp(nb2, &aj[jj2*nb2], &acj[counter*nb2]);
                             jac[counter] = i3;
                             counter ++;
                         }
                         else {
                             array_axpy(nb2, 1.0, &aj[jj2*nb2], &acj[Ps_marker[i3]*nb2]);
                             }
                         }
                     }
                     else {
                         for (jj3 = ip[i2]; jj3 < ip[i2+1]; jj3 ++) {
                             i3 = jp[jj3];
                             array_axpy(nb2, 1.0, &aj[jj2*nb2], &acj[Ps_marker[i3]*nb2]);
                             }
                         }
                     }
                 }
             }


    // setup coarse matrix B
    B->ROW=row; B->COL=col;
    B->IA=iac; B->JA=jac; B->val=acj;
    B->NNZ=B->IA[B->ROW]-B->IA[0];
    B->nb=A->nb;
    B->storage_manner = A->storage_manner;

    free(Ps_marker); Ps_marker = NULL;
}


/**
 * \fn dvector dbsr_getdiaginv(const dBSRmat *A)
 *
 * \brief Get D^{-1} of matrix A
 *
 * \param A   Pointer to the dBSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date   02/19/2013
 *
 * \note Works for general nb (Xiaozhe)
 */
dvector dbsr_getdiaginv(const dBSRmat *A)
{
    // members of A
    const INT     ROW = A->ROW;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    size = ROW*nb2;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;
    REAL         *val = A->val;

    dvector diaginv;

    INT i,k;

    // Variables for OpenMP
    SHORT nthreads = 1;

    // allocate memory
    diaginv.row = size;
    diaginv.val = (REAL *)calloc(size, sizeof(REAL));

    // get all the diagonal sub-blocks
    for (i = 0; i < ROW; ++i) {
        for (k = IA[i]; k < IA[i+1]; ++k) {
            if (JA[k] == i)
                memcpy(diaginv.val+i*nb2, val+k*nb2, nb2*sizeof(REAL));
        }
    }

    // compute the inverses of all the diagonal sub-blocks
    if (nb > 1) {
        for (i = 0; i < ROW; ++i) {
            ddense_inv_inplace(diaginv.val+i*nb2, nb);
        }
    }
    else {
        for (i = 0; i < ROW; ++i) {
            // zero-diagonal should be tested previously
            diaginv.val[i] = 1.0 / diaginv.val[i];
        }
    }

    return (diaginv);
}

/**
 * \fn static dCSRmat condenseBSR (const dBSRmat *A)
 *
 * \brief Form a dCSRmat matrix from a dBSRmat matrix: use the (1,1)-entry
 *
 * \param A    Pointer to the BSR format matrix
 *
 * \return     dCSRmat matrix if succeed, NULL if fail
 *
 * \author Xiaozhe Hu
 * \date   03/16/2012
 */
dCSRmat condenseBSR (const dBSRmat *A)
{
    // information about A
    const INT   ROW = A->ROW;
    const INT   COL = A->COL;
    const INT   NNZ = A->NNZ;
    const SHORT  nc = A->nb;
    const SHORT nc2 = nc*nc;
    const REAL  TOL = 1e-8;

    const REAL *val = A->val;
    const INT  *IA  = A->IA;
    const INT  *JA  = A->JA;

    // (1,1) block
    dCSRmat  P_csr = dcsr_create(ROW, COL, NNZ);
    REAL    *Pval  = P_csr.val;
    memcpy (P_csr.JA, JA, NNZ*sizeof(INT));
    memcpy (P_csr.IA, IA, (ROW+1)*sizeof(INT));

    INT i, j;

    for ( i=NNZ, j=NNZ*nc2-nc2 + (0*nc+0); i--; j-=nc2 ) Pval[i] = val[j];

    // compress CSR format
    dcsr_compress_inplace (&P_csr,TOL);

    // return P
    return P_csr;
}

/**
 * \fn static dCSRmat condenseBSRLinf (const dBSRmat *A)
 *
 * \brief Form a dCSRmat matrix from a dBSRmat matrix: use inf-norm of each block
 *
 * \param A    Pointer to the BSR format matrix
 *
 * \return     dCSRmat matrix if succeed, NULL if fail
 *
 * \author Xiaozhe Hu
 * \date   05/25/2014
 */
dCSRmat condenseBSRLinf (const dBSRmat *A)
{
    // information about A
    const INT   ROW = A->ROW;
    const INT   COL = A->COL;
    const INT   NNZ = A->NNZ;
    const SHORT  nc = A->nb;
    const SHORT nc2 = nc*nc;
    const REAL  TOL = 1e-8;

    const REAL *val = A->val;
    const INT  *IA  = A->IA;
    const INT  *JA  = A->JA;

    // CSR matrix
    dCSRmat  Acsr = dcsr_create (ROW, COL, NNZ);
    REAL    *Aval = Acsr.val;

    // get structure
    memcpy (Acsr.JA, JA, NNZ*sizeof(INT));
    memcpy (Acsr.IA, IA, (ROW+1)*sizeof(INT));

    INT i, j, k;
    INT row_start, row_end;

    for ( i=0; i<ROW; i++ ) {

        row_start = A->IA[i]; row_end = A->IA[i+1];

        for ( k = row_start; k < row_end; k++ ) {
            j = A->JA[k];
            Aval[k] = ddense_Linf (val+k*nc2, nc);
            if ( i != j ) Aval[k] = -Aval[k];
        }

    }

    // compress CSR format
    dcsr_compress_inplace(&Acsr,TOL);

    // return CSR matrix
    return Acsr;
}

/*********************************EOF***********************************/
