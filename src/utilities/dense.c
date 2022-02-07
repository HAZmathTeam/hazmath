/*! \file src/utilities/dense.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 5/13/15.
 *  Copyright 2016__HAZMATH__. All rights reserved.
 *
 *  \note dense matrix routines: allocation, multiplication, routines
 *                               to solve Ax=b, calculate inv(A),
 *                               transpose, axpy, qr, and other.
 *
 *  \note: Created by Xiaozhe Hu on 01/04/20.
 *  \note: modified by Xiaozhe Hu on 10/25/2016
 *  \note: modified on 20171012 (Ludmil)
 *  \note: modified on 20210807 (Ludmil)
 *
 *  \note: done cleanup for releasing -- Xiaozhe Hu 08/28/2021
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

/*********************************************************************
 *
 *  other dense matrix routines
 *
 **********************************************************************/
/*!
 * \fn ddense_solve_pivot(INT dopivot, INT n, REAL *A, REAL *b,
 *                        INT *p,REAL *piv)
 *
 *
 * \brief Solution of a linear system A*x = b with scaled partial
 *        pivoting. The right hand side is overwritten on output with
 *        the solution, that is x[] is the same as b[] on return.
 *
 * \param dopivot switch [0/1]; if idopivot=1 then do the LU
 *                decomposition of A. The L and U factors are stored
 *                in the lower and upper triangle of A respectively.
 *                if dopivot=0, then A is assumed to be already LU
 *                decomposed and this only performs the forward and
 *                backward substitutions.
 *
 * \param n    number of rows (and columns) of A.
 * \param A    the matrix as a one dimensional array by rows
 *
 */
INT ddense_solve_pivot(INT dopivot, INT n, REAL *A, REAL *b, INT *p,REAL *piv)
{
  INT nm1,i1,k1,pin,kswp,kp,i,j,k;
  REAL r,t,absaij;
  REAL *x=piv;
  if(dopivot) {
    for (i=0;i<n;i++){
      p[i]=i;
      piv[i] = fabs(A[i*n+0]);
      for (j=1;j<n;j++){
      absaij=fabs(A[i*n+j]);
      if(absaij > piv[i])
	piv[i] = absaij;
      }
      piv[i]=1./piv[i]; //here we need error stop if too small
    }
    nm1 = n-1;
    for (k = 0;k<nm1;k++){
      r = fabs(A[p[k]*n+k])*piv[p[k]];
      kp = k;
      for (i=k;i<n;i++){
	t = fabs(A[p[i]*n+k])*piv[p[i]];
	if (t > r){r = t; kp = i;}
      }
      kswp = p[kp]; p[kp] = p[k]; p[k] = kswp;
      k1 = k+1;
      for (i = k1;i<n;i++){
	pin=p[i]*n;
	A[pin+k] = A[pin+k]/A[p[k]*n+k];
	for(j = k1;j<n;j++){
	  A[pin+j] = A[pin+j]-A[pin+k]*A[p[k]*n+j];
	}
      }
    }
    /*end of decomposition; now solver part*/
  }
  x[0] = b[p[0]];
  for (i = 1;i<n;i++){
    x[i] = b[p[i]];
    for(j = 0;j<i;j++){
      x[i] = x[i]-A[p[i]*n+j]*x[j];
    }
  }
  nm1=n-1;
  x[nm1] = x[nm1]/A[p[nm1]*n+nm1];
  for (i = nm1;i>0;i--){
    i1=i-1;
    pin=p[i1]*n;
    for (j = i;j<n;j++){
      x[i1] = x[i1]-A[pin+j]*x[j];
    }
    x[i1] = x[i1]/A[pin+i1];
  }
  //  if(y) free(y);
  for(k=0;k<n;k++)b[k]=x[k];
  return 0;
}

/**************************************************************************/
/*
 * \fn SHORT ddense_lu(INT n, REAL *deta, REAL *A,INT *p,REAL *piv)
 *
 * \brief LU decomposition of an (nxn) matrix A. on output returns the
 *        value of the determinant of A as well. Uses scaled partial
 *        pivoting. The return is non-zero (1 or 2) if a diagonal element
 *        is very small.
 *
 * \param n    Number of rows and columns in A.
 * \param deta is the determinant of A
 * \param A    the matrix as a one dimensional array by rows
 * \param dopivot flag to indicate whether A is already decomposed in
 *                LU and we need to compute the determinant only;
 *
 */
SHORT ddense_lu(INT dopivot, INT n, REAL *deta, REAL *A,INT *p,REAL *piv)
{
  // return 0 if all is OK and returns 1 or 2 if there is a division by a very small number...<tol
  INT nm1,i1,k1,pin,kswp,kp,i,j,k;
  REAL det0,r,t,absaij,tol=1e-10;
  if(dopivot){
    for (i=0;i<n;i++){
      p[i]=i;
      piv[i] = fabs(A[i*n+0]);
      for (j=1;j<n;j++){
	absaij=fabs(A[i*n+j]);
	if(absaij > piv[i]) piv[i] = absaij;
      }
      if(fabs(piv[i])<tol) {
	return (SHORT )1;
      }
      piv[i]=1./piv[i]; //here we need error stop if too small
      //       fprintf(stderr,"\n*** i=%i; pivot=%g\n",i,piv[i]);
    }
    nm1 = n-1;
    for (k = 0;k<nm1;k++){
      r = fabs(A[p[k]*n+k])*piv[p[k]];
      kp = k;
      for (i=k;i<n;i++){
 	t = fabs(A[p[i]*n+k])*piv[p[i]];
 	if (t > r){r = t; kp = i;}
      }
      kswp = p[kp]; p[kp] = p[k]; p[k] = kswp;
      k1 = k+1;
      for (i = k1;i<n;i++){
 	pin=p[i]*n;
	if(fabs(A[p[k]*n+k])<tol) {
	  *deta=0e0;
	  return (SHORT )2;
	}
	A[pin+k] = A[pin+k]/A[p[k]*n+k];
	for(j = k1;j<n;j++){
	  A[pin+j] = A[pin+j]-A[pin+k]*A[p[k]*n+j];
	}
      }
    }
  }
  nm1=n-1;
  det0=A[p[nm1]*n+nm1];
  for (i = nm1;i>0;i--){
    i1=i-1;
    pin=p[i1]*n;
    det0 *= A[pin+i1];
  }
  *deta=det0;
  return (SHORT )0;
}
/**************************************************************************/
/*
 * \fn void ddense_inv(INT n, REAL *Ainv, REAL *A, void *wrk)
 *
 * \brief Inverse of a general (nxn) matrix A
 *
 * \param n    Number of rows
 * \param m    Number of columns
 * \param A    the matrix as a one dim. array by rows
 * \param Ainv the output inverse.
 *
 * \param wrk working array of size at least
 *            n*sizeof(INT)+n*(n+1)*sizeof(REAL)
 *
 */
void ddense_inv(REAL *Ainv,INT n, REAL *A, void *wrk)
{
  REAL *piv=(REAL *)wrk;  REAL *Awrk=piv+n;
  INT *p=(INT *)(wrk+n*(n+1)*sizeof(REAL));
  INT i,j,ji,ni;
  /* first transpose A beacuse we invert by rows;*/
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      ji=(n*j+i);
      Awrk[ji]=*A;
      A++;
    }
  }
  for(i=0;i<n;i++){
    ni=n*i;
    for(j=0;j<n;j++)
      Ainv[ni+j]=0e0;
    Ainv[ni+i]=1e0;
    // solve with rhs a row of the identity;
    ddense_solve_pivot((!i),n, Awrk, (Ainv+ni), p, piv);
  }
  return;
}
/**************************************************************************/
/*
 * \fn void ddense_abyb(const INT m, const INT p, REAL *c, REAL *a,
 *                   REAL *b, const INT n)
 *
 * \brief Computes c = a*b+c; for matrices
 *
 * \param m is the number of row in a and c
 * \param p is the number of columns in b and c
 * \param a is an m by n matrix
 * \param b is an n by p matrix
 * \param c is an m by p matrix.
 * \param n is the number of columns in a and number of rows in b.
 *
 */
void abybfull(const INT m, const INT p, REAL *c,	\
	      REAL *a, REAL *b, const INT n)
{
  /* matrices c = a*b+c; a is m by n, b is n by p; c is m by p */
  REAL cij;
  INT i,j,k,ij,ik,kj,in,ip;
  for(i=0;i<m;i++){
    in=i*n;
    ip=i*p;
    for(j=0;j<p;j++){
      cij=0.;
      for(k=0;k<n;k++){
	ik=in + k;
	kj=k*p + j;
	cij+=a[ik]*b[kj];
      }
      ij=ip+j;
      c[ij]+=cij;
    }
  }
  return;
}
/**************************************************************************/
/*
 * \fn void ddense_abyv(const INT m, REAL *y,REAL *a, REAL *x, const INT n)
 *
 * \brief Computes y = a*x+y; where x and y are vectors and a is a martix.
 *
 * \param m is the number of rows in a and the number of elements in x.
 * \param a is an m by n matrix
 * \param x is an n by 1 vector
 * \param y is an m by 1 vector
 * \param n is the number of columns in a and number of elements in x
 *
 */
void ddense_abyv(const INT m, REAL *y,REAL *a, REAL *x, const INT n)
{
  /* */
  REAL yi;
  INT i,j,in;
  for(i=0;i<m;i++){
    in=i*n; yi=0.;
    for(j=0;j<n;j++){
      yi+=a[in+j]*x[j];
    }
    y[i]+=yi;
  }
  return;
}
/**************************************************************************/
/*
 * \fn ddense_atbyv(const INT m, REAL *y,REAL *a, REAL *x, const INT n)
 *
 * \brief Computes y = transpose(a)*x+y;
 *
 * \param m is the number of rows in a (columns in transpose(a)).
 * \param a is an matrix m by n,
 * \param x is an m by 1 vector,
 * \param y is an n by 1 vector.
 *
 */
void ddense_atbyv(const INT m, REAL *y,REAL *a, REAL *x, const INT n)
{
  INT i,j,in;
  for(i=0;i<m;i++){
    in=i*n;
    for(j=0;j<n;j++){
      y[j]+=x[i]*a[in+j];
    }
  }
  return;
}
/**************************************************************************/

/**************************************************************************/
/*
 * \fn ddense_qr(const INT m, const INT n, REAL *A, REAL *Q, REAL *R)
 *
 * \brief QR decomposition of A: A=QR
 *
 * \param m     number of rows
 * \param n     number of columns
 * \param A     full matrix
 *
 * \return Q    orthogonal full matrix from QR decomposition
 * \return R    upper triangular matrix from QR decomposition
 */
void ddense_qr(const INT m, const INT n, REAL *A, REAL *Q, REAL *R)
{

  // local variables
  INT i,j;
  REAL *qj = (REAL *)calloc(m, sizeof(REAL));
  REAL *aj;

  // main loop
  for (j=0; j<n; j++)
  {

      // qj = aj
      aj = &A[j*m];
      array_cp(m, aj, qj);

      // inner loop
      for (i=0; i<j; i++)
      {

        // rij = qi'*qj
        R[i*n+j] = array_dotprod(m, &Q[i*m], qj);

        // qj = qj - rij*qi
        array_axpy(m, -R[i*n+j], &Q[i*m], qj);

      }

      // rjj = ||qj||_2
      R[j*n+j] = array_norm2(m, qj);

      // qj = qj/rjj
      array_cp(m, qj, &Q[j*m]);
      array_ax(m, 1./R[j*n+j], &Q[j*m]);

  }

  // free
  free (qj);

}

/**************************************************************************/
/*
 * \fn INT ddense_svd(INT m,INT n, REAL *A, REAL *U, REAL *VT, REAL* S,INT computeUV)
 *
 * \brief SVD decomposition of A in Full Format using LAPACK: A=U*S*VT
 *
 * \param m     number of rows
 * \param n     number of columns
 * \param A     full matrix (mxn)
 * \param computeUV 1 returns U and V^T matrices 0: only singular values
 *
 * \return U    orthogonal full matrix from SVD (mxm)
 * \return VT   orthogonal full matrix, V^T from SVD (nxn)
 * \return S    singular values (mxn diagonal)
 * \return      error message
 *
 * \note This calls LAPACK
 */
INT ddense_svd(INT m, INT n, REAL *A, REAL *U, REAL *VT, REAL* S,INT computeUV)
{
  INT info=-22;
#if WITH_LAPACK
  // We assume we are entering with a row-major full Matrix
  // but will convert for FORTRAN
  INT lda = m;
  INT ldu = m;
  INT ldvt = n;

  // Allocate a work array
  REAL* dwork;
  REAL workopt; // Dummy variable to figure out optimal workspace
  INT lwork=-1; // Indicates we're not solving yet

  // Convert to column-wise ordering
  r2c(m,n,sizeof(REAL),(void *)A);

  // Default to not compute U and V^T
  // Later can decide if you want the partial computations.
  char jobuv = 'N';
  if(computeUV) jobuv = 'A';

  // First call to allocate OPTIMAL workspace
  dgesvd_(&jobuv,&jobuv,&m,&n,A,&lda,S,U,&ldu,VT,&ldvt,&workopt,&lwork,&info);
  lwork = (INT) workopt;
  dwork = (REAL *) calloc(lwork,sizeof(REAL));

  // Compute svd
  dgesvd_(&jobuv,&jobuv,&m,&n,A,&lda,S,U,&ldu,VT,&ldvt,dwork,&lwork,&info);

  if(info){
    fprintf(stderr,"\nXXX: lapack error info during svd computations=%d\n",info);fflush(stderr);
    exit(16);
  }
  return info;
#else
  error_extlib(252, __FUNCTION__, "LAPACK");
  return info;
#endif
}

/**************************************************************************/
/*
 * \fn INT ddense_qr_lapack(INT m,INT n, REAL *A, REAL *Q, REAL *R,INT computeR)
 *
 * \brief QR decomposition of A in Full Format using LAPACK: A=QR
 *
 * \param m     number of rows
 * \param n     number of columns
 * \param A     full matrix (mxn)
 * \param computeR 1 returns R 0: only Q
 *
 * \return A   orthogonal full matrix, Q
 * \return R   upper triangular portion if requested
 * \return     error flag
 *
 * \note This calls LAPACK and A is destroyed!
 */
INT ddense_qr_lapack(INT m,INT n, REAL *A,REAL * Q,REAL *R,INT computeR)
{
  INT info=-22;
#if WITH_LAPACK

  // Loop counters
  INT i,j;

  // We assume we are entering with a row-major full Matrix
  // but will convert for FORTRAN
  INT lda = m;

  // Allocate a work array
  REAL* dwork;
  REAL workopt; // Dummy variable to figure out optimal workspace
  INT lwork=-1; // Indicates we're not solving yet

  // Convert to column-wise ordering
  r2c(m,n,sizeof(REAL),(void *)A);

  // Store reflectors to compute Q later
  REAL* tau = (REAL *) calloc(n,sizeof(REAL));

  // First call to allocate OPTIMAL workspace
  dgeqrf_(&m,&n,A,&lda,tau,&workopt,&lwork,&info);
  lwork = (INT) workopt;
  dwork = (REAL *) calloc(lwork,sizeof(REAL));

  // Compute QR
  dgeqrf_(&m,&n,A,&lda,tau,dwork,&lwork,&info);

  // Extract R if needed (back in row-major form)
  if(computeR) {
    for(i=0;i<m;i++) {
      for(j=i;j<n;j++) {
        R[i*n+j] = A[j*m+i];
      }
    }
  }

  // Now we need to compute Q from matrix
  // Form the leading n columns of Q explicitly
  dorgqr_(&m,&n,&n,A,&lda,tau,dwork,&lwork,&info);
  array_cp(m*n,A,Q);
  // Convert back to row-major form
  c2r(m,n,sizeof(REAL),(void *)Q);

  if(info){
    fprintf(stderr,"\nXXX: lapack error info during qr computations=%d\n",info);fflush(stderr);
    exit(16);
  }

  // Frees
  if(tau) free(tau);
  if(dwork) free(dwork);

  return info;

#else
  error_extlib(252, __FUNCTION__, "LAPACK");
  return info;
#endif
}

/*******************************************************************/
/* \fn  void c2r(const INT n, const INT m, const size_t sizeel, void *x)
 *
 * \note: converts 2d array stored by columns to a 2d array stored by
 *     rows (fortran to c); overwrites the input array x with the
 *     result. sizeel is the size of an element in the array in bytes.
 *
*/
void c2r(const INT n, const INT m, const size_t sizeel, void *x)
{
  INT i,j,ji,nms=n*m*sizeel;
  void *y=(void *)malloc(nms);
  memcpy(y,x,nms);
  for (i=0;i<n;i++){
    for (j=0;j<m;j++){
      ji=sizeel*(n*j+i);
      memcpy(x,(y+ji),sizeel);
      x+=sizeel;
    }
  }
  if(y) free(y);
  return;
}
/*******************************************************************/
/* \fn void r2c(const INT n, const INT m, const size_t sizeel, void *x)
 *
 * \note: converts 2d array stored by rows to a 2d array stored by
 *        columns (c to fortran). sizeel is the size of an element in
 *        the array in bytes
 *
*/
void r2c(const INT n, const INT m, const size_t sizeel, void *x)
{
  INT i,j,ji,nms=n*m*sizeel;
  void *y=(void *)malloc(nms);
  for (i=0;i<n;i++){
    for (j=0;j<m;j++){
      ji=sizeel*(n*j+i);
      memcpy((y+ji),x,sizeel);
      x+=sizeel;
    }
  }
  memcpy((x-nms),y,nms);
  if(y) free(y);
  return;
}
/*******************************************************************/
/*! \fn INT ddense_solve_pivot_l(INT dopivot,INT n,REAL16 *A,REAL16 *b,
 *			        INT *p, REAL16 *piv)
 *
 * \brief LONG_DOUBLE: Solution of a linear system A*x = b with scaled
 *        partial pivoting in LONG double. The right hand side is
 *        overwritten on output with the solution, that is x[] is the
 *        same as b[] on return.
 *
 * \param dopivot switch [0/1]; if idopivot=1 then do the LU
 *                 decomposition of A. The L and U factors are stored
 *                 in the lower and upper triangle of A respectively.
 *                 if dopivot=0, then A is assumed to be already LU
 *                 decomposed and this only performs the forward and
 *                 backward substitutions.
 *
 * \param n    number of rows (and columns) of A.
 *
 * \param A the matrix as a one dimensional array by rows (long
 *          double)
 *
 */
INT ddense_solve_pivot_l(INT dopivot,			\
			 INT n,				\
			 REAL16 *A,			\
			 REAL16 *b,			\
			 INT *p,			\
			 REAL16 *piv)
{
  INT nm1,i1,k1,pin,kswp,kp,i,j,k;
  REAL16 r,t,absaij;
  REAL16 *x=piv;
  if(dopivot){
    for (i=0;i<n;i++){
      p[i]=i;
      piv[i] = fabs(A[i*n+0]);
      for (j=1;j<n;j++){
	absaij=fabs(A[i*n+j]);
	if(absaij > piv[i])
	  piv[i] = absaij;
      }
      piv[i]=1./piv[i]; //here we need error stop if too small
    }
    nm1 = n-1;
    for (k = 0;k<nm1;k++){
      r = fabs(A[p[k]*n+k])*piv[p[k]];
      kp = k;
      for (i=k;i<n;i++){
	t = fabs(A[p[i]*n+k])*piv[p[i]];
	if (t > r){r = t; kp = i;}
      }
      kswp = p[kp]; p[kp] = p[k]; p[k] = kswp;
      k1 = k+1;
      for (i = k1;i<n;i++){
	pin=p[i]*n;
	A[pin+k] = A[pin+k]/A[p[k]*n+k];
	for(j = k1;j<n;j++){
	  A[pin+j] = A[pin+j]-A[pin+k]*A[p[k]*n+j];
	}
      }
    }
  }
  /*end of decomposition; now solver part*/
  x[0] = b[p[0]];
  for (i = 1;i<n;i++){
    x[i] = b[p[i]];
    for(j = 0;j<i;j++){
      x[i] = x[i]-A[p[i]*n+j]*x[j];
    }
  }
  nm1=n-1;
  x[nm1] = x[nm1]/A[p[nm1]*n+nm1];
  for (i = nm1;i>0;i--){
    i1=i-1;
    pin=p[i1]*n;
    for (j = i;j<n;j++){
      x[i1] = x[i1]-A[pin+j]*x[j];
    }
    x[i1] = x[i1]/A[pin+i1];
  }
  //  if(y) free(y);
  for(k=0;k<n;k++)b[k]=x[k];
  return 0;
}
/*******************************************************************/
/*! \fn INT zdense_solve_pivot_l(INT n,REAL16 *Ar, REAL16 *Ai,
 *                              REAL16 *br,REAL16 *bi)
 *
 * \brief COMPLEX LONG_DOUBLE: Solution of a linear system
 *        (Ar+i*Ai)*(xr+i*xi)=br+i*bi with scaled partial pivoting in
 *        LONG double, where Ar and Ai are the real and imaginary part
 *        of a matrix; br and bi are the real and imaginary part of
 *        the rhs. The right hand side is overwritten on output with
 *        the solution, that is, on return: xr[],xi[] occuppy the same
 *        memory as br[],bi[].
 *
 * \param n number of rows (and columns) of Ar and Ai.
 *
 * \param bi[] imaginary part of rhs cannot be null if Ai is not null;
 *             in such case bi[] must be allocated before entry here
 *             with n*sizeof(REAL16)
 *
 * \param both Ar an Ai the matrix as one dimensional arrays each
 *        stored by rows (long double)
 *
 */
INT zdense_solve_pivot_l(INT n,				\
			 REAL16 *Ar,			\
			 REAL16 *Ai,			\
			 REAL16 *br,			\
			 REAL16 *bi)
{
  INT k,j,n2=-10;
  REAL16 *A=NULL;
  REAL16 *b=NULL;
  REAL16 *piv=NULL;
  INT *perm=NULL;
  // Ai is null or bi is null; or both;
  if((Ai==NULL) && (bi==NULL)){
    piv=calloc(n,sizeof(REAL16));
    perm=calloc(n,sizeof(INT));
    ddense_solve_pivot_l(1,n,Ar,br,perm,piv);
    free(piv);
    free(perm);
    return 0;
  } else if(Ai==NULL){// bi is not null;
    piv=calloc(n,sizeof(REAL16));
    perm=calloc(n,sizeof(INT));
    ddense_solve_pivot_l(1,n,Ar,br,perm,piv);
    ddense_solve_pivot_l(1,n,Ar,bi,perm,piv);
    free(piv);
    free(perm);
    return 0;
  } else if(bi==NULL) {// this means Ai is not null and all is coupled, but bi is null, so we do an error solution stop:
    fprintf(stderr,"%%%% ******* ERROR: Ai[][] is not null but bi[] is null; please allocate bi[] before call to %s",__FUNCTION__);
    exit(4);
  } else { //neither Ai nor bi is null: compute
    n2=2*n;
    A=calloc(n2*n2,sizeof(REAL16));
    b=calloc(n2,sizeof(REAL16));
    memset(A,0,n2*n2*sizeof(REAL16));
    memset(b,0,n2*sizeof(REAL16));
    for (j=0;j<n;j++){
      for (k=0;k<n;k++){
	A[j*n2  + k]             =  Ar[j*n + k];//A11;
	A[j*n2  + k + n]         = -Ai[j*n + k];//A12;
	A[(j + n)*n2 + k]       =  Ai[j*n + k];//A21;
	A[(j + n)*n2  +  k + n] =  Ar[j*n + k];//A22;
      }
      b[j]     = br[j];
      b[j + n] = bi[j];
    }
    /* print_full_mat_l(n2,n2,A,"A"); */
    /* print_full_mat_l(n2,1,b,"b"); */
    piv=calloc(n2,sizeof(REAL16));
    perm=calloc(n2,sizeof(INT));
    ddense_solve_pivot_l(1,n2,A,b,perm,piv);
    //  if(y) free(y);
    for(k=0;k<n;k++) {
      br[k]=b[k];
      bi[k]=b[k+n];
    }
    free(perm);
    free(piv);
    free(A);
    free(b);
    return 0;
  }
}

/**
 * \fn void ddense_axm (REAL *a, const INT n, const REAL alpha)
 *
 * \brief Compute a = alpha*a (in place)
 *
 * \param a        Pointer to the REAL array which stands a n*n matrix
 * \param n        Dimension of the matrix
 * \param alpha    Scalar
 *
 * \author Xiaozhe Hu
 * \date   05/26/2014
 */
void ddense_axm (REAL       *a,
                 const INT   n,
                 const REAL  alpha)
{
    const INT  n2 = n*n;
    INT        i;

    for ( i = 0; i < n2; i++ ) a[i] *= alpha;

    return;
}

/**
 * \fn void ddense_add (const REAL *a, const REAL *b, const INT n,
 *                      const REAL alpha, const REAL beta, REAL *c)
 *
 * \brief Compute c = alpha*a + beta*b
 *
 * \param a        Pointer to the REAL array which stands a n*n matrix
 * \param b        Pointer to the REAL array which stands a n*n matrix
 * \param n        Dimension of the matrix
 * \param alpha    Scalar
 * \param beta     Scalar
 * \param c        Pointer to the REAL array which stands a n*n matrix
 *
 * \author Xiaozhe Hu
 * \date   05/26/2014
 */
void ddense_add (const REAL  *a,
                 const REAL  *b,
                 const INT    n,
                 const REAL   alpha,
                 const REAL   beta,
                 REAL        *c)
{
    const INT  n2 = n*n;
    INT        i;

    for ( i = 0; i < n2; i++ ) c[i] = alpha * a[i] + beta * b[i];

    return;
}


/**
 * \fn void ddense_mxv (const REAL *a, const REAL *b, REAL *c, const INT n)
 *
 * \brief Compute the product of a small full matrix a and a array b, stored in c
 *
 * \param a   Pointer to the REAL array which stands a n*n matrix
 * \param b   Pointer to the REAL array with length n
 * \param c   Pointer to the REAL array with length n
 * \param n   Dimension of the matrix
 *
 * \author Xiaozhe Hu
 * \date   04/21/2010
 */
void ddense_mxv (const REAL  *a,
                 const REAL  *b,
                 REAL        *c,
                 const INT    n)
{
    INT i,j,in=0;
    REAL temp;

    for (i=0; i<n; ++i, in+=n) {
        temp = 0.0;
        for (j=0; j<n; ++j) temp += a[in+j]*b[j];
        c[i]=temp;
    } // end for i

    return;
}


/**
 * \fn void ddense_mul (const REAL *a, const REAL *b, REAL *c, const INT n)
 *
 * \brief Compute the matrix product of two small full matrices a and b, stored in c
 *
 * \param a   Pointer to the REAL array which stands a n*n matrix
 * \param b   Pointer to the REAL array which stands a n*n matrix
 * \param c   Pointer to the REAL array which stands a n*n matrix
 * \param n   Dimension of the matrix
 *
 * \author Xiaozhe Hu
 * \date   04/21/2010
 */
void ddense_mul (const REAL  *a,
                 const REAL  *b,
                 REAL        *c,
                 const INT    n)
{
    const INT n2 = n*n;
    INT i,j,k;
    REAL temp;

    for (i=0; i<n2; i+=n) {
        for (j=0; j<n; ++j) {
            temp = 0.0;
            for (k=0; k<n; ++k) temp += a[i+k]*b[k*n+j];
            c[i+j] = temp;
        } // end for j
    } // end for i

    return;
}

/**
 * \fn void ddense_ypAx (const REAL *A, const REAL *x, REAL *y, const INT n)
 *
 * \brief Compute y := y + Ax, where 'A' is a n*n dense matrix
 *
 * \param A   Pointer to the n*n dense matrix
 * \param x   Pointer to the REAL array with length n
 * \param y   Pointer to the REAL array with length n
 * \param n   Dimension of the dense matrix
 *
 */
void ddense_ypAx (const REAL  *A,
                  const REAL  *x,
                  REAL        *y,
                  const INT    n)
{
    INT i,j,k;

    for ( k = i = 0; i < n; i++, k+=n ) {
        for ( j = 0; j < n; j++ ) {
            y[i] += A[k+j]*x[j];
        }
    }

    return;
}

/**
 * \fn void ddense_ymAx (const REAL *A, const REAL *x, REAL *y, const INT n)
 *
 * \brief Compute y := y - Ax, where 'A' is a n*n dense matrix
 *
 * \param A   Pointer to the n*n dense matrix
 * \param x   Pointer to the REAL array with length n
 * \param y   Pointer to the REAL array with length n
 * \param  n   the dimension of the dense matrix
 *
 * \author Xiaozhe Hu
 * \date   2010/10/25
 *
 */
void ddense_ymAx (const REAL  *A,
                  const REAL  *x,
                  REAL        *y,
                  const INT    n)
{
    INT i,j,k;

    for ( k = i = 0; i < n; i++, k+=n ) {
        for ( j = 0; j < n; j++ ) {
            y[i] -= A[k+j]*x[j];
        }
    }

    return;
}

/**
 * \fn void ddense_aAxpby (const REAL alpha, const REAL *A, const REAL *x,
 *                                 const REAL beta, REAL *y, const INT n)
 *
 * \brief Compute y:=alpha*A*x + beta*y
 *
 * \param alpha   REAL factor alpha
 * \param A       Pointer to the REAL array which stands for a n*n full matrix
 * \param x       Pointer to the REAL array with length n
 * \param beta    REAL factor beta
 * \param y       Pointer to the REAL array with length n
 * \param n       Length of array x and y
 *
 */
void ddense_aAxpby (const REAL   alpha,
                    const REAL  *A,
                    const REAL  *x,
                    const REAL   beta,
                    REAL        *y,
                    const INT    n)
{
    INT     i,j,k;
    REAL    tmp = 0.0;

    if ( alpha == 0 ) {
        for (i = 0; i < n; i ++) y[i] *= beta;
        return;
    }

    // y := (beta/alpha)y
    tmp = beta / alpha;
    if ( tmp != 1.0 ) {
        for (i = 0; i < n; i ++) y[i] *= tmp;
    }

    // y := y + A*x
    for ( k = i = 0; i < n; i++, k+=n ) {
        for (j = 0; j < n; j ++) {
            y[i] += A[k+j]*x[j];
        }
    }

    // y := alpha*y
    if ( alpha != 1.0 ) {
        for ( i = 0; i < n; i ++ ) y[i] *= alpha;
    }
}

/**
 * \fn void ddense_inv_inplace (REAL *a, const INT n)
 *
 * \brief Compute the inverse of a matrix using Gauss Elimination
 *
 * \param a   Pointer to the REAL array which stands a n*n matrix
 * \param n   Dimension of the matrix
 *
 * \author Xiaozhe Hu
 * \date   05/01/2010
 */
void ddense_inv_inplace (REAL      *a,
                         const INT  n)
{
    INT i,j,k,l,u,kn,in;
    REAL alinv;

    for (k=0; k<n; ++k) {

        kn = k*n;
        l  = kn+k;

        if (ABS(a[l]) < SMALLREAL) {
            printf("### HAZMATH ERROR: Diagonal entry is close to zero! ");
            printf("diag_%d = %.2e! [%s]\n", k, a[l], __FUNCTION__);
            exit(ERROR_SOLVER_EXIT);
        }
        alinv = 1.0/a[l];
        a[l] = alinv;

        for (j=0; j<k; ++j) {
            u = kn+j; a[u] *= alinv;
        }

        for (j=k+1; j<n; ++j) {
            u = kn+j; a[u] *= alinv;
        }

        for (i=0; i<k; ++i) {
            in = i*n;
            for (j=0; j<n; ++j)
                if (j!=k) {
                    u = in+j; a[u] -= a[in+k]*a[kn+j];
                } // end if (j!=k)
        }

        for (i=k+1; i<n; ++i) {
            in = i*n;
            for (j=0; j<n; ++j)
                if (j!=k) {
                    u = in+j; a[u] -= a[in+k]*a[kn+j];
                } // end if (j!=k)
        }

        for (i=0; i<k; ++i) {
            u=i*n+k; a[u] *= -alinv;
        }

        for (i=k+1; i<n; ++i) {
            u=i*n+k; a[u] *= -alinv;
        }

    } // end for (k=0; k<n; ++k)
}

/**
 * \fn REAL ddense_Linf (const REAL *A, const INT n )
 *
 * \brief Compute the L infinity norm of A
 *
 * \param A   Pointer to the n*n dense matrix
 * \param n   the dimension of the dense matrix
 *
 * \author Xiaozhe Hu
 * \date   05/26/2014
 */
REAL ddense_Linf (const REAL  *A,
                     const INT    n)
{

    REAL norm = 0.0, value;

    INT i,j;

    for ( i = 0; i < n; i++ ) {
        for ( value = 0.0, j = 0; j < n; j++ ) {
            value = value + ABS(A[i*n+j]);
        }
        norm = MAX(norm, value);
    }

    return norm;
}

/**
 * \fn void ddense_identity (REAL *a, INT n, INT n2)
 *
 * \brief Set a n*n full matrix to be a identity
 *
 * \param a      Pointer to the REAL vector which stands for a n*n full matrix
 * \param n      Size of full matrix
 * \param n2     Length of the REAL vector which stores the n*n full matrix
 *
 * \author Xiaozhe Hu
 * \date   2010/12/25
 */
void ddense_identity (REAL      *a,
                      const INT  n,
                      const INT  n2)
{
    memset(a, 0X0, n2*sizeof(REAL));

    INT l;
    for (l = 0; l < n; l ++) a[l*n+l] = 1.0;

}

/************************************** END ***************************************************/
/*! USED TO BE: \file src/utilities/solve_full.c
 *
 */
