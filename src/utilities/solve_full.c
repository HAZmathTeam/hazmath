/*! \file src/utilities/solve_full.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 5/13/15.
 *  Copyright 2016__HAZMATH__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/25/2016
 *  \note: modified on 10/12/2017 (ltz1)
 *
 *  \note routines to solve Ax=b, calculate inv(A), transpose, axpy,
 *        and qr.
 *
 */
#include "hazmath.h"
/**********************************************************************/
INT solve_pivot(INT dopivot, INT n, REAL *A, REAL *b, INT *p,REAL *piv)
{

/* \brief Solution of a linear system A*x = b with scaled partial
   pivoting. The right hand side is overwritten on output with the
   solution, that is x[] is the same as b[] on return.
 *
 * \param dopivot (switch [0/1]; if idopivot=1 then do the LU
 * decomposition of A. The L and U factors are stored in the lower and
 * upper triangle of A respectively.  if dopivot=0, then A is assumed
 * to be already LU decomposed and this only performs the forward and
 * backward substitutions.
 * \param n    number of rows (and columns) of A.
 * \param A    the matrix as a one dimensional array by rows
 *
 */
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
 * \fn SHORT lufull(INT n, REAL *deta, REAL *A,INT *p,REAL *piv)
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
SHORT lufull(INT dopivot, INT n, REAL *deta, REAL *A,INT *p,REAL *piv)
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
     //print_full_mat(n,n,A,"A");
     //fprintf(stderr,"\n*** ERR: zero diagonal in lufull; pivot=%e\n",A[p[k]*n+k]);
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
   //   print_full_mat(n,n,A,"A");
   //   fprintf(stderr,"\n*** det=%e\n",*deta);
   return (SHORT )0;
 }
/**************************************************************************/
/*
 * \fn void invfull(INT n, REAL *Ainv, REAL *A, void *wrk)
 *
 * \brief Inverts a general (nxn) matrix A
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
void invfull(REAL *Ainv,INT n, REAL *A, void *wrk)
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
  //  print_full_mat(n,n,Awrk,"At");
  for(i=0;i<n;i++){
    ni=n*i;
    for(j=0;j<n;j++)
      Ainv[ni+j]=0e0;
    Ainv[ni+i]=1e0;
    // solve with rhs a row of the identity;
    solve_pivot((!i),n, Awrk, (Ainv+ni), p, piv);
  }
  return;
}
/**************************************************************************/
/*
 * \fn void abybfull(const INT m, const INT p, REAL *c, REAL *a,
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
 * \fn void abyvfull(const INT m, REAL *y,REAL *a, REAL *x, const INT n)
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
void abyvfull(const INT m, REAL *y,REAL *a, REAL *x, const INT n)
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
 * \fn atbyvfull(const INT m, REAL *y,REAL *a, REAL *x, const INT n)
 *
 * \brief Computes y = transpose(a)*x+y;
 *
 * \param m is the number of rows in a (columns in transpose(a)).
 * \param a is an matrix m by n,
 * \param x is an m by 1 vector,
 * \param y is an n by 1 vector.
 *
 */
void atbyvfull(const INT m, REAL *y,REAL *a, REAL *x, const INT n)
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
 * \fn qr_full(const INT m, const INT n, REAL *A, REAL *Q, REAL *R)
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
void qr_full(const INT m, const INT n, REAL *A, REAL *Q, REAL *R)
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
 * \fn INT svd_full(INT m,INT n, REAL *A, REAL *U, REAL *VT, REAL* S,INT computeUV)
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
INT svd_full(INT m, INT n, REAL *A, REAL *U, REAL *VT, REAL* S,INT computeUV)
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
 * \fn INT qr_full_lapack(INT m,INT n, REAL *A, REAL *Q, REAL *R,INT computeR)
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
INT qr_full_lapack(INT m,INT n, REAL *A,REAL * Q,REAL *R,INT computeR)
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
/**********************************************************************/
INT solve_pivot_l(INT dopivot,			\
		  INT n,			\
		  REAL16 *A,			\
		  REAL16 *b,			\
		  INT *p,			\
		  REAL16 *piv)
{
/* \brief Solution of a linear system A*x = b with scaled partial
   pivoting in LONG double. The right hand side is overwritten on
   output with the solution, that is x[] is the same as b[] on return.
 *
 * \param dopivot (switch [0/1]; if idopivot=1 then do the LU
 * decomposition of A. The L and U factors are stored in the lower and
 * upper triangle of A respectively.  if dopivot=0, then A is assumed
 * to be already LU decomposed and this only performs the forward and
 * backward substitutions.
 * \param n    number of rows (and columns) of A.
 * \param A    the matrix as a one dimensional array by rows (long double)
 *
 */
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
