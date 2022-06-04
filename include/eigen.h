//
//
//  eigen.h
//
//
//  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20200529
//
//
#ifndef _eigen_h
#define _eigen_h
#endif

// LAPACK Routines

// Symmetric eigenvalue computation
void dsygv_( INT *itype, char *jobz, char *uplo, INT *n,	\
 	     REAL *a, INT *lda, REAL *b, INT *ldb, REAL *w,	\
	     REAL *work, INT *lwork, INT *info );

// Symmetric tridiagonal eigenvalur computation
void dsterf_(INT *N, REAL *D, REAL *E, INT *INFO);

// svd
void dgesvd_( char *joba, char *jobu, INT *m,	INT *n, \
       REAL *a, INT *lda, REAL *s, REAL *u, INT *ldu, REAL *vt,	\
       INT *ldvt, REAL *work, INT *lwork, INT *info );

// qr factorization
void dgeqrf_(INT *m, INT *n, REAL *a, INT *lda, REAL *tau, REAL* work, \
      INT* lwork, INT* info);

// extract Q from QR factorization
void dorgqr_(INT *m, INT *n1, INT *n2, REAL *a, INT *lda, REAL *tau, REAL* work, \
      INT* lwork, INT* info);
//
// nonsymmetric eigenvalue computations and SVD
//
void  dgeevx_( char *balanc, char *jobvl, char *jobvr, char *sense,	\
	       INT *n, REAL *a, INT *lda, REAL *wr, REAL *wi,		\
	       REAL *vl, INT *ldvl, REAL *vr, INT *ldvr,		\
	       INT *ilo, INT *ihi, REAL *scale, REAL *abnrm,		\
	       REAL *rconde, REAL *rcondv, REAL *work, INT *lwork,	\
	       INT *iwork, INT *info );

void  dgesvd_(char *jobu, char *jobvt, INT *m, INT *n,	\
	      REAL *a, INT *lda, REAL *s,		\
	      REAL *u,INT *ldu, REAL *vt, INT *ldvt,	\
	      REAL *work,INT *lwork, INT *info);
// this should be replaced by lapack.h if any exists in a standard install of lapack.
