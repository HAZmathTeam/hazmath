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
void dsygv_( INT *itype, char *jobz, char *uplo, INT *n,	\
 	     REAL *a, INT *lda, REAL *b, INT *ldb, REAL *w,	\
	     REAL *work, INT *lwork, INT *info );
// this should be replaced by lapack.h if any exists in a standard install of lapack. 
