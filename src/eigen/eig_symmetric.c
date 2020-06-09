
/*! \file src/eigen/eigen_symmetric.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 8/19/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 * \brief Eigenvalue computations using LAPACK
 *
 */

#include "hazmath.h"
INT eigsymm(dCSRmat *A,dCSRmat *B,REAL *evalues, REAL *evectors)
{
  INT info=-22;
#if WITH_LAPACK
  /*    
	evalues must have enough space for n reals; if evectors is a 
	NULL: compute only eigenvalues; otherwise evectors must be n*n 
	reals upon entering here. 
  */
  char uplo='U',jobz='N';
  INT n=A->row;
  INT itype=1,lwork=3*n+6; // itype=1 is to solve A*x =lambda*B*x
  REAL *af=NULL, *bf=NULL, *work=NULL;
  REAL *w=evalues;
  if(evectors){
    jobz='V';
    af=evectors;
  } else {
    af=calloc(n*n,sizeof(REAL));
  }
  bf=calloc(n*n+lwork,sizeof(REAL));
  dcsr2full(A,af);
  /* 
     if the matrix is not symmetric we need to do here c2r before we
     call the fortran to make the matrix stored "fortran friendly". We
     use c2r here just to be consistent
  */
  c2r(n,n,sizeof(REAL),(void *)af);
  /*B must be always symmetric as it should define an inner product  so no need of c2r; */
  dcsr2full(B,bf);
  work = bf + n*n;
  dsygv_( &itype, &jobz, &uplo, &n, af, &n, bf, &n, w, work,
  	  &lwork, &info );
  free(bf);
  if(!evectors)
    free(af);
  else{
    // fortran orders memory by columns, so we need to order by rows.
    c2r(n,n,sizeof(REAL),(void *)evectors);
  }
  if(info){
    fprintf(stderr,"\nXXX: lapack error info during eigenvalue computations=%d (n=%d)\n",info,n);fflush(stderr);
    exit(16);
  }
  return info;
#else
  error_extlib(252, __FUNCTION__, "LAPACK");
  return info;
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
