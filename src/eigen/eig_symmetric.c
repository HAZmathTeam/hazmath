
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
  //   evalues must be have enough space for n reals; if evectors is a
  //   NULL: compute only eigenvalues; otherwise evectors must be n*n
  //   reals upon entering here.
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
