/*! \file src/eigen/eigen_drivers.c
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
     if the matrix is not symmetric we need to do here col_2_row before we
     call the fortran to make the matrix stored "fortran friendly". We
     use col_2_row here just to be consistent
  */
  // TODO: Shouldn't this be row_2_col?  But it's symmetric anyway??
  col_2_row(n,n,sizeof(REAL),(void *)af);
  /*B must be always symmetric as it should define an inner product  so no need of col_2_row; */
  dcsr2full(B,bf);
  work = bf + n*n;
  dsygv_( &itype, &jobz, &uplo, &n, af, &n, bf, &n, w, work,
  	  &lwork, &info );
  free(bf);
  if(!evectors)
    free(af);
  else{
    // fortran orders memory by columns, so we need to order by rows.
    col_2_row(n,n,sizeof(REAL),(void *)evectors);
  }
  if(info){
    fprintf(stderr,"\nXXX: lapack error info during eigenvalue computations=%lld (n=%lld)\n",(long long )info,(long long )n);fflush(stderr);
    exit(16);
  }
  return info;
#else
  error_extlib(252, __FUNCTION__, "LAPACK");
  return info;
#endif
}
/***********************************************************************/
INT eiggeneral(INT n,REAL *a,REAL *wr, REAL *wi)
{
  INT info=-22;
#if WITH_LAPACK
  INT i,j;
  /*
    evalues must have enough space for 2*n reals;
    compute only eigenvalues; 
    IMPORTANT: enter here with a transposed instead of a. 
    a is destroyed on output. 
  */
  /* balanc can be 'N', 'P','S' or 'B'; use 'N' or 'B'*/
  char balanc='B',jobvl='N',jobvr='N',sense='N';
  INT lwork=n*(n+6)+16; //16 more....
  REAL *at=calloc(n*n+lwork + 3*n+2,sizeof(REAL));
  INT *iwork=calloc(2*n,sizeof(INT));
  REAL *work=at+n*n;
  REAL *scale = work + lwork;
  REAL *rconde=scale+n;
  REAL *rcondv=rconde+n;  
  REAL *vl=rcondv+n;
  REAL *vr=vl+1; // end is +1;
  INT lda=n,ldvl=1,ldvr=1,ilo,ihi;
  REAL abnrm;
  for(j=0;j<n;j++){
    for(i=0;i<n;i++){
      at[i*n+j]=a[j*n+i];
    }
  }
  /*  
      subroutine dgeevx( balanc, jobvl, jobvr, sense, n, a, lda, wr, wi,
      $                   vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm,
      $                   rconde, rcondv, work, lwork, iwork, info )
  */
  dgeevx_( &balanc, &jobvl, &jobvr, &sense,			\
	   &n, at, &lda, wr, wi,					\
	   vl, &ldvl, vr, &ldvr, &ilo, &ihi, scale, &abnrm,	\
	   rconde, rcondv, work, &lwork, iwork, &info);
  //  fprintf(stdout,"\n%%%s: INFO=%i;\n",__FUNCTION__,info);
  free(at);
  free(iwork);
  return info;
#else
  error_extlib(252, __FUNCTION__, "LAPACK");
  return info;
#endif
}
/**********************************************************************/
INT svdgeneral(INT nrow, INT ncol, REAL *a,REAL *smin, REAL *w)
{
  /* subroutine dgesvd	(	character 	jobu,
     character 	jobvt,
     integer 	m,
     integer 	n,
     double precision, dimension( lda, * ) 	a,
     integer 	lda,
     double precision, dimension( * ) 	s,
     double precision, dimension( ldu, * ) 	u,
     integer 	ldu,
     double precision, dimension( ldvt, * ) 	vt,
     integer 	ldvt,
     double precision, dimension( * ) 	work,
     integer 	lwork,
     integer 	info)		
  */
  INT info=-22;
#if WITH_LAPACK  
  char 	jobu='N',jobvt='N';
  INT m=nrow,n=ncol;
  INT i,j,memu,memvt,mnmin,mnmax;
  INT lda=m,ldu=m,ldvt=n,lwork=-10;
  REAL *allwork=NULL, *at=NULL, *s=NULL,*vt = NULL,*work = NULL, *u=NULL;
  //  LWORK=MAX(1,3*MIN(M,N) + MAX(M,N),5*MIN(M,N))
  /******************************************************************/
  // OUTPUT: w[] On output, w[] contains the RIGHT singular vector
  // corresponding to the minimal singular value, i.e. the last row in
  // vt, where a=u*s*vt) SIZE of w must be n at least.
  //
  // For different output you need to lookup the documentation of lapack. 
  //
  if(m<n){
    mnmin=m;
    mnmax=n;
  } else {
    mnmin=n;
    mnmax=m;
  }
  jobvt='S';
  ldvt=mnmin;
  memvt=ldvt*n; // either n*n if n is smaller or m*n if n is larger
  memu =1;
  ldu = 1;
  jobu='N'; //just in case
  lwork = 8*mnmin + mnmax;
  allwork=calloc(m*n + mnmin + memu + memvt + lwork,sizeof(REAL)); 
  at=allwork;
  s=at + m*n;  
  vt    = s  +  mnmin;
  u   = vt  +  memvt;
  work = u +  memu;
  REAL *memend=work+lwork;
  if((memend-at)>(m*n + mnmin + memu + memvt + lwork)) {
    fprintf(stdout,"\n%% HAZMATH WARNING ***  IN %s: REALLOC bc of insufficient memory (%lld REALs);\n",__FUNCTION__,(long long int)((memend-at)-(m*n + mnmin + memu + memvt + lwork)));
    allwork=realloc(allwork,(memend-at)*sizeof(REAL));
  }
  for(j=0;j<n;j++){
    for(i=0;i<m;i++){
      at[j*m+i]=a[i*n+j];// transpose to go to fortran; anyway a is
			 // destroyed after entering dgesvd so it is
			 // better to have a copy.
    }
  }
  dgesvd_(&jobu,&jobvt,&m,&n,			\
	  at,&lda,s,				\
	  u,&ldu,vt,&ldvt,			\
	  work,&lwork,&info);		
  *smin=s[mnmin-1];
  if(jobvt=='S'){
    // we need to get the last row of vt:
    i=mnmin-1;// last column of vt, as fortran stores this way
    for(j=0;j<n;j++) {
      w[j]=vt[i+j*mnmin];//
      //      fprintf(stdout,"\nvt000(%d)=%.18e;",j+1,vt[j*ldvt+ldvt-1]);// reading by rows
    }
  } else {
    fprintf(stderr,"\n\n*** ERROR: we should have either jobu=\'N\' or jobvt=\'S\'. Exiting.\n\n");
    exit(32);
  }
  free(allwork);
  if(info!=0){
    fprintf(stdout,"\n%% WARNING: IN %s: INFO IS NOT ZERO; INFO=%lld;\n",__FUNCTION__,(long long )info);
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
