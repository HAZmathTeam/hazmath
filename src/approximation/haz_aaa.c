/*! \file src/approximation/haz_aaa.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note containing programs for the AAA algorithm for rational approximation
 * 
 *  \note See: Yuji Nakatsukasa, Olivier S`ete, and Lloyd N.
 *             Trefethen.  The AAA algorithm for rational
 *             approximation.SIAM J. Sci.Comput., 40(3):A1494–A1522,
 *             2018. https://doi.org/10.1137/16M1106122
 *
 *  \note 20201025 (ltz)
 */
#include "hazmath.h"
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/**********************************************************************/
/*!
 * \fn static void get_res(const INT m,REAL16 *res,  REAL16 *z,	REAL16 *p,REAL16 *f,void *wrk)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
static void get_res(const INT m,		\
		    REAL16 *resr,		\
		    REAL16 *resi,		\
		    REAL16 *pr,			\
		    REAL16 *pi,			\
		    REAL16 *z,			\
		    REAL16 *f,			\
		    void *wrk)
{
  // wrk should have space m real16 and m integer. 
  //  given poles and right hand side, this finds the residues
  // Solves 
  INT j,k,m01=m-1;
  REAL16 zj,zpmod;
  REAL16 *apolr=(REAL16 *)calloc(m*m,sizeof(REAL16));
  REAL16 *apoli=(REAL16 *)calloc(m*m,sizeof(REAL16));
  //  REAL16 *piv = (REAL16 *)wrk;  
  //  INT *perm=(INT *)(wrk+2*m*sizeof(REAL16));
  for(j=0;j<m;j++){
    zj=z[j];
    for(k=0;k<m01;k++){
      zpmod=(zj-pr[k])*(zj-pr[k])+pi[k]*pi[k];
      apolr[j*m+k]=(zj-pr[k])/zpmod;
      apoli[j*m+k]=pi[k]/zpmod;
      /* fprintf(stdout,"\nzCpr(%d,%d)=%.18Le;",j+1,k+1,apolr[j*m+k]); */
      /* fprintf(stdout,"\nzCpi(%d,%d)=%.18Le;",j+1,k+1,apoli[j*m+k]); */
    }
    apolr[j*m+m01]=1e00;
    apoli[j*m+m01]=0e00;
    /* fprintf(stdout,"\nzCpr(%d,%d)=%.18Le;",j+1,m01+1,apolr[j*m+m01]); */
    /* fprintf(stdout,"\nzCpi(%d,%d)=%.18Le;",j+1,m01+1,apoli[j*m+m01]); */
  }
  //rhs
  for(j=0;j<m;j++) resr[j]=f[j];
  for(j=0;j<m;j++) resi[j]=0e0;
  zdense_solve_pivot_l(m,apolr,apoli,resr,resi);
  free(apolr);
  free(apoli);
  return;
}
/**********************************************************************/
/*!
 * \fn INT residues_poles(INT m,REAL *zd,REAL *wd,REAL *fd, 
 *                       REAL *resdr,REAL *resdi,REAL *poldr,REAL *poldi)
 *
 * \brief given nodes and weights computes residues and poles. 
 *
 * \param
 *
 * \return
 *
 * \note the res[m-1] is the coefficient corresponding to the constant
 *       shift.
 *
 */
INT residues_poles(INT m,				\
		   REAL *zd,				\
		   REAL *wd,				\
		   REAL *fd,				\
		   REAL *resdr,				\
		   REAL *resdi,				\
		   REAL *poldr,				\
		   REAL *poldi)
{
  INT m1=m+1,m01=m-1,mm1=m1*m1,i,j;
  //  REAL16 swp;
  // MEMORY
  size_t memall = 9*m1*sizeof(REAL16)+2*m1*sizeof(INT);// 2*m1 because we may have imaginary part;
  void *wrk=(void *)calloc(memall,sizeof(char));
  void *wrk0=wrk;
  REAL16 *z=(REAL16 *)wrk;
  wrk+=m1*sizeof(REAL16);
  REAL16 *w=(REAL16 *)wrk;
  wrk+=m1*sizeof(REAL16);
  REAL16 *f=(REAL16 *)wrk;
  wrk+=m1*sizeof(REAL16);
  REAL16 *polr=(REAL16 *)wrk;
  wrk+=m1*sizeof(REAL16);
  REAL16 *poli=(REAL16 *)wrk;
  wrk+=m1*sizeof(REAL16);
  REAL16 *resr=(REAL16 *)wrk;
  wrk+=m1*sizeof(REAL16);
  REAL16 *resi=(REAL16 *)wrk;
  wrk+=m1*sizeof(REAL16);
  REAL16 *pivot=(REAL16 *)wrk;
  wrk+=2*m1*sizeof(REAL16);
  //  INT *perm=(INT *)wrk;
  wrk+=2*m1*sizeof(INT);
  if((wrk-wrk0)>memall) {
    fprintf(stdout,"\nMemory error in %s: allocated=%ld .lt. needed=%ld\n",__FUNCTION__,memall,wrk-wrk0);
    exit(16);    
  }
  // END memory
  // grab the input params and convert them to higher precision
  for(j=0;j<m;j++){
    w[j]=(REAL16 )wd[j];
    z[j]=(REAL16 )zd[j];
    f[j]=(REAL16 )fd[j];
  }
  REAL16 *a_mat=(REAL16 *)calloc(m1*m1,sizeof(REAL16));
  REAL16 *b_mat=(REAL16 *)calloc(m1*m1,sizeof(REAL16));
  REAL16 *x000=(REAL16 *)calloc(m1*m1,sizeof(REAL16));
  REAL16 *alphar=(REAL16 *)calloc(m1,sizeof(REAL16));
  REAL16 *alphai=(REAL16 *)calloc(m1,sizeof(REAL16));
  REAL16 *beta=(REAL16 *)calloc(m1,sizeof(REAL16));
  INT *iter000=(INT *)calloc(m1,sizeof(INT));
  memset(a_mat,0,mm1*sizeof(REAL16));
  memset(b_mat,0,mm1*sizeof(REAL16));
  for(i=1;i<m1;i++){
    a_mat[i*m1+i]=z[i-1];
    a_mat[i*m1]=1e00;
    a_mat[i]=w[i-1];
    b_mat[i*m1+i]=1e0;
  }  
  /* for(i=0;i<m1;i++){ */
  /*   for(j=0;j<m1;j++){ */
  /*     fprintf(stdout,"\na_mat(%d,%d)=%.18Le;",i+1,j+1,a_mat[i*m1+j]);fflush(stdout); */
  /*   } */
  /* } */
  /* for(i=0;i<m1;i++){ */
  /*   for(j=0;j<m1;j++){ */
  /*     fprintf(stdout,"\nb_mat(%d,%d)=%.18Le;",i+1,j+1,b_mat[i*m1+j]);fflush(stdout); */
  /*   } */
  /* } */
  /* fflush(stdout); */
  INT wantx=0;
  REAL16 tol000,epsa,epsb;
  tol000=powl(2e0,-56e0);
  qzhes(m1,m1,a_mat,b_mat,wantx,x000);
  qzit(m1,m1,&a_mat,&b_mat,tol000,&epsa,&epsb,&iter000,wantx,&x000);
  qzval(m1,m1,
	&a_mat,&b_mat,
	epsa,epsb,
	&alphar,&alphai,&beta,
	wantx,&x000);
  /* for(j=0;j<m1;j++){ */
  /*   fprintf(stdout,"\nalphar(%d)=%.18Le;",j+1,alphar[j]); */
  /*   fprintf(stdout,"alphai(%d)=%.18Le;",j+1,alphai[j]); */
  /*   fprintf(stdout,"beta000(%d)=%.18Le;",j+1,beta[j]); */
  /* } */
  /* fprintf(stdout,"\ntol000=%Le,epsa=%Le,epsb=%Le\n",tol000,epsa,epsb);fflush(stdout); */
  // do two passes to remove the two eigenvalues with minimal beta:
  INT j1=0,j2;
  REAL16 betamin1=fabsl(beta[0]),betamin2;
  for(i=1;i<m1;i++){
    if(betamin1 > fabsl(beta[i])){
      betamin1=fabsl(beta[i]);
      j1=i;
    }
  }
  j2=-1;
  betamin2=1e20;
  for(i=0;i<m1;i++){
    if(i==j1) continue;
    if(betamin2 > fabsl(beta[i])){
      betamin2=fabsl(beta[i]);
      j2=i;
    }
  }
  // once j1 and j2 are known, we find the poles:
  //  INT imag_flag=0;
  j=0;
  for(i=1;i<m1;i++){
    if(i == j1 || i == j2) continue;
    // nonzero imaginary part:
    //    if(fabsl(alphai[i])>1e-16) {
    //      fprintf(stdout,"\n WARNING: pole[%d] haz nonzero imaginary part. Skipping it...\n");
    //    }
    polr[j]=alphar[i]/beta[i];
    if(alphai[i]!=0e0) {
      //      imag_flag=1;
      poli[j]=alphai[i]/beta[i];
    } else {
      poli[j]=0e0;
    }
    ++j;
  }
  //
  /* for(j=0;j<m01;j++) fprintf(stdout,"\npol0(%d)=%.16Le+sqrt(-1)*(%.16Le);",j+1,polr[j],poli[j]); */
  /* for(j=0;j<m;j++) fprintf(stdout,"\nz0(%d)=%.16Le;",j+1,z[j]); */
  /* for(j=0;j<m;j++) fprintf(stdout,"\nf0(%d)=%.16Le;",j+1,f[j]); */
  /* fprintf(stdout,"\n");fflush(stdout); */
  //
  free(a_mat);
  free(b_mat);
  free(x000);
  free(alphar);
  free(alphai);
  free(beta);
  free(iter000);
  ///////////////////////////////////////////////////////////////////////////////
  // get the residual  
  get_res(m,resr,resi,polr,poli,z,f,(void *)pivot);
  ///////////////////////////////////////////////////////////////////////////////
  /* for(j=0;j<m01;j++){ */
  /*   fprintf(stdout,"\npolr(%d)=%.18Le;poli(%d)=%.18e;",j+1,polr[j],j+1,poli[j]); */
  /* } */
  /* fprintf(stdout,"\n"); */
  /* for(j=0;j<m;j++){ */
  /*   fprintf(stdout,"\nresr(%d)=%.18Le;resi(%d)=%.18Le;",j+1,resr[j],j+1,resi[j]); */
  /* } */
  /* fprintf(stdout,"\n\n"); */
  ///////////////////double on output.
  //
  for(j=0;j<m;j++) resdr[j]=(REAL )resr[j];
  for(j=0;j<m01;j++) poldr[j]=(REAL )polr[j];
  for(j=0;j<m;j++) resdi[j]=(REAL )resi[j];
  for(j=0;j<m01;j++) poldi[j]=(REAL )poli[j];
  //
  free(wrk0);
  return 0;
}

/**********************************************************************/
/*!
 * \fn REAL get_rpzwf(INT numval,REAL16 *z, REAL16 *f, REAL **rpzwf,
 *                      INT *mmax_in,*INT *m_out,*REAL tolaaa,*INT
 *                      print_level)
 *
 *
 * \brief Uses AAA to compute residues, poles, nodes and weights using
 *        a rational fuunction approximation to "func()". 
 *
 * \param
 *
 * \return
 *
 * \note Uses long doubles as well as doubles. 
 *
 */
REAL get_rpzwf(INT numval,REAL16 *z, REAL16 *f,		\
	       REAL **rpzwf,				\
	       INT *mmax_in,				\
	       INT *m_out,				\
	       REAL tolaaa,
	       INT print_level)
{
  /*
   *  \param mmax_in is the max number of
   *  nodes taking part in the interpolation, mmax_in should be much
   *  smaller than z_in.row;
   *
   * \param m is the number of nodes in the final interpolation after
   *    tolerance tolaaa is achieved or mmax is reached.
   *   parameters for the function we are approximating.
   *   For example: s[2]*(x^s[0])+s[3]*(x**s[1])
   *   after a call to get_rpzwf, rpzwf[k] are pointers to the following
   *   arrays (each of m+1 elements but at most m are used:
   *
   *   rpzwf[0]->  real part of residues(resr[]);
   *   rpzwf[1]->  imaginary part of residues(resi[]);
   *   rpzwf[2]-> real part of poles (polr[]);
   *   rpzwf[3]-> imaginary part of poles (pol[]);
   *   rpzwf[4]-> nodes (z[]);
   *   rpzwf[5]-> weights (w[]);
   *   rpzwf[6]-> function values (f[]) 
   * (also as last entry rpzwf[0] contains the free term c[m-1]; 
   *
   *   the rational approximation is:
   ****   r(z)=res[m-1] + \sum_{i=0}^{m-2} res[i]/(z-pol[i]); 
   *
  */
  ////////////
  // rpzwf must contain 5 pointers: rpzwf[0:1] are the residues and
  // rpzwf[2:3] are the poles. **rpzwf must be allocated earlier, but
  // rpzwf[k], k=0:4 are allocated here.  func() is the function we want to
  // approximate and param are the parameters it may depend on
  INT m,i,j,k,krow,kmax,mmax;
  REAL16 r,wzkj,fnum,fden,rm,rmax,swp;
  REAL smin=-1e20;
  if(mmax_in[0]<2) mmax_in[0]=2;
  mmax=(INT )(numval/2);
  if(mmax_in[0] < mmax)
    mmax=mmax_in[0];
  else 
    mmax_in[0]=mmax;
  REAL16 tol=powl(2e0,-52e0);
  if(((REAL16 )tolaaa) > tol)
    tol=(REAL16 )tolaaa;
  INT mem=numval*mmax*sizeof(REAL) + mmax*sizeof(REAL);
  void *zfw_long=calloc(mem,sizeof(char));
  REAL *wd=(REAL *)zfw_long;
  zfw_long+=(mmax)*sizeof(REAL);
  REAL *x21d=(REAL *)zfw_long;
  zfw_long+=numval*mmax*sizeof(REAL);// end of it
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  // form f; take the mean of f; form the matrix;
  // set up z;  
  m=0;
  r=0e0;
  for(j=0;j<numval;++j){
    r+=f[j];
  }
  r/=(REAL16 )numval;
  REAL16 rmax_min=1e0/tol;
  //  INT m_minerr=mmax;
  while(1){
    rmax=-1.e00;
    kmax=-1;
    if(m==0){
      for(k=0;k<numval;k++){
	rm=fabsl(r-f[j]);
	if(rm>rmax){
	  rmax=rm;
	  kmax=k;
	}
      }
    } else {
      for(k=m;k<numval;k++){
	fnum=0e0;
	fden=0e0;
	for(j=0;j<m;j++){
	  wzkj=((REAL16 )wd[j])/(z[k]-z[j]);
	  fden+=wzkj;
	  fnum+=f[j]*wzkj;
	}
	rm=fabsl(fnum/fden-f[k]);
	if(rm>rmax){
	  rmax=rm;
	  kmax=k;
	}
      }
    }
    if(rmax_min>rmax){
      //      m_minerr=m;
      rmax_min=rmax;
    }
    if(print_level>2){
      fprintf(stdout,"\n%%%%iter=%d | rmax=%.20Le at %7d;",m,rmax,kmax); fflush(stdout);
    }
    if(rmax<tol || m>=(mmax-1)) break;
    swp=z[kmax];  z[kmax]=z[m];   z[m]=swp;
    swp=f[kmax];  f[kmax]=f[m];   f[m]=swp;
    m++;
    // form x21d (this needs improvement!:
    for(k=m;k<numval;k++){
      krow=k-m;
      for(j=0;j<m;j++){
	x21d[krow*m+j]=(REAL )((f[k]-f[j])/(z[k]-z[j]));
      }
    }
    INT info=svdgeneral(numval,m,x21d,&smin,wd);
    /////////    if(m>1) print_full_mat(m,1,wd,"wd");    
    if(info!=0){
      fprintf(stdout,"\n%% *** HAZMATH WARNING*** IN %s: SVD-INFO IS NOT ZERO; INFO=%d;\n",__FUNCTION__,info);
    }
  }
  //  fprintf(stdout,"\n\n%%%% MINerr=%.18Le; at iter=%d\n",rmax_min,m_minerr); fflush(stdout);
  m_out[0]=m;
  //
  // allocate space for results
  rpzwf[0]=calloc(7*(m+1), sizeof(REAL));
  rpzwf[1]=rpzwf[0] + m + 1;
  rpzwf[2]=rpzwf[1] + m + 1;
  rpzwf[3]=rpzwf[2] + m + 1;
  rpzwf[4]=rpzwf[3] + m + 1;
  rpzwf[5]=rpzwf[4] + m + 1;
  rpzwf[6]=rpzwf[5] + m + 1;
  for(i=0;i<m;i++) rpzwf[4][i]=z[i];
  memcpy(rpzwf[5],wd,m*sizeof(REAL));
  for(i=0;i<m;i++) rpzwf[6][i]=f[i];
  //   copy the functions values go last, poles go first, residues go second
  residues_poles(m,rpzwf[4],rpzwf[5],rpzwf[6],rpzwf[0],rpzwf[1],rpzwf[2],rpzwf[3]);
  //
  REAL rswpi,rswpr;
  INT m1=m-1;
  rswpr=rpzwf[0][m1];
  rswpi=rpzwf[1][m1];
  for(i=m1;i>0;i--){
    rpzwf[0][i]=rpzwf[0][i-1];
    rpzwf[1][i]=rpzwf[1][i-1];
  }
  rpzwf[0][0]=rswpr;
  rpzwf[1][0]=rswpi;

  return (REAL )rmax;
}
/**********************************************************************/
/*!
 * \fn static REAL16 **set_f_values(REAL16 (*func)(REAL16 , REAL16 , REAL16 , REAL16 , REAL16 ),
 *                                  REAL16 s,REAL16 t,REAL16 alpha,REAL16 beta,
 *                                  INT *numval_in, 
 *                                  REAL xmin_in, REAL xmax_in, INT print_level)
 *
 *
 * \brief generates two long double arrays with at most numval_in
 *        entries storing z[i] and func(z[i],s,t,alpha,beta). If the
 *        value of the function is too large, then such z[i] is
 *        removed (only finite values are stored).
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
REAL16 **set_f_values(REAL16 (*func)(REAL16 , REAL16 , REAL16 , REAL16 , REAL16 ), \
		      REAL16 s,REAL16 t,REAL16 alpha,REAL16 beta,	\
		      INT *numval_in,					\
		      REAL xmin_in,					\
		      REAL xmax_in,					\
		      INT print_level)
{
  INT numval=*numval_in;
  REAL16  fj,zj,tol;
  REAL16 xmin=(REAL16 )xmin_in,xmax=(REAL16 )xmax_in;
  tol=powl(2e0,-52e0); // small tolerance: 100/tol is considered infinity:)
  if(numval<4)numval=4;
  REAL16 **zf=malloc(2*sizeof(REAL16 *));
  zf[0]=calloc(numval,sizeof(REAL16));
  zf[1]=calloc(numval,sizeof(REAL16));
  REAL16 *z=zf[0],*f=zf[1];
  INT j,k;
  k=0;j=0;
  while(1){
    zj=xmin+((REAL16 )j)/((REAL16 )(numval-1))*(xmax-xmin);
    fj=func(zj,s,t,alpha,beta);
    if(fabsl(fj)<1e2/tol){
      f[k]=fj;
      z[k]=zj;
      k++;
    }
    j++;
    if(j>=numval) break;
  }
  if(k<numval){
    if(print_level>1){
      fprintf(stdout,"\n%%%% WARNING: some values of f were too big; removing %d of the values in z[]\n", numval-k);
    }
    numval=k;
    z=realloc(z,numval*sizeof(REAL16));
    f=realloc(f,numval*sizeof(REAL16));
  }
  *numval_in=numval;
  return zf;
}
/*EOF EOF*/
