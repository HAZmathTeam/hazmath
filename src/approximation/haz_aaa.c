/*! \file src/approximation/haz_aaa.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note containing programs for the AAA algorithm for rational approximation
 *  \note 20201025 (ltz)
 */
#include "hazmath.h"
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
		    REAL16 *res,		\
		    REAL16 *z,			\
		    REAL16 *p,			\
		    REAL16 *f,			\
		    void *wrk)
{
  // wrk should have space m real16 and m integer. 
  //  given poles and right hand side, this finds the residues
  // Solves 
  INT i,j,k,im,m01=m-1;
  REAL16 apoljk,zj,fj;
  REAL16 *apol=(REAL16 *)calloc(m*m,sizeof(REAL16));
  REAL16 *piv = (REAL16 *)wrk;  
  INT *perm=(INT *)(wrk+m*sizeof(REAL16));
  for(j=0;j<m;j++){
    zj=z[j];
    for(k=0;k<m01;k++){
      apol[j*m+k]=1e00/(zj-p[k]);
      //      fprintf(stdout,"\nzCp(%d,%d)=%.18Le;",j+1,k+1,apol[j*m+k]);
    }
    apol[j*m+m01]=1e00;
    //    fprintf(stdout,"\nzCp(%d,%d)=%.18Le;",j+1,m01+1,apol[j*m+m01]);
  }
  //rhs
  for(j=0;j<m;j++) res[j]=f[j];
  solve_pivot_l(1,m,apol,res,perm,piv);
  /* for(i=0;i<m;i++){ */
  /*   //    fprintf(stdout,"\npiv(%d)=%.18Le;",i+1,apol[i*m+i]); */
  /*   fprintf(stdout,"\nperm(%d)=%d;",i+1,perm[i]+1); */
  /* } */
  free(apol);
  return;
}
/**********************************************************************/
/*!
 * \fn INT residues_poles(INT m,REAL *zd,REAL *wd,REAL *fd, 
 *                       REAL *resd,REAL *pold)
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
		   REAL *resd,				\
		   REAL *pold)
{
  INT m1=m+1,m01=m-1,mm1=m1*m1,i,j,k;
  REAL16 swp;
  // MEMORY
  void *wrk=(void *)calloc(7*m1*sizeof(REAL16)+			\
			   2*m1*m1*sizeof(REAL16)+		\
			   m1*m1*sizeof(REAL)+			\
			   m1*sizeof(INT),			\
			   sizeof(char));
  REAL16 *z=(REAL16 *)wrk;
  wrk+=m1*sizeof(REAL16);
  REAL16 *w=(REAL16 *)wrk;
  wrk+=m1*sizeof(REAL16);
  REAL16 *f=(REAL16 *)wrk;
  wrk+=m1*sizeof(REAL16);
  REAL16 *pol=(REAL16 *)wrk;
  wrk+=m1*sizeof(REAL16);
  REAL16 *res=(REAL16 *)wrk;
  wrk+=m1*sizeof(REAL16);
  REAL16 *ewrk=(REAL16 *)wrk;
  wrk+=mm1*sizeof(REAL16);
  REAL16 *col_id=(REAL16 *)wrk;
  wrk+=m1*sizeof(REAL16);
  REAL *einv=(REAL *)wrk; // only REAL because it goes into the
			  // eigenvalue computation and lapack has no
			  // long doubles
  wrk+=mm1*sizeof(REAL);
  REAL16 *pivot=(REAL16 *)wrk;
  wrk+=m1*sizeof(REAL16);
  INT *perm=(INT *)wrk;
  wrk+=m1*sizeof(INT);
  // END memory
  // grab the input params and convert them to higher precision
  for(j=0;j<m;j++){
    w[j]=(REAL16 )wd[j];
    z[j]=(REAL16 )zd[j];
    f[j]=(REAL16 )fd[j];
    /* fprintf(stdout,"\nzin(%d)=%.18Le;",j+1,z[j]);  fflush(stdout); */
    /* fprintf(stdout,"\nwin(%d)=%.18Le;",j+1,w[j]); fflush(stdout); */
    /* fprintf(stdout,"\nfin(%d)=%.18Le;",j+1,f[j]); fflush(stdout); */
  }
  // swap z[1:m]=z[m:-1:1] and same for w.
  for(i=0;i<((INT ) m/2);i++){
    swp=z[i]; z[i]=z[m-i-1]; z[m-i-1]=swp;
    swp=w[i]; w[i]=w[m-i-1]; w[m-i-1]=swp;
  }
  //
  //  form the matrix
  // ewrk
  memset(ewrk,0,m1*m1*sizeof(REAL16));
  for(i=0;i<m;i++){
    ewrk[i*m1+i]=z[i];
    ewrk[i*m1+m1-1]=1e00;
    ewrk[m*m1+i]=w[i];
  }  
  /* for(i=0;i<m1;i++){ */
  /*   for(j=0;j<m1;j++){ */
  /*     fprintf(stdout,"\newrk(%d,%d)=%.18Le;",i+1,j+1,ewrk[i*m1+j]); */
  /*   } */
  /* } */
  for(i=0;i<m;i++){
    memset(col_id,0,m1*sizeof(REAL16));
    col_id[i]=1e00;
    if(i==0)
      solve_pivot_l(1,m1,ewrk,col_id,perm,pivot);
    else
      solve_pivot_l(0,m1,ewrk,col_id,perm,pivot);      
    for(j=0;j<m;j++){
      einv[i*m+j]=(REAL )col_id[j];
    }
  }
  /* for(i=0;i<m;i++){ */
  /*   for(j=0;j<m;j++){ */
  /*     fprintf(stdout,"\neinv(%d,%d)=%.18e;",i+1,j+1,einv[i*m+j]); */
  /*   } */
  /* } */
  // compute eigenvalues:
  // use ewrk as a working space here; we dont need it anymore.
  REAL *wr=(REAL *)ewrk;
  REAL *wi=(REAL *)ewrk+m1;
  eiggeneral(m,einv,wr,wi);
  // put the  the one with minimal abs to be the last
  i=0;swp=(REAL16 )wr[m01];
  for(j=(m-2);j>=0;j--){
    if(fabs(wr[j])<fabs((REAL )swp)){
      // swap j and m-1
      swp=(REAL16 ) wr[j];
      wr[j]=wr[m01];
      wr[m01]=(REAL )swp;
    }
  }
  // swap pol, ignore the last.
  for(j=0;j<m01;j++){
    pol[j]=1e00/((REAL16 )wr[m01-j-1]);
  }
  // permute z and w back.
  for(i=0;i<((INT ) m/2);i++){
    swp=z[i]; z[i]=z[m-i-1]; z[m-i-1]=swp;
    swp=w[i]; w[i]=w[m-i-1]; w[m-i-1]=swp;
  }
  // get the residuals
  get_res(m,res,z,pol,f,(void *)pivot);
  /* for(j=0;j<m01;j++){ */
  /*   fprintf(stdout,"\npol000(%d)=%.18Le;",j+1,pol[j]); */
  /* } */
  /* fprintf(stdout,"\n"); */
  /* for(j=0;j<m;j++){ */
  /*   fprintf(stdout,"\nres000(%d)=%.18Le;",j+1,res[j]); */
  /* } */
  /* fprintf(stdout,"\n\n"); */
  ///////////////////double on output.
  //
  for(j=0;j<m;j++) resd[j]=(REAL )res[j];
  for(j=0;j<m01;j++) pold[j]=(REAL )pol[j];
  //
  free((void *)z);
  return 0;
}
/**********************************************************************/
/*!
 * \fn REAL get_cpzwf(REAL16 (*func)(REAL16 x, void *param),void
                     *param, REAL **wzpc,INT *mbig_in,INT *mmax_in,
                     INT *m_out, REAL xmin_in, REAL xmax_in, REAL
                     tolaaa)
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
REAL get_cpzwf(REAL16 (*func)(REAL16 x, void *param),	\
	      void *param,				\
	      REAL **cpzwf,				\
	      INT *mbig_in,				\
	      INT *mmax_in,				\
	      INT *m_out,				\
	      REAL xmin_in,				\
	      REAL xmax_in,				\
	      REAL tolaaa,
	      INT print_level)
{
  /*
   *  \param mbig is the total number of points; mmax_in is the max number of
   *  nodes taking part in the interpolation, mmax_in should be much
   *  smaller than mbig;
   *
   * \param m is the number of nodes in the final interpolation after
   *    tolerance tolaaa is achieved or mmax is reached.
   *   parameters for the function we are approximating.
   *   For example: s[2]*(x^s[0])+s[3]*(x**s[1])
   *   after a call to get_cpzwf, cpzwf[k] are pointers to the following
   *   arrays (each of m+1 elements but at most m are used:
   *
   *   cpzwf[0]->  residues(res[]);
   *   cpzwf[1]-> poles (pol[]);
   *   cpzwf[2]-> nodes (z[]);
   *   cpzwf[3]-> weights (w[]);
   *   cpzwf[4]-> function values (f[]) 
   * (also as last entry cpzwf[0] contains the free term c[m-1]; 
   *
   *   the rational approximation is:
   ****   r(z)=res[m-1] + \sum_{i=0}^{m-2} res[i]/(z-pol[i]); 
   *
  */
  ////////////
  // cpzwf must contain 5 pointers: cpzwf[0] are the residues and
  // cpzwf[1] are the poles. **cpzwf must be allocated earlier, but
  // cpzwf[k], k=0:4 are allocated here.  func() is the function we want to
  // approximate and param are the parameters it may depend on
  INT m,i,j,k,krow,kmax,mmax;
  REAL16 r,wzkj,fnum,fden,rm,rmax,swp;
  REAL16  fj,zj,tol,xmin=(REAL16 )xmin_in,xmax=(REAL16 )xmax_in;
  REAL smin=-1e20;
  INT info=-22;
  INT mbig=mbig_in[0];
  if(mbig<4)mbig=4;
  if(mmax_in[0]<2) mmax_in[0]=2;
  mmax=(INT )(mbig/2);
  if(mmax_in[0] < mmax)
    mmax=mmax_in[0];
  else 
    mmax_in[0]=mmax;
  tol=powl(2e0,-42e0);
  if(((REAL16 )tolaaa) > tol)
    tol=(REAL16 )tolaaa;
  INT mem=mbig*mmax*sizeof(REAL) + mmax*sizeof(REAL)+2*mbig*sizeof(REAL16);
  void *zfw_long=calloc(mem,sizeof(char));
  REAL16 *z=(REAL16 *)zfw_long;
  zfw_long+=mbig*sizeof(REAL16);
  REAL *wd=(REAL *)zfw_long;
  zfw_long+=(mmax)*sizeof(REAL);
  REAL16 *f=(REAL16 *)zfw_long;
  zfw_long+=mbig*sizeof(REAL16);
  REAL *x21d=(REAL *)zfw_long;
  zfw_long+=mbig*mmax*sizeof(REAL);// end of it
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  // form f; take the mean of f; form the matrix;
  // set up z;  
  r=0e0;
  j=0;
  k=0;
  while(1){
    zj=xmin+((REAL16 )j)/((REAL16 )(mbig-1))*(xmax-xmin);
    fj=func(zj,(void *)param);
    if(fabsl(fj)<1e2/tol){
      r+=fj;
      f[k]=fj;
      z[k]=zj;
      k++;
    }
    j++;
    if(j>=mbig) break;
  }
  if(k<mbig){
    if(print_level>1){
      fprintf(stdout,"\n%%%% WARNING: some values of f were too big; removing %d of the values in z[]\n", mbig-k);
    }
    mbig=k;
    if(mmax>((INT )mbig/2))mmax=(INT )(mbig/2);
  }
  r/=(REAL16 )mbig;
  /* fprintf(stdout,"\nr=%.18Le;\n", r); */
  /* for(j=0;j<mbig;j++){ */
  /*     fprintf(stdout,"\nz(%d)=%.20Le;f(%d)=%.18Le;",j+1,z[j],j+1,f[j]); */
  /* } */
  /* fprintf(stdout,"\nmmax=%d\n",mmax); */
  m=0;  
  while(1){
    rmax=-1.e00;
    kmax=-1;
    if(m==0){
      for(k=0;k<mbig;k++){
	rm=fabsl(r-f[j]);
	if(rm>rmax){
	  rmax=rm;
	  kmax=k;
	}
      }
    } else {
      for(k=m;k<mbig;k++){
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
    if(print_level>2)
      fprintf(stdout,"\n%%%%iter=%d | rmax=%.20Le at %7d;",m,rmax,kmax);
    if(rmax<tol || m>=(mmax-1)) break;
    swp=z[kmax];  z[kmax]=z[m];   z[m]=swp;
    swp=f[kmax];  f[kmax]=f[m];   f[m]=swp;
    m++;
    // form x21d (this needs improvement!:
    for(k=m;k<mbig;k++){
      krow=k-m;
      for(j=0;j<m;j++){
	x21d[krow*m+j]=(REAL )((f[k]-f[j])/(z[k]-z[j]));
      }
    }
    /* for(k=m;k<mbig;k++){ */
    /*   krow=k-m; */
    /*   for(j=0;j<m;j++){ */
    /* 	fprintf(stdout,"\nX21(%d,%d)=%.15e",krow+1,j+1,x21d[krow*m+j]); */
    /*   } */
    /* } */
    info=svdgeneral(mbig,m,x21d,&smin,wd);
    /* fprintf(stdout,"\nsmin=%e\n",smin); */
    /* k=0; */
    /* for(i=0;i<m;i++) {  */
    /*   fprintf(stdout,"\nk=%d, wd(%d)=%.12e",k,i,wd[i]); */
    /*   if(fabs(wd[i])<1e-13) { */
    /* 	continue; */
    /*   } else { */
    /* 	wd[k]=wd[i]; */
    /* 	z[k]=z[i]; */
    /* 	f[k]=f[i]; */
    /* 	k++;	 */
    /*   } */
    /* } */
    /* m=k; */
  }
  /* fprintf(stdout,"\n%%%%OLDSTUFF==================\n"); */
  /* for(i=0;i<m;i++) { */
  /*   fprintf(stdout,"\nzold(%d)=%.16Le,wold(%d)=%.16e",i+1,z[i],i+1,wd[i]); */
  /* } */
  /* fprintf(stdout,"\n"); */
  // ALL DONE:
  m_out[0]=m;
  mbig_in[0]=mbig;
  // tricky stuff: realloc removing all REAL16 to REAL:
  cpzwf[0]=(REAL *)z;
  // we use that real requires less than real16:) hope this works:
  for(i=0;i<m;i++) *(cpzwf[0]+i)=z[i];
  cpzwf[1]=cpzwf[0] + m + 1;
  memcpy(cpzwf[1],wd,m*sizeof(REAL));
  cpzwf[2]=cpzwf[1] + m + 1;
  for(i=0;i<m;i++) *(cpzwf[2]+i)=f[i];
  // reallocate:
  cpzwf[0]=realloc(cpzwf[0],5*(m+1)*sizeof(REAL));
  cpzwf[3]=cpzwf[2] + m + 1;
  cpzwf[4]=cpzwf[3] + m + 1;
  memcpy(cpzwf[4],cpzwf[2],(m+1)*sizeof(REAL));
  memcpy(cpzwf[3],cpzwf[1],(m+1)*sizeof(REAL));
  memcpy(cpzwf[2],cpzwf[0],(m+1)*sizeof(REAL));
  //   copy the functions values go last, poles go first, residues go second
  residues_poles(m,cpzwf[2],cpzwf[3],cpzwf[4],cpzwf[0],cpzwf[1]);
  REAL rswp;
  INT m1=m-1;
  rswp=*(cpzwf[0]+m1);
  for(i=m1;i>0;i--)
    *(cpzwf[0]+i)=*(cpzwf[0]+i-1);
  *(cpzwf[0])=rswp;
  return (REAL )rmax;
}
/*EOF*/
