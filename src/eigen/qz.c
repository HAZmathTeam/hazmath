/*! \file src/fem/qz.c.c
*
* \brief Implements the the QZ algorithm for generalized eigenvalue problem.
*
*  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 2/17/15.
*  Copyright 2015__HAZMATH__. All rights reserved.
*
* \note The algorithm is described in: MR0345399 (49 10135) 65F15
*       Moler, C. B. ; Stewart, G. W. An algorithm for generalized
*       matrix eigenvalue problems.  SIAM J. Numer. Anal. 10 (1973),
*       241–256.
*       ForTran routines are found in: Moler, C. B. ;
*       Stewart, G. W.  An algorithm for the generalized matrix
*       eigenvalue problem Ax=\lambda Bx.  UT Austin tech report
*       N00014-67-A-0181-0023 issued jointly with Stanford techreport
*       STAN-CS-232-71
*
* \note ALL is in long double. 
*
* \note Ludmil 20210426
*/
/************************************************************************************/
#include "hazmath.h"
/************************************************************************************/
#ifndef SQRT
#define SQRT sqrtl
#endif
#ifndef FABS
#define FABS fabsl
#endif
/***************************************************************************/
/*
 * \fn static void del_row_col(INT ndp1,INT np1,REAL16 **a_io)
 *
 * \brief removes the first row and column of a matrix a(np1,np1). the
 *        result overwrites a and a is realloc-ed so on return it has
 *        the correct size.
 *
 * \note The algorithm is described in: MR0345399 (49 10135) 65F15
 *       Moler, C. B. ; Stewart, G. W. An algorithm for generalized
 *       matrix eigenvalue problems.  SIAM J. Numer. Anal. 10 (1973),
 *       241–256.
 *       ForTran routines are found in: Moler, C. B. ;
 *       Stewart, G. W.  An algorithm for the generalized matrix
 *       eigenvalue problem Ax=\lambda Bx.  UT Austin tech report
 *       N00014-67-A-0181-0023 issued jointly with Stanford techreport
 *       STAN-CS-232-71
 * \note Ludmil 20210426
 */
/*********************************************************************************/
static void del_row_col(INT ndp1,INT np1,REAL16 **a_io)
/*********************************************************************************/
{
  INT i,j,ip1,jp1,n=np1-1,nd=ndp1-1;
  REAL16 *a=a_io[0];
  for(i=0;i<nd;++i){
    ip1=i+1;
    for(j=0;j<n;++j){
      jp1=j+1;
      //      fprintf(stdout,"\n%% (%d,%d):%d <-- %d(%d,%d)",i,j,i*nd+j,ip1*np1+jp1,ip1,jp1);
      a[i*n+j]=a[ip1*np1+jp1];
    }
  }
  //  fprintf(stdout,"\n");
  a_io[0]=realloc(a,nd*n*sizeof(REAL16));
  return;
}
/*********************************************************************************/
/*
 * \fn static void add_row_col(INT nd,INT n,REAL16 **a_io)
 *
 * \brief ADDS a first row and column of zeros to a matrix a(nd by
 * n). the result keeps a as the a22 block; this is a hack to convert
 * from fortran to C.
 *
 * \note The algorithm is described in: MR0345399 (49 10135) 65F15
 *       Moler, C. B. ; Stewart, G. W. An algorithm for generalized
 *       matrix eigenvalue problems.  SIAM J. Numer. Anal. 10 (1973),
 *       241–256.
 *       ForTran routines are found in: Moler, C. B. ;
 *       Stewart, G. W.  An algorithm for the generalized matrix
 *       eigenvalue problem Ax=\lambda Bx.  UT Austin tech report
 *       N00014-67-A-0181-0023 issued jointly with Stanford techreport
 *       STAN-CS-232-71
 * \note Ludmil 20210426
   */
static void add_row_col(INT nd,INT n,REAL16 **a_io)
/*********************************************************************************/
{
  INT i,j,ip1,jp1,np1=n+1,ndp1=nd+1;
  a_io[0]=realloc(a_io[0],ndp1*np1*sizeof(REAL16));
  //
  REAL16 *a=a_io[0];
  for(i=(nd-1);i>=0;--i){
    ip1=i+1;
    for(j=(n-1);j>=0;--j){
      jp1=j+1;
      //      fprintf(stdout,"\n%%11:: (%d,%d):%d --> %d(%d,%d)",ip1,jp1,ip1*np1+jp1,i*nd+j,i,j);
      a[ip1*np1+jp1]=a[i*n+j];
    }
  }
  //  fprintf(stdout,"\n");
  for(i=0;i<ndp1;++i)
    a[i*np1+0]=0e0;
  for(j=0;j<np1;++j)
    a[0+j]=0e0;
  return;
}
/**********************************************************************/
/*
 * \fn static void househ3(REAL16 a1,REAL16 a2,REAL16 a3,
 *      		   REAL16 *u1, REAL16 *u2, REAL16 *u3,
 *		           REAL16 *v1, REAL16 *v2, REAL16 *v3)
 *
 * \brief  huseholder (a1,as2,a3)-->(*,0,0);
 *
 * \note The algorithm is described in: MR0345399 (49 10135) 65F15
 *       Moler, C. B. ; Stewart, G. W. An algorithm for generalized
 *       matrix eigenvalue problems.  SIAM J. Numer. Anal. 10 (1973),
 *       241–256.
 *       ForTran routines are found in: Moler, C. B. ;
 *       Stewart, G. W.  An algorithm for the generalized matrix
 *       eigenvalue problem Ax=\lambda Bx.  UT Austin tech report
 *       N00014-67-A-0181-0023 issued jointly with Stanford techreport
 *       STAN-CS-232-71
 * \note Ludmil 20210426
 */
/*****************************************************************/
static void househ3(REAL16 a1,REAL16 a2,REAL16 a3,		\
  		    REAL16 *u1, REAL16 *u2, REAL16 *u3,		\
		    REAL16 *v1, REAL16 *v2, REAL16 *v3)
{
  *u1=-1e19;  *u2=-1e19;  *u3=-1e19;
  *v1=-1e19;  *v2=-1e19;  *v3=-1e19;
  REAL16 r,s;
  if((a2==0e0)&&(a3==0e0)) {
    *u1=0e0;
  } else {
    s =FABS(a1) + FABS(a2) + FABS(a3);
    *u1 = a1/s;
    *u2 = a2/s;
    *u3 = a3/s;
    r =SQRT((*u1)*(*u1) + (*u2)*(*u2) + (*u3)*(*u3));
    if ((*u1)<0e0) r=-r;
    *v1 = -((*u1) + r)/r;
    *v2 = -(*u2)/r;
    *v3 = -(*u3)/r;
    *u1 = 1e0;
    *u2 = (*v2)/(*v1);
    *u3 = (*v3)/(*v1);
  }
  return;
}
/****************************************************************/
/*
 * \fn static void househ2(REAL16 a1,REAL16 a2,
 *                         REAL16 *u1, REAL16 *u2,
 *		           REAL16 *v1, REAL16 *v2)
 *
 * \brief   huseholder (a1,a2)-->(*,0);
 *
 * \note The algorithm is described in: MR0345399 (49 10135) 65F15
 *       Moler, C. B. ; Stewart, G. W. An algorithm for generalized
 *       matrix eigenvalue problems.  SIAM J. Numer. Anal. 10 (1973),
 *       241–256.
 *       ForTran routines are found in: Moler, C. B. ;
 *       Stewart, G. W.  An algorithm for the generalized matrix
 *       eigenvalue problem Ax=\lambda Bx.  UT Austin tech report
 *       N00014-67-A-0181-0023 issued jointly with Stanford techreport
 *       STAN-CS-232-71
 *
 * \note Ludmil 20210426
 *
*/
static void househ2(REAL16 a1,REAL16 a2,		\
		    REAL16 *u1, REAL16 *u2,		\
		    REAL16 *v1, REAL16 *v2)
{
  *u1=-1e18;  *u2=-1e18;
  *v1=-1e18;  *v2=-1e18;
  REAL16 r,s;
  if(a2==0e0) {
    *u1=0e0;
  } else {
    s =FABS(a1) + FABS(a2);
    *u1 = a1/s;
    *u2 = a2/s;
    r =SQRT((*u1)*(*u1) + (*u2)*(*u2));
    if ((*u1)<0e0) r=-r;
    *v1 = -((*u1) + r)/r;
    *v2 = -(*u2)/r;
    *u1 = 1e0;
    *u2 = (*v2)/(*v1);
  }
  return;
}
/****************************************************************/
/*
 * \fn static void chouseh2(REAL16 a1r,REAL16 a1i,
 *		     REAL16 a2r,REAL16 a2i,
 *		     REAL16 *c, REAL16 *sr, REAL16 *si)
 *
 * \brief complex huseholder (a1,a2,a3)-->(*,0,0);
 *
 * \note The algorithm is described in: MR0345399 (49 10135) 65F15
 *       Moler, C. B. ; Stewart, G. W. An algorithm for generalized
 *       matrix eigenvalue problems.  SIAM J. Numer. Anal. 10 (1973),
 *       241–256.
 *       ForTran routines are found in: Moler, C. B. ;
 *       Stewart, G. W.  An algorithm for the generalized matrix
 *       eigenvalue problem Ax=\lambda Bx.  UT Austin tech report
 *       N00014-67-A-0181-0023 issued jointly with Stanford techreport
 *       STAN-CS-232-71
 *
 * \note Ludmil 20210426
 */
static void chouseh2(REAL16 a1r,REAL16 a1i,		\
		     REAL16 a2r,REAL16 a2i,			\
		     REAL16 *c, REAL16 *sr, REAL16 *si)
{
  *c=-1e17;  *sr=-1e17;  *si=-1e17;
  REAL16 r;
  if((a2r==0e0)&&(a2i==0e0)) {
    *c=1e0;
    *sr=0e0;
    *si=0e0;
  } else if((a1r==0e0)&&(a1i==0e0)) {
    *c=0e0;
    *sr=1e0;
    *si=0e0;
  } else {
    r =SQRT(a1r*a1r + a1i*a1i);
    *c=r;
    *sr=(a1r*a2r + a1i*a2i)/r;
    *si=(a1r*a2i - a1i*a2r)/r;
    r=SQRT((*c)*(*c) + (*sr)*(*sr) + (*si)*(*si));
    *c/=r;  *sr/=r;    *si/=r;
  }
  return;
}
/***************************************************************************/
/*
 * \fn void xdivy(const REAL16 xr,const REAL16 xi,
 *		  const REAL16 yr,const REAL16 yi,
 *		  REAL16 *zr,REAL16 *zi)
 *
 * \brief complex division: compute (xr+i*xi)/(yr+i*yi);
 *
 * \note The algorithm is described in: MR0345399 (49 10135) 65F15
 *       Moler, C. B. ; Stewart, G. W. An algorithm for generalized
 *       matrix eigenvalue problems.  SIAM J. Numer. Anal. 10 (1973),
 *       241–256.
 *       ForTran routines are found in: Moler, C. B. ;
 *       Stewart, G. W.  An algorithm for the generalized matrix
 *       eigenvalue problem Ax=\lambda Bx.  UT Austin tech report
 *       N00014-67-A-0181-0023 issued jointly with Stanford techreport
 *       STAN-CS-232-71
 *
 * \note Ludmil 20210426
 */
void xdivy(const REAL16 xr,const REAL16 xi,	\
	   const REAL16 yr,const REAL16 yi,	\
	   REAL16 *zr,REAL16 *zi)
{
  REAL16 wr,wi,vr,vi,d;
  if(FABS(yr)>=FABS(yi)){
    wr=xr/yr;
    wi=xi/yr;
    vi=yi/yr;
    d=1.+vi*vi;
    *zr=(wr+wi*vi)/d;
    *zi=(wi-wr*vi)/d;
  } else {
    wr=xr/yi;
    wi=xi/yi;
    vr=yr/yi;
    d=1.+vr*vr;
    *zr=(wr*vr+wi)/d;
    *zi=(wi*vr-wr)/d;
  }
  return;
}
/****************************************************************************/
/*
 * \fn void qzhes(INT nd,INT n,REAL16 *a,REAL16 *b,const INT wantx, REAL16 *x)
 *
 * \brief transforms a in Hessenberg and b in upper triangular form.
 *
 * \note The algorithm is described in: MR0345399 (49 10135) 65F15
 *       Moler, C. B. ; Stewart, G. W. An algorithm for generalized
 *       matrix eigenvalue problems.  SIAM J. Numer. Anal. 10 (1973),
 *       241–256.
 *       ForTran routines are found in: Moler, C. B. ;
 *       Stewart, G. W.  An algorithm for the generalized matrix
 *       eigenvalue problem Ax=\lambda Bx.  UT Austin tech report
 *       N00014-67-A-0181-0023 issued jointly with Stanford techreport
 *       STAN-CS-232-71
 *
 * \note Ludmil 20210426
 */
/****************************************************************************/
void qzhes(INT nd,INT n,REAL16 *a,REAL16 *b,const INT wantx, REAL16 *x)
{
  INT i,j,k,l,l1,nm1,nm2,nk1,lb;
  REAL16 r,s,t,u1,u2,v1,v2,rho;
//initialize x, used to save transformations
  if (wantx){
    memset(x,0,nd*n*sizeof(REAL16));
    for(i=0;i<nd;++i)
      x[i*n+i]=1e0;
  }
  //reduce B to upper triangular...
  nm1=n-1;
  for(l=0;l<nm1;++l){
    l1=l+1;
    s=0e0;
    for(i=l1;i<n;++i){
      if(FABS(b[i*n+l])>s)
	s=FABS(b[i*n+l]);
    }
    if (s==0e0)
      continue;
    if(FABS(b[l*n+l])>s)
      s=FABS(b[l*n+l]);
    r=0e0;
    for(i=l;i<n;++i){
      b[i*n+l]=b[i*n+l]/s;
      r+=b[i*n+l]*b[i*n+l];
    }
    r=SQRT(r);
    if(b[l*n+l]<0e0)
      r=-r;
    b[l*n+l]+=r;
    rho=b[l*n+l]*r;
    for(j=l1;j<n;++j){
      t=0e0;
      for(i=l;i<n;++i){
	t+=b[i*n+l]*b[i*n+j];
      }
      t=-t/rho;
      for(i=l;i<n;++i){
	b[i*n+j]+=t*b[i*n+l];
      }
    }
    for(j=0;j<n;j++){
      t = 0e0;
      for(i=l;i<n;++i){
	t+= b[i*n+l]*a[i*n+j];
      }
      t=-t/rho;
      for(i=l;i<n;++i){
	a[i*n+j] += t*b[i*n+l];
      }
    }
    b[l*n+l]=-s*r;
    for(i=l1;i<n;++i){
      b[i*n+l] = 0e0;
    }
  }// 100
  //
  //reduce a to upper hessenberg, keep b triangular
  if (n<2)
    return; // go to 170
  nm2=n-2;
  for(k=0;k<nm2;++k){
    //k1 = k+1;
    nk1=nm2-k;// was n-k-1;
    for(lb=0;lb<nk1;++lb){
      l=nm2-lb;//was n-lb
      l1=l+1;
      //      fprintf(stdout,"\nn=%d,k=%d,k1=%d,lb=%d,l=%d,l1=%d",n,k,k1,lb,l,l1);
      /************************************************/
      househ2(a[l*n+k],a[l1*n+k],&u1,&u2,&v1,&v2);
      /************************************************/
      if(u1==1.){//do not go to 125
	for(j=k;j<n;++j){
	  t = a[l*n+j] + u2*a[l1*n+j];
	  a[l*n+j]+=t*v1;
	  a[l1*n+j]+=t*v2;
	}
	a[l1*n+k] = 0e0;
	for(j=l;j<n;++j){
	  t = b[l*n+j] + u2*b[l1*n+j];
	  b[l*n+j]+=t*v1;
	  b[l1*n+j]+=t*v2;
	}
      }
      //      fprintf(stdout,"\nmat='a';u1a=%.16e;u2a=%.16e;v1a=%.16e;v2a=%.16e;",u1,u2,v1,v2);fflush(stdout);
      //      fprintf(stdout,"\n%%Transforming b(%d,%d)=%.16e and
      //      b(%d,%d)=%.16e",l1,l1,b[l1*n+l1],l1,l,b[l1*n+l]);fflush(stdout);
      /************************************************/
      househ2(b[l1*n+l1],b[l1*n+l],&u1,&u2,&v1,&v2);
      /************************************************/
      //      fprintf(stdout,"\nmat='b';u1b=%.16e;u2b=%.16e;v1b=%.16e;v2b=%.16e;",u1,u2,v1,v2);fflush(stdout);
      if (u1!=1e0){
	continue;
      }
      //      fprintf(stdout,"\nmat='b';u1b=%.16e;u2b=%.16e;v1b=%.16e;v2b=%.16e;",u1,u2,v1,v2);fflush(stdout);
      for(i=0;i<=l1;++i){
	t = b[i*n+l1] + u2*b[i*n+l];
	b[i*n+l1]+= t*v1;
	b[i*n+l]+= t*v2;
      }
      b[l1*n+l] = 0e0;
      //      fprintf(stdout,"Qb=inv(By)*Bz;Qb=Qb.*(abs(Qb)>1e-12);chkb=norm(Qb*transpose(Qb)-eye(%d),1)",n);
      for(i=0;i<n;++i){
	t = a[i*n+l1] + u2*a[i*n+l];
	a[i*n+l1]+= t*v1;
	a[i*n+l]+= t*v2;
      }
//      fprintf(stdout,"Qa=inv(Ay)*Az;Qa=Qa.*(abs(Qa)>1e-12);chka=norm(Qa*transpose(Qa)-eye(%d),1)",n);
////      exit(22);
      if(wantx){
	for(i=0;i<n;++i){
	  t = x[i*n+l1] + u2*b[i*n+l];
	  x[i*n+l1]+= t*v1;
	  x[i*n+l]+= t*v2;
	}
      }
    }
  }
  return;
}
/****************************************************************************/
/*
 * \fn void qzit(INT nd,INT n,
 *	  REAL16 **a_io,REAL16 **b_io,
 *	  const REAL16 eps,
 *	  REAL16 *epsa,REAL16 *epsb,
 *	  INT **iter_io,
 *	  INT wantx,
 *	  REAL16 **x_io)
 *
 * \brief iteration to convert a with 2 by 2 blocks on the diagonal and keep b upper triangular
 *
 * \note The algorithm is described in: MR0345399 (49 10135) 65F15
 *       Moler, C. B. ; Stewart, G. W. An algorithm for generalized
 *       matrix eigenvalue problems.  SIAM J. Numer. Anal. 10 (1973),
 *       241–256.
 *       ForTran routines are found in: Moler, C. B. ;
 *       Stewart, G. W.  An algorithm for the generalized matrix
 *       eigenvalue problem Ax=\lambda Bx.  UT Austin tech report
 *       N00014-67-A-0181-0023 issued jointly with Stanford techreport
 *       STAN-CS-232-71
 *
 */
/********************************************************************************/
void qzit(INT nd,INT n,			\
	  REAL16 **a_io,REAL16 **b_io,			\
	  const REAL16 eps,				\
	  REAL16 *epsa,REAL16 *epsb,		\
	  INT **iter_io,				\
	  INT wantx,				\
	  REAL16 **x_io)
{
  //initialize iter and compute epsa,epsb
  REAL16 anorm,bnorm,t,ani,bni;
  REAL16 u1,u2,u3,v1,v2,v3,cc,old1,old2;
  REAL16 a10,a20,a30,a11,a12,a21,a22,a33,a34,a43,a44,	\
    b11,b12,b22,b33,b34,b44;
  INT i,j,k,l,l1,m,k1,k2,k3,mid,m0rn,l0r1,m1,km1,lb;
  //  INT mm,next;
  INT ndp1=nd+1,np1=n+1; //,ip1,jp1;
  add_row_col(nd,n,a_io);
  add_row_col(nd,n,b_io);
  add_row_col(nd,n,x_io);
  REAL16 *a=a_io[0], *b=b_io[0], *x=x_io[0];
  // set it up
  INT *iter=realloc(iter_io[0], np1*sizeof(INT));
  for(i=(nd-1);i>=0;--i)
    iter[i+1]=iter[i];
  anorm = 0e0;
  bnorm = 0e0;
  for(i=1;i<=n;++i){
    iter[i]=0;
    ani=0.;
    if (i!=0)
      ani =FABS(a[i*np1+(i-1)]);
    bni = 0e0;
    for(j=0;j<n;++j){
      ani += FABS(a[i*np1+j]);
      bni += FABS(b[i*np1+j]);
    }
    if (ani>anorm)
      anorm = ani;
    if (bni>bnorm)
      bnorm=bni;
  }
  (*epsa) = eps*anorm;
  (*epsb) = eps*bnorm;
  // fprintf(stdout,"\nanorm=%Le,bnorm=%Le,epsa=%Le,epsb=%Le\n",anorm,bnorm,*epsa,*epsb);
  //
  m=n;
 l000:
  if(m<=2) goto l190;
  for(lb=1;lb<=m;++lb){
    l = m+1-lb;
    if(l==1) goto l060;
    if(FABS(a[l*ndp1 + (l-1)])<=(*epsa)) goto l030;
  }
 l030:
  a[l*np1+(l-1)]=0e0;
  if (l<m-1) goto l060;
  m=l-1;
  goto l000;
  //check for small top of b
 l060:
  if (FABS(b[l*np1+l])>(*epsb)) goto l100;
  b[l*np1+l] =0e0;
  l1=l+1;
  househ2(a[l*np1+l],a[l1*np1+l],&u1,&u2,&v1,&v2);
  if (u1!=1e0) goto l080;
  for(j=l;j <=n;++j){
    t = a[l*np1+j] + u2*a[l1*np1+j];
    a[l*np1+j]+=t*v1;
    a[l1*np1+j]+=t*v2;
    t = b[l*np1+j] + u2*b[l1*np1+j];
    b[l*np1+j]+=t*v1;
    b[l1*np1+j]+=t*v2;
  }
 l080:
  l = l1;
  goto l030;
 l100:
  m1 = m - 1; // here m is .ge. 2
  l1 = l + 1;// l is at least 0
  cc= 0.75e0;
  iter[m]++;
  if(iter[m]==1) goto l105;
  if(a[m*np1+(m-1)]<cc*old1) goto l105;
  if(a[(m-1)*np1+(m-2)]<cc*old2) goto l105;
  if(iter[m]==10) goto l110;
  if(iter[m]>30) goto l180;
 l105:
  b11=b[l*np1+l];
  b22=b[l1*np1+l1];
  if (FABS(b22)<(*epsb)) b22 = (*epsb);
  b33=b[m1*np1+m1];
  if (FABS(b33)<(*epsb)) b33 = (*epsb);
  b44=b[m*np1+m];
  if (FABS(b44)<(*epsb)) b44 = (*epsb);
  a11 =a[l*np1+l]/b11;
  a12 =a[l*np1+l1]/b22;
  a21 =a[l1*np1+l]/b11;
  a22 =a[l1*np1+l1]/b22;
  //
  a33 =a[m1*np1+m1]/b33;
  a34 =a[m1*np1+m]/b44;
  a43 =a[m*np1+m1]/b33;
  a44 =a[m*np1+m]/b44;
  //
  b12 =b[l*np1+l1]/b22;
  b34 =b[m1*np1+m]/b44;
  a10=( (a33-a11)*(a44-a11) - a34*a43 + a43*b34*a11)/a21 + a12-a11*b12;
  a20= (a22-a11-a21*b12) -(a33-a11)-(a44-a11)+a43*b34;
  a30=a[(l+2)*np1+l1]/b22;
  goto l115;
  //      fprintf(stdout,"\na10=%.12e,a20=%.12e,a30=%.12e,a21=%f",a10,a20,a30,a21);fflush(stdout);
 l110:
  a10=0e0;
  a20=0e0;
  a30=1.1605e0;
  //
 l115:
  old1=FABS(a[m*np1+(m-1)]);
  old2=FABS(a[(m-1)*np1+(m-2)]);
  if(wantx){
    l0r1=1;
    m0rn=n;
  } else {
    l0r1=l;
    m0rn=m;
  }
  for(k=l;k<=m1;++k){
    mid=(INT )(k!=m1);
    k1=k+1;
    k2=k+2;
    k3=k+3;
    if(k3>m) k3=m;
    km1=k-1;
    if(km1<l) km1=l;
    if(k==l) househ3(a10,a20,a30,&u1,&u2,&u3,&v1,&v2,&v3);
    if(k>l && k<m1){
      househ3(a[k*np1+km1],a[k1*np1+km1],a[k2*np1+km1],	\
	      &u1,&u2,&u3,&v1,&v2,&v3);
    }
    if(k==m1){
      househ2(a[k*np1+km1],a[k1*np1+km1],		\
	      &u1,&u2,&v1,&v2);
    }
    if(u1!=1e0) goto l125;
    for(j=km1;j<=m0rn;++j){
      //	    fprintf(stdout,"\n%%3333j=%dk=%d,km1=%d,k1=%d,k2=%d,k3=%d",j,k,km1,k1,k2,k3);
      t=a[k*np1+j]+u2*a[k1*np1+j];
      if(mid) t+=u3*a[k2*np1+j];
      a[k*np1+j]+=t*v1;
      a[k1*np1+j]+=t*v2;
      if(mid) a[k2*np1+j]+=t*v3;
      //
      t=b[k*np1+j]+u2*b[k1*np1+j];
      if(mid) t+=u3*b[k2*np1+j];
      b[k*np1+j]+=t*v1;
      b[k1*np1+j]+=t*v2;
      if(mid) b[k2*np1+j]+=t*v3;
    }
    if(k==l) goto l125;
    a[k1*np1+(k-1)]=0e0;
    if(mid) a[k2*np1+(k-1)]=0e0;
  l125:
    if(k==m1) goto l140;
    househ3(b[k2*np1+k2],b[k2*np1+k1],b[k2*np1+k],	\
	    &u1,&u2,&u3,&v1,&v2,&v3);
    if(u1!=1e0) goto l140;
    for(i=l0r1;i<=k3;++i){
      t=a[i*np1+k2]+u2*a[i*np1+k1]+u3*a[i*np1+k];
      a[i*np1+k2]+=t*v1;
      a[i*np1+k1]+=t*v2;
      a[i*np1+k]+=t*v3;
      //
      t=b[i*np1+k2]+u2*b[i*np1+k1]+u3*b[i*np1+k];
      b[i*np1+k2]+=t*v1;
      b[i*np1+k1]+=t*v2;
      b[i*np1+k]+=t*v3;
    }
    b[k2*np1+k]=0e0;
    b[k2*np1+k1]=0e0;
    if(wantx){
      for(i=0;i<=n;++i){
	t=x[i*np1+k2]+u2*x[i*np1+k1]+u3*x[i*np1+k];
	x[i*np1+k2]+=t*v1;
	x[i*np1+k1]+=t*v2;
	x[i*np1+k]+=t*v3;
      }
    }
  l140:
    househ2(b[k1*np1+k1],b[k1*np1+k],&u1,&u2,&v1,&v2);
    if(u1!=1.) continue; //goto l160;
    for(i=l0r1;i<=k3;++i){
      t=a[i*np1+k1]+u2*a[i*np1+k];
      a[i*np1+k1]+=t*v1;
      a[i*np1+k]+=t*v2;
      //
      t=b[i*np1+k1]+u2*b[i*np1+k];
      b[i*np1+k1]+=t*v1;
      b[i*np1+k]+=t*v2;
      //
    }
    b[k1*np1+k]=0.;
    if(wantx){
      for(i=1;i<=n;++i){
	t=x[i*np1+k1]+u2*x[i*np1+k];
	x[i*np1+k1]+=t*v1;
	x[i*np1+k]+=t*v2;
      }
    }// wantx
  }
  goto l000;
 l180:
  for(i=0;i<m;++i){
    iter[i]=-1;
  }
 l190:
  del_row_col(ndp1,np1,&a);
  del_row_col(ndp1,np1,&b);
  del_row_col(ndp1,np1,&x);
  a_io[0]=a;  b_io[0]=b;  x_io[0]=x;
  for(i=0;i<n;++i)
    iter[i]=iter[i+1];
  iter_io[0]=realloc(iter,n*sizeof(INT));
  return;
}
/*********************************************************************************/
/*
 * \fn void qzval(INT nd,INT n,
 *	   REAL16 **a_io,REAL16 **b_io,
 *	   REAL16 epsa,REAL16 epsb,
 *	   REAL16 **alphar_io,REAL16 **alphai_io,REAL16 **beta_io,
 *	   INT wantx,REAL16 **x_io)
 *
 * \brief computes eigenvalues of the generalized eigenvalue problems
 *        a matrix that is block diagonal with 2x2 blocks on the
 *        diagonal and B is upper triangular.
 *
 * \note The algorithm is described in: MR0345399 (49 10135) 65F15
 *       Moler, C. B. ; Stewart, G. W. An algorithm for generalized
 *       matrix eigenvalue problems.  SIAM J. Numer. Anal. 10 (1973),
 *       241–256.
 *       ForTran routines are found in: Moler, C. B. ;
 *       Stewart, G. W.  An algorithm for the generalized matrix
 *       eigenvalue problem Ax=\lambda Bx.  UT Austin tech report
 *       N00014-67-A-0181-0023 issued jointly with Stanford techreport
 *       STAN-CS-232-71
 *
 */
/********************************************************************************/
void qzval(INT nd,INT n,
	   REAL16 **a_io,REAL16 **b_io,
	   REAL16 epsa,REAL16 epsb,
	   REAL16 **alphar_io,REAL16 **alphai_io,REAL16 **beta_io,
	   INT wantx,REAL16 **x_io)
{
  INT flip;
  //find eigenvalues of quasi-triangular matrices
  INT i,j,m,l,ndp1=nd+1,np1=n+1;//,ip1,jp1;
  REAL16  a11,a12,a21,a22,b11,b12,b22;
  REAL16  a11r,a12r,a21r,a22r;//,b11r,b12r,b22r;
  REAL16  a11i,a12i,a21i,a22i; // ,b11i,b12i,b22i;
  REAL16 u1,u2,v1,v2;
  REAL16 r,t,c,d,e,an,bn;
  REAL16 szr,szi,sqr,sqi,ssr,ssi,cz,cq;
  REAL16 ti,tr,er,ei,bdi,bdr;
  add_row_col(nd,n,a_io);
  add_row_col(nd,n,b_io);
  add_row_col(nd,n,x_io);
  // set it up values
  REAL16 *alphar=realloc(alphar_io[0], ndp1*sizeof(REAL16));
  REAL16 *alphai=realloc(alphai_io[0], ndp1*sizeof(REAL16));
  REAL16 *beta=realloc(beta_io[0], ndp1*sizeof(REAL16));
  for(i=(nd-1);i>=0;--i){
    alphar[i+1] = alphar[i];
    alphai[i+1] = alphai[i];
    beta[i+1]  = beta[i];
  }
  REAL16 *a=a_io[0], *b=b_io[0], *x=x_io[0];
  /***********************************************************************************/
  m = n;
  //  400 continue
 l200:
  if(m==1) goto l210;
  if(a[m*np1+(m-1)] != 0e0) goto l220;
  // one-by-one submatrix; one real root
 l210:
  alphar[m] = a[m*np1+m];
  beta[m]  = b[m*np1+m];
  alphai[m]=0e0;
  m--;
  goto l290;
  //two-by-two submatrix
 l220:
  l=m-1;
  if(FABS(b[l*np1+l]) > epsb) goto l225;
  b[l*np1+l]=0e0;
  househ2(a[l*np1+l],a[m*np1+l],&u1,&u2,&v1,&v2);
  goto l260;
 l225:
  if(FABS(b[m*np1+m]) > epsb) goto l230;
  b[m*np1+m]=0e0;
  househ2(a[m*np1+m],a[m*np1+l],&u1,&u2,&v1,&v2);
  bn=0e0;
  goto l235;
  //
 l230:
  an=FABS(a[l*np1+l]) + FABS(a[l*np1+m]) +FABS(a[m*np1+l]) + FABS(a[m*np1+m]);
  bn=FABS(b[l*np1+l]) + FABS(b[l*np1+m]) + FABS(b[m*np1+m]);
  a11=a[l*np1+l]/an;
  a12=a[l*np1+m]/an;
  a21=a[m*np1+l]/an;
  a22=a[m*np1+m]/an;
  //
  b11=b[l*np1+l]/bn;
  b12=b[l*np1+m]/bn;
  b22=b[m*np1+m]/bn;
  c=0.5e0*(a11*b22 + a22*b11-a21*b12);
  d = 0.5e0*(a22*b11-a11*b22-a21*b12);
  d = d*d + a21*b22*(a12*b11-a11*b12);
  if(d<0e0) goto l280;
  // two real roots:
  if(c<0e0)
    e=(c-SQRT(d))/(b11*b22);
  else
    e=(c+SQRT(d))/(b11*b22);
  a11-=e*b11;
  a12-=e*b12;
  a22-=e*b22;
  // reversed;
  flip=(INT )( (FABS(a11)+FABS(a12)) < (FABS(a21)+FABS(a22)) );
  if(flip)
    househ2(a22,a21,&u1,&u2,&v1,&v2);
  else
    househ2(a12,a11,&u1,&u2,&v1,&v2);
 l235:
  if(u1!=1e0) goto l250;
  for(i=1;i<=m;++i){
    t=a[i*np1+m]+u2*a[i*np1+l];
    a[i*np1+m]+=t*v1;
    a[i*np1+l]+=t*v2;
    //
    t=b[i*np1+m]+u2*b[i*np1+l];
    b[i*np1+m]+=t*v1;
    b[i*np1+l]+=t*v2;
  }
  if(wantx){
    for(i=1;i<=n;++i){
      t=x[i*np1+m]+u2*x[i*np1+l];
      x[i*np1+m]+=t*v1;
      x[i*np1+l]+=t*v2;
    }
  }
 l250:
  if(bn==0e0) goto l275;
  flip=(INT )(an < (FABS(e)*bn));
  if(flip)
    househ2(a[l*np1+l],a[m*np1+l],&u1,&u2,&v1,&v2);
  else
    househ2(b[l*np1+l],b[m*np1+l],&u1,&u2,&v1,&v2);
 l260:
  if(u1!=1e0) goto l275;
  for(j=l;j<=n;++j){
    t=a[l*np1+j]+u2*a[m*np1+j];
    a[l*np1+j]+=t*v1;
    a[m*np1+j]+=t*v2;
    //
    t=b[l*np1+j]+u2*b[m*np1+j];
    b[l*np1+j]+=t*v1;
    b[m*np1+j]+=t*v2;
  }
  //
 l275:
  a[m*np1+l]=0e0;
  b[m*np1+l]=0e0;
  alphar[l]=a[l*np1+l];
  alphar[m]=a[m*np1+m];
  alphai[l]=0e0;
  alphai[m]=0e0;
  beta[l]=b[l*np1+l];
  beta[m]=b[m*np1+m];
  m-=2;
  goto l290;
 l280:
  er=c/(b11*b22);
  ei=SQRT(-d)/(b11*b22);
  a11r=a11-er*b11;
  a11i=ei*b11;
  a12r=a12-er*b12;
  a12i=ei*b12;
  a21r=a21;
  a21i=0e0;
  a22r=a22-er*b22;
  a22i=ei*b22;
  flip=(INT )((FABS(a11r)+FABS(a11i)+FABS(a12r)+FABS(a12i))	\
	      < (FABS(a21r)+FABS(a22r)+FABS(a22i)));
  if(flip)
    chouseh2(a22r,a22i,-a21r,-a21i,&cz,&szr,&szi);
  else
    chouseh2(a12r,a12i,-a11r,-a11i,&cz,&szr,&szi);
  //
  flip=(INT )(an < ((FABS(er)+FABS(ei))*b22));
  if(flip)
    chouseh2(cz*a11 + szr*a12,szi*a12,cz*a21 + szr*a22,szi*a22,&cq,&sqr,&sqi);
  else
    chouseh2(cz*b11 + szr*b12,szi*b12,szr*b22,szi*b22,&cq,&sqr,&sqi);
  //
  ssr=sqr*szr + sqi*szi;
  ssi=sqr*szi - sqi*szr;
  tr=cq*cz*a11 + cq*szr*a12 + sqr*cz*a21 + ssr*a22;
  ti=cq*szi*a12 - sqi*cz*a21 + ssi*a22;
  bdr=cq*cz*b11 + cq*szr*b12 + ssr*b22;
  bdi=cq*szi*b12 + ssi*b22;
  r=SQRT(bdr*bdr + bdi*bdi);
  beta[l]=bn*r;
  alphar[l]=an*(tr*bdr+ti*bdi)/r;
  alphai[l]=an*(tr*bdi-ti*bdr)/r;
  //
  tr=ssr*a11 -sqr*cz*a12 - cq*szr*a21 + cq*cz*a22;
  ti= -ssi*a11 - sqi*cz*a12 + cq*szi*a21;
  bdr= ssr*b11 - sqr*cz*b12 + cq*cz*b22;
  bdi= -ssi*b11 - sqi*cz*b12;
  r=SQRT(bdr*bdr + bdi*bdi);
  beta[m]=bn*r;
  alphar[m]=an*(tr*bdr+ti*bdi)/r;
  alphai[m]=an*(tr*bdi-ti*bdr)/r;
  m-=2;
 l290:
  if(m>0) goto l200;
  del_row_col(ndp1,np1,&a);
  del_row_col(ndp1,np1,&b);
  del_row_col(ndp1,np1,&x);
  a_io[0]=a;  b_io[0]=b;  x_io[0]=x;
  for(i=0;i<nd;++i){
    alphar[i]=alphar[i+1];
    alphai[i]=alphai[i+1];
    beta[i]=beta[i+1];
  }
  alphar_io[0]=realloc(alphar,nd*sizeof(REAL16));
  alphai_io[0]=realloc(alphai,nd*sizeof(REAL16));
  beta_io[0]=realloc(beta,nd*sizeof(REAL16));
  return;
}
/*********************************************************************************/
/*
 *\fn void qzvec(INT nd, INT n,
 *	   REAL16 **a_io, REAL16 **b_io,
 *	   REAL16 epsa,REAL16 epsb,
 *	   REAL16 **alphar_io,REAL16 **alphai_io,
 *	   REAL16 **beta_io,
 *	   REAL16 **x_io)
 *
 * \brief DOES NOT WORK DO NOT USE: (generalized eigenvectors).
 *
 * \note The algorithm is described in: MR0345399 (49 10135) 65F15
 *       Moler, C. B. ; Stewart, G. W. An algorithm for generalized
 *       matrix eigenvalue problems.  SIAM J. Numer. Anal. 10 (1973),
 *       241–256.
 *       ForTran routines are found in: Moler, C. B. ;
 *       Stewart, G. W.  An algorithm for the generalized matrix
 *       eigenvalue problem Ax=\lambda Bx.  UT Austin tech report
 *       N00014-67-A-0181-0023 issued jointly with Stanford techreport
 *       STAN-CS-232-71
 *
 *   Ludmil 20210426
 */
/********************************************************************************/
void qzvec(INT nd, INT n,					\
	   REAL16 **a_io, REAL16 **b_io,			\
	   REAL16 epsa,REAL16 epsb,				\
	   REAL16 **alphar_io,REAL16 **alphai_io,		\
	   REAL16 **beta_io,					\
	   REAL16 **x_io)
{
  return;

  INT flip;//, wantx;
  //
  // find eigenvectors of quasi-triangular matrices
  // uses b for intermediate storage
  //
  INT m,mr,mi,l,l1,i,j,k,np1=n+1,ndp1=nd+1;
  REAL16 alpham,betam,sl,sk,d,dr,di,skr,ski,	\
    tkk,tkl,tlk,tll,almr,almi,tr,ti,slr,sli,	\
    tkkr,tkki,tklr,tkli,tlkr,tlki,tllr,tlli;
  REAL16 s,r,zr,zi;
  //
  // find eigenvectors of quasi-triangular matrices
  // uses b for intermediate storage
  //
  /*****************/
  add_row_col(nd,n,a_io);
  add_row_col(nd,n,b_io);
  add_row_col(nd,n,x_io);
  // set it up values
  REAL16 *alphar=realloc(alphar_io[0], ndp1*sizeof(REAL16));
  REAL16 *alphai=realloc(alphai_io[0], ndp1*sizeof(REAL16));
  REAL16 *beta=realloc(beta_io[0], ndp1*sizeof(REAL16));
  for(i=(nd-1);i>=0;--i){
    alphar[i+1] = alphar[i];
    alphai[i+1] = alphai[i];
    beta[i+1]  = beta[i];
  }
  REAL16 *a=a_io[0], *b=b_io[0], *x=x_io[0];
  /***********************************************************************************/
  m = n;
 l300:
  if(alphai[m]!=0e0) goto l350;
  alpham = alphar[m];
  betam = beta[m];
  if (fabs(alpham)<epsa)
    alpham = 0e0;
  if (abs(betam)<epsb)
    betam =0e0;
  b[m*np1+m]=1e0;
  l=m-1;
  if(l==0) goto l340;
 l310:
  l1=l+1;
  sl=0e0;
  for(j=l1;j<=m;++j)
    sl+=(betam*a[l*np1+j]-alpham*b[l*np1+j])*b[j*np1+m];
  if(l==1) goto l320;
  if(a[l*np1+(l-1)]!=0e0) goto l330;
 l320:
  d=betam*a[l*np1+l]-alpham*b[l*np1+l];
  if(d==0e0)
    d=0.5e0*(epsa+epsb);
  b[l*np1+m]=-sl/d;
  l--;
  goto l340;
 l330:
  k=l-1;
  sk=0e0;
  for(j=l1;j<=m;++j)
    sk+=(betam*a[k*np1+j]-alpham*b[k*np1+j])*b[j*np1+m];
  tkk=betam*a[k*np1+k]-alpham*b[k*np1+k];
  tkl=betam*a[k*np1+l]-alpham*b[k*np1+l];
  tlk=betam*a[l*np1+k];
  tll=betam*a[l*np1+l]-alpham*b[l*np1+l];
  d=tkk*tll-tkl*tlk;
  if(d==0e0)
    d=0.5*(epsa+epsb);
  b[l*np1+m]=(tlk*sk-tkk*sl)/d;
  flip=(INT )(fabs(tkk)<fabs(tlk));
  if(flip)
    b[k*np1+m]=-(sl + tll*b[l*np1+m])/tlk;
  else
    b[k*np1+m]=-(sk + tkl*b[l*np1+m])/tkk;
  l-=2;
 l340:
  if(l>0) goto l310;
  m--;
  goto l390;
 l350:
  // complex vector
  almr=alphar[m-1];
  almi=alphai[m-1];
  betam=beta[m-1];
  mr=m-1;
  mi=m;
  //
  b[(m-1)*np1+mr]=almi*b[m*np1+m]/(betam*a[m*np1+(m-1)]);
  b[(m-1)*np1+mi]=(betam*a[m*np1+m]-almr*b[m*np1+m])/(betam*a[m*np1+(m-1)]);
  b[m*np1+mr]=0e0;
  b[m*np1+mi]=-1e0;
  l=m-2;
  if(l==0) goto l385;
 l360:
  l1=l+1;
  slr=0e0;
  sli=0e0;
  for(j=l1;j<=m;j++){
    tr=betam*a[l*np1+j]-almr*b[l*np1+j];
    ti=-almi*b[l*np1+j];
    slr+=tr*b[j*np1+mr] - ti*b[j*np1+mi];
    sli+=tr*b[j*np1+mi] + ti*b[j*np1+mr];
  }
  if(l==1) goto l370;
  if(a[l*np1+(l-1)]!=0e0) goto l375;
 l370:
  dr=betam*a[l*np1+l]-almr*b[l*np1+l];
  di=-almi*b[l*np1+l];
  xdivy(-slr,-sli,dr,di,&b[l*np1+mr],&b[l*np1+mi]);
  l--;
  goto l385;
 l375:
  k=l-1;
  skr=0e0;
  ski=0e0;
  for(j=l1;j<=m;++j){
    tr=betam*a[k*np1+j]-almr*b[k*np1+j];
    ti=-almi*b[k*np1+j];
    skr+=tr*b[j*np1+mr] - ti*b[j*np1+mi];
    ski+=tr*b[j*np1+mi] + ti*b[j*np1+mr];
  }
  tkkr=betam*a[k*np1+k]-almr*b[k*np1+k];
  tkki=-almi*b[k*np1+k];
  //
  tklr=betam*a[k*np1+l]-almr*b[k*np1+l];
  tkli=-almi*b[k*np1+l];
  //
  tlkr=betam*a[l*np1+k];
  tlki=0e0;
  //
  tllr=betam*a[l*np1+l]-almr*b[l*np1+l];
  tlli=-almi*b[l*np1+l];
  //
  dr=tkkr*tllr - tkki*tlli - tklr*tlkr;
  di=tkkr*tlli + tkki*tllr - tkli*tlkr;
  if(dr==0e0 && di==0e0)
    dr=0.5e0*(epsa+epsb);
  xdivy(tlkr*skr-tkkr*slr+tkki*sli,		\
	tlkr*ski-tkkr*sli+tkki*slr,		\
	dr,di,&b[l*np1+mr],&b[l*np1+mi]);
  flip = (INT )((fabs(tkkr)+fabs(tkki))>=fabs(tlkr));
  if(flip){
    xdivy(-skr-tklr*b[l*np1+mr]+tkli*b[l*np1+mi],		\
	  -ski-tklr*b[l*np1+mi]-tkli*b[l*np1+mr],		\
	  tkkr, tkki,&b[k*np1+mr],&b[k*np1+mi]);
  }else{
    xdivy(-slr-tllr*b[l*np1+mr]+tlli*b[l*np1+mi],		\
	  -sli-tllr*b[l*np1+mi]-tlli*b[l*np1+mr],		\
	  tlkr, tlki,&b[k*np1+mr],&b[k*np1+mi]);
  }
  l-=2;
 l385:
  if(l>0) goto l360;
  m-=2;
 l390:
  if(m>0) goto l300;
  m=n;
 l400:
  for(i=1;i<=n;++i){
    s=0e0;
    for(j=1;j<=m;++j){
      s += x[i*np1+j]*b[j*np1+m];
    }
    x[i*np1+m]=s;
  }
  m--;
  if(m>0) goto l400;
  m=n;
 l430:
  s=0e0;
  if(alphai[m]!=0e0) goto l450;
  for(i=1;i<=n;++i){
    r=fabs(x[i*np1+m]);
    if(r<s) continue;
    s=r;
    d=x[i*np1+m];
  }
  for(i=1;i<=n;++i)
    x[i*np1+m]/=d;
  m--;
  goto l490;
 l450:
  for(i=1;i<=n;++i){
    r=x[i*np1+(m-1)]*x[i*np1+(m-1)] + x[i*np1+m]*x[i*np1+m];
    if(r<s) continue;
    s=r;
    dr=x[i*np1+(m-1)];
    di=x[i*np1+m];
  }
  for(i=1;i<=n;++i){
    xdivy(x[i*np1+(m-1)],x[i*np1+m],dr,di,&zr,&zi);
    x[i*np1+(m-1)]=zr;
    x[i*np1+m]=zi;
  }
  m-=2;
 l490:
  if(m>0) goto l430;
  del_row_col(ndp1,np1,&a);
  del_row_col(ndp1,np1,&b);
  del_row_col(ndp1,np1,&x);
  a_io[0]=a;  b_io[0]=b;  x_io[0]=x;
  for(i=0;i<nd;++i){
    alphar[i]=alphar[i+1];
    alphai[i]=alphai[i+1];
    beta[i]=beta[i+1];
  }
  alphar_io[0]=realloc(alphar,nd*sizeof(REAL16));
  alphai_io[0]=realloc(alphai,nd*sizeof(REAL16));
  beta_io[0]=realloc(beta,nd*sizeof(REAL16));
  return;
}
