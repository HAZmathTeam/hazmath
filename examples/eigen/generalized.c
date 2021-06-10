/*! \file examples/Solver/generalized.c
 *
 *  Created by Xiaozhe Hu on 01/01/19.
 *  Copyright 2019_HAZMATH__. All rights reserved.
 *
 * \brief This program read a matrices A and B and
 *         solves the generalized eigenvalue problem A*x=lambda*B*x using qz algorithm. 
 *
 * \note Ludmil 20210426
 *
 */

/************* HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
/*******************************************************************************/
/***********************************************************************************/
/*  \fn static void xprint_full_mat(const INT strti,const INT strtj,
 *                           const  INT n, const INT m, 
 *                           REAL16 *A,const char *varname)
 *
 *
 *  \note: prints a matrix A starting at (starti, startj) with (n)
 *         rows and (m) columns in matlab format e.g. 2 x 2 identity
 *         is printed as I2=[1. 0.;0. 1.]; if the varname= "I2"

 *  \note: most often strti=strtj=0 but if they are 1, then the matrix
 *         starts is 1:n; 1:m with the last two indices included.
 *
 */
static void xprint_full_mat(const INT strti,const INT strtj,			\
		     const  INT n, const INT m, REAL16 *A,const char *varname)
{
  INT nprt=1025,mprt=1025;
  if( (n<1) || (m<1) ) return;
  INT i,j,n1=n-1;
  if(n<=nprt) nprt=n;
  if(m<=mprt) mprt=m;
  if(varname==NULL){
    fprintf(stdout,"\nA=[");
  }else{
    fprintf(stdout,"\n%s=[",varname);
  }
  for (i = strti; i<(nprt+strti);++i){
    for(j = strtj;j<(mprt+strtj);++j){
      fprintf(stdout,"%25.20Le ", A[(m+strtj)*i+j]);
    }
    if(i!=n1+strti){
      fprintf(stdout,";");
    }else{
      fprintf(stdout,"];\n");
    }
  }
  return;
}
/****************************************************************************/
INT main()
{
  INT dummy,nr,nc,nd,n,i,j,k;  
  FILE *fp;
  REAL16 epsa=-1e20,epsb=-1e20,tol=1e-20; //tol is short for eps
  fp=fopen("a_matrix1.dat","r");
  fscanf(fp,"%d %d",&nr,&nc);
  nd=nr;
  n=nc;
  // nd=nc (number of cols; this is the leading dimension of a in C)
  REAL16 *a=(REAL16 *)calloc(nd*nd,sizeof(REAL16));
  REAL16 *b=(REAL16 *)calloc(nd*nd,sizeof(REAL16));
  REAL16 *x=(REAL16 *)calloc(nd*nd,sizeof(REAL16));
  REAL16 *alphar=(REAL16 *)calloc(nd,sizeof(REAL16));
  REAL16 *alphai=(REAL16 *)calloc(nd,sizeof(REAL16));
  REAL16 *beta=(REAL16 *)calloc(nd,sizeof(REAL16));
  INT *iter=(INT *)calloc(nd,sizeof(INT));
  // a
  for(i=0;i<nr;++i)
    for(j=0;j<nc;++j)
      fscanf(fp,"%Lg",(a+i*n+j));
  fscanf(fp,"%d",&dummy);
  // b
  for(i=0;i<nr;++i)
    for(j=0;j<nc;++j)
      fscanf(fp,"%Lg",(b + i*nc+j));
  // c
  INT wantx=1;
  xprint_full_mat(0,0,nd, nd, a,"Ao");
  xprint_full_mat(0,0,nd, nd, b,"Bo");
  qzhes(nd,nd,a,b,wantx,x);
  qzit(nd,nd,&a,&b,tol,&epsa,&epsb,&iter,wantx,&x);
  //  xprint_full_mat(0,0,nd, nd, a,"An");
  //  xprint_full_mat(0,0,nd, nd, b,"Bn");
  //  xprint_full_mat(0,0,nd, nd, x,"Xn");
  qzval(nd,nd,
	&a,&b,
	epsa,epsb,
	&alphar,&alphai,&beta,
	wantx,&x);
  //
  xprint_full_mat(0,0,nd, nd, a,"An");
  xprint_full_mat(0,0,nd, nd, b,"Bn");
  xprint_full_mat(0,0,nd, nd, x,"Xn");
  xprint_full_mat(0,0,nd,1, alphar,"alphar");
  xprint_full_mat(0,0,nd,1, alphai,"alphai");
  xprint_full_mat(0,0,nd,1, beta,  "beta");
  fprintf(stdout,"\n%s\n%s\n%s\n%s\n%s\n",					\
	  "er=alphar./beta;",						\
	  "ei=alphai./beta;",						\
	  "err_r=sort(real(eig(Ao,Bo)))-sort(er);",		\
	  "err_i=sort(imag(eig(Ao,Bo)))-sort(ei);",			\
	  "err=sqrt(err_r'*err_r+err_i'*err_i)");
  fprintf(stdout,"\n%s\n%s\n%s\n%s\n%s\n",					\
	  "er=alphar./beta;",						\
	  "ei=alphai./beta;",						\
	  "err_r=sort(real(eig(An,Bn)))-sort(er);",		\
	  "err_i=sort(imag(eig(An,Bn)))-sort(ei);",			\
	  "err=sqrt(err_r'*err_r+err_i'*err_i)");
  /* qzvec(nd,nd,&a,&b,				\ */
  /* 	epsa,epsb,				\ */
  /* 	&alphar,&alphai,&beta,			\ */
  /* 	&x); */
  /* xprint_full_mat(0,0,nd, nd, a,"Ae"); */
  /* xprint_full_mat(0,0,nd, nd, b,"Be"); */
  /* xprint_full_mat(0,0,nd, nd, x,"Xe"); */
  free(a);
  free(b);
  free(x);
  //
  free(alphar);
  free(alphai);
  free(beta);
  free(iter);
  return 0;
}
