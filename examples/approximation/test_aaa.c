/*! \file examples/approximation/test_aaa.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 2019/01/09.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief Computes approximations with rational functions using AAA
 *        algorithm
 *
 * \note This example shows how to use AAA method to compute rational
 *       approximation of a function f_to_approx.
 *
 *\note The function would return long double and its arguments are
 *      long double
 */

#include "hazmath.h"
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
#include "ra_examples.h"
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

INT main(int argc,char *argv[])
{
  REAL xmin_in=0e0,xmax_in=1e0;
  // parameters for the function we are approximating.
  // For example: s[2]*(x^s[0])+s[3]*(x**s[1])
  REAL16 s[4]={-0.5e0,0.5e0,1e0,0e0}; // s1, s2, alpha,beta
  INT print_level=1;

  fprintf(stderr,"\nUSAGE: %s<<EOF_FRAC\n s,t,alpha,beta,xmin,xmax\nEOF_FRAC\n",argv[0]);
  fprintf(stderr,"\nEXAMPLE:\n%s<<EOF_FRAC >frac_aaa.m\n %.2Lf %.2Lf %.2Lf %.2Lf %.2f %.2f\nEOF_FRAC\n", \
	  argv[0],s[0],s[1],s[2],s[3],xmin_in,xmax_in);
  INT k=fscanf(stdin,"%Lg %Lg %Lg %Lg %lg %lg",&s[0],&s[1],&s[2],&s[3],&xmin_in,&xmax_in);
  INT i;
  INT mbig=(1<<14)+1;
  INT mmax_in=(INT )(mbig/2),m=-22;
  REAL16 tolaaa=powl(2e0,-46e0);
  //
  INT mem,mem16,j,m1=m+1,m01=m-1,mm=m*m,mm1=m1*m1;
  // s[0]=s,s[1]=t,s[2]=alpha,s[3]=beta;
  /////////////////////////////////////////////
  //for debug: s[0]=0.5;  s[1]=-0.5;  s[2]=1e0;  s[3]=2e0; xmin_in=0e0;xmax_in=1e0;
  // s[0]=-0.50; s[1]=0.50;   s[2]=6.1031431e-03;   s[3]=1e0;   xmin_in=0e0;   xmax_in=1e0;
  /////////////////////////////////////////////////////////////////////////////////
  REAL **cpzwf=malloc(5*sizeof(REAL *));
  // after a call to get_cpzwf, cpzwf[k] are pointers to the following
  // arrays (each of m+1 elements but at most m are used:
  //
  // cpzwf[0]-> residues (res[]); (also as last entry contains the free
  //                               term c[m])
  // cpzwf[1]-> poles (pol[]);
  // cpzwf[2]-> nodes (z[]);
  // cpzwf[3]-> weights (w[]);
  // cpzwf[4]-> function values (f[])
  //
  // the rational approximation is:
  // r(z)=res[m-1] + \sum_{i=0}^{m-2} res[i]/(z-pol[i]);
  //
  REAL rmax=get_cpzwf(f_to_approx_l,				\
		      (void *)s,				\
		      cpzwf,					\
		      &mbig,&mmax_in,				\
		      &m,xmin_in,xmax_in,tolaaa,print_level);
  //rmax is the max error on the rest of z
  fprintf(stdout,"\n%%%% function [res,pol,z,w,f,er0,er2]=frac_aaa()\n");
  fprintf(stdout,"\n%% AUTO GENERATED\n");
  fprintf(stdout,"\n%%%%EXAMPLE(fractional):\n");
  fprintf(stdout,"\nx=linspace(%.8e,%.8e,%d);\nif(size(x,1)==1),x=x';end\n",xmin_in,xmax_in,1025);
  fprintf(stdout,"\nf_in = 1./(%.2Lf * x.^(%.1Lf) + %.2Lf * x.^(%.1Lf) );\n",s[2],s[0],s[3],s[1]);
  fprintf(stdout,"\n%%===============================================%%\n");
  fprintf(stdout,"\nbounds=[%.16e,%.16e];\n",xmin_in,xmax_in);
  fprintf(stdout,"\n\nm=%d; max_err=%.12e;\n",m,rmax);
  fprintf(stdout,"\nres=zeros(%d,1);pol=zeros(%d,1);\nz=zeros(%d,1);w=zeros(%d,1);f=zeros(%d,1);\n", \
	  m,m-1,m,m,m);
  fprintf(stdout,"\n%%===============================================%%\n");

  for(i=0;i<m;i++)fprintf(stdout,"\nres(%d)=%.16e;",i+1,*(cpzwf[0]+i));
  fprintf(stdout,"\n");
  for(i=0;i<m-1;i++) fprintf(stdout,"\npol(%d)=%.16e;",i+1,*(cpzwf[1]+i));
  fprintf(stdout,"\n");
  for(i=0;i<m;i++) fprintf(stdout,"\nz(%d)=%.16e;",i+1,*(cpzwf[2]+i));
  fprintf(stdout,"\n");
  for(i=0;i<m;i++)fprintf(stdout,"\nw(%d)=%.16e;",i+1,*(cpzwf[3]+i));
  fprintf(stdout,"\n");
  for(i=0;i<m;i++)fprintf(stdout,"\nf(%d)=%.16e;",i+1,*(cpzwf[4]+i));
  fprintf(stdout,"\n");  
  fprintf(stdout,"f_ra=1./(kron(x,ones(size(pol\')))-kron(ones(size(x)),pol\'));");
  fprintf(stdout,"f_ra=res(1)+f_ra*(res(2:%d));",m);
  fprintf(stdout,"\n%%%%%%%%\t\tfz = 1./(%.2Lf * z.^(%.1Lf) + %.2Lf * z.^(%.1Lf));\n",s[2],s[0],s[3],s[1]);
  fprintf(stdout,"\ner0=norm(f_in-f_ra);");
  fprintf(stdout,"\ner2=norm(f_in(2:length(x))-f_ra(2:length(x)));");
  fprintf(stdout,"\n%%%%return;end\n");
  fprintf(stdout,"\n%%===============================================%%\n\n");
  free(cpzwf[0]);// that is enough;
  free(cpzwf);
  return 0;
}
