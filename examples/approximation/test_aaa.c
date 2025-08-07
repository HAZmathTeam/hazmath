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
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*SOME SETTINGS BELOW*/
///////////////////////////////////////////////////////////////////////////////
#ifndef MAX_NUMBER_POLES
#define MAX_NUMBER_POLES 25
#endif
#ifndef TOL_POWER
#define TOL_POWER -46e0
#endif
#ifndef PRINT_LEVEL
#define PRINT_LEVEL 0
#endif
///////////////////////////////////////////////////////////////////////////////
INT main(int argc,char *argv[])
{
  REAL xmin_in=0e0,xmax_in=1e0;
  INT print_level=PRINT_LEVEL;
  // parameters for the function we are approximating.
  REAL16 s=0.6e0,t=-0.2e0,alpha=1e0,beta=1e0; // s, t, alpha,beta
  fprintf(stderr,"\nUSAGE: %s<<EOF_FRAC\n s t alpha beta xmin xmax\nEOF_FRAC\n",argv[0]);
  fprintf(stderr,"\nEXAMPLE:\n%s<<EOF_FRAC >frac_aaa.m\n %.2Lf %.2Lf %.2Lf %.2Lf %.2f %.2f\nEOF_FRAC\n", \
	  argv[0],s,t,alpha,beta,xmin_in,xmax_in);
  INT k=fscanf(stdin,"%Lg %Lg %Lg %Lg %lg %lg",&s,&t,&alpha,&beta,&xmin_in,&xmax_in);
  INT i;
  INT numval=(1<<14)+1;
  INT mmax_in=(INT )(numval/2),m;
  mmax_in=MAX_NUMBER_POLES;//
  mmax_in=40;
  REAL16 tolaaa=powl(2e0,TOL_POWER);
  //
  INT mem,mem16,j,m1=m+1,m01=m-1,mm=m*m,mm1=m1*m1;
  /////////////////////////////////////////////////////////////////////////////////
  REAL **rpzwf=malloc(7*sizeof(REAL *));
  // after a call to get_rpzwf, rpzwf[k] are pointers to the following
  // arrays (each of m+1 elements but at most m are used):
  //
  // rpzwf[0]-> real part of residues (resr[]); (also as last entry contains the free
  //                               term resr[m])
  // rpzwf[1]-> imaginary part of residues (resi[]); (also as last entry contains the free
  //                               term resi[m])
  // rpzwf[2]-> real(poles) (polr[]);
  // rpzwf[3]-> imag(poles) (poli[]);
  // rpzwf[4]-> nodes (z[]);
  // rpzwf[5]-> weights (w[]);
  // rpzwf[6]-> function values (f[])
  //
  // the rational approximation is:
  // r(z)=res[m-1] + \sum_{i=0}^{m-2} res[i]/(z-pol[i]);
  //
  REAL16 **zf=set_f_values(f_to_approx_l,		\
			   s,t,alpha,beta,		\
			   &numval,xmin_in,xmax_in,	\
			   print_level);
  //  for(i=0;i<numval;i++)fprintf(stdout,"\nz(%d)=%.16Le;f(%d)=%.16Le;",i+1,*(zf[0]+i),i+1,*(zf[1]+i));
  //  fprintf(stdout,"\n");
  REAL rmax=get_rpzwf(numval,zf[0],zf[1],	\
		      rpzwf,			\
		      &mmax_in, &m,		\
		      tolaaa,print_level);
  free(zf[0]);
  free(zf[1]);
  free(zf);
  //rmax is the max error on the rest of z
  fprintf(stdout,"\n%%%% function [res,pol,z,w,f,er0,er2]=frac_aaa()\n");
  fprintf(stdout,"\n%% AUTO GENERATED\n");
  fprintf(stdout,"\n%%%%EXAMPLE(fractional):\n");
  fprintf(stdout,"\nxmin=%.8e;xmax=%.8e;npts=%lld;h=(xmax-xmin)/(npts-1);\n",xmin_in,xmax_in,(long long )4097);
  fprintf(stdout,"\nx=linspace(xmin,xmax,npts);\nif(size(x,1)==1),x=x';end\n");
//  fprintf(stdout,"\nf_in = 1./(%.2Lf * x.^(%.1Lf) + %.2Lf * x.^(%.1Lf) );\n",alpha,s,beta,t);
  fprintf(stdout,"\nf_in = (%.2Lf * x.^(%.1Lf) + %.2Lf * x.^(%.1Lf) );\n",alpha,s,beta,t);
  fprintf(stdout,"\n%%===============================================%%\n");
  fprintf(stdout,"\nbounds=[%.16e,%.16e];\n",xmin_in,xmax_in);
  fprintf(stdout,"\n\nm=%lld; max_err=%.12e;\n",(long long )m,rmax);
  fprintf(stdout,"\nres=zeros(%lld,1);pol=zeros(%lld,1);\nz=zeros(%lld,1);w=zeros(%lld,1);f=zeros(%lld,1);", \
	  (long long )m,(long long )m-1,(long long )m,(long long )m,(long long )m);
  fprintf(stdout,"\n%%===============================================%%\n");
  for(i=0;i<m;i++)fprintf(stdout,"\nres(%lld)=%.16e + 1i*(%.16e);",(long long )(i+1),*(rpzwf[0]+i),*(rpzwf[1]+i));
  fprintf(stdout,"\n");
  for(i=0;i<m-1;i++) fprintf(stdout,"\npol(%lld)=%.16e + 1i*(%.16e);",(long long )(i+1),*(rpzwf[2]+i),*(rpzwf[3]+i));
  fprintf(stdout,"\n");
  for(i=0;i<m;i++) fprintf(stdout,"\nz(%lld)=%.16e;",(long long )(i+1),*(rpzwf[4]+i));
  fprintf(stdout,"\n");
  for(i=0;i<m;i++)fprintf(stdout,"\nw(%lld)=%.16e;",(long long )(i+1),*(rpzwf[5]+i));
  fprintf(stdout,"\n");
  for(i=0;i<m;i++)fprintf(stdout,"\nf(%lld)=%.16e;",(long long )(i+1),*(rpzwf[6]+i));
  fprintf(stdout,"\n");  
  fprintf(stdout,"f_ra=1./(kron(x,ones(size(pol\')))-kron(ones(size(x)),conj(pol\')));");
  fprintf(stdout,"f_ra=res(1)+f_ra*(res(2:%lld));",(long long )m);
//  fprintf(stdout,"\n%%%%%%%%\t\tf_of_z = 1./(%.2Lf * z.^(%.1Lf) + %.2Lf * z.^(%.1Lf));\n",alpha,s,beta,t);
  fprintf(stdout,"\n%%%%%%%%\t\tf_of_z = (%.2Lf * z.^(%.1Lf) + %.2Lf * z.^(%.1Lf));\n",alpha,s,beta,t);
  fprintf(stdout,"\nerL2=norm(f_in-f_ra,2)/norm(f_in,2);");
  fprintf(stdout,"\nerinf=norm(f_in-f_ra,inf)/norm(f_in,inf);");
  fprintf(stdout,"\n%%%%return;end\n");
  fprintf(stdout,"\n%%===============================================%%\n\n");
  free(rpzwf[0]);// that is enough;
  free(rpzwf);
  return 0;
}
