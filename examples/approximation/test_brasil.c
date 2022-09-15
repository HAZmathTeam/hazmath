/*! \file examples/approximation/test_brasil.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 2019/01/09.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief Computes approximations with rational functions using BRASIL
 *        algorithm
 *
 * \note See: C. Hofreither, "An algorithm for best rational
 *       approximation based on barycentric rational interpolation."
 *       Numerical Algorithms, 2021
 *       https://doi.org/10.1007/s11075-020-01042-0
 *
 * \author Created by Clemens Hofreither on 20210720.
 *
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
  REAL16 s[4]={-0.5e0,0.5e0,1e0,0.25e0}; // s1, s2, alpha,beta
  INT print_level=1;
  INT iter_brasil;
  fprintf(stderr,"\nUSAGE: %s<<EOF_FRAC\n s,t,alpha,beta,xmin,xmax\nEOF_FRAC\n",argv[0]);
  fprintf(stderr,"\nEXAMPLE:\n%s<<EOF_FRAC >frac_brasil.m\n %.2Lf %.2Lf %.2Lf %.2Lf %.2f %.2f\nEOF_FRAC\n", \
	  argv[0],s[0],s[1],s[2],s[3],xmin_in,xmax_in);
  INT k=fscanf(stdin,"%Lg %Lg %Lg %Lg %lg %lg",&s[0],&s[1],&s[2],&s[3],&xmin_in,&xmax_in);
  /////// BEST APPROXIMATION BY BRASIL ALGORITHM ///////
  REAL *rpzwf_brasil[7]; // for real and complex;
  const INT degree = 15;
  INT init_steps=128;
  INT  maxiter=1024;
  REAL   step_factor=1./16.;
  REAL   max_step_size=1./16.;
  REAL   tol_brasil=pow(2.,-14.);// < 1/16000
  if(fabs(xmin_in)<tol_brasil) 
    xmin_in=tol_brasil*((REAL )(1<<4));
  REAL rmax=get_rpzwf_brasil(old_f_to_approx_l, (void*)s, rpzwf_brasil,
			     xmin_in, xmax_in, degree,     // the remaining options can usually be kept at these defaults
			     init_steps, maxiter, step_factor, max_step_size, tol_brasil, &iter_brasil,print_level);  
  INT i;
  //rmax is the max error on the rest of z
  fprintf(stdout,"\n%%%%function [res,pol,z,w,f,er0,er2]=frac_brasil()\n");
  fprintf(stdout,"\n%% AUTO GENERATED\n");
  fprintf(stdout,"\n%%%%EXAMPLE(fractional):\n");
  fprintf(stdout,"\nx=linspace(%.8e,%.8e,%lld);\nif(size(x,1)==1),x=x';end\n",xmin_in,xmax_in,(long long )1025);
  //  fprintf(stdout,"\n%%f_in = @(x) ( %.2Lf * x.^(%.1Lf) + %.2Lf * x.^(%.1Lf) );\n",s[2],s[0],s[3],s[1]);
  fprintf(stdout,"\nf_in = %.2Lf * x.^(%.1Lf) + %.2Lf * x.^(%.1Lf) ;\n",s[2],s[0],s[3],s[1]);
  fprintf(stdout,"\n%%===============================================%%\n");
  fprintf(stdout,"\nbounds=[%.16e,%.16e];\n",xmin_in,xmax_in);
  fprintf(stdout,"\n\ndegree=%lld; num_iter=%lld; max_err=%.12e;\n",(long long )degree,(long long )iter_brasil,rmax);
  // pol is one dimension less when printing. 
  fprintf(stdout,"\nres=zeros(%lld,1);pol=zeros(%lld,1);\nz=zeros(%lld,1);w=zeros(%lld,1);f=zeros(%lld,1);\n", \
	  (long long )degree+1,(long long )degree,(long long )degree+1,(long long )degree+1,(long long )degree+1);
  fprintf(stdout,"\n%%===============================================%%\n");
  for(i=0;i<degree+1;i++) fprintf(stdout,"\nres(%lld)=%.16e;",(long long )i+1,*(rpzwf_brasil[0]+i));
  fprintf(stdout,"\n");
  for(i=0;i<degree;i++) fprintf(stdout,"\npol(%lld)=%.16e;",(long long )i+1,*(rpzwf_brasil[1]+i));
  fprintf(stdout,"\n");
  for(i=0;i<degree+1;i++) fprintf(stdout,"\nz(%lld)=%.16e;",(long long )i+1,*(rpzwf_brasil[2]+i));
  fprintf(stdout,"\n");
  for(i=0;i<degree+1;i++) fprintf(stdout,"\nw(%lld)=%.16e;",(long long )i+1,*(rpzwf_brasil[3]+i));
  fprintf(stdout,"\n");
  for(i=0;i<degree+1;i++) fprintf(stdout,"\nf(%lld)=%.16e;",(long long )i+1,*(rpzwf_brasil[4]+i));
  fprintf(stdout,"\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"f_ra=1./(kron(x,ones(size(pol\')))-kron(ones(size(x)),pol\'));");
  fprintf(stdout,"f_ra=res(1)+f_ra*(res(2:%lld));",(long long )degree+1);
  fprintf(stdout,"\n%%%%%%%%\t\tfz = %.2Lf * z.^(%.1Lf) + %.2Lf * z.^(%.1Lf);\n",s[2],s[0],s[3],s[1]);
  fprintf(stdout,"\ner0=norm(f_in-f_ra);");
  fprintf(stdout,"\ner2=norm(f_in(2:length(x))-f_ra(2:length(x)));");
  fprintf(stdout,"\n%%%%return;end\n");
  fprintf(stdout,"\n%%===============================================%%\n\n");
  free(rpzwf_brasil[0]);        // all 5 arrays are allocated as one big block
  return 0;
}
/* #if 0       // debug code */
/* static void bary_print(FILE *stream, const barycentric_t *bary) */
/* { */
/*   INT i; */
/*   fprintf(stream, "Barycentric function of degree %d:\n Nodes: [", bary->nn - 1); */
/*   for (i = 0; i < bary->nn; ++i) */
/*     fprintf(stream, "%f ", bary->z[i]); */
/*   fprintf(stream, "]\n Weights: ["); */
/*   for (i = 0; i < bary->nn; ++i) */
/*     fprintf(stream, "%f ", bary->w[i]); */
/*   fprintf(stream, "]\n Values: ["); */
/*   for (i = 0; i < bary->nn; ++i) */
/*     fprintf(stream, "%f ", bary->f[i]); */
/*   fprintf(stream, "]\n"); */
/* } */
/* static void test_bary() */
/* { */
/*   // set up rational function r(x) = (x + 2) / (x + 1) */
/*   barycentric_t bary; */
/*   bary_init(&bary, 2); */
/*   REAL nodes[] = { 0.0, 0.5, 1.0 }; */
/*   REAL values[] = { 2.0, 2.5 / 1.5 , 3.0/2.0 }; */
/*   REAL *L = calloc((bary->nn + 1) * (bary->nn + 1), sizeof(REAL)); */
/*   bary_interpolate(&bary, nodes, values,L); */
/*   bary_print(stdout, &bary); */
/*   for (int i = 0; i <= 10; ++i) { */
/*     REAL x = 0.1 * i; */
/*     printf("x = %f, r(x) = %f, f(x) = %f\n", x, bary_eval(&bary, x), (x + 2) / (x + 1)); */
/*   } */
/*   bary_free(&bary); */
/* } */

/* #endif      // end test code */
