/*! \file examples/approximation/aaa.c
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
#include "aaa_example.h"
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
INT main()
{
  REAL16 s[4]={0.5e0,-0.5e0,1e0,2e0}; // s1, s2, alpha,beta
  REAL xmin_in=0e0,xmax_in=1e0;
  ///////////////////////////////////////////////////////////////////////////
  //mbig is the total number of points; mmax_in is the max number of
  //nodes taking part in the interpolation, mmax is much smaller than
  //mbig;
  //
  // m is the number of nodes in the final interpolation after
  // tolerance tolaaa is achieved or mmax is reached.
  INT mbig=(1<<10)+1,mmax_in=2<<5+1,m=-22;
  REAL16 tolaaa=powl(2e0,-52e0);
  // parameters for the function we are approximating.
  //
  // For example: s[2]*(x^s[0])+s[3]*(x**s[1])
  s[0]=0.5;  s[1]=-0.5;  s[2]=2e0;  s[3]=1e0;
  //
  //
  REAL **wzpc=malloc(5*sizeof(REAL *));
  // after a call to get_wzpc, wzpc[k] are pointers to the following
  // arrays (each of m+1 elements but at most m are used:
  //
  // wz[0]-> nodes (z[]);
  // wz[1]-> weights (w[]);
  // wz[2]-> function values (f[]);
  // wz[3]-> poles (pol[]);
  // wz[4]-> residues (res[]) (also as last entry contains the free term c[m];
  //
  // the rational approximation is:
  // r(z)=res[m-1] + \sum_{i=0}^{m-2} res[i]/(z-pol[i]); 
  //
  REAL rmax=get_wzpc(f_to_approx_l,			\
		     (void *)s,				\
		     wzpc,				\
		     &mbig,&mmax_in,			\
		     &m,xmin_in,xmax_in,tolaaa);
  //rmax is the max error on the rest of z
  INT mem,mem16,i,j,k,m1=m+1,m01=m-1,mm=m*m,mm1=m1*m1;  
  fprintf(stdout,"\n\nm=%d; max_err=%.12e\n",m,rmax);
  
  for(i=0;i<m-1;i++)fprintf(stdout,"\nz(%d)=%.16e;",i+1,*(wzpc[0]+i));
  fprintf(stdout,"\n");
  for(i=0;i<m-1;i++) fprintf(stdout,"\nw(%d)=%.16e;",i+1,*(wzpc[1]+i));
  fprintf(stdout,"\n");
  for(i=0;i<m-1;i++) fprintf(stdout,"\nf(%d)=%.16e;",i+1,*(wzpc[2]+i));
  fprintf(stdout,"\n");
  for(i=0;i<m-1;i++)fprintf(stdout,"\npol(%d)=%.16e;",i+1,*(wzpc[3]+i));
  fprintf(stdout,"\n");
  for(i=0;i<m-1;i++)fprintf(stdout,"\nres(%d)=%.16e;",i+1,*(wzpc[4]+i));
  fprintf(stdout,"\n");
  free(wzpc[0]);// that is enough. the rest 1-4 are just shifts of wzpc[0] with m+1;
  free(wzpc);
  return 0;  
}
