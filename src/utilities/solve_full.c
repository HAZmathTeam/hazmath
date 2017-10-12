/*! \file src/utilities/solve_full.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 5/13/15.
 *  Copyright 2016__HAZMATH__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/25/2016
 *  \note: modified on 10/12/2017 (ltz1)
 *
 * \note Solves A x = b with scalled partial pivoting for general full
 * n,n A. if dopivot is true it does the decomposition and if dopivot
 * is false assumes A is already LU decomposed and only solves.
 * 
 */

#include "hazmath.h"

INT solve_pivot(INT dopivot, INT n, REAL *A, REAL *b, INT *p,REAL *piv)
{
/* solves A x = b using Gaussian elimination with scalled partial
   pivoting.*/
  INT nm1,i1,k1,pin,kswp,kp,i,j,k;
  REAL r,t,absaij;
  REAL *x=piv;
  if(dopivot) {
    for (i=0;i<n;i++){
      p[i]=i;
      piv[i] = fabs(A[i*n+0]);
      for (j=1;j<n;j++){
      absaij=fabs(A[i*n+j]);
      if(absaij > piv[i])
	piv[i] = absaij;
      }
      piv[i]=1./piv[i]; //here we need error stop if too small
    }
    nm1 = n-1;
    for (k = 0;k<nm1;k++){
      r = fabs(A[p[k]*n+k])*piv[p[k]];
      kp = k;
      for (i=k;i<n;i++){
	t = fabs(A[p[i]*n+k])*piv[p[i]]; 
	if (t > r){r = t; kp = i;}
      }
      kswp = p[kp]; p[kp] = p[k]; p[k] = kswp;
      k1 = k+1; 
      for (i = k1;i<n;i++){
	pin=p[i]*n;
	A[pin+k] = A[pin+k]/A[p[k]*n+k];
	for(j = k1;j<n;j++){
	  A[pin+j] = A[pin+j]-A[pin+k]*A[p[k]*n+j];
	}
      }
    }
    /*end of decomposition; now solver part*/
  }
  x[0] = b[p[0]];
  for (i = 1;i<n;i++){
    x[i] = b[p[i]];
    for(j = 0;j<i;j++){
      x[i] = x[i]-A[p[i]*n+j]*x[j];
    }
  }
  nm1=n-1;
  x[nm1] = x[nm1]/A[p[nm1]*n+nm1];
  for (i = nm1;i>0;i--){
    i1=i-1;
    pin=p[i1]*n;
    for (j = i;j<n;j++){
      x[i1] = x[i1]-A[pin+j]*x[j];
    }
    x[i1] = x[i1]/A[pin+i1];
  }
  //  if(y) free(y);
  for(k=0;k<n;k++)b[k]=x[k];
  return 0;
}
