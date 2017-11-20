/*! \file src/utilities/sorting.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 5/13/17.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *\note: sorting routines: straight insert sort. 
 */

#include "hazmath.h"
/********************************************************************/

/********************************************************************/
/*STRAIGHT INSERT sort */
void isi_sort(INT n, INT *a)
//implements STRAIGHT  INSERT sort for integers
{
  INT i,j,aj;
  for (j = 1; j < n; ++j) {
      aj = *(a + j);
      i = j - 1;
      while ((i >=0) && (aj<a[i])){
	  *(a + i + 1) = *(a + i);
	  --i;
      }
      *(a + i + 1) = aj; 
  }
  return;
}
void isi_sortp(const INT n, INT *a, INT *p, INT *invp)
{
//implements STRAIGHT INSERT sort for integers, returns prmutation and
//inverse permutation.
  INT i,j,aj,pj;
  for (i = 0; i < n; i++){p[i]=i;} 
  for (j = 1; j < n; j++){
    aj = *(a + j); pj = *(p + j);
    i = j - 1;
    while ((i >=0) && (aj<a[i])){
      *(a + i + 1) = *(a + i);  *(p+i+1)=*(p+i);
      --i;
    }
    *(a + i + 1) = aj; *(p+i+1)=pj;
  }
  for(i=0;i<n;i++){
    invp[p[i]]=i;
  }  
  for(i=0;i<n;i++){
    a[i]=a[invp[i]];
  }
  return;
}
void dsi_sortp(const INT n, REAL *a, INT *p, INT *invp)
//implements STRAIGHT INSERT sort for REALs, returns prmutation and
//inverse permutation.
{
  INT i,j,pj;
  REAL aj;
  for (i = 0; i < n; i++){p[i]=i;} 
  for (j = 1; j < n; j++){
    aj = *(a + j); pj = *(p + j);
    i = j - 1;
    while ((i >=0) && (aj<a[i])){
      *(a + i + 1) = *(a + i);  *(p+i+1)=*(p+i);
      --i;
    }
    *(a + i + 1) = aj; *(p+i+1)=pj;
  }
  for(i=0;i<n;i++){
    invp[p[i]]=i;
  }
  return;
}

