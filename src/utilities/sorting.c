/*! \file src/utilities/sorting.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 5/13/17.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note: sorting routines: straight insert sort.
 *  \note: done cleanup for releasing -- Xiaozhe Hu 08/28/2021
 */
/********************************************************************/
#include "hazmath.h"

/***************************************************************/
/*!
 * \fn void dlexsort(const INT nr, const INT nc,REAL *a, INT *p)
 *
 * \brief implements STRAIGHT INSERT sorting to order lexicographically nr
 names with nc components each.
 *
 * \param nr    Number of names
 * \param nc    Number of components
 * \param a     Pointer to the array (OUTPUT)
 * \param p     Pointer to the permutation
 *
 *
 */
void dlexsort(const INT nr, const INT nc,REAL *a,INT *p)
{
  /* REAL NUMBERS; identical below for ints. It can be combined later
     --ltz 20190625
     *
     implements STRAIGHT INSERT sorting to order lexicographically nr
     names with nc components each.  the array a is overwritten.  aj is
     a working double array with nc (this will be dim in our most
     common case) elements. used to sort coordinates of the vertices of
     the macroelements lexicographically.  on output, a is ordered, p
     is the permutation used to order a. The ogiginal a is recovered
     with inv permutation aorig[]=a[invp[i]] where invp[p[i]]=i;
  */
  INT i,j,k,pj;
  unsigned INT lt=0;
  REAL *aj=(REAL *)calloc(nc,sizeof(REAL));
  for (i = 0; i < nr; i++){p[i]=i;}
  for (j = 1; j < nr; ++j) {
    //    for(k=0;k<nc;k++)aj[k] = a[j*nc+k];
    memcpy(aj,(a+j*nc),nc*sizeof(REAL)); pj = *(p + j);
    i = j - 1;
    lt=0;
    for(k=0;k<nc;k++) {
      if(aj[k]>a[i*nc+k]){lt=0;break;}
      else if(aj[k]<a[i*nc+k]){lt=1;break;}
      else {continue;}
    }
    //    while ((i >=0) && (aj<a[i])){
    while ((i>=0) && (lt)){
      memcpy((a+(i+1)*nc),(a+i*nc),nc*sizeof(REAL)); *(p+i+1)=*(p+i);
      //      for(k=0;k<nc;k++) a[(i + 1)*nc+k] = a[i*nc+k];
      --i;
      if(i<0) break;
      lt=0;
      for(k=0;k<nc;k++) {
	if(aj[k]>a[i*nc+k]){lt=0;break;}
	else if(aj[k]<a[i*nc+k]){lt=1;break;}
	else{continue;}
      }
    }
    memcpy((a+(i+1)*nc),aj,nc*sizeof(REAL)); *(p+i+1)=pj;
  }
  if(aj) free(aj);
  return;
}

/*********************************************************************/
/*!
 * \fn void ilexsort(const INT nr, const INT nc, INT *a, INT *p)
 *
 * \brief implements STRAIGHT INSERT sorting to order lexicographically nr
 names with nc components each.
 *
 * \param nr    Number of names
 * \param nc    Number of components
 * \param a     Pointer to the integer array (OUTPUT)
 * \param p     Pointer to the permutation
 *
 *
 */
void ilexsort(const INT nr, const INT nc,INT *a,INT *p)
{
  /* INT NUMBERS; identical above for reals. It can be combined later into one.
     --ltz 20190625
     *
     implements STRAIGHT INSERT sorting to order lexicographically nr
     names with nc components each.  the array a is overwritten.  aj is
     a working INT array with nc (this will be dim in our most
     common case) elements. used to sort coordinates of the vertices of
     the macroelements lexicographically.  on output, a is ordered, p
     is the permutation used to order a. The ogiginal a is recovered
     with inv permutation aorig[]=a[invp[i]] where invp[p[i]]=i;
  */
  INT i,j,k,pj;
  unsigned INT lt=0;
  INT *aj=(INT *)calloc(nc,sizeof(INT));
  for (i = 0; i < nr; i++){p[i]=i;}
  for (j = 1; j < nr; ++j) {
    //    for(k=0;k<nc;k++)aj[k] = a[j*nc+k];
    memcpy(aj,(a+j*nc),nc*sizeof(INT)); pj = *(p + j);
    i = j - 1;
    lt=0;
    for(k=0;k<nc;k++) {
      if(aj[k]>a[i*nc+k]){lt=0;break;}
      else if(aj[k]<a[i*nc+k]){lt=1;break;}
      else {continue;}
    }
    //    while ((i >=0) && (aj<a[i])){
    while ((i>=0) && (lt)){
      memcpy((a+(i+1)*nc),(a+i*nc),nc*sizeof(INT)); *(p+i+1)=*(p+i);
      //      for(k=0;k<nc;k++) a[(i + 1)*nc+k] = a[i*nc+k];
      --i;
      if(i<0) break;
      lt=0;
      for(k=0;k<nc;k++) {
	if(aj[k]>a[i*nc+k]){lt=0;break;}
	else if(aj[k]<a[i*nc+k]){lt=1;break;}
	else{continue;}
      }
    }
    memcpy((a+(i+1)*nc),aj,nc*sizeof(INT)); *(p+i+1)=pj;
  }
  if(aj) free(aj);
  return;
}

/***************************************************************/
/*STRAIGHT INSERT sort */
/*!
 * \fn void isi_sort(INT n, INT *a)
 *
 * \brief implements STRAIGHT INSERT sorting
 *
 * \param n    Number of components
 * \param a    Pointer to the integer array (OUTPUT)
 *
 *
 */
void isi_sort(INT n, INT *a)
{
//implements STRAIGHT  INSERT sort for integers
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

/***************************************************************/
/*!
 * \fn void dsi_sort(INT n, REAL *a)
 *
 * \brief implements STRAIGHT INSERT sorting
 *
 * \param n    Number of components
 * \param a    Pointer to the REAL array (OUTPUT)
 *
 *
 */
void dsi_sort(INT n, REAL *a)
{
//implements STRAIGHT  INSERT sort for REAL
  INT i,j;
  REAL aj;
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

/***************************************************************/
/*!
 * \fn void qsi_sort(INT n, REAL16 *a)
 *
 * \brief implements STRAIGHT INSERT sorting
 *
 * \param n    Number of components
 * \param a    Pointer to the REAL16 array (OUTPUT)
 *
 *
 */
void qsi_sort(INT n, REAL16 *a)
{
//implements STRAIGHT  INSERT sort for REAL
  INT i,j;
  REAL16 aj;
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

/***************************************************************/
/*!
 * \fn void isi_sortp(const INT n, INT *a, INT *p, INT *invp)
 *
 * \brief //implements STRAIGHT INSERT sort for integers, returns prmutation and
 //inverse permutation.
 *
 * \param n     Number of components
 * \param a     Pointer to the integer array
 * \param p     Pointer to the permutation
 * \param invp  Pointer to the inverse permutation
 *
 *
 */
void isi_sortp(const INT n, INT *a, INT *p, INT *invp)
{
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
  return;
}

/***************************************************************/
/*!
 * \fn void dsi_sortp(const INT n, REAL *a, INT *p, INT *invp)
 *
 * \brief implements STRAIGHT INSERT sort for REALs, returns prmutation and
 //inverse permutation.
 *
 * \param n     Number of components
 * \param a     Pointer to the REAL array
 * \param p     Pointer to the permutation
 * \param invp  Pointer to the inverse permutation
 *
 *
 */
void dsi_sortp(const INT n, REAL *a, INT *p, INT *invp)
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
