/*! \file src/utilities/array.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 5/13/15.
 *  Copyright 2016__HAZMATH__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/25/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/27/2016
 *
 */

#include "hazmath.h"

/****************************************************************************************/
/*!
 * \fn INT darray_max(INT n, const REAL *u)
 *
 * \brief finds the max element in a linear array (vector); returns
 *        its index
 *
 * \param n (INPUT) size of the array
 * \param u array of REALs (INPUT)
 *
 *
 * \author Clemens Hofreither 20210620
 */
INT darray_max(INT n, const REAL *u)
{
  assert(n > 0);
  INT i_max = 0;
  REAL u_max = u[0];
  for (INT i = 1; i < n; ++i)
  {
    if (u[i] > u_max) {
      i_max = i;
      u_max = u[i];
    }
  }
  return i_max;
}

/****************************************************************************************/
/*!
 * \fn INT darray_min(INT n, const REAL *u)
 *
 * \brief finds the min element in a linear array (vector); returns
 *        its index
 *
 * \param n (INPUT) size of the array
 * \param u array of REALs (INPUT)
 *
 *
 * \author Clemens Hofreither 20210620
 */
// index of the minimum of an array
INT darray_min(INT n, const REAL *u)
{
  assert(n > 0);
  INT i_min = 0;
  REAL u_min = u[0];
  for (INT i = 1; i < n; ++i)
  {
    if (u[i] < u_min) {
      i_min = i;
      u_min = u[i];
    }
  }
  return i_min;
}

/****************************************************************************************/
/*!
 * \fn INT iarray_max(INT n, const INT *u)
 *
 * \brief finds the max element in a linear array (vector); returns
 *        its index
 *
 * \param n (INPUT) size of the array
 * \param u array of INTs (INPUT)
 *
 *
 * \author Clemens Hofreither 20210620
 */
INT iarray_max(INT n, const INT *u)
{
  assert(n > 0);
  INT i_max = 0;
  INT u_max = u[0];
  for (INT i = 1; i < n; ++i)
  {
    if (u[i] > u_max) {
      i_max = i;
      u_max = u[i];
    }
  }
  return i_max;
}

/****************************************************************************************/
/*!
 * \fn INT iarray_min(INT n, const INT *u)
 *
 * \brief finds the min element in a linear array (vector); returns
 *        its index
 *
 * \param n (INPUT) size of the array
 * \param u array of INTs (INPUT)
 *
 *
 * \author Clemens Hofreither 20210620
 */
// index of the minimum of an array
INT iarray_min(INT n, const INT *u)
{
  assert(n > 0);
  INT i_min = 0;
  INT u_min = u[0];
  for (INT i = 1; i < n; ++i)
  {
    if (u[i] < u_min) {
      i_min = i;
      u_min = u[i];
    }
  }
  return i_min;
}

/***********************************************************************************************/
/*!
 * \fn INT locate0(INT needle, INT *haystack, INT n)
 *
 * \brief finds an element in an array. Returns the index in the array
 *        where the element is found. If not found, returns (-1).
 *
 * \param haystack array of integers (INPUT)
 * \param needle element to find (INPUT)
 *
 */
INT locate0(INT needle, const INT *haystack, INT n)
{
  /*   */
  INT i;
  for (i=0;i<n;i++)
    if(needle==haystack[i])
      return i;
  return -1;
}

/***********************************************************************************************/
/*!
 * \fn INT iarrary_locate(INT needle, const INT *a, INT n)
 *
 * \brief finds an element in an array. Returns the index in the array
 *        where the element is found. If not found, returns (-1).
 *
 * \param needle    element to find (INPUT)
 * \param a         array of integers (INPUT)
 * \param n         length of the array (INPUT)
 *
 */
INT iarray_locate (INT needle,
                   const INT *a,
                   INT n)
{
  /*   */
  INT i;
  for (i=0; i<n; i++) {

    if(needle == a[i]) {
      return i;
    }

  }

  return -1;
}

/***********************************************************************************************/
/*!
 * \fn void array_set (const INT n, REAL *x, const REAL val)
 *
 * \brief Set value for the first n entries of a REAL array to be a given value
 *
 * \param n    Number of entries
 * \param x    Pointer to the REAL array (OUTPUT)
 * \param val  given value
 *
 */
void array_set (const INT n,
                REAL *x,
                const REAL val)
{
    INT i;

    if (val == 0.0) {
        memset(x, 0x0, sizeof(REAL)*n);
    }
    else {
        for (i=0; i<n; ++i) x[i] = val;
    }
}

/***********************************************************************************************/
/*!
 * \fn void iarray_set (const INT n, INT *x, const INT val)
 *
 * \brief Set value for the first n entries of a INT array to be a given integer
 *
 * \param n    Number of entries
 * \param x    Pointer to the INT array (OUTPUT)
 * \param val  given integer
 *
 */
void iarray_set (const INT n,
                 INT *x,
                 const INT val)
{
    INT i;

    if (val == 0) {
        memset(x, 0, sizeof(INT)*n);
    }
    else {
        for (i=0; i<n; ++i) x[i]=val;
    }
}

/***********************************************************************************************/
/*!
 * \fn void array_cp (const INT n, REAL *x, REAL *y)
 *
 * \brief Copy first n entries of a REAL array x to another REAL array y
 *
 * \param n    Number of entires
 * \param x    Pointer to the original REAL array
 * \param y    Pointer to the destination REAL array (OUTPUT)
 *
 */
void array_cp (const INT n,
               REAL *x,
               REAL *y)
{
    memcpy(y, x, n*sizeof(REAL));
}

/***********************************************************************************************/
/*!
 * \fn void iarray_cp (const INT n, INT *x, INT *y)
 *
 * \brief Copy first n entries of a INT array to another INT array y
 *
 * \param n    Number of entries
 * \param x    Pointer to the original INT array
 * \param y    Pointer to the destination INT array (OUTPUT)
 *
 */
void iarray_cp (const INT n,
                INT *x,
                INT *y)
{
    memcpy(y, x, n*sizeof(INT));
}

/***********************************************************************************************/
/*!
 * \fn void array_shuffle(const INT n, REAL *x)
 *
 * \brief shuffle a REAL array
 *
 * \param n    Number of entries
 * \param x    Pointer to the REAL array (OUTPUT)
 *
 */
void array_shuffle(const INT n,
                   REAL *x)
{
    if (n > 1) {
        INT i, j;
        REAL temp;

        for (i = 0; i < n-1; i++) {
          j = i + rand() / (RAND_MAX / (n - i) + 1);
          temp = x[j];
          x[j] = x[i];
          x[i] = temp;
        }

    }


}

/***********************************************************************************************/
/*!
 * \fn void iarray_shuffle(const INT n, INT *x)
 *
 * \brief shuffle an INT array
 *
 * \param n    Number of entries
 * \param x    Pointer to the INT array (OUTPUT)
 *
 */
void iarray_shuffle(const INT n,
                   INT *x)
{
    if (n > 1) {
        INT i, j, temp;

        for (i = 0; i < n-1; i++) {
          j = i + rand() / (RAND_MAX / (n - i) + 1);
          temp = x[j];
          x[j] = x[i];
          x[i] = temp;
        }

    }


}

/***********************************************************************************************/
/*!
 * \fn void array_ax (const INT n, const REAL a, REAL *x)
 *
 * \brief Compute x = a*x
 *
 * \param n    length of the array
 * \param a    Scalar REAL number
 * \param x    Pointer to a REAL array (OUTPUT)
 *
 * \note Result REAL array overwrites the original REAL array x
 */
void array_ax (const INT n,
               const REAL a,
               REAL *x)
{
    INT i;

    if (a != 1.0) {
        for (i=0; i<n; ++i) x[i] *= a;
    }
}

/***********************************************************************************************/
/*!
 * \fn void array_axpy (const INT n, const REAL a, REAL *x, REAL *y)
 *
 * \brief Compute y = a*x + y
 *
 * \param n    Length of the arrays
 * \param a    Scalar REAL number
 * \param x    Pointer to the REAL array x
 * \param y    Pointer to the REAL array y (OUTPUT)
 *
 *
 * \note .
 */
void array_axpy (const INT n,
                 const REAL a,
                 REAL *x,
                 REAL *y)
{
    INT i;

    if (a==1.0) {
        for (i=0; i<n; ++i) y[i] += x[i];
    }
    else if (a==-1.0) {
        for (i=0; i<n; ++i) y[i] -= x[i];
    }
    else {
        for (i=0; i<n; ++i) y[i] += a*x[i];
    }
}

/***********************************************************************************************/
/*!
 * \fn void array_axpyz(const INT n, const REAL a, REAL *x,
 *                                REAL *y, REAL *z)
 *
 * \brief Compute z = a*x + y
 *
 * \param n    Length of the arrays
 * \param a    Scalar REAL number
 * \param x    Pointer to the REAL array x
 * \param y    Pointer to the REAL array y
 * \param z    Pointer to the REAL array z (OUTPUT)
 *
 */
void array_axpyz (const INT n,
                  const REAL a,
                  REAL *x,
                  REAL *y,
                  REAL *z)
{
    INT i;

    if (a==1.0){
        for (i=0; i<n; ++i) z[i] = x[i]+y[i];
    }
    else if (a==-1.0){
        for (i=0; i<n; ++i) z[i] = y[i]-x[i];
    }
    else {
        for (i=0; i<n; ++i) z[i] = a*x[i]+y[i];
    }

}

/***********************************************************************************************/
/*!
 * \fn void array_axpby (const INT n, const REAL a, REAL *x, const REAL b, REAL *y)
 *
 * \brief y = a*x + b*y
 *
 * \param n    Length of the arrays
 * \param a    Scalar REAL number
 * \param x    Pointer to the REAL array x
 * \param b    Scalar REAL number
 * \param y    Pointer to the REAL array y (OUTPUT)
 *
 * \note y is reused to store the resulting array.
 */
void array_axpby (const INT n,
                  const REAL a,
                  REAL *x,
                  const REAL b,
                  REAL *y)
{
    INT i;

    for (i=0; i<n; ++i) y[i] = a*x[i]+b*y[i];

}

/***********************************************************************************************/
/*!
 * \fn REAL array_dotprod (const INT n, const REAL *x, const REAL *y)
 *
 * \brief Compute the inner product of two REAL arrays x and y
 *
 * \param n    Length of the arrays
 * \param x    Pointer to the REAL array x
 * \param y    Pointer to the REAL array y
 *
 * \return     Inner product (x,y)
 *
 */
REAL array_dotprod (const INT n,
                    const REAL * x,
                    const REAL * y)
{

    INT i;
    REAL value = 0.0;

    for (i=0; i<n; ++i) {
      value += x[i]*y[i];
    }

    return value;
}

/***********************************************************************************************/
/*!
 * \fn REAL array_norm1 (const INT n, const REAL * x)
 *
 * \brief Compute the l1 norm of a REAL array x
 *
 * \param n    Number of variables
 * \param x    Pointer to x
 *
 * \return     l1 norm of x
 *
 */
REAL array_norm1 (const INT n,
                  const REAL * x)
{
    INT i;
    REAL onenorm = 0.;

    for (i=0;i<n;++i) onenorm+=ABS(x[i]);

    return onenorm;
}

/***********************************************************************************************/
/*!
 * \fn REAL array_norm2 (const INT n, const REAL * x)
 *
 * \brief compute l2 norm of a REAL array x
 *
 * \param n    Length of the REAL array x
 * \param x    Pointer to the REAL array x
 *
 * \return     l2 norm of x
 *
 */
REAL array_norm2 (const INT n,
                  const REAL * x)
{
    INT i;
    REAL twonorm = 0.;

    for (i=0;i<n;++i) twonorm+=x[i]*x[i];

    return sqrt(twonorm);
}

/***********************************************************************************************/
/**
 * \fn REAL array_norminf (const INT n, const REAL * x)
 *
 * \brief compute infinity norm of a REAL array x
 *
 * \param n    Length of the REAL array x
 * \param x    Pointer to the REAL array x
 *
 * \return     infinity norm of x
 *
 */
REAL array_norminf (const INT n,
                    const REAL * x)
{
    INT i;
    REAL infnorm = 0.0;

    for (i=0;i<n;++i) infnorm=MAX(infnorm,ABS(x[i]));

    return infnorm;
}

/***********************************************************************************************/
/*!
 * \fn REAL array_normp (const INT n, const REAL * x)
 *
 * \brief compute lp norm of a REAL array x
 *
 * \param n    Length of the REAL array x
 * \param x    Pointer to the REAL array x
 *
 * \return     lp norm of x
 *
 */
REAL array_normp(const INT n,
                 const REAL * x,
                 const INT p)
{
    INT i;
    REAL pnorm = 0.;

    for (i=0;i<n;++i) pnorm+=ABS(pow(x[i],p));

    return pow(pnorm, 1.0/p);
}

/****************************************************************************************/
/*!
 * \fn det3D(REAL *mydet,REAL* vec1,REAL* vec2,REAL* vec3)
 *
 * \brief compute determinant of 3 3D arrays
 *
 * \param mydet   Determinant of the 3 arrays (OUTPUT)
 * \param vec1    Pointer to the first array
 * \param vec2    Pointer to the second array
 * \param vec2    Pointer to the third array
 *
 */
void det3D(REAL *mydet,
           REAL* vec1,
           REAL* vec2,
           REAL* vec3)
{
  /* gets determinant of 3 3D vectors */
  REAL dettmp;

  dettmp = vec1[0]*(vec2[1]*vec3[2]-vec2[2]*vec3[1]) - vec1[1]*(vec2[0]*vec3[2]-vec2[2]*vec3[0]) + vec1[2]*(vec2[0]*vec3[1]-vec2[1]*vec3[0]);

  *mydet = dettmp;

  return;
}
/*!
* \fn void cross_product(REAL u, REAL v, REAL cross, REAL mag)
*
* \brief Computes the cross product of two vectors in R3 and its magnitude
*
* \param u 			vector in R3
* \param v 			vector in R3
*
* \return cross		  u x v
* \return mag 		| u x v |
*
*/
void cross_product(REAL* u, REAL* v, REAL* cross, REAL* mag)
{
  cross[0] = u[1]*v[2] - u[2]*v[1];
  cross[1] = u[2]*v[0] - u[0]*v[2];
  cross[2] = u[0]*v[1] - u[1]*v[0];

  *mag = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);

  return;
}

/*!
* \fn void ld2d(REAL *dest,REAL16 *src, const INT n)
*
* \brief Copies a REAL16 array onto REAL array with the same number of elements;
*
* \param dest 			destination array (REAL   *)
* \param src 			source      array (REAL16 *)
*
*/
void ld2d(REAL *dest,REAL16 *src, const INT n)
{
  //  long double to double;
  INT j;
  for (j=0;j<n;++j) dest[j]=(REAL )src[j];
  return;
}
/*!
* \fn void d2ld(REAL16 *dest,REAL *src, const INT n)
*
* \brief Copies a REAL16 array onto REAL array with the same number of elements;
*
* \param dest 			destination array (REAL   *)
* \param src 			source      array (REAL16 *)
*
*/
void d2ld(REAL16 *dest,REAL *src, const INT n)
{
  //  double to long double;
  INT j;
  for (j=0;j<n;++j) dest[j]=(REAL16 )src[j];
  return;
}
/*!
 * \fn INT array_uniq(const INT n,INT *a)
 *
 * \brief in a sorted array, removes all repeated entries and returns
 *        the size of the array with no repetitions counted.
 *
 * \param n  the size of the array a
 *
 * \param a on input a is a pointer to a sorted integer array, on
 *          output the first j entries of a are the array with no
 *          repetitions.
 *
 *
 */
INT array_uniq(const INT n,INT *a)
{
  // remove duplicates from a sorted array; the array size is changed
  INT i,j;
  if(n<=1) return n;
  if((n==2) && (a[0]==a[1])){
    return 1;
  } else if (n==2){
    return 2;
  }
  j=n-1;
  for(i=j;i>0;--i)
    a[i]=a[i]-a[i-1];
  j=1;
  for(i=1;i<n;++i)
    if(a[i])
      a[j++]=a[i];
  for(i=1;i<j;++i)
    a[i]=a[i-1]+a[i];
  return j;
}
/*!
 * \fn INT array_cap(const INT n1i,INT *a1i,	\
 *	       const INT n2i, INT *a2i,		\
 *	       INT *acap, INT *wrk)
 *
 * \brief find the intersection of two sorted arrays with no repeating entries
 *
 * \param n1  the size of the array a1i
 * \param a1i a pointer to the first sorted integer array with no
 *            repeated entries
 * \param n2  the size of the array a2i
 * \param a2i a pointer to the second sorted integer array with no
 *            repeated entries
 *
 * \param acap array with enough space to hold the intersection.
 * \param wrk working array of at least n1i+n2i entries to copy the
 *            original arrays.
 *
 * \return 0 if the intersection is empty or the size of the intersection.
 *
 *
 *
 */
INT array_cap(const INT n1i,const INT *a1i,	\
	       const INT n2i, const INT *a2i,		\
	       INT *acap, INT *wrk)
{ // wrk must be of size (n1i+n2i)
  INT *a1=wrk;
  INT *a2=wrk+n1i;
  /*First array*/
  memcpy(a1,a1i,n1i*sizeof(INT));
  isi_sort(n1i,a1);
  INT n1=array_uniq(n1i,a1);
  /*Second array*/
  memcpy(a2,a2i,n2i*sizeof(INT));
  isi_sort(n2i,a2);
  INT n2=array_uniq(n2i,a2);
  //////////////////////////////////////////////////
  INT  ncap=0, i=0, j=0;
  while ((i < n1) && (j < n2)) {
    if (a1[i] < a2[j])
      i++;
    else if (a1[i] > a2[j]){
      j++;
    } else {
      acap[ncap]=a2[j];
      ncap++;
      j++;i++;
    }
  }
  if(ncap) {
    for(i=0;i<ncap;++i){
      wrk[i]=-1;
      for(j=0;j<n1i;j++){
	if(acap[i]!=a1i[j]) continue;
	wrk[i]=j;
      }
    }
  }
  /* if(ncap) */
  /*   acap=realloc(acap,(ncap)*sizeof(INT)); */
  /* else{ */
  /*   acap=NULL; */
  /* } */
  return ncap;
}
/*************************************  END  ******************************************************/
