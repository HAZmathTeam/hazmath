/*! \file src/utilities/array.c
 *
 *  Created by James Adler and Xiaozhe Hu on 5/13/15.
 *  Copyright 2016__HAZMATH__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/25/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/27/2016
 *
 */

#include "hazmath.h"

/***********************************************************************************************/
void array_null(REAL *x)
{
    /*!
     * \fn void array_null (REAL *x)
     *
     * \brief Initialize a REAL null array
     *
     * \param x    Null pointer to the vector (OUTPUT)
     *
     */
    
    x = NULL;
}

/***********************************************************************************************/
void iarray_null(INT *x)
{
    /*!
     * \fn void iarray_null (INT *x)
     *
     * \brief Initialize a INT null array
     *
     * \param x    Null pointer to the vector (OUTPUT)
     *
     */

    x = NULL;
}

/***********************************************************************************************/
void array_set (const INT n,
                REAL *x,
                const REAL val)
{
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
    
    INT i;
    
    if (val == 0.0) {
        memset(x, 0x0, sizeof(REAL)*n);
    }
    else {
        for (i=0; i<n; ++i) x[i] = val;
    }
}

/***********************************************************************************************/
void iarray_set (const INT n,
                 INT *x,
                 const INT val)
{
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
    
    INT i;
    
    if (val == 0) {
        memset(x, 0, sizeof(INT)*n);
    }
    else {
        for (i=0; i<n; ++i) x[i]=val;
    }
}

/***********************************************************************************************/
void array_cp (const INT n,
               REAL *x,
               REAL *y)
{
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
    
    memcpy(y, x, n*sizeof(REAL));
}

/***********************************************************************************************/
void iarray_cp (const INT n, 
                INT *x,
                INT *y)
{
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
    
    memcpy(y, x, n*sizeof(INT));
}

/***********************************************************************************************/
void array_ax (const INT n,
               const REAL a,
               REAL *x)
{
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
    
    INT i;
    
    if (a != 1.0) {
        for (i=0; i<n; ++i) x[i] *= a;
    }
}

/***********************************************************************************************/
void array_axpy (const INT n,
                 const REAL a,
                 REAL *x,
                 REAL *y)
{
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
void array_axpyz (const INT n,
                  const REAL a,
                  REAL *x,
                  REAL *y,
                  REAL *z)
{
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
void array_axpby (const INT n,
                  const REAL a,
                  REAL *x,
                  const REAL b,
                  REAL *y)
{
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
    
    INT i;

    for (i=0; i<n; ++i) y[i] = a*x[i]+b*y[i];
    
}

/***********************************************************************************************/
REAL array_dotprod (const INT n,
                    const REAL * x,
                    const REAL * y)
{
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
    
    INT i;
    REAL value = 0.0;
    
    for (i=0; i<n; ++i) value += x[i]*y[i];
    
    return value;
}

/***********************************************************************************************/
REAL array_norm1 (const INT n,
                  const REAL * x)
{
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
    
    INT i;
    REAL onenorm = 0.;
    
    for (i=0;i<n;++i) onenorm+=ABS(x[i]);
    
    return onenorm;
}

/***********************************************************************************************/
REAL array_norm2 (const INT n,
                  const REAL * x)
{
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
    
    INT i;
    REAL twonorm = 0.;
    
    for (i=0;i<n;++i) twonorm+=x[i]*x[i];
    
    return sqrt(twonorm);
}

/***********************************************************************************************/
REAL array_norminf (const INT n,
                    const REAL * x)
{
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
    
    INT i;
    REAL infnorm = 0.0;
    
    for (i=0;i<n;++i) infnorm=MAX(infnorm,ABS(x[i]));
    
    return infnorm;
}

/***********************************************************************************************/
REAL array_normp (const INT n,
                  const REAL * x,
                  const INT p)
{
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

    INT i;
    REAL pnorm = 0.;

    for (i=0;i<n;++i) pnorm+=ABS(pow(x[i],p));

    return pow(pnorm, 1.0/p);
}

/****************************************************************************************/
void det3D(REAL *mydet,
           REAL* vec1,
           REAL* vec2,
           REAL* vec3)
{
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

  /* gets determinant of 3 3D vectors */
  REAL dettmp;

  dettmp = vec1[0]*(vec2[1]*vec3[2]-vec2[2]*vec3[1]) - vec1[1]*(vec2[0]*vec3[2]-vec2[2]*vec3[0]) + vec1[2]*(vec2[0]*vec3[1]-vec2[1]*vec3[0]);

  *mydet = dettmp;

  return;
}

/*************************************  END  ******************************************************/
