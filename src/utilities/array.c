/*
 *  array.c
 *
 *  Created by James Adler and Xiaozhe Hu on 5/13/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

#include <math.h>

// Our Includes
#include "macro.h"
#include "grid.h"
#include "sparse.h"
#include "vec.h"
#include "functs.h"
#include "fem.h"



/***********************************************************************************************/
void array_null (REAL *x)
{
    
    /**
     * \fn void array_null (REAL *x)
     *
     * \brief Initialize an array
     *
     * \param x    Pointer to the vector
     *
     */
    
    x = NULL;
}

/***********************************************************************************************/
void array_set (const INT n,
                REAL *x,
                const REAL val)
{
    
    /**
     * \fn void array_set (const INT n, REAL *x, const REAL val)
     *
     * \brief Set initial value for an array to be x=val
     *
     * \param n    Number of variables
     * \param x    Pointer to the vector
     * \param val  Initial value for the REAL array
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
    /**
     * \fn void iarray_set (const INT n, INT *x, const INT val)
     *
     * \brief Set initial value for an array to be x=val
     *
     * \param n    Number of variables
     * \param x    Pointer to the vector
     * \param val  Initial value for the REAL array
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
    /**
     * \fn void array_cp (const INT n, REAL *x, REAL *y)
     *
     * \brief Copy an array to the other y=x
     *
     * \param n    Number of variables
     * \param x    Pointer to the original vector
     * \param y    Pointer to the destination vector
     *
     */
    
    memcpy(y, x, n*sizeof(REAL));
}

/***********************************************************************************************/
void iarray_cp (const INT n, 
                INT *x,
                INT *y)
{
    /**
     * \fn void iarray_cp (const INT n, INT *x, INT *y)
     *
     * \brief Copy an array to the other y=x
     *
     * \param n    Number of variables
     * \param x    Pointer to the original vector
     * \param y    Pointer to the destination vector
     *
     */
    
    memcpy(y, x, n*sizeof(INT));
}

/***********************************************************************************************/
void array_ax (const INT n,
               const REAL a,
               REAL *x)
{
    /**
     * \fn void array_ax (const INT n, const REAL a, REAL *x)
     *
     * \brief x = a*x
     *
     * \param n    Number of variables
     * \param a    Factor a
     * \param x    Pointer to x
     *
     * \note x is reused to store the resulting array.
     */
    
    INT i;
    
    if (a == 1.0) {
        
    }
    else {
        for (i=0; i<n; ++i) x[i] *= a;
    }
}

/***********************************************************************************************/
void array_axpy (const INT n,
                 const REAL a,
                 REAL *x,
                 REAL *y)
{
    /**
     * \fn void array_axpy (const INT n, const REAL a, REAL *x, REAL *y)
     *
     * \brief y = a*x + y
     *
     * \param n    Number of variables
     * \param a    Factor a
     * \param x    Pointer to x
     * \param y    Pointer to y
     *
     *
     * \note y is reused to store the resulting array.
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
    /**
     * \fn void array_axpyz(const INT n, const REAL a, REAL *x,
     *                                REAL *y, REAL *z)
     *
     * \brief z = a*x + y
     *
     * \param n    Number of variables
     * \param a    Factor a
     * \param x    Pointer to x
     * \param y    Pointer to y
     * \param z    Pointer to z
     *
     */
    
    INT i;
    
    for (i=0; i<n; ++i) z[i] = a*x[i]+y[i];
    
}

/***********************************************************************************************/
void array_axpby (const INT n,
                  const REAL a,
                  REAL *x,
                  const REAL b,
                  REAL *y)
{
    /**
     * \fn void array_axpby (const INT n, const REAL a, REAL *x,
     *                                 const REAL b, REAL *y)
     *
     * \brief y = a*x + b*y
     *
     * \param n    Number of variables
     * \param a    Factor a
     * \param x    Pointer to x
     * \param b    Factor b
     * \param y    Pointer to y
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
    /**
     * \fn REAL array_dotprod (const INT n, const REAL *x, const REAL *y)
     *
     * \brief Inner product of two arraies (x,y)
     *
     * \param n    Number of variables
     * \param x    Pointer to x
     * \param y    Pointer to y
     *
     * \return     Inner product (x,y)
     *
     */
    
    INT i;
    REAL value = 0.0;
    
    for ( i=0; i<n; ++i ) value += x[i]*y[i];
    
    return value;
}

/***********************************************************************************************/
REAL array_norm1 (const INT n,
                  const REAL * x)
{
    /**
     * \fn REAL array_norm1 (const INT n, const REAL * x)
     *
     * \brief L1 norm of array x
     *
     * \param n    Number of variables
     * \param x    Pointer to x
     *
     * \return     L1 norm of x
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
    /**
     * \fn REAL array_norm2 (const INT n, const REAL * x)
     *
     * \brief L2 norm of array x
     *
     * \param n    Number of variables
     * \param x    Pointer to x
     *
     * \return     L2 norm of x
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
     * \brief Linf norm of array x
     *
     * \param n    Number of variables
     * \param x    Pointer to x
     *
     * \return     L_inf norm of x
     *
     */
    
    INT i;
    REAL infnorm = 0.0;
    
    for (i=0;i<n;++i) infnorm=MAX(infnorm,ABS(x[i]));
    
    return infnorm;
}



