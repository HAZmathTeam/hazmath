/*! \file src/utilities/vec.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 5/13/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/29/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/30/2016 & 08/28/2021
 *
 */

#include "hazmath.h"

/***********************************************************************************************/
/*!
 * \fn INT dvec_isnan (dvector *u)
 *
 * \brief Check a dvector whether at least one entry is NAN
 *
 * \param u    Pointer to dvector
 *
 * \return     Return TRUE if there is an NAN
 *
 */
INT dvec_isnan (dvector *u)
{
    INT i;

    for ( i=0; i<u->row; i++ ) {
        if ( ISNAN(u->val[i]) ) return TRUE;
    }

    return FALSE;
}

/***********************************************************************************************/
/*!
 * \fn INT ivec_isnan (dvector *u)
 *
 * \brief Check a ivector whether at least one entry is NAN
 *
 * \param u    Pointer to ivector
 *
 * \return     Return TRUE if there is an NAN
 *
 */
INT ivec_isnan (ivector *u)
{
    INT i;

    for ( i=0; i<u->row; i++ ) {
        if ( ISNAN(u->val[i]) ) return TRUE;
    }

    return FALSE;
}

/***********************************************************************************************/
/*!
 * \fn dvector dvec_create (const INT m)
 *
 * \brief Create a dvector of given length
 *
 * \param m    length of the dvector
 *
 * \return u   The new dvector
 *
 */
dvector dvec_create (const INT m)
{
    dvector u;

    u.row = m;
    u.val = (REAL *)calloc(m,sizeof(REAL));

    return u;

}

/***********************************************************************************************/
/*!
 * \fn ivector ivec_create (const INT m)
 *
 * \brief Create an ivector of given length
 *
 * \param m   length of the ivector
 *
 * \return u  The new ivector
 *
 */
ivector ivec_create (const INT m)
{
    ivector u;

    u.row = m;
    u.val = (INT *)calloc(m,sizeof(INT));

    return u;
}

/***********************************************************************************************/
/*!
 * \fn void dvec_alloc (const INT m, dvector *u)
 *
 * \brief Allocate a dvector of given length
 *
 * \param m    length of the dvector
 * \param u    Pointer to dvector (OUTPUT)
 *
 */
void dvec_alloc (const INT m,
                 dvector *u)
{
    u->row = m;
    u->val = (REAL*)calloc(m,sizeof(REAL));

    return;
}

/***********************************************************************************************/
/*!
 * \fn void ivec_alloc (const INT m, ivector *u)
 *
 * \brief Allocate an ivector of given length
 *
 * \param m   length of the ivector
 * \param u   Pointer to ivector (OUTPUT)
 *
 */
void ivec_alloc (const INT m,
                 ivector *u)
{
    u->row = m;
    u->val = (INT*)calloc(m,sizeof(INT));

    return;
}

/***********************************************************************************************/
/*!
 * \fn void dvec_free (dvector *u)
 *
 * \brief Free the space of a dvector
 *
 * \param u   Pointer to dvector which needs to be deallocated
 *
 */
void dvec_free (dvector *u)
{
    if (u==NULL) return;

    free(u->val);
    u->row = 0; u->val = NULL;
}

/***********************************************************************************************/
/*!
 * \fn void ivec_free (ivector *u)
 *
 * \brief Free the space of an ivector
 *
 * \param u   Pointer to ivector which needs to be deallocated
 *
 * \note This function is same as dvec_free except input type.
 */
void ivec_free (ivector *u)
{
    if (u==NULL) return;

    free(u->val);
    u->row = 0; u->val = NULL;
}

/***********************************************************************************************/
/*!
 * \fn void dvec_null (dvector *u)
 *
 * \brief Initialize a dvector to NULL
 *
 * \param u   Pointer to the dvector
 *
 */
void dvec_null (dvector *u)
{
    u->row = 0; u->val = NULL;
}

/***********************************************************************************************/
/*!
 * \fn void ivec_null (ivector *u)
 *
 * \brief Initialize an ivector to NULL
 *
 * \param u   Pointer to the ivector
 *
 */
void ivec_null (ivector *u)
{
    u->row = 0; u->val = NULL;
}

/***********************************************************************************************/
/*!
 * \fn void dvec_rand (const INT n, dvector *x)
 *
 * \brief Generate random REAL vector in the range from 0 to 1
 *
 * \param n    Size of the vector
 * \param u    Pointer to dvector
 *
 * \note Use 1 as the seed of the random number generator, so the result dvector
 *       is actually repeatable.  -- Xiaozhe Hu
 *
 */
void dvec_rand (const INT n,
                dvector *u)
{
    INT s=1; srand(s);

    INT i;

    u->row = n;

    for (i=0; i<n; ++i){
        u->val[i] = rand()/(RAND_MAX+1.0);
    }
}

/***********************************************************************************************/
/*!
 * \fn void dvec_rand_true (const INT n, dvector *x)
 *
 * \brief Generate random REAL vector in the range from 0 to 1
 *
 * \param n    Size of the vector
 * \param u    Pointer to dvector
 *
 * \note Use time as the seed of the random number generator, so the result dvector
 *       is not repeatable.  -- Xiaozhe Hu
 *
 */
void dvec_rand_true (const INT n,
                     dvector *u)
{
    srand(time(0));

    INT i;

    u->row = n;

    for (i=0; i<n; ++i){
        u->val[i] = rand()/(RAND_MAX+1.0);
    }


}

/***********************************************************************************************/
/*!
 * \fn void dvec_set (INT n, dvector *x, REAL val)
 *
 * \brief Initialize dvector and set all the entries to be a given value
 *
 * \param n      Length of the dvector
 * \param u      Pointer to dvector (OUTPUT)
 * \param val    Given REAL value
 *
 * \note Memory space has to be allocated before calling this function -- Xiaozhe Hu
 *
 */
void dvec_set (INT n,
               dvector *u,
               REAL val)
{
    INT i;

    if (n>0) u->row=n;
    else n=u->row;

    if (val == 0.0) {
        memset(u->val, 0x0, sizeof(REAL)*n);
    }
    else {
        for (i=0; i<n; ++i) u->val[i]=val;
    }
}

/***********************************************************************************************/
/*!
 * \fn void ivec_set (INT n, ivector *u, const INT val)
 *
 * \brief Intialzie an ivector and set all the entries to be a given value
 *
 * \param  n    Length of the ivector
 * \param  u    Pointer to ivector
 * \param  val  Given integer value
 *
 * \note Memory space has to be allocated before calling this function -- Xiaozhe Hu
 *
 */
void ivec_set (INT n,
               ivector *u,
               const INT val)
{
    INT i;

    if (n>0) u->row=n;
    else n=u->row;

    for (i=0; i<n; ++i) u->val[i]=val;

}

/***********************************************************************************************/
/*!
 * \fn void dvec_cp (dvector *x, dvector *y)
 *
 * \brief Copy dvector x to dvector y
 *
 * \param x  Pointer to dvector
 * \param y  Pointer to dvector (OUTPUT) (MODIFIED)
 *
 */
void dvec_cp (dvector *x,
              dvector *y)
{
    y->row=x->row;
    memcpy(y->val,x->val,x->row*sizeof(REAL));

}

/***********************************************************************************************/
/*!
 * \fn void ivec_cp (ivector *x, ivector *y)
 *
 * \brief Copy ivector x to ivector y
 *
 * \param x  Pointer to ivector
 * \param y  Pointer to ivector (OUTPUT) (MODIFIED)
 *
 */
void ivec_cp (ivector *x,
              ivector *y)
{
    y->row=x->row;
    memcpy(y->val,x->val,x->row*sizeof(INT));

}


/***********************************************************************************************/
/*!
 * \fn REAL dvec_maxdiff (dvector *x, dvector *y)
 *
 * \brief Maximal difference of two dvector x and y (in absolute value)
 *
 * \param  x    Pointer to dvector
 * \param  y    Pointer to dvector
 *
 * \return      l-infinity norm of x-y
 *
 */
REAL dvec_maxdiff (dvector *x,
                   dvector *y)
{
    const INT length=x->row;
    REAL Linf=0., diffi=0.;

    INT i;

    for (i=0; i<length; ++i) {
        if ((diffi = ABS(x->val[i]-y->val[i])) > Linf) Linf = diffi;
    }

    return Linf;
}

/***********************************************************************************************/
/*!
 * \fn void dvec_ax (const REAL a, dvector *x)
 *
 * \brief Compute x = a*x
 *
 * \param a   REAL scalar number
 * \param x   Pointer to dvector x (OUTPUT)
 *
 */
void dvec_ax (const REAL a,
              dvector *x)
{
    INT i, m=x->row;

    if (a!=1.0){
        for (i=0; i<m; ++i) x->val[i] = a*x->val[i];
    }


}


/***********************************************************************************************/
/**
 * \fn void dvec_axpy (const REAL a, dvector *x, dvector *y)
 *
 * \brief Compute y = a*x + y
 *
 * \param a   REAL scalar number
 * \param x   Pointer to dvector x
 * \param y   Pointer to dvector y (OUTPUT)
 *
 */
void dvec_axpy (const REAL a,
                dvector *x,
                dvector *y)
{
    INT i, m=x->row;

    if ((y->row-m)!=0) {
        printf("### WARNING HAZMATH DANGER in function %s: Two vectors have different lengths!\n", __FUNCTION__);
        m = MIN(m, y->row);
        printf("Only first %lld entries will be computed!!\n", (long long )m);
    }

    if (a==1.0){
        for (i=0; i<m; ++i) y->val[i] += x->val[i];
    }
    else if(a==-1.0){
        for (i=0; i<m; ++i) y->val[i] -= x->val[i];
    }
    else {
        for (i=0; i<m; ++i) y->val[i] += a*x->val[i];
    }

}


/***********************************************************************************************/
/*!
 * \fn void dvec_axpyz (const REAL a, dvector *x, dvector *y, dvector *z)
 *
 * \brief Compute z = a*x + y, z is a third vector
 *
 * \param a   REAL factor a
 * \param x   Pointer to dvector x
 * \param y   Pointer to dvector y
 * \param z   Pointer to dvector z (OUTPUT)
 *
 * \note Memory of the dvector z shoule be allocated before calling this function -- Xiaozhe
 *
 */
void dvec_axpyz(const REAL a,
                dvector *x,
                dvector *y,
                dvector *z)
{
    INT m=x->row;

    if ((y->row-m)!=0) {
        printf("### WARNING HAZMATH DANGER in function %s: Two vectors have different lengths!\n", __FUNCTION__);
        m = MIN(m, y->row);
        printf("Only first %lld entries will be computed!!\n",(long long )m);
    }

    z->row = m;

    memcpy(z->val, y->val, m*sizeof(REAL));
    array_axpy(m, a, x->val, z->val);

}

/***********************************************************************************************/
/*!
 * \fn REAL dvec_dotprod (dvector *x, dvector *y)
 *
 * \brief Inner product of two dvectors x and y
 *
 * \param x   Pointer to dvector x
 * \param y   Pointer to dvector y
 *
 * \return Inner product of x and y: (x,y)
 *
 */
REAL dvec_dotprod (dvector *x,
                   dvector *y)
{
    REAL value=0;
    INT i, length=x->row;

    if ((y->row-length)!=0) {
        printf("### WARNING HAZMATH DANGER in function %s: Two vectors have different lengths!\n", __FUNCTION__);
        length = MIN(length, y->row);
        printf("Only first %lld entries will be computed!!\n", (long long )length);
    }

    for (i=0; i<length; ++i) value+=x->val[i]*y->val[i];

    return value;
}


/***********************************************************************************************/
/*!
 * \fn REAL dvec_relerr (dvector *x, dvector *y)
 *
 * \brief Relative differences between two dvector x and y (l2 norm)
 *
 * \param x   Pointer to dvector x
 * \param y   Pointer to dvector y
 *
 * \return relative differences: ||x-y||_2/||x||_2
 *
 */
REAL dvec_relerr (dvector *x,
                  dvector *y)
{
    REAL diff=0, temp=0;
    INT i;
    INT length=x->row;

    if (length!=y->row) {
        printf("### WARNING HAZMATH DANGER in function %s: Two vectors have different lengths!\n", __FUNCTION__);
        length = MIN(length, y->row);
        printf("Only first %lld entries will be computed!!\n", (long long )length);
    }

    for (i=0;i<length;++i) {
        temp += x->val[i]*x->val[i];
        diff += pow(x->val[i]-y->val[i],2);
    }

    return sqrt(diff/temp);
}

/***********************************************************************************************/
/*!
 * \fn REAL dvec_norm1 (dvector *x)
 *
 * \brief Compute l1 norm of a dvector x
 *
 * \param x   Pointer to dvector x
 *
 * \return l1 norm of x: ||x||_1
 *
 */
REAL dvec_norm1 (dvector *x)
{
    REAL onenorm=0;
    INT i;
    const INT length=x->row;

    for (i=0;i<length;++i) onenorm+=ABS(x->val[i]);

    return onenorm;
}

/***********************************************************************************************/
/*!
 * \fn REAL dvec_norm2 (dvector *x)
 *
 * \brief Compute l2 norm of dvector x
 *
 * \param x   Pointer to dvector x
 *
 * \return L2 norm of x: ||x||_2
 *
 */
REAL dvec_norm2 (dvector *x)
{
    REAL twonorm=0;
    INT i;
    const INT length=x->row;

    for (i=0;i<length;++i) twonorm+=x->val[i]*x->val[i];

    return sqrt(twonorm);
}

/***********************************************************************************************/
/*!
 * \fn REAL dvec_norminf (dvector *x)
 *
 * \brief Compute l_inf norm of dvector x
 *
 * \param x   Pointer to dvector x
 *
 * \return l_inf norm of x: ||x||_{inf}
 *
 */
REAL dvec_norminf (dvector *x)
{
    INT i;
    const INT length=x->row;
    register REAL infnorm=0;

    for (i=0;i<length;++i) infnorm=MAX(infnorm,ABS(x->val[i]));

    return infnorm;
}

/***********************************************************************************************/
/*!
 * \fn void dvec_orthog_const (dvector *x)
 *
 * \brief make x orthgonal to constant vector
 *
 * \param x   Pointer to dvector x
 *
 * \return x new x that is orthogonal to constant vector
 *
 */
void dvec_orthog_const(dvector *x)
{

  // local variables
  INT i;
  REAL sum, average;

  // compute sum
  sum = 0.0;
  for (i=0; i<x->row; i++) sum = sum+x->val[i];

  // compute average
  average = sum/x->row;

  // orthogonalize
  for (i=0; i<x->row; i++) x->val[i] = x->val[i]-average;

}

/***********************************************************************************************/
/*!
 * \fn void dvec_orthog(dvector *x, dvector *y)
 *
 * \brief make x orthgonal to y: x = x - (x'*y)*y;
 *
 * \param x   pointer to dvector x
 * \param y   pointer to dvector y (normalized)
 *
 * \return x  new x that is orthogonal to y
 *
 */
void dvec_orthog(dvector *x, dvector *y)
{

  // local variable
  REAL alpha = dvec_dotprod(x, y);

  // make it orthogonal
  dvec_axpy(-alpha, y, x);

}


/***********************************************************************************************/
/*!
 * \fn void dvec_inv (dvector *x)
 *
 * \brief Invert every element of dvector (0 -> 0 but produces a warning)
 *
 * \param x   Pointer to dvector x
 *
 * \return new x with x[i] = 1/x[i] for each i (if x[i] != 0)
 *
 */
void dvec_inv(dvector *x)
{
    // local variables
    INT i;
    const INT length=x->row;

    // invert all x_i != 0
    for (i=0;i<length;++i){
        if ( x->val[i] < SMALLREAL) printf("### WARNING HAZMATH DANGER in function %s: Cannot invert vector element of value zero!\n", __FUNCTION__);
        else x->val[i] = 1.0/x->val[i];
    }

}

/***********************************************************************************************/

/********************************** END *************************************************/
