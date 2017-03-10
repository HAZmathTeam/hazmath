/*! \file src/utilities/vec.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 5/13/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/29/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/30/2016
 *
 */

#include "hazmath.h"

/***********************************************************************************************/
INT dvec_isnan (dvector *u)
{
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
    
    INT i;
    
    for ( i=0; i<u->row; i++ ) {
        if ( ISNAN(u->val[i]) ) return TRUE;
    }

    return FALSE;
}

/***********************************************************************************************/
INT ivec_isnan (ivector *u)
{
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

    INT i;

    for ( i=0; i<u->row; i++ ) {
        if ( ISNAN(u->val[i]) ) return TRUE;
    }

    return FALSE;
}

/***********************************************************************************************/
dvector dvec_create (const INT m)
{
    
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
    
    dvector u;
    
    u.row = m;
    u.val = (REAL *)calloc(m,sizeof(REAL));
    
    return u;

}

/***********************************************************************************************/
ivector ivec_create (const INT m)
{    
    
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
    
    ivector u;
    
    u.row = m;
    u.val = (INT *)calloc(m,sizeof(INT));
    
    return u;
}

/***********************************************************************************************/
void dvec_alloc (const INT m,
                 dvector *u)
{    
    /*!
     * \fn void dvec_alloc (const INT m, dvector *u)
     *
     * \brief Allocate a dvector of given length
     *
     * \param m    length of the dvector
     * \param u    Pointer to dvector (OUTPUT)
     *
     */
    
    u->row = m;
    u->val = (REAL*)calloc(m,sizeof(REAL));
    
    return;
}

/***********************************************************************************************/
void ivec_alloc (const INT m,
                 ivector *u)
{
    /*!
     * \fn void ivec_alloc (const INT m, ivector *u)
     *
     * \brief Allocate an ivector of given length
     *
     * \param m   length of the ivector
     * \param u   Pointer to ivector (OUTPUT)
     *
     */
    
    u->row = m;
    u->val = (INT*)calloc(m,sizeof(INT));
    
    return;
}

/***********************************************************************************************/
void dvec_free (dvector *u)
{    
    /*!
     * \fn void dvec_free (dvector *u)
     *
     * \brief Free the space of a dvector
     *
     * \param u   Pointer to dvector which needs to be deallocated
     *
     */
    
    if (u==NULL) return;
    
    free(u->val);
    u->row = 0; u->val = NULL; 
}

/***********************************************************************************************/
void ivec_free (ivector *u)
{    
    /*!
     * \fn void ivec_free (ivector *u)
     *
     * \brief Free the space of an ivector
     *
     * \param u   Pointer to ivector which needs to be deallocated
     *
     * \note This function is same as dvec_free except input type.
     */
    
    if (u==NULL) return;
    
    free(u->val);
    u->row = 0; u->val = NULL; 
}

/***********************************************************************************************/
void dvec_null (dvector *u)
{
    
    /*!
     * \fn void dvec_null (dvector *u)
     *
     * \brief Initialize a dvector to NULL
     *
     * \param u   Pointer to the dvector
     *
     */
    
    u->row = 0; u->val = NULL;
}

/***********************************************************************************************/
void ivec_null (ivector *u)
{

    /*!
     * \fn void ivec_null (ivector *u)
     *
     * \brief Initialize an ivector to NULL
     *
     * \param u   Pointer to the ivector
     *
     */

    u->row = 0; u->val = NULL;
}

/***********************************************************************************************/
void dvec_rand (const INT n,
                dvector *u)
{
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
        
    INT s=1; srand(s);
    
    INT i;
    
    u->row = n;
    
    for (i=0; i<n; ++i){
        u->val[i] = rand()/(RAND_MAX+1.0);
    }
}

/***********************************************************************************************/
void dvec_rand_true (const INT n,
                     dvector *u)
{
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

    srand(time(0));

    INT i;

    u->row = n;

    for (i=0; i<n; ++i){
        u->val[i] = rand()/(RAND_MAX+1.0);
    }


}

/***********************************************************************************************/
void dvec_set (INT n,
               dvector *u,
               REAL val)
{
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
void ivec_set (INT n,
               ivector *u,
               const INT val)
{
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
    
    INT i;

    if (n>0) u->row=n;
    else n=u->row;
    
    for (i=0; i<n; ++i) u->val[i]=val;
    
}

/***********************************************************************************************/
void dvec_cp (dvector *x,
              dvector *y)
{
    /*!
     * \fn void dvec_cp (dvector *x, dvector *y)
     *
     * \brief Copy dvector x to dvector y
     *
     * \param x  Pointer to dvector
     * \param y  Pointer to dvector (OUTPUT) (MODIFIED)
     *
     */

    y->row=x->row;
    memcpy(y->val,x->val,x->row*sizeof(REAL));
    
}

/***********************************************************************************************/
void ivec_cp (ivector *x,
              ivector *y)
{
    /*!
     * \fn void ivec_cp (ivector *x, ivector *y)
     *
     * \brief Copy ivector x to ivector y
     *
     * \param x  Pointer to ivector
     * \param y  Pointer to ivector (OUTPUT) (MODIFIED)
     *
     */

    y->row=x->row;
    memcpy(y->val,x->val,x->row*sizeof(INT));

}


/***********************************************************************************************/
REAL dvec_maxdiff (dvector *x,
                   dvector *y)
{
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
    
    const INT length=x->row;
    REAL Linf=0., diffi=0.;
    
    INT i;
    
    for (i=0; i<length; ++i) {
        if ((diffi = ABS(x->val[i]-y->val[i])) > Linf) Linf = diffi;
    }
    
    return Linf;
}

/***********************************************************************************************/
void dvec_ax (const REAL a,
              dvector *x)
{
    /*!
     * \fn void dvec_ax (const REAL a, dvector *x)
     *
     * \brief Compute x = a*x
     *
     * \param a   REAL scalar number
     * \param x   Pointer to dvector x (OUTPUT)
     *
     */

    INT i, m=x->row;

    if (a!=1.0){
        for (i=0; i<m; ++i) x->val[i] = a*x->val[i];
    }


}


/***********************************************************************************************/
void dvec_axpy (const REAL a,
                dvector *x,
                dvector *y)
{
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
    
    INT i, m=x->row;
    
    if ((y->row-m)!=0) {
        printf("### WARNING HAZMATH DANGER in function %s: Two vectors have different lengths!\n", __FUNCTION__);
        m = MIN(m, y->row);
        printf("Only first %d entries will be computed!!\n", m);
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
void dvec_axpyz(const REAL a,
                dvector *x,
                dvector *y,
                dvector *z)
{
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

    INT m=x->row;
    
    if ((y->row-m)!=0) {
        printf("### WARNING HAZMATH DANGER in function %s: Two vectors have different lengths!\n", __FUNCTION__);
        m = MIN(m, y->row);
        printf("Only first %d entries will be computed!!\n", m);
    }
    
    z->row = m;

    memcpy(z->val, y->val, m*sizeof(REAL));
    array_axpy(m, a, x->val, z->val);

}

/***********************************************************************************************/
REAL dvec_dotprod (dvector *x,
                   dvector *y)
{
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
    
    REAL value=0;
    INT i, length=x->row;

    if ((y->row-length)!=0) {
        printf("### WARNING HAZMATH DANGER in function %s: Two vectors have different lengths!\n", __FUNCTION__);
        length = MIN(length, y->row);
        printf("Only first %d entries will be computed!!\n", length);
    }
    
    for (i=0; i<length; ++i) value+=x->val[i]*y->val[i];
    
    return value;
}


/***********************************************************************************************/
REAL dvec_relerr (dvector *x,
                  dvector *y)
{
    
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
    
    REAL diff=0, temp=0;
    INT i;
    INT length=x->row;
    
    if (length!=y->row) {
        printf("### WARNING HAZMATH DANGER in function %s: Two vectors have different lengths!\n", __FUNCTION__);
        length = MIN(length, y->row);
        printf("Only first %d entries will be computed!!\n", length);
    }
    
    for (i=0;i<length;++i) {
        temp += x->val[i]*x->val[i];
        diff += pow(x->val[i]-y->val[i],2);
    }
    
    return sqrt(diff/temp);
}

/***********************************************************************************************/
REAL dvec_norm1 (dvector *x)
{
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
    
    REAL onenorm=0;
    INT i;
    const INT length=x->row;
    
    for (i=0;i<length;++i) onenorm+=ABS(x->val[i]);
    
    return onenorm;
}

/***********************************************************************************************/
REAL dvec_norm2 (dvector *x)
{
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
    
    REAL twonorm=0;
    INT i;
    const INT length=x->row;
    
    for (i=0;i<length;++i) twonorm+=x->val[i]*x->val[i];
    
    return sqrt(twonorm);
}

/***********************************************************************************************/
REAL dvec_norminf (dvector *x)
{
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
    
    INT i;
    const INT length=x->row;
    register REAL infnorm=0;
    
    for (i=0;i<length;++i) infnorm=MAX(infnorm,ABS(x->val[i]));
    
    return infnorm;
}

/********************************** END *************************************************/

