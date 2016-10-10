/*! \file src/utilities/vec.c
 *
 *  Created by James Adler and Xiaozhe Hu on 5/13/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

#include "hazmat.h"

/***********************************************************************************************/
INT dvec_isnan (dvector *u)
{
    /**
     * \fn INT dvec_isnan (dvector *u)
     *
     * \brief Check a dvector whether there is NAN
     *
     * \param u    Pointer to dvector
     *
     * \return     Return TRUE if there is NAN
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
    
    /**
     * \fn dvector dvec_create (const INT m)
     *
     * \brief Create dvector data space of REAL type
     *
     * \param m    Number of rows
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
    
    /**
     * \fn ivector ivec_create (const INT m)
     *
     * \brief Create vector data space of INT type
     *
     * \param m   Number of rows
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
    /**
     * \fn void dvec_alloc (const INT m, dvector *u)
     *
     * \brief Create dvector data space of REAL type
     *
     * \param m    Number of rows
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
    /**
     * \fn void ivec_alloc (const INT m, ivector *u)
     *
     * \brief Create vector data space of INT type
     *
     * \param m   Number of rows
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
    /**
     * \fn void dvec_free (dvector *u)
     *
     * \brief Free vector data space of REAL type
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
    /**
     * \fn void ivec_free (ivector *u)
     *
     * \brief Free vector data space of INT type
     *
     * \param u   Pointer to ivector which needs to be deallocated
     *
     *
     * \note This function is same as dvec_free except input type.
     */
    
    if (u==NULL) return;
    
    free(u->val);
    u->row = 0; u->val = NULL; 
}

/***********************************************************************************************/
void dvec_null (dvector *x)
{
    
    /**
     * \fn void dvec_null (dvector *x)
     *
     * \brief Initialize dvector
     *
     * \param x   Pointer to dvector which needs to be initialized
     *
     */
    
    x->row = 0; x->val = NULL;
}

/***********************************************************************************************/
void dvec_rand (const INT n,
                dvector *x)
{
    /**
     * \fn void dvec_rand (const INT n, dvector *x)
     *
     * \brief Generate random REAL vector in the range from 0 to 1
     *
     * \param n    Size of the vector
     * \param x    Pointer to dvector
     *
     */
    
    const INT va=(REAL) 0;
    const INT vb=(REAL) n;
    
    INT s=1; srand(s);
    
    INT i,j;
    
    x->row = n;
    
    for (i=0; i<n; ++i){
        j = 1 + (INT) (((REAL)n)*rand()/(RAND_MAX+1.0));
        x->val[i] = (((REAL)j)-va)/(vb-va);
    }
}

/***********************************************************************************************/
void dvec_set (INT n,
               dvector *x,
               REAL val)
{
    /**
     * \fn void dvec_set (INT n, dvector *x, REAL val)
     *
     * \brief Initialize dvector x[i]=val for i=0:n-1
     *
     * \param n      Number of variables
     * \param x      Pointer to dvector
     * \param val    Initial value for the vector
     *
     */
    
    INT i;
    REAL *xpt=x->val;
    
    if (n>0) x->row=n; 
    else n=x->row;
   
    if (val == 0.0) {
        memset(xpt, 0x0, sizeof(REAL)*n);
    }
    
    else {
        for (i=0; i<n; ++i) xpt[i]=val;
    }
}

/***********************************************************************************************/
void ivec_set (const INT m,
               ivector *u)
{
    /**
     * \fn void ivec_set (const INT m, ivector *u)
     *
     * \brief Set ivector value to be m
     *
     * \param  m    Integer value of ivector
     * \param  u    Pointer to ivector (MODIFIED)
     *
     */
    
    INT i;
    INT n = u->row;
    
    for (i=0; i<n; ++i) u->val[i] = m;
    
}

/***********************************************************************************************/
void dvec_cp (dvector *x,
                   dvector *y)
{
    /**
     * \fn void dvec_cp (dvector *x, dvector *y)
     *
     * \brief Copy dvector x to dvector y
     *
     * \param x  Pointer to dvector
     * \param y  Pointer to dvector (MODIFIED)
     *
     */

    y->row=x->row;
    memcpy(y->val,x->val,x->row*sizeof(REAL));
    
}

/***********************************************************************************************/
REAL dvec_maxdiff (dvector *x,
                   dvector *y)
{
    /**
     * \fn REAL dvec_maxdiff (dvector *x, dvector *y)
     *
     * \brief Maximal difference of two dvector x and y
     *
     * \param  x    Pointer to dvector
     * \param  y    Pointer to dvector
     *
     * \return      Maximal norm of x-y
     *
     */
    
    const INT length=x->row;
    REAL Linf=0., diffi=0.;
    REAL *xpt=x->val, *ypt=y->val;
    
    INT i;
    
    for (i=0; i<length; ++i) {
        if ((diffi = ABS(xpt[i]-ypt[i])) > Linf) Linf = diffi;
    }
    
    return Linf;
}

/***********************************************************************************************/
void dvec_symdiagscale (dvector *b,
                             dvector *diag)
{
    
    /**
     * \fn void dvec_symdiagscale (dvector *b, dvector *diag)
     *
     * \brief Symmetric diagonal scaling D^{-1/2}b
     *
     * \param b       Pointer to dvector
     * \param diag    Pointer to dvector: the diagonal entries
     *
     */
    
    // information about dvector
    const INT n = b->row;
    REAL *val = b->val;
    
    // local variables
    INT i;
    
    if (diag->row != n) {
        printf("### ERROR: Sizes of diag = %d != dvector = %d!", diag->row, n);
    }
    
    for (i=0; i<n; ++i) val[i] = val[i]/sqrt(diag->val[i]);
    
    return;
}

/***********************************************************************************************/
void dvec_axpy (const REAL a,
                dvector *x,
                dvector *y)
{
    /**
     * \fn void dvec_axpy (const REAL a, dvector *x, dvector *y)
     *
     * \brief y = a*x + y
     *
     * \param a   REAL factor a
     * \param x   Pointer to dvector x
     * \param y   Pointer to dvector y
     *
     */
    
    INT i, m=x->row;
    REAL *xpt=x->val, *ypt=y->val;
    
    if ((y->row-m)!=0) {
        printf("### ERROR: Two vectors have different dimensions!\n");
    }
    
    for (i=0; i<m; ++i) ypt[i] += a*xpt[i];
    
}

/***********************************************************************************************/
void dvec_ax (const REAL a,
                dvector *x)
{
    /**
     * \fn void dvec_ax (const REAL a, dvector *x)
     *
     * \brief x = a*x
     *
     * \param a   REAL factor a
     * \param x   Pointer to dvector x
     *
     */
    
    INT i, m=x->row;
    REAL *xpt=x->val;
    
    
    for (i=0; i<m; ++i) xpt[i] = a*xpt[i];
    
}

/***********************************************************************************************/
void dvec_axpyz(const REAL a,
                dvector *x,
                dvector *y,
                dvector *z)
{
    /**
     * \fn void dvec_axpyz (const REAL a, dvector *x, dvector *y, dvector *z)
     *
     * \brief z = a*x + y, z is a third vector (z is cleared)
     *
     * \param a   REAL factor a
     * \param x   Pointer to dvector x
     * \param y   Pointer to dvector y
     * \param z   Pointer to dvector z
     *
     * note: something wrong with the function!! -- Xiaozhe 
     *
     */
    
    const INT m=x->row;
    REAL *xpt=x->val, *ypt=y->val, *zpt=z->val;
    
    if ((y->row-m)!=0) {
        printf("### ERROR: Two vectors have different dimensions!\n");
    }
    
    z->row = m;
    
    memcpy(ypt, zpt, m*sizeof(REAL));
    array_axpy(m, a, xpt, zpt);
}

/***********************************************************************************************/
REAL dvec_dotprod (dvector *x,
                   dvector *y)
{
    /**
     * \fn REAL dvec_dotprod (dvector *x, dvector *y)
     *
     * \brief Inner product of two vectors (x,y)
     *
     * \param x   Pointer to dvector x
     * \param y   Pointer to dvector y
     *
     * \return Inner product
     *
     */
    
    REAL value=0;
    INT i;
    const INT length=x->row;
    REAL *xpt=x->val, *ypt=y->val;
    
    for (i=0; i<length; ++i) value+=xpt[i]*ypt[i];
    
    return value;
}


/***********************************************************************************************/
REAL dvec_relerr (dvector *x,
                  dvector *y)
{
    
    /**
     * \fn REAL dvec_relerr (dvector *x, dvector *y)
     *
     * \brief Relative error of two dvector x and y
     *
     * \param x   Pointer to dvector x
     * \param y   Pointer to dvector y
     *
     * \return relative error ||x-y||/||x||
     *
     */
    
    REAL diff=0, temp=0;
    INT i;
    const INT length=x->row;
    REAL *xpt=x->val, *ypt=y->val;
    
    if (length!=y->row) {
        printf("### ERROR: Two vectors have different dimensions!\n");
    }
    
    for (i=0;i<length;++i) {
        temp += xpt[i]*xpt[i];
        diff += pow(xpt[i]-ypt[i],2);
    }
    
    return sqrt(diff/temp);
}

/***********************************************************************************************/
REAL dvec_norm1 (dvector *x)
{
    /**
     * \fn REAL dvec_norm1 (dvector *x)
     *
     * \brief L1 norm of dvector x
     *
     * \param x   Pointer to dvector x
     *
     * \return L1 norm of x
     *
     */
    
    REAL onenorm=0;
    INT i;
    const INT length=x->row;
    REAL *xpt=x->val;
    
    for (i=0;i<length;++i) onenorm+=ABS(xpt[i]);
    
    return onenorm;
}

/***********************************************************************************************/
REAL dvec_norm2 (dvector *x)
{
    /**
     * \fn REAL dvec_norm2 (dvector *x)
     *
     * \brief L2 norm of dvector x
     *
     * \param x   Pointer to dvector x
     *
     * \return L2 norm of x
     *
     */
    
    REAL twonorm=0;
    INT i;
    const INT length=x->row;
    REAL *xpt=x->val;
    
    for (i=0;i<length;++i) twonorm+=xpt[i]*xpt[i];
    
    return sqrt(twonorm);
}


/***********************************************************************************************/
REAL dvec_norminf (dvector *x)
{
    /**
     * \fn REAL dvec_norminf (dvector *x)
     *
     * \brief Linf norm of dvector x
     *
     * \param x   Pointer to dvector x
     *
     * \return L_inf norm of x
     *
     */
    
    INT i;
    const INT length=x->row;
    REAL *xpt=x->val;
    
    register REAL infnorm=0;
    
    for (i=0;i<length;++i) infnorm=MAX(infnorm,ABS(xpt[i]));
    
    return infnorm;
}

/****************************************************************************************/
void det3D(REAL *mydet,REAL* vec1,REAL* vec2,REAL* vec3)          
{
  /* gets determinant of 3 3D vectors */
  REAL dettmp;
  
  dettmp = vec1[0]*(vec2[1]*vec3[2]-vec2[2]*vec3[1]) - vec1[1]*(vec2[0]*vec3[2]-vec2[2]*vec3[0]) + vec1[2]*(vec2[0]*vec3[1]-vec2[1]*vec3[0]);

  *mydet = dettmp;

  return;
}
/****************************************************************************************/

