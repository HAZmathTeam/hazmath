/*
 *  precond.c
 *
 *  Created by James Adler and Xiaozhe Hu on 10/06/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

#include "hazmat.h"

#include "itsolver_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/***********************************************************************************************/
void precond_diag (REAL *r, 
                        REAL *z, 
                        void *data)
{
    /**
     * \fn void precond_diag (REAL *r, REAL *z, void *data)
     *
     * \brief Diagonal preconditioner z=inv(D)*r
     *
     * \param r     Pointer to the vector needs preconditioning
     * \param z     Pointer to preconditioned vector
     * \param data  Pointer to precondition data
     *
     * \author Xiaozhe Hu
     * \date   04/06/2010
     */
    
    dvector *diag=(dvector *)data;
    REAL *diagptr=diag->val;
    INT i, m=diag->row;    
    
    memcpy(z,r,m*sizeof(REAL));
    for (i=0;i<m;++i) {
        if (ABS(diag->val[i])>SMALLREAL) z[i]/=diagptr[i];
    }    
}

/***********************************************************************************************/
void precond_amg (REAL *r,
                       REAL *z,
                       void *data)
{
    /**
     * \fn void precond_amg (REAL *r, REAL *z, void *data)
     *
     * \brief AMG preconditioner
     *
     * \param r     Pointer to the vector needs preconditioning
     * \param z     Pointer to preconditioned vector
     * \param data  Pointer to precondition data
     *
     * \author Xiaozhe Hu
     * \date   12/30/2015
     */
    
    precond_data *pcdata=(precond_data *)data;
    const INT m=pcdata->mgl_data[0].A.row;
    const INT maxit=pcdata->maxit;
    INT i;
    
    AMG_param amgparam; param_amg_init(&amgparam);
    param_prec_to_amg(&amgparam,pcdata);
    
    AMG_data *mgl = pcdata->mgl_data;
    mgl->b.row=m; array_cp(m,r,mgl->b.val); // residual is an input
    mgl->x.row=m; dvec_set(m,&mgl->x,0.0);
    
    for (i=0;i<maxit;++i) mgcycle(mgl,&amgparam);
    
    array_cp(m,mgl->x.val,z);
}


/***********************************************************************************************/
void precond_amli (REAL *r,
                        REAL *z,
                        void *data)
{
    /**
     * \fn void precond_amli(REAL *r, REAL *z, void *data)
     *
     * \brief AMLI AMG preconditioner
     *
     * \param r     Pointer to the vector needs preconditioning
     * \param z     Pointer to preconditioned vector
     * \param data  Pointer to precondition data
     *
     * \author Xiaozhe Hu
     * \date   01/23/2011
     */
    
    precond_data *pcdata=(precond_data *)data;
    const INT m=pcdata->mgl_data[0].A.row;
    const INT maxit=pcdata->maxit;
    INT i;
    
    AMG_param amgparam; param_amg_init(&amgparam);
    param_prec_to_amg(&amgparam,pcdata);
    
    AMG_data *mgl = pcdata->mgl_data;
    mgl->b.row=m; array_cp(m,r,mgl->b.val); // residual is an input
    mgl->x.row=m; dvec_set(m,&mgl->x,0.0);
    
    for (i=0;i<maxit;++i) amli(mgl,&amgparam,0);
    
    array_cp(m,mgl->x.val,z);
}

/***********************************************************************************************/
void precond_nl_amli (REAL *r,
                           REAL *z,
                           void *data)
{
    /**
     * \fn void precond_nl_amli(REAL *r, REAL *z, void *data)
     *
     * \brief Nonlinear AMLI AMG preconditioner
     *
     * \param r     Pointer to the vector needs preconditioning
     * \param z     Pointer to preconditioned vector
     * \param data  Pointer to precondition data
     *
     * \author Xiaozhe Hu
     * \date   04/25/2011
     */
    
    precond_data *pcdata=(precond_data *)data;
    const INT m=pcdata->mgl_data[0].A.row;
    const INT maxit=pcdata->maxit;
    const SHORT num_levels = pcdata->max_levels;
    INT i;
    
    AMG_param amgparam; param_amg_init(&amgparam);
    param_prec_to_amg(&amgparam,pcdata);
    
    AMG_data *mgl = pcdata->mgl_data;
    mgl->b.row=m; array_cp(m,r,mgl->b.val); // residual is an input
    mgl->x.row=m; dvec_set(m,&mgl->x,0.0);
    
    for (i=0;i<maxit;++i) nl_amli(mgl, &amgparam, 0, num_levels);
    
    array_cp(m,mgl->x.val,z);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
