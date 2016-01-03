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
void precond_ilu (REAL *r,
                       REAL *z,
                       void *data)
{
    /**
     * \fn void precond_ilu (REAL *r, REAL *z, void *data)
     *
     * \brief ILU preconditioner
     *
     * \param r     Pointer to the vector needs preconditioning
     * \param z     Pointer to preconditioned vector
     * \param data  Pointer to precondition data
     *
     * \author Shiquan Zhang
     * \date   04/06/2010
     */
    
    ILU_data *iludata=(ILU_data *)data;
    const INT m=iludata->row, mm1=m-1, memneed=2*m;
    REAL *zz, *zr;
    
    if (iludata->nwork<memneed) goto MEMERR; // check this outside this subroutine!!
    
    zz = iludata->work;
    zr = iludata->work+m;
    array_cp(m, r, zr);
    
    {
        INT i, j, jj, begin_row, end_row, mm2=m-2;
        INT *ijlu=iludata->ijlu;
        REAL *lu=iludata->luval;
        
        // forward sweep: solve unit lower matrix equation L*zz=zr
        zz[0]=zr[0];
        
        for (i=1;i<=mm1;++i) {
            begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
            for (j=begin_row;j<=end_row;++j) {
                jj=ijlu[j];
                if (jj<i) zr[i]-=lu[j]*zz[jj];
                else break;
            }
            zz[i]=zr[i];
        }
        
        // backward sweep: solve upper matrix equation U*z=zz
        z[mm1]=zz[mm1]*lu[mm1];
        for (i=mm2;i>=0;i--) {
            begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
            for (j=end_row;j>=begin_row;j--) {
                jj=ijlu[j];
                if (jj>i) zz[i]-=lu[j]*z[jj];
                else break;
            }
            z[i]=zz[i]*lu[i];
        }
    }
    
    return;
    
MEMERR:
    printf("### ERROR: Need %d memory, only %d available!\n", memneed, iludata->nwork);
    exit(ERROR_ALLOC_MEM);
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
