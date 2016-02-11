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

/***********************************************************************************************/
void precond_hx_curl_additive (REAL *r,
                      REAL *z,
                      void *data)
{
    /**
     * \fn void precond_hx_curl_additive (REAL *r, REAL *z, void *data)
     *
     * \brief HX preconditioner for H(curl): additive version
     *
     * \param r     Pointer to the vector needs preconditioning
     * \param z     Pointer to preconditioned vector
     * \param data  Pointer to precondition data
     *
     * \author Xiaozhe Hu
     * \date   02/10/2016
     */
    
    HX_curl_data *hxcurldata=(HX_curl_data *)data;
    INT n = hxcurldata->A->row;
    SHORT smooth_iter = hxcurldata->smooth_iter;
    
    // make sure z is initialzied by zeros
    array_set(n, z, 0.0);
    
    // local variable
    dvector zz;
    zz.row = n; zz.val = z;
    dvector rr;
    rr.row = n; rr.val = r;
    
    SHORT maxit, i;
    
    // smoothing
    smoother_dcsr_sgs(&zz, hxcurldata->A, &rr, smooth_iter);

    // solve vector Laplacian
    AMG_param *amgparam_vgrad = hxcurldata->amgparam_vgrad;
    AMG_data *mgl_vgrad = hxcurldata->mgl_vgrad;
    maxit = amgparam_vgrad->maxit;
    
    mgl_vgrad->b.row = hxcurldata->A_vgrad->row;
    dcsr_mxv(hxcurldata->Pt_curl, r, mgl_vgrad->b.val);
    mgl_vgrad->x.row=hxcurldata->A_vgrad->row;
    dvec_set(hxcurldata->A_vgrad->row, &mgl_vgrad->x, 0.0);
    
    for (i=0;i<maxit;++i) mgcycle(mgl_vgrad, amgparam_vgrad);
    
    dcsr_aAxpy(1.0, hxcurldata->P_curl, mgl_vgrad->x.val, z);
    
    // solve scalar Laplacian
    AMG_param *amgparam_grad = hxcurldata->amgparam_grad;
    AMG_data *mgl_grad = hxcurldata->mgl_grad;
    maxit = amgparam_grad->maxit;
    
    mgl_grad->b.row = hxcurldata->A_grad->row;
    dcsr_mxv(hxcurldata->Gradt, r, mgl_grad->b.val);
    mgl_grad->x.row=hxcurldata->A_grad->row;
    dvec_set(hxcurldata->A_grad->row, &mgl_grad->x, 0.0);
    
    for (i=0;i<maxit;++i) mgcycle(mgl_grad, amgparam_grad);
    
    dcsr_aAxpy(1.0, hxcurldata->Grad, mgl_grad->x.val, z);
    
}

/***********************************************************************************************/
void precond_hx_curl_multiplicative (REAL *r,
                               REAL *z,
                               void *data)
{
    /**
     * \fn void precond_hx_curl_multiplicative (REAL *r, REAL *z, void *data)
     *
     * \brief HX preconditioner for H(curl): multiplicative version
     *
     * \param r     Pointer to the vector needs preconditioning
     * \param z     Pointer to preconditioned vector
     * \param data  Pointer to precondition data
     *
     * \author Xiaozhe Hu
     * \date   02/10/2016
     */
    
    HX_curl_data *hxcurldata=(HX_curl_data *)data;
    INT n = hxcurldata->A->row;
    SHORT smooth_iter = hxcurldata->smooth_iter;
    
    // backup r
    array_cp(n, r, hxcurldata->backup_r);
    
    // make sure z is initialzied by zeros
    array_set(n, z, 0.0);
    
    // local variable
    dvector zz;
    zz.row = n; zz.val = z;
    dvector rr;
    rr.row = n; rr.val = r;
    
    SHORT maxit, i;
    
    // smoothing
    smoother_dcsr_sgs(&zz, hxcurldata->A, &rr, smooth_iter);
    
    // update r
    dcsr_aAxpy(-1.0, hxcurldata->A, zz.val, rr.val);
    
    // solve vector Laplacian
    AMG_param *amgparam_vgrad = hxcurldata->amgparam_vgrad;
    AMG_data *mgl_vgrad = hxcurldata->mgl_vgrad;
    maxit = amgparam_vgrad->maxit;
    
    mgl_vgrad->b.row = hxcurldata->A_vgrad->row;
    dcsr_mxv(hxcurldata->Pt_curl, r, mgl_vgrad->b.val);
    mgl_vgrad->x.row=hxcurldata->A_vgrad->row;
    dvec_set(hxcurldata->A_vgrad->row, &mgl_vgrad->x, 0.0);
    
    for (i=0;i<maxit;++i) mgcycle(mgl_vgrad, amgparam_vgrad);
    
    dcsr_aAxpy(1.0, hxcurldata->P_curl, mgl_vgrad->x.val, z);
    
    // update r
    array_cp(n, hxcurldata->backup_r, r);
    dcsr_aAxpy(-1.0, hxcurldata->A, zz.val, rr.val);
    
    // solve scalar Laplacian
    AMG_param *amgparam_grad = hxcurldata->amgparam_grad;
    AMG_data *mgl_grad = hxcurldata->mgl_grad;
    maxit = amgparam_grad->maxit;
    
    mgl_grad->b.row = hxcurldata->A_grad->row;
    dcsr_mxv(hxcurldata->Gradt, r, mgl_grad->b.val);
    mgl_grad->x.row=hxcurldata->A_grad->row;
    dvec_set(hxcurldata->A_grad->row, &mgl_grad->x, 0.0);
    
    for (i=0;i<maxit;++i) mgcycle(mgl_grad, amgparam_grad);
    
    dcsr_aAxpy(1.0, hxcurldata->Grad, mgl_grad->x.val, z);
    
    // store r
    array_cp(n, hxcurldata->backup_r, r);
    
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
