/*! \file src/solver/precond.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 10/06/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note  Done cleanup for releasing -- Xiaozhe Hu 03/12/2017
 *
 *  \todo  Add block preconditoners for general size -- Xiaozhe Hu
 *  \todo  Add inexact block solvers for Stokes/NS equations -- Xiaozhe Hu
 *
 */

#include "hazmath.h"
#include "itsolver_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/***********************************************************************************************/
/**
 * \fn void precond_diag(REAL *r, REAL *z, void *data)
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
void precond_diag (REAL *r, 
                        REAL *z, 
                        void *data)
{   
    dvector *diag=(dvector *)data;
    REAL *diagptr=diag->val;
    INT i, m=diag->row;    
    
    memcpy(z,r,m*sizeof(REAL));
    for (i=0;i<m;++i) {
        if (ABS(diag->val[i])>SMALLREAL) z[i]/=diagptr[i];
    }    
}

/***********************************************************************************************/
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
void precond_amg(REAL *r,
                 REAL *z,
                 void *data)
{  
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
void precond_amli(REAL *r,
                  REAL *z,
                  void *data)
{  
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
void precond_nl_amli(REAL *r,
                     REAL *z,
                     void *data)
{    
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
void precond_hx_curl_additive(REAL *r,
                              REAL *z,
                              void *data)
{   
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
/**
 * \fn void precond_hx_curl_multiplicative(REAL *r, REAL *z, void *data)
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
void precond_hx_curl_multiplicative(REAL *r,
                                    REAL *z,
                                    void *data)
{   
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

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_2(REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved exactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   10/14/2016
 */
void precond_block_diag_2(REAL *r,
                          REAL *z,
                          void *data)
{
#if WITH_SUITESPARSE
  precond_block_data *precdata=(precond_block_data *)data;
  dCSRmat *A_diag = precdata->A_diag;
  dvector *tempr = &(precdata->r);
  
  const INT N0 = A_diag[0].row;
  const INT N1 = A_diag[1].row;
  const INT N = N0 + N1;
  
  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);
  
  // prepare
  //#if  WITH_UMFPACK
  void **LU_diag = precdata->LU_diag;
  dvector r0, r1, z0, z1;
  
  r0.row = N0; z0.row = N0;
  r1.row = N1; z1.row = N1;
  
  r0.val = r; r1.val = &(r[N0]);
  z0.val = z; z1.val = &(z[N0]);
  //#endif
  
  // Preconditioning A00 block
  //#if  WITH_UMFPACK
  /* use UMFPACK direct solver */
  umfpack_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);
  //#endif
  
  // Preconditioning A11 block
  //#if  WITH_UMFPACK
  /* use UMFPACK direct solver */
  umfpack_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);
  //#endif
  
  // restore r
  array_cp(N, tempr->val, r);
  
#endif
}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_2 (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (3x3 block matrix, each diagonal
 *        block is solved exactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   10/14/2016
 */
void precond_block_lower_2(REAL *r,
                           REAL *z,
                           void *data)
{ 
#if WITH_SUITESPARSE
  
  precond_block_data *precdata=(precond_block_data *)data;
  block_dCSRmat *A = precdata->Abcsr;
  dCSRmat *A_diag = precdata->A_diag;
  void **LU_diag = precdata->LU_diag;
  
  dvector *tempr = &(precdata->r);
  
  const INT N0 = A_diag[0].row;
  const INT N1 = A_diag[1].row;
  const INT N = N0 + N1;
  
  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);
  
  // prepare
  dvector r0, r1, z0, z1;
  
  r0.row = N0; z0.row = N0;
  r1.row = N1; z1.row = N1;
  
  r0.val = r; r1.val = &(r[N0]);
  z0.val = z; z1.val = &(z[N0]);
  
  // Preconditioning A00 block
  //#if  WITH_UMFPACK
  /* use UMFPACK direct solver */
  umfpack_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);
  //#endif
  
  // r1 = r1 - A2*z0
  dcsr_aAxpy(-1.0, A->blocks[2], z0.val, r1.val);
  
  // Preconditioning A11 block
  //#if  WITH_UMFPACK
  /* use UMFPACK direct solver */
  umfpack_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);
  //#endif
  
  // restore r
  array_cp(N, tempr->val, r);
  
#endif
  
}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_2 (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (2x2 block matrix, each diagonal
 *        block is solved exactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   10/14/2016
 */
void precond_block_upper_2(REAL *r,
                           REAL *z,
                           void *data)
{ 
#if WITH_SUITESPARSE
  
  precond_block_data *precdata=(precond_block_data *)data;
  block_dCSRmat *A = precdata->Abcsr;
  dCSRmat *A_diag = precdata->A_diag;
  void **LU_diag = precdata->LU_diag;
  
  dvector *tempr = &(precdata->r);
  
  const INT N0 = A_diag[0].row;
  const INT N1 = A_diag[1].row;
  const INT N = N0 + N1;
  
  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);
  
  // prepare
  dvector r0, r1, z0, z1;
  
  r0.row = N0; z0.row = N0;
  r1.row = N1; z1.row = N1;
  
  r0.val = r; r1.val = &(r[N0]);
  z0.val = z; z1.val = &(z[N0]);
  
  // Preconditioning A11 block
  //#if  WITH_UMFPACK
  /* use UMFPACK direct solver */
  umfpack_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);
  //#endif
  
  // r0 = r0 - A1*z1
  dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
  
  // Preconditioning A00 block
  //#if  WITH_UMFPACK
  /* use UMFPACK direct solver */
  umfpack_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);
  //#endif
  
  // restore r
  array_cp(N, tempr->val, r);
  
#endif
  
}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_3 (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (3x3 block matrix, each diagonal block
 *        is solved exactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   02/24/2014
 */
void precond_block_diag_3(REAL *r,
                          REAL *z,
                          void *data)
{   
#if WITH_SUITESPARSE
    precond_block_data *precdata=(precond_block_data *)data;
    dCSRmat *A_diag = precdata->A_diag;
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N = N0 + N1 + N2;
    
    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);
    
    // prepare
    void **LU_diag = precdata->LU_diag;
    dvector r0, r1, r2, z0, z1, z2;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    
    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    
    // Preconditioning A00 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);
    
    // Preconditioning A11 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);
    
    // Preconditioning A22 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[2], &r2, &z2, LU_diag[2], 0);
    
    // restore r
    array_cp(N, tempr->val, r);
   
#endif
}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_3 (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (3x3 block matrix, each diagonal
 *        block is solved exactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   02/24/2016
 */
void precond_block_lower_3(REAL *r,
                           REAL *z,
                           void *data)
{ 
#if WITH_SUITESPARSE
    
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_diag = precdata->A_diag;
    void **LU_diag = precdata->LU_diag;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N = N0 + N1 + N2;
    
    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, r2, z0, z1, z2;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    
    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    
    // Preconditioning A00 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);
    
    // r1 = r1 - A3*z0
    if (A->blocks[3] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[3], z0.val, r1.val);

    // Preconditioning A11 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);
    
    // r2 = r2 - A6*z0 - A7*z1
    if (A->blocks[6] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[6], z0.val, r2.val);
    if (A->blocks[7] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[7], z1.val, r2.val);
    
    // Preconditioning A22 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[2], &r2, &z2, LU_diag[2], 0);
    
    // restore r
    array_cp(N, tempr->val, r);
    
#endif
    
}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_3 (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (3x3 block matrix, each diagonal
 *        block is solved exactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   02/24/2016
 */
void precond_block_upper_3(REAL *r,
                           REAL *z,
                           void *data)
{ 
#if WITH_SUITESPARSE

    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_diag = precdata->A_diag;
    void **LU_diag = precdata->LU_diag;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N = N0 + N1 + N2;
    
    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, r2, z0, z1, z2;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    
    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    
    // Preconditioning A22 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[2], &r2, &z2, LU_diag[2], 0);
    
    // r1 = r1 - A5*z2
    if (A->blocks[5] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[5], z2.val, r1.val);
    
    // Preconditioning A11 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);
    
    // r0 = r0 - A1*z1 - A2*z2
    if (A->blocks[1] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
    if (A->blocks[2] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[2], z2.val, r0.val);
    
    // Preconditioning A00 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);
    
    // restore r
    array_cp(N, tempr->val, r);
    
#endif
    
}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_4 (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (4x4 block matrix, each diagonal block
 *        is solved exactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   01/20/2017
 */
void precond_block_diag_4(REAL *r,
                          REAL *z,
                          void *data)
{
#if WITH_SUITESPARSE
    precond_block_data *precdata=(precond_block_data *)data;
    dCSRmat *A_diag = precdata->A_diag;
    dvector *tempr = &(precdata->r);

    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N3 = A_diag[3].row;
    const INT N = N0 + N1 + N2 + N3;

    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);

    // prepare
    void **LU_diag = precdata->LU_diag;
    dvector r0, r1, r2, r3, z0, z1, z2, z3;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    r3.row = N3; z3.row = N3;

    r0.val = r;
    r1.val = &(r[N0]);
    r2.val = &(r[N0+N1]);
    r3.val = &(r[N0+N1+N2]);
    z0.val = z;
    z1.val = &(z[N0]);
    z2.val = &(z[N0+N1]);
    z3.val = &(z[N0+N1+N2]);

    // Preconditioning A00 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);

    // Preconditioning A11 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);

    // Preconditioning A22 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[2], &r2, &z2, LU_diag[2], 0);

    // Preconditioning A33 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[3], &r3, &z3, LU_diag[3], 0);

    // restore r
    array_cp(N, tempr->val, r);

#endif
}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_4 (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (4x4 block matrix, each diagonal
 *        block is solved exactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   01/28/2017
 *
 * A[0]  A[1]  A[2]  A[3]
 * A[4]  A[5]  A[6]  A[7]
 * A[8]  A[9]  A[10] A[11]
 * A[12] A[13] A[14] A[15]
 */
void precond_block_lower_4(REAL *r,
                           REAL *z,
                           void *data)
{
#if WITH_SUITESPARSE

    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_diag = precdata->A_diag;
    void **LU_diag = precdata->LU_diag;

    dvector *tempr = &(precdata->r);

    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N3 = A_diag[3].row;
    const INT N = N0 + N1 + N2 + N3;

    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);

    // prepare
    dvector r0, r1, r2, r3, z0, z1, z2, z3;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    r3.row = N3; z3.row = N3;

    r0.val = r;
    r1.val = &(r[N0]);
    r2.val = &(r[N0+N1]);
    r3.val = &(r[N0+N1+N2]);
    z0.val = z;
    z1.val = &(z[N0]);
    z2.val = &(z[N0+N1]);
    z3.val = &(z[N0+N1+N2]);

    // Preconditioning A00 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);

    // r1 = r1 - A4*z0
    if (A->blocks[4] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[4], z0.val, r1.val);

    // Preconditioning A11 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);

    // r2 = r2 - A8*z0 - A9*z1
    if (A->blocks[8] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[8], z0.val, r2.val);
    if (A->blocks[9] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[9], z1.val, r2.val);

    // Preconditioning A22 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[2], &r2, &z2, LU_diag[2], 0);

    // r3 = r3 - A12*z0 - A13*z1 - A14*z2
    if (A->blocks[12] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[12], z0.val, r3.val);
    if (A->blocks[13] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[13], z1.val, r3.val);
    if (A->blocks[14] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[14], z2.val, r3.val);

    // Preconditioning A33 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[3], &r3, &z3, LU_diag[3], 0);

    // restore r
    array_cp(N, tempr->val, r);

#endif

}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_4 (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (4x4 block matrix, each diagonal
 *        block is solved exactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   01/28/2017
 *
 * A[0]  A[1]  A[2]  A[3]
 * A[4]  A[5]  A[6]  A[7]
 * A[8]  A[9]  A[10] A[11]
 * A[12] A[13] A[14] A[15]
 */
void precond_block_upper_4(REAL *r,
                           REAL *z,
                           void *data)
{
#if WITH_SUITESPARSE

    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_diag = precdata->A_diag;
    void **LU_diag = precdata->LU_diag;

    dvector *tempr = &(precdata->r);

    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N3 = A_diag[3].row;
    const INT N = N0 + N1 + N2 + N3;

    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);

    // prepare
    dvector r0, r1, r2, r3, z0, z1, z2, z3;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    r3.row = N3; z3.row = N3;

    r0.val = r;
    r1.val = &(r[N0]);
    r2.val = &(r[N0+N1]);
    r3.val = &(r[N0+N1+N2]);
    z0.val = z;
    z1.val = &(z[N0]);
    z2.val = &(z[N0+N1]);
    z3.val = &(z[N0+N1+N2]);

    // Preconditioning A33 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[3], &r3, &z3, LU_diag[3], 0);

    // r2 = r2 - A11*z3
    if (A->blocks[11] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[11], z3.val, r2.val);

    // Preconditioning A22 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[2], &r2, &z2, LU_diag[2], 0);

    // r1 = r1 - A6*z2 - A7*z3
    if (A->blocks[6] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[6], z2.val, r1.val);
    if (A->blocks[7] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[7], z3.val, r1.val);

    // Preconditioning A11 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);

    // r0 = r0 - A1*z1 - A2*z2 - A3*z3
    if (A->blocks[1] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
    if (A->blocks[2] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[2], z2.val, r0.val);
    if (A->blocks[3] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[3], z3.val, r0.val);

    // Preconditioning A00 block
    /* use UMFPACK direct solver */
    umfpack_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);

    // restore r
    array_cp(N, tempr->val, r);

#endif

}


/*************** Special Preconditioners for Mixed Darcy Flow *********************************/

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_mixed_darcy (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   10/14/2016
 */
void precond_block_diag_mixed_darcy(REAL *r,
                                    REAL *z,
                                    void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);
  
  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;
  dvector *el_vol = precdata->el_vol;
  
  INT i;

  const INT N0 = A->blocks[0]->row;
  const INT N1 = A->blocks[2]->row;
  const INT N = N0 + N1;
  
  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);
  
  // prepare
  dvector r0, r1, z0, z1;
  
  r0.row = N0; z0.row = N0;
  r1.row = N1; z1.row = N1;
  
  r0.val = r; r1.val = &(r[N0]);
  z0.val = z; z1.val = &(z[N0]);
  //#endif
  
  // Preconditioning A00 block (flux)
  mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
  mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x,0.0);
  
  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
  array_cp(N0, mgl[0]->x.val, z0.val);
  
  // Preconditioning A11 block
  memcpy(z1.val,r1.val,N1*sizeof(REAL));
  for (i=0;i<N1;++i) {
    z1.val[i]/=el_vol->val[i];
  }
  
  // restore r
  array_cp(N, tempr->val, r);
  
}


/***********************************************************************************************/
/**
 * \fn void precond_block_lower_mixed_darcy (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   10/14/2016
 */
void precond_block_lower_mixed_darcy(REAL *r,
                                    REAL *z,
                                    void *data)
{  
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);
  
  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;
  dvector *el_vol = precdata->el_vol;
  
  INT i;
  
  const INT N0 = A->blocks[0]->row;
  const INT N1 = A->blocks[2]->row;
  const INT N = N0 + N1;
  
  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);
  
  // prepare
  dvector r0, r1, z0, z1;
  
  r0.row = N0; z0.row = N0;
  r1.row = N1; z1.row = N1;
  
  r0.val = r; r1.val = &(r[N0]);
  z0.val = z; z1.val = &(z[N0]);
  //#endif
  
  // Preconditioning A00 block (flux)
  mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
  mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x,0.0);
  
  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
  array_cp(N0, mgl[0]->x.val, z0.val);
  
  // r1 = r1 - A2*z0
  dcsr_aAxpy(-1.0, A->blocks[2], z0.val, r1.val);
  
  // Preconditioning A11 block
  memcpy(z1.val,r1.val,N1*sizeof(REAL));
  for (i=0;i<N1;++i) {
    z1.val[i]/=el_vol->val[i];
  }
  
  // restore r
  array_cp(N, tempr->val, r);
  
}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_mixed_darcy (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   10/14/2016
 */
void precond_block_upper_mixed_darcy(REAL *r,
                                     REAL *z,
                                     void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);
  
  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;
  dvector *el_vol = precdata->el_vol;
  
  INT i;
  
  const INT N0 = A->blocks[0]->row;
  const INT N1 = A->blocks[2]->row;
  const INT N = N0 + N1;
  
  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);
  
  // prepare
  dvector r0, r1, z0, z1;
  
  r0.row = N0; z0.row = N0;
  r1.row = N1; z1.row = N1;
  
  r0.val = r; r1.val = &(r[N0]);
  z0.val = z; z1.val = &(z[N0]);
  //#endif
  
  // Preconditioning A11 block
  memcpy(z1.val,r1.val,N1*sizeof(REAL));
  for (i=0;i<N1;++i) {
    z1.val[i]/=el_vol->val[i];
  }
  
  // r0 = r0 - A1*z1
  dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
  
  // Preconditioning A00 block (flux)
  mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
  mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x,0.0);
  
  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
  array_cp(N0, mgl[0]->x.val, z0.val);
  
  // restore r
  array_cp(N, tempr->val, r);
  
}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_mixed_darcy (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly by Krylov methods)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   10/14/2016
 */
void precond_block_diag_mixed_darcy_krylov(REAL *r,
                                           REAL *z,
                                           void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);
  
  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;
  dvector *el_vol = precdata->el_vol;
  
  INT i;
  
  const INT N0 = A->blocks[0]->row;
  const INT N1 = A->blocks[2]->row;
  const INT N = N0 + N1;
  
  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);
  
  // prepare
  dvector r0, r1, z0, z1;
  
  r0.row = N0; z0.row = N0;
  r1.row = N1; z1.row = N1;
  
  r0.val = r; r1.val = &(r[N0]);
  z0.val = z; z1.val = &(z[N0]);
  //#endif
  
  // Preconditioning A00 block (flux)
  precond_data pcdata_p;
  param_amg_to_prec(&pcdata_p,amgparam);
  pcdata_p.max_levels = mgl[0][0].num_levels;
  pcdata_p.mgl_data = mgl[0];
  
  precond pc_p; pc_p.data = &pcdata_p;
  pc_p.fct = precond_amg;
  
  dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc_p, 1e-3, 100, 100, 1, 1);

  
  // Preconditioning A11 block
  memcpy(z1.val,r1.val,N1*sizeof(REAL));
  for (i=0;i<N1;++i) {
    z1.val[i]/=el_vol->val[i];
  }
  
  // restore r
  array_cp(N, tempr->val, r);
  
}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_mixed_darcy (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly by Krylov methods)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   10/14/2016
 */
void precond_block_lower_mixed_darcy_krylov(REAL *r,
                                            REAL *z,
                                            void *data)
{ 
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);
  
  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;
  dvector *el_vol = precdata->el_vol;
  
  INT i;
  
  const INT N0 = A->blocks[0]->row;
  const INT N1 = A->blocks[2]->row;
  const INT N = N0 + N1;
  
  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);
  
  // prepare
  dvector r0, r1, z0, z1;
  
  r0.row = N0; z0.row = N0;
  r1.row = N1; z1.row = N1;
  
  r0.val = r; r1.val = &(r[N0]);
  z0.val = z; z1.val = &(z[N0]);
  //#endif
  
  // Preconditioning A00 block (flux)
  precond_data pcdata_p;
  param_amg_to_prec(&pcdata_p,amgparam);
  pcdata_p.max_levels = mgl[0][0].num_levels;
  pcdata_p.mgl_data = mgl[0];
  
  precond pc_p; pc_p.data = &pcdata_p;
  pc_p.fct = precond_amg;
  
  dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc_p, 1e-3, 100, 100, 1, 1);
  
  // r1 = r1 - A2*z0
  dcsr_aAxpy(-1.0, A->blocks[2], z0.val, r1.val);
  
  // Preconditioning A11 block
  memcpy(z1.val,r1.val,N1*sizeof(REAL));
  for (i=0;i<N1;++i) {
    z1.val[i]/=el_vol->val[i];
  }
  
  // restore r
  array_cp(N, tempr->val, r);
  
}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_mixed_darcy (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly by Krylov methods)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   10/14/2016
 */
void precond_block_upper_mixed_darcy_krylov(REAL *r,
                                            REAL *z,
                                            void *data)
{ 
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);
  
  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;
  dvector *el_vol = precdata->el_vol;
  
  INT i;
  
  const INT N0 = A->blocks[0]->row;
  const INT N1 = A->blocks[2]->row;
  const INT N = N0 + N1;
  
  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);
  
  // prepare
  dvector r0, r1, z0, z1;
  
  r0.row = N0; z0.row = N0;
  r1.row = N1; z1.row = N1;
  
  r0.val = r; r1.val = &(r[N0]);
  z0.val = z; z1.val = &(z[N0]);
  //#endif
  
  // Preconditioning A11 block
  memcpy(z1.val,r1.val,N1*sizeof(REAL));
  for (i=0;i<N1;++i) {
    z1.val[i]/=el_vol->val[i];
  }
  
  // r0 = r0 - A1*z1
  dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
  
  // Preconditioning A00 block (flux)
  precond_data pcdata_p;
  param_amg_to_prec(&pcdata_p,amgparam);
  pcdata_p.max_levels = mgl[0][0].num_levels;
  pcdata_p.mgl_data = mgl[0];
  
  precond pc_p; pc_p.data = &pcdata_p;
  pc_p.fct = precond_amg;
  
  dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc_p, 1e-3, 100, 100, 1, 1);

  
  // restore r
  array_cp(N, tempr->val, r);
  
}

/*************** Special Preconditioners for Biot 2-phase formulation **************************/
/***********************************************************************************************/
/**
 * \fn void precond_block_diag_biot_2phase (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   01/14/2017
 */
void precond_block_diag_biot_2phase(REAL *r,
                                    REAL *z,
                                    void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;

  INT i;

  const INT N0 = A->blocks[0]->row;
  const INT N1 = A->blocks[3]->row;
  const INT N = N0 + N1;

  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);

  // prepare
  dvector r0, r1, z0, z1;

  r0.row = N0; z0.row = N0;
  r1.row = N1; z1.row = N1;

  r0.val = r; r1.val = &(r[N0]);
  z0.val = z; z1.val = &(z[N0]);
  //#endif

  // Preconditioning A00 block (displacement)
  mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
  mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
  array_cp(N0, mgl[0]->x.val, z0.val);

  // Preconditioning A11 block (pressure)
  mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
  mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
  array_cp(N1, mgl[1]->x.val, z1.val);

  // restore r
  array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_biot_2phase (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   01/14/2017
 */
void precond_block_lower_biot_2phase(REAL *r,
                                     REAL *z,
                                     void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;

  INT i;

  const INT N0 = A->blocks[0]->row;
  const INT N1 = A->blocks[3]->row;
  const INT N = N0 + N1;

  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);

  // prepare
  dvector r0, r1, z0, z1;

  r0.row = N0; z0.row = N0;
  r1.row = N1; z1.row = N1;

  r0.val = r; r1.val = &(r[N0]);
  z0.val = z; z1.val = &(z[N0]);
  //#endif

  // Preconditioning A00 block (displacement)
  mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
  mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
  array_cp(N0, mgl[0]->x.val, z0.val);

  // r1 = r1 - A2*z0
  dcsr_aAxpy(-1.0, A->blocks[2], z0.val, r1.val);

  // Preconditioning A11 block
  mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
  mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
  array_cp(N1, mgl[1]->x.val, z1.val);

  // restore r
  array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_biot_2phase (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   01/14/2017
 */
void precond_block_upper_biot_2phase(REAL *r,
                                     REAL *z,
                                     void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;

  INT i;

  const INT N0 = A->blocks[0]->row;
  const INT N1 = A->blocks[3]->row;
  const INT N = N0 + N1;

  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);

  // prepare
  dvector r0, r1, z0, z1;

  r0.row = N0; z0.row = N0;
  r1.row = N1; z1.row = N1;

  r0.val = r; r1.val = &(r[N0]);
  z0.val = z; z1.val = &(z[N0]);
  //#endif

  // Preconditioning A11 block (pressure)
  mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
  mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
  array_cp(N1, mgl[1]->x.val, z1.val);

  // r0 = r0 - A1*z1
  dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);

  // Preconditioning A00 block (flux)
  mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
  mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
  array_cp(N0, mgl[0]->x.val, z0.val);

  // restore r
  array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_biot_2phase_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly by Krylov methods)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   01/14/2017
 */
void precond_block_diag_biot_2phase_krylov(REAL *r,
                                           REAL *z,
                                           void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;

  const INT N0 = A->blocks[0]->row;
  const INT N1 = A->blocks[3]->row;
  const INT N = N0 + N1;

  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);

  // prepare
  dvector r0, r1, z0, z1;

  r0.row = N0; z0.row = N0;
  r1.row = N1; z1.row = N1;

  r0.val = r; r1.val = &(r[N0]);
  z0.val = z; z1.val = &(z[N0]);
  //#endif

  // Preconditioning A00 block (displacement)
  precond_data pcdata_u;
  param_amg_to_prec(&pcdata_u,amgparam);
  pcdata_u.max_levels = mgl[0][0].num_levels;
  pcdata_u.mgl_data = mgl[0];

  precond pc_u; pc_u.data = &pcdata_u;
  pc_u.fct = precond_amg;

  dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc_u, 1e-3, 100, 100, 1, 1);


  // Preconditioning A11 block (pressure)
  precond_data pcdata_p;
  param_amg_to_prec(&pcdata_p,amgparam);
  pcdata_p.max_levels = mgl[1][0].num_levels;
  pcdata_p.mgl_data = mgl[1];

  precond pc_p; pc_p.data = &pcdata_p;
  pc_p.fct = precond_amg;

  dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc_p, 1e-3, 100, 100, 1, 1);

  // restore r
  array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_biot_2phase_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly by Krylov methods)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   01/14/2017
 */
void precond_block_lower_biot_2phase_krylov(REAL *r,
                                            REAL *z,
                                            void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;

  const INT N0 = A->blocks[0]->row;
  const INT N1 = A->blocks[3]->row;
  const INT N = N0 + N1;

  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);

  // prepare
  dvector r0, r1, z0, z1;

  r0.row = N0; z0.row = N0;
  r1.row = N1; z1.row = N1;

  r0.val = r; r1.val = &(r[N0]);
  z0.val = z; z1.val = &(z[N0]);
  //#endif

  // Preconditioning A00 block (displacement)
  precond_data pcdata_u;
  param_amg_to_prec(&pcdata_u,amgparam);
  pcdata_u.max_levels = mgl[0][0].num_levels;
  pcdata_u.mgl_data = mgl[0];

  precond pc_u; pc_u.data = &pcdata_u;
  pc_u.fct = precond_amg;

  dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc_u, 1e-3, 100, 100, 1, 1);

  // r1 = r1 - A2*z0
  dcsr_aAxpy(-1.0, A->blocks[2], z0.val, r1.val);

  // Preconditioning A11 block
  precond_data pcdata_p;
  param_amg_to_prec(&pcdata_p,amgparam);
  pcdata_p.max_levels = mgl[1][0].num_levels;
  pcdata_p.mgl_data = mgl[1];

  precond pc_p; pc_p.data = &pcdata_p;
  pc_p.fct = precond_amg;

  dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc_p, 1e-3, 100, 100, 1, 1);

  // restore r
  array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_biot_2phase_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly by Krylov methods)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   01/14/2017
 */
void precond_block_upper_biot_2phase_krylov(REAL *r,
                                            REAL *z,
                                            void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;

  const INT N0 = A->blocks[0]->row;
  const INT N1 = A->blocks[3]->row;
  const INT N = N0 + N1;

  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);

  // prepare
  dvector r0, r1, z0, z1;

  r0.row = N0; z0.row = N0;
  r1.row = N1; z1.row = N1;

  r0.val = r; r1.val = &(r[N0]);
  z0.val = z; z1.val = &(z[N0]);
  //#endif

  // Preconditioning A11 block (pressure)
  precond_data pcdata_p;
  param_amg_to_prec(&pcdata_p,amgparam);
  pcdata_p.max_levels = mgl[1][0].num_levels;
  pcdata_p.mgl_data = mgl[1];

  precond pc_p; pc_p.data = &pcdata_p;
  pc_p.fct = precond_amg;

  dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc_p, 1e-3, 100, 100, 1, 1);

  // r0 = r0 - A1*z1
  dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);

  // Preconditioning A00 block (displacement)
  precond_data pcdata_u;
  param_amg_to_prec(&pcdata_u,amgparam);
  pcdata_u.max_levels = mgl[0][0].num_levels;
  pcdata_u.mgl_data = mgl[0];

  precond pc_u; pc_u.data = &pcdata_u;
  pc_u.fct = precond_amg;

  dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc_u, 1e-3, 100, 100, 1, 1);

  // restore r
  array_cp(N, tempr->val, r);

}


/*************** Special Preconditioners for Maxwell equation **********************************/

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_maxwell (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (maxwell equation)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/12/2016
 */
void precond_block_diag_maxwell(REAL *r,
                                REAL *z,
                                void *data)
{   
    precond_block_data *precdata=(precond_block_data *)data;
    //block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_diag = precdata->A_diag;
    //void **LU_diag = precdata->LU_diag;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    HX_curl_data **hxcurldata = precdata->hxcurldata;
    
    INT i;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N = N0 + N1 + N2;
    
    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, r2, z0, z1, z2;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    
    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    
    // Preconditioning A00 block
    /* use AMG solver */
    /*
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x,0.0);
     
     for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
     array_cp(N0, mgl[0]->x.val, z0.val);
     */
    memcpy(z0.val,r0.val,N0*sizeof(REAL));
    for (i=0;i<N0;++i) {
        if (ABS(precdata->diag[0]->val[i])>SMALLREAL) z0.val[i]/=precdata->diag[0]->val[i];
    }
    
    
    // Preconditioning A11 block
    /* use HX preconditioner */
    precond_hx_curl_multiplicative(r1.val, z1.val, hxcurldata[1]);
    
    // Preconditioning A22 block
    /* use AMG solver */
    mgl[2]->b.row=N2; array_cp(N2, r2.val, mgl[2]->b.val); // residual is an input
    mgl[2]->x.row=N2; dvec_set(N2, &mgl[2]->x,0.0);
    
    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[2], amgparam);
    array_cp(N2, mgl[2]->x.val, z2.val);
    
    // restore r
    array_cp(N, tempr->val, r);
    
}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_maxwell (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (maxwell equation)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/12/2016
 */
void precond_block_lower_maxwell(REAL *r,
                                 REAL *z,
                                 void *data)
{  
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_diag = precdata->A_diag;
    //void **LU_diag = precdata->LU_diag;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    HX_curl_data **hxcurldata = precdata->hxcurldata;
    
    INT i;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N = N0 + N1 + N2;
    
    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, r2, z0, z1, z2;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    
    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    
    // Preconditioning A00 block
    /* use AMG solver */
    /*
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x,0.0);
    
    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
    array_cp(N0, mgl[0]->x.val, z0.val);
     */
    memcpy(z0.val,r0.val,N0*sizeof(REAL));
    for (i=0;i<N0;++i) {
        if (ABS(precdata->diag[0]->val[i])>SMALLREAL) z0.val[i]/=precdata->diag[0]->val[i];
    }
    
    // r1 = r1 - A3*z0
    dcsr_aAxpy(-1.0, A->blocks[3], z0.val, r1.val);
    
    // Preconditioning A11 block
    precond_hx_curl_multiplicative(r1.val, z1.val, hxcurldata[1]);
    
    // r2 = r2 - A6*z0 - A7*z1
    dcsr_aAxpy(-1.0, A->blocks[6], z0.val, r2.val);
    dcsr_aAxpy(-1.0, A->blocks[7], z1.val, r2.val);
    
    // Preconditioning A22 block
    /* use AMG solver */
    mgl[2]->b.row=N2; array_cp(N2, r2.val, mgl[2]->b.val); // residual is an input
    mgl[2]->x.row=N2; dvec_set(N2, &mgl[2]->x,0.0);
    
    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[2], amgparam);
    array_cp(N2, mgl[2]->x.val, z2.val);
    
    // restore r
    array_cp(N, tempr->val, r);
    
}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_maxwell (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (maxwell equation)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   02/25/2016
 */
void precond_block_upper_maxwell(REAL *r,
                                 REAL *z,
                                 void *data)
{   
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_diag = precdata->A_diag;
    //void **LU_diag = precdata->LU_diag;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    HX_curl_data **hxcurldata = precdata->hxcurldata;
    
    INT i;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N = N0 + N1 + N2;
    
    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, r2, z0, z1, z2;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    
    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    
    // Preconditioning A22 block
    /* use AMG solver */
    mgl[2]->b.row=N2; array_cp(N2, r2.val, mgl[2]->b.val); // residual is an input
    mgl[2]->x.row=N2; dvec_set(N2, &mgl[2]->x,0.0);
    
    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[2], amgparam);
    array_cp(N2, mgl[2]->x.val, z2.val);
    
    // r1 = r1 - A5*z2
    dcsr_aAxpy(-1.0, A->blocks[5], z2.val, r1.val);
    
    // Preconditioning A11 block
    precond_hx_curl_multiplicative(r1.val, z1.val, hxcurldata[1]);
    
    // r0 = r0 - A1*z1 - A2*z2
    dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
    dcsr_aAxpy(-1.0, A->blocks[2], z2.val, r0.val);
    
    // Preconditioning A00 block
    /* use AMG solver */
    /*
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x,0.0);
    
    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
    array_cp(N0, mgl[0]->x.val, z0.val);
     */
    memcpy(z0.val,r0.val,N0*sizeof(REAL));
    for (i=0;i<N0;++i) {
        if (ABS(precdata->diag[0]->val[i])>SMALLREAL) z0.val[i]/=precdata->diag[0]->val[i];
    }
    
    // restore r
    array_cp(N, tempr->val, r);
    
}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_maxwell_krylov (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (maxwell equation)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/12/2016
 */
void precond_block_diag_maxwell_krylov(REAL *r,
                                       REAL *z,
                                       void *data)
{   
    precond_block_data *precdata=(precond_block_data *)data;
    // block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_diag = precdata->A_diag;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    HX_curl_data **hxcurldata = precdata->hxcurldata;
    
    //INT i;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N = N0 + N1 + N2;
    
    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, r2, z0, z1, z2;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    
    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    
    // Preconditioning A00 block
    /* use AMG+Krylov solver */
    /*
    precond_data pcdata_B;
    param_amg_to_prec(&pcdata_B,amgparam);
    pcdata_B.max_levels = mgl[0][0].num_levels;
    pcdata_B.mgl_data = mgl[0];
    
    precond pc_B; pc_B.data = &pcdata_B;
    pc_B.fct = precond_amg;
    
    dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc_B, 1e-2, 100, 100, 1, 1);
     */
    precond pc_B; pc_B.data = precdata->diag[0];
    pc_B.fct = precond_diag;
    
    dcsr_pvfgmres(&A_diag[0], &r0, &z0, &pc_B, 1e-2, 100, 100, 1, 1);
    
    // Preconditioning A11 block
    /* use HX preconditioner+Krylov solver */
    precond pc_E; pc_E.data = hxcurldata[1];
    pc_E.fct = precond_hx_curl_multiplicative;
    dcsr_pvfgmres(&A_diag[1], &r1, &z1, &pc_E, 1e-2, 100, 100, 1, 1);

    // Preconditioning A22 block
    /* use AMG+Krylov solver */
    precond_data pcdata_p;
    param_amg_to_prec(&pcdata_p,amgparam);
    pcdata_p.max_levels = mgl[2][0].num_levels;
    pcdata_p.mgl_data = mgl[2];
    
    precond pc_p; pc_p.data = &pcdata_p;
    pc_p.fct = precond_amg;
    
    dcsr_pvfgmres(&mgl[2][0].A, &r2, &z2, &pc_p, 1e-2, 100, 100, 1, 1);
    
    // restore r
    array_cp(N, tempr->val, r);
    
}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_maxwell_krylov (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (maxwell equation)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   02/26/2016
 */
void precond_block_lower_maxwell_krylov(REAL *r,
                                        REAL *z,
                                        void *data)
{   
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_diag = precdata->A_diag;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    HX_curl_data **hxcurldata = precdata->hxcurldata;
    
    //    INT i;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N = N0 + N1 + N2;
    
    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, r2, z0, z1, z2;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    
    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    
    // Preconditioning A00 block
    /* use AMG+Krylov solver */
    /*
    precond_data pcdata_B;
    param_amg_to_prec(&pcdata_B,amgparam);
    pcdata_B.max_levels = mgl[0][0].num_levels;
    pcdata_B.mgl_data = mgl[0];
    
    precond pc_B; pc_B.data = &pcdata_B;
    pc_B.fct = precond_amg;
    
    dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc_B, 1e-2, 100, 100, 1, 1);
     */
    precond pc_B; pc_B.data = precdata->diag[0];
    pc_B.fct = precond_diag;
    
    dcsr_pvfgmres(&A_diag[0], &r0, &z0, &pc_B, 1e-2, 100, 100, 1, 1);
    
    // r1 = r1 - A3*z0
    dcsr_aAxpy(-1.0, A->blocks[3], z0.val, r1.val);
    
    // Preconditioning A11 block
    /* use HX preconditioner+Krylov solver */
    precond pc_E; pc_E.data = hxcurldata[1];
    pc_E.fct = precond_hx_curl_multiplicative;
    dcsr_pvfgmres(&A_diag[1], &r1, &z1, &pc_E, 1e-2, 100, 100, 1, 1);
    
    // r2 = r2 - A6*z0 - A7*z1
    dcsr_aAxpy(-1.0, A->blocks[6], z0.val, r2.val);
    dcsr_aAxpy(-1.0, A->blocks[7], z1.val, r2.val);

    // Preconditioning A22 block
    /* use AMG+Krylov solver */
    precond_data pcdata_p;
    param_amg_to_prec(&pcdata_p,amgparam);
    pcdata_p.max_levels = mgl[2][0].num_levels;
    pcdata_p.mgl_data = mgl[2];
    
    precond pc_p; pc_p.data = &pcdata_p;
    pc_p.fct = precond_amg;
    
    dcsr_pvfgmres(&mgl[2][0].A, &r2, &z2, &pc_p, 1e-2, 100, 100, 1, 1);
    
    // restore r
    array_cp(N, tempr->val, r);
    
}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_maxwell_krylov (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (maxwell equation)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   02/26/2016
 */
void precond_block_upper_maxwell_krylov(REAL *r,
                                        REAL *z,
                                        void *data)
{    
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_diag = precdata->A_diag;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    HX_curl_data **hxcurldata = precdata->hxcurldata;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N = N0 + N1 + N2;
    
    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, r2, z0, z1, z2;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    
    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    
    // Preconditioning A22 block
    /* use AMG+Krylov solver */
    precond_data pcdata_p;
    param_amg_to_prec(&pcdata_p,amgparam);
    pcdata_p.max_levels = mgl[2][0].num_levels;
    pcdata_p.mgl_data = mgl[2];
    
    precond pc_p; pc_p.data = &pcdata_p;
    pc_p.fct = precond_amg;
    
    dcsr_pvfgmres(&mgl[2][0].A, &r2, &z2, &pc_p, 1e-2, 100, 100, 1, 1);
    
    // r1 = r1 - A5*z2
    dcsr_aAxpy(-1.0, A->blocks[5], z2.val, r1.val);
    
    // Preconditioning A11 block
    precond pc_E; pc_E.data = hxcurldata[1];
    pc_E.fct = precond_hx_curl_multiplicative;
    dcsr_pvfgmres(&A_diag[1], &r1, &z1, &pc_E, 1e-2, 100, 100, 1, 1);
    
    // r0 = r0 - A1*z1 - A2*z2
    dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
    dcsr_aAxpy(-1.0, A->blocks[2], z2.val, r0.val);
    
    // Preconditioning A00 block
    /* use AMG+Krylov solver */
    /*
    precond_data pcdata_B;
    param_amg_to_prec(&pcdata_B,amgparam);
    pcdata_B.max_levels = mgl[0][0].num_levels;
    pcdata_B.mgl_data = mgl[0];
     
     precond pc_B; pc_B.data = &pcdata_B;
     pc_B.fct = precond_amg;
    
    dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc_B, 1e-2, 100, 100, 1, 1);
     */
    
    precond pc_B; pc_B.data = precdata->diag[0];
    pc_B.fct = precond_diag;
    
    dcsr_pvfgmres(&A_diag[0], &r0, &z0, &pc_B, 1e-2, 100, 100, 1, 1);
        
    // restore r
    array_cp(N, tempr->val, r);
    
}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_diag_maxwell (REAL *r, REAL *z, void *data)
 * \brief block diagonal/upper triangular preconditioning (maxwell equation)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/12/2016
 */
void precond_block_lower_diag_maxwell(REAL *r,
                                      REAL *z,
                                      void *data)
{   
    precond_block_data *precdata=(precond_block_data *)data;
    dCSRmat *A_diag = precdata->A_diag;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    HX_curl_data **hxcurldata = precdata->hxcurldata;
    
    //dCSRmat *G = precdata->G;
    //dCSRmat *K = precdata->K;
    dCSRmat *Gt = precdata->Gt;
    dCSRmat *Kt = precdata->Kt;
    
    INT i;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N = N0 + N1 + N2;
    
    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, r2, z0, z1, z2;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    
    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    
    // lower blocks
    
    // r1 = r1 + K^T * r0
    dcsr_aAxpy(1.0, Kt, r0.val, r1.val);
    
    // r2 = r2 + G^t * r1
    dcsr_aAxpy(1.0, Gt, r1.val, r2.val);
    
    // diagonal blocks
    
    // Preconditioning A00 block
    /* use AMG solver */
    /*
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x,0.0);
    
    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
    array_cp(N0, mgl[0]->x.val, z0.val);
     */
    memcpy(z0.val,r0.val,N0*sizeof(REAL));
    for (i=0;i<N0;++i) {
        if (ABS(precdata->diag[0]->val[i])>SMALLREAL) z0.val[i]/=precdata->diag[0]->val[i];
    }
    
    // Preconditioning A11 block
    /* use HX preconditioner */
    precond_hx_curl_multiplicative(r1.val, z1.val, hxcurldata[1]);
    
    // Preconditioning A22 block
    /* use AMG solver */
    mgl[2]->b.row=N2; array_cp(N2, r2.val, mgl[2]->b.val); // residual is an input
    mgl[2]->x.row=N2; dvec_set(N2, &mgl[2]->x,0.0);
    
    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[2], amgparam);
    array_cp(N2, mgl[2]->x.val, z2.val);
    
    // restore r
    array_cp(N, tempr->val, r);
    
}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_upper_maxwell_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal/upper triangular preconditioning (maxwell equation)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/12/2016
 */
void precond_block_diag_upper_maxwell(REAL *r,
                                      REAL *z,
                                      void *data)
{   
    precond_block_data *precdata=(precond_block_data *)data;
    dCSRmat *A_diag = precdata->A_diag;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    HX_curl_data **hxcurldata = precdata->hxcurldata;
    
    dCSRmat *G = precdata->G;
    dCSRmat *K = precdata->K;
    //dCSRmat *Gt = precdata->Gt;
    //dCSRmat *Kt = precdata->Kt;
    
    INT i;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N = N0 + N1 + N2;
    
    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, r2, z0, z1, z2;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    
    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    
    // diagonal blocks
    
    // Preconditioning A00 block
    /* use AMG solver */
    /*
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x,0.0);
    
    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
    array_cp(N0, mgl[0]->x.val, z0.val);
     */
    
    memcpy(z0.val,r0.val,N0*sizeof(REAL));
    for (i=0;i<N0;++i) {
        if (ABS(precdata->diag[0]->val[i])>SMALLREAL) z0.val[i]/=precdata->diag[0]->val[i];
    }
    
    // Preconditioning A11 block
    /* use HX preconditioner */
    precond_hx_curl_multiplicative(r1.val, z1.val, hxcurldata[1]);
    
    // Preconditioning A22 block
    /* use AMG solver */
    mgl[2]->b.row=N2; array_cp(N2, r2.val, mgl[2]->b.val); // residual is an input
    mgl[2]->x.row=N2; dvec_set(N2, &mgl[2]->x,0.0);
    
    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[2], amgparam);
    array_cp(N2, mgl[2]->x.val, z2.val);
    
    // upper blocks
    
    // z1 = z1 - G*z2
    dcsr_aAxpy(-1.0, G, z2.val, z1.val);
    
    // z0 = z0 - K*z1
    dcsr_aAxpy(-1.0, K, z1.val, z0.val);
    
    // restore r
    array_cp(N, tempr->val, r);
    
}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_diag_upper_maxwell_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal/upper triangular preconditioning (maxwell equation)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/12/2016
 */
void precond_block_lower_diag_upper_maxwell(REAL *r,
                                            REAL *z,
                                            void *data)
{   
    precond_block_data *precdata=(precond_block_data *)data;
    dCSRmat *A_diag = precdata->A_diag;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    HX_curl_data **hxcurldata = precdata->hxcurldata;
    
    dCSRmat *G = precdata->G;
    dCSRmat *K = precdata->K;
    dCSRmat *Gt = precdata->Gt;
    dCSRmat *Kt = precdata->Kt;
    
    INT i;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N = N0 + N1 + N2;
    
    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, r2, z0, z1, z2;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    
    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    
    // lower blocks
    
    // r1 = r1 + K^T * r0
    dcsr_aAxpy(1.0, Kt, r0.val, r1.val);
    
    // r2 = r2 + G^t * r1
    dcsr_aAxpy(1.0, Gt, r1.val, r2.val);
    
    // diagonal blocks
    
    // Preconditioning A00 block
    /* use AMG solver */
    /*
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x,0.0);
    
    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
    array_cp(N0, mgl[0]->x.val, z0.val);
     */
    
    memcpy(z0.val,r0.val,N0*sizeof(REAL));
    for (i=0;i<N0;++i) {
        if (ABS(precdata->diag[0]->val[i])>SMALLREAL) z0.val[i]/=precdata->diag[0]->val[i];
    }
    
    // Preconditioning A11 block
    /* use HX preconditioner */
    precond_hx_curl_multiplicative(r1.val, z1.val, hxcurldata[1]);
    
    // Preconditioning A22 block
    /* use AMG solver */
    mgl[2]->b.row=N2; array_cp(N2, r2.val, mgl[2]->b.val); // residual is an input
    mgl[2]->x.row=N2; dvec_set(N2, &mgl[2]->x,0.0);
    
    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[2], amgparam);
    array_cp(N2, mgl[2]->x.val, z2.val);
    
    // upper blocks
    
    // z1 = z1 - G*z2
    dcsr_aAxpy(-1.0, G, z2.val, z1.val);
    
    // z0 = z0 - K*z1
    dcsr_aAxpy(-1.0, K, z1.val, z0.val);
    
    // restore r
    array_cp(N, tempr->val, r);
    
}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_diag_maxwell_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal/upper triangular preconditioning (maxwell equation)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/12/2016
 */
void precond_block_lower_diag_maxwell_krylov(REAL *r,
                                             REAL *z,
                                             void *data)
{   
    precond_block_data *precdata=(precond_block_data *)data;
    dCSRmat *A_diag = precdata->A_diag;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    HX_curl_data **hxcurldata = precdata->hxcurldata;
    
    //dCSRmat *G = precdata->G;
    //dCSRmat *K = precdata->K;
    dCSRmat *Gt = precdata->Gt;
    dCSRmat *Kt = precdata->Kt;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N = N0 + N1 + N2;
    
    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, r2, z0, z1, z2;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    
    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    
    // lower blocks
    
    // r1 = r1 + K^T * r0
    dcsr_aAxpy(1.0, Kt, r0.val, r1.val);
    
    // r2 = r2 + G^t * r1
    dcsr_aAxpy(1.0, Gt, r1.val, r2.val);
    
    // diagonal blocks
    
    // Preconditioning A22 block
    /* use AMG+Krylov solver */
    //printf("solve p\n");
    
    precond_data pcdata_p;
    param_amg_to_prec(&pcdata_p,amgparam);
    pcdata_p.max_levels = mgl[2][0].num_levels;
    pcdata_p.mgl_data = mgl[2];
    
    precond pc_p; pc_p.data = &pcdata_p;
    pc_p.fct = precond_amg;

    dcsr_pvfgmres(&mgl[2][0].A, &r2, &z2, &pc_p, 1e-2, 100, 100, 1, 1);
    
    // Preconditioning A11 block
    /* use HX preconditioner+Krylov solver */
    //printf("solve E\n");
    
    precond pc_E; pc_E.data = hxcurldata[1];
    pc_E.fct = precond_hx_curl_multiplicative;
    dcsr_pvfgmres(&A_diag[1], &r1, &z1, &pc_E, 1e-2, 100, 100, 1, 1);
    
    // Preconditioning A00 block
    /* use AMG+Krylov solver */
    //printf("solve B\n");
    
    /*
    precond_data pcdata_B;
    param_amg_to_prec(&pcdata_B,amgparam);
    pcdata_B.max_levels = mgl[0][0].num_levels;
    pcdata_B.mgl_data = mgl[0];
    
    precond pc_B; pc_B.data = &pcdata_B;
    pc_B.fct = precond_amg;
     */
    
    precond pc_B; pc_B.data = precdata->diag[0];
    pc_B.fct = precond_diag;
    
    dcsr_pvfgmres(&A_diag[0], &r0, &z0, &pc_B, 1e-2, 100, 100, 1, 1);
    
    // restore r
    array_cp(N, tempr->val, r);
    
}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_upper_maxwell_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal/upper triangular preconditioning (maxwell equation)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/12/2016
 */
void precond_block_diag_upper_maxwell_krylov(REAL *r,
                                             REAL *z,
                                             void *data)
{   
    precond_block_data *precdata=(precond_block_data *)data;
    dCSRmat *A_diag = precdata->A_diag;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    HX_curl_data **hxcurldata = precdata->hxcurldata;
    
    dCSRmat *G = precdata->G;
    dCSRmat *K = precdata->K;
    //dCSRmat *Gt = precdata->Gt;
    //dCSRmat *Kt = precdata->Kt;

    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N = N0 + N1 + N2;
    
    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, r2, z0, z1, z2;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    
    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    
    // diagonal blocks
    
    // Preconditioning A22 block
    /* use AMG+Krylov solver */
    precond_data pcdata_p;
    param_amg_to_prec(&pcdata_p,amgparam);
    pcdata_p.max_levels = mgl[2][0].num_levels;
    pcdata_p.mgl_data = mgl[2];
    
    precond pc_p; pc_p.data = &pcdata_p;
    pc_p.fct = precond_amg;
    
    dcsr_pvfgmres(&mgl[2][0].A, &r2, &z2, &pc_p, 1e-2, 100, 100, 1, 1);
    
    
    // Preconditioning A11 block
    /* use HX preconditioner+Krylov solver */
    precond pc_E; pc_E.data = hxcurldata[1];
    pc_E.fct = precond_hx_curl_multiplicative;
    dcsr_pvfgmres(&A_diag[1], &r1, &z1, &pc_E, 1e-2, 100, 100, 1, 1);
    
    
    // Preconditioning A00 block
    /* use AMG+Krylov solver */
    /*
    precond_data pcdata_B;
    param_amg_to_prec(&pcdata_B,amgparam);
    pcdata_B.max_levels = mgl[0][0].num_levels;
    pcdata_B.mgl_data = mgl[0];
    
    precond pc_B; pc_B.data = &pcdata_B;
    pc_B.fct = precond_amg;
     
     dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc_B, 1e-2, 100, 100, 1, 1);
     */
    
    precond pc_B; pc_B.data = precdata->diag[0];
    pc_B.fct = precond_diag;
    
    dcsr_pvfgmres(&A_diag[0], &r0, &z0, &pc_B, 1e-2, 100, 100, 1, 1);

    
    // upper blocks
    
    // z1 = z1 - G*z2
    dcsr_aAxpy(-1.0, G, z2.val, z1.val);
    
    // z0 = z0 - K*z1
    dcsr_aAxpy(-1.0, K, z1.val, z0.val);
    
    // restore r
    array_cp(N, tempr->val, r);
    
}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_diag_upper_maxwell_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal/upper triangular preconditioning (maxwell equation)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/12/2016
 */
void precond_block_lower_diag_upper_maxwell_krylov(REAL *r,
                                                   REAL *z,
                                                   void *data)
{   
    precond_block_data *precdata=(precond_block_data *)data;
    dCSRmat *A_diag = precdata->A_diag;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    HX_curl_data **hxcurldata = precdata->hxcurldata;
    
    dCSRmat *G = precdata->G;
    dCSRmat *K = precdata->K;
    dCSRmat *Gt = precdata->Gt;
    dCSRmat *Kt = precdata->Kt;
    
    dvector *tempr = &(precdata->r);
    
    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N = N0 + N1 + N2;
    
    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);
    
    // prepare
    dvector r0, r1, r2, z0, z1, z2;
    
    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    
    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    
    // lower blocks
    
    // r1 = r1 + K^T * r0
    dcsr_aAxpy(1.0, Kt, r0.val, r1.val);
    
    // r2 = r2 + G^t * r1
    dcsr_aAxpy(1.0, Gt, r1.val, r2.val);
    
    // diagonal blocks
    
    // Preconditioning A22 block
    /* use AMG+Krylov solver */
    precond_data pcdata_p;
    param_amg_to_prec(&pcdata_p,amgparam);
    pcdata_p.max_levels = mgl[2][0].num_levels;
    pcdata_p.mgl_data = mgl[2];
    
    precond pc_p; pc_p.data = &pcdata_p;
    pc_p.fct = precond_amg;
    
    dcsr_pvfgmres(&mgl[2][0].A, &r2, &z2, &pc_p, 1e-2, 100, 100, 1, 1);
    
    // Preconditioning A11 block
    /* use HX preconditioner+Krylov solver */
    precond pc_E; pc_E.data = hxcurldata[1];
    pc_E.fct = precond_hx_curl_multiplicative;
    dcsr_pvfgmres(&A_diag[1], &r1, &z1, &pc_E, 1e-2, 100, 100, 1, 1);
    
    // Preconditioning A00 block
    /* use AMG+Krylov solver */
    /*
    precond_data pcdata_B;
    param_amg_to_prec(&pcdata_B,amgparam);
    pcdata_B.max_levels = mgl[0][0].num_levels;
    pcdata_B.mgl_data = mgl[0];
    
    precond pc_B; pc_B.data = &pcdata_B;
    pc_B.fct = precond_amg;
    
    dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc_B, 1e-2, 100, 100, 1, 1);
     */
    
    precond pc_B; pc_B.data = precdata->diag[0];
    pc_B.fct = precond_diag;
    
    dcsr_pvfgmres(&A_diag[0], &r0, &z0, &pc_B, 1e-2, 100, 100, 1, 1);
    
    // upper blocks
    
    // z1 = z1 - G*z2
    dcsr_aAxpy(-1.0, G, z2.val, z1.val);
    
    // z0 = z0 - K*z1
    dcsr_aAxpy(-1.0, K, z1.val, z0.val);
    
    // restore r
    array_cp(N, tempr->val, r);
    
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
