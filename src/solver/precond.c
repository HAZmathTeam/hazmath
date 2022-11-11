/*! \file src/solver/precond.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 10/06/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note  Done cleanup for releasing -- Xiaozhe Hu 03/12/2017 & 08/28/2021
 *
 *  \todo  Add block triangular preconditoners for general size -- Xiaozhe Hu
 *
 */

#include "hazmath.h"
/*! \file itsolver_util.inl
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 5/13/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note  Done cleanup for releasing -- Xiaozhe Hu 03/12/2017
 *
 */

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

//! Warning for residual false convergence
#define ITS_FACONV  printf("### HAZMATH WARNING: False convergence!\n")

//! Warning for solution close to zero
#define ITS_ZEROSOL printf("### HAZMATH WARNING: Iteration stopped due to the solution is almost zero! %s : %lld\n", __FUNCTION__,   (long long )__LINE__)

//! Warning for iteration restarted
#define ITS_RESTART printf("### HAZMATH WARNING: Iteration restarted due to stagnation! %s : %lld\n", __FUNCTION__,   (long long )__LINE__)

//! Warning for stagged iteration
#define ITS_STAGGED printf("### HAZMATH WARNING: Iteration stopped due to staggnation! %s : %lld\n", __FUNCTION__,   (long long )__LINE__)

//! Warning for tolerance practically close to zero
#define ITS_ZEROTOL printf("### HAZMATH WARNING: The tolerence might be too small! %s : %lld\n", __FUNCTION__,   (long long )__LINE__)

//! Warning for divided by zero
#define ITS_DIVZERO printf("### HAZMATH WARNING: Divided by zero! %s : %lld\n", __FUNCTION__,   (long long )__LINE__)

//! Warning for actual relative residual
#define ITS_REALRES(relres) printf("### HAZMATH WARNING: The actual relative residual = %e!\n",(relres))

//! Warning for computed relative residual
#define ITS_COMPRES(relres) printf("### HAZMATH WARNING: The computed relative residual = %e!\n",(relres))

//! Warning for too small sp
#define ITS_SMALLSP printf("### HAZMATH WARNING: sp is too small! %s : %lld\n", __FUNCTION__,   (long long )__LINE__)

//! Warning for restore previous iteration
#define ITS_RESTORE(iter) printf("### HAZMATH WARNING: Restore iteration %lld!\n",  (long long )(iter));

//! Output relative difference and residual
#define ITS_DIFFRES(reldiff,relres) printf("||u-u'|| = %e and the comp. rel. res. = %e.\n",(reldiff),(relres));

//! Output L2 norm of some variable
#define ITS_PUTNORM(name,value) printf("L2 norm of %s = %e.\n",(name),(value));

/***********************************************************************************************/
/**
 * \fn inline static void ITS_CHECK (const INT MaxIt, const REAL tol)
 * \brief Safeguard checks to prevent unexpected error for iterative solvers
 *
 * \param MaxIt   Maximal number of iterations
 * \param tol     Tolerance for convergence check
 *
 */
inline static void ITS_CHECK (const INT MaxIt, const REAL tol)
{
    if ( tol < SMALLREAL ) {
        printf("### HAZMATH WARNING: Convergence tolerance for iterative solver is too small!\n");
    }
    if ( MaxIt <= 0 ) {
        printf("### HAZMATH WARNING: Max number of iterations should be a POSITIVE integer!\n");
    }
}

/***********************************************************************************************/
/**
 * \fn inline static void ITS_FINAL (const INT iter, const INT MaxIt, const REAL relres)
 * \brief Print out final status of an iterative method
 *
 * \param iter    Number of iterations
 * \param MaxIt   Maximal number of iterations
 * \param relres  Relative residual
 *
 */
inline static void ITS_FINAL (const INT iter, const INT MaxIt, const REAL relres)
{
    if ( iter > MaxIt ) {
        printf("### HAZMATH WARNING: Max iter %lld reached with rel. resid. %e.\n",   (long long )MaxIt, relres);
    }
    else if ( iter >= 0 ) {
        printf("Number of iterations = %lld with relative residual %e.\n",   (long long )iter, relres);
    }
}
/***********************************************************************************************/
/**
 * \fn inline static void WARN_STATUS(const char *function_name,const char *call_to, const INT status)
 * \brief Print out a warning
 *
 * \param function_name    the name of calling function
 * \param call_to          the name of the function returning "status"
 * \param status           the status that needs to be reported.
 *
 */
inline static void WARN_STATUS(const char *function_name,const char *call_to, const INT status)
{
  fprintf(stderr,"\n\n%%%% ****WARNING in %s: status=%lld after exiting %s (WHILE SUCCESS .EQ. %lld)\n\n", \
	  function_name,  (long long )status,call_to,  (long long )SUCCESS);
}
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
 * \fn void precond_famg (REAL *r, REAL *z, void *data)
 *
 * \brief Fractional AMG preconditioner
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Ana Budisa
 * \date   2020-05-11
 */
void precond_famg(REAL *r,
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

    for (i=0;i<maxit;++i) fmgcycle(mgl,&amgparam);

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
 * \fn void precond_famli(REAL *r, REAL *z, void *data)
 *
 * \brief AMLI fractional AMG preconditioner
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Ana Budisa
 * \date   2020-05-11
 */
void precond_famli(REAL *r,
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

    for (i=0;i<maxit;++i) famli(mgl,&amgparam,0);

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
 * \fn void precond_amg_add (REAL *r, REAL *z, void *data)
 *
 * \brief AMG preconditioner (additive AMG)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   05/28/2020
 */
void precond_amg_add(REAL *r,
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

    for (i=0;i<maxit;++i) mgcycle_add(mgl,&amgparam);
    //mgcycle_add_update(mgl,&amgparam);

    array_cp(m,mgl->x.val,z);
}

/***********************************************************************************************/
/**
 * \fn void precond_famg_add (REAL *r, REAL *z, void *data)
 *
 * \brief fractional AMG preconditioner (additive AMG)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu, Ana Budisa
 * \date   2020-06-03
 *
 * Note: for 0<s<1, use famg(s); for -1<s<0, use famg((1+s)/2) A famg((1+s)/2)
 */
void precond_famg_add(REAL *r,
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

    if (pcdata->fpwr < 0 && pcdata->fpwr >= -1) {
        // temp work - TODO: change to array type
        dvector x1 = dvec_create(m);

        // famg preconditioner
        amgparam.fpwr = 0.5 * (1 + pcdata->fpwr);
        // solve famg((1+s)/2)
        mgl->b.row=m; array_cp(m,r,mgl->b.val); // residual is an input
        mgl->x.row=m; dvec_set(m,&mgl->x,0.0);
        for (i=0;i<maxit;++i) fmgcycle_add_update(mgl,&amgparam);

        // x1 = A x
        dcsr_mxv(&mgl->A, mgl->x.val, x1.val);

        // solve famg((1+s)/2)
        mgl->b.row=m; array_cp(m,x1.val,mgl->b.val); // x1 is an input
        mgl->x.row=m; dvec_set(m,&mgl->x,0.0);
        for (i=0;i<maxit;++i) fmgcycle_add_update(mgl,&amgparam);
    }
    else if (pcdata->fpwr >= 0 && pcdata->fpwr <= 1) {
        mgl->b.row=m; array_cp(m,r,mgl->b.val); // residual is an input
        mgl->x.row=m; dvec_set(m,&mgl->x,0.0);
        // call iterative solver
        for (i=0;i<maxit;++i) fmgcycle_add_update(mgl,&amgparam);
        //fmgcycle_add(mgl,&amgparam);
    }
    else {
        printf("\n !!! Fractionality s = %.2f is not within the limits -1 <= s <= 1 !!! \n", pcdata->fpwr);
        printf("Assuming s = 1.0 \n");
        amgparam.fpwr = 1.0;
        mgl->b.row=m; array_cp(m,r,mgl->b.val); // residual is an input
        mgl->x.row=m; dvec_set(m,&mgl->x,0.0);
        // call iterative solver
        for (i=0;i<maxit;++i) fmgcycle_add_update(mgl,&amgparam);
        //fmgcycle_add(mgl,&amgparam);
        return;
    }

    array_cp(m,mgl->x.val,z);
}


/***********************************************************************************************/
/**
 * \fn void precond_famg_add2 (REAL *r, REAL *z, void *data)
 *
 * \brief fractional AMG preconditioner (additive AMG)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Ana Budisa
 * \date   2020-07-06
 *
 * Note: for 0<s<1, use famg(A^s); for -1<s<0, use div famg(A_div^(1+s)) grad
 */
void precond_famg_add2(REAL *r,
                       REAL *z,
                       void *data)
{
    precond_data *pcdata=(precond_data *)data; // 0 - amg, 1 - grad, 2 - grad^T
    const INT m=pcdata[0].mgl_data[0].A.row;
    const INT maxit=pcdata[0].maxit;
    INT i;
    REAL fpwr = pcdata[0].fpwr;

    AMG_param amgparam; param_amg_init(&amgparam);
    param_prec_to_amg(&amgparam, &pcdata[0]);

    AMG_data *mgl = pcdata[0].mgl_data;

    if (fpwr < 0 && fpwr >= -1) {
        // pcdata[1].A is a pointer to Grad; pcdata[2].A is a pointer to Div;
        // mg_rhs = Grad r // (Grad r) is the input for amg
        dcsr_mxv(pcdata[1].A, r, mgl->b.val);

        // famg preconditioner for A_div^(1+s)
        amgparam.fpwr = fpwr + 1;
        // solve famg(A_div^(1+s))
        mgl->x.row=m; dvec_set(m,&mgl->x,0.0);
        for (i=0;i<maxit;++i) fmgcycle_add_update(mgl,&amgparam);
        // z = Div x
        dcsr_mxv(pcdata[2].A, mgl->x.val, z);
    }
    else if (fpwr >= 0 && fpwr <= 1) {
        mgl->b.row=m; array_cp(m,r,mgl->b.val); // residual is an input
        mgl->x.row=m; dvec_set(m,&mgl->x,0.0);
        // call iterative solver
        for (i=0;i<maxit;++i) fmgcycle_add_update(mgl,&amgparam);
        //fmgcycle_add(mgl,&amgparam);
        array_cp(m,mgl->x.val,z);
    }
    else {
        printf("\n !!! Fractionality s = %.2f is not within the limits -1 <= s <= 1 !!! \n", pcdata->fpwr);
        printf("Assuming s = 1.0 \n");
        amgparam.fpwr = 1.0;
        mgl->b.row=m; array_cp(m,r,mgl->b.val); // residual is an input
        mgl->x.row=m; dvec_set(m,&mgl->x,0.0);
        // call iterative solver
        for (i=0;i<maxit;++i) fmgcycle_add_update(mgl,&amgparam);
        //fmgcycle_add(mgl,&amgparam);
        array_cp(m,mgl->x.val,z);
    }

}


/***********************************************************************************************/
/**
 * \fn void precond_sum_famg_add (REAL *r, REAL *z, void *data)
 *
 * \brief fractional AMG preconditioner (additive AMG) for operators of type
 *        alpha * D^s + beta * D^(1+s)
 *        where D is a laplacian and s \in (-1,0)
 *        This consists of
 *        FAMG(s/2) * AMG(alpha * lump(M)^-1 + beta * lump(M)^-1 A lump(M)^-1) * FAMG(s/2)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Ana Budisa
 * \date   2020-06-18
 *
 */
void precond_sum_famg_add(REAL *r,
                          REAL *z,
                          void *data)
{
    // data[0] - FAMG; data[1] - AMG;
    precond_data *pcdata=(precond_data *)data;

    // first, FAMG(s/2)
    INT m1 = pcdata[0].mgl_data[0].A.row;
    REAL *x1 = (REAL *)calloc(m1, sizeof(REAL));
    precond_famg_add(r, x1, &pcdata[0]); // x1 is the result (preconditioned r)

    // second AMG(alpha * lump(M)^-1 + beta * lump(M)^-1 A lump(M)^-1)
    INT m2 = pcdata[1].mgl_data[0].A.row;
    REAL *x2 = (REAL *)calloc(m2, sizeof(REAL));
    precond_amg_add(x1, x2, &pcdata[1]); // x2 is the result (preconditioned x1)

    // third FAMG(s/2)
    precond_famg_add(x2, z, &pcdata[0]); // z is the result (preconditioned x2)

}


/***********************************************************************************************/
/**
 * \fn void precond_sum_famg_add2 (REAL *r, REAL *z, void *data)
 *
 * \brief fractional AMG preconditioner (additive AMG) for operators of type
 *        alpha * D^s + beta * D^(1+s)
 *        where D is a laplacian and s \in (-1,0)
 *        This consists of
 *        FAMG2(s/2) * AMG(alpha * lump(M)^-1 + beta * lump(M)^-1 A lump(M)^-1) * FAMG2(s/2)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Ana Budisa
 * \date   2020-06-18
 *
 */
void precond_sum_famg_add2(REAL *r,
                           REAL *z,
                           void *data)
{
    // data[0] - FAMG(Adiv^1+s/2); data[1] - Grad; data[2] - Grad^T; data[3] - AMG;
    precond_data *pcdata=(precond_data *)data;

    // first, Grad^T * FAMG(Adiv^1+s/2) * Grad
    INT m1 = pcdata[0].mgl_data[0].A.row;
    REAL *x1 = (REAL *)calloc(m1, sizeof(REAL));
    dvector temp1 = dvec_create(m1), temp2 = dvec_create(m1);
    //precond_famg_add2(r, x1, pcdata); // x1 is the result (preconditioned r)
    // direct solve instead
    dcsr_mxv(pcdata[1].A, r, temp1.val); // temp1 = Grad * r
    // direct solve Adiv^1+s/2 temp2 = temp1
    INT status;
    status = hazmath_solve(pcdata[0].A, &temp1, &temp2, pcdata[0].mgl_data[0].Numeric, pcdata[0].print_level);
    if(status) printf("Direct solve status: %lld \n",   (long long )status);
    dcsr_mxv(pcdata[2].A, temp2.val, x1); // x1 = Grad^T * temp2


    // second AMG(alpha * lump(M)^-1 + beta * lump(M)^-1 A lump(M)^-1)
    INT m2 = pcdata[3].mgl_data[0].A.row;
    REAL *x2 = (REAL *)calloc(m2, sizeof(REAL));

    precond_amg(x1, x2, &pcdata[3]); // x2 is the result (preconditioned x1)

    // third Grad^T * FAMG(Adiv^1+s/2) * Grad
    dvec_set(m1, &temp1, 0.0); dvec_set(m1, &temp2, 0.0);
    // precond_famg_add2(x2, z, pcdata); // z is the result (preconditioned x2)
    // direct solve instead
    dcsr_mxv(pcdata[1].A, x2, temp1.val); // temp1 = Grad * x2
    // direct solve Adiv^1+s/2 temp2 = temp1
    status = hazmath_solve(pcdata[0].A, &temp1, &temp2, pcdata[0].mgl_data[0].Numeric, pcdata[0].print_level);
    if(status) printf("Direct solve status: %lld \n",   (long long )status);
    dcsr_mxv(pcdata[2].A, temp2.val, z); // z = Grad^T * temp2
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
 * \fn void precond_hx_div_additive_2D (REAL *r, REAL *z, void *data)
 *
 * \brief HX preconditioner for H(div): additive version in 2D
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/03/2019
 */
void precond_hx_div_additive_2D(REAL *r,
                                REAL *z,
                                void *data)
{
    //printf("HX div additive precond 2D\n");
    HX_div_data *hxdivdata=(HX_div_data *)data;
    INT n = hxdivdata->A->row;
    SHORT smooth_iter = hxdivdata->smooth_iter;

    // make sure z is initialzied by zeros
    array_set(n, z, 0.0);

    // local variable
    dvector zz;
    zz.row = n; zz.val = z;
    dvector rr;
    rr.row = n; rr.val = r;

    SHORT maxit, i;

    //-----------
    // smoothing
    //-----------
    smoother_dcsr_sgs(&zz, hxdivdata->A, &rr, smooth_iter);
    //smoother_dcsr_jacobi(&zz, 0, n, 1, hxdivdata->A, &rr, smooth_iter);

    //----------------------------
    // solve div vector Laplacian
    //----------------------------
    AMG_param *amgparam_divgrad = hxdivdata->amgparam_divgrad;
    AMG_data *mgl_divgrad = hxdivdata->mgl_divgrad;
    maxit = amgparam_divgrad->maxit;

    mgl_divgrad->b.row = hxdivdata->A_divgrad->row;
    dcsr_mxv(hxdivdata->Pt_div, r, mgl_divgrad->b.val);
    mgl_divgrad->x.row = hxdivdata->A_divgrad->row;
    dvec_set(hxdivdata->A_divgrad->row, &mgl_divgrad->x, 0.0);

    for (i=0;i<maxit;++i) mgcycle(mgl_divgrad, amgparam_divgrad);

    //----------------------------
    // update solution (additive)
    //----------------------------
    dcsr_aAxpy(1.0, hxdivdata->P_div, mgl_divgrad->x.val, z);

    //------------------------
    // solve scalar Laplacian
    //------------------------
    AMG_param *amgparam_grad = hxdivdata->amgparam_grad;
    AMG_data *mgl_grad = hxdivdata->mgl_grad;
    maxit = amgparam_grad->maxit;

    mgl_grad->b.row = hxdivdata->Curlt->row;
    dcsr_mxv(hxdivdata->Curlt, r, mgl_grad->b.val);
    mgl_grad->x.row=hxdivdata->A_grad->row;
    dvec_set(hxdivdata->A_grad->row, &mgl_grad->x, 0.0);

    for (i=0;i<maxit;++i) mgcycle(mgl_grad, amgparam_grad);

    //---------------------------
    // update solution (additive)
    //---------------------------
    dcsr_aAxpy(1.0, hxdivdata->Curl, mgl_grad->x.val, z);

}

/***********************************************************************************************/
/**
 * \fn void precond_hx_div_multiplicative_2D (REAL *r, REAL *z, void *data)
 *
 * \brief HX preconditioner for H(div): multiplicative version in 2D
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/03/2019
 */
void precond_hx_div_multiplicative_2D(REAL *r,
                                      REAL *z,
                                      void *data)
{
    //printf("HX div multiplicative precond 2D\n");
    HX_div_data *hxdivdata=(HX_div_data *)data;
    INT n = hxdivdata->A->row;
    SHORT smooth_iter = hxdivdata->smooth_iter;

    // backup r
    array_cp(n, r, hxdivdata->backup_r);

    // make sure z is initialzied by zeros
    array_set(n, z, 0.0);

    // local variable
    dvector zz;
    zz.row = n; zz.val = z;
    dvector rr;
    rr.row = n; rr.val = r;

    //SHORT maxit;

    //-----------
    // smoothing
    //-----------
    smoother_dcsr_sgs(&zz, hxdivdata->A, &rr, smooth_iter);
    //smoother_dcsr_jacobi(&zz, 0, n, 1, hxdivdata->A, &rr, smooth_iter);

    //----------
    // update r
    //----------
    dcsr_aAxpy(-1.0, hxdivdata->A, zz.val, rr.val);

    //----------------------------
    // solve div vector Laplacian
    //----------------------------
    AMG_param *amgparam_divgrad = hxdivdata->amgparam_divgrad;
    AMG_data *mgl_divgrad = hxdivdata->mgl_divgrad;
    //maxit = amgparam_divgrad->maxit;

    /*
    // solve with AMG only
    mgl_divgrad->b.row = hxdivdata->A_divgrad->row;
    dcsr_mxv(hxdivdata->Pt_div, r, mgl_divgrad->b.val);
    mgl_divgrad->x.row = hxdivdata->A_divgrad->row;
    dvec_set(hxdivdata->A_divgrad->row, &mgl_divgrad->x, 0.0);

    for (i=0;i<100;++i) mgcycle(mgl_divgrad, amgparam_divgrad);
    //dcsr_pvfgmres(hxdivdata->A_divgrad, &mgl_divgrad->b, &mgl_divgrad->x, NULL, 1e-3, 1000, 1000, 1, 1);
    //directsolve_HAZ(hxdivdata->A_divgrad, &(mgl_divgrad->b), &(mgl_divgrad->x), 1);

    //-----------------
    // update solution
    //-----------------
    dcsr_aAxpy(1.0, hxdivdata->P_div, mgl_divgrad->x.val, z);
    */

    // solve wth AMG+krylov
    dvector r_divgrad;
    dvec_alloc(hxdivdata->A_divgrad->row, &r_divgrad);
    dcsr_mxv(hxdivdata->Pt_div, r, r_divgrad.val);

    dvector x_divgrad;
    dvec_alloc(hxdivdata->A_divgrad->row, &x_divgrad);
    dvec_set(hxdivdata->A_divgrad->row, &x_divgrad, 0.0);

    precond_data pcdata_divgrad;
    param_amg_to_prec(&pcdata_divgrad,amgparam_divgrad);
    precond pc_divgrad;
    pc_divgrad.fct = precond_amg;

    pcdata_divgrad.max_levels = mgl_divgrad->num_levels;
    pcdata_divgrad.mgl_data = mgl_divgrad;

    pc_divgrad.data = &pcdata_divgrad;

    dcsr_pvfgmres(hxdivdata->A_divgrad, &r_divgrad, &x_divgrad, &pc_divgrad, 1e-3, 100, 100, 1, 0);

    dcsr_aAxpy(1.0, hxdivdata->P_div, x_divgrad.val, z);

    dvec_free(&r_divgrad);
    dvec_free(&x_divgrad);

    //----------
    // update r
    //----------
    array_cp(n, hxdivdata->backup_r, r);
    dcsr_aAxpy(-1.0, hxdivdata->A, zz.val, rr.val);

    //------------------------
    // solve scalar Laplacian
    //------------------------
    AMG_param *amgparam_grad = hxdivdata->amgparam_grad;
    AMG_data *mgl_grad = hxdivdata->mgl_grad;
    //maxit = amgparam_grad->maxit;

    /*
    // solve with AMG only
    mgl_grad->b.row = hxdivdata->Curlt->row;
    dcsr_mxv(hxdivdata->Curlt, r, mgl_grad->b.val);
    mgl_grad->x.row=hxdivdata->A_grad->row;
    dvec_set(hxdivdata->A_grad->row, &mgl_grad->x, 0.0);

    for (i=0;i<100;++i) mgcycle(mgl_grad, amgparam_grad);
    //dcsr_pvfgmres(hxdivdata->A_curlgrad, &mgl_curlgrad->b, &mgl_curlgrad->x, NULL, 1e-3, 1000, 1000, 1, 1);
    //directsolve_HAZ(hxdivdata->A_curlgrad, &(mgl_curlgrad->b), &(mgl_curlgrad->x),1);

    //-----------------
    // update solution
    //-----------------
    dcsr_aAxpy(1.0, hxdivdata->Curl, mgl_grad->x.val, z);
    */

    // solve wth AMG+krylov
    dvector r_grad;
    dvec_alloc(hxdivdata->A_grad->row, &r_grad);
    dcsr_mxv(hxdivdata->Curlt, r, r_grad.val);

    dvector x_grad;
    dvec_alloc(hxdivdata->A_grad->row, &x_grad);
    dvec_set(hxdivdata->A_grad->row, &x_grad, 0.0);

    precond_data pcdata_grad;
    param_amg_to_prec(&pcdata_grad,amgparam_grad);
    precond pc_grad;
    pc_grad.fct = precond_amg;

    pcdata_grad.max_levels = mgl_grad->num_levels;
    pcdata_grad.mgl_data = mgl_grad;

    pc_grad.data = &pcdata_grad;

    dcsr_pvfgmres(hxdivdata->A_grad, &r_grad, &x_grad, &pc_grad, 1e-3, 100, 100, 1, 0);

    dcsr_aAxpy(1.0, hxdivdata->Curl, x_grad.val, z);

    // restore r
    array_cp(n, hxdivdata->backup_r, r);

}


/***********************************************************************************************/
/**
 * \fn void precond_hx_div_additive (REAL *r, REAL *z, void *data)
 *
 * \brief HX preconditioner for H(div): additive version
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Peter Ohm
 * \date   05/10/2018
 */
void precond_hx_div_additive(REAL *r,
                              REAL *z,
                              void *data)
{
    //printf("HX div additive precond\n");
    HX_div_data *hxdivdata=(HX_div_data *)data;
    INT n = hxdivdata->A->row;
    SHORT smooth_iter = hxdivdata->smooth_iter;

    // make sure z is initialzied by zeros
    array_set(n, z, 0.0);

    // local variable
    dvector zz;
    zz.row = n; zz.val = z;
    dvector rr;
    rr.row = n; rr.val = r;

    SHORT maxit, i;

    // smoothing
    smoother_dcsr_sgs(&zz, hxdivdata->A, &rr, smooth_iter);
    //smoother_dcsr_jacobi(&zz, 0, n, 1, hxdivdata->A, &rr, smooth_iter);

    // solve div vector Laplacian
    AMG_param *amgparam_divgrad = hxdivdata->amgparam_divgrad;
    AMG_data *mgl_divgrad = hxdivdata->mgl_divgrad;
    maxit = amgparam_divgrad->maxit;

    mgl_divgrad->b.row = hxdivdata->A_divgrad->row;
    dcsr_mxv(hxdivdata->Pt_div, r, mgl_divgrad->b.val);
    mgl_divgrad->x.row=hxdivdata->A_divgrad->row;
    dvec_set(hxdivdata->A_divgrad->row, &mgl_divgrad->x, 0.0);

    for (i=0;i<maxit;++i) mgcycle(mgl_divgrad, amgparam_divgrad);
    //dcsr_pvfgmres(hxdivdata->A_divgrad, &mgl_divgrad->b, &mgl_divgrad->x, NULL, 1e-3, 1000, 1000, 1, 1);
    //directsolve_HAZ(hxdivdata->A_divgrad, &(mgl_divgrad->b), &(mgl_divgrad->x), 1);

    dcsr_aAxpy(1.0, hxdivdata->P_div, mgl_divgrad->x.val, z);

    /*
    INT j;
    for(j=0;j<n;j++){
      if(z[j]!=z[j]){ printf("DIV z[%lld]=%f\n",j,z[j]);}
    }
    */

    // smoothing
    REAL *temp1 = (REAL*)calloc(hxdivdata->Curlt->row,sizeof(REAL));
    REAL *temp2 = (REAL*)calloc(hxdivdata->Curlt->row,sizeof(REAL));

    dvector Cz;
    Cz.row = hxdivdata->A_curl->row;
    Cz.val = temp1;// initial guess is zero

    dvector Cr;
    Cr.row = Cz.row;
    Cr.val = temp2;

    dcsr_mxv(hxdivdata->Curlt,r,Cr.val);
//    INT flag;
//    flag = directsolve_HAZ(hxdivdata->A_curl, &Cr, &Cz, 1);
//    printf("flag=%lld\n",flag);
    smoother_dcsr_sgs(&Cz, hxdivdata->A_curl, &Cr, smooth_iter);
//    //smoother_dcsr_jacobi(&Cz, 0, Cz.row, 1, hxdivdata->A_curl, &Cr, smooth_iter);
    dcsr_aAxpy(1.0,hxdivdata->Curl,Cz.val,z);

    // solve scalar Laplacian
    AMG_param *amgparam_curlgrad = hxdivdata->amgparam_curlgrad;
    AMG_data *mgl_curlgrad = hxdivdata->mgl_curlgrad;
    maxit = amgparam_curlgrad->maxit;

    REAL *temp = (REAL*)calloc(hxdivdata->Curlt->row,sizeof(REAL));
    dcsr_mxv(hxdivdata->Curlt, r, temp);
    mgl_curlgrad->b.row = hxdivdata->Pt_curl->row;
    dcsr_mxv(hxdivdata->Pt_curl, temp, mgl_curlgrad->b.val);
    dvec_set(hxdivdata->A_curlgrad->row, &mgl_curlgrad->x, 0.0);

    for (i=0;i<maxit;++i) mgcycle(mgl_curlgrad, amgparam_curlgrad);
    //dcsr_pvfgmres(hxdivdata->A_curlgrad, &mgl_curlgrad->b, &mgl_curlgrad->x, NULL, 1e-3, 1000, 1000, 1, 1);
    //directsolve_HAZ(hxdivdata->A_curlgrad, &(mgl_curlgrad->b), &(mgl_curlgrad->x),1);

    dcsr_mxv(hxdivdata->P_curl, mgl_curlgrad->x.val, temp);
    dcsr_aAxpy(1.0, hxdivdata->Curl, temp, z);
    /*
    for(j=0;j<n;j++){
      if(z[j]!=z[j]){ printf("z[%lld]=%f\n",j,z[j]);}
    }
    */

    // free
    free(temp);
    free(temp1);
    free(temp2);
}

/***********************************************************************************************/
/**
 * \fn void precond_hx_div_multiplicative (REAL *r, REAL *z, void *data)
 *
 * \brief HX preconditioner for H(div): multiplicative version
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Peter Ohm
 * \date   05/10/2018
 */
void precond_hx_div_multiplicative(REAL *r,
                              REAL *z,
                              void *data)
{
    //printf("Multiplicative\n");
    //--------------------------
    // preparation
    //--------------------------
    HX_div_data *hxdivdata=(HX_div_data *)data;
    INT n = hxdivdata->A->row;
    SHORT smooth_iter = hxdivdata->smooth_iter;

    // backup r
    array_cp(n, r, hxdivdata->backup_r);

    // make sure z is initialzied by zeros
    array_set(n, z, 0.0);

    // local variable
    dvector zz;
    zz.row = n; zz.val = z;
    dvector rr;
    rr.row = n; rr.val = r;

    SHORT maxit, i;

    //--------------------------
    // smoothing
    //--------------------------
    smoother_dcsr_sgs(&zz, hxdivdata->A, &rr, smooth_iter);
    //smoother_dcsr_jacobi(&zz, 0, n, 1, hxdivdata->A, &rr, smooth_iter);

    // update r
    dcsr_aAxpy(-1.0, hxdivdata->A, zz.val, rr.val);

    //----------------------------
    // solve div vector Laplacian
    //----------------------------
    //printf("div grad\n");
    AMG_param *amgparam_divgrad = hxdivdata->amgparam_divgrad;
    AMG_data *mgl_divgrad = hxdivdata->mgl_divgrad;
    maxit = amgparam_divgrad->maxit;

    // solve with AMG only
    mgl_divgrad->b.row = hxdivdata->A_divgrad->row;
    dcsr_mxv(hxdivdata->Pt_div, r, mgl_divgrad->b.val);

    mgl_divgrad->x.row=hxdivdata->A_divgrad->row;
    dvec_set(hxdivdata->A_divgrad->row, &mgl_divgrad->x, 0.0);

    for (i=0;i<maxit;++i) mgcycle(mgl_divgrad, amgparam_divgrad);
    //dcsr_pvfgmres(hxdivdata->A_divgrad, &mgl_divgrad->b, &mgl_divgrad->x, NULL, 1e-6, 1000, 1000, 1, 1);
    //directsolve_HAZ(hxdivdata->A_divgrad, &(mgl_divgrad->b), &(mgl_divgrad->x), 0);

    dcsr_aAxpy(1.0, hxdivdata->P_div, mgl_divgrad->x.val, z);

    /*
    // solve wth AMG+krylov
    dvector r_divgrad;
    dvec_alloc(hxdivdata->A_divgrad->row, &r_divgrad);
    dcsr_mxv(hxdivdata->Pt_div, r, r_divgrad.val);

    dvector x_divgrad;
    dvec_alloc(hxdivdata->A_divgrad->row, &x_divgrad);
    dvec_set(hxdivdata->A_divgrad->row, &x_divgrad, 0.0);

    precond_data pcdata_divgrad;
    param_amg_to_prec(&pcdata_divgrad,amgparam_divgrad);
    precond pc_divgrad;
    pc_divgrad.fct = precond_amg;

    pcdata_divgrad.max_levels = mgl_divgrad->num_levels;
    pcdata_divgrad.mgl_data = mgl_divgrad;

    pc_divgrad.data = &pcdata_divgrad;

    dcsr_pvfgmres(hxdivdata->A_divgrad, &r_divgrad, &x_divgrad, &pc_divgrad, 1e-3, 100, 100, 1, 0);

    dcsr_aAxpy(1.0, hxdivdata->P_div, x_divgrad.val, z);

    dvec_free(&r_divgrad);
    dvec_free(&x_divgrad);
    */

    // update r
    array_cp(n, hxdivdata->backup_r, r);
    dcsr_aAxpy(-1.0, hxdivdata->A, zz.val, rr.val);

    //--------------------------
    // smoothing
    //--------------------------
    dvector Cz;
    Cz.row = hxdivdata->A_curl->row;
    Cz.val = hxdivdata->w;
    dvec_set(Cz.row, &Cz, 0.0);

    dvector Cr;
    Cr.row = Cz.row;
    Cr.val = hxdivdata->w+Cz.row;

    dcsr_mxv(hxdivdata->Curlt,r,Cr.val);
    smoother_dcsr_sgs(&Cz, hxdivdata->A_curl, &Cr, smooth_iter);
    //smoother_dcsr_jacobi(&Cz, 0, Cz.row, 1, hxdivdata->A_curl, &Cr, smooth_iter);
    dcsr_aAxpy(1.0,hxdivdata->Curl,Cz.val,z);

    // update r
    array_cp(n, hxdivdata->backup_r, r);
    dcsr_aAxpy(-1.0, hxdivdata->A, zz.val, rr.val);

    //-----------------------------
    // solve vector curl Laplacian
    //-----------------------------
    //printf("curl grad\n");
    AMG_param *amgparam_curlgrad = hxdivdata->amgparam_curlgrad;
    AMG_data *mgl_curlgrad = hxdivdata->mgl_curlgrad;

    REAL *temp = hxdivdata->w;
    dcsr_mxv(hxdivdata->Curlt, r, temp);

    // solve with AMG only
    mgl_curlgrad->b.row = hxdivdata->Pt_curl->row;
    dcsr_mxv(hxdivdata->Pt_curl, temp, mgl_curlgrad->b.val);
    dvec_set(hxdivdata->A_curlgrad->row, &mgl_curlgrad->x, 0.0);

    for (i=0;i<maxit;++i) mgcycle(mgl_curlgrad, amgparam_curlgrad);
    //dcsr_pvfgmres(hxdivdata->A_curlgrad, &mgl_curlgrad->b, &mgl_curlgrad->x, NULL, 1e-6, 1000, 1000, 1, 1);
    //directsolve_HAZ(hxdivdata->A_curlgrad, &(mgl_curlgrad->b), &(mgl_curlgrad->x),0);

    dcsr_mxv(hxdivdata->P_curl, mgl_curlgrad->x.val, temp);
    dcsr_aAxpy(1.0, hxdivdata->Curl, temp, z);

    /*
    // solve wth AMG+krylov
    dvector r_curlgrad;
    dvec_alloc(hxdivdata->A_curlgrad->row, &r_curlgrad);
    dcsr_mxv(hxdivdata->Pt_curl, temp, r_curlgrad.val);

    dvector x_curlgrad;
    dvec_alloc(hxdivdata->A_curlgrad->row, &x_curlgrad);
    dvec_set(hxdivdata->A_curlgrad->row, &x_curlgrad, 0.0);

    precond_data pcdata_curlgrad;
    param_amg_to_prec(&pcdata_curlgrad,amgparam_curlgrad);
    precond pc_curlgrad;
    pc_curlgrad.fct = precond_amg;

    pcdata_curlgrad.max_levels = mgl_curlgrad->num_levels;
    pcdata_curlgrad.mgl_data = mgl_curlgrad;

    pc_curlgrad.data = &pcdata_curlgrad;

    dcsr_pvfgmres(hxdivdata->A_curlgrad, &r_curlgrad, &x_curlgrad, &pc_curlgrad, 1e-3, 100, 100, 1, 0);

    dcsr_mxv(hxdivdata->P_curl, x_curlgrad.val, temp);
    dcsr_aAxpy(1.0, hxdivdata->Curl, temp, z);

    dvec_free(&r_curlgrad);
    dvec_free(&x_curlgrad);
    */

    //--------------------------
    // store r
    //--------------------------
    //printf("restore r\n");
    array_cp(n, hxdivdata->backup_r, r);
    //printf("preconditioner done\n");

}

/***********************************************************************************************/
/**
 * \fn void precond_dbsr_diag (REAL *r, REAL *z, void *data)
 *
 * \brief Diagonal preconditioner z=inv(D)*r
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   10/26/2010
 *
 * \note Works for general nb (Xiaozhe)
 */
void precond_dbsr_diag (REAL *r,
                        REAL *z,
                        void *data)
{
    precond_diag_bsr * diag = (precond_diag_bsr *)data;
    const INT nb = diag->nb;

    REAL *diagptr = diag->diag.val;
    const INT nb2 = nb*nb;
    const INT m = diag->diag.row/nb2;
    INT i;

    for (i = 0; i < m; ++i) {
        ddense_mxv(&(diagptr[i*nb2]),&(r[i*nb]),&(z[i*nb]),nb);
    }



}

/***********************************************************************************************/
/**
 * \fn void precond_dbsr_amg(REAL *r, REAL *z, void *data)
 *
 * \brief AMG preconditioner
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   08/07/2011
 */
void precond_dbsr_amg(REAL *r,
                      REAL *z,
                      void *data)
{
    precond_data_bsr *predata=(precond_data_bsr *)data;
    const INT row=predata->mgl_data[0].A.ROW;
    const INT nb = predata->mgl_data[0].A.nb;
    const INT maxit=predata->maxit;
    const INT m = row*nb;

	INT i;

    AMG_param amgparam; param_amg_init(&amgparam);
    amgparam.cycle_type = predata->cycle_type;
    amgparam.smoother   = predata->smoother;
    amgparam.presmooth_iter  = predata->presmooth_iter;
    amgparam.postsmooth_iter = predata->postsmooth_iter;
    amgparam.relaxation = predata->relaxation;
    amgparam.coarse_scaling = predata->coarse_scaling;
    amgparam.tentative_smooth = predata->tentative_smooth;

    AMG_data_bsr *mgl = predata->mgl_data;
    mgl->b.row=m; array_cp(m,r,mgl->b.val); // residual is an input
    mgl->x.row=m; dvec_set(m,&mgl->x,0.0);

    for ( i=maxit; i--; ) mgcycle_bsr(mgl,&amgparam);

    array_cp(m,mgl->x.val,z);
}


/***********************************************************************************************/
/**
 * \fn void precond_bdcsr_amg(REAL *r, REAL *z, void *data)
 *
 * \brief AMG preconditioner
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/16/2022
 */
void precond_bdcsr_amg(REAL *r,
                       REAL *z,
                       void *data)
{
    precond_data_bdcsr *predata=(precond_data_bdcsr *)data;
    //const INT brow=predata->mgl_data[0].A.brow;
  //  const INT bcol=predata->mgl_data[0].A.bcol;
    const INT maxit=predata->maxit;
    const INT total_row = predata->total_row;
    const INT total_col = predata->total_col;

	INT i;

    AMG_param amgparam; param_amg_init(&amgparam);
    amgparam.cycle_type = predata->cycle_type;
    amgparam.smoother   = predata->smoother;
    amgparam.presmooth_iter  = predata->presmooth_iter;
    amgparam.postsmooth_iter = predata->postsmooth_iter;
    amgparam.relaxation = predata->relaxation;
    amgparam.coarse_solver = predata->coarse_solver;
    amgparam.coarse_scaling = predata->coarse_scaling;
    amgparam.tentative_smooth = predata->tentative_smooth;

    AMG_data_bdcsr *mgl = predata->mgl_data;
    mgl->b.row=total_row; array_cp(total_row, r, mgl->b.val); // residual is an input
    mgl->x.row=total_col; dvec_set(total_col, &mgl->x, 0.0);

    for ( i=maxit; i--; ) mgcycle_bdcsr(mgl,&amgparam);

    array_cp(total_col, mgl->x.val, z);
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
  //#if WITH_SUITESPARSE
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
  hazmath_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);
  //#endif

  // Preconditioning A11 block
  //#if  WITH_UMFPACK
  /* use UMFPACK direct solver */
  hazmath_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);
  //#endif

  // restore r
  array_cp(N, tempr->val, r);

  //#endif
}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_2_amg (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved by AMG)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   08/08/2018
 */
void precond_block_diag_2_amg(REAL *r,
                          REAL *z,
                          void *data)
{
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
    INT i;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, z0, z1;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;

    r0.val = r; r1.val = &(r[N0]);
    z0.val = z; z1.val = &(z[N0]);

    // Preconditioning A00 block
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
    array_cp(N0, mgl[0]->x.val, z0.val);

    // Preconditioning A11 block
    mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
    mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
    array_cp(N1, mgl[1]->x.val, z1.val);

    // restore r
    array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_2_amg_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved by AMG preconditioned krylov method)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   08/08/2018
 */
void precond_block_diag_2_amg_krylov(REAL *r,
                          REAL *z,
                          void *data)
{
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
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, z0, z1;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;

    r0.val = r; r1.val = &(r[N0]);
    z0.val = z; z1.val = &(z[N0]);

    precond_data pcdata;
    param_amg_to_prec(&pcdata,amgparam);
    precond pc;
    pc.fct = precond_amg;

    // Preconditioning A00 block
    pcdata.max_levels = mgl[0][0].num_levels;
    pcdata.mgl_data = mgl[0];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc, 1e-2, 100, 100, 1, 1);

    // Preconditioning A11 block
    pcdata.max_levels = mgl[1][0].num_levels;
    pcdata.mgl_data = mgl[1];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc, 1e-2, 100, 100, 1, 1);

    // restore r
    array_cp(N, tempr->val, r);

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
  //#if WITH_SUITESPARSE

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
  hazmath_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);
  //#endif

  // r1 = r1 - A2*z0
  dcsr_aAxpy(-1.0, A->blocks[2], z0.val, r1.val);

  // Preconditioning A11 block
  //#if  WITH_UMFPACK
  /* use UMFPACK direct solver */
  hazmath_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);
  //#endif

  // restore r
  array_cp(N, tempr->val, r);

  //#endif

}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_2_amg (REAL *r, REAL *z, void *data)
 * \brief block lower diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved by AMG)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   08/08/2018
 */
void precond_block_lower_2_amg(REAL *r,
                          REAL *z,
                          void *data)
{
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_diag = precdata->A_diag;
    dvector *tempr = &(precdata->r);

    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N = N0 + N1;

    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);

    // prepare
    INT i;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, z0, z1;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;

    r0.val = r; r1.val = &(r[N0]);
    z0.val = z; z1.val = &(z[N0]);

    // Preconditioning A00 block
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
    array_cp(N0, mgl[0]->x.val, z0.val);

    // r1 = r1 - A2*z0
    if (A->blocks[2] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[2], z0.val, r1.val);

    // Preconditioning A11 block
    mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
    mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
    array_cp(N1, mgl[1]->x.val, z1.val);

    // restore r
    array_cp(N, tempr->val, r);

}


/***********************************************************************************************/
/**
 * \fn void precond_block_lower_2_amg_krylov (REAL *r, REAL *z, void *data)
 * \brief block lower diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved by AMG preconditioned krylov method)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   08/08/2018
 */
void precond_block_lower_2_amg_krylov(REAL *r,
                          REAL *z,
                          void *data)
{
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_diag = precdata->A_diag;
    dvector *tempr = &(precdata->r);

    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    //const INT N2 = A_diag[2].row;
    const INT N = N0 + N1;

    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);

    // prepare
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, z0, z1;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;

    r0.val = r; r1.val = &(r[N0]);
    z0.val = z; z1.val = &(z[N0]);

    precond_data pcdata;
    param_amg_to_prec(&pcdata,amgparam);
    precond pc;
    pc.fct = precond_amg;

    // Preconditioning A00 block
    pcdata.max_levels = mgl[0][0].num_levels;
    pcdata.mgl_data = mgl[0];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc, 1e-2, 100, 100, 1, 1);

    // r1 = r1 - A2*z0
    if (A->blocks[2] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[2], z0.val, r1.val);

    // Preconditioning A11 block
    pcdata.max_levels = mgl[1][0].num_levels;
    pcdata.mgl_data = mgl[1];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc, 1e-2, 100, 100, 1, 1);

    // restore r
    array_cp(N, tempr->val, r);

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
  //#if WITH_SUITESPARSE

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
  hazmath_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);
  //#endif

  // r0 = r0 - A1*z1
  dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);

  // Preconditioning A00 block
  //#if  WITH_UMFPACK
  /* use UMFPACK direct solver */
  hazmath_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);
  //#endif

  // restore r
  array_cp(N, tempr->val, r);

  //#endif

}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_2_amg (REAL *r, REAL *z, void *data)
 * \brief block upper diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved by AMG)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   08/08/2018
 */
void precond_block_upper_2_amg(REAL *r,
                          REAL *z,
                          void *data)
{
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_diag = precdata->A_diag;
    dvector *tempr = &(precdata->r);

    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N = N0 + N1;

    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);

    // prepare
    INT i;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, z0, z1;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;

    r0.val = r; r1.val = &(r[N0]);
    z0.val = z; z1.val = &(z[N0]);

    // Preconditioning A11 block
    mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
    mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
    array_cp(N1, mgl[1]->x.val, z1.val);

    // r0 = r0 - A1*z1
    if (A->blocks[1] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);

    // Preconditioning A00 block
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
    array_cp(N0, mgl[0]->x.val, z0.val);

    // restore r
    array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_2_amg_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved by AMG preconditioned krylov method)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   08/08/2018
 */
void precond_block_upper_2_amg_krylov(REAL *r,
                          REAL *z,
                          void *data)
{
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dCSRmat *A_diag = precdata->A_diag;
    dvector *tempr = &(precdata->r);

    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N = N0 + N1;

    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);

    // prepare
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, z0, z1;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;

    r0.val = r; r1.val = &(r[N0]);
    z0.val = z; z1.val = &(z[N0]);

    precond_data pcdata;
    param_amg_to_prec(&pcdata,amgparam);
    precond pc;
    pc.fct = precond_amg;

    // Preconditioning A11 block
    pcdata.max_levels = mgl[1][0].num_levels;
    pcdata.mgl_data = mgl[1];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc, 1e-3, 100, 100, 1, 1);

    // r0 = r0 - A1*z1
    if (A->blocks[1] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);

    // Preconditioning A00 block
    pcdata.max_levels = mgl[0][0].num_levels;
    pcdata.mgl_data = mgl[0];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc, 1e-3, 100, 100, 1, 1);

    // restore r
    array_cp(N, tempr->val, r);

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
  //#if WITH_SUITESPARSE
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
    hazmath_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);

    // Preconditioning A11 block
    /* use UMFPACK direct solver */
    hazmath_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);

    // Preconditioning A22 block
    /* use UMFPACK direct solver */
    hazmath_solve(&A_diag[2], &r2, &z2, LU_diag[2], 0);

    // restore r
    array_cp(N, tempr->val, r);

    //#endif
}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_3_amg (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (3x3 block matrix, each diagonal block
 *        is solved by AMG)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   05/16/2018
 */
void precond_block_diag_3_amg(REAL *r,
                          REAL *z,
                          void *data)
{
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
    INT i;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, r2, z0, z1, z2;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;

    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);

    // Preconditioning A00 block
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
    array_cp(N0, mgl[0]->x.val, z0.val);

    // Preconditioning A11 block
    mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
    mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
    array_cp(N1, mgl[1]->x.val, z1.val);

    // Preconditioning A22 block
    mgl[2]->b.row=N2; array_cp(N2, r2.val, mgl[2]->b.val); // residual is an input
    mgl[2]->x.row=N2; dvec_set(N2, &mgl[2]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[2], amgparam);
    array_cp(N2, mgl[2]->x.val, z2.val);

    // restore r
    array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_3_amg_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (3x3 block matrix, each diagonal block
 *        is solved by AMG preconditioned krylov method)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   05/16/2018
 */
void precond_block_diag_3_amg_krylov(REAL *r,
                          REAL *z,
                          void *data)
{
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
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, r2, z0, z1, z2;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;

    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);

    precond_data pcdata;
    param_amg_to_prec(&pcdata,amgparam);
    precond pc;
    pc.fct = precond_amg;

    // Preconditioning A00 block
    pcdata.max_levels = mgl[0][0].num_levels;
    pcdata.mgl_data = mgl[0];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc, 1e-3, 100, 100, 1, 1);

    // Preconditioning A11 block
    pcdata.max_levels = mgl[1][0].num_levels;
    pcdata.mgl_data = mgl[1];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc, 1e-3, 100, 100, 1, 1);

    // Preconditioning A22 block
    pcdata.max_levels = mgl[2][0].num_levels;
    pcdata.mgl_data = mgl[2];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[2][0].A, &r2, &z2, &pc, 1e-3, 100, 100, 1, 1);

    // restore r
    array_cp(N, tempr->val, r);

}


/***********************************************************************************************/
/**
 * \fn void precond_block_lower_3 (REAL *r, REAL *z, void *data)
 * \brief block lower triangular preconditioning (3x3 block matrix, each diagonal
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
  //#if WITH_SUITESPARSE

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
    hazmath_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);

    // r1 = r1 - A3*z0
    if (A->blocks[3] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[3], z0.val, r1.val);

    // Preconditioning A11 block
    /* use UMFPACK direct solver */
    hazmath_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);

    // r2 = r2 - A6*z0 - A7*z1
    if (A->blocks[6] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[6], z0.val, r2.val);
    if (A->blocks[7] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[7], z1.val, r2.val);

    // Preconditioning A22 block
    /* use UMFPACK direct solver */
    hazmath_solve(&A_diag[2], &r2, &z2, LU_diag[2], 0);

    // restore r
    array_cp(N, tempr->val, r);

    //#endif

}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_3_amg (REAL *r, REAL *z, void *data)
 * \brief block lower diagonal preconditioning (3x3 block matrix, each diagonal block
 *        is solved by AMG)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   05/16/2018
 */
void precond_block_lower_3_amg(REAL *r,
                          REAL *z,
                          void *data)
{
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
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
    INT i;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, r2, z0, z1, z2;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;

    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);

    // Preconditioning A00 block
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
    array_cp(N0, mgl[0]->x.val, z0.val);

    // r1 = r1 - A3*z0
    if (A->blocks[3] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[3], z0.val, r1.val);

    // Preconditioning A11 block
    mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
    mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
    array_cp(N1, mgl[1]->x.val, z1.val);

    // r2 = r2 - A6*z0 - A7*z1
    if (A->blocks[6] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[6], z0.val, r2.val);
    if (A->blocks[7] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[7], z1.val, r2.val);

    // Preconditioning A22 block
    mgl[2]->b.row=N2; array_cp(N2, r2.val, mgl[2]->b.val); // residual is an input
    mgl[2]->x.row=N2; dvec_set(N2, &mgl[2]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[2], amgparam);
    array_cp(N2, mgl[2]->x.val, z2.val);

    // restore r
    array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_3_amg_krylov (REAL *r, REAL *z, void *data)
 * \brief block lower diagonal preconditioning (3x3 block matrix, each diagonal block
 *        is solved by AMG preconditioned krylov method)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   05/16/2018
 */
void precond_block_lower_3_amg_krylov(REAL *r,
                          REAL *z,
                          void *data)
{
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
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
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, r2, z0, z1, z2;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;

    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);

    precond_data pcdata;
    param_amg_to_prec(&pcdata,amgparam);
    precond pc;
    pc.fct = precond_amg;

    // Preconditioning A00 block
    pcdata.max_levels = mgl[0][0].num_levels;
    pcdata.mgl_data = mgl[0];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc, 1e-3, 100, 100, 1, 1);

    // r1 = r1 - A3*z0
    if (A->blocks[3] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[3], z0.val, r1.val);

    // Preconditioning A11 block
    pcdata.max_levels = mgl[1][0].num_levels;
    pcdata.mgl_data = mgl[1];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc, 1e-3, 100, 100, 1, 1);

    // r2 = r2 - A6*z0 - A7*z1
    if (A->blocks[6] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[6], z0.val, r2.val);
    if (A->blocks[7] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[7], z1.val, r2.val);

    // Preconditioning A22 block
    pcdata.max_levels = mgl[2][0].num_levels;
    pcdata.mgl_data = mgl[2];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[2][0].A, &r2, &z2, &pc, 1e-3, 100, 100, 1, 1);

    // restore r
    array_cp(N, tempr->val, r);

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
  //#if WITH_SUITESPARSE

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
    hazmath_solve(&A_diag[2], &r2, &z2, LU_diag[2], 0);

    // r1 = r1 - A5*z2
    if (A->blocks[5] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[5], z2.val, r1.val);

    // Preconditioning A11 block
    /* use UMFPACK direct solver */
    hazmath_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);

    // r0 = r0 - A1*z1 - A2*z2
    if (A->blocks[1] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
    if (A->blocks[2] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[2], z2.val, r0.val);

    // Preconditioning A00 block
    /* use UMFPACK direct solver */
    hazmath_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);

    // restore r
    array_cp(N, tempr->val, r);

    //#endif

}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_3_amg (REAL *r, REAL *z, void *data)
 * \brief block upper diagonal preconditioning (3x3 block matrix, each diagonal block
 *        is solved by AMG)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   05/16/2018
 */
void precond_block_upper_3_amg(REAL *r,
                          REAL *z,
                          void *data)
{
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
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
    INT i;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, r2, z0, z1, z2;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;

    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);

    // Preconditioning A22 block
    mgl[2]->b.row=N2; array_cp(N2, r2.val, mgl[2]->b.val); // residual is an input
    mgl[2]->x.row=N2; dvec_set(N2, &mgl[2]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[2], amgparam);
    array_cp(N2, mgl[2]->x.val, z2.val);

    // r1 = r1 - A5*z2
    if (A->blocks[5] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[5], z2.val, r1.val);

    // Preconditioning A11 block
    mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
    mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
    array_cp(N1, mgl[1]->x.val, z1.val);

    // r0 = r0 - A1*z1 - A2*z2
    if (A->blocks[1] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
    if (A->blocks[2] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[2], z2.val, r0.val);

    // Preconditioning A00 block
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
    array_cp(N0, mgl[0]->x.val, z0.val);

    // restore r
    array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_3_amg_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (3x3 block matrix, each diagonal block
 *        is solved by AMG preconditioned krylov method)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   05/16/2018
 */
void precond_block_upper_3_amg_krylov(REAL *r,
                          REAL *z,
                          void *data)
{
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
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
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, r2, z0, z1, z2;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;

    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);

    precond_data pcdata;
    param_amg_to_prec(&pcdata,amgparam);
    precond pc;
    pc.fct = precond_amg;

    // Preconditioning A22 block
    pcdata.max_levels = mgl[2][0].num_levels;
    pcdata.mgl_data = mgl[2];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[2][0].A, &r2, &z2, &pc, 1e-3, 100, 100, 1, 1);

    // r1 = r1 - A5*z2
    if (A->blocks[5] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[5], z2.val, r1.val);

    // Preconditioning A11 block
    pcdata.max_levels = mgl[1][0].num_levels;
    pcdata.mgl_data = mgl[1];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc, 1e-3, 100, 100, 1, 1);

    // r0 = r0 - A1*z1 - A2*z2
    if (A->blocks[1] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
    if (A->blocks[2] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[2], z2.val, r0.val);

    // Preconditioning A00 block
    pcdata.max_levels = mgl[0][0].num_levels;
    pcdata.mgl_data = mgl[0];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc, 1e-3, 100, 100, 1, 1);

    // restore r
    array_cp(N, tempr->val, r);

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
  //#if WITH_SUITESPARSE
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
    hazmath_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);

    // Preconditioning A11 block
    /* use UMFPACK direct solver */
    hazmath_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);

    // Preconditioning A22 block
    /* use UMFPACK direct solver */
    hazmath_solve(&A_diag[2], &r2, &z2, LU_diag[2], 0);

    // Preconditioning A33 block
    /* use UMFPACK direct solver */
    hazmath_solve(&A_diag[3], &r3, &z3, LU_diag[3], 0);

    // restore r
    array_cp(N, tempr->val, r);

    //#endif
}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_4_amg (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (4x4 block matrix, each diagonal block
 *        is solved by AMG)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   08/23/2021
 */
void precond_block_diag_4_amg(REAL *r,
                              REAL *z,
                              void *data)
{
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
    INT i;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
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
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
    array_cp(N0, mgl[0]->x.val, z0.val);

    // Preconditioning A11 block
    mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
    mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
    array_cp(N1, mgl[1]->x.val, z1.val);

    // Preconditioning A22 block
    mgl[2]->b.row=N2; array_cp(N2, r2.val, mgl[2]->b.val); // residual is an input
    mgl[2]->x.row=N2; dvec_set(N2, &mgl[2]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[2], amgparam);
    array_cp(N2, mgl[2]->x.val, z2.val);

    // Preconditioning A33 block
    mgl[3]->b.row=N3; array_cp(N3, r3.val, mgl[3]->b.val); // residual is an input
    mgl[3]->x.row=N3; dvec_set(N3, &mgl[3]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[3], amgparam);
    array_cp(N3, mgl[3]->x.val, z3.val);

    // restore r
    array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_4_amg_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (4x4 block matrix, each diagonal block
 *        is solved by AMG preconditioned Krylov method)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   08/23/2021
 */
void precond_block_diag_4_amg_krylov(REAL *r,
                                     REAL *z,
                                     void *data)
{
  //#if WITH_SUITESPARSE
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
    //    INT i;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
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

    precond_data pcdata;
    param_amg_to_prec(&pcdata,amgparam);
    precond pc;
    pc.fct = precond_amg;

    // Preconditioning A00 block
    pcdata.max_levels = mgl[0][0].num_levels;
    pcdata.mgl_data = mgl[0];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc, 1e-3, 100, 100, 1, 1);

    // Preconditioning A11 block
    pcdata.max_levels = mgl[1][0].num_levels;
    pcdata.mgl_data = mgl[1];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc, 1e-3, 100, 100, 1, 1);

    // Preconditioning A22 block
    pcdata.max_levels = mgl[2][0].num_levels;
    pcdata.mgl_data = mgl[2];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[2][0].A, &r2, &z2, &pc, 1e-3, 100, 100, 1, 1);

    // Preconditioning A33 block
    pcdata.max_levels = mgl[3][0].num_levels;
    pcdata.mgl_data = mgl[3];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[3][0].A, &r3, &z3, &pc, 1e-3, 100, 100, 1, 1);

    // restore r
    array_cp(N, tempr->val, r);

    //#endif
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
  //#if WITH_SUITESPARSE

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
    hazmath_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);

    // r1 = r1 - A4*z0
    if (A->blocks[4] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[4], z0.val, r1.val);

    // Preconditioning A11 block
    /* use UMFPACK direct solver */
    hazmath_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);

    // r2 = r2 - A8*z0 - A9*z1
    if (A->blocks[8] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[8], z0.val, r2.val);
    if (A->blocks[9] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[9], z1.val, r2.val);

    // Preconditioning A22 block
    /* use UMFPACK direct solver */
    hazmath_solve(&A_diag[2], &r2, &z2, LU_diag[2], 0);

    // r3 = r3 - A12*z0 - A13*z1 - A14*z2
    if (A->blocks[12] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[12], z0.val, r3.val);
    if (A->blocks[13] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[13], z1.val, r3.val);
    if (A->blocks[14] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[14], z2.val, r3.val);

    // Preconditioning A33 block
    /* use UMFPACK direct solver */
    hazmath_solve(&A_diag[3], &r3, &z3, LU_diag[3], 0);

    // restore r
    array_cp(N, tempr->val, r);

    //#endif

}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_4_amg (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (4x4 block matrix, each diagonal
 *        block is solved by AMG)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   08/23/2021
 *
 * A[0]  A[1]  A[2]  A[3]
 * A[4]  A[5]  A[6]  A[7]
 * A[8]  A[9]  A[10] A[11]
 * A[12] A[13] A[14] A[15]
 */
void precond_block_lower_4_amg(REAL *r,
                               REAL *z,
                               void *data)
{

    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
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
    INT i;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
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
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
    array_cp(N0, mgl[0]->x.val, z0.val);

    // r1 = r1 - A4*z0
    if (A->blocks[4] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[4], z0.val, r1.val);

    // Preconditioning A11 block
    mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
    mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
    array_cp(N1, mgl[1]->x.val, z1.val);

    // r2 = r2 - A8*z0 - A9*z1
    if (A->blocks[8] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[8], z0.val, r2.val);
    if (A->blocks[9] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[9], z1.val, r2.val);

    // Preconditioning A22 block
    mgl[2]->b.row=N2; array_cp(N2, r2.val, mgl[2]->b.val); // residual is an input
    mgl[2]->x.row=N2; dvec_set(N2, &mgl[2]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[2], amgparam);
    array_cp(N2, mgl[2]->x.val, z2.val);

    // r3 = r3 - A12*z0 - A13*z1 - A14*z2
    if (A->blocks[12] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[12], z0.val, r3.val);
    if (A->blocks[13] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[13], z1.val, r3.val);
    if (A->blocks[14] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[14], z2.val, r3.val);

    // Preconditioning A33 block
    mgl[3]->b.row=N3; array_cp(N3, r3.val, mgl[3]->b.val); // residual is an input
    mgl[3]->x.row=N3; dvec_set(N3, &mgl[3]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[3], amgparam);
    array_cp(N3, mgl[3]->x.val, z3.val);

    // restore r
    array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_4_amg_krylov (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (4x4 block matrix, each diagonal
 *        block is solved by AMG preconditioned Krylov method)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   08/23/2021
 *
 * A[0]  A[1]  A[2]  A[3]
 * A[4]  A[5]  A[6]  A[7]
 * A[8]  A[9]  A[10] A[11]
 * A[12] A[13] A[14] A[15]
 */
void precond_block_lower_4_amg_krylov(REAL *r,
                                      REAL *z,
                                      void *data)
{

    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
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
    //    INT i;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
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

    precond_data pcdata;
    param_amg_to_prec(&pcdata,amgparam);
    precond pc;
    pc.fct = precond_amg;

    // Preconditioning A00 block
    pcdata.max_levels = mgl[0][0].num_levels;
    pcdata.mgl_data = mgl[0];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc, 1e-3, 100, 100, 1, 1);

    // r1 = r1 - A4*z0
    if (A->blocks[4] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[4], z0.val, r1.val);

    // Preconditioning A11 block
    pcdata.max_levels = mgl[1][0].num_levels;
    pcdata.mgl_data = mgl[1];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc, 1e-3, 100, 100, 1, 1);

    // r2 = r2 - A8*z0 - A9*z1
    if (A->blocks[8] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[8], z0.val, r2.val);
    if (A->blocks[9] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[9], z1.val, r2.val);

    // Preconditioning A22 block
    pcdata.max_levels = mgl[2][0].num_levels;
    pcdata.mgl_data = mgl[2];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[2][0].A, &r2, &z2, &pc, 1e-3, 100, 100, 1, 1);

    // r3 = r3 - A12*z0 - A13*z1 - A14*z2
    if (A->blocks[12] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[12], z0.val, r3.val);
    if (A->blocks[13] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[13], z1.val, r3.val);
    if (A->blocks[14] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[14], z2.val, r3.val);

    // Preconditioning A33 block
    pcdata.max_levels = mgl[3][0].num_levels;
    pcdata.mgl_data = mgl[3];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[3][0].A, &r3, &z3, &pc, 1e-3, 100, 100, 1, 1);

    // restore r
    array_cp(N, tempr->val, r);

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
  //#if WITH_SUITESPARSE

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
    hazmath_solve(&A_diag[3], &r3, &z3, LU_diag[3], 0);

    // r2 = r2 - A11*z3
    if (A->blocks[11] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[11], z3.val, r2.val);

    // Preconditioning A22 block
    /* use UMFPACK direct solver */
    hazmath_solve(&A_diag[2], &r2, &z2, LU_diag[2], 0);

    // r1 = r1 - A6*z2 - A7*z3
    if (A->blocks[6] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[6], z2.val, r1.val);
    if (A->blocks[7] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[7], z3.val, r1.val);

    // Preconditioning A11 block
    /* use UMFPACK direct solver */
    hazmath_solve(&A_diag[1], &r1, &z1, LU_diag[1], 0);

    // r0 = r0 - A1*z1 - A2*z2 - A3*z3
    if (A->blocks[1] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
    if (A->blocks[2] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[2], z2.val, r0.val);
    if (A->blocks[3] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[3], z3.val, r0.val);

    // Preconditioning A00 block
    /* use UMFPACK direct solver */
    hazmath_solve(&A_diag[0], &r0, &z0, LU_diag[0], 0);

    // restore r
    array_cp(N, tempr->val, r);

    //#endif

}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_4_amg (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (4x4 block matrix, each diagonal
 *        block is solved by AMG)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   08/23/2021
 *
 * A[0]  A[1]  A[2]  A[3]
 * A[4]  A[5]  A[6]  A[7]
 * A[8]  A[9]  A[10] A[11]
 * A[12] A[13] A[14] A[15]
 */
void precond_block_upper_4_amg(REAL *r,
                               REAL *z,
                               void *data)
{

    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
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
    INT i;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
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
    mgl[3]->b.row=N3; array_cp(N3, r3.val, mgl[3]->b.val); // residual is an input
    mgl[3]->x.row=N3; dvec_set(N3, &mgl[3]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[3], amgparam);
    array_cp(N3, mgl[3]->x.val, z3.val);

    // r2 = r2 - A11*z3
    if (A->blocks[11] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[11], z3.val, r2.val);

    // Preconditioning A22 block
    mgl[2]->b.row=N2; array_cp(N2, r2.val, mgl[2]->b.val); // residual is an input
    mgl[2]->x.row=N2; dvec_set(N2, &mgl[2]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[2], amgparam);
    array_cp(N2, mgl[2]->x.val, z2.val);

    // r1 = r1 - A6*z2 - A7*z3
    if (A->blocks[6] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[6], z2.val, r1.val);
    if (A->blocks[7] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[7], z3.val, r1.val);

    // Preconditioning A11 block
    mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
    mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
    array_cp(N1, mgl[1]->x.val, z1.val);

    // r0 = r0 - A1*z1 - A2*z2 - A3*z3
    if (A->blocks[1] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
    if (A->blocks[2] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[2], z2.val, r0.val);
    if (A->blocks[3] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[3], z3.val, r0.val);

    // Preconditioning A00 block
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
    array_cp(N0, mgl[0]->x.val, z0.val);

    // restore r
    array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_4_amg_krylov (REAL *r, REAL *z, void *data)
 * \brief block upper triangular preconditioning (4x4 block matrix, each diagonal
 *        block is solved by AMG preconditioned Krylov method)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   08/23/2021
 *
 * A[0]  A[1]  A[2]  A[3]
 * A[4]  A[5]  A[6]  A[7]
 * A[8]  A[9]  A[10] A[11]
 * A[12] A[13] A[14] A[15]
 */
void precond_block_upper_4_amg_krylov(REAL *r,
                                      REAL *z,
                                      void *data)
{

    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
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
    //    INT i;
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
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

    precond_data pcdata;
    param_amg_to_prec(&pcdata,amgparam);
    precond pc;
    pc.fct = precond_amg;

    // Preconditioning A33 block
    pcdata.max_levels = mgl[3][0].num_levels;
    pcdata.mgl_data = mgl[3];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[3][0].A, &r3, &z3, &pc, 1e-3, 100, 100, 1, 1);

    // r2 = r2 - A11*z3
    if (A->blocks[11] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[11], z3.val, r2.val);

    // Preconditioning A22 block
    pcdata.max_levels = mgl[2][0].num_levels;
    pcdata.mgl_data = mgl[2];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[2][0].A, &r2, &z2, &pc, 1e-3, 100, 100, 1, 1);

    // r1 = r1 - A6*z2 - A7*z3
    if (A->blocks[6] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[6], z2.val, r1.val);
    if (A->blocks[7] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[7], z3.val, r1.val);

    // Preconditioning A11 block
    pcdata.max_levels = mgl[1][0].num_levels;
    pcdata.mgl_data = mgl[1];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc, 1e-3, 100, 100, 1, 1);

    // r0 = r0 - A1*z1 - A2*z2 - A3*z3
    if (A->blocks[1] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
    if (A->blocks[2] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[2], z2.val, r0.val);
    if (A->blocks[3] != NULL)
        dcsr_aAxpy(-1.0, A->blocks[3], z3.val, r0.val);

    // Preconditioning A00 block
    pcdata.max_levels = mgl[0][0].num_levels;
    pcdata.mgl_data = mgl[0];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc, 1e-3, 100, 100, 1, 1);

    // restore r
    array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_5_amg_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (5x5 block matrix, each diagonal block
 *        is solved by AMG preconditioned krylov method)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   06/19/2020
 */
void precond_block_diag_5_amg_krylov(REAL *r,
                          REAL *z,
                          void *data)
{
    precond_block_data *precdata=(precond_block_data *)data;
    dCSRmat *A_diag = precdata->A_diag;
    dvector *tempr = &(precdata->r);

    const INT N0 = A_diag[0].row;
    const INT N1 = A_diag[1].row;
    const INT N2 = A_diag[2].row;
    const INT N3 = A_diag[3].row;
    const INT N4 = A_diag[4].row;
    const INT N = N0 + N1 + N2 + N3 + N4;

    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);

    // prepare
    AMG_param *amgparam = precdata->amgparam;
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, r2, r3, r4, z0, z1, z2, z3, z4;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;
    r2.row = N2; z2.row = N2;
    r3.row = N3; z3.row = N3;
    r4.row = N4; z4.row = N4;

    r0.val = r; r1.val = &(r[N0]); r2.val = &(r[N0+N1]);
    r3.val = &(r[N0+N1+N2]); r4.val = &(r[N0+N1+N2+N3]);
    z0.val = z; z1.val = &(z[N0]); z2.val = &(z[N0+N1]);
    z3.val = &(z[N0+N1+N2]); z4.val = &(z[N0+N1+N2+N3]);

    precond_data pcdata;
    param_amg_to_prec(&pcdata,amgparam);
    precond pc;
    pc.fct = precond_amg;

    // Preconditioning A00 block
    //printf("block 0\n");
    pcdata.max_levels = mgl[0][0].num_levels;
    pcdata.mgl_data = mgl[0];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc, 1e-3, 100, 100, 1, 1);

    // Preconditioning A11 block
    //printf("block 1\n");
    pcdata.max_levels = mgl[1][0].num_levels;
    pcdata.mgl_data = mgl[1];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc, 1e-3, 100, 100, 1, 1);

    // Preconditioning A22 block
    //printf("block 2\n");
    pcdata.max_levels = mgl[2][0].num_levels;
    pcdata.mgl_data = mgl[2];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[2][0].A, &r2, &z2, &pc, 1e-3, 100, 100, 1, 1);

    // Preconditioning A33 block
    //printf("block 3\n");
    pcdata.max_levels = mgl[3][0].num_levels;
    pcdata.mgl_data = mgl[3];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[3][0].A, &r3, &z3, &pc, 1e-3, 100, 100, 1, 1);

    // Preconditioning A44 block
    //printf("block 4\n");
    pcdata.max_levels = mgl[4][0].num_levels;
    pcdata.mgl_data = mgl[4];

    pc.data = &pcdata;

    dcsr_pvfgmres(&mgl[4][0].A, &r4, &z4, &pc, 1e-3, 100, 100, 1, 1);

    // restore r
    array_cp(N, tempr->val, r);

}


/***********************************************************************************************/
/**
 * \fn void precond_block_diag(REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (nxn block matrix, each diagonal block
 *        is solved exactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   04/05/2017
 */
void precond_block_diag(REAL *r,
                          REAL *z,
                          void *data)
{
  //#if WITH_SUITESPARSE
  precond_block_data *precdata=(precond_block_data *)data;
  dCSRmat *A_diag = precdata->A_diag;
  dvector *tempr = &(precdata->r);
  block_dCSRmat *A = precdata->Abcsr;

  const INT nb = A->brow;
  const INT N = tempr->row;

  INT i, istart = 0;;

  //const INT N0 = A_diag[0].row;
  //const INT N1 = A_diag[1].row;

  // back up r, setup z;
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);

  // prepare
  void **LU_diag = precdata->LU_diag;
  dvector ri, zi;

  for (i=0; i<nb; i++)
  {

      // get dvector for solutions and right hand sides
      ri.row = A_diag[i].row;
      zi.row = A_diag[i].row;

      ri.val = &(r[istart]);
      zi.val = &(z[istart]);

      // solve
      hazmath_solve(&A_diag[i], &ri, &zi, LU_diag[i], 0);

      // update istart
      istart = istart + A_diag[i].row;

  }

  // restore r
  array_cp(N, tempr->val, r);

  //#endif
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
 * \fn void precond_block_lower_mixed_darcy (REAL *r, REAL *z, void *data)
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

/***********************************************************************************************/
///**
// * \fn void precond_block_diag_mixed_darcy_HX(REAL *r, REAL *z, void *data)
// * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
// *        is solved inexactly by Krylov methods) (Use HX preconditioner)
// *
// * \param r     Pointer to the vector needs preconditioning
// * \param z     Pointer to preconditioned vector
// * \param data  Pointer to precondition data
// *
// * \author Xiaozhe Hu
// * \date   03/03/2019
// */
//void precond_block_diag_mixed_darcy_krylov_HX(REAL *r,
//                                              REAL *z,
//                                              void *data)
//{
//  precond_block_data *precdata=(precond_block_data *)data;
//  dvector *tempr = &(precdata->r);
//
//  block_dCSRmat *A = precdata->Abcsr;
//  AMG_param *amgparam = precdata->amgparam;
//  AMG_data **mgl = precdata->mgl;
//  dvector *el_vol = precdata->el_vol;
//  HX_div_data **hxdivdata = precdata->hxdivdata;
//
//  INT i;
//
//  const INT N0 = A->blocks[0]->row;
//  const INT N1 = A->blocks[2]->row;
//  const INT N = N0 + N1;
//
//  // back up r, setup z;
//  array_cp(N, r, tempr->val);
//  array_set(N, z, 0.0);
//
//  // prepare
//  dvector r0, r1, z0, z1;
//
//  r0.row = N0; z0.row = N0;
//  r1.row = N1; z1.row = N1;
//
//  r0.val = r; r1.val = &(r[N0]);
//  z0.val = z; z1.val = &(z[N0]);
//  //#endif
//
//  // Preconditioning A00 block (flux) using HX preconditioner
//  precond pc_flux; pc_flux.data = hxdivdata[0];
//  if (hxdivdata[0]->P_curl == NULL)
//  {
//    //pc_flux.fct = precond_hx_div_additive_2D;
//    pc_flux.fct = precond_hx_div_multiplicative_2D;
//  }
//
//  dcsr_pvfgmres(hxdivdata[0]->A, &r0, &z0, &pc_flux, 1e-3, 100, 100, 1, 1);
//
//  // Preconditioning A11 block
//  memcpy(z1.val,r1.val,N1*sizeof(REAL));
//  for (i=0;i<N1;++i) {
//    z1.val[i]/=el_vol->val[i];
//  }
//
//  // restore r
//  array_cp(N, tempr->val, r);
//
//}
//
///***********************************************************************************************/
///**
// * \fn void precond_block_lower_mixed_darcy_HX(REAL *r, REAL *z, void *data)
// * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
// *        is solved inexactly by Krylov methods) (Use HX preconditioner)
// *
// * \param r     Pointer to the vector needs preconditioning
// * \param z     Pointer to preconditioned vector
// * \param data  Pointer to precondition data
// *
// * \author Xiaozhe Hu
// * \date   03/03/2019
// */
//void precond_block_lower_mixed_darcy_krylov_HX(REAL *r,
//                                               REAL *z,
//                                               void *data)
//{
//  precond_block_data *precdata=(precond_block_data *)data;
//  dvector *tempr = &(precdata->r);
//
//  block_dCSRmat *A = precdata->Abcsr;
//  AMG_param *amgparam = precdata->amgparam;
//  AMG_data **mgl = precdata->mgl;
//  dvector *el_vol = precdata->el_vol;
//  HX_div_data **hxdivdata = precdata->hxdivdata;
//
//  INT i;
//
//  const INT N0 = A->blocks[0]->row;
//  const INT N1 = A->blocks[2]->row;
//  const INT N = N0 + N1;
//
//  // back up r, setup z;
//  array_cp(N, r, tempr->val);
//  array_set(N, z, 0.0);
//
//  // prepare
//  dvector r0, r1, z0, z1;
//
//  r0.row = N0; z0.row = N0;
//  r1.row = N1; z1.row = N1;
//
//  r0.val = r; r1.val = &(r[N0]);
//  z0.val = z; z1.val = &(z[N0]);
//  //#endif
//
//  // Preconditioning A00 block (flux)
//  precond pc_flux; pc_flux.data = hxdivdata[0];
//  if (hxdivdata[0]->P_curl == NULL)
//  {
//    //pc_flux.fct = precond_hx_div_additive_2D;
//    pc_flux.fct = precond_hx_div_multiplicative_2D;
//  }
//
//  dcsr_pvfgmres(hxdivdata[0]->A, &r0, &z0, &pc_flux, 1e-3, 100, 100, 1, 1);
//
//  // r1 = r1 - A2*z0
//  dcsr_aAxpy(-1.0, A->blocks[2], z0.val, r1.val);
//
//  // Preconditioning A11 block
//  memcpy(z1.val,r1.val,N1*sizeof(REAL));
//  for (i=0;i<N1;++i) {
//    z1.val[i]/=el_vol->val[i];
//  }
//
//  // restore r
//  array_cp(N, tempr->val, r);
//
//}
//
///***********************************************************************************************/
///**
// * \fn void precond_block_upper_mixed_darcy_HX(REAL *r, REAL *z, void *data)
// * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
// *        is solved inexactly by Krylov methods) (Use HX preconditioner)
// *
// * \param r     Pointer to the vector needs preconditioning
// * \param z     Pointer to preconditioned vector
// * \param data  Pointer to precondition data
// *
// * \author Xiaozhe Hu
// * \date   03/03/2019
// */
//void precond_block_upper_mixed_darcy_krylov_HX(REAL *r,
//                                               REAL *z,
//                                               void *data)
//{
//  precond_block_data *precdata=(precond_block_data *)data;
//  dvector *tempr = &(precdata->r);
//
//  block_dCSRmat *A = precdata->Abcsr;
//  AMG_param *amgparam = precdata->amgparam;
//  AMG_data **mgl = precdata->mgl;
//  dvector *el_vol = precdata->el_vol;
//  HX_div_data **hxdivdata = precdata->hxdivdata;
//
//  INT i;
//
//  const INT N0 = A->blocks[0]->row;
//  const INT N1 = A->blocks[2]->row;
//  const INT N = N0 + N1;
//
//  // back up r, setup z;
//  array_cp(N, r, tempr->val);
//  array_set(N, z, 0.0);
//
//  // prepare
//  dvector r0, r1, z0, z1;
//
//  r0.row = N0; z0.row = N0;
//  r1.row = N1; z1.row = N1;
//
//  r0.val = r; r1.val = &(r[N0]);
//  z0.val = z; z1.val = &(z[N0]);
//  //#endif
//
//  // Preconditioning A11 block
//  memcpy(z1.val,r1.val,N1*sizeof(REAL));
//  for (i=0;i<N1;++i) {
//    z1.val[i]/=el_vol->val[i];
//  }
//
//  // r0 = r0 - A1*z1
//  dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);
//
//  // Preconditioning A00 block (flux)
//  precond pc_flux; pc_flux.data = hxdivdata[0];
//  if (hxdivdata[0]->P_curl == NULL)
//  {
//    //pc_flux.fct = precond_hx_div_additive_2D;
//    pc_flux.fct = precond_hx_div_multiplicative_2D;
//  }
//
//  dcsr_pvfgmres(hxdivdata[0]->A, &r0, &z0, &pc_flux, 1e-3, 100, 100, 1, 1);
//
//
//  // restore r
//  array_cp(N, tempr->val, r);
//
//}


/***********************************************************************************************/
/**
 * \fn void precond_block_diag_mixed_darcy_lap (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   09/28/2017
 */
void precond_block_diag_mixed_darcy_lap(REAL *r,
                                    REAL *z,
                                    void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;
  //dvector *el_vol = precdata->el_vol;

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
  mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
  mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
  array_cp(N1, mgl[1]->x.val, z1.val);

  // restore r
  array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_diag_mixed_darcy_HX(REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly by Krylov methods) (Use HX preconditioner)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/03/2019
 */
void precond_block_diag_mixed_darcy_krylov_HX(REAL *r,
                                              REAL *z,
                                              void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  dvector *el_vol = precdata->el_vol;
  HX_div_data **hxdivdata = precdata->hxdivdata;

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

  // Preconditioning A00 block (flux) using HX preconditioner
  precond pc_flux; pc_flux.data = hxdivdata[0];
  // 2D case
  if (hxdivdata[0]->P_curl == NULL)
  {
    //pc_flux.fct = precond_hx_div_additive_2D;
    pc_flux.fct = precond_hx_div_multiplicative_2D;
  }
  // 3D case
  else
  {
    //pc_fluc.fct = precond_hx_div_additive;
    pc_flux.fct = precond_hx_div_multiplicative;
  }

  dcsr_pvfgmres(hxdivdata[0]->A, &r0, &z0, &pc_flux, 1e-2, 100, 100, 1, 1);

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
 * \fn void precond_block_lower_mixed_darcy_HX(REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly by Krylov methods) (Use HX preconditioner)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/03/2019
 */
void precond_block_lower_mixed_darcy_krylov_HX(REAL *r,
                                               REAL *z,
                                               void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  dvector *el_vol = precdata->el_vol;
  HX_div_data **hxdivdata = precdata->hxdivdata;

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
  precond pc_flux; pc_flux.data = hxdivdata[0];
  if (hxdivdata[0]->P_curl == NULL)
  {
    //pc_flux.fct = precond_hx_div_additive_2D;
    pc_flux.fct = precond_hx_div_multiplicative_2D;
  }
  // 3D case
  else
  {
    //pc_fluc.fct = precond_hx_div_additive;
    pc_flux.fct = precond_hx_div_multiplicative;
  }

  dcsr_pvfgmres(hxdivdata[0]->A, &r0, &z0, &pc_flux, 1e-2, 100, 100, 1, 1);

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
 * \fn void precond_block_upper_mixed_darcy_HX(REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly by Krylov methods) (Use HX preconditioner)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/03/2019
 */
void precond_block_upper_mixed_darcy_krylov_HX(REAL *r,
                                               REAL *z,
                                               void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  dvector *el_vol = precdata->el_vol;
  HX_div_data **hxdivdata = precdata->hxdivdata;

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
  precond pc_flux; pc_flux.data = hxdivdata[0];
  if (hxdivdata[0]->P_curl == NULL)
  {
    //pc_flux.fct = precond_hx_div_additive_2D;
    pc_flux.fct = precond_hx_div_multiplicative_2D;
  }
  // 3D case
  else
  {
    //pc_fluc.fct = precond_hx_div_additive;
    pc_flux.fct = precond_hx_div_multiplicative;
  }

  dcsr_pvfgmres(hxdivdata[0]->A, &r0, &z0, &pc_flux, 1e-2, 100, 100, 1, 1);


  // restore r
  array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_mixed_darcy_lap (REAL *r, REAL *z, void *data)
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
void precond_block_lower_mixed_darcy_lap(REAL *r,
                                    REAL *z,
                                    void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;
  //dvector *el_vol = precdata->el_vol;

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
  mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
  mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
  array_cp(N1, mgl[1]->x.val, z1.val);

  // restore r
  array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_mixed_darcy_lap (REAL *r, REAL *z, void *data)
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
void precond_block_upper_mixed_darcy_lap(REAL *r,
                                     REAL *z,
                                     void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;
  //dvector *el_vol = precdata->el_vol;

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
 * \fn void precond_block_ilu_mixed_darcy_lap (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly by Krylov methods)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   09/29/2017
 */
void precond_block_ilu_mixed_darcy_lap(REAL *r,
                                       REAL *z,
                                       void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;
  //dvector *el_vol = precdata->el_vol;

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

  // Preconditioning A11 block (head)
  mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
  mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
  array_cp(N1, mgl[1]->x.val, z1.val);

  // z1 = -z1
  dvec_ax(-1.0, &z1);

  //r0 = A00*z0 - A1*z1
  dcsr_mxv(A->blocks[0], z0.val, r0.val); // r0 = A00*z0
  dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);

  // Preconditioning A00 block again
  mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
  mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
  array_cp(N0, mgl[0]->x.val, z0.val);

  // restore r
  array_cp(N, tempr->val, r);

}


/***********************************************************************************************/
/**
 * \fn void precond_block_diag_mixed_darcy_lap_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly by Krylov methods)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   09/27/2017
 */
void precond_block_diag_mixed_darcy_lap_krylov(REAL *r,
                                           REAL *z,
                                           void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;
  //dvector *el_vol = precdata->el_vol;

  //INT i;

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

  // Preconditioning A11 block (head)
  precond_data pcdata_h;
  param_amg_to_prec(&pcdata_h,amgparam);
  pcdata_h.max_levels = mgl[1][0].num_levels;
  pcdata_h.mgl_data = mgl[1];

  precond pc_h; pc_h.data = &pcdata_h;
  pc_h.fct = precond_amg;

  dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc_h, 1e-3, 100, 100, 1, 1);

  // restore r
  array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_lower_mixed_darcy_lap_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly by Krylov methods)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   09/27/2017
 */
void precond_block_lower_mixed_darcy_lap_krylov(REAL *r,
                                            REAL *z,
                                            void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;
  //dvector *el_vol = precdata->el_vol;

  //INT i;

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

  // Preconditioning A11 block (head)
  precond_data pcdata_h;
  param_amg_to_prec(&pcdata_h,amgparam);
  pcdata_h.max_levels = mgl[1][0].num_levels;
  pcdata_h.mgl_data = mgl[1];

  precond pc_h; pc_h.data = &pcdata_h;
  pc_h.fct = precond_amg;

  dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc_h, 1e-3, 100, 100, 1, 1);

  // restore r
  array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_upper_mixed_darcy_lap_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly by Krylov methods)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   09/27/2017
 */
void precond_block_upper_mixed_darcy_lap_krylov(REAL *r,
                                            REAL *z,
                                            void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;
  //dvector *el_vol = precdata->el_vol;

  //INT i;

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

  // Preconditioning A11 block (head)
  precond_data pcdata_h;
  param_amg_to_prec(&pcdata_h,amgparam);
  pcdata_h.max_levels = mgl[1][0].num_levels;
  pcdata_h.mgl_data = mgl[1];

  precond pc_h; pc_h.data = &pcdata_h;
  pc_h.fct = precond_amg;

  dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc_h, 1e-3, 100, 100, 1, 1);

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

/***********************************************************************************************/
/**
 * \fn void precond_block_ilu_mixed_darcy_lap_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly by Krylov methods)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   09/28/2017
 */
void precond_block_ilu_mixed_darcy_lap_krylov(REAL *r,
                                            REAL *z,
                                            void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;
  //dvector *el_vol = precdata->el_vol;

  //INT i;

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

  // Preconditioning A11 block (head)
  precond_data pcdata_h;
  param_amg_to_prec(&pcdata_h,amgparam);
  pcdata_h.max_levels = mgl[1][0].num_levels;
  pcdata_h.mgl_data = mgl[1];

  precond pc_h; pc_h.data = &pcdata_h;
  pc_h.fct = precond_amg;

  dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc_h, 1e-3, 100, 100, 1, 1);

  // z1 = -z1
  dvec_ax(-1.0, &z1);

  //r0 = A00*z0 - A1*z1
  dcsr_mxv(A->blocks[0], z0.val, r0.val); // r0 = A00*z0
  dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);

  // Preconditioning A00 block again
  array_set(N0, z0.val, 0.0);
  dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc_p, 1e-3, 100, 100, 1, 1);

  // restore r
  array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block_ilu_mixed_darcy_graph_lap_krylov (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (2x2 block matrix, each diagonal block
 *        is solved inexactly by Krylov methods)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   09/28/2017
 */
void precond_block_ilu_mixed_darcy_graph_lap_krylov(REAL *r,
                                            REAL *z,
                                            void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  block_dCSRmat *A = precdata->Abcsr;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;
  dvector **diag = precdata->diag;
  //dvector *el_vol = precdata->el_vol;

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
  memcpy(z0.val,r0.val,N0*sizeof(REAL));
  for (i=0;i<N0;++i) {
      if (ABS(diag[0]->val[i])>SMALLREAL) z0.val[i]/=diag[0]->val[i];
  }

  // r1 = r1 - A2*z0
  dcsr_aAxpy(-1.0, A->blocks[2], z0.val, r1.val);

  // Preconditioning A11 block (head)
  precond_data pcdata_h;
  param_amg_to_prec(&pcdata_h,amgparam);
  pcdata_h.max_levels = mgl[1][0].num_levels;
  pcdata_h.mgl_data = mgl[1];

  precond pc_h; pc_h.data = &pcdata_h;
  pc_h.fct = precond_amg;

  dcsr_pvfgmres(&mgl[1][0].A, &r1, &z1, &pc_h, 1e-3, 100, 100, 1, 1);

  // z1 = -z1
  dvec_ax(-1.0, &z1);

  //r0 = A00*z0 - A1*z1
  dcsr_mxv(A->blocks[0], z0.val, r0.val); // r0 = A00*z0
  dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);

  // Preconditioning A00 block again
  /*
  array_set(N0, z0.val, 0.0);
  dcsr_pvfgmres(&mgl[0][0].A, &r0, &z0, &pc_p, 1e-3, 100, 100, 1, 1);
  */
  memcpy(z0.val,r0.val,N0*sizeof(REAL));
  for (i=0;i<N0;++i) {
      if (ABS(diag[0]->val[i])>SMALLREAL) z0.val[i]/=diag[0]->val[i];
  }

  // restore r
  array_cp(N, tempr->val, r);

}







/*
\note all the functions below needs to be moved -- Xiaozhe Hu
*/
/*************** Special Preconditioners for 3d-1d **********************************/
/**
 * \fn void precond_block_diag_3d1d (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning (3x3 block matrix, each diagonal block
 *        is solved inexactly)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Ana Budisa
 * \date   2020-08-24
 */
void precond_block_diag_3d1d(REAL *r,
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
  const INT N1 = A->blocks[4]->row;
  const INT N2 = A->blocks[8]->row;
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
  //#endif

  // Preconditioning A00 block (3D) - AMG on laplacian
  mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
  mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
  array_cp(N0, mgl[0]->x.val, z0.val);

  // Preconditioning A11 block (1D) - AMG on laplacian
  mgl[1]->b.row=N1; array_cp(N1, r1.val, mgl[1]->b.val); // residual is an input
  mgl[1]->x.row=N1; dvec_set(N1, &mgl[1]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[1], amgparam);
  array_cp(N1, mgl[1]->x.val, z1.val);

  // Preconditioning A11 block (multiplier) - FAMG SUM Fracs
  precond_sum_famg_add2(r2.val, z2.val, precdata);

  // restore r
  array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block2_babuska_diag (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning for a standard babuska problem
 *        (2x2 block matrix, first block is solved by AMG, second block is solved by
 *        rational approximation)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Ana Budisa
 * \date   2020-10-19
 *
 * \note modified by Xiaozhe on 12/26/2020
 */
void precond_block2_babuska_diag(REAL *r,
                                 REAL *z,
                                 void *data)
{
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dvector *tempr = &(precdata->r);
    INT status;// = SUCCESS;

    const INT N0 = A->blocks[0]->row;
    const INT N1 = A->blocks[3]->row;
    const INT N = N0 + N1;

    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);

    // prepare
    INT i; //not used:, j;
    AMG_param *amgparam = precdata->amgparam; // array!!
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, z0, z1;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;

    r0.val = r; r1.val = &(r[N0]);
    z0.val = z; z1.val = &(z[N0]);

    //--------------------------------------------------------
    // Step 1: Preconditioning A00 block
    //--------------------------------------------------------
    // mgl[0] is pointer to mgl data for the first diagonal block
    /*
    // apply AMG only to the first diagonal block
    {
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], &(amgparam[0]));
    array_cp(N0, mgl[0]->x.val, z0.val);
    }*/

    // apply AMG + Krylov to the first diagonal block
    {
      precond pc00;
      pc00.fct = precond_amg;
      precond_data pcdata00;
      param_amg_to_prec(&pcdata00, &(amgparam[0]));
      pc00.data = &pcdata00;

      pcdata00.max_levels = mgl[0]->num_levels;
      pcdata00.mgl_data = mgl[0];

      // solve
      //status = dcsr_pvfgmres(&(mgl[0][0].A), &r0, &z0, &pc00, 1e-6, 100, 100, 1, 1);
      status = dcsr_pcg(&(mgl[0][0].A), &r0, &z0, &pc00, 1e-12, 100, 1, 1);
      if(status<SUCCESS)
	WARN_STATUS(__FUNCTION__,"dcsr_pcg(...)",status);
    }

    // direct solve the first diagonal block
    //directsolve_HAZ(&(mgl[0][0].A), &r0, &z0, 1);
    //--------------------------------------------------------

    //--------------------------------------------------------
    // Step 2: Preconditioning A11 block
    //--------------------------------------------------------

    /*----------------------------------------*/
    // get scaled mass matrix
    dCSRmat *scaled_M = precdata->scaled_M;
    dvector *diag_scaled_M = precdata->diag_scaled_M;

    // get scaled alpha and beta
    REAL scaled_alpha = precdata->scaled_alpha;
    REAL scaled_beta  = precdata->scaled_beta;

    // get poles and residues
    //not used:    dvector *poles = precdata->poles;
    dvector *residues = precdata->residues;

    INT k = residues->row;
    /*----------------------------------------*/

    /*----------------------------------------*/
    /* set up preconditioners */
    /*----------------------------------------*/
    // pc for scaled mass matrix
    precond pc_scaled_M;
    pc_scaled_M.data = diag_scaled_M;
    pc_scaled_M.fct  = precond_diag;

    // pc for krylov for shifted Laplacians
    precond pc_frac_A;
    pc_frac_A.fct = precond_amg;
    precond_data pcdata;
    param_amg_to_prec(&pcdata, &(amgparam[1]));
    pc_frac_A.data = &pcdata;
    /*----------------------------------------*/


    /*----------------------------------------*/
    /* main loop of applying rational approximation */
    /*----------------------------------------*/
    // scaling r1
    if (scaled_alpha > scaled_beta)
    {
      dvec_ax(1./scaled_alpha, &r1);
    }
    else
    {
      dvec_ax(1./scaled_beta, &r1);
    }

    // z1 = residues(0)*(scaled_M\scaled_r1)
    status = dcsr_pcg(scaled_M, &r1, &z1, &pc_scaled_M, 1e-12, 100, 1, 1);


    dvector update = dvec_create(N1);
    array_ax(N1, residues->val[0], z1.val);

    for(i = 1; i < k; ++i) {

        mgl[i]->b.row = N1; array_cp(N1, r1.val, mgl[i]->b.val); // residual is an input
        mgl[i]->x.row = N1; dvec_set(N1, &mgl[i]->x, 0.0);

        // set precond data and param
        pcdata.max_levels = mgl[i]->num_levels;
        pcdata.mgl_data = mgl[i];

        // set update to zero
        dvec_set(update.row, &update, 0.0);

        // solve
        //status = dcsr_pvfgmres(&(mgl[i][0].A), &r1, &update, &pc_frac_A, 1e-6, 100, 100, 1, 1);
        status = dcsr_pcg(&(mgl[i][0].A), &r1, &update, &pc_frac_A, 1e-12, 100, 1, 1);

        // z = z + residues[i+1]*update
        array_axpy(N1, residues->val[i], update.val, z1.val);

    }

    // direct solve A11 block
    //directsolve_HAZ(As, &r1, &z1, 1);
    //--------------------------------------------------------

    // restore r
    array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block2_babuska_lower (REAL *r, REAL *z, void *data)
 * \brief block lower triangular preconditioning for a standard babuska problem
 *        (2x2 block matrix, first block is solved by AMG, second block is solved by
 *        rational approximation)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Ana Budisa
 * \date   2020-10-19
 *
 * \note modified by Xiaozhe on 12/27/2020
 */
void precond_block2_babuska_lower(REAL *r,
                                  REAL *z,
                                  void *data)
{
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dvector *tempr = &(precdata->r);
    INT status;// = SUCCESS;

    const INT N0 = A->blocks[0]->row;
    const INT N1 = A->blocks[3]->row;
    const INT N = N0 + N1;

    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);

    // prepare
    INT i;// not used:, j;
    AMG_param *amgparam = precdata->amgparam; // array!!
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, z0, z1;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;

    r0.val = r; r1.val = &(r[N0]);
    z0.val = z; z1.val = &(z[N0]);

    //--------------------------------------------------------
    // Step 1: Preconditioning A00 block
    //--------------------------------------------------------
    // mgl[0] is pointer to mgl data for the first diagonal block
    /*
    // apply AMG only to the first diagonal block
    {
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], &(amgparam[0]));
    array_cp(N0, mgl[0]->x.val, z0.val);
    }*/

    // apply AMG + Krylov to the first diagonal block
    {
      precond pc00;
      pc00.fct = precond_amg;
      precond_data pcdata00;
      param_amg_to_prec(&pcdata00, &(amgparam[0]));
      pc00.data = &pcdata00;

      pcdata00.max_levels = mgl[0]->num_levels;
      pcdata00.mgl_data = mgl[0];

      // solve
      //status = dcsr_pvfgmres(&(mgl[0][0].A), &r0, &z0, &pc00, 1e-6, 100, 100, 1, 1);
      status = dcsr_pcg(&(mgl[0][0].A), &r0, &z0, &pc00, 1e-12, 100, 1, 1);
      if(status<SUCCESS)
	WARN_STATUS(__FUNCTION__,"dcsr_pcg(...)",status);
    }

    // direct solve the first diagonal block
    //directsolve_HAZ(&(mgl[0][0].A), &r0, &z0, 1);
    //--------------------------------------------------------

    //--------------------------------------------------------
    // Step 2: update r1
    //--------------------------------------------------------
    // r1 = r1 - A2*z0
    dcsr_aAxpy(-1.0, A->blocks[2], z0.val, r1.val);

    //--------------------------------------------------------
    // Step 3: Preconditioning A11 block
    //--------------------------------------------------------

    /*----------------------------------------*/
    // get scaled mass matrix
    dCSRmat *scaled_M = precdata->scaled_M;
    dvector *diag_scaled_M = precdata->diag_scaled_M;

    // get scaled alpha and beta
    REAL scaled_alpha = precdata->scaled_alpha;
    REAL scaled_beta  = precdata->scaled_beta;

    // get poles and residues
    //not used:    dvector *poles = precdata->poles;
    dvector *residues = precdata->residues;

    INT k = residues->row;
    /*----------------------------------------*/

    /*----------------------------------------*/
    /* set up preconditioners */
    /*----------------------------------------*/
    // pc for scaled mass matrix
    precond pc_scaled_M;
    pc_scaled_M.data = diag_scaled_M;
    pc_scaled_M.fct  = precond_diag;

    // pc for krylov for shifted Laplacians
    precond pc_frac_A;
    pc_frac_A.fct = precond_amg;
    precond_data pcdata;
    param_amg_to_prec(&pcdata, &(amgparam[1]));
    pc_frac_A.data = &pcdata;
    /*----------------------------------------*/


    /*----------------------------------------*/
    /* main loop of applying rational approximation */
    /*----------------------------------------*/
    // scaling r1
    if (scaled_alpha > scaled_beta)
    {
      dvec_ax(1./scaled_alpha, &r1);
    }
    else
    {
      dvec_ax(1./scaled_beta, &r1);
    }

    // z1 = residues(0)*(scaled_M\scaled_r1)
    status = dcsr_pcg(scaled_M, &r1, &z1, &pc_scaled_M, 1e-12, 100, 1, 1);


    dvector update = dvec_create(N1);
    array_ax(N1, residues->val[0], z1.val);

    for(i = 1; i < k; ++i) {

        mgl[i]->b.row = N1; array_cp(N1, r1.val, mgl[i]->b.val); // residual is an input
        mgl[i]->x.row = N1; dvec_set(N1, &mgl[i]->x, 0.0);

        // set precond data and param
        pcdata.max_levels = mgl[i]->num_levels;
        pcdata.mgl_data = mgl[i];

        // set update to zero
        dvec_set(update.row, &update, 0.0);

        // solve
        //status = dcsr_pvfgmres(&(mgl[i][0].A), &r1, &update, &pc_frac_A, 1e-6, 100, 100, 1, 1);
        status = dcsr_pcg(&(mgl[i][0].A), &r1, &update, &pc_frac_A, 1e-12, 100, 1, 1);

        // z = z + residues[i+1]*update
        array_axpy(N1, residues->val[i], update.val, z1.val);

    }

    // direct solve A11 block
    //directsolve_HAZ(As, &r1, &z1, 1);
    //--------------------------------------------------------

    // restore r
    array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_block2_babuska_upper(REAL *r, REAL *z, void *data)
 * \brief block preconditioning for a standard babuska problem
 *        (2x2 block matrix, first block is solved by AMG, second block is solved by
 *        rational approximation)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   2020-12-27
 *
 */
void precond_block2_babuska_upper(REAL *r,
                                  REAL *z,
                                  void *data)
{
    precond_block_data *precdata=(precond_block_data *)data;
    block_dCSRmat *A = precdata->Abcsr;
    dvector *tempr = &(precdata->r);
    INT status;// = SUCCESS;

    const INT N0 = A->blocks[0]->row;
    const INT N1 = A->blocks[3]->row;
    const INT N = N0 + N1;

    // back up r, setup z;
    array_cp(N, r, tempr->val);
    array_set(N, z, 0.0);

    // prepare
    INT i;//not used:, j;
    AMG_param *amgparam = precdata->amgparam; // array!!
    AMG_data **mgl = precdata->mgl;
    dvector r0, r1, z0, z1;

    r0.row = N0; z0.row = N0;
    r1.row = N1; z1.row = N1;

    r0.val = r; r1.val = &(r[N0]);
    z0.val = z; z1.val = &(z[N0]);

    //--------------------------------------------------------
    // Step 1: Preconditioning A11 block
    //--------------------------------------------------------

    /*----------------------------------------*/
    // get scaled mass matrix
    dCSRmat *scaled_M = precdata->scaled_M;
    dvector *diag_scaled_M = precdata->diag_scaled_M;

    // get scaled alpha and beta
    REAL scaled_alpha = precdata->scaled_alpha;
    REAL scaled_beta  = precdata->scaled_beta;

    // get poles and residues
    //    dvector *poles = precdata->poles;
    dvector *residues = precdata->residues;

    INT k = residues->row;
    /*----------------------------------------*/

    /*----------------------------------------*/
    /* set up preconditioners */
    /*----------------------------------------*/
    // pc for scaled mass matrix
    precond pc_scaled_M;
    pc_scaled_M.data = diag_scaled_M;
    pc_scaled_M.fct  = precond_diag;

    // pc for krylov for shifted Laplacians
    precond pc_frac_A;
    pc_frac_A.fct = precond_amg;
    precond_data pcdata;
    param_amg_to_prec(&pcdata, &(amgparam[1]));
    pc_frac_A.data = &pcdata;
    /*----------------------------------------*/

    /*----------------------------------------*/
    /* main loop of applying rational approximation */
    /*----------------------------------------*/
    // scaling r1
    if (scaled_alpha > scaled_beta)
    {
      dvec_ax(1./scaled_alpha, &r1);
    }
    else
    {
      dvec_ax(1./scaled_beta, &r1);
    }

    // z1 = residues(0)*(scaled_M\scaled_r1)
    status = dcsr_pcg(scaled_M, &r1, &z1, &pc_scaled_M, 1e-6, 100, 1, 1);
    if(status<SUCCESS)
      WARN_STATUS(__FUNCTION__,"dcsr_pcg(scaled_M,...)",status);

    dvector update = dvec_create(N1);
    array_ax(N1, residues->val[0], z1.val);

    for(i = 1; i < k; ++i) {

        mgl[i]->b.row = N1; array_cp(N1, r1.val, mgl[i]->b.val); // residual is an input
        mgl[i]->x.row = N1; dvec_set(N1, &mgl[i]->x, 0.0);

        // set precond data and param
        pcdata.max_levels = mgl[i]->num_levels;
        pcdata.mgl_data = mgl[i];

        // set update to zero
        dvec_set(update.row, &update, 0.0);

        // solve
        //status = dcsr_pvfgmres(&(mgl[i][0].A), &r1, &update, &pc_frac_A, 1e-6, 100, 100, 1, 1);
        status = dcsr_pcg(&(mgl[i][0].A), &r1, &update, &pc_frac_A, 1e-6, 100, 1, 1);
	if(status<SUCCESS)
	  WARN_STATUS(__FUNCTION__,"dcsr_pcg(...)",status);

        // z = z + residues[i+1]*update
        array_axpy(N1, residues->val[i], update.val, z1.val);

    }

    // direct solve A11 block
    //directsolve_HAZ(As, &r1, &z1, 1);
    //--------------------------------------------------------

    //--------------------------------------------------------
    // Step 2: update r0
    //--------------------------------------------------------
    // r0 = r0 - A1*z1
    dcsr_aAxpy(-1.0, A->blocks[1], z1.val, r0.val);

    //--------------------------------------------------------
    // Step 3: Preconditioning A00 block
    //--------------------------------------------------------
    // mgl[0] is pointer to mgl data for the first diagonal block
    /*
    // apply AMG only to the first diagonal block
    {
    mgl[0]->b.row=N0; array_cp(N0, r0.val, mgl[0]->b.val); // residual is an input
    mgl[0]->x.row=N0; dvec_set(N0, &mgl[0]->x, 0.0);

    for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], &(amgparam[0]));
    array_cp(N0, mgl[0]->x.val, z0.val);
    }*/

    // apply AMG + Krylov to the first diagonal block
    {
    precond pc00;
    pc00.fct = precond_amg;
    precond_data pcdata00;
    param_amg_to_prec(&pcdata00, &(amgparam[0]));
    pc00.data = &pcdata00;

    pcdata00.max_levels = mgl[0]->num_levels;
    pcdata00.mgl_data = mgl[0];

    // solve
    //status = dcsr_pvfgmres(&(mgl[0][0].A), &r0, &z0, &pc00, 1e-6, 100, 100, 1, 1);
    status = dcsr_pcg(&(mgl[0][0].A), &r0, &z0, &pc00, 1e-6, 100, 1, 1);
    if(status<SUCCESS)
      WARN_STATUS(__FUNCTION__,"dcsr_pcg(...)",status);
    }

    // direct solve the first diagonal block
    //directsolve_HAZ(&(mgl[0][0].A), &r0, &z0, 1);
    //--------------------------------------------------------



    // restore r
    array_cp(N, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_ra_fenics (REAL *r, REAL *z, void *data)
 * \brief preconditioning a fractional problem with rational approximation solver
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precond_ra_data
 *
 * \author Ana Budisa
 * \date   2021-02-03
 *
 * \note 2021-10-28 updated to include complex valued poles an residues
 *                  by using algorithm Ludmil provided;
 * \note 2021-10-28 here we solve for both pole=a+ib and pole=a-ib even though
 *                  the solving process is the same; FIXME
 */
void precond_ra_fenics(REAL *r, REAL *z, void *data)
{
    // fprintf(stdout, " \n ---- Call from python ---- \n"); fflush(stdout);
    // local variables
  INT status;// = SUCCESS;
    precond_ra_data *precdata=(precond_ra_data *)data;
    AMG_data **mgl = precdata->mgl; // count from 0
    AMG_param *amgparam = precdata->amgparam; // this is not an array anymore

    INT n = precdata->scaled_A->col; // general size of the problem
    INT i;

    // back up r, setup z
    /* array_set(n, z, 0.0); */
    /* dvector z_vec = dvec_create(n); memcpy(z_vec.val, z, n*sizeof(REAL)); */
    dvector z_vec, r_vec;
    z_vec.row = n;
    z_vec.val = z;
    array_set(n, z_vec.val, 0e0);
    r_vec.row = n;
    r_vec.val = r;
    // printf("Norm residual (n=%lld)before copy: %.5e\n", n, array_norm2(n, r));
    /* dvector r_vec = dvec_create(n); // we probably can just have r_vec.val=r; */
    /* fprintf(stderr,"\nr_vec %ld: %lld\n",sizeof(r_vec.val)/sizeof(REAL),r_vec.row); */
    /* if(r_vec.val == NULL){ */
    /*   fprintf(stderr,"\nCOULD NOT ALLOCATE r_vec\n"); */
    /*   fflush(stderr); */
    /* } */
    /* memcpy(r_vec.val, r, n*sizeof(REAL)); */
    /*----------------------------------------*/
    // get scaled mass matrix
    dCSRmat *scaled_M = precdata->scaled_M;
    dvector *diag_scaled_M = precdata->diag_scaled_M;

    // get scaled alpha and beta
    REAL scaled_alpha = precdata->scaled_alpha;
    REAL scaled_beta  = precdata->scaled_beta;

    // get poles and residues
    dvector *poles = precdata->poles;
    dvector *residues = precdata->residues;

    // number of poles
    INT npoles = (INT)poles->row/2; // because we have imag parts appended after real parts

    /*----------------------------------------*/

    /*----------------------------------------*/
    /* set up preconditioners */
    /*----------------------------------------*/
    // pc for scaled mass matrix
    precond pc_scaled_M;
    pc_scaled_M.data = diag_scaled_M;
    pc_scaled_M.fct  = precond_diag;

    // pc for krylov for shifted Laplacians
    precond pc_frac_A;
    pc_frac_A.fct = precond_amg;
    precond_data pcdata;
    param_amg_to_prec(&pcdata, amgparam);
    pc_frac_A.data = &pcdata;

    /*----------------------------------------*/

    /*----------------------------------------*/
    /* main loop of applying rational approximation */
    /*----------------------------------------*/
    // scaling r
    if (scaled_alpha > scaled_beta)
    {
      dvec_ax(1./scaled_alpha, &r_vec);
    }
    else
    {
      dvec_ax(1./scaled_beta, &r_vec);
    }

    // z = residues(0)*(scaled_M\scaled_r1)
    // NOTE: here we assume imag(residues(0)) = 0
    if(fabs(residues->val[0]) > 0.) {
        status = dcsr_pcg(scaled_M, &r_vec, &z_vec, &pc_scaled_M, 1e-6, 100, 1, 0);
	if(status<SUCCESS)
	  WARN_STATUS(__FUNCTION__,"dcsr_pcg(...)",status);
        array_ax(n, residues->val[0], z_vec.val);
    }

    dvector update = dvec_create(n);
    /* dvector u000 = dvec_create(n); */
    // INT solver_flag,jjj;
    /* printf("\nNumber of poles: %lld\n", npoles); */
    // INT count = 0;

    for(i = 0; i < npoles; ++i) {

        // fprintf(stdout, "We are at pole %lld with value %.10f + i %.10f \n", i, poles->val[i], poles->val[i+npoles]); fflush(stdout);
        if(fabs(poles->val[i+npoles]) > 0.) {
            // then we have a nonzero imag part of that pole and we do the 2x2 block algorithm
            /* solve
            [(A - Re(pole)*I),       Im(pole)*I ] [upd ] = [r]
            [   - Im(pole)*I ,  (A - Re(pole)*I)] [iupd] = [0]
            */
            REAL p_im = poles->val[i+npoles];
            dvector iupdate = dvec_create(n); // imag part of the update
            INT K = 5; // number of loops for this algorithm
            INT k;

            // set initial updates to zero
            dvec_set(update.row, &update, 0.0);
            dvec_set(iupdate.row, &iupdate, 0.0);

            // create right hand sides for inv(A - dI);
            // rhs_1 = r - Im(pole) * iupdate; rhs_2 = Im(pole) * update
            dvector rhs1 = dvec_create(n);
            dvector rhs2 = dvec_create(n);

            // amg data for this pole
            pcdata.max_levels = mgl[i]->num_levels;
            pcdata.mgl_data = mgl[i];

            for(k = 0; k < K; ++k) {
                // set right hand sides
                dvec_set(rhs1.row, &rhs1, 0.0);
                dvec_set(rhs2.row, &rhs2, 0.0);
                dvec_axpyz(-p_im, &iupdate, &r_vec, &rhs1);
                dvec_axpy(p_im, &update, &rhs2);

                // (1) solve (A - Re(pole)*I) update = rhs1
                // set amg data
                mgl[i]->b.row = n; array_cp(n, rhs1.val, mgl[i]->b.val);
                mgl[i]->x.row = n; dvec_set(n, &mgl[i]->x, 0.0);

                status = dcsr_pcg(&(mgl[i][0].A), &rhs1, &update, &pc_frac_A, 1e-6, 100, 1, 0);
		        if(status<SUCCESS)
		            WARN_STATUS(__FUNCTION__,"dcsr_pcg((1) solve...)",status);
                // (2) solve (A - Re(pole)*I) iupdate = rhs2
                // set amg data
                mgl[i]->b.row = n; array_cp(n, rhs2.val, mgl[i]->b.val);
                mgl[i]->x.row = n; dvec_set(n, &mgl[i]->x, 0.0);

                status = dcsr_pcg(&(mgl[i][0].A), &rhs2, &iupdate, &pc_frac_A, 1e-6, 100, 1, 0);
		        if(status<SUCCESS)
		            WARN_STATUS(__FUNCTION__,"dcsr_pcg((2) solve...)",status);
            }

            // update next increment
            // first check if Im(residue) > 0
            if(fabs(residues->val[(npoles+1)+i+1]) > 0.) {
                // z = z + residues[i+1]*update - residues[npoles+1+i+1]*iupdate
                array_axpy(n, 2*residues->val[i+1], update.val, z_vec.val);
                array_axpy(n, -2*residues->val[(npoles+1)+i+1], iupdate.val, z_vec.val);
            }
            else {
                // z = z + residues[i+1]*update
                array_axpy(n, 2*residues->val[i+1], update.val, z_vec.val);
            }

            // free memory
            dvec_free(&iupdate);
            dvec_free(&rhs1);
            dvec_free(&rhs2);

        }
        else{
            // else we do the standard algorithm to solve (D - dI) * update = r
            mgl[i]->b.row = n; array_cp(n, r_vec.val, mgl[i]->b.val); // residual is an input
            mgl[i]->x.row = n; dvec_set(n, &mgl[i]->x, 0.0);

            // set precond data and param
            pcdata.max_levels = mgl[i]->num_levels;
            pcdata.mgl_data = mgl[i];

            // set update to zero
            dvec_set(update.row, &update, 0.0);

            // solve
            // printf("\tPole %lld, norm of r = %e\n", i, dvec_norm2(&r_vec));
            // status = dcsr_pvfgmres(&(mgl[i][0].A), &r1, &update, &pc_frac_A, 1e-6, 100, 100, 1, 1);
            status = dcsr_pcg(&(mgl[i][0].A), &r_vec, &update, &pc_frac_A, 1e-6, 100, 1, 0);
	    if(status<SUCCESS)
	      WARN_STATUS(__FUNCTION__,"dcsr_pcg(...)",status);
            /* void *numeric=NULL;  // prepare for direct solve.
            numeric=factorize_HAZ(&(mgl[i][0].A),0);
            solver_flag=(INT )solve_HAZ(&(mgl[i][0].A),	\
                              &r_vec,		\
                              &update,		\
                              numeric,
                              0);
            free(numeric); */
            // if(status > 0) count += status;
            // printf("\tPole %lld, norm of update = %e\n", i, dvec_norm2(&update));
            // z = z + residues[i+1]*update
            array_axpy(n, residues->val[i+1], update.val, z_vec.val);
            /* for(jjj=0;jjj<z_vec.row;jjj++) z[jjj]=z_vec.val[jjj]; */
        }
    }

    // if(count) printf("Inner solver took total of %lld iterations. \n", count);

    // cleanup
    // UNSCALLING r
    if (scaled_alpha > scaled_beta) {
      dvec_ax(scaled_alpha, &r_vec);
    }
    else {
      dvec_ax(scaled_beta, &r_vec);
    }
    dvec_free(&update);
    //    dvec_free(&r_vec);
    /* dvec_free(&z_vec); */
    return;
}

/***********************************************************************************************/
/**
 * \fn void precond_bdcsr_metric_amg_exact(REAL *r, REAL *z, void *data)
 *
 * \brief metric AMG preconditioner
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/16/2022
 *
 * \note direct solver on interface + AMG
 *
 */
void precond_bdcsr_metric_amg_exact(REAL *r,
                       REAL *z,
                       void *data)
{
    precond_data_bdcsr *predata=(precond_data_bdcsr *)data;
    // const INT brow=predata->mgl_data[0].A.brow;
    // const INT bcol=predata->mgl_data[0].A.bcol;
    const INT maxit=predata->maxit;
    const INT total_row = predata->total_row;
    const INT total_col = predata->total_col;
    dvector *tempr = &(predata->r);

    // back up r, setup z;
    array_cp(total_row, r, tempr->val);
    array_set(total_row, z, 0.0);

    // local variables
    INT i;

    // Schwarz method on the interface part
    // Schwarz_param *schwarz_param = predata->schwarz_param;
    //#if WITH_SUITESPARSE
    dvector rr, zz;

    rr.row = predata->A->blocks[3]->row; rr.val = r+predata->A->blocks[0]->row;
    zz.row = predata->A->blocks[3]->col; zz.val = z+predata->A->blocks[0]->row;
    Schwarz_data *schwarz_data = predata->schwarz_data;
    void **LU_data = predata->LU_data;
    hazmath_solve(&schwarz_data->A, &rr, &zz, LU_data[0], 0);
    //directsolve_HAZ(&schwarz_data->A, &rr, &zz, 1);
    //    #else
    //    error_extlib(257, __FUNCTION__, "SuiteSparse");
    //#endif


    // AMG solve on the whole matrix
    AMG_param amgparam; param_amg_init(&amgparam);
    amgparam.cycle_type = predata->cycle_type;
    amgparam.smoother   = predata->smoother;
    amgparam.presmooth_iter  = predata->presmooth_iter;
    amgparam.postsmooth_iter = predata->postsmooth_iter;
    amgparam.relaxation = predata->relaxation;
    amgparam.coarse_scaling = predata->coarse_scaling;
    amgparam.tentative_smooth = predata->tentative_smooth;

    AMG_data_bdcsr *mgl = predata->mgl_data;
    mgl->b.row=total_row; array_cp(total_row, r, mgl->b.val); // residual is an input
    mgl->x.row=total_col; array_cp(total_row, z, mgl->x.val); //dvec_set(total_col, &mgl->x, 0.0);

    for ( i=maxit; i--; ) mgcycle_bdcsr(mgl,&amgparam);

    array_cp(total_col, mgl->x.val, z);

    // restore residual
    array_cp(total_row, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_bdcsr_metric_amg_exact_additive(REAL *r, REAL *z, void *data)
 *
 * \brief metric AMG preconditioner (additive version)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/16/2022
 *
 * \note Additive version: direct solver on interface + AMG
 *
 */
void precond_bdcsr_metric_amg_exact_additive(REAL *r,
                       REAL *z,
                       void *data)
{
    precond_data_bdcsr *predata=(precond_data_bdcsr *)data;
    // const INT brow=predata->mgl_data[0].A.brow;
    // const INT bcol=predata->mgl_data[0].A.bcol;
    const INT maxit=predata->maxit;
    const INT total_row = predata->total_row;
    const INT total_col = predata->total_col;
    dvector *tempr = &(predata->r);

    // back up r, setup z;
    array_cp(total_row, r, tempr->val);
    array_set(total_row, z, 0.0);

    // local variables
	INT i;

    //smoother_dcsr_Schwarz_forward(schwarz_data, schwarz_param, &zz, &rr);
    //smoother_dcsr_Schwarz_backward(schwarz_data, schwarz_param, &zz, &rr);
	//#if WITH_SUITESPARSE
    dvector rr, zz;

    rr.row = predata->A->blocks[3]->row; rr.val = r+predata->A->blocks[0]->row;
    zz.row = predata->A->blocks[3]->col; zz.val = z+predata->A->blocks[0]->row;

    // Schwarz method on the interface part
    // Schwarz_param *schwarz_param = predata->schwarz_param;
    Schwarz_data *schwarz_data = predata->schwarz_data;
    void **LU_data = predata->LU_data;
    hazmath_solve(&schwarz_data->A, &rr, &zz, LU_data[0], 0);
    //directsolve_HAZ(&schwarz_data->A, &rr, &zz, 1);
    //#else
    //    error_extlib(257, __FUNCTION__, "SuiteSparse");
    //#endif

    // update residual
    //bdcsr_aAxpy(-1.0, predata->A, z, r);

    // AMG solve on the whole matrix
    AMG_param amgparam; param_amg_init(&amgparam);
    amgparam.cycle_type = predata->cycle_type;
    amgparam.smoother   = predata->smoother;
    amgparam.presmooth_iter  = predata->presmooth_iter;
    amgparam.postsmooth_iter = predata->postsmooth_iter;
    amgparam.relaxation = predata->relaxation;
    amgparam.coarse_scaling = predata->coarse_scaling;
    amgparam.tentative_smooth = predata->tentative_smooth;

    AMG_data_bdcsr *mgl = predata->mgl_data;
    mgl->b.row=total_row; array_cp(total_row, r, mgl->b.val); // residual is an input
    mgl->x.row=total_col; dvec_set(total_col, &mgl->x, 0.0);

    for ( i=maxit; i--; ) mgcycle_bdcsr(mgl,&amgparam);

    // update solution (additive)
    array_axpy(total_col, 1.0, mgl->x.val, z);
    //array_cp(total_col, mgl->x.val, z);

    // restore residual
    array_cp(total_row, tempr->val, r);

}


/***********************************************************************************************/
/**
 * \fn void precond_bdcsr_metric_amg(REAL *r, REAL *z, void *data)
 *
 * \brief metric AMG preconditioner (nonsymmeric multiplicative version)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/16/2022
 *
 * \note symmetric Schwarz on interface + AMG (overall multiplicaive and nonsymmetric)
 *
 */
void precond_bdcsr_metric_amg(REAL *r,
                       REAL *z,
                       void *data)
{
    precond_data_bdcsr *predata=(precond_data_bdcsr *)data;
    // const INT brow=predata->mgl_data[0].A.brow;
    // const INT bcol=predata->mgl_data[0].A.bcol;
    const INT maxit=predata->maxit;
    const INT total_row = predata->total_row;
    const INT total_col = predata->total_col;
    dvector *tempr = &(predata->r);

    // back up r, setup z;
    array_cp(total_row, r, tempr->val);
    array_set(total_row, z, 0.0);

    // local variables
	INT i;
    dvector rr, zz;

    rr.row = predata->A->blocks[3]->row; rr.val = r+predata->A->blocks[0]->row;
    zz.row = predata->A->blocks[3]->col; zz.val = z+predata->A->blocks[0]->row;

    // Schwarz method on the interface part
    Schwarz_param *schwarz_param = predata->schwarz_param;
    Schwarz_data *schwarz_data = predata->schwarz_data;
    smoother_dcsr_Schwarz_forward(schwarz_data, schwarz_param, &zz, &rr);
    smoother_dcsr_Schwarz_backward(schwarz_data, schwarz_param, &zz, &rr);
    //directsolve_HAZ(&schwarz_data->A, &rr, &zz, 1);

    // AMG solve on the whole matrix
    AMG_param amgparam; param_amg_init(&amgparam);
    amgparam.cycle_type = predata->cycle_type;
    amgparam.smoother   = predata->smoother;
    amgparam.presmooth_iter  = predata->presmooth_iter;
    amgparam.postsmooth_iter = predata->postsmooth_iter;
    amgparam.relaxation = predata->relaxation;
    amgparam.coarse_scaling = predata->coarse_scaling;
    amgparam.tentative_smooth = predata->tentative_smooth;

    AMG_data_bdcsr *mgl = predata->mgl_data;
    mgl->b.row=total_row; array_cp(total_row, r, mgl->b.val); // residual is an input
    mgl->x.row=total_col; array_cp(total_row, z, mgl->x.val); //dvec_set(total_col, &mgl->x, 0.0);

    for ( i=maxit; i--; ) mgcycle_bdcsr(mgl,&amgparam);

    array_cp(total_col, mgl->x.val, z);

    // restore residual
    array_cp(total_row, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_bdcsr_metric_amg_additive(REAL *r, REAL *z, void *data)
 *
 * \brief metric AMG preconditioner (additive version)
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/16/2022
 *
 * \note symmetric Schwarz on interface + AMG (overall additive version)
 *
 */
void precond_bdcsr_metric_amg_additive(REAL *r,
                       REAL *z,
                       void *data)
{
    precond_data_bdcsr *predata=(precond_data_bdcsr *)data;
    // const INT brow=predata->mgl_data[0].A.brow;
    // const INT bcol=predata->mgl_data[0].A.bcol;
    const INT maxit=predata->maxit;
    const INT total_row = predata->total_row;
    const INT total_col = predata->total_col;
    dvector *tempr = &(predata->r);

    // back up r, setup z;
    array_cp(total_row, r, tempr->val);
    array_set(total_row, z, 0.0);

    // local variables
	INT i;
    dvector rr, zz;

    rr.row = predata->A->blocks[3]->row; rr.val = r+predata->A->blocks[0]->row;
    zz.row = predata->A->blocks[3]->col; zz.val = z+predata->A->blocks[0]->row;

    // Schwarz method on the interface part
    Schwarz_param *schwarz_param = predata->schwarz_param;
    Schwarz_data *schwarz_data = predata->schwarz_data;
    smoother_dcsr_Schwarz_forward(schwarz_data, schwarz_param, &zz, &rr);
    smoother_dcsr_Schwarz_backward(schwarz_data, schwarz_param, &zz, &rr);
    //directsolve_HAZ(&schwarz_data->A, &rr, &zz, 1);

    // AMG solve on the whole matrix
    AMG_param amgparam; param_amg_init(&amgparam);
    amgparam.cycle_type = predata->cycle_type;
    amgparam.smoother   = predata->smoother;
    amgparam.presmooth_iter  = predata->presmooth_iter;
    amgparam.postsmooth_iter = predata->postsmooth_iter;
    amgparam.relaxation = predata->relaxation;
    amgparam.coarse_scaling = predata->coarse_scaling;
    amgparam.tentative_smooth = predata->tentative_smooth;

    AMG_data_bdcsr *mgl = predata->mgl_data;
    mgl->b.row=total_row; array_cp(total_row, r, mgl->b.val); // residual is an input
    mgl->x.row=total_col; dvec_set(total_col, &mgl->x, 0.0);
    // array_cp(total_row, z, mgl->x.val);

    for ( i=maxit; i--; ) mgcycle_bdcsr(mgl,&amgparam);

    // update solution
    array_axpy(total_col, 1.0, mgl->x.val, z);

    // restore residual
    array_cp(total_row, tempr->val, r);

}

/***********************************************************************************************/
/**
 * \fn void precond_bdcsr_metric_amg_symmetric(REAL *r, REAL *z, void *data)
 *
 * \brief metric AMG preconditioner
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/16/2022
 *
 * \note symmetric Schwarz on interface + AMG
 *
 * \note edited 2022-11-11 by ana (add permutation)
 */
void precond_bdcsr_metric_amg_symmetric(REAL *r,
                       REAL *z,
                       void *data)
{
    precond_data_bdcsr *predata=(precond_data_bdcsr *)data;
    // const INT brow=predata->mgl_data[0].A.brow;
    // const INT bcol=predata->mgl_data[0].A.bcol;
    const INT maxit=predata->maxit;
    const INT total_row = predata->total_row;
    const INT total_col = predata->total_col;
    dvector *tempr = &(predata->r);

    // back up r, setup z;
    array_cp(total_row, r, tempr->val);
    array_set(total_row, z, 0.0);

    // local variables
	INT i;
    dvector rr, zz;

    // permute residual (IN PLACE, USING THE BACKUP!)
    if(&(predata->perm)){
        for(i = 0; i < total_row; ++i) r[i] = tempr->val[predata->perm.val[i]];
    }

    rr.row = predata->A->blocks[3]->row; rr.val = r+predata->A->blocks[0]->row;
    zz.row = predata->A->blocks[3]->col; zz.val = z+predata->A->blocks[0]->row;

    // Schwarz method on the interface part
    Schwarz_param *schwarz_param = predata->schwarz_param;
    Schwarz_data *schwarz_data = predata->schwarz_data;
    smoother_dcsr_Schwarz_forward(schwarz_data, schwarz_param, &zz, &rr);
    //directsolve_HAZ(&schwarz_data->A, &rr, &zz, 1);

    // AMG solve on the whole matrix
    AMG_param amgparam; param_amg_init(&amgparam);
    amgparam.cycle_type = predata->cycle_type;
    amgparam.smoother   = predata->smoother;
    amgparam.presmooth_iter  = predata->presmooth_iter;
    amgparam.postsmooth_iter = predata->postsmooth_iter;
    amgparam.relaxation = predata->relaxation;
    amgparam.coarse_scaling = predata->coarse_scaling;
    amgparam.tentative_smooth = predata->tentative_smooth;

    AMG_data_bdcsr *mgl = predata->mgl_data;
    mgl->b.row=total_row; array_cp(total_row, r, mgl->b.val); // residual is an input
    mgl->x.row=total_col; array_cp(total_row, z, mgl->x.val);
    //dvec_set(total_col, &mgl->x, 0.0);

    for ( i=maxit; i--; ) mgcycle_bdcsr(mgl,&amgparam);

    // update residual
    bdcsr_aAxpy(-1.0, predata->A, mgl->x.val, r);

    // setup z
    dvec_set(zz.row, &zz, 0.0);

    // Schwarz method on the interface part
    smoother_dcsr_Schwarz_backward(schwarz_data, schwarz_param, &zz, &rr);

    // update solution
    array_axpy(total_col, 1.0, mgl->x.val, z);

    // permute back the solution (USING precond_data_bdcsr temp work variable)
    if(&(predata->perm)) {
        predata->w = (REAL*)calloc(total_row, sizeof(REAL*));
        array_cp(total_row, z, predata->w);
        for(i = 0; i < total_row; ++i) z[predata->perm.val[i]] = predata->w[i];
        free(predata->w);
    }

    // restore residual
    array_cp(total_row, tempr->val, r);

}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
