/*! \file src/solver/mgcycle.c
 *
 *  Abstract multigrid cycle
 *  Routines for algebraic multigrid cycles
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 12/25/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 * \note   Done cleanup for releasing -- Xiaozhe Hu 03/12/2017 & 08/28/2021
 *
 * \todo Combine pre- and post-smoothing and use a flag to make sure the symmetry whenever is necessary -- Xiaozhe Hu
 *
 */

#include "hazmath.h"

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/
/***********************************************************************************************/
/**
 * \fn static void coarse_itsolver(dCSRmat *A, dvector *b, dvector *x,
 *                                 const REAL ctol, const SHORT prt_lvl)
 *
 * \brief Iterative on the coarset level
 *
 * \param  A         pointer to matrix data
 * \param  b         pointer to rhs data
 * \param  x         pointer to sol data
 * \param  ctol      tolerance for the coarsest level
 * \param  prt_lvl   level of output
 *
 */
static void coarse_itsolver(dCSRmat *A,
                            dvector *b,
                            dvector *x,
                            const REAL ctol,
                            const SHORT prt_lvl)
{
    const INT n = A->row;
    const INT maxit = 20*n;

    INT status = dcsr_pcg(A, b, x, NULL, ctol, maxit, 1, 0);

    // If CG fails to converge, use GMRES as another safe net
    if ( status < 0 ) {
        status = dcsr_pvgmres(A, b, x, NULL, ctol, maxit, 30, 1, 0);
    }

    if ( status < 0 && prt_lvl >= PRINT_MORE ) {
        printf("### HAZMATH WARNING: Coarse level solver failed to converge!\n");
    }
}

/***********************************************************************************************/
/**
 * \fn static void dcsr_presmoothing (const SHORT smoother, dCSRmat *A,
 *                                         dvector *b, dvector *x,
 *                                         const INT nsweeps, const INT istart,
 *                                         const INT iend, const INT istep,
 *                                         const REAL relax)
 *
 * \brief  Pre-smoothing
 *
 * \param  smoother  type of smoother
 * \param  A         pointer to matrix data
 * \param  b         pointer to rhs data
 * \param  x         pointer to sol data
 * \param  nsweeps   number of smoothing sweeps
 * \param  istart    starting index
 * \param  iend      ending index
 * \param  istep     step size
 * \param  relax     relaxation parameter for SOR-type smoothers
 *
 */
static void dcsr_presmoothing(SHORT smoother,
                              dCSRmat *A,
                              dvector *b,
                              dvector *x,
                              const INT nsweeps,
                              const INT istart,
                              const INT iend,
                              const INT istep,
                              const REAL relax)
{

    switch (smoother) {

        case SMOOTHER_GS:
            smoother_dcsr_gs(x, istart, iend, istep, A, b, nsweeps);
            break;

        case SMOOTHER_SGS:
            smoother_dcsr_sgs(x, A, b, nsweeps);
            break;

        case SMOOTHER_JACOBI:
            smoother_dcsr_jacobi(x, istart, iend, istep, A, b, nsweeps);
            break;

        case SMOOTHER_L1DIAG:
            smoother_dcsr_L1diag(x, istart, iend, istep, A, b, nsweeps);
            break;

        case SMOOTHER_SOR:
            smoother_dcsr_sor(x, istart, iend, istep, A, b, nsweeps, relax);
            break;

        case SMOOTHER_SSOR:
            smoother_dcsr_sor(x, istart, iend, istep, A, b, nsweeps, relax);
            smoother_dcsr_sor(x, iend, istart,-istep, A, b, nsweeps, relax);
            break;

        case SMOOTHER_GSOR:
            smoother_dcsr_gs (x, istart, iend, istep, A, b, nsweeps);
            smoother_dcsr_sor(x, iend, istart,-istep, A, b, nsweeps, relax);
            break;

        case SMOOTHER_SGSOR:
            smoother_dcsr_gs (x, istart, iend, istep, A, b, nsweeps);
            smoother_dcsr_gs (x, iend, istart,-istep, A, b, nsweeps);
            smoother_dcsr_sor(x, istart, iend, istep, A, b, nsweeps, relax);
            smoother_dcsr_sor(x, iend, istart,-istep, A, b, nsweeps, relax);
            break;

        case SMOOTHER_CG:
            dcsr_pcg(A, b, x, NULL, 1e-3, nsweeps, 1, PRINT_NONE);
            break;

        case SMOOTHER_USERDEF:
            /*if (smoother_function) {
                printf("Smoother function is not NULL! Let's try to use it. \n");
                smoother_function(x, istart, iend, istep, A, b, nsweeps, relax);
                printf("Smoother function ran successfully! But the correct result is not guaranteed. \n");
            }
            else {
                printf("### ERROR: User defined smoother is actually not defined (NULL ptr)!\n");
                check_error(ERROR_AMG_SMOOTH_TYPE, __FUNCTION__);
            } */
            /*printf("Smoother function is not NULL! Let's try to use it. \n");
            smmv->data = smoother_data_alloc(1);
            fenics_smoother_data_setup(istart, iend, istep, nsweeps, relax, A, smmv->data);

            py_callback_eval(b->val, x->val, smmv);
            printf("Smoother function ran successfully! But the correct result is not guaranteed. \n");
            smoother_data_free(smmv->data);
            */

            printf("Running GS just in case. \n");
            smoother_dcsr_gs(x, iend, istart, istep, A, b, nsweeps);
            break;

        default:
            printf("### ERROR: Wrong smoother type %d!\n", smoother);
            check_error(ERROR_INPUT_PAR, __FUNCTION__);
    }
}

/***********************************************************************************************/
/**
 * \fn static void dcsr_postsmoothing (const SHORT smoother, dCSRmat *A,
 *                                          dvector *b, dvector *x,
 *                                          const INT nsweeps, const INT istart,
 *                                          const INT iend, const INT istep,
 *                                          const REAL relax)
 *
 * \brief  Post-smoothing
 *
 * \param  smoother  type of smoother
 * \param  A         pointer to matrix data
 * \param  b         pointer to rhs data
 * \param  x         pointer to sol data
 * \param  nsweeps   number of smoothing sweeps
 * \param  istart    starting index
 * \param  iend      ending index
 * \param  istep     step size
 * \param  relax     relaxation parameter for SOR-type smoothers
 *
 */
static void dcsr_postsmoothing(SHORT smoother,
                               dCSRmat *A,
                               dvector *b,
                               dvector *x,
                               const INT nsweeps,
                               const INT istart,
                               const INT iend,
                               const INT istep,
                               const REAL relax)
{
    // const short smoother = smmv->type;

    switch (smoother) {

        case SMOOTHER_GS:
            smoother_dcsr_gs(x, iend, istart, istep, A, b, nsweeps);
            break;

        case SMOOTHER_SGS:
            smoother_dcsr_sgs(x, A, b, nsweeps);
            break;

        case SMOOTHER_JACOBI:
            smoother_dcsr_jacobi(x, iend, istart, istep, A, b, nsweeps);
            break;

        case SMOOTHER_L1DIAG:
            smoother_dcsr_L1diag(x, iend, istart, istep, A, b, nsweeps);
            break;

        case SMOOTHER_SOR:
            smoother_dcsr_sor(x, iend, istart, istep, A, b, nsweeps, relax);
            break;

        case SMOOTHER_SSOR:
            smoother_dcsr_sor(x, istart, iend, -istep, A, b, nsweeps, relax);
            smoother_dcsr_sor(x, iend, istart,  istep, A, b, nsweeps, relax);
            break;

        case SMOOTHER_GSOR:
            smoother_dcsr_sor(x, istart, iend, -istep, A, b, nsweeps, relax);
            smoother_dcsr_gs (x, iend, istart,  istep, A, b, nsweeps);
            break;

        case SMOOTHER_SGSOR:
            smoother_dcsr_sor(x, istart, iend, -istep, A, b, nsweeps, relax);
            smoother_dcsr_sor(x, iend, istart,  istep, A, b, nsweeps, relax);
            smoother_dcsr_gs (x, istart, iend, -istep, A, b, nsweeps);
            smoother_dcsr_gs (x, iend, istart,  istep, A, b, nsweeps);
            break;

        case SMOOTHER_CG:
            dcsr_pcg(A, b, x, NULL, 1e-3, nsweeps, 1, PRINT_NONE);
            break;

        case SMOOTHER_USERDEF:
            /*if (smoother_function) {
                printf("Smoother function is not NULL! Let's try to use it. \n");
                smoother_function(x, istart, iend, istep, A, b, nsweeps, relax);
                printf("Smoother function ran successfully! But the correct result is not guaranteed. \n");
                printf("Running GS just in case. \n");
                smoother_dcsr_gs(x, iend, istart, istep, A, b, nsweeps);
            }
            else {
                printf("### ERROR: User defined smoother is actually not defined (NULL ptr)!\n");
                check_error(ERROR_AMG_SMOOTH_TYPE, __FUNCTION__);
            } */
            /*printf("Smoother function is not NULL! Let's try to use it. \n");
            smmv->data = smoother_data_alloc(1);
            fenics_smoother_data_setup(istart, iend, istep, nsweeps, relax, A, smmv->data);

            // python_callback(b->val, x->val, smmv);
            printf("Smoother function ran successfully! But the correct result is not guaranteed. \n");
            smoother_data_free(smmv->data);
            */

            printf("Running GS just in case. \n");
            smoother_dcsr_gs(x, iend, istart, istep, A, b, nsweeps);
            break;

        default:
            printf("### ERROR: Wrong smoother type %d!\n", smoother);
            check_error(ERROR_INPUT_PAR, __FUNCTION__);
    }
}




/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/**
 * \fn void mgcycle (AMG_data *mgl, AMG_param *param)
 *
 * \brief Solve Ax=b with non-recursive multigrid cycle (V- and W-cycle)
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \author Xiaozhe Hu
 * \date   12/25/2015
 *
 */
void mgcycle(AMG_data *mgl,
             AMG_param *param)
{
    const SHORT  prtlvl = param->print_level;
    const SHORT  amg_type = param->AMG_type;
    const SHORT  smoother = param->smoother;
    const SHORT  cycle_type = param->cycle_type;
    const SHORT  coarse_solver = param->coarse_solver;
    const SHORT  nl = mgl[0].num_levels;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol * 1e-2;
    //smoother_matvec *smmv = (smoother_matvec*)mgl->wdata;

    // Schwarz parameters
    Schwarz_param swzparam;

    // local variables
    REAL alpha = 1.0;
    INT  num_lvl[MAX_AMG_LVL] = {0}, l = 0;

ForwardSweep:
    while ( l < nl-1 ) {

        num_lvl[l]++;

        // pre-smoothing with Schwarz method
        if ( l < mgl->Schwarz_levels ) {
            swzparam.Schwarz_blksolver = mgl[l].Schwarz.blk_solver;
            switch (mgl[l].Schwarz.Schwarz_type) {
                case SCHWARZ_SYMMETRIC:
                    smoother_dcsr_Schwarz_forward(&mgl[l].Schwarz, &swzparam,
                                                  &mgl[l].x, &mgl[l].b);
                    smoother_dcsr_Schwarz_backward(&mgl[l].Schwarz, &swzparam,
                                                   &mgl[l].x, &mgl[l].b);
                    break;
                default:
                    smoother_dcsr_Schwarz_forward(&mgl[l].Schwarz, &swzparam,
                                                  &mgl[l].x, &mgl[l].b);
                    break;
            }
        }
        else
        { // pre-smoothing with standard smoothers
          dcsr_presmoothing(smoother, &mgl[l].A, &mgl[l].b, &mgl[l].x,
                            param->presmooth_iter, 0, mgl[l].A.row-1, 1,
                            relax);
        }

        // form residual r = b - A x
        array_cp(mgl[l].A.row, mgl[l].b.val, mgl[l].w.val);
        dcsr_aAxpy(-1.0,&mgl[l].A, mgl[l].x.val, mgl[l].w.val);

        // restriction r1 = R*r0
        switch ( amg_type ) {
            case UA_AMG:
                dcsr_mxv_agg(&mgl[l].R, mgl[l].w.val, mgl[l+1].b.val);
                break;
            default:
                dcsr_mxv(&mgl[l].R, mgl[l].w.val, mgl[l+1].b.val);
                break;
        }

        // prepare for the next level
        ++l; dvec_set(mgl[l].A.row, &mgl[l].x, 0.0);

    }

    // If AMG only has one level or we have arrived at the coarsest level,
    // call the coarse space solver:
    switch ( coarse_solver ) {

#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
            // use UMFPACK direct solver on the coarsest level
            umfpack_solve(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, mgl[nl-1].Numeric, 0);
            break;
        }
#endif
        default:
            // use iterative solver on the coarsest level
            coarse_itsolver(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, tol, prtlvl);
            break;

    }

    // BackwardSweep:
    while ( l > 0 ) {

        --l;

        // find the optimal scaling factor alpha
        if ( param->coarse_scaling == ON ) {
            alpha = array_dotprod(mgl[l+1].A.row, mgl[l+1].x.val, mgl[l+1].b.val)
                  / dcsr_vmv(&mgl[l+1].A, mgl[l+1].x.val, mgl[l+1].x.val);
            alpha = MIN(alpha, 1.0); // Add this for safety! --Chensong on 10/04/2014
        }

        // prolongation u = u + alpha*P*e1
        switch ( amg_type ) {
            case UA_AMG:
                dcsr_aAxpy_agg(alpha, &mgl[l].P, mgl[l+1].x.val, mgl[l].x.val);
                break;
            default:
                dcsr_aAxpy(alpha, &mgl[l].P, mgl[l+1].x.val, mgl[l].x.val);
                break;
        }

        // post-smoothing with Schwarz method
        if ( l < mgl->Schwarz_levels ) {
          swzparam.Schwarz_blksolver = mgl[l].Schwarz.blk_solver;
            switch (mgl[l].Schwarz.Schwarz_type) {
                case SCHWARZ_SYMMETRIC:
                    smoother_dcsr_Schwarz_backward(&mgl[l].Schwarz, &swzparam,
                                                   &mgl[l].x, &mgl[l].b);
                    smoother_dcsr_Schwarz_forward(&mgl[l].Schwarz, &swzparam,
                                                  &mgl[l].x, &mgl[l].b);
                    break;
                default:
                    smoother_dcsr_Schwarz_backward(&mgl[l].Schwarz, &swzparam,
                                                   &mgl[l].x, &mgl[l].b);
                    break;
            }
        }
        else
        { // post-smoothing with standard methods
          dcsr_postsmoothing(smoother, &mgl[l].A, &mgl[l].b, &mgl[l].x,
                             param->postsmooth_iter, 0, mgl[l].A.row-1, -1,
                             relax);
        }

        if ( num_lvl[l] < cycle_type ) break;
        else num_lvl[l] = 0;
    }

    if ( l > 0 ) goto ForwardSweep;


}


/**
 * \fn void amli (AMG_data *mgl, AMG_param *param, INT level)
 *
 * \brief Solve Ax=b with recursive AMLI-cycle
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 * \param level  Current level
 *
 * \author Xiaozhe Hu
 * \date   01/23/2011
 *
 * \note AMLI polynomial computed by the best approximation of 1/x.
 *       Refer to Johannes K. Kraus, Panayot S. Vassilevski, Ludmil T. Zikatanov,
 *       "Polynomial of best uniform approximation to $x^{-1}$ and smoothing in
 *        two-level methods", 2013.
 *
 */
void amli(AMG_data *mgl,
          AMG_param *param,
          INT level)
{
    const SHORT  amg_type=param->AMG_type;
    const SHORT  prtlvl = param->print_level;
    const SHORT  smoother = param->smoother;
    const SHORT  coarse_solver = param->coarse_solver;
    const SHORT  degree= param->amli_degree;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol*1e-4;

    // local variables
    REAL   alpha  = 1.0;
    REAL * coef   = param->amli_coef;

    dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x;   // fine level b and x
    dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x

    dCSRmat *A0 = &mgl[level].A;   // fine level matrix
    dCSRmat *A1 = &mgl[level+1].A; // coarse level matrix

    const INT m0 = A0->row, m1 = A1->row;

    REAL     *r        = mgl[level].w.val;      // work array for residual
    REAL     *r1       = mgl[level+1].w.val+m1; // work array for residual

    if ( prtlvl >= PRINT_MOST )
        printf("AMLI level %d, smoother %d.\n", level, smoother);

    if ( level < mgl[level].num_levels-1 ) {

        // presmoothing
        dcsr_presmoothing(smoother,A0,b0,e0,param->presmooth_iter,
                          0,m0-1,1,relax);

        // form residual r = b - A x
        array_cp(m0,b0->val,r);
        dcsr_aAxpy(-1.0,A0,e0->val,r);

        // restriction r1 = R*r0
        switch (amg_type) {
            case UA_AMG:
                dcsr_mxv_agg(&mgl[level].R, r, b1->val); break;
            default:
                dcsr_mxv(&mgl[level].R, r, b1->val); break;
        }

        // coarse grid correction
        {
            INT i;

            array_cp(m1,b1->val,r1);

            for ( i=1; i<=degree; i++ ) {
                dvec_set(m1,e1,0.0);
                amli(mgl, param, level+1);

                // b1 = (coef[degree-i]/coef[degree])*r1 + A1*e1;
                // First, compute b1 = A1*e1
                dcsr_mxv(A1, e1->val, b1->val);
                // Then, compute b1 = b1 + (coef[degree-i]/coef[degree])*r1
                array_axpy(m1, coef[degree-i]/coef[degree], r1, b1->val);
            }

            dvec_set(m1,e1,0.0);
            amli(mgl, param, level+1);
        }

        // find the optimal scaling factor alpha
        array_ax(m1, coef[degree], e1->val);
        if ( param->coarse_scaling == ON ) {
            alpha = array_dotprod(m1, e1->val, r1)
            / dcsr_vmv(A1, e1->val, e1->val);
            alpha = MIN(alpha, 1.0);
        }

        // prolongation e0 = e0 + alpha * P * e1
        switch (amg_type) {
            case UA_AMG:
                dcsr_aAxpy_agg(alpha, &mgl[level].P, e1->val, e0->val);
                break;
            default:
                dcsr_aAxpy(alpha, &mgl[level].P, e1->val, e0->val);
                break;
        }

        // postsmoothing
        dcsr_postsmoothing(smoother,A0,b0,e0,param->postsmooth_iter,
                           0,m0-1,-1,relax);

    }

    else { // coarsest level solver

        switch (coarse_solver) {

#if WITH_SUITESPARSE
            case SOLVER_UMFPACK:
                // use UMFPACK direct solver on the coarsest level //
                umfpack_solve(A0, b0, e0, mgl[level].Numeric, 0);
                break;
#endif

            default:
                /* use iterative solver on the coarsest level */
                coarse_itsolver(A0, b0, e0, tol, prtlvl);

        }

    }

}

/**
 * \fn void nl_amli(AMG_data *mgl, AMG_param *param, INT level, INT num_levels)
 *
 * \brief Solve Ax=b with recursive nonlinear AMLI-cycle
 *
 * \param mgl         Pointer to AMG_data data
 * \param param       Pointer to AMG parameters
 * \param level       Current level
 * \param num_levels  Total number of levels
 *
 * \author Xiaozhe Hu
 * \date   04/06/2010
 *
 * \note Refer to Xiazhe Hu, Panayot S. Vassilevski, Jinchao Xu
 *       "Comparative Convergence Analysis of Nonlinear AMLI-cycle Multigrid", 2013.
 *
 */
void nl_amli (AMG_data *mgl,
              AMG_param *param,
              INT level,
              INT num_levels)
{
    const SHORT  amg_type=param->AMG_type;
    const SHORT  prtlvl = param->print_level;
    const SHORT  smoother = param->smoother;
    const SHORT  coarse_solver = param->coarse_solver;
    const REAL   relax = param->relaxation;
    const INT    maxit = param->amli_degree+1;
    const REAL   tol = param->tol*1e-4;

    dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x;   // fine level b and x
    dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x

    dCSRmat *A0 = &mgl[level].A;   // fine level matrix
    dCSRmat *A1 = &mgl[level+1].A; // coarse level matrix

    const INT m0 = A0->row, m1 = A1->row;

    REAL    *r = mgl[level].w.val;      // work array for residual

    dvector uH, bH;  // for coarse level correction
    uH.row = m1; uH.val = mgl[level+1].w.val + m1;
    bH.row = m1; bH.val = mgl[level+1].w.val + 2*m1;

    if ( prtlvl >= PRINT_MOST )
        printf("Nonlinear AMLI level %d, smoother %d.\n", num_levels, smoother);

    if ( level < num_levels-1 ) {

        // presmoothing
        dcsr_presmoothing(smoother,A0,b0,e0,param->presmooth_iter,
                          0,m0-1,1,relax);

        // form residual r = b - A x
        array_cp(m0,b0->val,r);
        dcsr_aAxpy(-1.0,A0,e0->val,r);

        // restriction r1 = R*r0
        switch (amg_type) {
            case UA_AMG:
                dcsr_mxv_agg(&mgl[level].R, r, b1->val);
                break;
            default:
                dcsr_mxv(&mgl[level].R, r, b1->val);
                break;
        }

        // call nonlinear AMLI-cycle recursively
        {
            dvec_set(m1,e1,0.0);

            // The coarsest problem is solved exactly.
            // No need to call krylov method on second coarest level
            if ( level == num_levels-2 ) {
                nl_amli(&mgl[level+1], param, 0, num_levels-1);
            }
            else { // recursively call preconditioned Krylov method on coarse grid

                precond_data pcdata;

                param_amg_to_prec(&pcdata, param);
                pcdata.maxit = 1;
                pcdata.max_levels = num_levels-1;
                pcdata.mgl_data = &mgl[level+1];

                precond pc;
                pc.data = &pcdata;
                pc.fct = precond_nl_amli;

                array_cp(m1, b1->val, bH.val);
                array_cp(m1, e1->val, uH.val);

                switch (param->nl_amli_krylov_type) {
                    case SOLVER_GCG: // Use GCG
                        dcsr_pgcg(A1,&bH,&uH,&pc,tol,maxit,1,PRINT_NONE);
                        break;
                    default: // Use FGMRES
                        dcsr_pvfgmres(A1,&bH,&uH,&pc,tol,maxit,30,1,PRINT_NONE);
                        break;
                }

                array_cp(m1, bH.val, b1->val);
                array_cp(m1, uH.val, e1->val);
            }

        }

        // prolongation e0 = e0 + P*e1
        switch (amg_type) {
            case UA_AMG:
                dcsr_aAxpy_agg(1.0, &mgl[level].P, e1->val, e0->val);
                break;
            default:
                dcsr_aAxpy(1.0, &mgl[level].P, e1->val, e0->val);
                break;
        }

        // postsmoothing
        dcsr_postsmoothing(smoother,A0,b0,e0,param->postsmooth_iter,
                           0,m0-1,-1,relax);

    }

    else { // coarsest level solver

        switch (coarse_solver) {

#if WITH_SUITESPARSE
            case SOLVER_UMFPACK:
                // use UMFPACK direct solver on the coarsest level //
                umfpack_solve(A0, b0, e0, mgl[level].Numeric, 0);
                break;
#endif

            default:
                /* use iterative solver on the coarsest level */
                coarse_itsolver(A0, b0, e0, tol, prtlvl);

        }

    }

}

/**
 * \fn void mgcycle_add(AMG_data *mgl, AMG_param *param)
 *
 * \brief Solve Ax=b with additive multigrid cycle
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \author Xiaozhe Hu
 * \date   05/28/2020
 *
 */
void mgcycle_add(AMG_data *mgl,
                 AMG_param *param)
{
    const SHORT  prtlvl = param->print_level;
    const SHORT  amg_type = param->AMG_type;
    const SHORT  smoother = param->smoother;
    const SHORT  coarse_solver = param->coarse_solver;
    const SHORT  nl = mgl[0].num_levels;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol * 1e-2;

    // Schwarz parameters
    Schwarz_param swzparam;

    // local variables
    REAL alpha = 1.0;
    INT l = 0;

    // compute the residual on the finest level
    array_cp(mgl[0].A.row, mgl[0].b.val, mgl[0].w.val);
    dcsr_aAxpy(-1.0,&mgl[0].A, mgl[0].x.val, mgl[0].w.val);

    // main loop
    while ( l < nl-1 ) {

        // pre-smoothing with Schwarz method
        if ( l < mgl->Schwarz_levels ) {
            swzparam.Schwarz_blksolver = mgl[l].Schwarz.blk_solver;
            switch (mgl[l].Schwarz.Schwarz_type) {
                case SCHWARZ_SYMMETRIC:
                    smoother_dcsr_Schwarz_forward(&mgl[l].Schwarz, &swzparam,
                                                  &mgl[l].x, &mgl[l].b);
                    smoother_dcsr_Schwarz_backward(&mgl[l].Schwarz, &swzparam,
                                                   &mgl[l].x, &mgl[l].b);
                    break;
                default:
                    smoother_dcsr_Schwarz_forward(&mgl[l].Schwarz, &swzparam,
                                                  &mgl[l].x, &mgl[l].b);
                    break;
            }
        }
        else
        { // pre-smoothing with standard smoothers
          dcsr_presmoothing(smoother, &mgl[l].A, &mgl[l].b, &mgl[l].x,
                            param->presmooth_iter, 0, mgl[l].A.row-1, 1,
                            relax);
        }

        // restriction rH = R*rh (restrict residual, not the right-hand-side)
        if (l==0) { // on the finest level, residual is stored in w, not b.
          switch ( amg_type ) {
              case UA_AMG:
                  dcsr_mxv_agg(&mgl[l].R, mgl[l].w.val, mgl[l+1].b.val);
                  break;
              default:
                  dcsr_mxv(&mgl[l].R, mgl[l].w.val, mgl[l+1].b.val);
                  break;
          }
        }
        else { // on coarse levels, residual is stored in b.
          switch ( amg_type ) {
              case UA_AMG:
                  dcsr_mxv_agg(&mgl[l].R, mgl[l].b.val, mgl[l+1].b.val);
                  break;
              default:
                  dcsr_mxv(&mgl[l].R, mgl[l].b.val, mgl[l+1].b.val);
                  break;
          }
        }

        // prepare for the next level
        ++l; dvec_set(mgl[l].A.row, &mgl[l].x, 0.0);

    }

    // If AMG only has one level or we have arrived at the coarsest level,
    // call the coarse space solver:
    switch ( coarse_solver ) {

#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
            // use UMFPACK direct solver on the coarsest level
            umfpack_solve(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, mgl[nl-1].Numeric, 0);
            break;
        }
#endif
        default:
            // use iterative solver on the coarsest level
            coarse_itsolver(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, tol, prtlvl);
            break;

    }

    // BackwardSweep (update solution only, no postsmoothing)
    while ( l > 0 ) {

        --l;

        // find the optimal scaling factor alpha
        if ( param->coarse_scaling == ON ) {
            alpha = array_dotprod(mgl[l+1].A.row, mgl[l+1].x.val, mgl[l+1].b.val)
                  / dcsr_vmv(&mgl[l+1].A, mgl[l+1].x.val, mgl[l+1].x.val);
            alpha = MIN(alpha, 1.0);
        }

        // prolongation u = u + alpha*P*e1
        switch ( amg_type ) {
            case UA_AMG:
                dcsr_aAxpy_agg(alpha, &mgl[l].P, mgl[l+1].x.val, mgl[l].x.val);
                break;
            default:
                dcsr_aAxpy(alpha, &mgl[l].P, mgl[l+1].x.val, mgl[l].x.val);
                break;
        }

    }

}


/**
 * \fn void mgcycle_add_update(AMG_data *mgl, AMG_param *param)
 *
 * \brief Solve Ae=r with additive multigrid cycle
 * \note This subroutine assumes that the input right hand side is residual
 *       and the output solutio is the update.  Therefore, the initial guess has to zero
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \author Xiaozhe Hu
 * \date   05/28/2020
 *
 */
void mgcycle_add_update(AMG_data *mgl,
                 AMG_param *param)
{
    const SHORT  prtlvl = param->print_level;
    const SHORT  amg_type = param->AMG_type;
    const SHORT  smoother = param->smoother;
    const SHORT  coarse_solver = param->coarse_solver;
    const SHORT  nl = mgl[0].num_levels;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol * 1e-2;

    // Schwarz parameters
    Schwarz_param swzparam;

    // local variables
    REAL alpha = 1.0;
    INT l = 0;

    // make sure the initial guess is zero
    dvec_set(mgl[0].A.row, &mgl[0].x, 0.0);

    // main loop
    while ( l < nl-1 ) {

        // pre-smoothing with Schwarz method
        if ( l < mgl->Schwarz_levels ) {
            swzparam.Schwarz_blksolver = mgl[l].Schwarz.blk_solver;
            switch (mgl[l].Schwarz.Schwarz_type) {
                case SCHWARZ_SYMMETRIC:
                    smoother_dcsr_Schwarz_forward(&mgl[l].Schwarz, &swzparam,
                                                  &mgl[l].x, &mgl[l].b);
                    smoother_dcsr_Schwarz_backward(&mgl[l].Schwarz, &swzparam,
                                                   &mgl[l].x, &mgl[l].b);
                    break;
                default:
                    smoother_dcsr_Schwarz_forward(&mgl[l].Schwarz, &swzparam,
                                                  &mgl[l].x, &mgl[l].b);
                    break;
            }
        }
        else
        { // pre-smoothing with standard smoothers
          dcsr_presmoothing(smoother, &mgl[l].A, &mgl[l].b, &mgl[l].x,
                            param->presmooth_iter, 0, mgl[l].A.row-1, 1,
                            relax);
        }

        // restriction rH = R*rh (restrict residual)
        switch ( amg_type ) {
            case UA_AMG:
                dcsr_mxv_agg(&mgl[l].R, mgl[l].b.val, mgl[l+1].b.val);
                break;
            default:
                dcsr_mxv(&mgl[l].R, mgl[l].b.val, mgl[l+1].b.val);
                break;
        }

        // prepare for the next level
        ++l; dvec_set(mgl[l].A.row, &mgl[l].x, 0.0);

    }

    // If AMG only has one level or we have arrived at the coarsest level,
    // call the coarse space solver:
    switch ( coarse_solver ) {

#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
            // use UMFPACK direct solver on the coarsest level
            umfpack_solve(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, mgl[nl-1].Numeric, 0);
            break;
        }
#endif
        default:
            // use iterative solver on the coarsest level
            coarse_itsolver(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, tol, prtlvl);
            break;

    }

    // BackwardSweep (update solution only, no postsmoothing)
    while ( l > 0 ) {

        --l;

        // find the optimal scaling factor alpha
        if ( param->coarse_scaling == ON ) {
            alpha = array_dotprod(mgl[l+1].A.row, mgl[l+1].x.val, mgl[l+1].b.val)
                  / dcsr_vmv(&mgl[l+1].A, mgl[l+1].x.val, mgl[l+1].x.val);
            alpha = MIN(alpha, 1.0);
        }

        // prolongation u = u + alpha*P*e1
        switch ( amg_type ) {
            case UA_AMG:
                dcsr_aAxpy_agg(alpha, &mgl[l].P, mgl[l+1].x.val, mgl[l].x.val);
                break;
            default:
                dcsr_aAxpy(alpha, &mgl[l].P, mgl[l+1].x.val, mgl[l].x.val);
                break;
        }

    }

}


/**
 * \fn void cascadic_eigen(AMG_data *mgl, AMG_param *param, INT level, INT num_levels)
 *
 * \brief Solve Ax=\lambda x using cascadic multigrid
 *
 * \param mgl         Pointer to AMG_data data
 * \param param       Pointer to AMG parameters
 * \param level       Current level
 * \param num_levels  Total number of levels
 *
 * \author Xiaozhe Hu
 * \date   06/25/2019
 *
 * \note Refer to John Urschel, Xiaozhe Hu, Jinchao Xu, and Ludmil Zikatanov
 *       "A Cascadic Multigrid Algorithm for Computing the Fiedler Vector of Graph Laplacians", 2015.
 *
 */
void cascadic_eigen(AMG_data *mgl,
                    AMG_param *param,
                    INT level,
                    const INT num_eigen)
{

  // local variables
  INT i;

  dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x;   // fine level b and x
  //dvector *b1 = &mgl[level+1].b, // coarse level b
  dvector *e1 = &mgl[level+1].x; // coarse level x

  dCSRmat *A0 = &mgl[level].A;   // fine level matrix
  dCSRmat *A1 = &mgl[level+1].A; // coarse level matrix

  const INT n0 = A0->row, n1 = A1->row;

  INT nsmooth = pow(2, level-1)*param->postsmooth_iter;

  if ( level < mgl[level].num_levels-1 )
  {

    // coarse grid correction
    cascadic_eigen(mgl, param, level+1, num_eigen);

    // prolongation
    dvec_set(e0->row, e0, 0.0);
    for (i=0; i<num_eigen; i++)
    {
      dcsr_aAxpy(1.0, &mgl[level].P, &(e1->val[i*n1]), &(e0->val[i*n0]));
    }

    // set zero right hand side
    dvec_set(A0->row, b0, 0.0);

    // postsmoothing
    smoother_dcsr_sgs_graph_eigen(e0, A0, b0, nsmooth, num_eigen);

  }
  else // coarsest grid
  {

    // random initial guess
    dvec_rand(num_eigen*A0->row, e0);

    // set zero right hand side
    dvec_set(A0->row, b0, 0.0);

    // solve eigenvalue problem on coarsest grid
    smoother_dcsr_sgs_graph_eigen(e0, A0, b0, 10*A0->row, num_eigen);

  }

}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
