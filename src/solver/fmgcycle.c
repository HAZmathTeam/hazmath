/*! \file src/solver/fmgcycle.c
 *
 *  Abstract (fractional) multigrid cycle
 *  Routines for algebraic multigrid cycles with fractional smoothers
 *
 *  Created by Ana Budisa 2020-05-08 from original mgcycle.c (by HAZ people)
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *
 *
 */

#include "hazmath.h"

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/
/**
 * \fn static void coarse_fitsolver(dCSRmat *A, dvector *b, dvector *x,
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
static void coarse_fitsolver(dCSRmat *A,
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

/**
 * \fn static void coarse_fracinv(dCSRmat *A, dvector *b, dvector *x, dCSRmat *M,
 *                                const REAL ctol, const SHORT prt_lvl)
 *
 * \brief Inverse of A^s via a generalized eigenvalue solver
 *
 * \param  A         pointer to stiffness matrix data
 * \param  b         pointer to rhs data
 * \param  x         pointer to sol data
 * \param  M         pointer to mass matrix data
 * \param  ctol      tolerance for the coarsest level
 * \param  prt_lvl   level of output
 *
 */
static void coarse_fracinv(dCSRmat *A,
                           dvector *b,
                           dvector *x,
                           dCSRmat *M,
                           const REAL power)
{
    const INT n = A->row;
    INT i;

    // first initialize eigenvalues and vectors
    REAL *lambda = calloc(n, sizeof(REAL));
    REAL *v = calloc(n*n, sizeof(REAL));
    // set x to zero
    dvec_set(n, x, 0.0);

    // call eigensolver: A v_i = lambda_i M v_i,  i = 1, ..., ncol
    // then A V = M V Lambda
    eigsymm(A, M, lambda, v);

    // x = A^-s b = V Lambda^-s V^T b
    REAL *temp = calloc(n, sizeof(REAL)); // temp = 0
    ddense_atbyv(n, temp, v, b->val, n); // temp = V^T b + temp
    for(i = 0; i < n; ++i) temp[i] *= pow(lambda[i], -power); // temp[i] = temp[i] * lambda_i^-s
    ddense_abyv(n, x->val, v, temp, n); // x = V * temp + x (if x = 0)

}

/**
 * \fn static void dcsr_fpresmoothing (const SHORT smoother, dCSRmat *A,
 *                                         dvector *b, dvector *x, dCSRmat *M,
 *                                         const REAL p,
 *                                         const INT nsweeps, const INT istart,
 *                                         const INT iend, const INT istep,
 *                                         const REAL relax)
 *
 * \brief  Fractional Pre-smoothing (now only Jacobi, GS and SGS)
 *
 * \param  smoother  type of smoother
 * \param  A         pointer to matrix data
 * \param  b         pointer to rhs data
 * \param  x         pointer to sol data
 * \param  M         pointer to mass matrix data
 * \param  p         fractional exponent
 * \param  nsweeps   number of smoothing sweeps
 * \param  istart    starting index
 * \param  iend      ending index
 * \param  istep     step size
 * \param  relax     relaxation parameter for SOR-type smoothers
 *
 */
static void dcsr_fpresmoothing(const SHORT smoother,
                               dCSRmat *A,
                               dvector *b,
                               dvector *x,
                               dCSRmat *M,
                               const REAL p,
                               const INT nsweeps,
                               const INT istart,
                               const INT iend,
                               const INT istep,
                               const REAL relax)
{
    switch (smoother) {

        case SMOOTHER_FGS:
            smoother_dcsr_fgs(x, istart, iend, istep, A, b, M, p, nsweeps);
            break;

        case SMOOTHER_FSGS:
            smoother_dcsr_fsgs(x, A, b, M, p, nsweeps);
            break;

        case SMOOTHER_FJACOBI:
            smoother_dcsr_fjacobi(x, istart, iend, istep, A, b, M, p, relax, nsweeps);
            break;

        default:
            printf("### ERROR: Wrong smoother type %lld!\n", (long long )smoother);
            check_error(ERROR_INPUT_PAR, __FUNCTION__);
    }
}

/**
 * \fn static void dcsr_fpostsmoothing (const SHORT smoother, dCSRmat *A,
 *                                          dvector *b, dvector *x, dCSRmat *M,
 *                                          const REAL p,
 *                                          const INT nsweeps, const INT istart,
 *                                          const INT iend, const INT istep,
 *                                          const REAL relax)
 *
 * \brief  Fractional Post-smoothing
 *
 * \param  smoother  type of smoother
 * \param  A         pointer to matrix data
 * \param  b         pointer to rhs data
 * \param  x         pointer to sol data
 * \param  M         pointer to mass matrix data
 * \param  p         fractional exponent
 * \param  nsweeps   number of smoothing sweeps
 * \param  istart    starting index
 * \param  iend      ending index
 * \param  istep     step size
 * \param  relax     relaxation parameter for SOR-type smoothers
 *
 */
static void dcsr_fpostsmoothing(const SHORT smoother,
                                dCSRmat *A,
                                dvector *b,
                                dvector *x,
                                dCSRmat *M,
                                const REAL p,
                                const INT nsweeps,
                                const INT istart,
                                const INT iend,
                                const INT istep,
                                const REAL relax)
{
    switch (smoother) {

        case SMOOTHER_FGS:
            smoother_dcsr_fgs(x, iend, istart, istep, A, b, M, p, nsweeps);
            break;

        case SMOOTHER_FSGS:
            smoother_dcsr_fsgs(x, A, b, M, p, nsweeps);
            break;

        case SMOOTHER_FJACOBI:
            smoother_dcsr_fjacobi(x, iend, istart, istep, A, b, M, p, relax, nsweeps);
            break;

        default:
            printf("### ERROR: Wrong smoother type %lld!\n", (long long )smoother);
            check_error(ERROR_INPUT_PAR, __FUNCTION__);
    }
}


/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/**
 * \fn void fmgcycle (AMG_data *mgl, AMG_param *param)
 *
 * \brief Solve Ax=b with non-recursive multigrid cycle (V- and W-cycle) with
 * fractional smoothers
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \author Ana Budisa (disclaimer: edited from original mgcycle)
 * \date   2020-05-08
 *
 */
void fmgcycle(AMG_data *mgl,
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
    const REAL   power = param->fpwr; // fractional exponent
    
    // Schwarz parameters
    Schwarz_param swzparam;

    // local variables
    REAL alpha = 1.0;
    INT  num_lvl[MAX_AMG_LVL] = {0}, l = 0,sch_type=SCHWARZ_SYMMETRIC_LOCAL;

ForwardSweep:
    while ( l < nl-1 ) {

        num_lvl[l]++;
        if ( l < mgl->Schwarz_levels ) {
            swzparam.Schwarz_blksolver = mgl[l].Schwarz.blk_solver;
            /* switch (mgl[l].Schwarz.Schwarz_type) { */
            /*     case SCHWARZ_SYMMETRIC: */
	    smoother_dcsr_Schwarz(&mgl[l].Schwarz,&mgl[l].x, &mgl[l].b,param->presmooth_iter);
	    /* smoother_dcsr_Schwarz_forward(&mgl[l].Schwarz, &swzparam, */
	    /*                               &mgl[l].x, &mgl[l].b); */
	    /* smoother_dcsr_Schwarz_backward(&mgl[l].Schwarz, &swzparam, */
	    /*                                &mgl[l].x, &mgl[l].b); */
            /*         break; */
            /*     default: */
            /*         smoother_dcsr_Schwarz_forward(&mgl[l].Schwarz, &swzparam, */
            /*                                       &mgl[l].x, &mgl[l].b); */
            /*         break; */
            /* } */
        }
        else
        { // pre-smoothing with standard smoothers
          dcsr_fpresmoothing(smoother, &mgl[l].A, &mgl[l].b, &mgl[l].x, &mgl[l].M,
                             power, param->presmooth_iter, 0, mgl[l].A.row-1, 1,
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

      //#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
            // use UMFPACK direct solver on the coarsest level
            hazmath_solve(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, mgl[nl-1].Numeric, 0);
            break;
        }
	  //#endif
        default:
            // use iterative solver on the coarsest level
            coarse_fitsolver(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, tol,
            prtlvl);
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
	  sch_type=mgl[l].Schwarz.Schwarz_type;
	  switch(sch_type){
	  case SCHWARZ_FORWARD:
	    mgl[l].Schwarz.Schwarz_type=SCHWARZ_BACKWARD;
	    smoother_dcsr_Schwarz(&mgl[l].Schwarz,&mgl[l].x, &mgl[l].b,param->postsmooth_iter);
	    mgl[l].Schwarz.Schwarz_type=sch_type;
	    break;
	  case SCHWARZ_FORWARD_LOCAL:
	    mgl[l].Schwarz.Schwarz_type=SCHWARZ_BACKWARD_LOCAL;
	    smoother_dcsr_Schwarz(&mgl[l].Schwarz,&mgl[l].x, &mgl[l].b,param->postsmooth_iter);
	    mgl[l].Schwarz.Schwarz_type=sch_type;
	    break;
	  case SCHWARZ_BACKWARD:
	    mgl[l].Schwarz.Schwarz_type=SCHWARZ_FORWARD;
	    smoother_dcsr_Schwarz(&mgl[l].Schwarz,&mgl[l].x, &mgl[l].b,param->postsmooth_iter);
	    mgl[l].Schwarz.Schwarz_type=sch_type;
	    break;
	  case SCHWARZ_BACKWARD_LOCAL:
	    mgl[l].Schwarz.Schwarz_type=SCHWARZ_FORWARD_LOCAL;
	    smoother_dcsr_Schwarz(&mgl[l].Schwarz,&mgl[l].x, &mgl[l].b,param->postsmooth_iter);
	    mgl[l].Schwarz.Schwarz_type=sch_type;
	    break;
	  default: //symmetric or symmetric local stays the same
	    smoother_dcsr_Schwarz(&mgl[l].Schwarz,&mgl[l].x, &mgl[l].b,param->postsmooth_iter);
	    break;
	  }
	  /* switch (mgl[l].Schwarz.Schwarz_type) { */
	  /*     case SCHWARZ_SYMMETRIC: */
	  /*         smoother_dcsr_Schwarz_backward(&mgl[l].Schwarz, &swzparam, */
	  /*                                        &mgl[l].x, &mgl[l].b); */
	  /*         smoother_dcsr_Schwarz_forward(&mgl[l].Schwarz, &swzparam, */
	  /*                                       &mgl[l].x, &mgl[l].b); */
	  /*         break; */
	  /*     default: */
	  /*         smoother_dcsr_Schwarz_backward(&mgl[l].Schwarz, &swzparam, */
	  /*                                        &mgl[l].x, &mgl[l].b); */
	  /*         break; */
	  /* } */
        }  else { // post-smoothing with standard methods
          dcsr_fpostsmoothing(smoother, &mgl[l].A, &mgl[l].b, &mgl[l].x, &mgl[l].M,power, param->postsmooth_iter, 0, mgl[l].A.row-1,-1, relax);
        }

        if ( num_lvl[l] < cycle_type ) break;
        else num_lvl[l] = 0;
    }

    if ( l > 0 ) goto ForwardSweep;


}

/**
 * \fn void fmgcycle_add_update(AMG_data *mgl, AMG_param *param)
 *
 * \brief Solve Ae=r with additive multigrid cycle with fractional smoothers
 * \note This subroutine assumes that the input right hand side is residual
 *       and the output solution is the update.  Therefore, the initial guess has to be zero
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \author Xiaozhe Hu, Ana Budisa
 * \date   2020-06-03
 *
 */
void fmgcycle_add_update(AMG_data *mgl,
                         AMG_param *param)
{
  // not used:  const SHORT  prtlvl = param->print_level;
    const SHORT  amg_type = param->AMG_type;
    const SHORT  smoother = param->smoother;
    const SHORT  coarse_solver = param->coarse_solver;
    const SHORT  nl = mgl[0].num_levels;
    const REAL   relax = param->relaxation;
    // not used:  const REAL   tol = param->tol * 1e-2;
    const REAL   power = param->fpwr; // fractional exponent

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
	    smoother_dcsr_Schwarz(&mgl[l].Schwarz,&mgl[l].x, &mgl[l].b,param->presmooth_iter);
            /* switch (mgl[l].Schwarz.Schwarz_type) { */
            /*     case SCHWARZ_SYMMETRIC: */
            /*         smoother_dcsr_Schwarz_forward(&mgl[l].Schwarz, &swzparam, */
            /*                                       &mgl[l].x, &mgl[l].b); */
            /*         smoother_dcsr_Schwarz_backward(&mgl[l].Schwarz, &swzparam, */
            /*                                        &mgl[l].x, &mgl[l].b); */
            /*         break; */
            /*     default: */
            /*         smoother_dcsr_Schwarz_forward(&mgl[l].Schwarz, &swzparam, */
            /*                                       &mgl[l].x, &mgl[l].b); */
            /*         break; */
            /* } */
        } else { // pre-smoothing with standard smoothers
          dcsr_fpresmoothing(smoother, &mgl[l].A, &mgl[l].b, &mgl[l].x, &mgl[l].M,
                             power, param->presmooth_iter, 0, mgl[l].A.row-1, 1,
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

      //#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
            // use UMFPACK direct solver on the coarsest level
            hazmath_solve(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, mgl[nl-1].Numeric, 0);
            break;
        }
	  //#endif
        default:
            // use eigensolver to approximate coarse A^-s
            coarse_fracinv(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, &mgl[nl-1].M, power);
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
        // alpha = pow(0.5, l);
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
 * \fn void famli (AMG_data *mgl, AMG_param *param, INT level)
 *
 * \brief Solve Ax=b with recursive AMLI-cycle with fractional power
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 * \param level  Current level
 *
 * \author Ana Budisa
 * \date   2020-05-08
 *
 * \note AMLI polynomial computed by the best approximation of 1/x.
 *       Refer to Johannes K. Kraus, Panayot S. Vassilevski, Ludmil T. Zikatanov,
 *       "Polynomial of best uniform approximation to $x^{-1}$ and smoothing in
 *        two-level methods", 2013.
 *
 */
void famli(AMG_data *mgl,
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
    const REAL   power = param->fpwr; // fractional exponent

    // local variables
    REAL   alpha  = 1.0;
    REAL * coef   = param->amli_coef;

    // fine level b, x and M
    dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x;
    // coarse level b, x and M
    dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x;

    dCSRmat *A0 = &mgl[level].A;   // fine level matrix
    dCSRmat *A1 = &mgl[level+1].A; // coarse level matrix

    dCSRmat *M0 = &mgl[level].M; // fine level mass mass matrix

    const INT m0 = A0->row, m1 = A1->row;

    REAL     *r        = mgl[level].w.val;      // work array for residual
    REAL     *r1       = mgl[level+1].w.val+m1; // work array for residual

    if ( prtlvl >= PRINT_MOST )
        printf("AMLI level %lld, smoother %lld.\n", (long long )level, (long long )smoother);

    if ( level < mgl[level].num_levels-1 ) {

        // presmoothing
        dcsr_fpresmoothing(smoother, A0, b0, e0, M0, power,
                           param->presmooth_iter, 0, m0-1, 1, relax);

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
        dcsr_fpostsmoothing(smoother, A0, b0, e0, M0, power,
                            param->postsmooth_iter, 0, m0-1, -1, relax);

    }

    else { // coarsest level solver

        switch (coarse_solver) {

	  //#if WITH_SUITESPARSE
            case SOLVER_UMFPACK:
                // use UMFPACK direct solver on the coarsest level //
                hazmath_solve(A0, b0, e0, mgl[level].Numeric, 0);
                break;
		//#endif

            default:
                /* use iterative solver on the coarsest level */
                coarse_fitsolver(A0, b0, e0, tol, prtlvl);

        }

    }

}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
