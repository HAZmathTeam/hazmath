/*! src/solver/amg_solve.c
 *
 *  Algebraic multigrid iterations with fractional smoothers: SOLVE phase.
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 12/24/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 * \note   Variation on amg_solve.c, but with functions from fmgcycle.c
 *         (Ana Budisa, 2020-05-14)
 *
 */

#include "hazmath.h"
/*! \file itsolver_util.inl
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 5/13/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *
 */

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn inline static void ITS_FINAL_ (const INT iter, const INT MaxIt, const REAL relres)
 * \brief Print out final status of an iterative method
 *
 * \param iter    Number of iterations
 * \param MaxIt   Maximal number of iterations
 * \param relres  Relative residual
 *
 */
inline static void ITS_FINAL_ (const INT iter, const INT MaxIt, const REAL relres)
{
    if ( iter > MaxIt ) {
        printf("### HAZMATH WARNING: Max iter %lld reached with rel. resid. %e.\n", (long long )MaxIt, relres);
    }
    else if ( iter >= 0 ) {
        printf("Number of iterations = %lld with relative residual %e.\n", (long long )iter, relres);
    }
}

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/**
 * \fn INT famg_solve (AMG_data *mgl, AMG_param *param)
 *
 * \brief Solve phase for AMG method (as standard alone iterative solver)
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       Iteration number if converges; ERROR otherwise.
 *
 */
INT famg_solve(AMG_data *mgl,
               AMG_param *param)
{
    dCSRmat      *ptrA = &mgl[0].A;
    dvector      *b = &mgl[0].b, *x = &mgl[0].x, *r = &mgl[0].w;

    const SHORT   prtlvl = param->print_level;
    const INT     MaxIt  = param->maxit;
    const REAL    tol    = param->tol;
    const REAL    sumb   = dvec_norm2(b);

    // local variables
    REAL  solve_start, solve_end;
    REAL  relres1 = BIGREAL, absres0 = sumb, absres, factor;
    INT   iter = 0;

    get_time(&solve_start);

    // Print iteration information if needed
    print_itsolver_info(prtlvl, STOP_REL_RES, iter, 1.0, sumb, 0.0);

    // Main loop
    while ( (++iter <= MaxIt) & (sumb > SMALLREAL) ) {

        // Call one multigrid cycle
        fmgcycle(mgl, param);

        // Form residual r = b - A*x
        dvec_cp(b, r);
        dcsr_aAxpy(-1.0, ptrA, x->val, r->val);

        // Compute norms of r and convergence factor
        absres  = dvec_norm2(r);
        relres1 = absres/sumb;
        factor  = absres/absres0;
        absres0 = absres;

        // Print iteration information if needed
        print_itsolver_info(prtlvl, STOP_REL_RES, iter, relres1, absres, factor);

        // Check convergence
        if ( relres1 < tol ) break;
    }

    if ( prtlvl > PRINT_NONE ) {
        ITS_FINAL_(iter, MaxIt, relres1);
        get_time(&solve_end);
        print_cputime("FAMG solve",solve_end - solve_start);
    }


    return iter;
}


/**
 * \fn INT amg_solve_famli (AMG_data *mgl, AMG_param *param)
 *
 * \brief Solve phase for AMG method using AMLI-cycle (as standard alone iterative solver)
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       Iteration number if converges; ERROR otherwise.
 *
 * \author Ana Budisa
 * \date   2020-05-14
 *
 * \note AMLI polynomial computed by the best approximation of 1/x.
 *       Refer to Johannes K. Kraus, Panayot S. Vassilevski,
 *       Ludmil T. Zikatanov, "Polynomial of best uniform approximation to $x^{-1}$
 *       and smoothing in two-level methods", 2013.
 *
 */
INT amg_solve_famli (AMG_data *mgl,
                     AMG_param *param)
{
    dCSRmat     *ptrA = &mgl[0].A;
    dvector     *b = &mgl[0].b, *x = &mgl[0].x, *r = &mgl[0].w;

    const INT    MaxIt  = param->maxit;
    const SHORT  prtlvl = param->print_level;
    const REAL   tol    = param->tol;
    const REAL   sumb   = dvec_norm2(b); // L2norm(b)

    // local variables
    REAL         solve_start, solve_end, solve_time;
    REAL         relres1 = BIGREAL, absres0 = sumb, absres, factor;
    INT          iter = 0;

    get_time(&solve_start);

    // Print iteration information if needed
    print_itsolver_info(prtlvl, STOP_REL_RES, iter, 1.0, sumb, 0.0);

    // MG solver here
    while ( (++iter <= MaxIt) & (sumb > SMALLREAL) ) {

        // Call one AMLI cycle
        famli(mgl, param, 0);

        // Form residual r = b-A*x
        dvec_cp(b, r);
        dcsr_aAxpy(-1.0, ptrA, x->val, r->val);

        // Compute norms of r and convergence factor
        absres  = dvec_norm2(r);
        relres1 = absres/sumb;
        factor  = absres/absres0;
        absres0 = absres;

        // Print iteration information if needed
        print_itsolver_info(prtlvl, STOP_REL_RES, iter, relres1, absres, factor);

        // Check convergence
        if ( relres1 < tol ) break;
    }

    if ( prtlvl > PRINT_NONE ) {
        ITS_FINAL_(iter, MaxIt, relres1);
        get_time(&solve_end);
        solve_time = solve_end - solve_start;
        print_cputime("FAMLI solve", solve_time);
    }

    return iter;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
