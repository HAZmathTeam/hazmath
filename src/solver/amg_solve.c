/*! src/solver/amg_solve.c
 *
 *  Algebraic multigrid iterations: SOLVE phase.
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 12/24/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 * \note   Done cleanup for releasing -- Xiaozhe Hu 03/12/2017 & 08/27/2021
 *
 */

#include "hazmath.h"

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

//! Warning for residual false convergence
#define ITS_FACONV  printf("### HAZMATH WARNING: False convergence!\n")

//! Warning for solution close to zero
#define ITS_ZEROSOL printf("### HAZMATH WARNING: Iteration stopped due to the solution is almost zero! %s : %d\n", __FUNCTION__, __LINE__)

//! Warning for iteration restarted
#define ITS_RESTART printf("### HAZMATH WARNING: Iteration restarted due to stagnation! %s : %d\n", __FUNCTION__, __LINE__)

//! Warning for stagged iteration
#define ITS_STAGGED printf("### HAZMATH WARNING: Iteration stopped due to staggnation! %s : %d\n", __FUNCTION__, __LINE__)

//! Warning for tolerance practically close to zero
#define ITS_ZEROTOL printf("### HAZMATH WARNING: The tolerence might be too small! %s : %d\n", __FUNCTION__, __LINE__)

//! Warning for divided by zero
#define ITS_DIVZERO printf("### HAZMATH WARNING: Divided by zero! %s : %d\n", __FUNCTION__, __LINE__)

//! Warning for actual relative residual
#define ITS_REALRES(relres) printf("### HAZMATH WARNING: The actual relative residual = %e!\n",(relres))

//! Warning for computed relative residual
#define ITS_COMPRES(relres) printf("### HAZMATH WARNING: The computed relative residual = %e!\n",(relres))

//! Warning for too small sp
#define ITS_SMALLSP printf("### HAZMATH WARNING: sp is too small! %s : %d\n", __FUNCTION__, __LINE__)

//! Warning for restore previous iteration
#define ITS_RESTORE(iter) printf("### HAZMATH WARNING: Restore iteration %d!\n",(iter));

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
        printf("### HAZMATH WARNING: Max iter %d reached with rel. resid. %e.\n", MaxIt, relres);
    }
    else if ( iter >= 0 ) {
      printf("Num_iter(amg_solve.c) = %d with relative residual %e.\n", iter, relres);
    }
}

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/***********************************************************************************************/
/**
 * \fn INT amg_solve (AMG_data *mgl, AMG_param *param)
 *
 * \brief Solve phase for AMG method (as standard alone iterative solver)
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       Iteration number if converges; ERROR otherwise.
 *
 */
INT amg_solve(AMG_data *mgl,
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
        mgcycle(mgl, param);

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
        ITS_FINAL(iter, MaxIt, relres1);
        get_time(&solve_end);
        print_cputime("AMG solve",solve_end - solve_start);
    }


    return iter;
}

/***********************************************************************************************/
/**
 * \fn INT amg_solve_amli (AMG_data *mgl, AMG_param *param)
 *
 * \brief Solve phase for AMG method using AMLI-cycle (as standard alone iterative solver)
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   01/23/2011
 *
 * \note AMLI polynomial computed by the best approximation of 1/x.
 *       Refer to Johannes K. Kraus, Panayot S. Vassilevski,
 *       Ludmil T. Zikatanov, "Polynomial of best uniform approximation to $x^{-1}$
 *       and smoothing in two-level methods", 2013.
 *
 */
INT amg_solve_amli (AMG_data *mgl,
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
        amli(mgl, param, 0);

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
        ITS_FINAL(iter, MaxIt, relres1);
        get_time(&solve_end);
        solve_time = solve_end - solve_start;
        print_cputime("AMLI solve", solve_time);
    }

    return iter;
}

/***********************************************************************************************/
/**
 * \fn INT amg_solve_nl_amli(AMG_data *mgl, AMG_param *param)
 *
 * \brief Solve phase for AMG method using nonlinear AMLI-cycle (as standard alone iterative solver)
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   04/30/2011
 *
 * \note Nonlinear AMLI-cycle.
 *       Refer to Xiazhe Hu, Panayot S. Vassilevski, Jinchao Xu
 *       "Comparative Convergence Analysis of Nonlinear AMLI-cycle Multigrid", 2013.
 *
 */
INT amg_solve_nl_amli(AMG_data *mgl,
                      AMG_param *param)
{
    dCSRmat      *ptrA = &mgl[0].A;
    dvector      *b = &mgl[0].b, *x = &mgl[0].x, *r = &mgl[0].w;

    const INT     MaxIt  = param->maxit;
    const SHORT   prtlvl = param->print_level;
    const REAL    tol    = param->tol;
    const REAL    sumb   = dvec_norm2(b); // L2norm(b)

    // local variables
    REAL          solve_start, solve_end;
    REAL          relres1 = BIGREAL, absres0 = BIGREAL, absres, factor;
    INT           iter = 0;

    get_time(&solve_start);

    // Print iteration information if needed
    print_itsolver_info(prtlvl, STOP_REL_RES, iter, 1.0, sumb, 0.0);

    while ( (++iter <= MaxIt) & (sumb > SMALLREAL) ) // MG solver here
    {
        // Call nonlinear AMLI-cycle
        nl_amli(mgl, param, 0, mgl[0].num_levels);

        // Computer r = b-A*x
        dvec_cp(b, r);
        dcsr_aAxpy(-1.0, ptrA, x->val, r->val);

        absres  = dvec_norm2(r); // residual ||r||
        relres1 = absres/sumb;       // relative residual ||r||/||b||
        factor  = absres/absres0;    // contraction factor
        absres0 = absres;


        // output iteration information if needed
        print_itsolver_info(prtlvl, STOP_REL_RES, iter, relres1, absres, factor);

        if ( relres1 < tol ) break; // early exit condition
    }

    if ( prtlvl > PRINT_NONE ) {
        ITS_FINAL(iter, MaxIt, relres1);
        get_time(&solve_end);
        print_cputime("Nonlinear AMLI solve", solve_end - solve_start);
    }

    return iter;
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
