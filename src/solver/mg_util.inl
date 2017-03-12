/*! \file mg_util.inl
 *
 *  Routines for algebraic multigrid cycles
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 12/25/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note  Done cleanup for releasing -- Xiaozhe Hu 03/12/2017
 *
 *  \todo Combine pre- and post-smoothing and use a flag to make sure the symmetry whenever is necessary -- Xiaozhe Hu
 *
 */

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/
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
static void dcsr_presmoothing(const SHORT smoother,
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

        default:
            printf("### ERROR: Wrong smoother type %d!\n", smoother);
            check_error(ERROR_INPUT_PAR, __FUNCTION__);
    }
}

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
static void dcsr_postsmoothing(const SHORT smoother,
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

        default:
            printf("### ERROR: Wrong smoother type %d!\n", smoother);
            check_error(ERROR_INPUT_PAR, __FUNCTION__);
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
