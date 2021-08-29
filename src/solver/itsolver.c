/*! \file src/solver/itsolver.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 10/06/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note  Done cleanup for releasing -- Xiaozhe Hu 03/12/2017 & 08/28/2021
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
        printf("### HAZMATH WARNING: Max iter %d reached with rel. resid. %e.\n", MaxIt, relres);
    }
    else if ( iter >= 0 ) {
      //        printf("Number of iterations = %d with relative residual %e.\n", iter, relres);
	printf("Num_iter(itsolver.c) = %d with relative residual %e.\n", iter, relres);
    }
}

/***********************************************************************************************/
/**
 * \fn inline static REAL16 frac_inv(REAL16 x, void *param)
 * \brief inverse of the
 *
 * \param iter    Number of iterations
 * \param MaxIt   Maximal number of iterations
 * \param relres  Relative residual
 *
 */
inline static REAL16 frac_inv(REAL16 x, void *param)
{
  REAL16 *s,s1,s2,alpha,beta;
  if(param!=NULL){
    s=(REAL16 *)param;
    s1=s[0];
    s2=s[1];
    alpha=s[2];
    beta=s[3];
  } else {
    s1=0.5e0;
    s2=-0.5e0;
    alpha=1e0;
    beta=2e0;
  }
  return 1./(alpha*powl(x,s1)+beta*powl(x,s2));
}
/**/

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/********************************************************************************************/
// general itrative solver for different format
/********************************************************************************************/
/**
 * \fn INT solver_dcsr_linear_itsolver (dCSRmat *A, dvector *b, dvector *x,
 *                                    precond *pc, linear_itsolver_param *itparam)
 *
 * \brief Solve Ax=b by preconditioned Krylov methods for CSR matrices
 *
 * \param A        Pointer to the coeff matrix in dCSRmat format
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param pc       Pointer to the preconditioning action
 * \param itparam  Pointer to parameters for lienar iterative solvers
 *
 * \return         Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   10/06/2015
 *
 * \note This is an abstract interface for iterative methods.
 */
INT solver_dcsr_linear_itsolver(dCSRmat *A,
                                dvector *b,
                                dvector *x,
                                precond *pc,
                                linear_itsolver_param *itparam)
{
    const SHORT prtlvl        = itparam->linear_print_level;
    const SHORT itsolver_type = itparam->linear_itsolver_type;
    const SHORT stop_type     = itparam->linear_stop_type;
    const SHORT restart       = itparam->linear_restart;
    const INT   MaxIt         = itparam->linear_maxit;
    const REAL  tol           = itparam->linear_tol;

    /* Local Variables */
    REAL solver_start, solver_end, solver_duration;
    INT iter;

    get_time(&solver_start);

    /* Safe-guard checks on parameters */
    ITS_CHECK ( MaxIt, tol );

    /* Choose a desirable Krylov iterative solver */
    switch ( itsolver_type ) {
        case 1:
            if ( prtlvl > PRINT_NONE ) {
                printf("**********************************************************\n");
                printf(" --> using Conjugate Gradient Method:\n");
            }
            iter = dcsr_pcg(A, b, x, pc, tol, MaxIt, stop_type, prtlvl);
            break;

        case 2:
            if ( prtlvl > PRINT_NONE ) {
                printf("**********************************************************\n");
                printf(" --> using MINRES Method:\n");
            }
            iter = dcsr_pminres(A, b, x, pc, tol, MaxIt, stop_type, prtlvl);
            printf(" NOTHING IMPLEMENTED FOR MINRES\n");
            break;

        case 3:
            if ( prtlvl > PRINT_NONE )  {
                printf("**********************************************************\n");
                printf(" --> using GMRES Method:\n");
            }
            iter = dcsr_pvgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, prtlvl);
            break;

        case 4:
            if ( prtlvl > PRINT_NONE ) {
                printf("**********************************************************\n");
                printf(" --> using Flexible GMRES Method:\n");
            }
            iter = dcsr_pvfgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, prtlvl);
            break;

        default:
            printf("### ERROR: Unknown itertive solver type %d!\n", itsolver_type);
            return ERROR_SOLVER_TYPE;

    }

    if ( (prtlvl >= PRINT_SOME) && (iter >= 0) ) {
        get_time(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Iterative method", solver_duration);
        printf("**********************************************************\n");
    }


    return iter;
}

/********************************************************************************************/
/**
 * \fn INT solver_bdcsr_linear_itsolver (block_dCSRmat *A, dvector *b, dvector *x,
 *                                     precond *pc, linear_itsolver_param *itparam)
 *
 * \brief Solve Ax = b by standard Krylov methods
 *
 * \param A        Pointer to the coeff matrix in block_dCSRmat format
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param pc       Pointer to the preconditioning action
 * \param itparam  Pointer to parameters for iterative solvers
 *
 * \return         Iteration number if converges; ERROR otherwise.
 *
 *
 * \author Xiaozhe Hu
 * \date   02/17/2016
 */
INT solver_bdcsr_linear_itsolver(block_dCSRmat *A,
                                 dvector *b,
                                 dvector *x,
                                 precond *pc,
                                 linear_itsolver_param *itparam)
{

    const SHORT prtlvl =        itparam->linear_print_level;
    const SHORT itsolver_type = itparam->linear_itsolver_type;
    const SHORT stop_type =     itparam->linear_stop_type;
    const SHORT restart =       itparam->linear_restart;
    const INT   MaxIt =         itparam->linear_maxit;
    const REAL  tol =           itparam->linear_tol;

    REAL  solver_start, solver_end, solver_duration;
    INT   iter = ERROR_SOLVER_TYPE;

    get_time(&solver_start);

    /* Safe-guard checks on parameters */
    ITS_CHECK ( MaxIt, tol );

    switch (itsolver_type) {

        case SOLVER_CG:
            if ( prtlvl > PRINT_NONE ) {
                printf("**********************************************************\n");
                printf(" --> using Conjugate Gradient Method (Block CSR):\n");
            }
            iter = bdcsr_pcg(A, b, x, pc, tol, MaxIt, stop_type, prtlvl);
            break;

        case SOLVER_MinRes:
            if ( prtlvl > PRINT_NONE ) {
                printf("**********************************************************\n");
                printf(" --> using MINRES Method (Block CSR):\n");
            }
            iter = bdcsr_pminres(A, b, x, pc, tol, MaxIt, stop_type, prtlvl);
            break;

        case SOLVER_VGMRES:
            if ( prtlvl > PRINT_NONE ) {
                printf("**********************************************************\n");
                printf(" --> using GMRES Method (Block CSR):\n");
            }
            iter = bdcsr_pvgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, prtlvl);
            break;

        case SOLVER_VFGMRES:
            if ( prtlvl > PRINT_NONE ) {
                printf("**********************************************************\n");
                printf(" --> using Flexible GMRES Method (Block CSR):\n");
            }
            iter = bdcsr_pvfgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, prtlvl);
            break;

        default:
            printf("### ERROR: Unknown itertive solver type %d!\n", itsolver_type);

    }

    if ( (prtlvl >= PRINT_MIN) && (iter >= 0) ) {
        get_time(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Iterative method", solver_duration);
        printf("**********************************************************\n");
    }

    return iter;
}

/********************************************************************************************/
/**
 * \fn INT solver_general_linear_itsolver (matvec *mxv, dvector *b, dvector *x,
 *                                    precond *pc, linear_itsolver_param *itparam)
 *
 * \brief Solve Ax=b by preconditioned Krylov methods (matrix-free version)
 *
 * \param mxv      Pointer to the function that performes matrix vector multiplication
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param pc       Pointer to the preconditioning action
 * \param itparam  Pointer to parameters for lienar iterative solvers
 *
 * \return         Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   05/31/2015
 *
 * \note This is an abstract interface for iterative methods.
 */
INT solver_general_linear_itsolver(matvec *mxv,
                                   dvector *b,
                                   dvector *x,
                                   precond *pc,
                                   linear_itsolver_param *itparam)
{
    const SHORT prtlvl        = itparam->linear_print_level;
    const SHORT itsolver_type = itparam->linear_itsolver_type;
    const SHORT stop_type     = itparam->linear_stop_type;
    const SHORT restart       = itparam->linear_restart;
    const INT   MaxIt         = itparam->linear_maxit;
    const REAL  tol           = itparam->linear_tol;

    /* Local Variables */
    REAL solver_start, solver_end, solver_duration;
    INT iter;

    get_time(&solver_start);

    /* Safe-guard checks on parameters */
    ITS_CHECK ( MaxIt, tol );

    /* Choose a desirable Krylov iterative solver */
    switch ( itsolver_type ) {
        case 1:
            if ( prtlvl > PRINT_NONE ) {
                printf("**********************************************************\n");
                printf(" --> using Conjugate Gradient Method:\n");
            }
            iter = general_pcg(mxv, b, x, pc, tol, MaxIt, stop_type, prtlvl);
            break;

        case 2:
            if ( prtlvl > PRINT_NONE ) {
                printf("**********************************************************\n");
                printf(" --> using MINRES Method:\n");
            }
            iter = general_pminres(mxv, b, x, pc, tol, MaxIt, stop_type, prtlvl);
            printf(" NOTHING IMPLEMENTED FOR MINRES\n");
            break;

        case 3:
            if ( prtlvl > PRINT_NONE )  {
                printf("**********************************************************\n");
                printf(" --> using GMRES Method:\n");
            }
            iter = general_pvgmres(mxv, b, x, pc, tol, MaxIt, restart, stop_type, prtlvl);
            break;

        case 4:
            if ( prtlvl > PRINT_NONE ) {
                printf("**********************************************************\n");
                printf(" --> using Flexible GMRES Method:\n");
            }

            iter = general_pvfgmres(mxv, b, x, pc, tol, MaxIt, restart, stop_type, prtlvl);
            break;

        default:
            printf("### ERROR: Unknown itertive solver type %d!\n", itsolver_type);
            return ERROR_SOLVER_TYPE;

    }

    if ( (prtlvl >= PRINT_SOME) && (iter >= 0) ) {
        get_time(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Iterative method", solver_duration);
        printf("**********************************************************\n");
    }


    return iter;
}


/********************************************************************************************/
// AMG method for CSR format
/********************************************************************************************/
/**
 * \fn void linear_solver_amg (dCSRmat *A, dvector *b, dvector *x,
 *                           AMG_param *param)
 *
 * \brief Solve Ax = b by algebraic multigrid methods
 *
 * \param A      Pointer to dCSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param x      Pointer to dvector: the unknowns
 * \param param  Pointer to AMG_param: AMG parameters
 *
 * \author Xiaozhe Hu
 * \date   12/25/2015
 *
 *
 */
INT linear_solver_amg(dCSRmat *A,
                      dvector *b,
                      dvector *x,
                      AMG_param *param)
{
    const SHORT   max_levels  = param->max_levels;
    const SHORT   prtlvl      = param->print_level;
    const SHORT   amg_type    = param->AMG_type;
    const SHORT   cycle_type  = param->cycle_type;
    const INT     nnz = A->nnz, m = A->row, n = A->col;

    // local variables
    INT      status = SUCCESS;
    AMG_data *    mgl = amg_data_create(max_levels);
    REAL          AMG_start, AMG_end;

    if ( prtlvl > PRINT_NONE ) get_time(&AMG_start);

    // check matrix data
    if ( m != n ) {
        printf("### ERROR: A is not a square matrix!\n");
        check_error(ERROR_MAT_SIZE, __FUNCTION__);
    }

    if ( nnz <= 0 ) {
        printf("### ERROR: A has no nonzero entries!\n");
        check_error(ERROR_MAT_SIZE, __FUNCTION__);
    }

    // Step 0: initialize mgl[0] with A, b and x
    mgl[0].A = dcsr_create(m, n, nnz);
    dcsr_cp(A, &mgl[0].A);

    mgl[0].b = dvec_create(n);
    dvec_cp(b, &mgl[0].b);

    mgl[0].x = dvec_create(n);
    dvec_cp(x, &mgl[0].x);

    // Step 1: AMG setup phase
    switch (amg_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl, param);
        break;

        case SA_AMG: // Smoothed Aggregation AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = amg_setup_sa(mgl, param);
        break;

        default: // Unsmoothed Aggregation AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl, param);
        break;

    }

    // Step 2: AMG solve phase
    if ( status == SUCCESS ) { // call a multilevel cycle

        switch (cycle_type) {

            case AMLI_CYCLE: // AMLI-cycle
                status = amg_solve_amli(mgl, param); break;

            case NL_AMLI_CYCLE: // Nonlinear AMLI-cycle
                status = amg_solve_nl_amli(mgl, param); break;

            default: // V,W-cycles (determined by param)
                status = amg_solve(mgl, param); break;

        }

        dvec_cp(&mgl[0].x, x);

    }

    else { // call a backup solver

        if ( prtlvl > PRINT_MIN ) {
            printf("### HAZMATH WARNING: AMG setup failed!\n");
            printf("### HAZMATH WARNING: Use a backup solver instead.\n");
        }
        status = dcsr_pvgmres(A, b, x, NULL, param->tol, param->maxit,
                              20, 1, prtlvl);

    }

    // clean-up memory
    amg_data_free(mgl, param);

    // print out CPU time if needed
    if ( prtlvl > PRINT_NONE ) {
        get_time(&AMG_end);
        print_cputime("AMG totally", AMG_end - AMG_start);
        printf("**********************************************************\n");
    }

    return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_frac_rational_approx (dCSRmat *A, dvector *b, dvector *x, dCSRmat *M,
 *                                             AMG_param *amgparam, linear_itsolver_param *itparam,
 *                                             dvector *poles, devector *residues)
 *
 * \brief Solve (alpha*A^s + beta*A^t)x = b by rational approximations
 *
 * \note fractional matrix A^s is defined as A^s = M U Lambda^s U^T M, where A U = M U Lambda
 * \note rational approximation to (alpha*z^x + beta*z^t)^{-1} should have been computed already
 *
 * \param A      Pointer to dCSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param x      Pointer to dvector: the unknowns
 * \param M      Pointer to dCSRmat: the mass matrix
 * \param param  Pointer to AMG_param: AMG parameters
 * \param
 * \param poles  Pointer to dvector: poles
 *
 * \note If M = NULL, we assume that M = I  (this does not work for python ... - Xiaozhe)
 *
 * \author Xiaozhe Hu
 * \date   2020-09-27
 *
 * \note revised by Xiaozhe Hu on 12/26/2020
 *
 *
 */
INT linear_solver_frac_rational_approx(dCSRmat *A,
                                       dvector *b,
                                       dvector *x,
                                       dCSRmat *M,
                                       AMG_param *amgparam,
                                       linear_itsolver_param *itparam,
                                       dvector *poles,
                                       dvector *residues)
{

  // local variables
  INT k = poles->row;
  INT status = SUCCESS;
  INT i;
  dvector update = dvec_create(x->row);
  dCSRmat shiftA;
  dCSRmat I;

  if (M=NULL)
  {
    I = dcsr_create_identity_matrix(A->row, 0);
  }

  // apply rational approximation
  // x = residues(0)*(M\b)
  if (M=NULL)
  {
    dvec_cp(b, x);
  }
  else
  {
    status = linear_solver_dcsr_krylov_amg(M, b, x, itparam, amgparam);
  }
  dvec_ax(residues->val[0], x);

  // main loop
  INT count = 0;
  for (i=0; i<k; i++)
  {
    // set update to zero
    dvec_set(update.row, &update, 0.0);

    // form (A - poles[i]*M)
    if (M=NULL)
    {
      dcsr_add(A, 1.0, &I, -poles->val[i], &shiftA);
    }
    else
    {
      dcsr_add(A, 1.0, M, -poles->val[i], &shiftA);
    }

    // solve (A - poles[i]*M)e=f
    status = linear_solver_dcsr_krylov_amg(&shiftA, b, &update, itparam, amgparam);

    // x = x + residues[i+1]*update
    dvec_axpy(residues->val[i+1], &update, x);

    // free shiftA
    dcsr_free(&shiftA);

    // iter
    if(status > 0) count += status;
  }

  // cleanup
  dvec_free(&update);
  if (M=NULL) dcsr_free(&I);

  return count;

}

/********************************************************************************************/
// preconditioned Krylov methods for CSR format
/********************************************************************************************/
/**
 * \fn INT linear_solver_dcsr_krylov (dCSRmat *A, dvector *b, dvector *x,
 *                                  linear_itsolver_param *itparam)
 *
 * \brief Solve Ax=b by standard Krylov methods for CSR matrices
 *
 * \param A        Pointer to the coeff matrix in dCSRmat format
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param itparam  Pointer to parameters for linear iterative solvers
 *
 * \return         Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   10/06/2015
 */
INT linear_solver_dcsr_krylov(dCSRmat *A,
                             dvector *b,
                             dvector *x,
                             linear_itsolver_param *itparam)
{
    const SHORT prtlvl = itparam->linear_print_level;

    /* Local Variables */
    INT      status = SUCCESS;
    REAL     solver_start, solver_end, solver_duration;

    get_time(&solver_start);

    status = solver_dcsr_linear_itsolver(A,b,x,NULL,itparam);

    if ( prtlvl >= PRINT_MIN ) {
        get_time(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Krylov method totally", solver_duration);
        printf("**********************************************************\n");
    }

    return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_dcsr_krylov_diag (dCSRmat *A, dvector *b, dvector *x,
 *                                       linear_itsolver_param *itparam)
 *
 * \brief Solve Ax=b by diagonal preconditioned Krylov methods
 *
 * \param A        Pointer to the coeff matrix in dCSRmat format
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param itparam  Pointer to parameters for linear iterative solvers
 *
 * \return         Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   10/06/2015
 */
INT linear_solver_dcsr_krylov_diag(dCSRmat *A,
                                   dvector *b,
                                   dvector *x,
                                   linear_itsolver_param *itparam)
{
    const SHORT prtlvl = itparam->linear_print_level;

    /* Local Variables */
    INT       status = SUCCESS;
    REAL      solver_start, solver_end, solver_duration;

    get_time(&solver_start);

    // setup preconditioner
    dvector diag; dcsr_getdiag(0,A,&diag);

    precond pc;
    pc.data = &diag;
    pc.fct  = precond_diag;

    // call iterative solver
    status = solver_dcsr_linear_itsolver(A,b,x,&pc,itparam);

    if ( prtlvl >= PRINT_MIN ) {
        get_time(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Diag_Krylov method totally", solver_duration);
        printf("**********************************************************\n");

    }

    dvec_free(&diag);

    return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_dcsr_krylov_amg (dCSRmat *A, dvector *b, dvector *x,
 *                                      linear_itsolver_param *itparam, AMG_param *amgparam)
 *
 * \brief Solve Ax=b by AMG preconditioned Krylov methods
 *
 * \param A         Pointer to the coeff matrix in dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG methods
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   10/06/2015
 */
INT linear_solver_dcsr_krylov_amg(dCSRmat *A,
                                  dvector *b,
                                  dvector *x,
                                  linear_itsolver_param *itparam,
                                  AMG_param *amgparam)
{
    const SHORT prtlvl = itparam->linear_print_level;
    const SHORT max_levels = amgparam->max_levels;
    const INT nnz = A->nnz, m = A->row, n = A->col;

    /* Local Variables */
    INT      status = SUCCESS;
    REAL     solver_start, solver_end, solver_duration;

    get_time(&solver_start);

    // initialize A, b, x for mgl[0]
    AMG_data *mgl=amg_data_create(max_levels);
    mgl[0].A=dcsr_create(m,n,nnz); dcsr_cp(A,&mgl[0].A);
    mgl[0].b=dvec_create(n); mgl[0].x=dvec_create(n);

    // setup preconditioner
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl, amgparam);
        break;

        case SA_AMG: // Smoothed Aggregation AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = amg_setup_sa(mgl, amgparam);
        break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl, amgparam);
        break;

    }

    if (status < 0) goto FINISHED;

    // setup preconditioner
    precond_data pcdata;
    param_amg_to_prec(&pcdata,amgparam);
    pcdata.max_levels = mgl[0].num_levels;
    pcdata.mgl_data = mgl;

    precond pc; pc.data = &pcdata;

    switch (amgparam->cycle_type) {

        case AMLI_CYCLE: // AMLI cycle
            pc.fct = precond_amli;
            break;

        case NL_AMLI_CYCLE: // Nonlinear AMLI AMG
            pc.fct = precond_nl_amli;
            break;

        case ADD_CYCLE: // additive cycle
            pc.fct = precond_amg_add;
            break;

        default: // V,W-Cycle AMG
            pc.fct = precond_amg;
            break;

    }

    // call iterative solver
    status = solver_dcsr_linear_itsolver(A, b, x, &pc, itparam);

    if ( prtlvl >= PRINT_MIN ) {
        get_time(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("AMG_Krylov method totally", solver_duration);
        fprintf(stdout,"**********************************************************\n");
    }

FINISHED:
    amg_data_free(mgl, amgparam);free(mgl);
    return status;
}


/********************************************************************************************/
/**
 * \fn INT linear_solver_dcsr_krylov_famg (dCSRmat *A_frac, dvector *bb, dvector *x, dCSRmat *M, dCSRmat *A,
 *                                      linear_itsolver_param *itparam, AMG_param *amgparam)
 *
 * \brief Solve Ax=b by AMG preconditioned Krylov methods
 *
 * \param A_frac    Pointer to the coeff matrix in dCSRmat format
 * \param bb         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param M         Pointer to the mass matrix diagonal in dCSRmat format
 * \param A         Pointer to the stiff matrix in dCSRmat format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG methods
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Ana Budisa
 * \date   2020-05-14  //
 *
 */
INT linear_solver_dcsr_krylov_famg(dCSRmat *A_frac,
                                   dvector *bb,
                                   dvector *x,
                                   dCSRmat *M,
                                   dCSRmat *A,
                                   linear_itsolver_param *itparam,
                                   AMG_param *amgparam)
{
    const SHORT prtlvl = itparam->linear_print_level;
    const SHORT max_levels = amgparam->max_levels;
    const INT nnz_A = A->nnz, m_A = A->row, n_A = A->col;
    const INT nnz_M = M->nnz, m_M = M->row, n_M = M->col;

    /* Local Variables */
    INT      status = SUCCESS;
    REAL     solver_start, solver_end, solver_duration;

    get_time(&solver_start);

    // initialize A, b, x, M for mgl[0]
    AMG_data *mgl = amg_data_create(max_levels);

    mgl[0].A      = dcsr_create(m_A, n_A, nnz_A); dcsr_cp(A, &mgl[0].A);

    mgl[0].M      = dcsr_create(m_M, n_M, nnz_M); dcsr_cp(M, &mgl[0].M);

    mgl[0].b      = dvec_create(m_A);
    mgl[0].x      = dvec_create(n_A);
    // randomize input
    // dvec_rand_true(n_A, &mgl[0].x);

    // setup preconditioner
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = famg_setup_ua(mgl, amgparam);
            printf("AMG status: %d \n", status);
        break;

        case SA_AMG: // Smoothed Aggregation AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = famg_setup_sa(mgl, amgparam);
            printf("AMG status: %d \n", status);
        break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = famg_setup_ua(mgl, amgparam);
            printf("AMG status: %d \n", status);
        break;

    }

    if (status < 0) goto FINISHED;

    // setup preconditioner
    precond_data pcdata;
    param_amg_to_prec(&pcdata, amgparam);
    pcdata.max_levels = mgl[0].num_levels;
    pcdata.mgl_data = mgl;

    precond pc; pc.data = &pcdata;

    switch (amgparam->cycle_type) {

        case AMLI_CYCLE: // AMLI cycle
            pc.fct = precond_famli;
        break;

        default: // V,W-Cycle AMG
            pc.fct = precond_famg_add;
        break;

    }

    // call iterative solver (status <=> iter)
    status = solver_dcsr_linear_itsolver(A_frac, bb, x, &pc, itparam);

    if ( prtlvl >= PRINT_MIN ) {
        get_time(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("AMG_Krylov method totally", solver_duration);
        printf("**********************************************************\n");
    }

FINISHED:
    amg_data_free(mgl, amgparam);
    return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_dcsr_krylov_famg (dCSRmat *A_frac, dvector *bb, dvector *x, dCSRmat *M, dCSRmat *A,
 *                                      linear_itsolver_param *itparam, AMG_param *amgparam)
 *
 * \brief Solve Ax=b by AMG preconditioned Krylov methods
 *
 * \param A_frac    Pointer to the coeff matrix in dCSRmat format
 * \param bb         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param M         Pointer to the mass matrix diagonal in dCSRmat format
 * \param A         Pointer to the stiff matrix in dCSRmat format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG methods
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Ana Budisa
 * \date   2020-07-07  //
 *
 */
INT linear_solver_dcsr_krylov_famg2(dCSRmat *A_frac,
                                    dvector *bb,
                                    dvector *x,
                                    dCSRmat *M,
                                    dCSRmat *A,
                                    dCSRmat *Grad,
                                    linear_itsolver_param *itparam,
                                    AMG_param *amgparam)
{
    const SHORT prtlvl = itparam->linear_print_level;
    const SHORT max_levels = amgparam->max_levels;
    const INT nnz_A = A->nnz, m_A = A->row, n_A = A->col;
    const INT nnz_M = M->nnz, m_M = M->row, n_M = M->col;
    const INT n_Af = A_frac->col;

    /* Local Variables */
    INT      status = SUCCESS;
    REAL     solver_start, solver_end, solver_duration;

    get_time(&solver_start);

    // initialize A, b, x, M for mgl[0]
    AMG_data *mgl = amg_data_create(max_levels);

    mgl[0].A      = dcsr_create(m_A, n_A, nnz_A); dcsr_cp(A, &mgl[0].A);

    mgl[0].M      = dcsr_create(m_M, n_M, nnz_M); dcsr_cp(M, &mgl[0].M);

    mgl[0].b      = dvec_create(m_A);
    mgl[0].x      = dvec_create(n_A);

    // randomize input (but orthogonal to constants!)
    dvec_rand_true(n_Af, x);
    dvec_orthog_const(x);

    // setup preconditioner
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = famg_setup_ua(mgl, amgparam);
            printf("AMG status: %d \n", status);
        break;

        case SA_AMG: // Smoothed Aggregation AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = famg_setup_sa(mgl, amgparam);
            printf("AMG status: %d \n", status);
        break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = famg_setup_ua(mgl, amgparam);
            printf("AMG status: %d \n", status);
        break;

    }

    if (status < 0) goto FINISHED;

    // setup preconditioners: 0 - FAMG(Adiv^s),
    //                        1 - Grad
    //                        2 - Grad^T
    precond_data pcdata0, pcdata1, pcdata2;

    param_amg_to_prec(&pcdata0, amgparam);
    pcdata0.max_levels = mgl[0].num_levels;
    pcdata0.mgl_data = mgl;

    precond_data_null(&pcdata1); precond_data_null(&pcdata2);
    pcdata1.A = dcsr_create_p(Grad->row, Grad->col, Grad->nnz);
    dcsr_cp(Grad, pcdata1.A);
    pcdata2.A = dcsr_create_p(Grad->col, Grad->row, Grad->nnz);
    dcsr_trans(Grad, pcdata2.A);

    // INT i;
    // for(i = 0; i < Grad->nnz; ++i) printf("%.5f \t", Grad->val[i]);
    // printf("\n");

    // array of preconditioners
    precond pc;
    pc.data = (precond_data*)malloc(3*sizeof(precond_data));
    ((precond_data*)pc.data)[0] = pcdata0; ((precond_data*)pc.data)[1] = pcdata1;
    ((precond_data*)pc.data)[2] = pcdata2;

    pc.fct = precond_famg_add2;

    // call iterative solver (status <=> iter)
    status = solver_dcsr_linear_itsolver(A_frac, bb, x, &pc, itparam);

    if ( prtlvl >= PRINT_MIN ) {
        get_time(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("AMG_Krylov method totally", solver_duration);
        printf("**********************************************************\n");
    }

FINISHED:
    amg_data_free(mgl, amgparam);
    // precond_free(&pc)??
    return status;
}


/********************************************************************************************/
/**
 * \fn INT linear_solver_dcsr_krylov_famg_sum (dCSRmat *A_frac, dvector *bb, dvector *x, dCSRmat *M, dCSRmat *A,
 *                                             const REAL falpha, const REAL beta, linear_itsolver_param *itparam, AMG_param *amgparam)
 *
 * \brief Solve A_frac x = b by FAMG preconditioned Krylov methods where
 *              A_frac = falpha * A^s + fbeta * A^(1+s), s in (-1, 0)
 *
 * \param A_frac    Pointer to the coeff matrix in dCSRmat format
 * \param bb        Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param M         Pointer to the mass matrix diagonal in dCSRmat format
 * \param A         Pointer to the stiff matrix in dCSRmat format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG methods
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Ana Budisa
 * \date   2020-06-22  //
 *
 */
INT linear_solver_dcsr_krylov_famg_sum(dCSRmat *A_frac,
                                       dvector *bb,
                                       dvector *x,
                                       dCSRmat *M,
                                       dCSRmat *A,
                                       const REAL falpha,
                                       const REAL fbeta,
                                       linear_itsolver_param *itparam,
                                       AMG_param *amgparam,
                                       AMG_param *famgparam)
{
    const SHORT prtlvl = itparam->linear_print_level;
    const SHORT max_levels_famg = famgparam->max_levels;
    const SHORT max_levels_amg = amgparam->max_levels;
    const INT nnz_A = A->nnz, m_A = A->row, n_A = A->col;
    const INT nnz_M = M->nnz, m_M = M->row, n_M = M->col;

    /* Local Variables */
    INT      status = SUCCESS;
    REAL     solver_start, solver_end, solver_duration;

    get_time(&solver_start);

    // initialize A, b, x, M for fractional mgl[0]
    AMG_data *fmgl = amg_data_create(max_levels_famg);

    fmgl[0].A      = dcsr_create(m_A, n_A, nnz_A); dcsr_cp(A, &fmgl[0].A);
    fmgl[0].M      = dcsr_create(m_M, n_M, nnz_M); dcsr_cp(M, &fmgl[0].M);
    fmgl[0].b      = dvec_create(m_A);
    fmgl[0].x      = dvec_create(n_A);
    // randomize input
    // dvec_rand_true(n_A, &mgl[0].x);

    // initialize A, b, x for standard mgl[0]
    AMG_data *mgl  = amg_data_create(max_levels_amg);

    mgl[0].A      = dcsr_create(m_A, n_A, nnz_A);
    mgl[0].b      = dvec_create(m_A);
    mgl[0].x      = dvec_create(n_A);

    // create alpha * lump(M)^-1 + beta * lump(M)^-1 A lump(M)^-1
    // lump(M)^-1
    dvector M_lump = dvec_create(m_M), x1 = dvec_create(n_M);
    dvec_set(n_M, &x1, 1.0);
    dcsr_mxv(M, x1.val, M_lump.val);
    dvec_inv(&M_lump); // M^-1 as vector
    dCSRmat M_lump_mat= dcsr_create_diagonal_matrix(&M_lump); // M^-1 as matrix

    // make MAM
    dCSRmat MAM = dcsr_create(m_A, n_A, nnz_A), MAM_t;
    dcsr_cp(A, &MAM); // MAM = A
    dcsr_row_scale(&MAM, &M_lump); // MAM = M^-1 A
    dcsr_trans(&MAM, &MAM_t); // MAM_t = A^T M^-1
    dcsr_row_scale(&MAM_t, &M_lump); // MAM_t = M^-1 A^T M^-1
    dcsr_free(&MAM); // clean MAM TODO: is this necessary?? maybe realloc?
    dcsr_trans(&MAM_t, &MAM); // MAM = M^-1 A M^-1

    // make alpha * lump(M)^-1 + beta * lump(M)^-1 A lump(M)^-1
    dcsr_add(&M_lump_mat, falpha, &MAM, fbeta, &mgl[0].A); // save matrix as finest amg level A

    // setup AMG
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl, amgparam);
            printf("AMG status: %d \n", status);

        break;

        case SA_AMG: // Smoothed Aggregation AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = amg_setup_sa(mgl, amgparam);
            printf("AMG status: %d \n", status);
        break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl, amgparam);
            printf("AMG status: %d \n", status);

        break;

    }

    // setup FAMG
    switch (famgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA FAMG ...\n");
            status = famg_setup_ua(fmgl, famgparam);
            printf("FAMG status: %d \n", status);
        break;

        case SA_AMG: // Smoothed Aggregation AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = famg_setup_sa(fmgl, famgparam);
            printf("FAMG status: %d \n", status);
        break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA FAMG ...\n");
            status = famg_setup_ua(fmgl, famgparam);
            printf("FAMG status: %d \n", status);
        break;

    }

    if (status < 0) goto FINISHED;

    // setup preconditioners: 1 - FAMG(s/2),
    //                        2 - AMG(alpha * lump(M)^-1 + beta * lump(M)^-1 A lump(M)^-1)
    precond_data pcdata1, pcdata2;

    param_amg_to_prec(&pcdata1, famgparam);
    param_amg_to_prec(&pcdata2, amgparam);

    pcdata1.max_levels = fmgl[0].num_levels;
    pcdata1.mgl_data = fmgl;

    pcdata2.max_levels = mgl[0].num_levels;
    pcdata2.mgl_data = mgl;

    // array of preconditioners
    precond pc;
    pc.data = (precond_data*)malloc(2*sizeof(precond_data));
    ((precond_data*)pc.data)[0] = pcdata1; ((precond_data*)pc.data)[1] = pcdata2;
    //pc.data = &pcdata;
    pc.fct = precond_sum_famg_add;

    // call iterative solver (status <=> iter)
    status = solver_dcsr_linear_itsolver(A_frac, bb, x, &pc, itparam);

    if ( prtlvl >= PRINT_MIN ) {
        get_time(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("AMG_Krylov method totally", solver_duration);
        printf("**********************************************************\n");
    }

FINISHED:
    amg_data_free(fmgl, famgparam);
    amg_data_free(mgl, amgparam);
    dcsr_free(&MAM); dcsr_free(&MAM_t); dcsr_free(&M_lump_mat);
    dvec_free(&M_lump); dvec_free(&x1);
    // precond_free(&pc); ??

    fprintf(stdout,"------------------- Leaving hazmath \n"); fflush(stdout);

    return status;
}


/********************************************************************************************/
/**
 * \fn INT linear_solver_dcsr_krylov_famg_sum2 (dCSRmat *A_frac, dvector *bb, dvector *x, dCSRmat *M, dCSRmat *A,
 *                                             const REAL falpha, const REAL beta, linear_itsolver_param *itparam, AMG_param *amgparam)
 *
 * \brief Solve A_frac x = b by FAMG preconditioned Krylov methods where
 *              A_frac = falpha * A^s + fbeta * A^(1+s), s in (-1, 0)
 *
 * \param A_frac    Pointer to the coeff matrix in dCSRmat format
 * \param bb        Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param M         Pointer to the mass matrix diagonal in dCSRmat format
 * \param A         Pointer to the stiff matrix in dCSRmat format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG methods
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Ana Budisa
 * \date   2020-07-07  //
 *
 */
INT linear_solver_dcsr_krylov_famg_sum2(dCSRmat *A_frac,
                                        dvector *bb,
                                        dvector *x,
                                        dCSRmat *MS,
                                        dCSRmat *AS,
                                        dCSRmat *Mdiv,
                                        dCSRmat *Adiv,
                                        dCSRmat *Adivfrac,
                                        dCSRmat *Grad,
                                        const REAL falpha,
                                        const REAL fbeta,
                                        linear_itsolver_param *itparam,
                                        AMG_param *amgparam,
                                        AMG_param *famgparam)
{
    const SHORT prtlvl = itparam->linear_print_level;
    const SHORT max_levels_famg = famgparam->max_levels;
    const SHORT max_levels_amg = amgparam->max_levels;

    const INT nnz_A = AS->nnz, m_A = AS->row, n_A = AS->col;
    //not used:    const INT nnz_M = MS->nnz,
    const INT m_M = MS->row, n_M = MS->col;

    const INT nnz_Adiv = Adiv->nnz, m_Adiv = Adiv->row, n_Adiv = Adiv->col;
    const INT nnz_Mdiv = Mdiv->nnz, m_Mdiv = Mdiv->row, n_Mdiv = Mdiv->col;

    /* Local Variables */
    INT      status = SUCCESS;
    REAL     solver_start, solver_end, solver_duration;

    get_time(&solver_start);

    // initialize A, b, x, M for fractional mgl[0]
    AMG_data *fmgl = amg_data_create(max_levels_famg);

    fmgl[0].A      = dcsr_create(m_Adiv, n_Adiv, nnz_Adiv); dcsr_cp(Adiv, &fmgl[0].A);
    fmgl[0].M      = dcsr_create(m_Mdiv, n_Mdiv, nnz_Mdiv); dcsr_cp(Mdiv, &fmgl[0].M);
    fmgl[0].b      = dvec_create(m_Adiv);
    fmgl[0].x      = dvec_create(n_Adiv);
    fmgl[0].Numeric = umfpack_factorize(Adivfrac, prtlvl); // LU factorization of Adiv^1+s/2 for direct solve

    // initialize A, b, x for standard mgl[0]
    AMG_data *mgl  = amg_data_create(max_levels_amg);

    mgl[0].A      = dcsr_create(m_A, n_A, nnz_A);
    mgl[0].b      = dvec_create(m_A);
    mgl[0].x      = dvec_create(n_A);

    // create alpha * lump(M)^-1 + beta * lump(M)^-1 A lump(M)^-1
    // lump(M)^-1
    dvector M_lump = dvec_create(m_M), x1 = dvec_create(n_M);
    dvec_set(n_M, &x1, 1.0);
    dcsr_mxv(MS, x1.val, M_lump.val);
    dvec_inv(&M_lump); // M^-1 as vector
    dCSRmat M_lump_mat= dcsr_create_diagonal_matrix(&M_lump); // M^-1 as matrix

    // make MAM
    dCSRmat MAM = dcsr_create(m_A, n_A, nnz_A), MAM_t;
    dcsr_cp(AS, &MAM); // MAM = A
    dcsr_row_scale(&MAM, &M_lump); // MAM = M^-1 A
    dcsr_trans(&MAM, &MAM_t); // MAM_t = A^T M^-1
    dcsr_row_scale(&MAM_t, &M_lump); // MAM_t = M^-1 A^T M^-1
    dcsr_free(&MAM); // clean MAM TODO: is this necessary?? maybe realloc?
    dcsr_trans(&MAM_t, &MAM); // MAM = M^-1 A M^-1

    // make alpha * lump(M)^-1 + beta * lump(M)^-1 A lump(M)^-1
    dcsr_add(&M_lump_mat, falpha, &MAM, fbeta, &mgl[0].A); // save matrix as finest amg level A

    // setup AMG
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl, amgparam);
            printf("AMG status: %d \n", status);

        break;

        case SA_AMG: // Smoothed Aggregation AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = amg_setup_sa(mgl, amgparam);
            printf("AMG status: %d \n", status);
        break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl, amgparam);
            printf("AMG status: %d \n", status);

        break;

    }

    // setup FAMG
    switch (famgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA FAMG ...\n");
            status = famg_setup_ua(fmgl, famgparam);
            printf("FAMG status: %d \n", status);
        break;

        case SA_AMG: // Smoothed Aggregation AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = famg_setup_sa(fmgl, famgparam);
            printf("FAMG status: %d \n", status);
        break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA FAMG ...\n");
            status = famg_setup_ua(fmgl, famgparam);
            printf("FAMG status: %d \n", status);
        break;

    }

    if (status < 0) goto FINISHED;

    // setup preconditioners: 0-2 - FAMG(Adiv^1+s/2),
    //                        3   - AMG(alpha * lump(M)^-1 + beta * lump(M)^-1 A lump(M)^-1)
    precond_data pcdata0, pcdata1, pcdata2, pcdata3;

    param_amg_to_prec(&pcdata0, famgparam);
    pcdata0.max_levels = fmgl[0].num_levels;
    pcdata0.mgl_data = fmgl;
    dcsr_alloc(Adivfrac->row, Adivfrac->col, Adivfrac->nnz, pcdata0.A);
    dcsr_cp(Adivfrac, pcdata0.A);// for direct solve!

    precond_data_null(&pcdata1); precond_data_null(&pcdata2);
    pcdata1.A = dcsr_create_p(Grad->row, Grad->col, Grad->nnz);
    dcsr_cp(Grad, pcdata1.A);
    pcdata2.A = dcsr_create_p(Grad->col, Grad->row, Grad->nnz);
    dcsr_trans(Grad, pcdata2.A);

    param_amg_to_prec(&pcdata3, amgparam);
    pcdata3.max_levels = mgl[0].num_levels;
    pcdata3.mgl_data = mgl;

    // array of preconditioners
    precond pc;
    pc.data = (precond_data*)malloc(4*sizeof(precond_data));
    ((precond_data*)pc.data)[0] = pcdata0; ((precond_data*)pc.data)[1] = pcdata1;
    ((precond_data*)pc.data)[2] = pcdata2; ((precond_data*)pc.data)[3] = pcdata3;
    //pc.data = &pcdata;
    pc.fct = precond_sum_famg_add2;

    // randomize first guess
    //dvec_rand_true(A_frac->col, x);
    //dvec_orthog_const(x);
    // set first guess to zero
    dvec_set(A_frac->col, x, 0.0);

    // call iterative solver (status <=> iter)
    status = solver_dcsr_linear_itsolver(A_frac, bb, x, &pc, itparam);

    if ( prtlvl >= PRINT_MIN ) {
        get_time(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("AMG_Krylov method totally", solver_duration);
        printf("**********************************************************\n");
    }

FINISHED:
    amg_data_free(fmgl, famgparam);
    amg_data_free(mgl, amgparam);
    dcsr_free(&MAM); dcsr_free(&MAM_t); dcsr_free(&M_lump_mat);
    dvec_free(&M_lump); dvec_free(&x1);
    // precond_free(&pc); ??

    fprintf(stdout,"------------------- Leaving hazmath \n"); fflush(stdout);

    return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_dcsr_krylov_hx_curl (dCSRmat *A, dvector *b, dvector *x,
 *                                      linear_itsolver_param *itparam, AMG_param *amgparam,
 *                                      dCSRmat P_curl, dCSRmat Grad)
 *
 * \brief Solve Ax=b by HX (H(curl)) preconditioned Krylov methods
 *
 * \param A         Pointer to the coeff matrix in dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG methods
 * \param P_curl    Pointer to the Pi_curl interpolation in dCSRmat format
 * \param Grad      Pointer to the Grad operator in dCSRmat format
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   02/10/2016
 */
INT linear_solver_dcsr_krylov_hx_curl(dCSRmat *A,
                                      dvector *b,
                                      dvector *x,
                                      linear_itsolver_param *itparam,
                                      AMG_param *amgparam,
                                      dCSRmat *P_curl,
                                      dCSRmat *Grad)
{
    const SHORT prtlvl = itparam->linear_print_level;
    const SHORT max_levels = amgparam->max_levels;

    /*------------------------*/
    /* Local Variables */
    /*------------------------*/
    INT      status = SUCCESS;
    REAL     solver_start, solver_end, solver_duration;

    get_time(&solver_start);

    /*------------------------*/
    /* setup vector Laplacian */
    /*------------------------*/

    // get transpose of P_curl
    dCSRmat Pt_curl;
    dcsr_trans(P_curl, &Pt_curl);

    // get A_vgrad
    dCSRmat A_vgrad;
    dcsr_rap(&Pt_curl, A, P_curl, &A_vgrad);

    // initialize A, b, x for mgl_vgrad[0]
    AMG_data *mgl_vgrad = amg_data_create(max_levels);
    mgl_vgrad[0].A=dcsr_create(A_vgrad.row,A_vgrad.col,A_vgrad.nnz);
    dcsr_cp(&A_vgrad, &mgl_vgrad[0].A);
    mgl_vgrad[0].b=dvec_create(A_vgrad.col);
    mgl_vgrad[0].x=dvec_create(A_vgrad.row);

    // setup AMG for vector Laplacian
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_vgrad, amgparam); break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = amg_setup_sa(mgl_vgrad, amgparam); break;

        default: // Classical AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_vgrad, amgparam); break;

    }

    if (status < 0) goto FINISHED;

    /*------------------------*/
    /* setup scalar Laplacian */
    /*------------------------*/
    // get transpose of Grad
    dCSRmat Gradt;
    dcsr_trans(Grad, &Gradt);

    // get A_grad
    dCSRmat A_grad;
    dcsr_rap(&Gradt, A, Grad, &A_grad);

    // initialize A, b, x for mgl_grad[0]
    AMG_data *mgl_grad = amg_data_create(max_levels);
    mgl_grad[0].A=dcsr_create(A_grad.row,A_grad.col,A_grad.nnz);
    dcsr_cp(&A_grad, &mgl_grad[0].A);
    mgl_grad[0].b=dvec_create(A_grad.col);
    mgl_grad[0].x=dvec_create(A_grad.row);

    // setup AMG for scalar Laplacian
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_grad, amgparam);
        break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = amg_setup_ua(mgl_grad, amgparam);
        break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_grad, amgparam);
        break;

    }

    if (status < 0) goto FINISHED;

    /*------------------------*/
    // setup preconditioner
    HX_curl_data hxcurldata;

    hxcurldata.A = A;

    hxcurldata.smooth_type = 1;
    hxcurldata.smooth_iter = itparam->HX_smooth_iter;

    hxcurldata.P_curl = P_curl;
    hxcurldata.Pt_curl = &Pt_curl;
    hxcurldata.A_vgrad = &A_vgrad;
    hxcurldata.amgparam_vgrad = amgparam;
    hxcurldata.mgl_vgrad = mgl_vgrad;

    hxcurldata.Grad = Grad;
    hxcurldata.Gradt = &Gradt;
    hxcurldata.A_grad = &A_grad;
    hxcurldata.amgparam_grad = amgparam;
    hxcurldata.mgl_grad = mgl_grad;

    hxcurldata.backup_r = (REAL*)calloc(A->row, sizeof(REAL));
    hxcurldata.w = (REAL*)calloc(A->row, sizeof(REAL));

    precond pc; pc.data = &hxcurldata;
    switch (itparam->linear_precond_type) {

        case PREC_HX_CURL_A: //additive HX preconditioner
            pc.fct = precond_hx_curl_additive;
            break;

        default:  // multiplicative HX preconditioner
            pc.fct = precond_hx_curl_multiplicative;
            break;

    }

    // call iterative solver
    status = solver_dcsr_linear_itsolver(A, b, x, &pc, itparam);

    if ( prtlvl >= PRINT_MIN ) {
        get_time(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("HX_curl_Krylov method totally", solver_duration);
        printf("**********************************************************\n");
    }

FINISHED:
    HX_curl_data_free(&hxcurldata, FALSE);

    return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_dcsr_krylov_hx_div (dCSRmat *A, dvector *b, dvector *x,
 *                                      linear_itsolver_param *itparam, AMG_param *amgparam,
 *                                      dCSRmat P_curl, dCSRmat Curl)
 *
 * \brief Solve Ax=b by HX (H(div)) preconditioned Krylov methods
 *
 * \param A         Pointer to the coeff matrix in dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG methods
 * \param P_curl    Pointer to the Pi_curl interpolation in dCSRmat format
 * \param P_div     Pointer to the Pi_div interpolation in dCSRmat format
 * \param Curl      Pointer to the Curl operator in dCSRmat format
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \note  P_curl is only needed for 3D H(div) problem (==NULL for 2D problem)
 * \note  Curl is rotated gradient in 2D
 *
 * \author Xiaozhe Hu
 * \date   02/10/2016
 */
INT linear_solver_dcsr_krylov_hx_div(dCSRmat *A,
                                     dvector *b,
                                     dvector *x,
                                     linear_itsolver_param *itparam,
                                     AMG_param *amgparam,
                                     dCSRmat *P_curl,
                                     dCSRmat *P_div,
                                     dCSRmat *Curl)
{

    const SHORT prtlvl = itparam->linear_print_level;
    //amgparam->max_levels = 2;
    const SHORT max_levels = amgparam->max_levels;

    /*------------------------*/
    /* Local Variables */
    /*------------------------*/
    INT      status = SUCCESS;
    REAL     solver_start, solver_end, solver_duration;

    get_time(&solver_start);

    /*------------------------*/
    /* setup Laplacians */
    /*------------------------*/
    // get transpose of P_curl
    dCSRmat Pt_curl;
    if (P_curl == NULL) // 2D problem
    {
      dcsr_null(&Pt_curl);
    }
    else // 3D problem
    {
      dcsr_trans(P_curl, &Pt_curl);
    }
    // get transpose of P_div
    dCSRmat Pt_div;
    dcsr_trans(P_div, &Pt_div);
    // get transpose of Curl
    dCSRmat Curlt;
    dcsr_trans(Curl, &Curlt);

    // get A_curl (in 3D)
    dCSRmat A_curl;
    if (P_curl == NULL) // 2D
    {
      dcsr_null(&A_curl);
    }
    else // 3D
    {
      dcsr_rap(&Curlt, A, Curl, &A_curl);
    }
    // get A_grad (in 2D)
    dCSRmat A_grad;
    if (P_curl == NULL) //2D
    {
      dcsr_rap(&Curlt, A, Curl, &A_grad);
    }
    else // 3D
    {
      dcsr_null(&A_grad);
    }
    // get A_curlgrad
    dCSRmat A_curlgrad;
    if (P_curl == NULL) // 2D problem
    {
      dcsr_null(&A_curlgrad);
    }
    else // 3D problem
    {
      dcsr_rap(&Pt_curl, &A_curl, P_curl, &A_curlgrad);
    }
    // get A_divgrad
    dCSRmat A_divgrad;
    dcsr_rap(&Pt_div, A, P_div, &A_divgrad);

    /*------------------------*/
    /* setup AMG */
    /*------------------------*/
    AMG_data *mgl_grad = amg_data_create(max_levels);
    AMG_data *mgl_curlgrad = amg_data_create(max_levels);

    if (P_curl == NULL) // 2D
    {
      // initialize A, b, x for mgl_grad[0]
      mgl_grad[0].A=dcsr_create(A_grad.row,A_grad.col,A_grad.nnz);
      dcsr_cp(&A_grad, &mgl_grad[0].A);
      mgl_grad[0].b=dvec_create(A_grad.col);
      mgl_grad[0].x=dvec_create(A_grad.row);

      // setup AMG for Laplacian
      switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
        status = amg_setup_ua(mgl_grad, amgparam); break;

        case SA_AMG: // Smoothed Aggregation AMG
        status = amg_setup_sa(mgl_grad, amgparam); break;

        default: // Unsmoothed Aggregation AMG
        status = amg_setup_ua(mgl_grad, amgparam); break;

      }

      if (status < 0) goto FINISHED;

    }
    else // 3D
    {
      // initialize A, b, x for mgl_curlgrad[0]
      mgl_curlgrad[0].A=dcsr_create(A_curlgrad.row,A_curlgrad.col,A_curlgrad.nnz);
      dcsr_cp(&A_curlgrad, &mgl_curlgrad[0].A);
      mgl_curlgrad[0].b=dvec_create(A_curlgrad.col);
      mgl_curlgrad[0].x=dvec_create(A_curlgrad.row);

      // setup AMG for vector Laplacian
      switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
        status = amg_setup_ua(mgl_curlgrad, amgparam); break;

        case SA_AMG: // Smoothed Aggregation AMG
        status = amg_setup_sa(mgl_curlgrad, amgparam); break;

        default: // Unsmoothed Aggregation AMG
        status = amg_setup_ua(mgl_curlgrad, amgparam); break;

      }

      if (status < 0) goto FINISHED;

    }

    // initialize A, b, x for mgl_divgrad[0]
    AMG_data *mgl_divgrad = amg_data_create(max_levels);
    mgl_divgrad[0].A=dcsr_create(A_divgrad.row,A_divgrad.col,A_divgrad.nnz);
    dcsr_cp(&A_divgrad, &mgl_divgrad[0].A);
    mgl_divgrad[0].b=dvec_create(A_divgrad.col);
    mgl_divgrad[0].x=dvec_create(A_divgrad.row);

    // setup AMG for vector Laplacian
    switch (amgparam->AMG_type) {

      case UA_AMG: // Unsmoothed Aggregation AMG
      status = amg_setup_ua(mgl_divgrad, amgparam); break;

      case SA_AMG: // Smoothed Aggregation AMG
      status = amg_setup_sa(mgl_divgrad, amgparam); break;

      default: // Unsmoothed Aggregation AMG
      status = amg_setup_ua(mgl_divgrad, amgparam); break;

    }

    if (status < 0) goto FINISHED;

    /*------------------------*/
    // setup preconditioner
    HX_div_data hxdivdata;

    hxdivdata.A = A;

    hxdivdata.smooth_type = 1;
    hxdivdata.smooth_iter = itparam->HX_smooth_iter;

    hxdivdata.P_curl = P_curl;
    hxdivdata.Pt_curl = &Pt_curl;
    hxdivdata.P_div = P_div;
    hxdivdata.Pt_div = &Pt_div;
    hxdivdata.Curl = Curl;
    hxdivdata.Curlt = &Curlt;

    hxdivdata.A_curlgrad = &A_curlgrad;
    hxdivdata.A_divgrad = &A_divgrad;
    hxdivdata.A_curl = &A_curl;
    hxdivdata.A_grad = &A_grad;
    hxdivdata.amgparam_curlgrad = amgparam;
    hxdivdata.mgl_curlgrad = mgl_curlgrad;
    hxdivdata.amgparam_divgrad = amgparam;
    hxdivdata.mgl_divgrad = mgl_divgrad;
    hxdivdata.amgparam_grad = amgparam;
    hxdivdata.mgl_grad = mgl_divgrad;

    hxdivdata.backup_r = (REAL*)calloc(A->row, sizeof(REAL));
    hxdivdata.w = (REAL*)calloc(2*(A_curl.row)+A->row, sizeof(REAL));

    precond pc; pc.data = &hxdivdata;
    switch (itparam->linear_precond_type) {

        case PREC_HX_DIV_A: //additive HX preconditioner
            if (P_curl == NULL) // 2D
            {
              pc.fct = precond_hx_div_additive_2D;
            }
            else // 3D
            {
              pc.fct = precond_hx_div_additive;
            }
            break;

        default:  // multiplicative HX preconditioner
            if (P_curl == NULL) // 2D
            {
              pc.fct = precond_hx_div_multiplicative_2D;
            }
            else // 3D
            {
              pc.fct = precond_hx_div_multiplicative;
            }
            break;

    }

    // call iterative solver
    status = solver_dcsr_linear_itsolver(A, b, x, &pc, itparam);

    if ( prtlvl >= PRINT_MIN ) {
        get_time(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("HX_curl_Krylov method totally", solver_duration);
        printf("**********************************************************\n");
    }

FINISHED:
    HX_div_data_free(&hxdivdata, FALSE);

    return status;
}

/********************************************************************************************/
// preconditioned Krylov methods for block CSR format
/********************************************************************************************/
/**
 * \fn INT linear_solver_bdcsr_krylov (block_dCSRmat *A, dvector *b, dvector *x,
 *                                   linear_itsolver_param *itparam)
 *
 * \brief Solve Ax = b by standard Krylov methods
 *
 * \param A         Pointer to the coeff matrix in block_dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   07/18/2010
 */
INT linear_solver_bdcsr_krylov(block_dCSRmat *A,
                               dvector *b,
                               dvector *x,
                               linear_itsolver_param *itparam)
{
    const SHORT prtlvl = itparam->linear_print_level;

    INT status = SUCCESS;
    REAL solver_start, solver_end, solver_duration;

    // solver part
    get_time(&solver_start);

    status = solver_bdcsr_linear_itsolver(A,b,x,NULL,itparam);

    get_time(&solver_end);

    solver_duration = solver_end - solver_start;

    if ( prtlvl >= PRINT_MIN ){
        print_cputime("Krylov method totally", solver_duration);
        printf("**********************************************************\n");
    }


    return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_bdcsr_krylov_block_2 (block_dCSRmat *A, dvector *b, dvector *x,
 *                                           itsolver_param *itparam,
 *                                           AMG_param *amgparam, dCSRmat *A_diag)
 *
 * \brief Solve Ax = b by preconditioned Krylov methods
 *
 * \note  Use block preconditioners where the diagonal blocks of the preconditioner is stored in A_diag
 * \note  Each diagonal block is solved exactly or approximately by AMG/AMG+krylov methods
 *
 * \param A         Pointer to the coeff matrix in block_dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG solvers
 * \param A_diag    Digonal blocks of the block preconditioner
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   10/14/2016
 *
 * \note only works for 2 by 2 block dCSRmat problems!! -- Xiaozhe Hu
 */
INT linear_solver_bdcsr_krylov_block_2(block_dCSRmat *A,
                                       dvector *b,
                                       dvector *x,
                                       linear_itsolver_param *itparam,
                                       AMG_param *amgparam,
                                       dCSRmat *A_diag)
{
  const SHORT prtlvl = itparam->linear_print_level;
  const SHORT precond_type = itparam->linear_precond_type;

  INT i;
  INT status = SUCCESS;
  REAL setup_start, setup_end, setup_duration;
  REAL solver_start, solver_end, solver_duration;

#if WITH_SUITESPARSE
    void **LU_diag = (void **)calloc(2, sizeof(void *));
#else
    error_extlib(257, __FUNCTION__, "SuiteSparse");
#endif

  SHORT max_levels;
  if (amgparam) max_levels = amgparam->max_levels;
  AMG_data **mgl = (AMG_data **)calloc(2, sizeof(AMG_data *));

  /* setup preconditioner */
  get_time(&setup_start);

  if (precond_type > 0 && precond_type < 20) {
  /* diagonal blocks are solved exactly */
#if WITH_SUITESPARSE
    // Need to sort the diagonal blocks for UMFPACK format
    dCSRmat A_tran;

    for (i=0; i<2; i++){

        dcsr_trans(&A_diag[i], &A_tran);
        dcsr_cp(&A_tran, &A_diag[i]);

        if ( prtlvl > PRINT_NONE ) printf("Factorization for %d-th diagonal block: \n", i);
        LU_diag[i] = umfpack_factorize(&A_diag[i], prtlvl);

        dcsr_free(&A_tran);

    }
#endif
  }
  else {

      for (i=0; i<2; i++){
          /* set AMG for diagonal blocks */
          mgl[i] = amg_data_create(max_levels);
          dcsr_alloc(A_diag[i].row, A_diag[i].row, A_diag[i].nnz, &mgl[i][0].A);
          dcsr_cp(&(A_diag[i]), &mgl[i][0].A);
          mgl[i][0].b=dvec_create(A_diag[i].row);
          mgl[i][0].x=dvec_create(A_diag[i].row);

          switch (amgparam->AMG_type) {

              case UA_AMG: // Unsmoothed Aggregation AMG
                  if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                  status = amg_setup_ua(mgl[i], amgparam);
                  break;

              case SA_AMG: // Smoothed Aggregation AMG
                  if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
                  status = amg_setup_sa(mgl[i], amgparam);
                  break;

              default: // UA AMG
                  if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                  status = amg_setup_ua(mgl[i], amgparam);
                  break;

          }

      }

  }


  precond_block_data precdata;
  precond_block_data_null(&precdata);

  precdata.Abcsr = A;

  precdata.A_diag = A_diag;
  precdata.r = dvec_create(b->row);
  if (amgparam) precdata.amgparam = amgparam;

  if (precond_type > 0 && precond_type < 20) {
  /* diagonal blocks are solved exactly */
#if WITH_SUITESPARSE
      precdata.LU_diag = LU_diag;
#endif
  }
  else {
      precdata.mgl = mgl;
  }

  precdata.hxcurldata = NULL;
  precdata.hxdivdata = NULL;

  precond prec; prec.data = &precdata;

  switch (precond_type)
  {
    case 10:
      prec.fct = precond_block_diag_2;
      break;

    case 11:
      prec.fct = precond_block_lower_2;
      break;

    case 12:
      prec.fct = precond_block_upper_2;
      break;

    case 20:
      prec.fct = precond_block_diag_2_amg;
      break;

    case 21:
      prec.fct = precond_block_lower_2_amg;
      break;

    case 22:
      prec.fct = precond_block_upper_2_amg;
      break;

    case 30:
      prec.fct = precond_block_diag_2_amg_krylov;
      break;

    case 31:
      prec.fct = precond_block_lower_2_amg_krylov;
      break;

    case 32:
      prec.fct = precond_block_upper_2_amg_krylov;
      break;

    default:
      break;
  }


  if ( prtlvl >= PRINT_MIN ) {
    get_time(&setup_end);
    setup_duration = setup_end - setup_start;
    print_cputime("Setup totally", setup_duration);
  }

  // solver part
  get_time(&solver_start);

  status=solver_bdcsr_linear_itsolver(A,b,x, &prec,itparam);

  get_time(&solver_end);

  solver_duration = solver_end - solver_start;

  if ( prtlvl >= PRINT_MIN ) {
    print_cputime("Krylov method totally", solver_duration);
    printf("**********************************************************\n");
  }

  // clean
  precond_block_data_free(&precdata, 2, TRUE);

  return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_bdcsr_krylov_block_3 (block_dCSRmat *A, dvector *b, dvector *x,
 *                                           itsolver_param *itparam,
 *                                           AMG_param *amgparam, dCSRmat *A_diag)
 *
 * \brief Solve Ax = b by preconditioned Krylov methods
 *
 * \note  Use block preconditioners where the diagonal blocks of the preconditioner is stored in A_diag
 * \note  Each diagonal block is solved exactly or approximately by AMG/AMG+krylov methods
 *
 * \param A         Pointer to the coeff matrix in block_dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG solvers
 * \param A_diag    Digonal blocks of A
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   02/24/2014
 *
 * \note only works for 3 by 3 block dCSRmat problems!! -- Xiaozhe Hu
 */
INT linear_solver_bdcsr_krylov_block_3(block_dCSRmat *A,
                                       dvector *b,
                                       dvector *x,
                                       linear_itsolver_param *itparam,
                                       AMG_param *amgparam,
                                       dCSRmat *A_diag)
{

    const SHORT prtlvl = itparam->linear_print_level;
    const SHORT precond_type = itparam->linear_precond_type;

//#if WITH_SUITESPARSE
    INT i;
//#endif
    INT status = SUCCESS;
    REAL setup_start, setup_end, setup_duration;
    REAL solver_start, solver_end, solver_duration;

#if WITH_SUITESPARSE
    void **LU_diag = (void **)calloc(3, sizeof(void *));
//#else
//    error_extlib(257, __FUNCTION__, "SuiteSparse");
#endif


    SHORT max_levels;
    if (amgparam) max_levels = amgparam->max_levels;
    AMG_data **mgl = (AMG_data **)calloc(3, sizeof(AMG_data *));

    /* setup preconditioner */
    get_time(&setup_start);

    if (precond_type > 0 && precond_type < 20) {
    /* diagonal blocks are solved exactly */
#if WITH_SUITESPARSE
        // Need to sort the diagonal blocks for UMFPACK format
        dCSRmat A_tran;

        for (i=0; i<3; i++){

            dcsr_trans(&A_diag[i], &A_tran);
            dcsr_cp(&A_tran, &A_diag[i]);

            if ( prtlvl > PRINT_NONE ) printf("Factorization for %d-th diagonal block:\n", i);
            LU_diag[i] = umfpack_factorize(&A_diag[i], prtlvl);

            dcsr_free(&A_tran);

        }

#endif
    }
    else {

        for (i=0; i<3; i++){

            /* set AMG for diagonal blocks */
            mgl[i] = amg_data_create(max_levels);
            dcsr_alloc(A_diag[i].row, A_diag[i].row, A_diag[i].nnz, &mgl[i][0].A);
            dcsr_cp(&(A_diag[i]), &mgl[i][0].A);
            mgl[i][0].b=dvec_create(A_diag[i].row);
            mgl[i][0].x=dvec_create(A_diag[i].row);

            switch (amgparam->AMG_type) {

                case UA_AMG: // Unsmoothed Aggregation AMG
                    if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                    status = amg_setup_ua(mgl[i], amgparam);
                    break;

                case SA_AMG: // Smoothed Aggregation AMG
                    if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
                    status = amg_setup_sa(mgl[i], amgparam);
                    break;

                default: // UA AMG
                    if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                    status = amg_setup_ua(mgl[i], amgparam);
                    break;

            }

        }

    }

    precond_block_data precdata;
    precond_block_data_null(&precdata);

    precdata.Abcsr = A;

    precdata.A_diag = A_diag;
    precdata.r = dvec_create(b->row);
    if (amgparam) precdata.amgparam = amgparam;


    if (precond_type > 0 && precond_type < 20) {
    /* diagonal blocks are solved exactly */
#if WITH_SUITESPARSE
        precdata.LU_diag = LU_diag;
#endif
    }
    else {
        precdata.mgl = mgl;
    }

    precond prec; prec.data = &precdata;

    switch (precond_type)
    {
        case 10:
            prec.fct = precond_block_diag_3;
            break;

        case 11:
            prec.fct = precond_block_lower_3;
            break;

        case 12:
            prec.fct = precond_block_upper_3;
            break;

        case 20:
            prec.fct = precond_block_diag_3_amg;
            break;

        case 21:
            prec.fct = precond_block_lower_3_amg;
            break;

        case 22:
            prec.fct = precond_block_upper_3_amg;
            break;

        case 30:
            prec.fct = precond_block_diag_3_amg_krylov;
            break;

        case 31:
            prec.fct = precond_block_lower_3_amg_krylov;
            break;

        case 32:
            prec.fct = precond_block_upper_3_amg_krylov;
            break;

        default:
            break;
    }


    if ( prtlvl >= PRINT_MIN ) {
        get_time(&setup_end);
        setup_duration = setup_end - setup_start;
        print_cputime("Setup totally", setup_duration);
    }

    // solver part
    get_time(&solver_start);

    status=solver_bdcsr_linear_itsolver(A,b,x, &prec,itparam);

    get_time(&solver_end);

    solver_duration = solver_end - solver_start;

    if ( prtlvl >= PRINT_MIN ) {
        print_cputime("Krylov method totally", solver_duration);
        printf("**********************************************************\n");
    }

    // clean
    precond_block_data_free(&precdata, 3, TRUE);

    return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_bdcsr_krylov_block_4 (block_dCSRmat *A, dvector *b, dvector *x,
 *                                           itsolver_param *itparam,
 *                                           AMG_param *amgparam, dCSRmat *A_diag)
 *
 * \brief Solve Ax = b by preconditioned Krylov methods
 *
 * \note  Use block preconditioners where the diagonal blocks of the preconditioner is stored in A_diag
 * \note  Each diagonal block is solved exactly or approximately by AMG/AMG+krylov methods
 *
 * \param A         Pointer to the coeff matrix in block_dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG solvers
 * \param A_diag    Digonal blocks of A
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   01/20/2017
 *
 * \note only works for 4 by 4 block dCSRmat problems!! -- Xiaozhe Hu
 */
INT linear_solver_bdcsr_krylov_block_4(block_dCSRmat *A,
                                       dvector *b,
                                       dvector *x,
                                       linear_itsolver_param *itparam,
                                       AMG_param *amgparam,
                                       dCSRmat *A_diag)
{
    const SHORT prtlvl = itparam->linear_print_level;
    const SHORT precond_type = itparam->linear_precond_type;

    INT i;
    INT status = SUCCESS;
    REAL setup_start, setup_end, setup_duration;
    REAL solver_start, solver_end, solver_duration;

#if WITH_SUITESPARSE
    void **LU_diag = (void **)calloc(4, sizeof(void *));
#else
    error_extlib(258, __FUNCTION__, "SuiteSparse");
#endif

    SHORT max_levels;
    if (amgparam) max_levels = amgparam->max_levels;
    AMG_data **mgl = (AMG_data **)calloc(4, sizeof(AMG_data *));

    /* setup preconditioner */
    get_time(&setup_start);

    if (precond_type > 0 && precond_type < 20) {
    /* diagonal blocks are solved exactly */
#if WITH_SUITESPARSE
        // Need to sort the diagonal blocks for UMFPACK format
        dCSRmat A_tran;

        for (i=0; i<4; i++){

            dcsr_trans(&A_diag[i], &A_tran);
            dcsr_cp(&A_tran, &A_diag[i]);

            if ( prtlvl > PRINT_NONE ) printf("Factorization for %d-th diagonal block:\n", i);
            LU_diag[i] = umfpack_factorize(&A_diag[i], prtlvl);

            dcsr_free(&A_tran);


        }
#endif
    }
    else {

      for (i=0; i<4; i++){

          /* set AMG for diagonal blocks */
          mgl[i] = amg_data_create(max_levels);
          dcsr_alloc(A_diag[i].row, A_diag[i].row, A_diag[i].nnz, &mgl[i][0].A);
          dcsr_cp(&(A_diag[i]), &mgl[i][0].A);
          mgl[i][0].b=dvec_create(A_diag[i].row);
          mgl[i][0].x=dvec_create(A_diag[i].row);

          switch (amgparam->AMG_type) {

              case UA_AMG: // Unsmoothed Aggregation AMG
                  if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                  status = amg_setup_ua(mgl[i], amgparam);
                  break;

              case SA_AMG: // Smoothed Aggregation AMG
                  if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
                  status = amg_setup_sa(mgl[i], amgparam);
                  break;

              default: // UA AMG
                  if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                  status = amg_setup_ua(mgl[i], amgparam);
                  break;

          }

          if (status < 0) {
            printf("### HAZMATH WARNING: AMG setup failed!\n");
            goto FINISHED;
          }

      }

    }

    precond_block_data precdata;
    precond_block_data_null(&precdata);

    precdata.Abcsr = A;
    precdata.A_diag = A_diag;
    precdata.r = dvec_create(b->row);
    if (amgparam) precdata.amgparam = amgparam;

    if (precond_type > 0 && precond_type < 20) {
    /* diagonal blocks are solved exactly */
#if WITH_SUITESPARSE
        precdata.LU_diag = LU_diag;
#endif
    }
    else {
      precdata.mgl = mgl;
    }

    precond prec; prec.data = &precdata;

    switch (precond_type)
    {
        case 10:
            prec.fct = precond_block_diag_4;
            break;

        case 11:
            prec.fct = precond_block_lower_4;
            break;

        case 12:
            prec.fct = precond_block_upper_4;
            break;

        case 20:
            prec.fct = precond_block_diag_4_amg;
            break;

        case 21:
            prec.fct = precond_block_lower_4_amg;
            break;

        case 22:
            prec.fct = precond_block_upper_4_amg;
            break;

        case 30:
            prec.fct = precond_block_diag_4_amg_krylov;
            break;

        case 31:
            prec.fct = precond_block_lower_4_amg_krylov;
            break;

        case 32:
            prec.fct = precond_block_upper_4_amg_krylov;
            break;

        default:
            break;
    }

    if ( prtlvl >= PRINT_MIN ) {
        get_time(&setup_end);
        setup_duration = setup_end - setup_start;
        print_cputime("Setup totally", setup_duration);
    }

    // solver part
    get_time(&solver_start);

    status=solver_bdcsr_linear_itsolver(A,b,x, &prec,itparam);

    get_time(&solver_end);

    solver_duration = solver_end - solver_start;

    if ( prtlvl >= PRINT_MIN ) {
        print_cputime("Krylov method totally", solver_duration);
        printf("**********************************************************\n");
    }

FINISHED:
    // clean
    precond_block_data_free(&precdata, 4, TRUE);

    return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_bdcsr_krylov_block_5 (block_dCSRmat *A, dvector *b, dvector *x,
 *                                           itsolver_param *itparam,
 *                                           AMG_param *amgparam, dCSRmat *A_diag)
 *
 * \brief Solve Ax = b by preconditioned Krylov methods
 *
 * \note  Use block preconditioners where the diagonal blocks of the preconditioner is stored in A_diag
 * \note  Each diagonal block is solved exactly or approximately by AMG/AMG+krylov methods
 *
 * \param A         Pointer to the coeff matrix in block_dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG solvers
 * \param A_diag    Digonal blocks of A
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   06/19/2020
 *
 * \note only works for 5 by 5 block dCSRmat problems!! -- Xiaozhe Hu
 * \note only block diagonal preconditioners have been implemented!! -- Xiaozhe Hu
 */
INT linear_solver_bdcsr_krylov_block_5(block_dCSRmat *A,
                                       dvector *b,
                                       dvector *x,
                                       linear_itsolver_param *itparam,
                                       AMG_param *amgparam,
                                       dCSRmat *A_diag)
{
    const SHORT prtlvl = itparam->linear_print_level;
    const SHORT precond_type = itparam->linear_precond_type;

    INT i;
    INT status = SUCCESS;
    REAL setup_start, setup_end, setup_duration;
    REAL solver_start, solver_end, solver_duration;

#if WITH_SUITESPARSE
    void **LU_diag = (void **)calloc(5, sizeof(void *));
//#else
//    error_extlib(257, __FUNCTION__, "SuiteSparse");
#endif

    SHORT max_levels;
    if (amgparam) max_levels = amgparam->max_levels;
    AMG_data **mgl = (AMG_data **)calloc(5, sizeof(AMG_data *));

    /* setup preconditioner */
    get_time(&setup_start);

    if (precond_type > 0 && precond_type < 20) {
    /* diagonal blocks are solved exactly */
#if WITH_SUITESPARSE
        // Need to sort the diagonal blocks for UMFPACK format
        dCSRmat A_tran;

        for (i=0; i<5; i++){

            dcsr_trans(&A_diag[i], &A_tran);
            dcsr_cp(&A_tran, &A_diag[i]);

            if ( prtlvl > PRINT_NONE ) printf("Factorization for %d-th diagonal block:\n", i);
            LU_diag[i] = umfpack_factorize(&A_diag[i], prtlvl);

            dcsr_free(&A_tran);

        }

#endif
    }
    else {

        for (i=0; i<5; i++){

            /* set AMG for diagonal blocks */
            mgl[i] = amg_data_create(max_levels);
            dcsr_alloc(A_diag[i].row, A_diag[i].row, A_diag[i].nnz, &mgl[i][0].A);
            dcsr_cp(&(A_diag[i]), &mgl[i][0].A);
            mgl[i][0].b=dvec_create(A_diag[i].row);
            mgl[i][0].x=dvec_create(A_diag[i].row);

            switch (amgparam->AMG_type) {

                case UA_AMG: // Unsmoothed Aggregation AMG
                    if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                    status = amg_setup_ua(mgl[i], amgparam);
                    break;

                case SA_AMG: // Smoothed Aggregation AMG
                    if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
                    status = amg_setup_sa(mgl[i], amgparam);
                    break;

                default: // UA AMG
                    if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                    status = amg_setup_ua(mgl[i], amgparam);
                    break;

            }

        }

    }

    precond_block_data precdata;
    precond_block_data_null(&precdata);

    precdata.Abcsr = A;

    precdata.A_diag = A_diag;
    precdata.r = dvec_create(b->row);
    if (amgparam) precdata.amgparam = amgparam;


    if (precond_type > 0 && precond_type < 20) {
    /* diagonal blocks are solved exactly */
#if WITH_SUITESPARSE
        precdata.LU_diag = LU_diag;
#endif
    }
    else {
        precdata.mgl = mgl;
    }

    precond prec; prec.data = &precdata;

    switch (precond_type)
    {

        case 10:
            prec.fct = precond_block_diag;
            break;
        /*
        case 11:
            prec.fct = precond_block_lower_5;
            break;

        case 12:
            prec.fct = precond_block_upper_5;
            break;

        case 20:
            prec.fct = precond_block_diag_5_amg;
            break;

        case 21:
            prec.fct = precond_block_lower_5_amg;
            break;

        case 22:
            prec.fct = precond_block_upper_5_amg;
            break;
       */

        case 30:
            prec.fct = precond_block_diag_5_amg_krylov;
            break;

        //case 31:
        //    prec.fct = precond_block_lower_5_amg_krylov;
        //    break;

        //case 32:
        //    prec.fct = precond_block_upper_5_amg_krylov;
        //    break;

        default:
            break;
    }


    if ( prtlvl >= PRINT_MIN ) {
        get_time(&setup_end);
        setup_duration = setup_end - setup_start;
        print_cputime("Setup totally", setup_duration);
    }

    // solver part
    get_time(&solver_start);

    status=solver_bdcsr_linear_itsolver(A,b,x, &prec,itparam);

    get_time(&solver_end);

    solver_duration = solver_end - solver_start;

    if ( prtlvl >= PRINT_MIN ) {
        print_cputime("Krylov method totally", solver_duration);
        printf("**********************************************************\n");
    }

    // clean
    precond_block_data_free(&precdata, 3, TRUE);

    return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_bdcsr_krylov_block (block_dCSRmat *A, dvector *b, dvector *x,
 *                                           itsolver_param *itparam,
 *                                           AMG_param *amgparam, dCSRmat *A_diag)
 *
 * \brief Solve Ax = b by preconditioned Krylov methods
 *
 * \note  Use block preconditioners where the diagonal blocks of the preconditioner is stored in A_diag
 * \note  Each diagonal block is solved exactly
 *
 * \param A         Pointer to the coeff matrix in block_dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG solvers
 * \param A_diag    Digonal blocks of A
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   04/05/2017
 *
 * \note  works for general block dCSRmat problems!! -- Xiaozhe Hu
 * \note  only block diagonal preconditioner has been implemented!! -- Xiaozhe Hu
 *
 */
INT linear_solver_bdcsr_krylov_block(block_dCSRmat *A,
                                       dvector *b,
                                       dvector *x,
                                       linear_itsolver_param *itparam,
                                       AMG_param *amgparam,
                                       dCSRmat *A_diag)
{
  const SHORT prtlvl = itparam->linear_print_level;
  const SHORT precond_type = itparam->linear_precond_type;

  const INT nb = A->brow;

#if WITH_SUITESPARSE
  INT i;
#endif
  INT status = SUCCESS;
  REAL setup_start, setup_end, setup_duration;
  REAL solver_start, solver_end, solver_duration;

#if WITH_SUITESPARSE
  void **LU_diag = (void **)calloc(nb, sizeof(void *));
#else
    error_extlib(256, __FUNCTION__, "SuiteSparse");
#endif

  /* setup preconditioner */
  get_time(&setup_start);

  /* diagonal blocks are solved exactly */
#if WITH_SUITESPARSE
  // Need to sort the diagonal blocks for UMFPACK format
  dCSRmat A_tran;

  for (i=0; i<nb; i++){

    dcsr_trans(&A_diag[i], &A_tran);
    dcsr_cp(&A_tran, &A_diag[i]);

    if ( prtlvl > PRINT_NONE ) printf("Factorization for %d-th diagonal block:\n", i);
    LU_diag[i] = umfpack_factorize(&A_diag[i], prtlvl);

    dcsr_free(&A_tran);

  }
#endif

  precond_block_data precdata;
  precond_block_data_null(&precdata);

  precdata.Abcsr = A;

  precdata.A_diag = A_diag;
  precdata.r = dvec_create(b->row);

#if WITH_SUITESPARSE
  precdata.LU_diag = LU_diag;
#endif

  precond prec; prec.data = &precdata;


  switch (precond_type)
  {
    case 10:
      prec.fct = precond_block_diag;
      break;

    /*
    case 11:
      prec.fct = precond_block_lower;
      break;

    case 12:
      prec.fct = precond_block_upper;
      break;
   */

    default:
      prec.fct = precond_block_diag;
      break;
  }


  if ( prtlvl >= PRINT_MIN ) {
    get_time(&setup_end);
    setup_duration = setup_end - setup_start;
    print_cputime("Setup totally", setup_duration);
  }

  // solver part
  get_time(&solver_start);

  status=solver_bdcsr_linear_itsolver(A,b,x, &prec,itparam);

  get_time(&solver_end);

  solver_duration = solver_end - solver_start;

  if ( prtlvl >= PRINT_MIN ) {
    print_cputime("Krylov method totally", solver_duration);
    printf("**********************************************************\n");
  }

  // clean
  precond_block_data_free(&precdata, nb, TRUE);

  return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_bdcsr_krylov_mixed_darcy (block_dCSRmat *A, dvector *b, dvector *x,
 *                                           itsolver_param *itparam,
 *                                           AMG_param *amgparam, dCSRmat *A_diag)
 *
 * \brief Solve Ax = b by standard Krylov methods
 *
 * \param A         Pointer to the coeff matrix in block_dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG solvers
 * \param A_diag    Digonal blocks of A
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   10/14/2016
 *
 * \note only works for 2 by 2 block dCSRmat problems!! -- Xiaozhe Hu
 * \note this function needs to be moved!! -- Xiaozhe Hu
 */
INT linear_solver_bdcsr_krylov_mixed_darcy(block_dCSRmat *A,
                                           dvector *b,
                                           dvector *x,
                                           linear_itsolver_param *itparam,
                                           AMG_param *amgparam,
                                           dCSRmat *P_div,
                                           dCSRmat *Curl,
                                           dCSRmat *P_curl,
                                           dvector *el_vol)
{
  // variables
  const SHORT prtlvl = itparam->linear_print_level;
  const SHORT precond_type = itparam->linear_precond_type;

  INT status = SUCCESS;
  REAL setup_start, setup_end, setup_duration;
  REAL solver_start, solver_end, solver_duration;

  const SHORT max_levels = amgparam->max_levels;
  INT i, n, np;

  // initilize data for preconditioners of each diagnal block
  AMG_data **mgl = (AMG_data **)calloc(2, sizeof(AMG_data *));
  for (i=0; i<2; i++) mgl[i]=NULL;
  dvector **diag = (dvector **)calloc(2, sizeof(dvector *));
  for (i=0; i<2; i++) diag[i]=NULL;
  HX_div_data **hxdivdata = (HX_div_data **)calloc(2, sizeof(HX_div_data *));
  for (i=0; i<2; i++) hxdivdata[i] = NULL;

  // data for argumented Lagrange type preconditioner
  dCSRmat BTB;

  // data for HX preconditioner
  dCSRmat A_div, A_curl, A_grad, A_divgrad, A_curlgrad;
  dCSRmat Pt_div, Curlt, Pt_curl;
  AMG_data *mgl_grad = amg_data_create(max_levels);
  AMG_data *mgl_curlgrad = amg_data_create(max_levels);
  AMG_data *mgl_divgrad = amg_data_create(max_levels);

  /* setup preconditioner */
  get_time(&setup_start);

  /*----------------------*/
  /* Use argumented Lagrange type preconditioner */
  /*----------------------*/
  if (precond_type < 30)
  {
      /* set AMG for the flux block */
      mgl[0] = amg_data_create(max_levels);
      n = A->blocks[0]->row;
      np = A->blocks[3]->row;

      dCSRmat invMp = dcsr_create(np,np,np);
      for (i=0;i<np;i++)
      {
          invMp.IA[i] = i;
          invMp.JA[i] = i;
          if (el_vol->val[i] > SMALLREAL) invMp.val[i]   = 1.0/(el_vol->val[i]);
          else invMp.val[i] = 1.0;
      }
      invMp.IA[np] = np;

      //dcsr_mxm(A->blocks[1], A->blocks[2], &BTB);
      dcsr_rap(A->blocks[1], &invMp, A->blocks[2], &BTB);
      dcsr_add(&BTB, 1.0, A->blocks[0], 1.0, &mgl[0][0].A);

      mgl[0][0].b=dvec_create(n); mgl[0][0].x=dvec_create(n);

      switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
          status = amg_setup_ua(mgl[0], amgparam);
          break;

       case SA_AMG: // Smoothed Aggregation AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
          status = amg_setup_sa(mgl[0], amgparam);
          break;

        default: // Unsmoothed Aggregation AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
          status = amg_setup_ua(mgl[0], amgparam);
          break;

      }

      dcsr_free(&BTB);
      //dcsr_free(&invMp);

  }
  else if (precond_type > 29 && precond_type < 50)
  {
    /*  HX preconditioner for the flux block */
    hxdivdata[0] = (HX_div_data *)calloc(1, sizeof(HX_div_data));

    //-----------------------
    // form divdiv block
    //-----------------------
    n = A->blocks[0]->row;
    np = A->blocks[3]->row;

    dCSRmat invMp = dcsr_create(np,np,np);
    for (i=0;i<np;i++)
    {
        invMp.IA[i] = i;
        invMp.JA[i] = i;
        if (el_vol->val[i] > SMALLREAL) invMp.val[i]   = 1.0/(el_vol->val[i]);
        else invMp.val[i] = 1.0;
    }
    invMp.IA[np] = np;

    //dcsr_mxm(A->blocks[1], A->blocks[2], &BTB);
    dcsr_rap(A->blocks[1], &invMp, A->blocks[2], &BTB);
    //dcsr_add(&BTB, 1, A->blocks[0], 1.0, &A_div);
    dcsr_add(&BTB, itparam->AL_scaling_param, A->blocks[0], 1.0, &A_div);

    // free
    dcsr_free(&BTB);
    dcsr_free(&invMp);

    //-----------------------
    // Setup HX div preconditioner
    //-----------------------
    // get transpose of P_div
    //dCSRmat Pt_div;
    dcsr_trans(P_div, &Pt_div);
    // get transpose of Curl
    //dCSRmat Curlt;
    dcsr_trans(Curl, &Curlt);

    // 2D case: P_curl = NULL
    if (P_curl == NULL)
    {
      // get A_grad (scalar H(grad))
      dcsr_rap(&Curlt, &A_div, Curl, &A_grad);
      // get A_divgrad (vector H(grad))
      dcsr_rap(&Pt_div, &A_div, P_div, &A_divgrad);

      // initialize A, b, x for mgl_grad[0]
      //AMG_data *mgl_grad = amg_data_create(max_levels);
      mgl_grad[0].A = dcsr_create(A_grad.row, A_grad.col, A_grad.nnz);
      dcsr_cp(&A_grad, &mgl_grad[0].A);
      mgl_grad[0].b=dvec_create(A_grad.col);
      mgl_grad[0].x=dvec_create(A_grad.row);

      // setup AMG for A_grad (scalar H(grad))
      switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            status = amg_setup_ua(mgl_grad, amgparam); break;

        default: // Classical AMG
            status = amg_setup_ua(mgl_grad, amgparam); break;

      }
      if (status < 0) goto FINISHED;

      // initialize A, b, x for mgl_divgrad[0]
      //AMG_data *mgl_divgrad = amg_data_create(max_levels);
      mgl_divgrad[0].A=dcsr_create(A_divgrad.row,A_divgrad.col,A_divgrad.nnz);
      dcsr_cp(&A_divgrad, &mgl_divgrad[0].A);
      mgl_divgrad[0].b=dvec_create(A_divgrad.col);
      mgl_divgrad[0].x=dvec_create(A_divgrad.row);

      // setup AMG for vector Laplacian
      switch (amgparam->AMG_type) {

          case UA_AMG: // Unsmoothed Aggregation AMG
              status = amg_setup_ua(mgl_divgrad, amgparam); break;

          case SA_AMG: // Smoothed Aggregation AMG
                  status = amg_setup_sa(mgl_divgrad, amgparam); break;

          default: // UA AMG
              status = amg_setup_ua(mgl_divgrad, amgparam); break;

      }
      if (status < 0) goto FINISHED;

      // free amg data for curl Grad
      free(mgl_curlgrad); mgl_curlgrad = NULL;
    }
    // 3D case: P_curl exists
    else
    {
      // get transpose of P_curl
      dcsr_trans(P_curl, &Pt_curl);

      // get A_curl
      dcsr_rap(&Curlt, &A_div, Curl, &A_curl);
      // get A_curlgrad (vector H(grad) from curl)
      dcsr_rap(&Pt_curl, &A_curl, P_curl, &A_curlgrad);
      // get A_divgrad (vector H(grad) from div)
      dcsr_rap(&Pt_div, &A_div, P_div, &A_divgrad);
      //dcsr_write_dcoo("A_divgrad_0.dat", &A_divgrad);

      // initialize A, b, x for mgl_curlgrad[0]
      mgl_curlgrad[0].A=dcsr_create(A_curlgrad.row,A_curlgrad.col,A_curlgrad.nnz);
      dcsr_cp(&A_curlgrad, &mgl_curlgrad[0].A);
      mgl_curlgrad[0].b=dvec_create(A_curlgrad.col);
      mgl_curlgrad[0].x=dvec_create(A_curlgrad.row);

      // setup AMG for vector Laplacian
      switch (amgparam->AMG_type) {

          case UA_AMG: // Unsmoothed Aggregation AMG
              status = amg_setup_ua(mgl_curlgrad, amgparam); break;

          case SA_AMG: // Smoothed Aggregation AMG
              status = amg_setup_sa(mgl_curlgrad, amgparam); break;

          default: // UA AMG
              status = amg_setup_ua(mgl_curlgrad, amgparam); break;

      }
      if (status < 0) goto FINISHED;

      // initialize A, b, x for mgl_divgrad[0]
      mgl_divgrad[0].A=dcsr_create(A_divgrad.row,A_divgrad.col,A_divgrad.nnz);
      dcsr_cp(&A_divgrad, &mgl_divgrad[0].A);
      mgl_divgrad[0].b=dvec_create(A_divgrad.col);
      mgl_divgrad[0].x=dvec_create(A_divgrad.row);

      // setup AMG for vector Laplacian
      switch (amgparam->AMG_type) {

          case UA_AMG: // Unsmoothed Aggregation AMG
              status = amg_setup_ua(mgl_divgrad, amgparam); break;

          case SA_AMG: // Smoothed Aggregation AMG
              status = amg_setup_ua(mgl_divgrad, amgparam); break;

          default: // UA AMG
              status = amg_setup_ua(mgl_divgrad, amgparam); break;

      }
      if (status < 0) goto FINISHED;

      // free amg data for grad
      free(mgl_grad); mgl_grad = NULL;
    }

    //------------------------
    // assign HX preconditioner data
    //------------------------
    hxdivdata[0]->A = &A_div;
    hxdivdata[0]->smooth_type = 1;
    hxdivdata[0]->smooth_iter = itparam->HX_smooth_iter;

    // 2D case: P_curl = NULL
    if (P_curl == NULL)
    {
      hxdivdata[0]->P_curl = NULL;
      hxdivdata[0]->Pt_curl = NULL;
      hxdivdata[0]->P_div = P_div;
      hxdivdata[0]->Pt_div = &Pt_div;
      hxdivdata[0]->Curl = Curl;
      hxdivdata[0]->Curlt = &Curlt;
      hxdivdata[0]->A_curlgrad = NULL;
      hxdivdata[0]->amgparam_curlgrad = NULL;
      hxdivdata[0]->mgl_curlgrad = NULL;
      hxdivdata[0]->A_divgrad = &A_divgrad;
      hxdivdata[0]->amgparam_divgrad = amgparam;
      hxdivdata[0]->mgl_divgrad = mgl_divgrad;
      hxdivdata[0]->A_curl = NULL;
      hxdivdata[0]->A_grad = &A_grad;
      hxdivdata[0]->amgparam_grad = amgparam;
      hxdivdata[0]->mgl_grad = mgl_grad;

      hxdivdata[0]->backup_r = (REAL*)calloc(A_div.row, sizeof(REAL));
      hxdivdata[0]->w = (REAL*)calloc(A_div.row, sizeof(REAL));
    }
    // 3D case : P_curl exists
    else
    {
      hxdivdata[0]->P_curl = P_curl;
      hxdivdata[0]->Pt_curl = &Pt_curl;
      hxdivdata[0]->P_div = P_div;
      hxdivdata[0]->Pt_div = &Pt_div;
      hxdivdata[0]->Curl = Curl;
      hxdivdata[0]->Curlt = &Curlt;
      hxdivdata[0]->A_curlgrad = &A_curlgrad;
      hxdivdata[0]->amgparam_curlgrad = amgparam;
      hxdivdata[0]->mgl_curlgrad = mgl_curlgrad;
      hxdivdata[0]->A_divgrad = &A_divgrad;
      hxdivdata[0]->amgparam_divgrad = amgparam;
      hxdivdata[0]->mgl_divgrad = mgl_divgrad;
      hxdivdata[0]->A_curl = &A_curl;
      hxdivdata[0]->A_grad = NULL;
      hxdivdata[0]->amgparam_grad = amgparam;
      hxdivdata[0]->mgl_grad = NULL;

      hxdivdata[0]->backup_r = (REAL*)calloc(A_div.row, sizeof(REAL));
      hxdivdata[0]->w = (REAL*)calloc(2*A_curl.row, sizeof(REAL));
    }

  }
  /*----------------------*/
  /* Use Pressure Poission type preconditioner */
  /*----------------------*/
  else if (precond_type > 49 && precond_type < 70)
  {
      /* set AMG for the flux block */
      mgl[0] = amg_data_create(max_levels);
      n = A->blocks[0]->row;
      mgl[0][0].A = dcsr_create(n, n, A->blocks[0]->nnz);
      dcsr_cp(A->blocks[0], &mgl[0][0].A);
      mgl[0][0].b=dvec_create(n); mgl[0][0].x=dvec_create(n);

      switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
          status = amg_setup_ua(mgl[0], amgparam);
          break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = amg_setup_sa(mgl[0], amgparam);
            break;

        default: // Unsmoothed Aggregation AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
          status = amg_setup_ua(mgl[0], amgparam);
          break;

      }

      /* set AMG for the pressure block */
      mgl[1] = amg_data_create(max_levels);
      dvector diag_M;
      dCSRmat invM = dcsr_create(n,n,n);

      dcsr_getdiag(n,A->blocks[0],&diag_M);
      for (i=0;i<n;i++)
      {
          invM.IA[i] = i;
          invM.JA[i] = i;
          if (diag_M.val[i] > SMALLREAL) invM.val[i]   = 1.0/diag_M.val[i];
          else invM.val[i] = 1.0;
      }
      invM.IA[n] = n;

      //dcsr_rap(A->blocks[2], &invM, A->blocks[1], &mgl[1][0].A);
      dcsr_rap(A->blocks[2], &invM, A->blocks[1], &BTB);
      dcsr_add(&BTB, 1.0, A->blocks[3], -1.0, &mgl[1][0].A);
      mgl[1][0].b=dvec_create(A->blocks[2]->row); mgl[1][0].x=dvec_create(A->blocks[2]->row);

      switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
          status = amg_setup_ua(mgl[1], amgparam);
          break;

        case SA_AMG: // Smoothed Aggregation AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
          status = amg_setup_sa(mgl[1], amgparam);
          break;

        default: // Unsmoothed Aggregation AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
          status = amg_setup_ua(mgl[1], amgparam);
          break;

      }

      dcsr_free(&BTB);
      dcsr_free(&invM);
      dvec_free(&diag_M);

  }
  else if (precond_type >= 70)
  {
      /* set diagonal preconditioner for the flux block */
      diag[0] = (dvector *)calloc(1, sizeof(dvector));
      dcsr_getdiag(0, A->blocks[0], diag[0]);

      /* set AMG for the pressure block */
      mgl[1] = amg_data_create(max_levels);
      n = A->blocks[0]->row;

      dCSRmat invM = dcsr_create(n,n,n);

      for (i=0;i<n;i++)
      {
          invM.IA[i] = i;
          invM.JA[i] = i;
          if (diag[0]->val[i] > SMALLREAL) invM.val[i]   = 1.0/diag[0]->val[i];
          else invM.val[i] = 1.0;
      }
      invM.IA[n] = n;
      dcsr_rap(A->blocks[2], &invM, A->blocks[1], &BTB);
      dcsr_add(&BTB, 1.0, A->blocks[3], -1.0, &mgl[1][0].A);
      mgl[1][0].b=dvec_create(A->blocks[2]->row); mgl[1][0].x=dvec_create(A->blocks[2]->row);

      switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
          status = amg_setup_ua(mgl[1], amgparam);
          break;

        case SA_AMG: // Smoothed Aggregation AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
          status = amg_setup_sa(mgl[1], amgparam);
          break;

        default: // Unsmoothed Aggregation AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
          status = amg_setup_ua(mgl[1], amgparam);
          break;

      }

      dcsr_free(&BTB);
      dcsr_free(&invM);

  }
  else
  {
      check_error(ERROR_SOLVER_PRECTYPE, __FUNCTION__);
  }

  /*----------------------------*/
  /* Setup  preconditioner data */
  /*----------------------------*/
  // overall
  precond_block_data precdata;
  precond_block_data_null(&precdata);

  precdata.Abcsr = A;
  precdata.r = dvec_create(b->row);
  precdata.amgparam = amgparam;
  precdata.mgl = mgl;
  precdata.diag = diag;

  if (precond_type > 29 && precond_type < 50)
  {
    // scale el_vol
    dvec_ax(1./itparam->AL_scaling_param, el_vol);
    precdata.el_vol = el_vol;

    precdata.hxdivdata = hxdivdata;
  }
  else
  {
    precdata.el_vol = el_vol;
  }

  precond prec; prec.data = &precdata;

  switch (precond_type)
  {
    case 10:
      prec.fct = precond_block_diag_mixed_darcy;
      break;

    case 11:
      prec.fct = precond_block_lower_mixed_darcy;
      break;

    case 12:
      prec.fct = precond_block_upper_mixed_darcy;
      break;

    case 20:
      prec.fct = precond_block_diag_mixed_darcy_krylov;
      break;

    case 21:
      prec.fct = precond_block_lower_mixed_darcy_krylov;
      break;

    case 22:
      prec.fct = precond_block_upper_mixed_darcy_krylov;
      break;

    case 40:
      prec.fct = precond_block_diag_mixed_darcy_krylov_HX;
      break;

    case 41:
      prec.fct = precond_block_lower_mixed_darcy_krylov_HX;
      break;

    case 42:
      prec.fct = precond_block_upper_mixed_darcy_krylov_HX;
      break;

    case 50:
      prec.fct = precond_block_diag_mixed_darcy_lap;
      break;

    case 51:
      prec.fct = precond_block_lower_mixed_darcy_lap;
      break;

    case 52:
      prec.fct = precond_block_upper_mixed_darcy_lap;
      break;

    case 53:
      prec.fct = precond_block_ilu_mixed_darcy_lap;
      break;

    case 60:
      prec.fct = precond_block_diag_mixed_darcy_lap_krylov;
      break;

    case 61:
      prec.fct = precond_block_lower_mixed_darcy_lap_krylov;
      break;

    case 62:
      prec.fct = precond_block_upper_mixed_darcy_lap_krylov;
      break;

    case 63:
      prec.fct = precond_block_ilu_mixed_darcy_lap_krylov;
      break;

    case 73:
      prec.fct = precond_block_ilu_mixed_darcy_graph_lap_krylov;
      break;


    default:
      break;
  }

  if ( prtlvl >= PRINT_MIN ) {
    get_time(&setup_end);
    setup_duration = setup_end - setup_start;
    print_cputime("Setup totally", setup_duration);
  }

  // solver part
  get_time(&solver_start);

  status=solver_bdcsr_linear_itsolver(A,b,x, &prec,itparam);
  //status=solver_bdcsr_linear_itsolver(A,b,x, NULL,itparam);

  get_time(&solver_end);

  solver_duration = solver_end - solver_start;

  if ( prtlvl >= PRINT_MIN ) {
    print_cputime("Krylov method totally", solver_duration);
    printf("**********************************************************\n");
  }

FINISHED:
  // clean
  precond_block_data_free(&precdata, 2, FALSE);
  if (precond_type > 29 && precond_type < 50)
  {
    dcsr_free(&A_div);
    dvec_ax(itparam->AL_scaling_param, el_vol);
  }

  return status;
}


/********************************************************************************************/
/*!
 * \fn INT linear_solver_bdcsr_babuska_block_2 (block_dCSRmat *A, dvector *b, dvector *x,
 *                                              itsolver_param *itparam,
 *                                              AMG_param *amgparam, dCSRmat *A_diag)
 *
 * \brief Solve Ax = b by standard Krylov methods
 *
 * \param A         Pointer to the coeff matrix in block_dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG solvers * \param A_diag    Diagonal blocks of A
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Ana Budisa
 * \date   2020-10-19
 *
 * \note this function needs to be moved!! -- Xiaozhe Hu
 *
 */
INT linear_solver_bdcsr_babuska_block_2(block_dCSRmat *A,
                                        dvector *b,
                                        dvector *x,
                                        dCSRmat *AS,
                                        dCSRmat *MS,
                                        linear_itsolver_param *itparam,
                                        AMG_param *amgparam1,
                                        AMG_param *amgparam2,
                                        REAL s_frac_power,
                                        REAL t_frac_power,
                                        REAL alpha,
                                        REAL beta,
                                        REAL scaling_a,
                                        REAL scaling_m
                                        )
{

    const SHORT prtlvl = itparam->linear_print_level;

    INT i;
    INT status = SUCCESS;
    REAL setup_start, setup_end, setup_duration;
    REAL solver_start, solver_end, solver_duration;

    SHORT max_levels1, max_levels2;
    if(amgparam1) max_levels1 = amgparam1->max_levels;
    if(amgparam2) max_levels2 = amgparam2->max_levels;

    /*------------------------------------------------*/
    /* setup preconditioner */
    /*------------------------------------------------*/
    get_time(&setup_start);

    INT n0 = A->blocks[0]->col;
    INT ns = AS->col;

    //------------------------------------------------
    // compute the rational approximation
    //------------------------------------------------
    // poles and residues of rational approximation
    dvector poles;
    dvector residues;

    // scaling parameters for rational approximation
    REAL scaled_alpha = alpha, scaled_beta = beta;
    // scale alpha = alpha*sa^(-s)*sm^(s-1)
    scaled_alpha = alpha*pow(scaling_a, -s_frac_power)*pow(scaling_m, s_frac_power-1.);
    // scale beta = beta*sa^(-t)*sm^(t-1)
    scaled_beta  = beta*pow(scaling_a, -t_frac_power)*pow(scaling_m, t_frac_power-1.);

    // parameters used in the function
    REAL16 func_param[4];
    func_param[0] = s_frac_power;
    func_param[1] = t_frac_power;
    if (scaled_alpha > scaled_beta)
    {
      func_param[2] = 1.;
      func_param[3] = scaled_beta/scaled_alpha;
    }
    else
    {
      func_param[2] = scaled_alpha/scaled_beta;
      func_param[3] = 1.;
    }

    // parameters used in the AAA algorithm
    REAL xmin_in=0.e0, xmax_in=1.e0;  // interval for x
    INT mbig=(1<<14)+1;  // initial number of points on the interval [x_min, x_max]
    INT mmax_in=(INT )(mbig/2);  // maximal final number of pole + 1
    REAL16 AAA_tol=powl(2e0,-52e0);  // tolerance of the AAA algorithm
    INT k=-22; // k is the number of nodes in the final interpolation after tolerance is achieved or mmax is reached.
    INT print_level=0; // print level for AAA

    // output of the AAA algorithm
    REAL **rpnwf=malloc(5*sizeof(REAL *));  // output of the AAA algorithm.  It contains residues, poles, nodes, weights, function values

    // compute the rational approximation using AAA algorithms
    //    REAL err_max=get_cpzwf(frac_inv, (void *)func_param,	rpnwf, &mbig, &mmax_in, &k, xmin_in, xmax_in, AAA_tol, print_level);
    get_cpzwf(frac_inv, (void *)func_param,rpnwf, &mbig, &mmax_in, &k, xmin_in, xmax_in, AAA_tol, print_level);

    // assign poles and residules
    dvec_alloc(k,  &residues);
    dvec_alloc(k-1, &poles);
    array_cp(k, rpnwf[0], residues.val);
    array_cp(k-1, rpnwf[1], poles.val);
    //------------------------------------------------

    //------------------------------------------------
    // setup AMG preconditioners for the first block
    // and shifted matrices used in rational approximation
    //------------------------------------------------
    INT npoles = poles.row;
    AMG_data **mgl = (AMG_data **)calloc(npoles+1, sizeof(AMG_data *));

    /* first amg data is to set up AMG for the first block */
    mgl[0] = amg_data_create(max_levels1);
    dcsr_alloc(n0, n0, A->blocks[0]->nnz, &(mgl[0][0].A));
    dcsr_cp(A->blocks[0], &(mgl[0][0].A));
    mgl[0][0].b = dvec_create(n0);
    mgl[0][0].x = dvec_create(n0);

    switch (amgparam1->AMG_type) {

      case UA_AMG: // Unsmoothed Aggregation AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
          status = amg_setup_ua(mgl[0], amgparam1);
          break;

      case SA_AMG: // Smoothed Aggregation AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
          status = amg_setup_sa(mgl[0], amgparam1);
          break;

      default: // UA AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
          status = amg_setup_ua(mgl[0], amgparam1);
          break;

    }

    /* get the scaled mass matrix  */
    // scaling mass matrix
    dCSRmat scaled_M;
    dcsr_alloc(ns, ns, MS->nnz, &scaled_M);
    dcsr_cp(MS, &scaled_M);
    dcsr_axm(&scaled_M, scaling_m);

    // get diagonal entries of the scaled mass matrix
    dvector diag_scaled_M;
    dcsr_getdiag(0, &scaled_M, &diag_scaled_M);

    /* second amg data is to set up all AMG for shifted laplacians */
    // assemble all amg data for all shifted laplacians (scaling_a*A - poles[i] * scaling_m*M)
    //dCSRmat IS = dcsr_create_identity_matrix(ns, 0);
    for(i = 1; i < npoles+1; ++i) {
        mgl[i] = amg_data_create(max_levels2);
        dcsr_alloc(ns, ns, 0, &(mgl[i][0].A));
        dcsr_add(AS, scaling_a, MS, -poles.val[i-1]*scaling_m, &(mgl[i][0].A));
        //dcsr_alloc(ns, ns, MS->nnz, &(mgl[i][0].M)); // TODO: edit this row and one below
        //dcsr_cp(MS, &(mgl[i][0].M));
        mgl[i][0].b = dvec_create(ns);
        mgl[i][0].x = dvec_create(ns);

        switch (amgparam2->AMG_type) {

          case UA_AMG: // Unsmoothed Aggregation AMG
              if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
              status = amg_setup_ua(mgl[i], amgparam2);
              break;

          case SA_AMG: // Smoothed Aggregation AMG
              if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
              status = amg_setup_sa(mgl[i], amgparam2);
              break;

          default: // UA AMG
              if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
              status = amg_setup_ua(mgl[i], amgparam2);
              break;

        }

    }

    //------------------------------------------------
    // setup preconditioner data
    //------------------------------------------------
    precond_block_data precdata;
    precond_block_data_null(&precdata);

    // precdata.Abcsr = A;
    precdata.Abcsr = (block_dCSRmat*)calloc(1, sizeof(block_dCSRmat));
    bdcsr_alloc(A->brow, A->bcol, precdata.Abcsr);
    bdcsr_cp(A, precdata.Abcsr);

    precdata.r = dvec_create(b->row);
    if (amgparam1 && amgparam2) {
        precdata.amgparam = (AMG_param *)calloc(2, sizeof(AMG_param));
        param_amg_cp(amgparam1, &(precdata.amgparam[0]));
        param_amg_cp(amgparam2, &(precdata.amgparam[1]));
    }

    precdata.mgl = mgl;

    // save scaled Mass matrix
    precdata.scaled_M = &scaled_M;
    precdata.diag_scaled_M = &diag_scaled_M;

    // save scaled alpha and beta
    precdata.scaled_alpha = scaled_alpha;
    precdata.scaled_beta = scaled_beta;

    // save poles and residues
    precdata.poles = &poles;
    precdata.residues = &residues;

    precond pc;
    pc.data = &precdata;

    switch (itparam->linear_precond_type)
    {

      case 20:
        pc.fct = precond_block2_babuska_diag;
        break;

      case 21:
        pc.fct = precond_block2_babuska_lower;
        break;

      case 22:
        pc.fct = precond_block2_babuska_upper;
        break;

      default:
        break;
    }

    if ( prtlvl >= PRINT_MIN ) {
        get_time(&setup_end);
        setup_duration = setup_end - setup_start;
        print_cputime("Setup totally", setup_duration);
    }

    //------------------------------------------------
    // solver part
    //------------------------------------------------
    get_time(&solver_start);

    status = solver_bdcsr_linear_itsolver(A, b, x, &pc, itparam);

    get_time(&solver_end);

    solver_duration = solver_end - solver_start;

    if ( prtlvl >= PRINT_MIN ) {
        print_cputime("Krylov method totally", solver_duration);
        printf("**********************************************************\n");
    }

    // -----------------------------------------------
    // clean
    // -----------------------------------------------
    // free rational approximation part
    //if (func_param) free(func_param);
    if (rpnwf[0]) free(rpnwf[0]);
    if (rpnwf) free(rpnwf);
    dvec_free(&poles);
    dvec_free(&residues);

    // free data for preconditioner
    /*
    dcsr_free(&scaled_M);
    dvec_free(&diag_scaled_M);

    bdcsr_free(precdata.Abcsr);
    //dvec_free(&precdata.r);

    amg_data_free(precdata.mgl[0], &(precdata.amgparam[0]));
    for(i = 1; i < npoles+1; ++i) amg_data_free(precdata.mgl[i], &(precdata.amgparam[1])); //free(precdata.mgl[0]);
    if (precdata.mgl) free(precdata.mgl);
    precdata.mgl = NULL;

    // free block data
    precond_block_data_free(&precdata, 2, FALSE);
    */

    return status;
}



/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
