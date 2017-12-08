/*! \file src/solver/itsolver.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 10/06/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note  Done cleanup for releasing -- Xiaozhe Hu 03/12/2017
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
        printf("Number of iterations = %d with relative residual %e.\n", iter, relres);
    }
}

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
 * \fn INT solver_general_linear_itsolver (dCSRmat *A, dvector *b, dvector *x,
 *                                    precond *pc, linear_itsolver_param *itparam)
 *
 * \brief Solve Ax=b by preconditioned Krylov methods (matrix-free version)
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
    
    //dvec_write("diag.dat", &diag);
    
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
        printf("**********************************************************\n");
    }
    
FINISHED:
    amg_data_free(mgl, amgparam);
    
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
            status = amg_setup_ua(mgl_vgrad, amgparam); break;
            
        default: // Classical AMG
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

#if WITH_SUITESPARSE  
  INT i;
#endif
  INT status = SUCCESS;
  REAL setup_start, setup_end, setup_duration;
  REAL solver_start, solver_end, solver_duration;
    
#if WITH_SUITESPARSE
  void **LU_diag = (void **)calloc(2, sizeof(void *));
#else
    error_extlib(256, __FUNCTION__, "SuiteSparse");
#endif
  
  /* setup preconditioner */
  get_time(&setup_start);
  
  /* diagonal blocks are solved exactly */
#if WITH_SUITESPARSE
  // Need to sort the diagonal blocks for UMFPACK format
  dCSRmat A_tran;
  
  for (i=0; i<2; i++){
    
    dcsr_trans(&A_diag[i], &A_tran);
    dcsr_cp(&A_tran, &A_diag[i]);
    
    printf("Factorization for %d-th diagnol: \n", i);
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
      prec.fct = precond_block_diag_2;
      break;
      
    case 11:
      prec.fct = precond_block_lower_2;
      break;
      
    case 12:
      prec.fct = precond_block_upper_2;
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
  precond_block_data_free(&precdata, 2);

  return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_bdcsr_krylov_block_3 (block_dCSRmat *A, dvector *b, dvector *x,
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
    
#if WITH_SUITESPARSE  
    INT i;
#endif
    INT status = SUCCESS;
    REAL setup_start, setup_end, setup_duration;
    REAL solver_start, solver_end, solver_duration;
        
#if WITH_SUITESPARSE
    void **LU_diag = (void **)calloc(3, sizeof(void *));
#else
    error_extlib(257, __FUNCTION__, "SuiteSparse");
#endif
    
    /* setup preconditioner */
    get_time(&setup_start);
    
    /* diagonal blocks are solved exactly */        
#if WITH_SUITESPARSE
        // Need to sort the diagonal blocks for UMFPACK format
        dCSRmat A_tran;
        
        for (i=0; i<3; i++){
            
            dcsr_trans(&A_diag[i], &A_tran);
            dcsr_cp(&A_tran, &A_diag[i]);
            
            printf("Factorization for %d-th diagnol: \n", i);
            LU_diag[i] = umfpack_factorize(&A_diag[i], prtlvl);

            dcsr_free(&A_tran);
            
        }   
#endif
            
    precond_block_data precdata;
    precond_block_data_null(&precdata);

    precdata.Abcsr = A;

    precdata.A_diag = A_diag;
    precdata.r = dvec_create(b->row);
    
    /* diagonal blocks are solved exactly */
#if WITH_SUITESPARSE
        precdata.LU_diag = LU_diag;
#endif
    
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
    precond_block_data_free(&precdata, 3);
    
    return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_bdcsr_krylov_block_4 (block_dCSRmat *A, dvector *b, dvector *x,
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

#if WITH_SUITESPARSE  
    INT i;
#endif
    INT status = SUCCESS;
    REAL setup_start, setup_end, setup_duration;
    REAL solver_start, solver_end, solver_duration;

#if WITH_SUITESPARSE
    void **LU_diag = (void **)calloc(4, sizeof(void *));
#else
    error_extlib(258, __FUNCTION__, "SuiteSparse");
#endif

    /* setup preconditioner */
    get_time(&setup_start);

    /* diagonal blocks are solved exactly */
#if WITH_SUITESPARSE
        // Need to sort the diagonal blocks for UMFPACK format
        dCSRmat A_tran;

        for (i=0; i<4; i++){

            dcsr_trans(&A_diag[i], &A_tran);
            dcsr_cp(&A_tran, &A_diag[i]);

            printf("Factorization for %d-th diagnol: \n", i);
            LU_diag[i] = umfpack_factorize(&A_diag[i], prtlvl);

            dcsr_free(&A_tran);


        }


#endif

    precond_block_data precdata;
    precond_block_data_null(&precdata);

    precdata.Abcsr = A;
    precdata.A_diag = A_diag;
    precdata.r = dvec_create(b->row);

    /* diagonal blocks are solved exactly */
#if WITH_SUITESPARSE
        precdata.LU_diag = LU_diag;
#endif

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
    precond_block_data_free(&precdata, 4);

    return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_bdcsr_krylov_block (block_dCSRmat *A, dvector *b, dvector *x,
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
 * \date   04/05/2017
 *
 * \note  works for general block dCSRmat problems
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

    printf("Factorization for %d-th diagnol: \n", i);
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
  precond_block_data_free(&precdata, nb);

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
 */
INT linear_solver_bdcsr_krylov_mixed_darcy(block_dCSRmat *A,
                                           dvector *b,
                                           dvector *x,
                                           linear_itsolver_param *itparam,
                                           AMG_param *amgparam,
                                           dvector *el_vol)
{  
  const SHORT prtlvl = itparam->linear_print_level;
  const SHORT precond_type = itparam->linear_precond_type;
  
  INT status = SUCCESS;
  REAL setup_start, setup_end, setup_duration;
  REAL solver_start, solver_end, solver_duration;
  
  const SHORT max_levels = amgparam->max_levels;
  INT i, n;
  
  AMG_data **mgl = (AMG_data **)calloc(2, sizeof(AMG_data *));
  for (i=0; i<2; i++) mgl[i]=NULL;
  dvector **diag = (dvector **)calloc(2, sizeof(dvector *));
  for (i=0; i<2; i++) diag[i]=NULL;

  dCSRmat BTB;

  /* setup preconditioner */
  get_time(&setup_start);

  /*----------------------*/
  /* Use argumented Lagrange type preconditioner */
  /*----------------------*/
  if (precond_type < 39)
  {
      /* set AMG for the flux block */
      mgl[0] = amg_data_create(max_levels);
      n = A->blocks[0]->row;
      dcsr_mxm(A->blocks[1],A->blocks[2], &BTB);      
      dcsr_add(&BTB, 1000.0, A->blocks[0], 1.0, &mgl[0][0].A);
      mgl[0][0].b=dvec_create(n); mgl[0][0].x=dvec_create(n);
  
      switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
          status = amg_setup_ua(mgl[0], amgparam);
          break;

        default: // Unsmoothed Aggregation AMG
          if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
          status = amg_setup_ua(mgl[0], amgparam);
          break;

      }

      dcsr_free(&BTB);

  }
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

  precond_block_data precdata;
  precond_block_data_null(&precdata);

  precdata.Abcsr = A;
  precdata.r = dvec_create(b->row);
  precdata.amgparam = amgparam;
  precdata.mgl = mgl;
  precdata.diag = diag;
  precdata.el_vol = el_vol;

  precond prec; prec.data = &precdata;
  
  switch (precond_type)
  {
    case 20:
      prec.fct = precond_block_diag_mixed_darcy;
      break;
      
    case 21:
      prec.fct = precond_block_lower_mixed_darcy;
      break;
      
    case 22:
      prec.fct = precond_block_upper_mixed_darcy;
      break;
      
    case 30:
      prec.fct = precond_block_diag_mixed_darcy_krylov;
      break;
      
    case 31:
      prec.fct = precond_block_lower_mixed_darcy_krylov;
      break;
      
    case 32:
      prec.fct = precond_block_upper_mixed_darcy_krylov;
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

  // clean
  precond_block_data_free(&precdata, 2);

  return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_bdcsr_krylov_biot_2phase (block_dCSRmat *A, dvector *b, dvector *x,
 *                                           itsolver_param *itparam,
 *                                           AMG_param *amgparam, dCSRmat *Mp)
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
 * \date   01/14/2017
 *
 * \note only works for 2 by 2 block dCSRmat problems!! -- Xiaozhe Hu
 */
INT linear_solver_bdcsr_krylov_biot_2phase(block_dCSRmat *A,
                                           dvector *b,
                                           dvector *x,
                                           linear_itsolver_param *itparam,
                                           AMG_param *amgparam,
                                           dCSRmat *Mp)
{
  const SHORT prtlvl = itparam->linear_print_level;
  const SHORT precond_type = itparam->linear_precond_type;

  INT status = SUCCESS;
  REAL setup_start, setup_end, setup_duration;
  REAL solver_start, solver_end, solver_duration;

  const SHORT max_levels = amgparam->max_levels;

  AMG_data **mgl = (AMG_data **)calloc(2, sizeof(AMG_data *));

  /* setup preconditioner */
  get_time(&setup_start);

  /* set AMG for the displacement block */
  mgl[0] = amg_data_create(max_levels);
  dcsr_alloc(A->blocks[0]->row, A->blocks[0]->row, A->blocks[0]->nnz, &mgl[0][0].A);
  dcsr_cp(A->blocks[0], &mgl[0][0].A);
  mgl[0][0].b=dvec_create(A->blocks[0]->row);
  mgl[0][0].x=dvec_create(A->blocks[0]->row);

  switch (amgparam->AMG_type) {

    case UA_AMG: // Unsmoothed Aggregation AMG
      if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
      status = amg_setup_ua(mgl[0], amgparam);
      break;

    default: // UA AMG
      if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
      status = amg_setup_ua(mgl[0], amgparam);
    break;

  }

  /* set AMG for the presssure block */
  mgl[1] = amg_data_create(max_levels);
  dcsr_add(A->blocks[3], 1.0, Mp, 1.0, &mgl[1][0].A);
  mgl[1][0].b=dvec_create(A->blocks[3]->row);
  mgl[1][0].x=dvec_create(A->blocks[3]->row);

  switch (amgparam->AMG_type) {

    case UA_AMG: // Unsmoothed Aggregation AMG
      if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
      status = amg_setup_ua(mgl[1], amgparam);
    break;

    default: // UA AMG
      if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
      status = amg_setup_ua(mgl[1], amgparam); break;
  }

  /* set the whole preconditioner data */
  precond_block_data precdata;
  precond_block_data_null(&precdata);

  precdata.Abcsr = A;
  precdata.r = dvec_create(b->row);
  precdata.amgparam = amgparam;
  precdata.mgl = mgl;

  precond prec; prec.data = &precdata;

  switch (precond_type)
  {
    case 20:
      prec.fct = precond_block_diag_biot_2phase;
      break;

    case 21:
      prec.fct = precond_block_lower_biot_2phase;
      break;

    case 22:
      prec.fct = precond_block_upper_biot_2phase;
      break;

    case 30:
      prec.fct = precond_block_diag_biot_2phase_krylov;
      break;

    case 31:
      prec.fct = precond_block_lower_biot_2phase_krylov;
      break;

    case 32:
      prec.fct = precond_block_upper_biot_2phase_krylov;
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

  // clean
  precond_block_data_free(&precdata, 2);

  return status;
}


/********************************************************************************************/
/**
 * \fn INT linear_solver_bdcsr_krylov_maxwell (block_dCSRmat *A, dvector *b, dvector *x,
 *                                           itsolver_param *itparam,
 *                                           AMG_param *amgparam, dCSRmat *A_diag,
 *                                           dCSRmat *P_curl, dCSRmat *Grad)
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
 * \date   02/25/2014
 *
 * \note Xiaozhe: this is a special solver only for time dependent Maxwell problem invloving
 *       magnetic filed B, electrical field E, and pressure p
 *
 * \note Xiaozhe: A = [A_BB  A_BE  A_Bp
 *                     A_EB  A_EE  A_Ep
 *                     A_pB  A_pE  A_pp]
 */
INT linear_solver_bdcsr_krylov_maxwell(block_dCSRmat *A,
                                       dvector *b,
                                       dvector *x,
                                       linear_itsolver_param *itparam,
                                       AMG_param *amgparam,
                                       dCSRmat *A_diag,
                                       dCSRmat *P_curl,
                                       dCSRmat *Grad,
                                       dCSRmat *Gb,
                                       dCSRmat *Kb,
                                       dCSRmat *Gtb,
                                       dCSRmat *Ktb)
{    
    const SHORT prtlvl = itparam->linear_print_level;
    const SHORT precond_type = itparam->linear_precond_type;
    

    INT status = SUCCESS;
    REAL setup_start, setup_end, setup_duration;
    REAL solver_start, solver_end, solver_duration;
    
    const SHORT max_levels = amgparam->max_levels;
    INT m, n, nnz;
    
    void **LU_diag = (void **)calloc(3, sizeof(void *));
    dvector **diag = (dvector **)calloc(3, sizeof(dvector *));
    AMG_data **mgl = (AMG_data **)calloc(3, sizeof(AMG_data *));
    HX_curl_data **hxcurldata = (HX_curl_data **)calloc(3, sizeof(HX_curl_data *));
        
    dCSRmat Pt_curl;
    dCSRmat Gradt;
    dCSRmat A_vgrad;
    dCSRmat A_grad;
#if WITH_SUITESPARSE
    dCSRmat A_tran;
#endif
    
    AMG_data *mgl_vgrad;
    AMG_data *mgl_grad;
    
    /*------------------------*/
    /* setup preconditioner */
    /*------------------------*/
    get_time(&setup_start);
    
    /*--------------------------------------------------------------------- */
    /* magnetic filed block A_BB */
    /*--------------------------------------------------------------------- */
    if ( prtlvl >= PRINT_MIN ){
        printf("\n*******************************************\n");
        printf("Setup for magnetic diagonal block: \n");
        printf("*******************************************\n");
    }
    
    if (precond_type < 20){
        
#if WITH_SUITESPARSE
        dcsr_trans(&A_diag[0], &A_tran);
        dcsr_cp(&A_tran, &A_diag[0]);
        if ( prtlvl >= PRINT_MIN )
            printf("Factorization for magnetic diagonal block: \n");
        LU_diag[0] = umfpack_factorize(&A_diag[0], prtlvl);
#else
    error_extlib(259, __FUNCTION__, "SuiteSparse");
#endif
        
    }
    else {
        
        // AMG for A_BB
        /*
        mgl[0] = amg_data_create(max_levels);
        m = A_diag[0].row; n = A_diag[0].col; nnz = A_diag[0].nnz;
        mgl[0][0].A=dcsr_create(m,n,nnz); dcsr_cp(&A_diag[0],&mgl[0][0].A);
        mgl[0][0].b=dvec_create(n); mgl[0][0].x=dvec_create(n);
        
        switch (amgparam->AMG_type) {
            case UA_AMG: // Unsmoothed Aggregation AMG
                status = amg_setup_ua(mgl[0], amgparam); break;
            default: // Classical AMG
                status = amg_setup_c(mgl[0], amgparam); break;
        }
     
        if (status < 0) goto FINISHED;
         */
        
        // diagonal precondition for A_BB
        diag[0] = (dvector *)calloc(1, sizeof(dvector));
        dcsr_getdiag(0, &A_diag[0], diag[0]);
        
    }
    
    /*--------------------------------------------------------------------- */
    /* electrical field A_EE */
    /*--------------------------------------------------------------------- */
    if ( prtlvl >= PRINT_MIN ) {
        printf("\n*******************************************\n");
        printf("Setup for electrical diagonal block: \n");
        printf("*******************************************\n");
    }
    
    if (precond_type < 20){
        
#if WITH_SUITESPARSE
        // direct solver for A_EE
        dcsr_trans(&A_diag[1], &A_tran);
        dcsr_cp(&A_tran, &A_diag[1]);
        if ( prtlvl >= PRINT_MIN )
            printf("Factorization for electrical diagonal block: \n");
        LU_diag[1] = umfpack_factorize(&A_diag[1], prtlvl);
#endif
        
    }
    else{
    
        // HX for A_EE
        /*------------------------*/
        /* setup vector Laplacian */
        /*------------------------*/
        // get transpose of P_curl
        dcsr_trans(P_curl, &Pt_curl);
        
        // get A_vgrad
        dcsr_rap(&Pt_curl, &A_diag[1], P_curl, &A_vgrad);
        
        // initialize A, b, x for mgl_vgrad[0]
        mgl_vgrad = amg_data_create(max_levels);
        mgl_vgrad[0].A=dcsr_create(A_vgrad.row,A_vgrad.col,A_vgrad.nnz);
        dcsr_cp(&A_vgrad, &mgl_vgrad[0].A);
        mgl_vgrad[0].b=dvec_create(A_vgrad.col);
        mgl_vgrad[0].x=dvec_create(A_vgrad.row);
        
        // setup AMG for vector Laplacian
        switch (amgparam->AMG_type) {
                
            case UA_AMG: // Unsmoothed Aggregation AMG
                if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                status = amg_setup_ua(mgl_vgrad, amgparam);
            break;
                
            default: // UA AMG
                if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                status = amg_setup_ua(mgl_vgrad, amgparam);
            break;
                
        }
        
        if (status < 0) goto FINISHED;
        
        /*------------------------*/
        /* setup scalar Laplacian */
        /*------------------------*/
        
        // get transpose of Grad
        dcsr_trans(Grad, &Gradt);
        
        // get A_grad
        dcsr_rap(&Gradt, &A_diag[1], Grad, &A_grad);
        
        // initialize A, b, x for mgl_grad[0]
        mgl_grad = amg_data_create(max_levels);
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
                
            default: // UA AMG
                if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                status = amg_setup_ua(mgl_grad, amgparam);
            break;
                
        }
        
        if (status < 0) goto FINISHED;
        
        // setup HX data
        hxcurldata[1] = (HX_curl_data *)calloc(1, sizeof(HX_curl_data));

        hxcurldata[1]->A = &A_diag[1];
        
        hxcurldata[1]->smooth_type = 1;
        hxcurldata[1]->smooth_iter = itparam->HX_smooth_iter;
        
        hxcurldata[1]->P_curl = P_curl;
        hxcurldata[1]->Pt_curl = &Pt_curl;
        hxcurldata[1]->A_vgrad = &A_vgrad;
        hxcurldata[1]->amgparam_vgrad = amgparam;
        hxcurldata[1]->mgl_vgrad = mgl_vgrad;
        
        hxcurldata[1]->Grad = Grad;
        hxcurldata[1]->Gradt = &Gradt;
        hxcurldata[1]->A_grad = &A_grad;
        hxcurldata[1]->amgparam_grad = amgparam;
        hxcurldata[1]->mgl_grad = mgl_grad;
        
        hxcurldata[1]->backup_r = (REAL*)calloc(A_diag[1].row, sizeof(REAL));
        hxcurldata[1]->w = (REAL*)calloc(A_diag[1].row, sizeof(REAL));
        
    }
    
    /*--------------------------------------------------------------------- */
    /* pressue field A_pp */
    /*--------------------------------------------------------------------- */
    if ( prtlvl >= PRINT_MIN ) {
        printf("\n*******************************************\n");
        printf("Setup for pressure diagonal block: \n");
        printf("*******************************************\n");
    }
    
    if (precond_type < 20){
        
#if WITH_SUITESPARSE
        // direct solver for A_pp
        dcsr_trans(&A_diag[2], &A_tran);
        dcsr_cp(&A_tran, &A_diag[2]);
        if ( prtlvl >= PRINT_MIN )
        printf("Factorization for pressure diagonal block: \n");
        LU_diag[2] = umfpack_factorize(&A_diag[2], prtlvl);
#endif
        
    }
    else {
        
        // AMG for A_pp
        mgl[2] = amg_data_create(max_levels);
        m = A_diag[2].row; n = A_diag[2].col; nnz = A_diag[2].nnz;
        mgl[2][0].A=dcsr_create(m,n,nnz); dcsr_cp(&A_diag[2],&mgl[2][0].A);
        mgl[2][0].b=dvec_create(n); mgl[2][0].x=dvec_create(n);
    
        switch (amgparam->AMG_type) {

            case UA_AMG: // Unsmoothed Aggregation AMG
                if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                status = amg_setup_ua(mgl[2], amgparam);
            break;

            default: // UA AMG
                if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                status = amg_setup_ua(mgl[2], amgparam);
            break;

        }
    }
    
    /*--------------------------------------------------------------------- */
    // setup precond data
    /*--------------------------------------------------------------------- */
    precond_block_data precdata;
    precond_block_data_null(&precdata);

    precdata.Abcsr = A;

    precdata.A_diag = A_diag;
    precdata.r = dvec_create(b->row);
    
    precdata.G = Gb;
    precdata.K = Kb;
    precdata.Gt = Gtb;
    precdata.Kt = Ktb;
    
    precdata.LU_diag = LU_diag;
    precdata.diag = diag;
    precdata.amgparam = amgparam;
    precdata.mgl = mgl;
    precdata.hxcurldata = hxcurldata;
    
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
            prec.fct = precond_block_diag_maxwell;
            break;
            
        case 21:
            prec.fct = precond_block_lower_maxwell;
            break;
     
        case 22:
            prec.fct = precond_block_upper_maxwell;
            break;
            
        case 30:
            prec.fct = precond_block_diag_maxwell_krylov;
            break;
            
        case 31:
            prec.fct = precond_block_lower_maxwell_krylov;
            break;
            
        case 32:
            prec.fct = precond_block_upper_maxwell_krylov;
            break;
            
        case 41:
            prec.fct = precond_block_lower_diag_maxwell;
            break;
            
        case 42:
            prec.fct = precond_block_diag_upper_maxwell;
            break;
            
        case 43:
            prec.fct = precond_block_lower_diag_upper_maxwell;
            break;
            
        case 51:
            prec.fct = precond_block_lower_diag_maxwell_krylov;
            break;
            
        case 52:
            prec.fct = precond_block_diag_upper_maxwell_krylov;
            break;
            
        case 53:
            prec.fct = precond_block_lower_diag_upper_maxwell_krylov;
            break;
     
        default:
            check_error(ERROR_SOLVER_PRECTYPE, __FUNCTION__);
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
        print_cputime("Krylov method for Maxwell totally", solver_duration);
        printf("**********************************************************\n");
    }

FINISHED:
    // clean
    precond_block_data_free(&precdata,3);
    
    return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_bdcsr_krylov_bubble_stokes (block_dCSRmat *A, dvector *b, dvector *x,
 *                                           itsolver_param *itparam,
 *                                           AMG_param *amgparam, dCSRmat *Mp)
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
 * \date   01/14/2017
 *
 * \note only works for 2 by 2 block dCSRmat problems!! -- Xiaozhe Hu
 */
INT linear_solver_bdcsr_krylov_bubble_stokes(block_dCSRmat *A,
                                           dvector *b,
                                           dvector *x,
                                           linear_itsolver_param *itparam,
                                           AMG_param *amgparam,
                                           dCSRmat *Mp)
{
  const SHORT prtlvl = itparam->linear_print_level;
  const SHORT precond_type = itparam->linear_precond_type;

  INT status = SUCCESS;
  REAL setup_start, setup_end, setup_duration;
  REAL solver_start, solver_end, solver_duration;

  const SHORT max_levels = amgparam->max_levels;

  AMG_data **mgl = (AMG_data **)calloc(2, sizeof(AMG_data *));

  /* setup preconditioner */
  get_time(&setup_start);

  /* set AMG for the displacement block */
  mgl[0] = amg_data_create(max_levels);
  dcsr_alloc(A->blocks[0]->row, A->blocks[0]->row, A->blocks[0]->nnz, &mgl[0][0].A);
  dcsr_cp(A->blocks[0], &mgl[0][0].A);
  mgl[0][0].b=dvec_create(A->blocks[0]->row);
  mgl[0][0].x=dvec_create(A->blocks[0]->row);

  switch (amgparam->AMG_type) {

    case UA_AMG: // Unsmoothed Aggregation AMG
      if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
      status = amg_setup_ua(mgl[0], amgparam);
      break;

    default: // UA AMG
      if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
      status = amg_setup_ua(mgl[0], amgparam);
    break;

  }

  /* set AMG for the presssure block */
  mgl[1] = amg_data_create(max_levels);
  dcsr_alloc(Mp->row,Mp->row,Mp->nnz, &mgl[1][0].A);
  dcsr_cp(Mp, &mgl[1][0].A);
  mgl[1][0].b=dvec_create(Mp->row);
  mgl[1][0].x=dvec_create(Mp->row);

  switch (amgparam->AMG_type) {

    case UA_AMG: // Unsmoothed Aggregation AMG
      if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
      status = amg_setup_ua(mgl[1], amgparam);
    break;

    default: // UA AMG
      if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
      status = amg_setup_ua(mgl[1], amgparam); break;
  }

  /* set the whole preconditioner data */
  precond_block_data precdata;
  precond_block_data_null(&precdata);

  precdata.Abcsr = A;
  precdata.r = dvec_create(b->row);
  precdata.amgparam = amgparam;
  precdata.mgl = mgl;

  precond prec; prec.data = &precdata;


  switch (precond_type)
  {
    case 20:
      prec.fct = precond_block_diag_bubble_stokes;
      break;

    default:
      prec.fct = precond_block_diag_bubble_stokes;
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

  // clean
  precond_block_data_free(&precdata, 2);

  return status;
}

/********************************************************************************************/
/**
 * \fn INT linear_solver_bdcsr_biot_bubble (block_dCSRmat *A, 
 *                                                   dvector *b, 
 *                                                   dvector *x,
 *                                                   itsolver_param *itparam,
 *                                                   AMG_param *amgparam,
 *                                                   dCSRmat *Mp)
 *
 * \brief Solve Ax = b by standard Krylov methods or AMG
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
 * \date   01/14/2017
 *
 * \note only works for 3 by 3 block dCSRmat problems!! -- Peter Ohm
 */
INT linear_solver_bdcsr_biot_bubble(block_dCSRmat *A,
                                           dvector *b,
                                           dvector *x,
                                           linear_itsolver_param *itparam,
                                           AMG_param *amgparam,
                                           dCSRmat * A_diag)
{
  const SHORT prtlvl = itparam->linear_print_level;
  const SHORT precond_type = itparam->linear_precond_type;

  INT status = SUCCESS;
  REAL setup_start, setup_end, setup_duration;
  REAL solver_start, solver_end, solver_duration;

  const SHORT max_levels = amgparam->max_levels;

  AMG_data **mgl = (AMG_data **)calloc(3, sizeof(AMG_data *));

  /* setup preconditioner */
  get_time(&setup_start);

  /* set AMG for the displacement block */
  mgl[0] = amg_data_create(max_levels);
  dcsr_alloc(A_diag[0].row, A_diag[0].row, A_diag[0].nnz, &mgl[0][0].A);
  dcsr_cp(&(A_diag[0]), &mgl[0][0].A);
  mgl[0][0].b=dvec_create(A_diag[0].row);
  mgl[0][0].x=dvec_create(A_diag[0].row);

  switch (amgparam->AMG_type) {

    case UA_AMG: // Unsmoothed Aggregation AMG
      if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
      status = amg_setup_ua(mgl[0], amgparam);
      break;

    default: // UA AMG
      if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
      status = amg_setup_ua(mgl[0], amgparam);
    break;

  }
  printf("finished displacement block\n");

  /* set AMG for the darcy block */
  mgl[1] = amg_data_create(max_levels);
  dcsr_alloc(A_diag[1].row, A_diag[1].row, A_diag[1].nnz, &mgl[1][0].A);
  dcsr_cp(&(A_diag[1]), &mgl[1][0].A);
  mgl[1][0].b=dvec_create(A_diag[1].row);
  mgl[1][0].x=dvec_create(A_diag[1].row);

  switch (amgparam->AMG_type) {

    case UA_AMG: // Unsmoothed Aggregation AMG
      if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
      status = amg_setup_ua(mgl[1], amgparam);
    break;

    default: // UA AMG
      if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
      status = amg_setup_ua(mgl[1], amgparam); break;
  }

  /* set AMG for the presssure block */
  mgl[2] = amg_data_create(max_levels);
  dcsr_alloc(A_diag[2].row, A_diag[2].row, A_diag[2].nnz, &mgl[2][0].A);
  dcsr_cp(&(A_diag[2]), &mgl[2][0].A);
  mgl[2][0].b=dvec_create(A_diag[2].row);
  mgl[2][0].x=dvec_create(A_diag[2].row);

  //switch (amgparam->AMG_type) {

  //  case UA_AMG: // Unsmoothed Aggregation AMG
  //    if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
  //    status = amg_setup_ua(mgl[2], amgparam);
  //  break;

  //  default: // UA AMG
  //    if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
  //    status = amg_setup_ua(mgl[2], amgparam); break;
  //}

  /* set the whole preconditioner data */
  precond_block_data precdata;
  precond_block_data_null(&precdata);

  precdata.Abcsr = A;
  precdata.A_diag = A_diag;
  precdata.r = dvec_create(b->row);
  precdata.amgparam = amgparam;
  precdata.mgl = mgl;

  precond prec; prec.data = &precdata;


  switch (precond_type)
  {
    default:
      prec.fct = precond_block_diag_biot_bubble_krylov;
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

  // clean
  precond_block_data_free(&precdata, 2);

  return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
