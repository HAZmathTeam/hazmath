/*
 *  itsolver.c
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

INT solver_dcsr_linear_itsolver (dCSRmat *A,
                               dvector *b,
                               dvector *x,
                               precond *pc,
                               linear_itsolver_param *itparam)
{
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
    
    const SHORT prtlvl        = itparam->linear_print_level;
    const SHORT itsolver_type = itparam->linear_itsolver_type;
    const SHORT stop_type     = itparam->linear_stop_type;
    const SHORT restart       = itparam->linear_restart;
    const INT   MaxIt         = itparam->linear_maxit;
    const REAL  tol           = itparam->linear_tol;
    
    /* Local Variables */
    REAL solver_start, solver_end, solver_duration;
    INT iter;
    
    gettime(&solver_start);
    
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
        gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Iterative method", solver_duration);
    }
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return iter;
}


INT linear_solver_amg (dCSRmat *A,
                      dvector *b,
                      dvector *x,
                      AMG_param *param)
{
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
     * \note Refer to "Multigrid"
     *       by U. Trottenberg, C. W. Oosterlee and A. Schuller
     *       Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben)
     *       Academic Press Inc., San Diego, CA, 2001.
     *
     */
    
    const SHORT   max_levels  = param->max_levels;
    const SHORT   prtlvl      = param->print_level;
    const SHORT   amg_type    = param->AMG_type;
    const SHORT   cycle_type  = param->cycle_type;
    const INT     nnz = A->nnz, m = A->row, n = A->col;
    
    // local variables
    INT      status = SUCCESS;
    AMG_data *    mgl = amg_data_create(max_levels);
    REAL          AMG_start, AMG_end;
    
    if ( prtlvl > PRINT_NONE ) gettime(&AMG_start);
    
    // check matrix data
    if ( m != n ) {
        printf("### ERROR: A is not a square matrix!\n");
        chkerr(ERROR_MAT_SIZE, __FUNCTION__);
    }
    
    if ( nnz <= 0 ) {
        printf("### ERROR: A has no nonzero entries!\n");
        chkerr(ERROR_MAT_SIZE, __FUNCTION__);
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
          
            /*
        case SA_AMG: // Smoothed Aggregation AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\nCalling SA AMG ...\n");
            status = amg_setup_sa(mgl, param); break;
             */
            
        case UA_AMG: // Unsmoothed Aggregation AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\nCalling UA AMG ...\n");
            status = amg_setup_ua(mgl, param); break;
            
        default: // Classical AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\nCalling classical AMG ...\n");
            status = amg_setup_c(mgl, param); break;
            
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
            printf("### WARNING: AMG setup failed!\n");
            printf("### WARNING: Use a backup solver instead.\n");
        }
        status = dcsr_pvgmres (A, b, x, NULL, param->tol, param->maxit,
                                  20, 1, prtlvl);
        
    }
    
    // clean-up memory
    amg_data_free(mgl, param);
    
    // print out CPU time if needed
    if ( prtlvl > PRINT_NONE ) {
        gettime(&AMG_end);
        print_cputime("AMG totally", AMG_end - AMG_start);
        printf("**********************************************************\n");
    }
    
    return status;
}

INT linear_solver_dcsr_krylov (dCSRmat *A,
                             dvector *b,
                             dvector *x,
                             linear_itsolver_param *itparam)
{
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
    const SHORT prtlvl = itparam->linear_print_level;
    
    /* Local Variables */
    INT      status = SUCCESS;
    REAL     solver_start, solver_end, solver_duration;
    
    gettime(&solver_start);
    
    status = solver_dcsr_linear_itsolver(A,b,x,NULL,itparam);
    
    if ( prtlvl >= PRINT_MIN ) {
        gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Krylov method totally", solver_duration);
        printf("**********************************************************");
    }
    
    return status;
}

INT linear_solver_dcsr_krylov_diag (dCSRmat *A,
                                  dvector *b,
                                  dvector *x,
                                  linear_itsolver_param *itparam)
{
    
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
    const SHORT prtlvl = itparam->linear_print_level;
    
    /* Local Variables */
    INT       status = SUCCESS;
    REAL      solver_start, solver_end, solver_duration;
    
    gettime(&solver_start);
    
    // setup preconditioner
    dvector diag; dcsr_getdiag(0,A,&diag);
    
    precond pc;
    pc.data = &diag;
    pc.fct  = precond_diag;
    
    // call iterative solver
    status = solver_dcsr_linear_itsolver(A,b,x,&pc,itparam);
    
    if ( prtlvl >= PRINT_MIN ) {
        gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Diag_Krylov method totally", solver_duration);
        printf("**********************************************************\n");

    }
    
    dvec_free(&diag);
    
    return status;
}


INT linear_solver_dcsr_krylov_ilu (dCSRmat *A,
                                 dvector *b,
                                 dvector *x,
                                 linear_itsolver_param *itparam,
                                 ILU_param *iluparam)
{
    /**
     * \fn INT linear_solver_dcsr_krylov_ilu (dCSRmat *A, dvector *b, dvector *x,
     *                                      itsolver_param *itparam, ILU_param *iluparam)
     *
     * \brief Solve Ax=b by ILUs preconditioned Krylov methods
     *
     * \param A         Pointer to the coeff matrix in dCSRmat format
     * \param b         Pointer to the right hand side in dvector format
     * \param x         Pointer to the approx solution in dvector format
     * \param itparam   Pointer to parameters for iterative solvers
     * \param iluparam  Pointer to parameters for ILU
     *
     * \return          Iteration number if converges; ERROR otherwise.
     *
     * \author Chensong Zhang, Shiquan Zhang
     * \date   09/25/2009
     */
    
    const SHORT prtlvl = itparam->linear_print_level;
    
    /* Local Variables */
    INT      status = SUCCESS;
    REAL     solver_start, solver_end, solver_duration;
    
    gettime(&solver_start);
    
    // ILU setup for whole matrix
    ILU_data LU;
    if ( (status = ilu_dcsr_setup(A,&LU,iluparam)) < 0 ) goto FINISHED;
    
    // check iludata
    if ( (status = mem_iludata_check(&LU)) < 0 ) goto FINISHED;
    
    // set preconditioner
    precond pc;
    pc.data = &LU;
    pc.fct  = precond_ilu;
    
    // call iterative solver
    status = solver_dcsr_linear_itsolver(A,b,x,&pc,itparam);
    
    if ( prtlvl >= PRINT_MIN ) {
        gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        
        switch (iluparam->ILU_type) {
            case ILUt:
                print_cputime("ILUt_Krylov method totally", solver_duration);
                break;
            case ILUtp:
                print_cputime("ILUtp_Krylov method totally", solver_duration);
                break;
            default: // ILUk
                print_cputime("ILUk_Krylov method totally", solver_duration);
                break;
        }
    }
    
FINISHED:
    ilu_data_free(&LU);
    
    return status;
}


INT linear_solver_dcsr_krylov_amg (dCSRmat *A,
                                 dvector *b,
                                 dvector *x,
                                 linear_itsolver_param *itparam,
                                 AMG_param *amgparam)
{
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
    
    const SHORT prtlvl = itparam->linear_print_level;
    const SHORT max_levels = amgparam->max_levels;
    const INT nnz = A->nnz, m = A->row, n = A->col;
    
    /* Local Variables */
    INT      status = SUCCESS;
    REAL     solver_start, solver_end, solver_duration;
    
    gettime(&solver_start);
    
    // initialize A, b, x for mgl[0]
    AMG_data *mgl=amg_data_create(max_levels);
    mgl[0].A=dcsr_create(m,n,nnz); dcsr_cp(A,&mgl[0].A);
    mgl[0].b=dvec_create(n); mgl[0].x=dvec_create(n);
    
    
    // setup preconditioner
    switch (amgparam->AMG_type) {
            
        //case SA_AMG: // Smoothed Aggregation AMG
        //    status = amg_setup_sa(mgl, amgparam); break;
            
        case UA_AMG: // Unsmoothed Aggregation AMG
            status = amg_setup_ua(mgl, amgparam); break;
            
        default: // Classical AMG
            status = amg_setup_c(mgl, amgparam); break;
            
    }
    
    if (status < 0) goto FINISHED;
    
    // setup preconditioner
    precond_data pcdata;
    param_amg_to_prec(&pcdata,amgparam);
    pcdata.max_levels = mgl[0].num_levels;
    pcdata.mgl_data = mgl;
    
    precond pc; pc.data = &pcdata;
    //pc.fct = precond_amg;
    
    //if (itparam->precond_type == PREC_FMG) {
    //    pc.fct = fasp_precond_famg; // Full AMG
    //}
    //else {
        switch (amgparam->cycle_type) {
            case AMLI_CYCLE: // AMLI cycle
                pc.fct = precond_amli; break;
            case NL_AMLI_CYCLE: // Nonlinear AMLI AMG
                pc.fct = precond_nl_amli; break;
            default: // V,W-Cycle AMG
                pc.fct = precond_amg; break;
        }
    //}
    
    // call iterative solver
    status = solver_dcsr_linear_itsolver(A, b, x, &pc, itparam);
    
    if ( prtlvl >= PRINT_MIN ) {
        gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("AMG_Krylov method totally", solver_duration);
        printf("**********************************************************\n");
    }
    
FINISHED:
    amg_data_free(mgl, amgparam);
    

    return status;
}

INT linear_solver_dcsr_krylov_hx_curl (dCSRmat *A,
                                   dvector *b,
                                   dvector *x,
                                   linear_itsolver_param *itparam,
                                   AMG_param *amgparam,
                                   dCSRmat *P_curl,
                                   dCSRmat *Grad)
{
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
    
    const SHORT prtlvl = itparam->linear_print_level;
    const SHORT max_levels = amgparam->max_levels;
    const INT nnz = A->nnz, m = A->row, n = A->col;
    
    /*------------------------*/
    /* Local Variables */
    /*------------------------*/
    INT      status = SUCCESS;
    REAL     solver_start, solver_end, solver_duration;
    
    gettime(&solver_start);
    
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
            status = amg_setup_c(mgl_vgrad, amgparam); break;
            
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
            status = amg_setup_ua(mgl_grad, amgparam); break;
            
        default: // Classical AMG
            status = amg_setup_c(mgl_grad, amgparam); break;
            
    }
    
    if (status < 0) goto FINISHED;
    
    /*------------------------*/
    // setup preconditioner
    HX_curl_data hxcurldata;
    
    hxcurldata.A = A;
    
    hxcurldata.smooth_type = 1;
    hxcurldata.smooth_iter = 1;
    
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
    
    precond pc; pc.data = &hxcurldata;
    switch (itparam->linear_precond_type) {
            
        case PREC_HX_CURL_A:
            pc.fct = precond_hx_curl_additive;
            break;
        
        default:
            pc.fct = precond_hx_curl_multiplicative;
            break;

    }
    
    // call iterative solver
    status = solver_dcsr_linear_itsolver(A, b, x, &pc, itparam);
    
    if ( prtlvl >= PRINT_MIN ) {
        gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("HX_curl_Krylov method totally", solver_duration);
        printf("**********************************************************\n");
    }
    
FINISHED:
    amg_data_free(mgl_vgrad, amgparam);
    amg_data_free(mgl_grad, amgparam);
    
    //todo: clean memory
    
    return status;
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
