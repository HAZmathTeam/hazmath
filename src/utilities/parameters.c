/*! \file src/utilities/parameters.c
 *
 *  Created by James Adler and Xiaozhe Hu on 10/06/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

#include "hazmat.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
void param_input_init (input_param *iniparam)
{
    
    //----------------
    // output flags
    //----------------
    iniparam->print_level              = PRINT_SOME;
    iniparam->output_type              = 0;
    
    //----------------
    // files
    //----------------
    strcpy(iniparam->workdir,"./");
    strcpy(iniparam->inifile,"./input.dat");
    strcpy(iniparam->gridfile,"../grids/2D/unitSQ_hp125.dat");
    
    //--------------------------
    // finite element parameters
    //--------------------------
    // general parameter
    iniparam->dim                      = 2;
    iniparam->nquad                    = 2;
    
    // parameters for H(D) equations
    iniparam->FE_type                  = 1;
    
    // paramters for Stokes/NS equations
    iniparam->FE_type_velocity         = 2;
    iniparam->FE_type_pressure         = 0;
    
    //----------------------------
    // time steppng paramters
    //----------------------------
    iniparam->time_step_type           = 0;
    iniparam->time_steps               = 2;
    iniparam->time_step_size           = 0.01;
    
    //----------------------------
    // nonlinear solver parameters
    //----------------------------
    iniparam->nonlinear_itsolver_maxit = 5;
    iniparam->nonlinear_itsolver_tol   = 1e-8;
    
    //-------------------------
    // linear solver parameters
    //-------------------------
    // Iterative solver
    iniparam->linear_itsolver_type     = SOLVER_VGMRES;
    iniparam->linear_precond_type      = PREC_NULL;
    iniparam->linear_stop_type         = STOP_REL_RES;
    
    // Solver parameters
    iniparam->linear_itsolver_tol      = 1e-6;
    iniparam->linear_itsolver_maxit    = 500;
    iniparam->linear_restart           = 25;
    
    // ILU method parameters
    iniparam->ILU_type                 = ILUt;
    iniparam->ILU_lfil                 = 0;
    iniparam->ILU_droptol              = 0.01;
    iniparam->ILU_relax                = 0;
    iniparam->ILU_permtol              = 0.0;
    
    // AMG method parameters
    iniparam->AMG_type                 = CLASSIC_AMG;
    iniparam->AMG_levels               = 20;
    iniparam->AMG_cycle_type           = V_CYCLE;
    iniparam->AMG_smoother             = SMOOTHER_GS;
    iniparam->AMG_smooth_order         = CF_ORDER;
    iniparam->AMG_presmooth_iter       = 1;
    iniparam->AMG_postsmooth_iter      = 1;
    iniparam->AMG_relaxation           = 1.0;
    iniparam->AMG_coarse_dof           = 500;
    iniparam->AMG_coarse_solver        = 0;
    iniparam->AMG_tol                  = 1e-6;
    iniparam->AMG_maxit                = 1;
    iniparam->AMG_ILU_levels           = 0;
    iniparam->AMG_coarse_scaling       = OFF;
    iniparam->AMG_amli_degree          = 1;
    iniparam->AMG_nl_amli_krylov_type  = 2;
    
    // Classical AMG specific
    iniparam->AMG_coarsening_type      = 1;
    iniparam->AMG_interpolation_type   = 1;
    iniparam->AMG_max_row_sum          = 0.9;
    iniparam->AMG_strong_threshold     = 0.3;
    iniparam->AMG_truncation_threshold = 0.2;
    iniparam->AMG_aggressive_level     = 0;
    iniparam->AMG_aggressive_path      = 1;
    
    // Aggregation AMG specific
    iniparam->AMG_aggregation_type     = VMB;
    iniparam->AMG_quality_bound        = 8.0;
    iniparam->AMG_pair_number          = 2;
    iniparam->AMG_strong_coupled       = 0.08;
    iniparam->AMG_max_aggregation      = 20;

    // HX Preconditioner
    iniparam->HX_smooth_iter           = 1;
    
}


void param_amg_init (AMG_param *amgparam)
{
    /**
     * \fn void param_amg_init (AMG_param *amgparam)
     *
     * \brief Initialize AMG parameters
     *
     * \param amgparam    Parameters for AMG
     *
     * \author Xiaozhe Hu
     * \date   10/06/2015
     */
    
    // General AMG parameters
    amgparam->AMG_type             = CLASSIC_AMG;
    amgparam->print_level          = PRINT_NONE;
    amgparam->maxit                = 1;
    amgparam->tol                  = 1e-6;
    amgparam->max_levels           = 20;
    amgparam->coarse_dof           = 100;
    amgparam->cycle_type           = V_CYCLE;
    amgparam->smoother             = SMOOTHER_GS;
    amgparam->smooth_order         = NO_ORDER;
    amgparam->presmooth_iter       = 1;
    amgparam->postsmooth_iter      = 1;
    amgparam->coarse_solver        = SOLVER_DEFAULT;
    amgparam->relaxation           = 1.0;
    amgparam->polynomial_degree    = 3;
    amgparam->coarse_scaling       = OFF;
    amgparam->amli_degree          = 2;
    amgparam->amli_coef            = NULL;
    amgparam->nl_amli_krylov_type  = SOLVER_GCG;
    
    // Classical AMG specific
    amgparam->coarsening_type      = COARSE_C;
    amgparam->interpolation_type   = INTERP_STD;
    amgparam->max_row_sum          = 0.9;
    amgparam->strong_threshold     = 0.3;
    amgparam->truncation_threshold = 0.2;
    amgparam->aggressive_level     = 0;
    amgparam->aggressive_path      = 1;
    
    // Aggregation AMG specific
    amgparam->aggregation_type     = PAIRWISE;
    amgparam->quality_bound        = 10.0;
    amgparam->pair_number          = 2;
    amgparam->strong_coupled       = 0.25;
    amgparam->max_aggregation      = 20;
    amgparam->tentative_smooth     = 0.67; // important for SA
    amgparam->smooth_filter        = ON;
    
    // ILU smoother parameters
    amgparam->ILU_type             = ILUk;
    amgparam->ILU_levels           = 0;
    amgparam->ILU_lfil             = 0;
    amgparam->ILU_droptol          = 0.001;
    amgparam->ILU_relax            = 0;
    
    // Schwarz smoother parameters
    amgparam->Schwarz_levels       = 0; // how many levels will use Schwarz smoother
    amgparam->Schwarz_mmsize       = 200;
    amgparam->Schwarz_maxlvl       = 3; // block size -- all vertices at distance .le. this
    amgparam->Schwarz_type         = 1;
    amgparam->Schwarz_blksolver    = SOLVER_DEFAULT;
}

void param_ilu_init (ILU_param *iluparam)
{
    /**
     * \fn void param_ilu_init (ILU_param *iluparam)
     *
     * \brief Initialize ILU parameters
     *
     * \param iluparam  Parameters for ILU
     *
     */
    
    iluparam->print_level  = PRINT_NONE;
    iluparam->ILU_type     = ILUt;
    iluparam->ILU_lfil     = 2;
    iluparam->ILU_droptol  = 0.001;
    iluparam->ILU_relax    = 0;
    iluparam->ILU_permtol  = 0.01;
}

void param_linear_solver_init (linear_itsolver_param *itsparam)
{
    itsparam->linear_print_level   = 0;
    itsparam->linear_itsolver_type = SOLVER_VFGMRES;
    itsparam->linear_precond_type  = PREC_NULL;
    itsparam->linear_stop_type     = STOP_REL_RES;
    itsparam->linear_maxit         = 500;
    itsparam->linear_restart       = 500;
    itsparam->linear_tol           = 1e-6;
    
    // HX preconditioner
    itsparam->HX_smooth_iter       = 3;
}

void param_ilu_set (ILU_param *iluparam,
                         input_param *iniparam)
{
    /**
     * \fn void param_ilu_set (ILU_param *iluparam, input_param *iniparam)
     *
     * \brief Set ILU_param with INPUT
     *
     * \param iluparam    Parameters for ILU
     * \param iniparam     Input parameters
     *
     */
    
    iluparam->print_level = iniparam->print_level;
    iluparam->ILU_type    = iniparam->ILU_type;
    iluparam->ILU_lfil    = iniparam->ILU_lfil;
    iluparam->ILU_droptol = iniparam->ILU_droptol;
    iluparam->ILU_relax   = iniparam->ILU_relax;
    iluparam->ILU_permtol = iniparam->ILU_permtol;
}

void param_solver_set (linear_itsolver_param *itsparam,
                            input_param *iniparam)
{
    /**
     * \fn void param_linear_solver_set (linear_itsolver_param *itsparam, input_param *iniparam)
     *
     * \brief Set itsolver_param with INPUT
     *
     * \param itsparam   Parameters for linear iterative solvers
     * \param iniparam    Input parameters
     *
     */
    
    itsparam->linear_print_level    = iniparam->print_level;
    itsparam->linear_itsolver_type  = iniparam->linear_itsolver_type;
    itsparam->linear_stop_type      = iniparam->linear_stop_type;
    itsparam->linear_restart        = iniparam->linear_restart;
    itsparam->linear_precond_type   = iniparam->linear_precond_type;
    
    if ( itsparam->linear_itsolver_type == SOLVER_AMG ) {
        itsparam->linear_tol   = iniparam->AMG_tol;
        itsparam->linear_maxit = iniparam->AMG_maxit;
    }
    else {
        itsparam->linear_tol   = iniparam->linear_itsolver_tol;
        itsparam->linear_maxit = iniparam->linear_itsolver_maxit;
    }
    
    itsparam->HX_smooth_iter = iniparam->HX_smooth_iter;
    
}

void param_amg_set (AMG_param *param,
                    input_param *iniparam)
{
    /**
     * \fn void param_amg_set (AMG_param *param, input_param *iniparam)
     *
     * \brief Set AMG_param from INPUT
     *
     * \param param     Parameters for AMG
     * \param iniparam   Input parameters
     *
     */
    
    param->AMG_type    = iniparam->AMG_type;
    param->print_level = iniparam->print_level;
    
    if (iniparam->linear_itsolver_type == SOLVER_AMG) {
        param->maxit = iniparam->linear_itsolver_maxit;
        param->tol   = iniparam->linear_itsolver_tol;
    }
    /*
    else if (iniparam->linear_itsolver_type == SOLVER_FMG) {
        param->maxit = iniparam->linear_itsolver_maxit;
        param->tol   = iniparam->linear_itsolver_tol;
    }
     */
    else {
        param->maxit = iniparam->AMG_maxit;
        param->tol   = iniparam->AMG_tol;
    }
    
    param->max_levels           = iniparam->AMG_levels;
    param->cycle_type           = iniparam->AMG_cycle_type;
    param->smoother             = iniparam->AMG_smoother;
    param->smooth_order         = iniparam->AMG_smooth_order;
    param->relaxation           = iniparam->AMG_relaxation;
    param->coarse_solver        = iniparam->AMG_coarse_solver;
    //param->polynomial_degree    = iniparam->AMG_polynomial_degree;
    param->presmooth_iter       = iniparam->AMG_presmooth_iter;
    param->postsmooth_iter      = iniparam->AMG_postsmooth_iter;
    param->coarse_solver        = iniparam->AMG_coarse_solver;
    param->coarse_dof           = iniparam->AMG_coarse_dof;
    param->coarse_scaling       = iniparam->AMG_coarse_scaling;
    param->amli_degree          = iniparam->AMG_amli_degree;
    param->amli_coef            = NULL;
    param->nl_amli_krylov_type  = iniparam->AMG_nl_amli_krylov_type;
    
    param->coarsening_type      = iniparam->AMG_coarsening_type;
    param->interpolation_type   = iniparam->AMG_interpolation_type;
    param->strong_threshold     = iniparam->AMG_strong_threshold;
    param->truncation_threshold = iniparam->AMG_truncation_threshold;
    param->max_row_sum          = iniparam->AMG_max_row_sum;
    param->aggressive_level     = iniparam->AMG_aggressive_level;
    param->aggressive_path      = iniparam->AMG_aggressive_path;
    
    param->aggregation_type     = iniparam->AMG_aggregation_type;
    param->pair_number          = iniparam->AMG_pair_number;
    param->quality_bound        = iniparam->AMG_quality_bound;
    param->strong_coupled       = iniparam->AMG_strong_coupled;
    param->max_aggregation      = iniparam->AMG_max_aggregation;
    //param->tentative_smooth     = iniparam->AMG_tentative_smooth;
    //param->smooth_filter        = iniparam->AMG_smooth_filter;
    
    param->ILU_levels           = iniparam->AMG_ILU_levels;
    param->ILU_type             = iniparam->ILU_type;
    param->ILU_lfil             = iniparam->ILU_lfil;
    param->ILU_droptol          = iniparam->ILU_droptol;
    param->ILU_relax            = iniparam->ILU_relax;
    param->ILU_permtol          = iniparam->ILU_permtol;
    
    /*
    param->Schwarz_levels       = iniparam->AMG_Schwarz_levels;
    param->Schwarz_mmsize       = iniparam->Schwarz_mmsize;
    param->Schwarz_maxlvl       = iniparam->Schwarz_maxlvl;
    param->Schwarz_type         = iniparam->Schwarz_type;
    */
}

void param_linear_solver_print (linear_itsolver_param *param)
{
    /**
     * \fn void param_linear_solver_print (linear_itsolver_param *param)
     *
     * \brief Print out itsolver parameters
     *
     * \param param    Paramters for linear iterative solvers
     *
     * \author Xiaozhe Hu
     * \date   10/06/2015
     */
    
    if ( param ) {
        
        printf("\n       Parameters in linear_itsolver_param\n");
        printf("-----------------------------------------------\n");
        
        printf("Solver print level:                %d\n", param->linear_print_level);
        printf("Solver type:                       %d\n", param->linear_itsolver_type);
        printf("Solver precond type:               %d\n", param->linear_precond_type);
        printf("Solver max num of iter:            %d\n", param->linear_maxit);
        printf("Solver tolerance:                  %.2e\n", param->linear_tol);
        printf("Solver stopping type:              %d\n", param->linear_stop_type);
        printf("Solver restart number:             %d\n", param->linear_restart);
        
        printf("HX smooth num of iter:             %d\n", param->HX_smooth_iter);
    
        printf("-----------------------------------------------\n\n");
        
    }
    else {
        printf("### WARNING: param has not been set!\n");
    }
}

void param_ilu_print (ILU_param *param)
{
    /**
     * \fn void param_ilu_print (ILU_param *param)
     *
     * \brief Print out ILU parameters
     *
     * \param param    Parameters for ILU
     *
     */
    
    if ( param ) {
        
        printf("\n       Parameters in ILU_param\n");
        printf("-----------------------------------------------\n");
        printf("ILU print level:                   %d\n",   param->print_level);
        printf("ILU type:                          %d\n",   param->ILU_type);
        printf("ILU level of fill-in:              %d\n",   param->ILU_lfil);
        printf("ILU relaxation factor:             %.4f\n", param->ILU_relax);
        printf("ILU drop tolerance:                %.2e\n", param->ILU_droptol);
        printf("ILU permutation tolerance:         %.2e\n", param->ILU_permtol);
        printf("-----------------------------------------------\n\n");
        
    }
    else {
        printf("### WARNING: param has not been set!\n");
    }
}

void param_amg_print (AMG_param *param)
{
    
    /**
     * \fn void param_amg_print (AMG_param *param)
     *
     * \brief Print out AMG parameters
     *
     * \param param   Parameters for AMG
     *
     * \author Xiaozhe Hu
     * \date   10/06/2015
     */
    
    if ( param ) {
        
        printf("\n       Parameters in AMG_param\n");
        printf("-----------------------------------------------\n");
        
        printf("AMG print level:                   %d\n", param->print_level);
        printf("AMG max num of iter:               %d\n", param->maxit);
        printf("AMG type:                          %d\n", param->AMG_type);
        printf("AMG tolerance:                     %.2e\n", param->tol);
        printf("AMG max levels:                    %d\n", param->max_levels);
        printf("AMG cycle type:                    %d\n", param->cycle_type);
        printf("AMG coarse solver type:            %d\n", param->coarse_solver);
        printf("AMG scaling of coarse correction:  %d\n", param->coarse_scaling);
        printf("AMG smoother type:                 %d\n", param->smoother);
        printf("AMG smoother order:                %d\n", param->smooth_order);
        printf("AMG num of presmoothing:           %d\n", param->presmooth_iter);
        printf("AMG num of postsmoothing:          %d\n", param->postsmooth_iter);
        
        if ( param->smoother == SMOOTHER_SOR  ||
            param->smoother == SMOOTHER_SSOR ||
            param->smoother == SMOOTHER_GSOR ||
            param->smoother == SMOOTHER_SGSOR ) {
            printf("AMG relax factor:                  %.4f\n", param->relaxation);
        }
        
        /*
        if ( param->smoother == SMOOTHER_POLY ) {
            printf("AMG polynomial smoother degree:    %d\n", param->polynomial_degree);
        }
         */
        
        if ( param->cycle_type == AMLI_CYCLE ) {
            printf("AMG AMLI degree of polynomial:     %d\n", param->amli_degree);
        }
        
        if ( param->cycle_type == NL_AMLI_CYCLE ) {
            printf("AMG Nonlinear AMLI Krylov type:    %d\n", param->nl_amli_krylov_type);
        }
        
        switch (param->AMG_type) {
            case CLASSIC_AMG:
                printf("AMG coarsening type:               %d\n", param->coarsening_type);
                printf("AMG interpolation type:            %d\n", param->interpolation_type);
                printf("AMG dof on coarsest grid:          %d\n", param->coarse_dof);
                printf("AMG strong threshold:              %.4f\n", param->strong_threshold);
                printf("AMG truncation threshold:          %.4f\n", param->truncation_threshold);
                printf("AMG max row sum:                   %.4f\n", param->max_row_sum);
                printf("AMG aggressive levels:             %d\n", param->aggressive_level);
                printf("AMG aggressive path:               %d\n", param->aggressive_path);
                break;
                
            default: // SA_AMG or UA_AMG
                printf("Aggregation type:                  %d\n", param->aggregation_type);
                if ( param->aggregation_type == PAIRWISE ) {
                    printf("Aggregation number of pairs:       %d\n", param->pair_number);
                    printf("Aggregation quality bound:         %.2f\n", param->quality_bound);
                }
                if ( param->aggregation_type == VMB ) {
                    printf("Aggregation AMG strong coupling:   %.4f\n", param->strong_coupled);
                    printf("Aggregation AMG max aggregation:   %d\n", param->max_aggregation);
                    printf("Aggregation AMG tentative smooth:  %.4f\n", param->tentative_smooth);
                    printf("Aggregation AMG smooth filter:     %d\n", param->smooth_filter);
                }
                break;
        }
        
        if (param->ILU_levels>0) {
            printf("AMG ILU smoother level:            %d\n", param->ILU_levels);
            printf("AMG ILU type:                      %d\n", param->ILU_type);
            printf("AMG ILU level of fill-in:          %d\n", param->ILU_lfil);
            printf("AMG ILU drop tol:                  %e\n", param->ILU_droptol);
            printf("AMG ILU relaxation:                %f\n", param->ILU_relax);
        }
        
        if (param->Schwarz_levels>0){
            printf("AMG Schwarz smoother level:        %d\n", param->Schwarz_levels);
            printf("AMG Schwarz type:                  %d\n", param->Schwarz_type);
            printf("AMG Schwarz forming block level:   %d\n", param->Schwarz_maxlvl);
            printf("AMG Schwarz maximal block size:    %d\n", param->Schwarz_mmsize);
        }
        
        printf("-----------------------------------------------\n\n");
        
    }
    else {
        printf("### WARNING: param has not been set!\n");
    } // end if (param)
    
}

void param_amg_to_prec (precond_data *pcdata,
                             AMG_param *amgparam)
{
    /**
     * \fn void param_amg_to_prec (precond_data *pcdata, AMG_param *amgparam)
     *
     * \brief Set precond_data with AMG_param
     *
     * \param pcdata      Preconditioning data structure
     * \param amgparam    Parameters for AMG
     *
     */
    
    pcdata->AMG_type            = amgparam->AMG_type;
    pcdata->print_level         = amgparam->print_level;
    pcdata->maxit               = amgparam->maxit;
    pcdata->max_levels          = amgparam->max_levels;
    pcdata->tol                 = amgparam->tol;
    pcdata->cycle_type          = amgparam->cycle_type;
    pcdata->smoother            = amgparam->smoother;
    pcdata->smooth_order        = amgparam->smooth_order;
    pcdata->presmooth_iter      = amgparam->presmooth_iter;
    pcdata->postsmooth_iter     = amgparam->postsmooth_iter;
    pcdata->coarsening_type     = amgparam->coarsening_type;
    pcdata->coarse_solver       = amgparam->coarse_solver;
    pcdata->relaxation          = amgparam->relaxation;
    pcdata->polynomial_degree   = amgparam->polynomial_degree;
    pcdata->coarse_scaling      = amgparam->coarse_scaling;
    pcdata->amli_degree         = amgparam->amli_degree;
    pcdata->amli_coef           = amgparam->amli_coef;
    pcdata->nl_amli_krylov_type = amgparam->nl_amli_krylov_type;
    
}

void param_prec_to_amg (AMG_param *amgparam,
                             precond_data *pcdata)
{
    /**
     * \fn void fasp_param_prec_to_amg (AMG_param *amgparam, precond_data *pcdata)
     *
     * \brief Set AMG_param with precond_data
     *
     * \param amgparam    Parameters for AMG
     * \param pcdata      Preconditioning data structure
     *
     */
    
    amgparam->AMG_type            = pcdata->AMG_type;
    amgparam->print_level         = pcdata->print_level;
    amgparam->cycle_type          = pcdata->cycle_type;
    amgparam->smoother            = pcdata->smoother;
    amgparam->smooth_order        = pcdata->smooth_order;
    amgparam->presmooth_iter      = pcdata->presmooth_iter;
    amgparam->postsmooth_iter     = pcdata->postsmooth_iter;
    amgparam->relaxation          = pcdata->relaxation;
    amgparam->polynomial_degree   = pcdata->polynomial_degree;
    amgparam->coarse_solver       = pcdata->coarse_solver;
    amgparam->coarse_scaling      = pcdata->coarse_scaling;
    amgparam->amli_degree         = pcdata->amli_degree;
    amgparam->amli_coef           = pcdata->amli_coef;
    amgparam->nl_amli_krylov_type = pcdata->nl_amli_krylov_type;
    amgparam->ILU_levels          = pcdata->mgl_data->ILU_levels;
}

void amg_amli_coef (const REAL lambda_max,
                         const REAL lambda_min,
                         const INT degree,
                         REAL *coef)
{
    
    /**
     * \fn void amg_amli_coef (const REAL lambda_max, const REAL lambda_min,
     *                              const INT degree, REAL *coef)
     *
     * \brief Compute the coefficients of the polynomial used by AMLI-cycle
     *
     * \param lambda_max  Maximal lambda
     * \param lambda_min  Minimal lambda
     * \param degree      Degree of polynomial approximation
     * \param coef        Coefficient of AMLI (output)
     *
     * \author Xiaozhe Hu
     * \date   01/23/2011
     */
    
    const REAL mu0 = 1.0/lambda_max, mu1 = 1.0/lambda_min;
    const REAL c = (sqrt(mu0)+sqrt(mu1))*(sqrt(mu0)+sqrt(mu1));
    const REAL a = (4*mu0*mu1)/(c);
    
    const REAL kappa = lambda_max/lambda_min; // condition number
    const REAL delta = (sqrt(kappa) - 1.0)/(sqrt(kappa)+1.0);
    const REAL b = delta*delta;
    
    if (degree == 0) {
        coef[0] = 0.5*(mu0+mu1);
    }
    
    else if (degree == 1) {
        coef[0] = 0.5*c;
        coef[1] = -1.0*mu0*mu1;
    }
    
    else if (degree > 1) {
        INT i;
        
        // allocate memory
        REAL *work = (REAL *)calloc(2*degree-1, sizeof(REAL));
        REAL *coef_k, *coef_km1;
        coef_k = work; coef_km1 = work+degree;
        
        // get q_k
        amg_amli_coef(lambda_max, lambda_min, degree-1, coef_k);
        // get q_km1
        amg_amli_coef(lambda_max, lambda_min, degree-2, coef_km1);
        
        // get coef
        coef[0] = a - b*coef_km1[0] + (1+b)*coef_k[0];
        
        for (i=1; i<degree-1; i++) {
            coef[i] = -b*coef_km1[i] + (1+b)*coef_k[i] - a*coef_k[i-1];
        }
        
        coef[degree-1] = (1+b)*coef_k[degree-1] - a*coef_k[degree-2];
        
        coef[degree] = -a*coef_k[degree-1];
        
        // clean memory
        if (work) free(work);
    }
    
    else {
        printf("### ERROR: Wrong AMLI degree %d!\n", degree);
        check_error(ERROR_INPUT_PAR, __FUNCTION__);
    }
    
    return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
