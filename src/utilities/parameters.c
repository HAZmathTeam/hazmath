/*! \file src/utilities/parameters.c
 *
 *  Created by James Adler and Xiaozhe Hu on 10/06/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/27/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/28/2016
 *
 */

#include "hazmat.h"

/*************************************************************************************/
void param_input_init (input_param *inparam)
{
    /*!
     * \fn void param_input_init (input_param *inparam)
     *
     * \brief Initialize input parameters
     *
     * \param inparam    Pointer to input_param structure
     *
     */

    //----------------
    // output flags
    //----------------
    inparam->print_level              = PRINT_SOME;

    //----------------
    // files
    //----------------
    strcpy(inparam->inifile,"./input.dat");
    strcpy(inparam->gridfile,"../grids/2D/unitSQ_hp125.dat");
    strcpy(inparam->output_dir,"./output/");
    
    //--------------------------
    // finite element parameters
    //--------------------------
    // general parameter
    inparam->nquad                    = 2;
    
    // parameters for H(D) equations
    inparam->FE_type                  = 1;
    
    //----------------------------
    // time steppng paramters
    //----------------------------
    inparam->time_step_type           = 0;
    inparam->time_steps               = 2;
    inparam->time_step_size           = 0.01;
    
    //----------------------------
    // nonlinear solver parameters
    //----------------------------
    inparam->nonlinear_itsolver_type    = 0;
    inparam->nonlinear_itsolver_maxit   = 6;
    inparam->nonlinear_itsolver_tol     = 1e-6;
    inparam->nonlinear_itsolver_toltype	= 0;
    
    //-------------------------
    // linear solver parameters
    //-------------------------
    // Iterative solver
    inparam->linear_itsolver_type     = SOLVER_VGMRES;
    inparam->linear_precond_type      = PREC_NULL;
    inparam->linear_stop_type         = STOP_REL_RES;
    
    // Solver parameters
    inparam->linear_itsolver_tol      = 1e-6;
    inparam->linear_itsolver_maxit    = 500;
    inparam->linear_restart           = 25;
    
    // ILU method parameters
    inparam->ILU_type                 = ILUt;
    inparam->ILU_lfil                 = 0;
    inparam->ILU_droptol              = 0.01;
    inparam->ILU_relax                = 0;
    inparam->ILU_permtol              = 0.0;
    
    // AMG method parameters
    inparam->AMG_type                 = UA_AMG;
    inparam->AMG_levels               = 10;
    inparam->AMG_cycle_type           = V_CYCLE;
    inparam->AMG_smoother             = SMOOTHER_GS;
    inparam->AMG_smooth_order         = NO_ORDER;
    inparam->AMG_presmooth_iter       = 1;
    inparam->AMG_postsmooth_iter      = 1;
    inparam->AMG_polynomial_degree    = 2;
    inparam->AMG_relaxation           = 1.2;
    inparam->AMG_coarse_dof           = 100;
    inparam->AMG_coarse_solver        = SOLVER_DEFAULT;
    inparam->AMG_tol                  = 1e-6;
    inparam->AMG_maxit                = 1;
    inparam->AMG_ILU_levels           = 0;
    inparam->AMG_coarse_scaling       = OFF;
    inparam->AMG_amli_degree          = 1;
    inparam->AMG_nl_amli_krylov_type  = 2;
    inparam->AMG_Schwarz_levels       = 0;
    
    // Classical AMG specific
    inparam->AMG_coarsening_type      = 1;
    inparam->AMG_interpolation_type   = 1;
    inparam->AMG_max_row_sum          = 0.9;
    inparam->AMG_strong_threshold     = 0.3;
    inparam->AMG_truncation_threshold = 0.2;
    inparam->AMG_aggressive_level     = 0;
    inparam->AMG_aggressive_path      = 1;
    
    // Aggregation AMG specific
    inparam->AMG_aggregation_type     = VMB;
    inparam->AMG_quality_bound        = 8.0;
    inparam->AMG_pair_number          = 2;
    inparam->AMG_strong_coupled       = 0.00;
    inparam->AMG_max_aggregation      = 20;

    // HX Preconditioner
    inparam->HX_smooth_iter           = 1;
    
}

/*************************************************************************************/
void param_amg_init (AMG_param *amgparam)
{
    /*!
     * \fn void param_amg_init (AMG_param *amgparam)
     *
     * \brief Initialize AMG parameters
     *
     * \param amgparam    Pointer to AMG_param structure
     *
     * \author Xiaozhe Hu
     * \date   10/06/2015
     *
     */
    
    // General AMG parameters
    amgparam->AMG_type             = UA_AMG;
    amgparam->print_level          = PRINT_NONE;
    amgparam->maxit                = 1;
    amgparam->tol                  = 1e-6;
    amgparam->max_levels           = 10;
    amgparam->coarse_dof           = 100;
    amgparam->cycle_type           = V_CYCLE;
    amgparam->smoother             = SMOOTHER_GS;
    amgparam->smooth_order         = NO_ORDER;
    amgparam->presmooth_iter       = 1;
    amgparam->postsmooth_iter      = 1;
    amgparam->coarse_solver        = SOLVER_DEFAULT;
    amgparam->relaxation           = 1.0;
    amgparam->polynomial_degree    = 2;
    amgparam->coarse_scaling       = OFF;
    amgparam->amli_degree          = 2;
    amgparam->amli_coef            = NULL;
    amgparam->nl_amli_krylov_type  = SOLVER_VFGMRES;
    
    // Classical AMG specific
    amgparam->coarsening_type      = COARSE_C;
    amgparam->interpolation_type   = INTERP_STD;
    amgparam->max_row_sum          = 0.9;
    amgparam->strong_threshold     = 0.3;
    amgparam->truncation_threshold = 0.2;
    amgparam->aggressive_level     = 0;
    amgparam->aggressive_path      = 1;
    
    // Aggregation AMG specific
    amgparam->aggregation_type     = VMB;
    amgparam->quality_bound        = 8.0;
    amgparam->pair_number          = 2;
    amgparam->strong_coupled       = 0.00;
    amgparam->max_aggregation      = 20;
    
    // ILU smoother parameters
    amgparam->ILU_type             = ILUt;
    amgparam->ILU_levels           = 0;
    amgparam->ILU_lfil             = 0;
    amgparam->ILU_droptol          = 0.01;
    amgparam->ILU_relax            = 0;
    
    // Schwarz smoother parameters
    amgparam->Schwarz_levels       = 0; // how many levels will use Schwarz smoother
    amgparam->Schwarz_mmsize       = 200;
    amgparam->Schwarz_maxlvl       = 3; // block size -- all vertices at distance .le. this
    amgparam->Schwarz_type         = 1;
    amgparam->Schwarz_blksolver    = SOLVER_DEFAULT;
}

/*************************************************************************************/
void param_ilu_init (ILU_param *iluparam)
{
    /*!
     * \fn void param_ilu_init (ILU_param *iluparam)
     *
     * \brief Initialize ILU parameters
     *
     * \param iluparam  Pointer to the LIU_param structure
     *
     */
    
    iluparam->print_level  = PRINT_NONE;
    iluparam->ILU_type     = ILUt;
    iluparam->ILU_lfil     = 0;
    iluparam->ILU_droptol  = 0.01;
    iluparam->ILU_relax    = 0;
    iluparam->ILU_permtol  = 0.01;

}

/*************************************************************************************/
void param_linear_solver_init (linear_itsolver_param *itsparam)
{
    /*!
     * \fn void param_linear_solver_init (linear_itsolver_param *itsparam)
     *
     * \brief Initialize linear iterative solver parameters
     *
     * \param itsparam  Pointer to the linear_itsolver_param structure
     *
     */

    itsparam->linear_print_level   = 0;
    itsparam->linear_itsolver_type = SOLVER_VFGMRES;
    itsparam->linear_precond_type  = PREC_NULL;
    itsparam->linear_stop_type     = STOP_REL_RES;
    itsparam->linear_maxit         = 500;
    itsparam->linear_restart       = 100;
    itsparam->linear_tol           = 1e-6;
    
    // HX preconditioner
    itsparam->HX_smooth_iter       = 1;

}

/*************************************************************************************/
void param_ilu_set (ILU_param *iluparam,
                    input_param *inparam)
{
    /*!
     * \fn void param_ilu_set (ILU_param *iluparam, input_param *iniparam)
     *
     * \brief Set ILU_param using input_param
     *
     * \param iluparam     Pointer to ILU_param structure
     * \param iniparam     Pointer to input_param structure
     *
     */
    
    iluparam->print_level = inparam->print_level;
    iluparam->ILU_type    = inparam->ILU_type;
    iluparam->ILU_lfil    = inparam->ILU_lfil;
    iluparam->ILU_droptol = inparam->ILU_droptol;
    iluparam->ILU_relax   = inparam->ILU_relax;
    iluparam->ILU_permtol = inparam->ILU_permtol;

}

/*************************************************************************************/
void param_linear_solver_set (linear_itsolver_param *itsparam,
                       input_param *inparam)
{
    /*!
     * \fn void param_linear_solver_set (linear_itsolver_param *itsparam, input_param *iniparam)
     *
     * \brief Set linear_itsolver_param using input_param
     *
     * \param itsparam   Pointer to linear_itsovler_param structure
     * \param inparam    Pointer to input_param structure
     *
     */
    
    itsparam->linear_print_level    = inparam->print_level;
    itsparam->linear_itsolver_type  = inparam->linear_itsolver_type;
    itsparam->linear_stop_type      = inparam->linear_stop_type;
    itsparam->linear_restart        = inparam->linear_restart;
    itsparam->linear_precond_type   = inparam->linear_precond_type;
    
    if ( itsparam->linear_itsolver_type == SOLVER_AMG ) {
        itsparam->linear_tol   = inparam->AMG_tol;
        itsparam->linear_maxit = inparam->AMG_maxit;
    }
    else {
        itsparam->linear_tol   = inparam->linear_itsolver_tol;
        itsparam->linear_maxit = inparam->linear_itsolver_maxit;
    }
    
    itsparam->HX_smooth_iter = inparam->HX_smooth_iter;
    
}

/*************************************************************************************/
void param_amg_set (AMG_param *amgparam,
                    input_param *inparam)
{
    /*!
     * \fn void param_amg_set (AMG_param *param, input_param *inparam)
     *
     * \brief Set AMG_param using input_param
     *
     * \param amgparam   Pointer to the AMG_param structure
     * \param inparam    Pointer to the input_param structure
     *
     */
    
    amgparam->AMG_type    = inparam->AMG_type;
    amgparam->print_level = inparam->print_level;
    
    if (inparam->linear_itsolver_type == SOLVER_AMG) {
        amgparam->maxit = inparam->linear_itsolver_maxit;
        amgparam->tol   = inparam->linear_itsolver_tol;
    }
    else {
        amgparam->maxit = inparam->AMG_maxit;
        amgparam->tol   = inparam->AMG_tol;
    }
    
    amgparam->max_levels           = inparam->AMG_levels;
    amgparam->cycle_type           = inparam->AMG_cycle_type;
    amgparam->smoother             = inparam->AMG_smoother;
    amgparam->smooth_order         = inparam->AMG_smooth_order;
    amgparam->relaxation           = inparam->AMG_relaxation;
    amgparam->coarse_solver        = inparam->AMG_coarse_solver;
    amgparam->polynomial_degree    = inparam->AMG_polynomial_degree;
    amgparam->presmooth_iter       = inparam->AMG_presmooth_iter;
    amgparam->postsmooth_iter      = inparam->AMG_postsmooth_iter;
    amgparam->coarse_solver        = inparam->AMG_coarse_solver;
    amgparam->coarse_dof           = inparam->AMG_coarse_dof;
    amgparam->coarse_scaling       = inparam->AMG_coarse_scaling;
    amgparam->amli_degree          = inparam->AMG_amli_degree;
    amgparam->amli_coef            = NULL;
    amgparam->nl_amli_krylov_type  = inparam->AMG_nl_amli_krylov_type;
    
    amgparam->coarsening_type      = inparam->AMG_coarsening_type;
    amgparam->interpolation_type   = inparam->AMG_interpolation_type;
    amgparam->strong_threshold     = inparam->AMG_strong_threshold;
    amgparam->truncation_threshold = inparam->AMG_truncation_threshold;
    amgparam->max_row_sum          = inparam->AMG_max_row_sum;
    amgparam->aggressive_level     = inparam->AMG_aggressive_level;
    amgparam->aggressive_path      = inparam->AMG_aggressive_path;
    
    amgparam->aggregation_type     = inparam->AMG_aggregation_type;
    amgparam->pair_number          = inparam->AMG_pair_number;
    amgparam->quality_bound        = inparam->AMG_quality_bound;
    amgparam->strong_coupled       = inparam->AMG_strong_coupled;
    amgparam->max_aggregation      = inparam->AMG_max_aggregation;

    amgparam->ILU_levels           = inparam->AMG_ILU_levels;
    amgparam->ILU_type             = inparam->ILU_type;
    amgparam->ILU_lfil             = inparam->ILU_lfil;
    amgparam->ILU_droptol          = inparam->ILU_droptol;
    amgparam->ILU_relax            = inparam->ILU_relax;
    amgparam->ILU_permtol          = inparam->ILU_permtol;

    amgparam->Schwarz_levels       = inparam->AMG_Schwarz_levels;
    amgparam->Schwarz_mmsize       = inparam->Schwarz_mmsize;
    amgparam->Schwarz_maxlvl       = inparam->Schwarz_maxlvl;
    amgparam->Schwarz_type         = inparam->Schwarz_type;

}

/*************************************************************************************/
void param_linear_solver_print (linear_itsolver_param *itsparam)
{
    /*!
     * \fn void param_linear_solver_print (linear_itsolver_param *itsparam)
     *
     * \brief Print out linear iterative solver parameters
     *
     * \param param  Pointer to the lienar_itsolver_parame structure
     *
     */
    
    if ( itsparam ) {
        
        printf("\n     Parameters in linear_itsolver_param     \n");
        printf("-----------------------------------------------\n");
        
        printf("Solver print level:                %d\n", itsparam->linear_print_level);
        printf("Solver type:                       %d\n", itsparam->linear_itsolver_type);
        printf("Solver precond type:               %d\n", itsparam->linear_precond_type);
        printf("Solver max num of iter:            %d\n", itsparam->linear_maxit);
        printf("Solver tolerance:                  %.2e\n", itsparam->linear_tol);
        printf("Solver stopping type:              %d\n", itsparam->linear_stop_type);
        printf("Solver restart number:             %d\n", itsparam->linear_restart);
        
        if ( (itsparam->linear_precond_type == PREC_HX_CURL_A) || (itsparam->linear_precond_type == PREC_HX_CURL_M) )
            printf("HX precond number of smooth:       %d\n", itsparam->HX_smooth_iter);
    
        printf("-----------------------------------------------\n\n");
        
    }
    else {
        printf("### WARNING HAZMAT DANGER: linear solver parameters have not been set!!! \n");
        printf("Set linear solver parameters to default !!!\n");
        param_linear_solver_init(itsparam);
    }
}

/*************************************************************************************/
void param_ilu_print (ILU_param *iluparam)
{
    /*!
     * \fn void param_ilu_print (ILU_param *iluparam)
     *
     * \brief Print out ILU solver or smoother parameters
     *
     * \param param    Pointer to the ILU_param structure
     *
     */
    
    if ( iluparam ) {
        
        printf("\n       Parameters in ILU_param\n");
        printf("-----------------------------------------------\n");
        printf("ILU print level:                   %d\n",   iluparam->print_level);
        printf("ILU type:                          %d\n",   iluparam->ILU_type);
        printf("ILU level of fill-in:              %d\n",   iluparam->ILU_lfil);
        printf("ILU relaxation factor:             %.4f\n", iluparam->ILU_relax);
        printf("ILU drop tolerance:                %.2e\n", iluparam->ILU_droptol);
        printf("ILU permutation tolerance:         %.2e\n", iluparam->ILU_permtol);
        printf("-----------------------------------------------\n\n");
        
    }
    else {
        printf("### WARNING HAZMAT DANGER: ILU parameters have not been set!\n");
        printf("Set ILU parameters to default !!!\n");
        param_ilu_init(iluparam);
    }
}

/*************************************************************************************/
void param_amg_print (AMG_param *amgparam)
{
    
    /*!
     * \fn void param_amg_print (AMG_param *amgparam)
     *
     * \brief Print out AMG parameters
     *
     * \param param   Pointer to AMG_param structure
     *
     */
    
    if ( amgparam ) {
        
        printf("\n       Parameters in AMG_param\n");
        printf("-----------------------------------------------\n");
        
        printf("AMG print level:                   %d\n", amgparam->print_level);
        printf("AMG max num of iter:               %d\n", amgparam->maxit);
        printf("AMG type:                          %d\n", amgparam->AMG_type);
        printf("AMG tolerance:                     %.2e\n", amgparam->tol);
        printf("AMG max levels:                    %d\n", amgparam->max_levels);
        printf("AMG cycle type:                    %d\n", amgparam->cycle_type);
        printf("AMG coarse dof:                    %d\n", amgparam->coarse_dof);
        printf("AMG coarse solver type:            %d\n", amgparam->coarse_solver);
        printf("AMG scaling of coarse correction:  %d\n", amgparam->coarse_scaling);
        printf("AMG smoother type:                 %d\n", amgparam->smoother);
        printf("AMG smoother order:                %d\n", amgparam->smooth_order);
        printf("AMG num of presmoothing:           %d\n", amgparam->presmooth_iter);
        printf("AMG num of postsmoothing:          %d\n", amgparam->postsmooth_iter);
        
        if ( amgparam->smoother == SMOOTHER_SOR  ||
             amgparam->smoother == SMOOTHER_SSOR ||
             amgparam->smoother == SMOOTHER_GSOR ||
             amgparam->smoother == SMOOTHER_SGSOR ) {
            printf("AMG relax factor:                  %.4f\n", amgparam->relaxation);
        }
        
        /*
        if ( param->smoother == SMOOTHER_POLY ) {
            printf("AMG polynomial smoother degree:    %d\n", param->polynomial_degree);
        }
         */
        
        if ( amgparam->cycle_type == AMLI_CYCLE ) {
            printf("AMG AMLI degree of polynomial:     %d\n", amgparam->amli_degree);
        }
        
        if ( amgparam->cycle_type == NL_AMLI_CYCLE ) {
            printf("AMG Nonlinear AMLI Krylov type:    %d\n", amgparam->nl_amli_krylov_type);
        }
        
        switch (amgparam->AMG_type) {
            case CLASSIC_AMG:
                printf("AMG coarsening type:               %d\n", amgparam->coarsening_type);
                printf("AMG interpolation type:            %d\n", amgparam->interpolation_type);
                printf("AMG dof on coarsest grid:          %d\n", amgparam->coarse_dof);
                printf("AMG strong threshold:              %.4f\n", amgparam->strong_threshold);
                printf("AMG truncation threshold:          %.4f\n", amgparam->truncation_threshold);
                printf("AMG max row sum:                   %.4f\n", amgparam->max_row_sum);
                printf("AMG aggressive levels:             %d\n", amgparam->aggressive_level);
                printf("AMG aggressive path:               %d\n", amgparam->aggressive_path);
                break;
                
            default: // UA_AMG
                printf("Aggregation type:                  %d\n", amgparam->aggregation_type);
                if ( amgparam->aggregation_type == PAIRWISE ) {
                    printf("Aggregation number of pairs:       %d\n", amgparam->pair_number);
                    printf("Aggregation quality bound:         %.2f\n", amgparam->quality_bound);
                }
                if ( amgparam->aggregation_type == VMB ) {
                    printf("Aggregation AMG strong coupling:   %.4f\n", amgparam->strong_coupled);
                    printf("Aggregation AMG max aggregation:   %d\n", amgparam->max_aggregation);
                }
                break;
        }
        
        if ( amgparam->ILU_levels>0 ) {
            printf("AMG ILU smoother level:            %d\n", amgparam->ILU_levels);
            printf("AMG ILU type:                      %d\n", amgparam->ILU_type);
            printf("AMG ILU level of fill-in:          %d\n", amgparam->ILU_lfil);
            printf("AMG ILU drop tol:                  %e\n", amgparam->ILU_droptol);
            printf("AMG ILU relaxation:                %f\n", amgparam->ILU_relax);
        }
        
        if ( amgparam->Schwarz_levels>0 ){
            printf("AMG Schwarz smoother level:        %d\n", amgparam->Schwarz_levels);
            printf("AMG Schwarz type:                  %d\n", amgparam->Schwarz_type);
            printf("AMG Schwarz forming block level:   %d\n", amgparam->Schwarz_maxlvl);
            printf("AMG Schwarz maximal block size:    %d\n", amgparam->Schwarz_mmsize);
        }
        
        printf("-----------------------------------------------\n\n");
        
    }
    else {
        printf("### WARNING HAZMAT DANGER: AMG parameters have not been set!\n");
        printf("Set AMG parameters to default !!!\n");
        param_amg_init(amgparam);
    } // end if (param)
    
}

/*************************************************************************************/
void param_amg_to_prec (precond_data *pcdata,
                        AMG_param *amgparam)
{
    /*!
     * \fn void param_amg_to_prec (precond_data *pcdata, AMG_param *amgparam)
     *
     * \brief Set parameters in precond_data using AMG parameters
     *
     * \param pcdata      Pointer to the precond_data structure
     * \param amgparam    Pointer to the AMG_param structure
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

/*************************************************************************************/
void param_prec_to_amg (AMG_param *amgparam,
                        precond_data *pcdata)
{
    /*!
     * \fn void fasp_param_prec_to_amg (AMG_param *amgparam, precond_data *pcdata)
     *
     * \brief Set AMG parameters using parameters in precond_data
     *
     * \param amgparam    Pointer to the AMG_param structuree
     * \param pcdata      Pointer to the precond_data structure
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

/*************************************************************************************/
void amg_amli_coef (const REAL lambda_max,
                    const REAL lambda_min,
                    const INT degree,
                    REAL *coef)
{
    
    /*!
     * \fn void amg_amli_coef (const REAL lambda_max, const REAL lambda_min,
     *                              const INT degree, REAL *coef)
     *
     * \brief Compute the coefficients of the polynomial used by AMLI-cycle
     *
     * \param lambda_max  Maximal eigenvalue (estimated)
     * \param lambda_min  Minimal eigenvalue (estimated)
     * \param degree      Degree of polynomial used in AMLI-cycle
     * \param coef        Pointer to the coefficients of AMLI-cycle (output)
     *
     * \note Best polynomial approximation of 1/x is used here -- Xiaozhe Hu
     * \todo Other polynomials -- Xiaozhe Hu
     *
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
        printf("### ERROR HAZMAT DANGER: Wrong AMLI degree %d!\n", degree);
        check_error(ERROR_INPUT_PAR, __FUNCTION__);
    }
    
    return;
}

/******************************** END ************************************************/
