/*! \file src/utilities/parameters.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 10/06/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/27/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/28/2016 & 08/28/2021
 *
 */

#include "hazmath.h"

/*************************************************************************************/
/*!
 * \fn input_param *param_input_initP()
 *
 * \brief Initialize input parameters
 *
 * \return inparam    Pointer to input_param structure
 *
 */
input_param *param_input_init_p()
{
  input_param *inparam=malloc(1*sizeof(input_param));
  param_input_init (inparam);
  return inparam;
}

/*************************************************************************************/
/*!
 * \fn void param_input_init (input_param *inparam)
 *
 * \brief Initialize input parameters
 *
 * \param inparam    Pointer to input_param structure
 *
 * \note added frac. exponent (Ana Budisa, 2020-05-13)
 */
void param_input_init (input_param *inparam)
{
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
    // mesh parameters
    //--------------------------
    inparam->read_mesh_from_file=0;
    inparam->spatial_dim=2;
    inparam->refinement_type=11;
    inparam->refinement_levels=2;
    inparam->boundary_codes=1;

    //--------------------------
    // finite element parameters
    //--------------------------
    // general parameter
    inparam->nquad                    = 2;

    // parameters for H(D) equations
    inparam->FE_type                  = 1;
    inparam->Mass_lump                  = 0;

    //----------------------------
    // time steppng paramters
    //----------------------------
    inparam->time_start           = 0.0;
    inparam->time_step_type           = 0;
    inparam->time_steps               = 0;
    inparam->time_step_size           = 0.01;
    inparam->rhs_time_dep           = 1;

    //----------------------------
    // nonlinear solver parameters
    //----------------------------
    inparam->nonlinear_itsolver_type    = 0;
    inparam->nonlinear_itsolver_maxit   = 0;
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

    // AMG method parameters
    inparam->AMG_type                 = UA_AMG;
    inparam->AMG_levels               = 10;
    inparam->AMG_cycle_type           = W_CYCLE;
    inparam->AMG_smoother             = SMOOTHER_GS;
    inparam->AMG_presmooth_iter       = 1;
    inparam->AMG_postsmooth_iter      = 1;
    inparam->AMG_polynomial_degree    = 2;
    inparam->AMG_relaxation           = 1.2;
    inparam->AMG_coarse_dof           = 200;
    inparam->AMG_coarse_solver        = SOLVER_DEFAULT;
    inparam->AMG_tol                  = 1e-6;
    inparam->AMG_maxit                = 1;
    inparam->AMG_Schwarz_levels       = 0;
    inparam->AMG_coarse_scaling       = OFF;
    inparam->AMG_amli_degree          = 1;
    inparam->AMG_nl_amli_krylov_type  = 2;
    inparam->AMG_fpwr                 = 1.0;

    // Aggregation AMG parameters
    inparam->AMG_aggregation_type     = HEC;
    inparam->AMG_strong_coupled       = 0.04;
    inparam->AMG_max_aggregation      = 20;

    inparam->AMG_tentative_smooth     = 0.67;
    inparam->AMG_smooth_filter        = ON;

    // Schwarz method parameters
    inparam->Schwarz_mmsize           = 200;
    inparam->Schwarz_maxlvl           = 2;
    inparam->Schwarz_type             = 1;
    inparam->Schwarz_blksolver        = SOLVER_DEFAULT;

    // HX Preconditioner
    inparam->HX_smooth_iter           = 1;

    // BSR Preconditioner
    inparam->BSR_alpha           = 1.;
    inparam->BSR_omega           = 1.;
    return;
}

/*************************************************************************************/
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
 * \note added frac. exponent (Ana Budisa, 2020-05-13)
 * \note added function pointer to user def smoother (Ana Budisa, 2021-04-27)
 */
void param_amg_init (AMG_param *amgparam)
{
    // General AMG parameters
    amgparam->AMG_type             = UA_AMG;
    amgparam->print_level          = PRINT_NONE;
    amgparam->maxit                = 1;
    amgparam->tol                  = 1e-6;
    amgparam->max_levels           = 10;
    amgparam->coarse_dof           = 200;
    amgparam->cycle_type           = W_CYCLE;
    amgparam->smoother             = SMOOTHER_GS;
    amgparam->presmooth_iter       = 1;
    amgparam->postsmooth_iter      = 1;
    amgparam->coarse_solver        = SOLVER_DEFAULT;
    amgparam->relaxation           = 1.0;
    amgparam->polynomial_degree    = 2;
    amgparam->coarse_scaling       = OFF;
    amgparam->amli_degree          = 2;
    amgparam->amli_coef            = NULL;
    amgparam->nl_amli_krylov_type  = SOLVER_VFGMRES;
    amgparam->fpwr                 = 1.0;

    // Aggregation AMG parameters
    amgparam->aggregation_type     = HEC;
    amgparam->strong_coupled       = 0.04;
    amgparam->max_aggregation      = 20;

    amgparam->tentative_smooth     = 0.67;
    amgparam->smooth_filter        = ON;

    // Schwarz smoother parameters
    amgparam->Schwarz_levels       = 1; // how many levels will use Schwarz smoother
    amgparam->Schwarz_mmsize       = 200;
    amgparam->Schwarz_maxlvl       = -1; // blocksize -- vertices with smaller distance
    amgparam->Schwarz_type         = 1;
    amgparam->Schwarz_blksolver    = SOLVER_UMFPACK;

    // Other smoother param
    amgparam->HAZDIR     = NULL;
    amgparam->Schwarz_on_blk     = NULL;
    amgparam->Schwarz_patch_type = NULL;
    amgparam->damping_param        = 1.0;
    amgparam->BSR_alpha            = -1000.;
    amgparam->BSR_omega            = -1000.;

    // user def smoother
    // amgparam->smoother_function = NULL;
}

/*************************************************************************************/
/*!
 * \fn void param_Schwarz_init (Schwarz_param *schparam)
 *
 * \brief Initialize Schwarz parameters
 *
 * \param schparam    Parameters for Schwarz method
 *
 * \author Xiaozhe Hu
 * \date   05/22/2012
 *
 */
void param_Schwarz_init (Schwarz_param *schparam)
{
    schparam->print_level       = PRINT_NONE;
    schparam->Schwarz_type      = 3;
    schparam->Schwarz_maxlvl    = 1;
    schparam->Schwarz_mmsize    = 200;
    schparam->Schwarz_blksolver = 0;
}


/*************************************************************************************/
/*!
 * \fn void param_linear_solver_init (linear_itsolver_param *itsparam)
 *
 * \brief Initialize linear iterative solver parameters
 *
 * \param itsparam  Pointer to the linear_itsolver_param structure
 *
 */
void param_linear_solver_init (linear_itsolver_param *itsparam)
{
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
/*!
 * \fn void param_linear_solver_set (linear_itsolver_param *itsparam, input_param *iniparam)
 *
 * \brief Set linear_itsolver_param using input_param
 *
 * \param itsparam   Pointer to linear_itsovler_param structure
 * \param inparam    Pointer to input_param structure
 *
 */
void param_linear_solver_set (linear_itsolver_param *itsparam,
                       input_param *inparam)
{
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
/*!
 * \fn void param_amg_set (AMG_param *param, input_param *inparam)
 *
 * \brief Set AMG_param using input_param
 *
 * \param amgparam   Pointer to the AMG_param structure
 * \param inparam    Pointer to the input_param structure
 *
 * \note added frac. exponent (Ana Budisa, 2020-05-13)
 */
void param_amg_set (AMG_param *amgparam,
                    input_param *inparam)
{
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
    amgparam->fpwr                 = inparam->AMG_fpwr;

    amgparam->aggregation_type     = inparam->AMG_aggregation_type;
    amgparam->strong_coupled       = inparam->AMG_strong_coupled;
    amgparam->max_aggregation      = inparam->AMG_max_aggregation;

    amgparam->tentative_smooth     = inparam->AMG_tentative_smooth;
    amgparam->smooth_filter        = inparam->AMG_smooth_filter;

    amgparam->Schwarz_levels       = inparam->AMG_Schwarz_levels;
    amgparam->Schwarz_mmsize       = inparam->Schwarz_mmsize;
    amgparam->Schwarz_maxlvl       = inparam->Schwarz_maxlvl;
    amgparam->Schwarz_type         = inparam->Schwarz_type;
    amgparam->Schwarz_blksolver    = inparam->Schwarz_blksolver;
    amgparam->Schwarz_on_blk     = NULL;
    amgparam->HAZDIR     = NULL;
    amgparam->Schwarz_on_blk     = NULL;
    amgparam->Schwarz_patch_type = NULL;
    amgparam->damping_param=0.;
    amgparam->BSR_alpha            = inparam->BSR_alpha;
    amgparam->BSR_omega            = inparam->BSR_omega;

}

/*************************************************************************************/
/*!
 * \fn void param_amg_cp (AMG_param *amgparam1, AMG_param *amgparam2)
 *
 * \brief Copy AMG_param amgparam1 to amgparam2
 *
 * \param amgparam1   Pointer to the AMG_param structure
 * \param amgparam2   Pointer to the AMG_param structure
 *
 * \note maybe the order of pointers should be switched, but I don't know what is
 *       the convention in hazmath  -- Ana
 * \note the order is okay -- Xiaozhe
 */
void param_amg_cp (AMG_param *amgparam1,
                   AMG_param *amgparam2)
{
    amgparam2->AMG_type    = amgparam1->AMG_type;
    amgparam2->print_level = amgparam1->print_level;
    amgparam2->maxit = amgparam1->maxit;
    amgparam2->tol   = amgparam1->tol;

    amgparam2->max_levels           = amgparam1->max_levels;
    amgparam2->cycle_type           = amgparam1->cycle_type;
    amgparam2->smoother             = amgparam1->smoother;
    amgparam2->relaxation           = amgparam1->relaxation;
    amgparam2->coarse_solver        = amgparam1->coarse_solver;
    amgparam2->polynomial_degree    = amgparam1->polynomial_degree;
    amgparam2->presmooth_iter       = amgparam1->presmooth_iter;
    amgparam2->postsmooth_iter      = amgparam1->postsmooth_iter;
    amgparam2->coarse_solver        = amgparam1->coarse_solver;
    amgparam2->coarse_dof           = amgparam1->coarse_dof;
    amgparam2->coarse_scaling       = amgparam1->coarse_scaling;
    amgparam2->amli_degree          = amgparam1->amli_degree;

    if(amgparam1->amli_coef) array_cp(amgparam1->amli_degree + 1, amgparam1->amli_coef, amgparam2->amli_coef);

    amgparam2->nl_amli_krylov_type  = amgparam1->nl_amli_krylov_type;
    amgparam2->fpwr                 = amgparam1->fpwr;

    amgparam2->aggregation_type     = amgparam1->aggregation_type;
    amgparam2->strong_coupled       = amgparam1->strong_coupled;
    amgparam2->max_aggregation      = amgparam1->max_aggregation;

    amgparam2->tentative_smooth     = amgparam1->tentative_smooth;
    amgparam2->smooth_filter        = amgparam1->smooth_filter;

    amgparam2->Schwarz_levels       = amgparam1->Schwarz_levels;
    amgparam2->Schwarz_mmsize       = amgparam1->Schwarz_mmsize;
    amgparam2->Schwarz_maxlvl       = amgparam1->Schwarz_maxlvl;
    amgparam2->Schwarz_type         = amgparam1->Schwarz_type;
    amgparam2->Schwarz_blksolver    = amgparam1->Schwarz_blksolver;

    //if(amgparam1->Schwarz_on_blk) iarray_cp(size, amgparam1->Schwarz_on_blk, amgparam2->Schwarz_on_blk);
    //if(amgparam1->Schwarz_patch_type) iarray_cp(size, amgparam1->Schwarz_patch_type, amgparam2->Schwarz_patch_type);
    //if(amgparam1->HAZDIR){
    //    amgparam2->HAZDIR = malloc(strlen(amgparam1->HAZDIR) + 1);
    //    strcpy(amgparam2->HAZDIR, amgparam1->HAZDIR);
    // }
    amgparam2->Schwarz_on_blk = NULL;
    amgparam2->Schwarz_patch_type = NULL;
    amgparam2->HAZDIR = NULL;

    amgparam2->damping_param        = amgparam1->damping_param;
    amgparam2->BSR_alpha            = amgparam1->BSR_alpha;
    amgparam2->BSR_omega            = amgparam1->BSR_omega;

}

/*************************************************************************************/
/**
 * \fn void param_Schwarz_set (Schwarz_param *schparam, input_param *iniparam)
 *
 * \brief Set Schwarz_param with INPUT
 *
 * \param schparam    Parameters for Schwarz method
 * \param iniparam     Input parameters
 *
 * \author Xiaozhe Hu
 * \date   05/22/2012
 */
void param_Schwarz_set (Schwarz_param *schparam,
                             input_param *inparam)
{
    schparam->print_level       = inparam->print_level;
    schparam->Schwarz_type      = inparam->Schwarz_type;
    schparam->Schwarz_maxlvl    = inparam->Schwarz_maxlvl;
    schparam->Schwarz_mmsize    = inparam->Schwarz_mmsize;
    schparam->Schwarz_blksolver = inparam->Schwarz_blksolver;
}

/*************************************************************************************/
/*!
 * \fn void param_linear_solver_print (linear_itsolver_param *itsparam)
 *
 * \brief Print out linear iterative solver parameters
 *
 * \param param  Pointer to the lienar_itsolver_parame structure
 *
 */
void param_linear_solver_print (linear_itsolver_param *itsparam)
{
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
        printf("### WARNING HAZMATH DANGER: linear solver parameters have not been set!!! \n");
        printf("Set linear solver parameters to default !!!\n");
        param_linear_solver_init(itsparam);
    }
}

/*************************************************************************************/
/*!
 * \fn void param_amg_print (AMG_param *amgparam)
 *
 * \brief Print out AMG parameters
 *
 * \param param   Pointer to AMG_param structure
 *
 * \note added frac. exponent (Ana Budisa, 2020-05-13)
 */
void param_amg_print (AMG_param *amgparam)
{
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
        printf("AMG num of presmoothing:           %d\n", amgparam->presmooth_iter);
        printf("AMG num of postsmoothing:          %d\n", amgparam->postsmooth_iter);

        if ( amgparam->smoother == SMOOTHER_SOR  ||
             amgparam->smoother == SMOOTHER_SSOR ||
             amgparam->smoother == SMOOTHER_GSOR ||
             amgparam->smoother == SMOOTHER_SGSOR ) {
            printf("AMG relax factor:                  %.4f\n", amgparam->relaxation);
        }

        if ( amgparam->smoother == SMOOTHER_FJACOBI ||
             amgparam->smoother == SMOOTHER_FGS     ||
             amgparam->smoother == SMOOTHER_FSGS    ) {
            printf("AMG fractional exponent:           %.4f\n", amgparam->fpwr);
        }

        if ( amgparam->cycle_type == AMLI_CYCLE ) {
            printf("AMG AMLI degree of polynomial:     %d\n", amgparam->amli_degree);
        }

        if ( amgparam->cycle_type == NL_AMLI_CYCLE ) {
            printf("AMG Nonlinear AMLI Krylov type:    %d\n", amgparam->nl_amli_krylov_type);
        }

        switch (amgparam->AMG_type) {

            case SA_AMG:
                printf("Aggregation type:                  %d\n", amgparam->aggregation_type);
                printf("Aggregation AMG strong coupling:   %.4f\n", amgparam->strong_coupled);
                printf("Aggregation AMG max aggregation:   %d\n", amgparam->max_aggregation);
                printf("SA AMG tentative smooth parameter: %.4f\n", amgparam->tentative_smooth);
                printf("SA AMG smooth filter:              %d\n", amgparam->smooth_filter);

            default: // UA_AMG
                printf("Aggregation type:                  %d\n", amgparam->aggregation_type);
                printf("Aggregation AMG strong coupling:   %.4f\n", amgparam->strong_coupled);
                printf("Aggregation AMG max aggregation:   %d\n", amgparam->max_aggregation);
                break;
        }

        if (amgparam->Schwarz_levels>0){
            printf("AMG Schwarz smoother level:        %d\n", amgparam->Schwarz_levels);
            printf("AMG Schwarz type:                  %d\n", amgparam->Schwarz_type);
            printf("AMG Schwarz forming block level:   %d\n", amgparam->Schwarz_maxlvl);
            printf("AMG Schwarz maximal block size:    %d\n", amgparam->Schwarz_mmsize);
        }

        printf("-----------------------------------------------\n\n");

    }
    else {
        printf("### WARNING HAZMATH DANGER: AMG parameters have not been set!\n");
        printf("Set AMG parameters to default !!!\n");
        param_amg_init(amgparam);
    } // end if (param)

}

/*************************************************************************************/
/**
 * \fn void param_Schwarz_print (Schwarz_param *param)
 *
 * \brief Print out Schwarz parameters
 *
 * \param param    Parameters for Schwarz
 *
 * \author Xiaozhe Hu
 * \date   05/22/2012
 */
void param_Schwarz_print (Schwarz_param *schparam)
{
    if ( schparam ) {

        printf("\n       Parameters in Schwarz_param\n");
        printf("-----------------------------------------------\n");
        printf("Schwarz print level:               %d\n",   schparam->print_level);
        printf("Schwarz type:                      %d\n",   schparam->Schwarz_type);
        printf("Schwarz forming block level:       %d\n",   schparam->Schwarz_maxlvl);
        printf("Schwarz maximal block size:        %d\n",   schparam->Schwarz_mmsize);
        printf("Schwarz block solver type:         %d\n",   schparam->Schwarz_blksolver);
        printf("-----------------------------------------------\n\n");

    }
    else {
        printf("### WARNING: param has not been set!\n");
    }
}

/*************************************************************************************/
/*!
 * \fn void param_amg_to_prec (precond_data *pcdata, AMG_param *amgparam)
 *
 * \brief Set parameters in precond_data using AMG parameters
 *
 * \param pcdata      Pointer to the precond_data structure
 * \param amgparam    Pointer to the AMG_param structure
 *
 * \note added frac. exponent (Ana Budisa, 2020-05-13)
 * \note added function pointer to user defined smoother (Ana Budisa, 2021-04-27)
 */
void param_amg_to_prec (precond_data *pcdata,
                        AMG_param *amgparam)
{
    pcdata->AMG_type            = amgparam->AMG_type;
    pcdata->print_level         = amgparam->print_level;
    pcdata->maxit               = amgparam->maxit;
    pcdata->max_levels          = amgparam->max_levels;
    pcdata->tol                 = amgparam->tol;
    pcdata->cycle_type          = amgparam->cycle_type;
    pcdata->smoother            = amgparam->smoother;
    pcdata->presmooth_iter      = amgparam->presmooth_iter;
    pcdata->postsmooth_iter     = amgparam->postsmooth_iter;
    pcdata->coarse_solver       = amgparam->coarse_solver;
    pcdata->relaxation          = amgparam->relaxation;
    pcdata->polynomial_degree   = amgparam->polynomial_degree;
    pcdata->coarse_scaling      = amgparam->coarse_scaling;
    pcdata->amli_degree         = amgparam->amli_degree;
    pcdata->amli_coef           = amgparam->amli_coef;
    pcdata->nl_amli_krylov_type = amgparam->nl_amli_krylov_type;
    pcdata->fpwr                = amgparam->fpwr;
}

/*************************************************************************************/
/*!
 * \fn void param_prec_to_amg (AMG_param *amgparam, precond_data *pcdata)
 *
 * \brief Set AMG parameters using parameters in precond_data
 *
 * \param amgparam    Pointer to the AMG_param structuree
 * \param pcdata      Pointer to the precond_data structure
 *
 * \note added frac. exponent (Ana Budisa, 2020-05-13)
 * \note added function pointer to user defined smoother (Ana Budisa, 2021-04-27)
 */
void param_prec_to_amg (AMG_param *amgparam,
                        precond_data *pcdata)
{
    amgparam->AMG_type            = pcdata->AMG_type;
    amgparam->print_level         = pcdata->print_level;
    amgparam->cycle_type          = pcdata->cycle_type;
    amgparam->smoother            = pcdata->smoother;
    amgparam->presmooth_iter      = pcdata->presmooth_iter;
    amgparam->postsmooth_iter     = pcdata->postsmooth_iter;
    amgparam->relaxation          = pcdata->relaxation;
    amgparam->polynomial_degree   = pcdata->polynomial_degree;
    amgparam->coarse_solver       = pcdata->coarse_solver;
    amgparam->coarse_scaling      = pcdata->coarse_scaling;
    amgparam->amli_degree         = pcdata->amli_degree;
    amgparam->amli_coef           = pcdata->amli_coef;
    amgparam->nl_amli_krylov_type = pcdata->nl_amli_krylov_type;
    amgparam->fpwr                = pcdata->fpwr;
}

/*************************************************************************************/
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
void amg_amli_coef (const REAL lambda_max,
                    const REAL lambda_min,
                    const INT degree,
                    REAL *coef)
{
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
        printf("### ERROR HAZMATH DANGER: Wrong AMLI degree %d!\n", degree);
        check_error(ERROR_INPUT_PAR, __FUNCTION__);
    }

    return;
}

/******************************** END ************************************************/
