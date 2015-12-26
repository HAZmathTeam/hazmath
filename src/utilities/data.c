/*
 *  data.c
 *
 *  Created by James Adler and Xiaozhe Hu on 12/23/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

#include "hazmat.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/***********************************************************************************************/
void precond_data_null (precond_data *pcdata)
{
    /**
     * \fn void precond_data_null (precond_data *pcdata)
     *
     * \brief Initialize precond_data
     *
     * \param pcdata   Preconditioning data structure
     *
     * \author Chensong Zhang & Xiaozhe Hu
     * \date   12/23/2015
     */
    
    pcdata->AMG_type            = UA_AMG;
    pcdata->print_level         = PRINT_NONE;
    pcdata->maxit               = 500;
    pcdata->max_levels          = 20;
    pcdata->tol                 = 1e-8;
    pcdata->cycle_type          = V_CYCLE;
    pcdata->smoother            = SMOOTHER_GS;
    pcdata->smooth_order        = CF_ORDER;
    pcdata->presmooth_iter      = 1;
    pcdata->postsmooth_iter     = 1;
    pcdata->relaxation          = 1.1;
    pcdata->coarsening_type     = 1;
    pcdata->coarse_scaling      = ON;
    pcdata->amli_degree         = 1;
    pcdata->nl_amli_krylov_type = SOLVER_GCG;
}

/***********************************************************************************************/
AMG_data * amg_data_create (SHORT max_levels)
{
    /**
     * \fn AMG_data * amg_data_create (SHORT max_levels)
     *
     * \brief Create and initialize AMG_data
     *
     * \param max_levels   Max number of levels allowed
     *
     * \return Pointer to the AMG_data data structure
     *
     * \author Chensong Zhang & Xiaozhe
     * \date   12/24/2015
     */
    
    max_levels = MAX(1, max_levels); // at least allocate one level
    
    AMG_data *mgl = (AMG_data *)calloc(max_levels,sizeof(AMG_data));
    
    INT i;
    for ( i=0; i<max_levels; ++i ) {
        mgl[i].max_levels = max_levels;
        mgl[i].num_levels = 0;
        mgl[i].near_kernel_dim = 0;
        mgl[i].near_kernel_basis = NULL;
        mgl[i].cycle_type = 0;
    }
    
    return(mgl);
}

/***********************************************************************************************/
void amg_data_free (AMG_data *mgl,
                         AMG_param *param)
{
    /**
     * \fn void amg_data_free (AMG_data *mgl, AMG_param *param)
     *
     * \brief Free AMG_data data memeory space
     *
     * \param mgl    Pointer to the AMG_data
     * \param param  Pointer to AMG parameters
     *
     * \author Chensong Zhang & Xiaozhe Hu
     * \date   12/24/2015
     *
     */
    
    const INT max_levels = MAX(1,mgl[0].num_levels);
    
    INT i;
    
    for (i=0; i<max_levels; ++i) {
        dcsr_free(&mgl[i].A);
        dcsr_free(&mgl[i].P);
        dcsr_free(&mgl[i].R);
        dvec_free(&mgl[i].b);
        dvec_free(&mgl[i].x);
        dvec_free(&mgl[i].w);
        ivec_free(&mgl[i].cfmark);
    }
    
    for (i=0; i<mgl->near_kernel_dim; ++i) {
        if (mgl->near_kernel_basis[i]) free(mgl->near_kernel_basis[i]);
        mgl->near_kernel_basis[i] = NULL;
    }
    
    // Clean direct solver data in necessary
    switch (param->coarse_solver) {
            
        default: // Do nothing!
            break;
    }
    
    free(mgl->near_kernel_basis); mgl->near_kernel_basis = NULL;
    
    free(mgl); mgl = NULL;
    
    if (param != NULL) {
        if ( param->cycle_type == AMLI_CYCLE ) free(param->amli_coef);
    }
}

/***********************************************************************************************/
void precond_null (precond *pcdata)
{
    /**
     * \fn void precond_null (precond *pcdata)
     *
     * \brief Initialize precond data
     *
     * \param pcdata   Pointer to precond
     *
     * \author Chensong Zhang & Xiaozhe Hu
     * \date   12/23/2015
     */
    
    pcdata->data = NULL;
    pcdata->fct  = NULL;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
