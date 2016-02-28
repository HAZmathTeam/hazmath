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
void ilu_data_null (ILU_data *ILUdata)
{
    /**
     * \fn void ilu_data_null (ILU_data *ILUdata)
     *
     * \brief Initialize ILU data
     *
     * \param ILUdata   Pointer to ILU_data
     *
     * \author Chensong Zhang
     * \date   2010/03/23
     */
    
    ILUdata->row = ILUdata->col = ILUdata->nzlu = 0;
    ILUdata->ijlu = NULL; ILUdata->luval = NULL;
}

/***********************************************************************************************/
void ilu_data_alloc (const INT iwk,
                          const INT nwork,
                          ILU_data *iludata)
{
    /**
     * \fn void ilu_data_alloc (const INT iwk, const INT nwork, ILU_data *iludata)
     *
     * \brief Allocate workspace for ILU factorization
     *
     * \param iwk       Size of the index array
     * \param nwork     Size of the work array
     * \param iludata   Pointer to the ILU_data
     *
     * \author Chensong Zhang
     * \date   2010/04/06
     */

    iludata->ijlu=(INT*)calloc(iwk, sizeof(INT));
    
    iludata->luval=(REAL*)calloc(iwk, sizeof(REAL));
    
    iludata->work=(REAL*)calloc(nwork, sizeof(REAL));
    
    return;
}

/***********************************************************************************************/
void ilu_data_free (ILU_data *ILUdata)
{
    /**
     * \fn void ilu_data_free (ILU_data *ILUdata)
     *
     * \brief Create ILU_data sturcture
     *
     * \param ILUdata   Pointer to ILU_data
     *
     * \author Chensong Zhang
     * \date   2010/04/03
     */
    
    if (ILUdata==NULL) return;
    
    free(ILUdata->ijlu);  ILUdata->ijlu  = NULL;
    free(ILUdata->luval); ILUdata->luval = NULL;
    free(ILUdata->work);  ILUdata->work  = NULL;
    
    ILUdata->row = ILUdata->col = ILUdata->nzlu = ILUdata ->nwork = ILUdata->nb = 0;
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
        ilu_data_free(&mgl[i].LU);
    }
    
    for (i=0; i<mgl->near_kernel_dim; ++i) {
        if (mgl->near_kernel_basis[i]) free(mgl->near_kernel_basis[i]);
        mgl->near_kernel_basis[i] = NULL;
    }
    
    // Clean direct solver data in necessary
    switch (param->coarse_solver) {
            
#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
            free(mgl[max_levels-1].Numeric);
            break;
        }
#endif
            
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
void HX_curl_data_free (HX_curl_data *hxcurldata, SHORT flag)
{
    /**
     * \fn void HX_curl_data_free (HX_curl_data *hxcurldata)
     *
     * \brief Free HX_curl_data data memeory space
     *
     * \param hxcurldata    Pointer to the HX_curl_data
     * \param flag          A flag: 0 means A, P_curl, and Grad will be reused | 1 means free everything
     *
     * \author Xiaozhe Hu
     * \date   02/11/2016
     *
     *
     */
    
    if (flag == TRUE) {
        dcsr_free(hxcurldata->A);
        dcsr_free(hxcurldata->P_curl);
        dcsr_free(hxcurldata->Grad);
    }
    
    dcsr_free(hxcurldata->Pt_curl);
    dcsr_free(hxcurldata->A_vgrad);
    
    amg_data_free(hxcurldata->mgl_vgrad, hxcurldata->amgparam_vgrad);
    
    dcsr_free(hxcurldata->Gradt);
    dcsr_free(hxcurldata->A_grad);
    
    amg_data_free(hxcurldata->mgl_grad, hxcurldata->amgparam_grad);
    
    if (hxcurldata->backup_r) free(hxcurldata->backup_r);
    
    if (hxcurldata->w) free(hxcurldata->w);
        

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
