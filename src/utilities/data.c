/*! \file src/utilities/data.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 12/23/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note modified by Xiaozhe Hu 10/27/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/27/2016 & 08/28/2021
 *
 */

#include "hazmath.h"

/***********************************************************************************************/
/*!
 * \fn void precond_data_null (precond_data *pcdata)
 *
 * \brief Initialize precond_data (pointers are set to NULL) (OUTPUT)
 *
 * \param pcdata   Preconditioning data structure
 *
 */
void precond_data_null (precond_data *pcdata)
{
    pcdata->AMG_type            = UA_AMG;
    pcdata->print_level         = PRINT_MIN;
    pcdata->maxit               = 100;
    pcdata->max_levels          = 20;
    pcdata->tol                 = 1e-8;
    pcdata->cycle_type          = V_CYCLE;
    pcdata->smoother            = SMOOTHER_GS;
    pcdata->presmooth_iter      = 1;
    pcdata->postsmooth_iter     = 1;
    pcdata->relaxation          = 1.2;
    pcdata->polynomial_degree   = 2;
    pcdata->coarse_solver       = SOLVER_UMFPACK;
    pcdata->coarse_scaling      = OFF;
    pcdata->amli_degree         = 2;
    pcdata->nl_amli_krylov_type = SOLVER_VFGMRES;
    pcdata->fpwr                = 1.0;

    pcdata->amli_coef           = NULL;
    pcdata->mgl_data            = NULL;
    pcdata->A                   = NULL;

    pcdata->A_nk                = NULL;
    pcdata->P_nk                = NULL;
    pcdata->R_nk                = NULL;

    pcdata->r                   = NULL;
    pcdata->w                   = NULL;

}


/***********************************************************************************************/
/*!
 * \fn AMG_data * amg_data_create (SHORT max_levels)
 *
 * \brief Create AMG_data structure (but all values are 0 and pointers point to NULL)
 *
 * \param max_levels   Max number of levels allowed
 *
 * \return Pointer to the AMG_data structure
 *
 */
AMG_data *amg_data_create(SHORT max_levels)
{
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
/*!
 * \fn void amg_data_free(AMG_data *mgl, AMG_param *param)
 *
 * \brief Free AMG_data structure
 *
 * \param mgl    Pointer to the AMG_data (OUTPUT)
 * \param param  Pointer to AMG parameters
 *
 *
 */
void amg_data_free(AMG_data *mgl,
                   AMG_param *param)
{
    const INT max_levels = MAX(1,mgl[0].num_levels);

    INT i;

    for (i=0; i<max_levels; ++i) {
        dcsr_free(&mgl[i].A);
        dcsr_free(&mgl[i].P);
        dcsr_free(&mgl[i].R);
        dcsr_free(&mgl[i].M);
        dvec_free(&mgl[i].b);
        dvec_free(&mgl[i].x);
        dvec_free(&mgl[i].w);

        // free Schwarz data
        if ( i < param->Schwarz_levels ) {
            schwarz_data_free(&mgl[i].Schwarz);
        }
    }

    for (i=0; i<mgl->near_kernel_dim; ++i) {
        if (mgl->near_kernel_basis[i]) free(mgl->near_kernel_basis[i]);
        mgl->near_kernel_basis[i] = NULL;
    }

    // Clean direct solver data if necessary
    switch (param->coarse_solver) {

      //#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
	  hazmath_free_numeric(&(mgl[max_levels-1].Numeric));
	  break;
        }
	  //#endif

        default: // Do nothing!
            break;
    }

    free(mgl->near_kernel_basis);
    mgl->near_kernel_basis = NULL;

    if (param != NULL) {
        if ( param->cycle_type == AMLI_CYCLE )
            free(param->amli_coef);
    }

    free(mgl);


}

/**
 * \fn AMG_data_bsr * amg_data_bsr_create (SHORT max_levels)
 *
 * \brief Create and initialize AMG_data data sturcture for AMG/SAMG (BSR format)
 *
 * \param max_levels   Max number of levels allowed
 *
 * \return Pointer to the AMG_data data structure
 *
 * \author Xiaozhe Hu
 * \date   08/07/2011
 */
AMG_data_bsr * amg_data_bsr_create (SHORT max_levels)
{
    max_levels = MAX(1, max_levels); // at least allocate one level

    AMG_data_bsr *mgl = (AMG_data_bsr *)calloc(max_levels,sizeof(AMG_data_bsr));

    INT i;
    for (i=0; i<max_levels; ++i) {
        mgl[i].max_levels = max_levels;
        mgl[i].num_levels = 0;
        mgl[i].near_kernel_dim = 0;
        mgl[i].near_kernel_basis = NULL;
        mgl[i].A_nk = NULL;
        mgl[i].P_nk = NULL;
        mgl[i].R_nk = NULL;
    }

    return(mgl);
}

/**
 * \fn void amg_data_bsr_free (AMG_data_bsr *mgl)
 *
 * \brief Free AMG_data_bsr data memeory space
 *
 * \param mgl  Pointer to the AMG_data_bsr
 *
 * \author Xiaozhe Hu
 * \date   2013/02/13
 *
 */
void amg_data_bsr_free (AMG_data_bsr *mgl)
{
    const INT max_levels = MAX(1,mgl[0].num_levels);

    INT i;

    for ( i = 0; i < max_levels; ++i ) {

        dbsr_free(&mgl[i].A);
        if ( max_levels > 1 ) {
            dbsr_free(&mgl[i].P);
            dbsr_free(&mgl[i].R);
        }
        dvec_free(&mgl[i].b);
        dvec_free(&mgl[i].x);
        dvec_free(&mgl[i].diaginv);
        dvec_free(&mgl[i].diaginv_SS);
        dcsr_free(&mgl[i].Ac);

        dcsr_free(&mgl[i].PP);
        dbsr_free(&mgl[i].SS);
        dvec_free(&mgl[i].diaginv_SS);
        dvec_free(&mgl[i].w);
        ivec_free(&mgl[i].cfmark);

        free(mgl[i].pw); mgl[i].pw = NULL;
        free(mgl[i].sw); mgl[i].sw = NULL;
    }

    for ( i = 0; i < mgl->near_kernel_dim; ++i ) {
        free(mgl->near_kernel_basis[i]); mgl->near_kernel_basis[i] = NULL;
    }
    free(mgl->near_kernel_basis); mgl->near_kernel_basis = NULL;
    free(mgl); mgl = NULL;
}

/***********************************************************************************************/
/*!
 * \fn AMG_data * amg_data_bdcsr_create (SHORT max_levels)
 *
 * \brief Create AMG_data_bdcsr structure (but all values are 0 and pointers point to NULL)
 *
 * \param max_levels   Max number of levels allowed
 *
 * \return Pointer to the AMG_data structure
 *
 */
AMG_data_bdcsr *amg_data_bdcsr_create(SHORT max_levels)
{
    max_levels = MAX(1, max_levels); // at least allocate one level

    AMG_data_bdcsr *mgl = (AMG_data_bdcsr *)calloc(max_levels,sizeof(AMG_data_bdcsr));

    INT i;
    for ( i=0; i<max_levels; ++i ) {
        mgl[i].max_levels = max_levels;
        mgl[i].num_levels = 0;
        mgl[i].near_kernel_dim = 0;
        mgl[i].near_kernel_basis = NULL;
        mgl[i].cycle_type = 0;
    }

    mgl[0].A_gamma = NULL;

    return(mgl);
}

/**
 * \fn void amg_data_bdcsr_free (AMG_data_bsr *mgl, AMG_param *param)
 *
 * \brief Free AMG_data_bdcsr data memeory space
 *
 * \param mgl       Pointer to the AMG_data_bdcsr
 * \param param     Pointer to the AMG parameters
 *
 * \author Xiaozhe Hu
 * \date   03/16/2022
 *
 */
void amg_data_bdcsr_free (AMG_data_bdcsr *mgl,
                          AMG_param *param)
{
    const INT max_levels = MAX(1,mgl[0].num_levels);

    INT i, j;
    INT brow = mgl[0].A.brow;

    //printf("in free\n");

    for ( i = 0; i < max_levels; ++i ) {

        //printf("i = %d\n", i);

        bdcsr_free(&mgl[i].A);
        //printf("done free A\n");
        if ( max_levels > 1 ) {
            bdcsr_free(&mgl[i].P);
            bdcsr_free(&mgl[i].R);
        }
        //printf("done free P and R\n");
        dvec_free(&mgl[i].b);
        //printf("done free b\n");
        dvec_free(&mgl[i].x);
        //printf("done free x\n");
        dcsr_free(&mgl[i].Ac);
        //printf("done free Ac\n");
        dvec_free(&mgl[i].w);
        //printf("done free w\n");

        for (j=0; j<brow; j++) dcsr_free(&mgl[i].A_diag[j]);
        if (mgl[i].A_diag) free(mgl[i].A_diag);

    }


    // clean AMG data for diaognal blocks
    for (i=0; i<brow; i++){
        amg_data_free(mgl[0].mgl_diag[i], param);
        if (mgl[0].mgl_diag[i]) free(mgl[0].mgl_diag[i]);
    }
    if (mgl[0].mgl_diag) free(mgl[0].mgl_diag);

    // Clean direct solver data if necessary
    switch (param->coarse_solver) {

      //#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
	  hazmath_free_numeric(&(mgl[max_levels-1].Numeric));
            break;
        }
	  //#endif

        default: // Do nothing!
            break;
    }

    //printf("free kernel\n");

    for ( i = 0; i < mgl->near_kernel_dim; ++i ) {
        free(mgl->near_kernel_basis[i]); mgl->near_kernel_basis[i] = NULL;
    }
    free(mgl->near_kernel_basis); mgl->near_kernel_basis = NULL;

    //printf("done free kernel, start free interface \n");

    bdcsr_free(mgl[0].A_gamma);
    dbsr_free(&mgl[0].A_gamma_bsr);
    dvec_free(&mgl[0].A_gamma_diaginv);

    //printf("done free interface\n");

    free(mgl); mgl = NULL;
}


/***********************************************************************************************/
/*!
 * \fn void HX_curl_data_null(HX_curl_data *hxcurldata)
 *
 * \brief Initalize HX_curl_data structure (set values to 0 and pointers to NULL) (OUTPUT)
 *
 * \param hxcurldata    Pointer to the HX_curl_data structure
 *
 */
void HX_curl_data_null (HX_curl_data *hxcurldata)
{
    hxcurldata->A               = NULL;

    hxcurldata->smooth_type     = 0;
    hxcurldata->smooth_iter     = 0;

    hxcurldata->P_curl          = NULL;
    hxcurldata->Pt_curl         = NULL;
    hxcurldata->A_vgrad         = NULL;

    hxcurldata->amgparam_vgrad  = NULL;
    hxcurldata->mgl_vgrad       = NULL;

    hxcurldata->Grad            = NULL;
    hxcurldata->Gradt           = NULL;
    hxcurldata->A_grad          = NULL;

    hxcurldata->amgparam_grad   = NULL;
    hxcurldata->mgl_grad        = NULL;

    hxcurldata->backup_r        = NULL;
    hxcurldata->w               = NULL;

}

/***********************************************************************************************/
/*!
 * \fn void HX_curl_data_free (HX_curl_data *hxcurldata, SHORT flag)
 *
 * \brief Free HX_curl_data structure (set values to 0 and pointers to NULL)
 *
 * \param hxcurldata    Pointer to the HX_curl_data structure (OUTPUT)
 * \param flag          flag of whether the date will be reused:
 *                      flag = False - A, P_curl, and Grad will be reused
 *                      flag = TRUE  - free everything
 *
 */
void HX_curl_data_free (HX_curl_data *hxcurldata,
                        SHORT flag)
{
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
/*!
 * \fn void HX_div_data_null(HX_div_data *hxdivdata)
 *
 * \brief Initalize HX_div_data structure (set values to 0 and pointers to NULL) (OUTPUT)
 *
 * \param hxdivdata    Pointer to the HX_div_data structure
 *
 */
void HX_div_data_null (HX_div_data *hxdivdata)
{
    hxdivdata->A                = NULL;

    hxdivdata->smooth_type      = 0;
    hxdivdata->smooth_iter      = 0;

    hxdivdata->P_curl           = NULL;
    hxdivdata->Pt_curl          = NULL;
    hxdivdata->P_div            = NULL;
    hxdivdata->Pt_div           = NULL;
    hxdivdata->A_curlgrad       = NULL;
    hxdivdata->A_divgrad        = NULL;

    hxdivdata->amgparam_curlgrad  = NULL;
    hxdivdata->mgl_curlgrad       = NULL;
    hxdivdata->amgparam_divgrad  = NULL;
    hxdivdata->mgl_divgrad       = NULL;

    hxdivdata->Curl            = NULL;
    hxdivdata->Curlt           = NULL;
    hxdivdata->A_grad          = NULL;
    hxdivdata->A_curl          = NULL;

    hxdivdata->amgparam_grad   = NULL;
    hxdivdata->mgl_grad        = NULL;

    hxdivdata->backup_r        = NULL;
    hxdivdata->w               = NULL;

}

/***********************************************************************************************/
/*!
 * \fn void HX_div_data_free (HX_div_data *hxdivdata, SHORT flag)
 *
 * \brief Free HX_div_data structure (set values to 0 and pointers to NULL)
 *
 * \param hxcurldata    Pointer to the HX_curl_data structure (OUTPUT)
 * \param flag          flag of whether the date will be reused:
 *                      flag = False - A, P_curl, and Grad will be reused
 *                      flag = TRUE  - free everything
 *
 */
void HX_div_data_free (HX_div_data *hxdivdata,
                        SHORT flag)
{
    if (flag == TRUE) {
        dcsr_free(hxdivdata->A);
        dcsr_free(hxdivdata->P_curl);
        dcsr_free(hxdivdata->P_div);
        dcsr_free(hxdivdata->Curl);
    }

    dcsr_free(hxdivdata->Pt_curl);
    dcsr_free(hxdivdata->Pt_div);
    dcsr_free(hxdivdata->Curlt);
    dcsr_free(hxdivdata->A_curl);
    dcsr_free(hxdivdata->A_curlgrad);
    dcsr_free(hxdivdata->A_divgrad);

    if (hxdivdata->mgl_curlgrad) amg_data_free(hxdivdata->mgl_curlgrad, hxdivdata->amgparam_curlgrad);
    if (hxdivdata->mgl_divgrad) amg_data_free(hxdivdata->mgl_divgrad, hxdivdata->amgparam_divgrad);

    if (hxdivdata->backup_r) free(hxdivdata->backup_r);
    if (hxdivdata->w) free(hxdivdata->w);

}
/***********************************************************************/
/**
 * \fn void schwarz_data_init (Schwarz_data *schwarzdata)
 * \brief initialize Schwarz data memory space
 *
 * \param schwarzdata      Pointer to the Schwarz_data for Schwarz methods
 *
 * \author Ludmil
 * \date   20221213
 */
void schwarz_data_init(Schwarz_data *Schwarz)
{
  Schwarz->A=dcsr_create(0,0,0); // NULL dcsr matrix.
  Schwarz->nblk=0;
  Schwarz->iblock=0;
  Schwarz->jblock=0;
  Schwarz->Schwarz_type=3; //this should be set.
  Schwarz->blk_solver=SOLVER_UMFPACK; //direct solve of all blocks.
  Schwarz->memt=0;
  Schwarz->mask=NULL;
  Schwarz->maxbs=0;
  Schwarz->blk_data=NULL;
  Schwarz->rhsloc1=dvec_create(0);
  Schwarz->xloc1=dvec_create(0);
  Schwarz->numeric=NULL;
  Schwarz->swzparam=NULL;
}
/***********************************************************************/
/**
 * \fn void schwarz_data_free (Schwarz_data *schwarzdata)
 * \brief Free Schwarz data memeory space
 *
 * \param swzdata      Pointer to the Schwarz_data for Schwarz methods
 *
 * \author Xiaozhe Hu
 * \date   2010/04/06
 * \modified 20221213 --ltz
 */
void schwarz_data_free(Schwarz_data *schwarzdata)
{
  INT i;
  if ( schwarzdata == NULL ) return; // There is nothing to do!
  //
  if(&(schwarzdata->A)!=NULL)
    dcsr_free(&schwarzdata->A);
  if(schwarzdata->Schwarz_type==SCHWARZ_FORWARD || schwarzdata->Schwarz_type==SCHWARZ_BACKWARD || schwarzdata->Schwarz_type==SCHWARZ_SYMMETRIC){
    for ( i=0; i<schwarzdata->nblk; ++i ) {
      dcsr_free(&((schwarzdata->blk_data)[i]));
      if (schwarzdata->blk_solver == SOLVER_UMFPACK){
	//#if WITH_SUITESPARSE
	if (schwarzdata->numeric[i])
	  hazmath_free_numeric(&(schwarzdata->numeric[i]));
	//#endif
      }
    }
  } else {
    // only one matrix then:
    dcsr_free(&((schwarzdata->blk_data)[0]));
    if (schwarzdata->blk_solver == SOLVER_UMFPACK){
      //#if WITH_SUITESPARSE
      if (schwarzdata->numeric[0])
	  hazmath_free_numeric(&(schwarzdata->numeric[0]));
      //#endif
    }
  }
  schwarzdata->nblk = 0;
  if (schwarzdata->blk_data) free(schwarzdata->blk_data);
  schwarzdata->blk_data = NULL;
  //
  if (schwarzdata->blk_solver == SOLVER_UMFPACK){
    //#if WITH_SUITESPARSE
    if (schwarzdata->numeric) free(schwarzdata->numeric);
    schwarzdata->numeric = NULL;
    //#endif
  }
  if (schwarzdata->iblock) free(schwarzdata->iblock);
  schwarzdata->iblock = NULL;
  //
  if (schwarzdata->jblock) free(schwarzdata->jblock);
  schwarzdata->jblock = NULL;
  //
  dvec_free(&schwarzdata->rhsloc1);
  //
  dvec_free(&schwarzdata->xloc1);
  //
  schwarzdata->memt = 0;
  if (schwarzdata->mask) free(schwarzdata->mask);
  schwarzdata->mask = NULL;
  //
  /* if (schwarzdata->maxa) free(schwarzdata->maxa); */
  /* schwarzdata->maxa = NULL; */
  return;
}

/***********************************************************************************************/
/**
 * \fn void precond_null(precond *pcdata)
 *
 * \brief Initialize precond data (set pointers to NULL)
 *
 * \param pcdata   Pointer to precond
 *
 */
void precond_null(precond *pcdata)
{
    pcdata->data = NULL;
    pcdata->fct  = NULL;
}

/***********************************************************************************************/
/**
 * \fn void precond_block_data_null(precond_block_data *precdata)
 *
 * \brief Initialize precond block data (set pointers to NULL)
 *
 * \param precdata   Pointer to precond block data
 *
 */
void precond_block_data_null(precond_block_data *precdata)
{

    precdata->Abcsr = NULL;

    precdata->A_diag = NULL;
    precdata->diag = NULL;

#if WITH_SUITESPARSE
    precdata->LU_diag = NULL;
#endif

    precdata->mgl = NULL;
    precdata->amgparam = NULL;

    precdata->hxcurldata = NULL;
    precdata->hxdivdata = NULL;

    precdata->el_vol = NULL;

    precdata->G = NULL;
    precdata->K = NULL;
    precdata->Gt = NULL;
    precdata->Kt = NULL;

    precdata->scaled_M = NULL;
    precdata->diag_scaled_M = NULL;
    precdata->poles = NULL;
    precdata->residues = NULL;

}

/***********************************************************************************************/
/*!
 * \fn void precond_block_data_free(precond_block_data *precdata,SHORT flag)
 *
 * \brief Free precond_block_data structure (set values to 0 and pointers to NULL)
 *
 * \param precdata      Pointer to the precond_block_data structure (OUTPUT)
 * \param nb            number of blocks
 *
 */
void precond_block_data_free(precond_block_data *precdata,
                             const INT nb,
                             SHORT flag)
{
    INT i;

    for (i=0; i<nb; i++)
    {

        if(precdata->diag) {
           if(precdata->diag[i]) dvec_free(precdata->diag[i]);
        }

        if(precdata->mgl) {
            if(precdata->mgl[i])
            {
              amg_data_free(precdata->mgl[i], &precdata->amgparam[i]);
              free(precdata->mgl[i]);
            }
        }

        if(precdata->hxcurldata) {
            if(precdata->hxcurldata[i])
            {
              HX_curl_data_free(precdata->hxcurldata[i],flag);
              free(precdata->hxcurldata[i]);
            }
        }

        if(precdata->hxdivdata) {
            if(precdata->hxdivdata[i])
            {
              HX_div_data_free(precdata->hxdivdata[i],flag);
              free(precdata->hxdivdata[i]);
            }
        }

    }

    if(precdata->diag) free(precdata->diag);
    if(precdata->mgl) free(precdata->mgl);
    if(precdata->hxcurldata) free(precdata->hxcurldata);
    if(precdata->hxdivdata)  free(precdata->hxdivdata);

    //#if WITH_SUITESPARSE
    for (i=0; i<nb; i++)
    {
        if(precdata->LU_diag){
	  if(precdata->LU_diag[i])
	    hazmath_free_numeric(&(precdata->LU_diag[i]));
        }
    }
    if(precdata->LU_diag) free(precdata->LU_diag);
    //#endif

    dvec_free(&precdata->r);

    if (flag == TRUE)
    {
      if(precdata->G)  dcsr_free(precdata->G);
      if(precdata->K)  dcsr_free(precdata->K);
      if(precdata->Gt) dcsr_free(precdata->Gt);
      if(precdata->Kt) dcsr_free(precdata->Kt);
    }

    if(precdata->scaled_M) free(precdata->scaled_M);
    if(precdata->diag_scaled_M) free(precdata->diag_scaled_M);
    if(precdata->poles) free(precdata->poles);
    if(precdata->residues) free(precdata->residues);

    return;
}

/***********************************************************************************************/
/*!
 * \fn void amli_coef_free(AMG_param *param)
 *
 * \brief Free amli coefficients in AMG_param structure
 *
 * \param param  Pointer to AMG parameters
 *
 *
 */
void amli_coef_free(AMG_param *param)
{
    if (param != NULL) {
        if ( param->cycle_type == AMLI_CYCLE )
            free(param->amli_coef);
    }

}

/***********************************************************************************************/
/*!
 * \fn void precond_ra_data_free(precond_ra_data *precdata)
 *
 * \brief Free precond_ra_data structure
 *
 * \param precdata      Pointer to the precond_ra_data structure
 *
 */
void precond_ra_data_free(precond_ra_data *precdata)
{

    INT np = precdata->poles->row;
    INT i;

    for (i = 0; i < np; i++)
    {
        if(precdata->mgl) {
            if(precdata->mgl[i])
            {
              amg_data_free(precdata->mgl[i], NULL);
              free(precdata->mgl[i]);
            }
        }
    }
    if(precdata->mgl) free(precdata->mgl);

    if(precdata->amgparam) amli_coef_free(precdata->amgparam);

#if WITH_SUITESPARSE
    for (i = 0; i < np; i++)
    {
        if(precdata->LU_diag){
	  if(precdata->LU_diag[i])
	    hazmath_free_numeric(&(precdata->LU_diag[i]));
        }
    }
    if(precdata->LU_diag) free(precdata->LU_diag);
#endif

    if(precdata->scaled_M)  dcsr_free(precdata->scaled_M);
    if(precdata->scaled_A)  dcsr_free(precdata->scaled_A);

    if(precdata->diag_scaled_M) dvec_free(precdata->diag_scaled_M);
    if(precdata->poles) dvec_free(precdata->poles);
    if(precdata->residues) dvec_free(precdata->residues);
    if(precdata->r) dvec_free(precdata->r);

    if(precdata->w) free(precdata->w);

    return;
}

/***********************************************************************************************/
/*!
 * \fn void precond_data_free(precond_ra_data *precdata)
 *
 * \brief Free precond_data structure
 *
 * \param precdata      Pointer to the precond_data structure
 *
 */
void precond_data_free(precond_data *precdata)
{
    if(precdata->amli_coef) free(precdata->amli_coef);

    if(precdata->mgl_data) {
        amg_data_free(precdata->mgl_data, NULL);
        free(precdata->mgl_data);
    }

    if(precdata->A)  dcsr_free(precdata->A);
    if(precdata->A_nk)  dcsr_free(precdata->A_nk);
    if(precdata->P_nk)  dcsr_free(precdata->P_nk);
    if(precdata->R_nk)  dcsr_free(precdata->R_nk);

    if(precdata->r) dvec_free(precdata->r);

    if(precdata->w) free(precdata->w);

    return;
}


/***********************************************************************************************/
/*!
 * \fn void precond_data_bdcsr_free(precond_data_bdcsr *precdata)
 *
 * \brief Free precond_data_bdcsr structure
 *
 * \param precdata      Pointer to the precond_data_bdcsr structure
 *
 */
void precond_data_bdcsr_free(precond_data_bdcsr *precdata)
{
    if(precdata->amli_coef) free(precdata->amli_coef);

    if(precdata->mgl_data) {
        amg_data_bdcsr_free(precdata->mgl_data, NULL);
        free(precdata->mgl_data);
    }

    if(precdata->schwarz_data) schwarz_data_free(precdata->schwarz_data);
    if(precdata->LU_data){
        if(precdata->LU_data[0]) hazmath_free_numeric(&precdata->LU_data[0]);
    }
    if(precdata->LU_data) free(precdata->LU_data);

    if(precdata->A)  bdcsr_free(precdata->A);

    if(&(precdata->r)) dvec_free(&(precdata->r));
    if(precdata->w) free(precdata->w);
    if(&(precdata->perm)) ivec_free(&(precdata->perm));

    return;
}

/*************************************  END  ***************************************************/
