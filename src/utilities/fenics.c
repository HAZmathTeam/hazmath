/*! \file src/utilities/fenics.c
 *
 *  Created by Ana Budisa on 2020-11-25.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note: some useful functions so python/fenics can stop bothering me with
 *         lost pointers and seg faults
 *  \note: make sure all functions here return either void, int, or pointers
 *         (helps with the implementation in python's ctypes)
 *
 */

#include "hazmath.h"

/***********************************************************************************************/
/*!
 * \fn INT fenics_amg_data_setup (dCSRmat *A, AMG_data *mgl, AMG_param *amgparam)
 *
 * \brief Setup AMG_data structure and build levels
 *
 * \param A         Pointer to dCSRmat on coarsest level
 * \param mgl       Pointer to AMG_data
 * \param amgparam  Pointer to AMG_param
 *
 * \return          INT (1 if successful setup, 0 else)
 *
 * \author          Ana Budisa
 * \date            2020-11-25
 */
INT fenics_amg_data_setup(dCSRmat *A,
                          AMG_data *mgl,
                          AMG_param *amgparam)
{
    const SHORT prtlvl = amgparam->print_level;
    const SHORT max_levels = amgparam->max_levels;
    const INT m = A->row, n = A->col, nnz = A->nnz;
    INT status = SUCCESS;

    // Create enough mg levels
    // mgl = amg_data_create(max_levels);
    // initialize coarsest level
    mgl[0].A = dcsr_create(m, n, nnz); dcsr_cp(A, &mgl[0].A);
    mgl[0].b = dvec_create(m);
    mgl[0].x = dvec_create(n);

    // Build levels
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl, amgparam);
        break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = amg_setup_sa(mgl, amgparam);
        break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl, amgparam);
        break;

    }

    return((status < 0) ? 0 : 1);

}

/***********************************************************************************************/
/*!
 * \fn INT fenics_famg_data_setup (dCSRmat *A, dCSRmat *M, AMG_data *mgl, AMG_param *amgparam)
 *
 * \brief Setup AMG_data structure and build levels
 *
 * \param A         Pointer to dCSRmat stiff matrix on coarsest level
 * \param M         Pointer to dCSRmat mass matrix on coarses level
 * \param mgl       Pointer to AMG_data
 * \param amgparam  Pointer to AMG_param
 *
 * \return          INT (1 if successful setup, 0 else)
 *
 * \author          Ana Budisa
 * \date            2020-11-27
 */
INT fenics_famg_data_setup(dCSRmat *A,
                           dCSRmat *M,
                           AMG_data *mgl,
                           AMG_param *amgparam)
{
    const SHORT prtlvl = amgparam->print_level;
    const SHORT max_levels = amgparam->max_levels;
    const INT m = A->row, n = A->col, nnz = A->nnz, nnz_M = M->nnz;
    INT status = SUCCESS;

    // Create enough mg levels
    // mgl = amg_data_create(max_levels);
    // initialize coarsest level
    mgl[0].A = dcsr_create(m, n, nnz); dcsr_cp(A, &mgl[0].A);
    mgl[0].M = dcsr_create(m, n, nnz_M); dcsr_cp(M, &mgl[0].M);
    mgl[0].b = dvec_create(m);
    mgl[0].x = dvec_create(n);

    // Build levels
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = famg_setup_ua(mgl, amgparam);
        break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = famg_setup_sa(mgl, amgparam);
        break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = famg_setup_ua(mgl, amgparam);
        break;

    }

    return((status < 0) ? 0 : 1);

}


/***********************************************************************************************/
/*!
 * \fn void fenics_precond_data_setup (AMG_data *mgl, AMG_param *amgparam, precond_data *pcdata)
 *
 * \brief Setup precond_data structure from AMG_data and AMG_param
 *
 * \param mgl       Pointer to AMG_data
 * \param amgparam  Pointer to AMG_param
 * \param pcdata    Pointer to precond_data
 *
 * \return          INT (1 if successful setup, 0 else)
 *
 * \author          Ana Budisa
 * \date            2020-11-25
 */
void fenics_precond_data_setup(AMG_data *mgl,
                              AMG_param *amgparam,
                              precond_data *pcdata)
{
    // Initialize precond_data
    // precond_data_null(pcdata);

    // Setup parameters
    param_amg_to_prec(pcdata, amgparam);

    // Setup AMG levels
    pcdata->max_levels = mgl[0].num_levels;
    pcdata->mgl_data = mgl;

}

/*!
 * \fn precond_data *precond_data_alloc (SHORT max_size)
 *
 * \brief Allocate precond_data array of length max_size; each component
 *        is initialized
 *
 * \param max_size  Size of precond_data array (usually 1)
 *
 * \return pcdata   precond_data* (callocated)
 *
 * \author          Ana Budisa
 * \date            2020-11-26
 */
precond_data *precond_data_alloc(SHORT max_size)
{
    max_size = MAX(1, max_size);

    precond_data *pcdata = (precond_data*)calloc(max_size, sizeof(precond_data));

    INT i;
    for(i = 0; i < max_size; ++i) precond_data_null(&(pcdata[i]));

    return(pcdata);
}

/*!
 * \fn AMG_param *amg_param_alloc (SHORT max_size)
 *
 * \brief Allocate AMG_param array of length max_size; each component
 *        is initialized
 *
 * \param max_size      Size of AMG_param array (usually 1)
 *
 * \return amgparam     AMG_param* (callocated)
 *
 * \author          Ana Budisa
 * \date            2020-11-26
 */
AMG_param *amg_param_alloc(SHORT max_size)
{
    max_size = MAX(1, max_size);

    AMG_param *amgparam = (AMG_param*)calloc(max_size, sizeof(AMG_param));

    INT i;
    for(i = 0; i < max_size; ++i)  param_amg_init(&(amgparam[i]));

    return(amgparam);
}


void save_poles_residues(INT k, REAL *poles, REAL *residues, precond_data *pcdata)
{
    pcdata->r = (dvector*)calloc(1, sizeof(dvector));
    pcdata->r->row = 2 * k + 1;
    pcdata->r->val = (REAL *)calloc(2 * k + 1, sizeof(REAL));
    array_cp(k, poles, pcdata->r->val);
    array_cp(k + 1, residues, &(pcdata->r->val[k]));
}