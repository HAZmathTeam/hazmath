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

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

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
/*--      Public Functions      --*/
/*---------------------------------*/

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

/*!
 * \fn void save_poles_residues(INT k, REAL *poles, REAL *residues, precond_data *pcdata)
 *
 * \brief Save poles & residues from rational approximation to precond_data
 *
 * \param k         Number of poles
 * \param poles     Array of poles
 * \param residues  Array of residues
 * \param pcdata    Preconditioner data
 *
 * \return
 *
 * \author          --
 * \date            --
 */
void save_poles_residues(INT k, REAL *poles, REAL *residues, precond_data *pcdata)
{
    pcdata->r = (dvector*)calloc(1, sizeof(dvector));
    pcdata->r->row = 2 * k + 1;
    pcdata->r->val = (REAL *)calloc(2 * k + 1, sizeof(REAL));
    array_cp(k, poles, pcdata->r->val);
    array_cp(k + 1, residues, &(pcdata->r->val[k]));
}

/*!
 * \fn precond_data *precond_block_data_alloc (SHORT max_size)
 *
 * \brief Allocate precond_block_data array of length max_size; each component
 *        is initialized
 *
 * \param max_size  Size of precond_block_data array (usually 1)
 *
 * \return pcdata   pointer to precond_block_data (callocated)
 *
 * \author          Ana Budisa
 * \date            2021-01-21
 */
precond_block_data *precond_block_data_alloc(SHORT max_size)
{
    max_size = MAX(1, max_size);

    precond_block_data *pcdata = (precond_block_data*)calloc(max_size, sizeof(precond_block_data));

    INT i;
    for(i = 0; i < max_size; ++i) precond_block_data_null(&(pcdata[i]));

    return(pcdata);
}

/***********************************************************************************************/
/*!
 * \fn void fenics_precond_block_data_setup (AMG_data **mgl, AMG_param *amgparam, precond_block_data *pcdata)
 *
 * \brief Setup precond_block_data structure from AMG_data and AMG_param
 *
 * \param mgl       Pointer to pointer to AMG_data
 * \param amgparam  Pointer to AMG_param
 * \param pcdata    Pointer to precond_block_data
 *
 * \return          --
 *
 * \author          Ana Budisa
 * \date            2021-01-21
 */
void fenics_precond_block_data_setup(block_dCSRmat *A,
                                     AMG_data **mgl,
                                     AMG_param *amgparam,
                                     dCSRmat *scaled_M,
                                     dvector *diag_scaled_M,
                                     REAL scaled_alpha,
                                     REAL scaled_beta,
                                     dvector *poles,
                                     dvector *residues,
                                     precond_block_data *pcdata)
{
    // Initialize precond_data
    // precond_data_null(pcdata);

    // Setup matrix data
    pcdata->Abcsr = (block_dCSRmat*)calloc(1, sizeof(block_dCSRmat));
    bdcsr_alloc(A->brow, A->bcol, pcdata->Abcsr);
    bdcsr_cp(A, pcdata->Abcsr);

    // Setup parameters
    pcdata->amgparam = amgparam;

    // Setup scaled matrices and parameters
    pcdata->scaled_M = scaled_M;
    pcdata->diag_scaled_M = diag_scaled_M;
    pcdata->scaled_alpha = scaled_alpha;
    pcdata->scaled_beta = scaled_beta;

    // Poles and residues from rational approximation
    pcdata->poles = poles;
    pcdata->residues = residues;

    // Setup AMG levels
    pcdata->mgl = mgl;

}

/***********************************************************************************************/
/*!
 * \fn INT fenics_ra_setup (dCSRmat *A, dCSRmat *M, REAL s_frac_power, REAL t_frac_power,
 *                          REAL alpha, REAL beta, REAL scaling_a, REAL scaling_m,
 *                          AMG_param *amgparam)
 *
 * \brief Setup AMG_data structure and build levels
 *
 * \param A             Pointer to dCSRmat stiff matrix
 * \param M             Pointer to dCSRmat mass matrix
 * \param s_frac_power  "First" fractional power
 * \param t_frac_power  "Second" fractional power
 * \param alpha         "First" weight
 * \param beta          "Second" weight
 * \param scaling a     Scaling on matrix A (usually ||A||_inf)
 * \param scaling m     Scaling on matrix M (usually ||M||_inf)
 * \param amgparam      Pointer to AMG_param
 * \param pcdata        Pointer to precond_block_data
 *
 * \return          INT (1 if successful setup, 0 else)
 *
 * \author          Ana Budisa
 * \date            2021-01-21
 */
INT fenics_ra_setup(dCSRmat *A,
                    dCSRmat *M,
                    REAL s_frac_power,
                    REAL t_frac_power,
                    REAL alpha,
                    REAL beta,
                    REAL scaling_a,
                    REAL scaling_m,
                    AMG_param *amgparam,
                    precond_block_data *pcdata)
{
    // AMG_param *amgparam1 = &(amgparam[1]);
    const SHORT prtlvl = amgparam->print_level;
    const SHORT max_levels = amgparam->max_levels;
    const INT m = A->row, n = A->col, nnz = A->nnz, nnz_M = M->nnz;
    INT status = SUCCESS;
    INT i;

    //------------------------------------------------
    // compute the rational approximation
    //------------------------------------------------
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

    /* AAA algorithm for the rational approximation */
    // parameters used in the AAA algorithm
    REAL xmin_in=0.e0, xmax_in=1.e0;  // interval for x
    INT mbig=(1<<14)+1;  // initial number of points on the interval [x_min, x_max]
    INT mmax_in=(INT )(mbig/2);  // maximal final number of pole + 1
    REAL16 AAA_tol=powl(2e0,-52e0);  // tolerance of the AAA algorithm
    INT k=-22; // k is the number of nodes in the final interpolation after tolerance is achieved or mmax is reached.
    INT print_level=0; // print level for AAA

    // output of the AAA algorithm.  It contains residues, poles, nodes, weights, function values
    REAL **rpnwf=malloc(5*sizeof(REAL *));

    // compute the rational approximation using AAA algorithms
    REAL err_max=get_cpzwf(frac_inv, (void *)func_param,	rpnwf, &mbig, &mmax_in, &k, xmin_in, xmax_in, AAA_tol, print_level);

    // assign poles and residules
    pcdata->residues = dvec_create_p(k);
    pcdata->poles = dvec_create_p(k-1);
    array_cp(k, rpnwf[0], pcdata->residues->val);
    array_cp(k-1, rpnwf[1], pcdata->poles->val);


    /* --------------------------------------------- */

    // scaling mass matrix
    pcdata->scaled_M = dcsr_create_p(m, n, nnz_M);
    dcsr_cp(M, pcdata->scaled_M);
    dcsr_axm(pcdata->scaled_M, scaling_m);

    // get diagonal entries of the scaled mass matrix
    pcdata->diag_scaled_M = dvec_create_p(n);
    dcsr_getdiag(0, pcdata->scaled_M, pcdata->diag_scaled_M);


    //------------------------------------------------
    // Set up all AMG for shifted laplacians
    //------------------------------------------------
    INT npoles = k-1;
    pcdata->mgl = (AMG_data **)calloc(npoles, sizeof(AMG_data *));
    // assemble amg data for all shifted laplacians:
    // (scaling_a*A - poles[i] * scaling_m*M)
    for(i = 0; i < npoles; ++i) {
        pcdata->mgl[i] = amg_data_create(max_levels);
        dcsr_alloc(n, n, 0, &(pcdata->mgl[i][0].A));
        dcsr_add(A, scaling_a, M, -pcdata->poles->val[i]*scaling_m, &(pcdata->mgl[i][0].A));
        pcdata->mgl[i][0].b = dvec_create(n);
        pcdata->mgl[i][0].x = dvec_create(n);

        switch (amgparam->AMG_type) {

          case UA_AMG: // Unsmoothed Aggregation AMG
              if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
              status = amg_setup_ua(pcdata->mgl[i], amgparam);
              break;

          case SA_AMG: // Smoothed Aggregation AMG
              if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
              status = amg_setup_sa(pcdata->mgl[i], amgparam);
              break;

          default: // UA AMG
              if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
              status = amg_setup_ua(pcdata->mgl[i], amgparam);
              break;

        }

        if(status < 0)
        {
            printf("Unsuccessful AMG setup at pole %d with status = %d\n", i, status);
            return 0;
        }
    }

    //------------------------------------------------
    // setup preconditioner data
    //------------------------------------------------
    pcdata->amgparam = amgparam;

    // save scaled alpha and beta
    pcdata->scaled_alpha = scaled_alpha;
    pcdata->scaled_beta = scaled_beta;

    // clean
    // for(i = 0; i < 5; ++i){
    //     if (rpnwf[i]) free(rpnwf[i]);
    // }
    // if (rpnwf) free(rpnwf);

    return 1;
}


void precond_block_data_print(precond_block_data *precdata, SHORT flag)
{
    INT i;
    INT nb = precdata->poles->row;
    printf("poles: row = %d\n", precdata->poles->row);

    for (i=0; i<nb; i++)
    {

        if(precdata->mgl){
            if(precdata->mgl[i]){
                AMG_data *mgl_i = precdata->mgl[i];
                printf("Pole %d\n", i);
                if(&(mgl_i->A)) printf("A: row = %d, col = %d, nnz = %d\n", mgl_i->A.row, mgl_i->A.col, mgl_i->A.nnz);
                if(&(mgl_i->P)) printf("P: row = %d, col = %d, nnz = %d\n", mgl_i->P.row, mgl_i->P.col, mgl_i->P.nnz);
                if(&(mgl_i->R)) printf("R: row = %d, col = %d, nnz = %d\n", mgl_i->R.row, mgl_i->R.col, mgl_i->R.nnz);
                if(&(mgl_i->M)) printf("M: row = %d, col = %d, nnz = %d\n", mgl_i->M.row, mgl_i->M.col, mgl_i->M.nnz);

            }
        }

    }

    if(precdata->amgparam) {
            param_amg_print(precdata->amgparam);
    }

    if (flag == TRUE)
    {
      if(precdata->G)  printf("G: row = %d, col = %d, nnz = %d\n", precdata->G->row, precdata->G->col, precdata->G->nnz);
      if(precdata->K)  printf("K: row = %d, col = %d, nnz = %d\n", precdata->K->row, precdata->K->col, precdata->K->nnz);
      if(precdata->Gt) printf("Gt: row = %d, col = %d, nnz = %d\n", precdata->Gt->row, precdata->Gt->col, precdata->Gt->nnz);
      if(precdata->Kt) printf("Kt: row = %d, col = %d, nnz = %d\n", precdata->Kt->row, precdata->Kt->col, precdata->Kt->nnz);
    }

    if(precdata->scaled_M) printf("scaled M: row = %d, col = %d, nnz = %d\n", precdata->scaled_M->row, precdata->scaled_M->col, precdata->scaled_M->nnz);
    if(precdata->diag_scaled_M) printf("diag M: row = %d\n", precdata->diag_scaled_M->row);
    if(precdata->poles) printf("poles: row = %d\n", precdata->poles->row);
    if(precdata->residues) printf("residues: row = %d\n", precdata->residues->row);

    return;
}

REAL *array_calloc(const INT n)
{
    REAL *ar = (REAL*)calloc(n, sizeof(REAL));

    return ar;
}