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
/***********************************************************************************************/
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
  REAL16 *s,s1,s2,alpha,beta; //,f123;
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
  //  fprintf(stdout,"\nf-param: s=%Le,t=%Le,alpha=%Le,beta=%Le\n",s1,s2,alpha,beta);
  //  fflush(stdout);
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
    //not used:    const SHORT max_levels = amgparam->max_levels;
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
	  if ( prtlvl > PRINT_NONE ) fprintf(stdout,"\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl, amgparam);
        break;

        case SA_AMG: // Smoothed Aggregation AMG
	  if ( prtlvl > PRINT_NONE ) fprintf(stdout,"\n Calling SA AMG ...\n");
            status = amg_setup_sa(mgl, amgparam);
        break;

        default: // Unsmoothed Aggregation AMG
	  if ( prtlvl > PRINT_NONE ) fprintf(stdout,"\n Calling UA AMG ...\n");
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
    //not used: const SHORT max_levels = amgparam->max_levels;
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
	  if ( prtlvl > PRINT_NONE ) fprintf(stdout,"\n Calling UA AMG ...\n");
            status = famg_setup_ua(mgl, amgparam);
        break;

        case SA_AMG: // Smoothed Aggregation AMG
	  if ( prtlvl > PRINT_NONE ) fprintf(stdout,"\n Calling SA AMG ...\n");
            status = famg_setup_sa(mgl, amgparam);
        break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) fprintf(stdout,"\n Calling UA AMG ...\n");
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


/***********************************************************************************************/
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


/***********************************************************************************************/
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


/***********************************************************************************************/
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


/***********************************************************************************************/
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
    //    const INT m = A->row, n = A->col, nnz = A->nnz, nnz_M = M->nnz;
    const INT m = A->row, n = A->col, nnz_M = M->nnz;
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
    fprintf(stdout,"\nf-param: scaled_alpha=%e,scaled_beta=%e,func_param[2]=%Le,func_param[3]=%Le\n",scaled_alpha,scaled_beta,func_param[2],func_param[3]);
    fflush(stdout);
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
    //    REAL err_max=get_cpzwf(frac_inv, (void *)func_param,rpnwf, &mbig, &mmax_in, &k, xmin_in, xmax_in, AAA_tol, print_level);
    get_cpzwf(frac_inv, (void *)func_param,rpnwf, &mbig, &mmax_in, &k, xmin_in, xmax_in, AAA_tol, print_level);
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
              if ( prtlvl > PRINT_NONE ) fprintf(stdout,"\n Calling UA AMG ...\n");
              status = amg_setup_ua(pcdata->mgl[i], amgparam);
              break;

          case SA_AMG: // Smoothed Aggregation AMG
              if ( prtlvl > PRINT_NONE ) fprintf(stdout,"\n Calling SA AMG ...\n");
              status = amg_setup_sa(pcdata->mgl[i], amgparam);
              break;

          default: // UA AMG
              if ( prtlvl > PRINT_NONE ) fprintf(stdout,"\n Calling UA AMG ...\n");
              status = amg_setup_ua(pcdata->mgl[i], amgparam);
              break;

        }

        if(status < 0)
        {
            fprintf(stdout,"Unsuccessful AMG setup at pole %d with status = %d\n", i, status);
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


/***********************************************************************************************/
/*!
 * \fn void precond_ra_data_null (precond_ra_data *pcdata)
 *
 * \brief Initialize (with NULL) data of precond_ra_data
 *
 * \param pcdata    Pointer to precond_ra_data
 *
 * \return
 *
 * \author          Ana Budisa
 * \date            2021-03-02
 */
void precond_ra_data_null(precond_ra_data *precdata)
{

#if WITH_SUITESPARSE
    precdata->LU_diag = NULL;
#endif

    precdata->mgl = NULL;
    precdata->amgparam = NULL;

    precdata->scaled_A = NULL;
    precdata->scaled_M = NULL;
    precdata->diag_scaled_M = NULL;
    precdata->poles = NULL;
    precdata->residues = NULL;

    precdata->r = NULL;
    precdata->w = NULL;

}

/*!
 * \fn precond_data *precond_ra_data_alloc (SHORT max_size)
 *
 * \brief Allocate precond_ra_data array of length max_size; each component
 *        is initialized
 *
 * \param max_size  Size of precond_ra_data array (usually 1)
 *
 * \return pcdata   Pointer to precond_ra_data (callocated)
 *
 * \author          Ana Budisa
 * \date            2021-02-03
 */
precond_ra_data *precond_ra_data_alloc(SHORT max_size)
{
    max_size = MAX(1, max_size);

    precond_ra_data *pcdata = (precond_ra_data*)calloc(max_size, sizeof(precond_ra_data));

    INT i;
    for(i = 0; i < max_size; ++i) precond_ra_data_null(&(pcdata[i]));

    return(pcdata);
}


/***********************************************************************************************/
/*!
 * \fn INT fenics_precond_ra_setup (dCSRmat *A, dCSRmat *M, REAL s_frac_power,
                                    REAL t_frac_power, REAL alpha, REAL beta,
                                    REAL scaling_a, REAL scaling_m, AMG_param *amgparam)
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
 * \date            2021-02-03
 */
INT fenics_precond_ra_data_setup(dCSRmat *A,
                                 dCSRmat *M,
                                 REAL s_frac_power,
                                 REAL t_frac_power,
                                 REAL alpha,
                                 REAL beta,
                                 REAL scaling_a,
                                 REAL scaling_m,
                                 AMG_param *amgparam,
                                 precond_ra_data *pcdata)
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
    REAL16 AAA_tol=powl(2e0,-46e0);  // tolerance of the AAA algorithm
    INT k=-22; // k is the number of nodes in the final interpolation after tolerance is achieved or mmax is reached.
    INT print_level=0; // print level for AAA

    // output of the AAA algorithm.  It contains residues, poles, nodes, weights, function values
    REAL **rpnwf=malloc(5*sizeof(REAL *));

    // compute the rational approximation using AAA algorithms
    //    REAL err_max=get_cpzwf(frac_inv, (void *)func_param, rpnwf, &mbig, &mmax_in, &k, xmin_in, xmax_in, AAA_tol, print_level);
    get_cpzwf(frac_inv, (void *)func_param, rpnwf, &mbig, &mmax_in, &k, xmin_in, xmax_in, AAA_tol, print_level);
    if(rpnwf == NULL) {
      fprintf(stderr,"\nUnsuccessful AAA computation of rational approximation\n");
      fflush(stderr);
      return 0;
    }

    // assign poles and residules
    pcdata->residues = dvec_create_p(k);
    pcdata->poles = dvec_create_p(k-1);
    array_cp(k, rpnwf[0], pcdata->residues->val);
    array_cp(k-1, rpnwf[1], pcdata->poles->val);
    //////// pcdata->residues->val[0]+=2e-15;
    REAL polez;
    for(i = 0; i < k-1; ++i) {
      polez=pcdata->poles->val[i];
      if(polez>0e0){
	fprintf(stderr,"\n%%%%%% *** HAZMATH WARNING*** Positive pole in function=%s", \
		__FUNCTION__);
	fprintf(stdout,"\n%%%%%%  0 < pole(%d)=%.16e\n", i, polez);
	break;
      }
    }
    /*if(prtlvl>5){
      fprintf(stdout,"\n%%%%%%params:");
      fprintf(stdout,"\ns=%.2Le;t=%.2Le;alpha=%.8Le;beta=%.8Le;",	\
	      func_param[0],func_param[1],func_param[2],func_param[3]);
      fprintf(stdout,"\n%%%%%%  POLES:\n");
      for(i = 0; i < k-1; ++i)
	fprintf(stdout,"pole(%d)=%.16e;\n", i+1, pcdata->poles->val[i]);
      fprintf(stdout,"\n%%%%%%  RESIDUES:\n");
      for(i = 0; i < k; ++i)
        fprintf(stdout,"res(%d)=%.16e;\n", i+1, pcdata->residues->val[i]);
    }
    fprintf(stdout,"\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");*/
    /* --------------------------------------------- */
    // scaling stiffness matrix
    pcdata->scaled_A = dcsr_create_p(m, n, nnz);
    dcsr_cp(A, pcdata->scaled_A);
    dcsr_axm(pcdata->scaled_A, scaling_a);

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
    fprintf(stdout,"\nNumber of poles: %d\n", npoles);
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
              if ( prtlvl > PRINT_NONE ) fprintf(stdout,"\n Calling UA AMG ...\n");
              status = amg_setup_ua(pcdata->mgl[i], amgparam);
              break;

          case SA_AMG: // Smoothed Aggregation AMG
              if ( prtlvl > PRINT_NONE ) fprintf(stdout,"\n Calling SA AMG ...\n");
              status = amg_setup_sa(pcdata->mgl[i], amgparam);
              break;

          default: // UA AMG
              if ( prtlvl > PRINT_NONE ) fprintf(stdout,"\n Calling UA AMG ...\n");
              status = amg_setup_ua(pcdata->mgl[i], amgparam);
              break;

        }

        if(status < 0)
        {
            fprintf(stdout,"Unsuccessful AMG setup at pole %d with status = %d\n", i, status);
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
    pcdata->s_power = s_frac_power;
    pcdata->t_power = t_frac_power;

    return 1;
}


/***********************************************************************************************/
/*!
 * \fn INT precond_block_data_print(precond_block_data *precdata, SHORT flag)
 *
 * \brief Print precond_block_data
 *
 * \param precdata      Pointer to precond_block_data
 * \param flag          Flag for Maxwell data
 *
 * \return
 *
 * \author          Ana Budisa
 * \date            2021-01-21
 */
void precond_block_data_print(precond_block_data *precdata, SHORT flag)
{
    INT i;
    INT nb = precdata->poles->row;
    fprintf(stdout,"poles: row = %d\n", precdata->poles->row);

    for (i=0; i<nb; i++)
    {

        if(precdata->mgl){
            if(precdata->mgl[i]){
                AMG_data *mgl_i = precdata->mgl[i];
                fprintf(stdout,"Pole %d\n", i);
                if(&(mgl_i->A)) fprintf(stdout,"A: row = %d, col = %d, nnz = %d\n", mgl_i->A.row, mgl_i->A.col, mgl_i->A.nnz);
                if(&(mgl_i->P)) fprintf(stdout,"P: row = %d, col = %d, nnz = %d\n", mgl_i->P.row, mgl_i->P.col, mgl_i->P.nnz);
                if(&(mgl_i->R)) fprintf(stdout,"R: row = %d, col = %d, nnz = %d\n", mgl_i->R.row, mgl_i->R.col, mgl_i->R.nnz);
                if(&(mgl_i->M)) fprintf(stdout,"M: row = %d, col = %d, nnz = %d\n", mgl_i->M.row, mgl_i->M.col, mgl_i->M.nnz);

            }
        }

    }

    if(precdata->amgparam) {
            param_amg_print(precdata->amgparam);
    }

    if (flag == TRUE)
    {
      if(precdata->G)  fprintf(stdout,"G: row = %d, col = %d, nnz = %d\n", precdata->G->row, precdata->G->col, precdata->G->nnz);
      if(precdata->K)  fprintf(stdout,"K: row = %d, col = %d, nnz = %d\n", precdata->K->row, precdata->K->col, precdata->K->nnz);
      if(precdata->Gt) fprintf(stdout,"Gt: row = %d, col = %d, nnz = %d\n", precdata->Gt->row, precdata->Gt->col, precdata->Gt->nnz);
      if(precdata->Kt) fprintf(stdout,"Kt: row = %d, col = %d, nnz = %d\n", precdata->Kt->row, precdata->Kt->col, precdata->Kt->nnz);
    }

    if(precdata->scaled_M) fprintf(stdout,"scaled M: row = %d, col = %d, nnz = %d\n", precdata->scaled_M->row, precdata->scaled_M->col, precdata->scaled_M->nnz);
    if(precdata->diag_scaled_M) fprintf(stdout,"diag M: row = %d\n", precdata->diag_scaled_M->row);
    if(precdata->poles) fprintf(stdout,"poles: row = %d\n", precdata->poles->row);
    if(precdata->residues) fprintf(stdout,"residues: row = %d\n", precdata->residues->row);

    return;
}


/***********************************************************************************************/
/*!
 * \fn REAL *array_calloc(const INT n)
 *
 * \brief Callocate an array of n doubles
 *
 * \param n      Array size
 *
 * \return       Pointer to array (REAL*)
 *
 * \author          Ana Budisa
 * \date            2021-02-02
 */
REAL *array_calloc(const INT n)
{
    REAL *ar = (REAL*)calloc(n, sizeof(REAL));

    return ar;
}


/***********************************************************************************************/
/*!
 * \fn dCSRmat *dcsr_calloc(const INT n)
 *
 * \brief Callocate an array of n dCSRmat matrices
 *
 * \param n      Array size
 *
 * \return       Pointer to array (dCSRmat*)
 *
 * \author          Ana Budisa
 * \date            2021-03-01
 */
dCSRmat *dcsr_calloc(const INT n)
{
    dCSRmat *ar = (dCSRmat*)calloc(n, sizeof(dCSRmat));

    return ar;
}


/***********************************************************************************************/
/*!
 * \fn void stupid_append_function(const INT k, dCSRmat *A_n, dCSRmat *A_diag_ptr)
 *
 * \brief
 *
 * \param n
 *
 * \return
 *
 * \author          Ana Budisa
 * \date            2021-03-01
 */
void stupid_append_function(const INT k, dCSRmat *A_n, dCSRmat *A_diag_ptr)
{
    INT m = A_n->row, n = A_n->col, nnz = A_n->nnz;

    dcsr_alloc(m, n, nnz, &(A_diag_ptr[k]));

    dcsr_cp(A_n, &(A_diag_ptr[k]));
}


/***********************************************************************************************/
/*!
 * \fn HX_curl_data *HX_curl_data_alloc (SHORT max_size)
 *
 * \brief Allocate HX_curl_data array of length max_size; each component
 *        is initialized
 *
 * \param max_size  Size of HX_curl_data array (usually 1)
 *
 * \return pcdata   Pointer to HX_curl_data (callocated)
 *
 * \author          Ana Budisa
 * \date            2021-03-18
 */
HX_curl_data *HX_curl_data_alloc(SHORT max_size)
{
    max_size = MAX(1, max_size);

    HX_curl_data *hxcurldata = (HX_curl_data*)calloc(max_size, sizeof(HX_curl_data));

    INT i;
    for(i = 0; i < max_size; ++i) HX_curl_data_null(&(hxcurldata[i]));

    return(hxcurldata);
}


/***********************************************************************************************/
/*!
 * \fn void fenics_HX_curl_data_setup(dCSRmat *Acurl, dCSRmat *Pcurl, dCSRmat *Grad,
 *                                    AMG_param *amgparam, HX_curl_data *hxcurldata)
 *
 * \brief Setup HX_curl_data structure for HX curl precond
 *
 * \param Acurl         Pointer to dCSRmat curl-curl matrix
 * \param Pcurl         Pointer to dCSRmat edges-to-nodes projection matrix
 * \param Grad          Pointer to dCSRmat grad operator matrix
 * \param amgparam      Pointer to AMG_param
 * \param hxcurldata    Pointer to HX_curl_data
 *
 * \return          INT (1 if successful setup, 0 else)
 *
 * \author          Ana Budisa
 * \date            2021-03-18
 * TODO: CHANGE LOCAL MATRICES INTO POINTERS!!
 */
int fenics_HX_curl_data_setup(dCSRmat *Acurl,
                              dCSRmat *Pcurl,
                              dCSRmat *Grad,
                              AMG_param *amgparam,
                              HX_curl_data *hxcurldata)
{
    const SHORT prtlvl = amgparam->print_level;
    const SHORT max_levels = amgparam->max_levels;

    /*------------------------*/
    /* Local Variables */
    /*------------------------*/
    INT      status = SUCCESS;

    /*------------------------*/
    /* setup vector Laplacian */
    /*------------------------*/

    // get transpose of P_curl
    dCSRmat Pt_curl;
    dcsr_trans(Pcurl, &Pt_curl);

    // get A_vgrad
    dCSRmat A_vgrad;
    dcsr_rap(&Pt_curl, Acurl, Pcurl, &A_vgrad);

    // initialize A, b, x for mgl_vgrad[0]
    AMG_data *mgl_vgrad = amg_data_create(max_levels);
    mgl_vgrad[0].A = dcsr_create(A_vgrad.row, A_vgrad.col, A_vgrad.nnz);
    dcsr_cp(&A_vgrad, &mgl_vgrad[0].A);
    mgl_vgrad[0].b = dvec_create(A_vgrad.row);
    mgl_vgrad[0].x = dvec_create(A_vgrad.col);

    // setup AMG for vector Laplacian
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_vgrad, amgparam); break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = amg_setup_sa(mgl_vgrad, amgparam); break;

        default: // Classical AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_vgrad, amgparam); break;

    }

    if(status < 0)
    {
        fprintf(stdout,"Unsuccessful AMG setup for vector Laplacian with status = %d\n", status);
        return 0;
    }

    /*------------------------*/
    /* setup scalar Laplacian */
    /*------------------------*/

    // get transpose of Grad
    dCSRmat Gradt;
    dcsr_trans(Grad, &Gradt);

    // get A_grad
    dCSRmat A_grad;
    dcsr_rap(&Gradt, Acurl, Grad, &A_grad);

    // initialize A, b, x for mgl_grad[0]
    AMG_data *mgl_grad = amg_data_create(max_levels);
    mgl_grad[0].A = dcsr_create(A_grad.row, A_grad.col, A_grad.nnz);
    dcsr_cp(&A_grad, &mgl_grad[0].A);
    mgl_grad[0].b = dvec_create(A_grad.row);
    mgl_grad[0].x = dvec_create(A_grad.col);

    // setup AMG for scalar Laplacian
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_grad, amgparam);
        break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = amg_setup_ua(mgl_grad, amgparam);
        break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_grad, amgparam);
        break;

    }

    if(status < 0)
    {
        fprintf(stdout,"Unsuccessful AMG setup for scalar Laplacian with status = %d\n", status);
        return 0;
    }

    /*------------------------*/
    hxcurldata->A = Acurl;

    hxcurldata->smooth_type = 1;
    hxcurldata->smooth_iter = 1;

    hxcurldata->P_curl = Pcurl;
    hxcurldata->Pt_curl = &Pt_curl;
    hxcurldata->A_vgrad = &A_vgrad;
    hxcurldata->amgparam_vgrad = amgparam;
    hxcurldata->mgl_vgrad = mgl_vgrad;

    hxcurldata->Grad = Grad;
    hxcurldata->Gradt = &Gradt;
    hxcurldata->A_grad = &A_grad;
    hxcurldata->amgparam_grad = amgparam;
    hxcurldata->mgl_grad = mgl_grad;

    hxcurldata->backup_r = (REAL*)calloc(Acurl->row, sizeof(REAL));
    hxcurldata->w = (REAL*)calloc(Acurl->row, sizeof(REAL));

    return 1;
}


/***********************************************************************************************/
/*!
 * \fn HX_div_data *HX_div_data_alloc (SHORT max_size)
 *
 * \brief Allocate HX_div_data array of length max_size; each component
 *        is initialized
 *
 * \param max_size  Size of HX_div_data array (usually 1)
 *
 * \return pcdata   Pointer to HX_div_data (callocated)
 *
 * \author          Ana Budisa
 * \date            2021-03-18
 */
HX_div_data *HX_div_data_alloc(SHORT max_size)
{
    max_size = MAX(1, max_size);

    HX_div_data *hxdivdata = (HX_div_data*)calloc(max_size, sizeof(HX_div_data));

    INT i;
    for(i = 0; i < max_size; ++i) HX_div_data_null(&(hxdivdata[i]));

    return(hxdivdata);
}


/***********************************************************************************************/
/*!
 * \fn void fenics_HX_div_data_3D_setup(dCSRmat *Adiv, dCSRmat *P_div, dCSRmat *Curl,
 *                                      dCSRmat *P_curl, AMG_param *amgparam,
 *                                      HX_div_data *hxdivdata)
 *
 * \brief Setup HX_div_data structure for HX div precond in 3D
 *
 * \param Adiv          Pointer to dCSRmat div-div matrix
 * \param P_div         Pointer to dCSRmat faces-to-nodes projection matrix
 * \param Curl          Pointer to dCSRmat curl operator matrix
 * \param P_curl        Pointer to dCSRmat edges-to-nodes projection matrix
 * \param amgparam      Pointer to AMG_param
 * \param hxdivdata     Pointer to HX_div_data
 *
 * \return          INT (1 if successful setup, 0 else)
 *
 * \author          Ana Budisa
 * \date            2021-03-18
 * TODO: CHANGE LOCAL MATRICES INTO POINTERS!!
 */
int fenics_HX_div_data_3D_setup(dCSRmat *Adiv,
                                dCSRmat *P_div,
                                dCSRmat *Curl,
                                dCSRmat *P_curl,
                                AMG_param *amgparam,
                                HX_div_data *hxdivdata)
{
    const SHORT prtlvl = amgparam->print_level;
    const SHORT max_levels = amgparam->max_levels;

    /*------------------------*/
    /* Local Variables */
    /*------------------------*/
    INT      status = SUCCESS;

    /*------------------------*/
    /* setup vector Laplacian */
    /*------------------------*/
    // get transpose of P_curl
    dCSRmat Pt_curl;
    dcsr_trans(P_curl, &Pt_curl);
    // get transpose of P_div
    dCSRmat Pt_div;
    dcsr_trans(P_div, &Pt_div);
    // get transpose of Curl
    dCSRmat Curlt;
    dcsr_trans(Curl, &Curlt);

    // get A_curl
    dCSRmat A_curl;
    dcsr_rap(&Curlt, Adiv, Curl, &A_curl);
    // get A_curlgrad
    dCSRmat A_curlgrad;
    dcsr_rap(&Pt_curl, &A_curl, P_curl, &A_curlgrad);
    // get A_divgrad
    dCSRmat A_divgrad;
    dcsr_rap(&Pt_div, Adiv, P_div, &A_divgrad);

    // initialize A, b, x for mgl_curlgrad[0]
    AMG_data *mgl_curlgrad = amg_data_create(max_levels);
    mgl_curlgrad[0].A = dcsr_create(A_curlgrad.row, A_curlgrad.col, A_curlgrad.nnz);
    dcsr_cp(&A_curlgrad, &mgl_curlgrad[0].A);
    mgl_curlgrad[0].b = dvec_create(A_curlgrad.row);
    mgl_curlgrad[0].x = dvec_create(A_curlgrad.col);

    // setup AMG for vector Laplacian
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_curlgrad, amgparam); break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_sa(mgl_curlgrad, amgparam); break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_curlgrad, amgparam); break;

    }

    if(status < 0)
    {
        fprintf(stdout,"Unsuccessful AMG setup for curlgrad Laplacian with status = %d\n", status);
        return 0;
    }

    // initialize A, b, x for mgl_divgrad[0]
    AMG_data *mgl_divgrad = amg_data_create(max_levels);
    mgl_divgrad[0].A = dcsr_create(A_divgrad.row, A_divgrad.col, A_divgrad.nnz);
    dcsr_cp(&A_divgrad, &mgl_divgrad[0].A);
    mgl_divgrad[0].b = dvec_create(A_divgrad.row);
    mgl_divgrad[0].x = dvec_create(A_divgrad.col);

    // setup AMG for vector Laplacian
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_divgrad, amgparam); break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_sa(mgl_divgrad, amgparam); break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_divgrad, amgparam); break;

    }

    if(status < 0)
    {
        fprintf(stdout,"Unsuccessful AMG setup for divgrad Laplacian with status = %d\n", status);
        return 0;
    }

    /*------------------------*/
    // setup preconditioner
    hxdivdata->A = Adiv;

    hxdivdata->smooth_type = 1;
    hxdivdata->smooth_iter = 1;

    hxdivdata->P_curl = P_curl;
    hxdivdata->Pt_curl = &Pt_curl;
    hxdivdata->P_div = P_div;
    hxdivdata->Pt_div = &Pt_div;
    hxdivdata->Curl = Curl;
    hxdivdata->Curlt = &Curlt;
    hxdivdata->A_curlgrad = &A_curlgrad;
    hxdivdata->A_divgrad = &A_divgrad;
    hxdivdata->amgparam_curlgrad = amgparam;
    hxdivdata->mgl_curlgrad = mgl_curlgrad;
    hxdivdata->amgparam_divgrad = amgparam;
    hxdivdata->mgl_divgrad = mgl_divgrad;

    hxdivdata->A_curl = &A_curl;
    hxdivdata->A_grad = NULL;
    hxdivdata->amgparam_grad = amgparam;
    hxdivdata->mgl_grad = NULL;

    hxdivdata->backup_r = (REAL*)calloc(Adiv->row, sizeof(REAL));
    hxdivdata->w = (REAL*)calloc(2*(A_curl.row), sizeof(REAL));

    return 1;
}


/***********************************************************************************************/
/*!
 * \fn void fenics_HX_div_data_2D_setup(dCSRmat *Adiv, dCSRmat *P_div, dCSRmat *Curl,
 *                                      AMG_param *amgparam, HX_div_data *hxdivdata)
 *
 * \brief Setup HX_div_data structure for HX div precond in 2D
 *
 * \param Adiv          Pointer to dCSRmat div-div matrix
 * \param P_div         Pointer to dCSRmat faces-to-nodes projection matrix
 * \param Curl          Pointer to dCSRmat curl operator matrix
 * \param amgparam      Pointer to AMG_param
 * \param hxdivdata     Pointer to HX_div_data
 *
 * \return          INT (1 if successful setup, 0 else)
 *
 * \author          Ana Budisa
 * \date            2021-03-18
 * TODO: CHANGE LOCAL MATRICES INTO POINTERS!!
 */
int fenics_HX_div_data_2D_setup(dCSRmat *Adiv,
                                dCSRmat *P_div,
                                dCSRmat *Curl,
                                AMG_param *amgparam,
                                HX_div_data *hxdivdata)
{
    const SHORT prtlvl = amgparam->print_level;
    const SHORT max_levels = amgparam->max_levels;

    /*------------------------*/
    /* Local Variables */
    /*------------------------*/
    INT status = SUCCESS;

    /*------------------------*/
    /* setup vector Laplacian */
    /*------------------------*/
    // get transpose of P_div
    dCSRmat Pt_div;
    dcsr_trans(P_div, &Pt_div);
    // get transpose of Curl
    dCSRmat Curlt;
    dcsr_trans(Curl, &Curlt);

    // get A_grad
    dCSRmat A_grad;
    dcsr_rap(&Curlt, Adiv, Curl, &A_grad);

    // get A_divgrad
    dCSRmat A_divgrad;
    dcsr_rap(&Pt_div, Adiv, P_div, &A_divgrad);

    // initialize A, b, x for mgl_grad[0]
    AMG_data *mgl_grad = amg_data_create(max_levels);
    mgl_grad[0].A = dcsr_create(A_grad.row, A_grad.col, A_grad.nnz);
    dcsr_cp(&A_grad, &mgl_grad[0].A);
    mgl_grad[0].b = dvec_create(A_grad.row);
    mgl_grad[0].x = dvec_create(A_grad.col);

    // setup AMG for scalar Laplacian
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_grad, amgparam); break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_sa(mgl_grad, amgparam); break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_grad, amgparam); break;

    }

    if(status < 0)
    {
        fprintf(stdout,"Unsuccessful AMG setup for grad Laplacian with status = %d\n", status);
        return 0;
    }

    // initialize A, b, x for mgl_divgrad[0]
    AMG_data *mgl_divgrad = amg_data_create(max_levels);
    mgl_divgrad[0].A = dcsr_create(A_divgrad.row, A_divgrad.col, A_divgrad.nnz);
    dcsr_cp(&A_divgrad, &mgl_divgrad[0].A);
    mgl_divgrad[0].b = dvec_create(A_divgrad.row);
    mgl_divgrad[0].x = dvec_create(A_divgrad.col);

    // setup AMG for vector Laplacian
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_divgrad, amgparam); break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_sa(mgl_divgrad, amgparam); break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_divgrad, amgparam); break;

    }

    if(status < 0)
    {
        fprintf(stdout,"Unsuccessful AMG setup for divgrad Laplacian with status = %d\n", status);
        return 0;
    }

    /*------------------------*/
    // setup preconditioner
    hxdivdata->A = Adiv;

    hxdivdata->smooth_type = 1;
    hxdivdata->smooth_iter = 1;

    hxdivdata->P_curl = NULL;
    hxdivdata->Pt_curl = NULL;
    hxdivdata->P_div = P_div;
    hxdivdata->Pt_div = &Pt_div;
    hxdivdata->Curl = Curl;
    hxdivdata->Curlt = &Curlt;
    hxdivdata->A_curlgrad = NULL;
    hxdivdata->A_divgrad = &A_divgrad;
    hxdivdata->amgparam_curlgrad = NULL;
    hxdivdata->mgl_curlgrad = NULL;
    hxdivdata->amgparam_divgrad = amgparam;
    hxdivdata->mgl_divgrad = mgl_divgrad;
    hxdivdata->A_curl = NULL;
    hxdivdata->A_grad = &A_grad;
    hxdivdata->amgparam_grad = amgparam;
    hxdivdata->mgl_grad = mgl_grad;

    hxdivdata->backup_r = (REAL*)calloc(Adiv->row, sizeof(REAL));
    hxdivdata->w = (REAL*)calloc(Adiv->row, sizeof(REAL));

    return 1;
}


/***********************************************************************************************/
/*!
 * \fn void smoother_data_null (smoother_data *smdata)
 *
 * \brief Initialize smoother_data (pointers are set to NULL)
 *
 * \param smdata   Smoother data structure
 *
 * \author         Ana Budisa
 * \date           2021-04-29
 */
void smoother_data_null (smoother_data *smdata)
{
    smdata->A = NULL;

}

/***********************************************************************************************/
/*!
 * \fn void smoother_data_free (smoother_data *smdata)
 *
 * \brief Free smoother_data
 *
 * \param smdata   Smoother data structure
 *
 * \author         Ana Budisa
 * \date           2021-04-29
 */
void smoother_data_free (smoother_data *smdata)
{
    if(smdata->A) dcsr_free(smdata->A);

}

/***********************************************************************************************/
/*!
 * \fn smoother_data *smoother_data_alloc (SHORT max_size)
 *
 * \brief Allocate smoother_data array of length max_size; each component
 *        is initialized
 *
 * \param max_size  Size of smoother_data array (usually 1)
 *
 * \return SMdata   Pointer to smoother_data (callocated)
 *
 * \author          Ana Budisa
 * \date            2021-04-29
 */
smoother_data *smoother_data_alloc(SHORT max_size)
{
    max_size = MAX(1, max_size);

    smoother_data *smdata = (smoother_data*)calloc(max_size, sizeof(smoother_data));

    INT i;
    for(i = 0; i < max_size; ++i) smoother_data_null(&(smdata[i]));

    return(smdata);
}


/***********************************************************************************************/
/*!
 * \fn void fenics_smoother_data_setup (AMG_data *mgl, AMG_param *amgparam, precond_data *pcdata)
 *
 * \brief Setup smoother_data structure from AMG_data and AMG_param
 *
 * \param istart    Start index of type INT
 * \param iend      End index of type INT
 * \param istep     Step size of type INT
 * \param nsweeps   Number of smoother iterations of type INT
 * \param relax     Relaxation parameter of type REAL
 * \param A         Pointer to dCSRmat (coefficient matrix for the smoother)
 * \param smdata    Pointer to smoother_data
 *
 * \author          Ana Budisa
 * \date            2021-04-29
 */
void fenics_smoother_data_setup(INT istart, INT iend, INT istep, INT nsweeps, REAL relax, dCSRmat *A, smoother_data *smdata)
{
    // Setup parameters
    smdata->istart = istart;
    smdata->iend = iend;
    smdata->istep = istep;
    smdata->nsweeps = nsweeps;
    smdata->relax = relax;

    smdata->A = dcsr_create_p(A->row, A->col, A->nnz);
    dcsr_cp(A, smdata->A);

}


/***********************************************************************************************/
/*!
 * \fn void smoother_matvec_null (smoother_matvec *smmv)
 *
 * \brief Initialize smoother_matvec (pointers are set to NULL)
 *
 * \param smmv   Smoother matvec structure
 *
 * \author       Ana Budisa
 * \date         2021-04-29
 */
void smoother_matvec_null (smoother_matvec *smmv)
{
    smmv->data = NULL;
    smmv->fct = NULL;

}

/***********************************************************************************************/
/*!
 * \fn smoother_matvec *smoother_matvec_alloc (SHORT max_size)
 *
 * \brief Allocate smoother_matvec array of length max_size; each component
 *        is initialized
 *
 * \param max_size  Size of smoother_matvec array (usually 1)
 *
 * \return smmv     Pointer to smoother_matvec (callocated)
 *
 * \author          Ana Budisa
 * \date            2021-04-29
 */
smoother_matvec *smoother_matvec_alloc(SHORT max_size)
{
    max_size = MAX(1, max_size);

    smoother_matvec *smmv = (smoother_matvec*)calloc(max_size, sizeof(smoother_matvec));

    INT i;
    for(i = 0; i < max_size; ++i) smoother_matvec_null(&(smmv[i]));

    return(smmv);
}


/***********************************************************************************************/
