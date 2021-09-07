//#include "math.h"
#include "hazmath.h"

/// A HACK(LUDMIL)

#ifndef NPY_INTP
#define NPY_INTP long
#endif

/// END A HACK(LUDMIL)

//#include "helper.hidden"

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
  REAL16 *s,s1,s2,alpha,beta;//,f123;
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
  // fprintf(stdout,"\nf-param: s=%Le,t=%Le,alpha=%Le,beta=%Le\n",s1,s2,alpha,beta);
  // fflush(stdout);
  return 1./(alpha*powl(x,s1)+beta*powl(x,s2));
}

/**/

/*---------------------------------*/
/*--      Public Functions      --*/
/*---------------------------------*/

dCSRmat* create_matrix(double *A, int nnz, int *ja, int nnz2, int *ia, int n)
{
  dCSRmat* mat;  
  mat = (dCSRmat *) calloc(1, sizeof(dCSRmat));
  /* for now not copy arrays, make test that produce seg. fault  */ 
  mat->row = n-1;
  mat->col = n-1;
  mat->nnz = nnz;
  mat->IA = (INT *)calloc(n, sizeof(INT)); // ia is of length nrow + 1
  mat->JA = (INT *)calloc(nnz, sizeof(INT)); // ja is of length nnz
  mat->val = (REAL *)calloc(nnz, sizeof(REAL)); // val is of lenght nnz
  iarray_cp(n, ia, mat->IA);
  iarray_cp(nnz, ja, mat->JA);
  array_cp(nnz, A, mat->val);

  return mat; 
}

dCSRmat* create_matrix2(double *A, int nnz, int *ja, int nnz2, int *ia, int n, int ncol)
{
  dCSRmat* mat;
  mat = (dCSRmat *) calloc(1, sizeof(dCSRmat));
  /* for now not copy arrays, make test that produce seg. fault  */
  mat->row = n-1;
  mat->col = ncol;
  mat->nnz = nnz;
  mat->IA = (INT *)calloc(n, sizeof(INT)); // ia is of length nrow + 1
  mat->JA = (INT *)calloc(nnz, sizeof(INT)); // ja is of length nnz
  mat->val = (REAL *)calloc(nnz, sizeof(REAL)); // val is of lenght nnz
  iarray_cp(n, ia, mat->IA);
  iarray_cp(nnz, ja, mat->JA);
  array_cp(nnz, A, mat->val);

  return mat;
}


dvector* create_dvector(double *x, int n)
{
  /* for now not copy arrays, make test that produce seg. fault  */ 
  dvector* vec; 
  vec = (dvector* ) calloc(2, sizeof(dvector)); 
  vec->row = n; 
  vec->val = x; 
  return vec; 
}

input_param* create_input_param()
{
  input_param* in_param = (input_param *)calloc(2, sizeof(input_param));
  param_input_init(in_param);
  param_input("./input.dat", in_param);
  return in_param; 
}

AMG_param* create_AMG_param(input_param *in_param)
{
  AMG_param* amg_param = (AMG_param *)calloc(2, sizeof(AMG_param));
  param_amg_init(amg_param); 
  param_amg_set(amg_param, in_param); 
  return amg_param; 
}

linear_itsolver_param* create_linear_itsolver_param(input_param *in_param)
{
  linear_itsolver_param* lin_param = (linear_itsolver_param*)calloc(2, sizeof(linear_itsolver_param));
  param_linear_solver_init(lin_param); 
  param_linear_solver_set(lin_param, in_param); 
  return lin_param; 
}

precond* set_precond(void *data, void (*foo)(REAL*, REAL*, void*))
{
    precond *pc = (precond*)calloc(1, sizeof(precond));

    pc->data = data;
    pc->fct = foo;

    return pc;
}

precond* create_precond(dCSRmat *A, AMG_param *amgparam)
{

    precond *pc = (precond*)calloc(1, sizeof(precond));

    precond_data *pcdata = (precond_data*)calloc(1, sizeof(precond_data));

    const INT nnz = A->nnz, m = A->row, n = A->col;
    short prtlvl = amgparam->print_level;
    INT      status = SUCCESS;

    AMG_data *mgl = amg_data_create(amgparam->max_levels);
    mgl[0].A = dcsr_create(m,n,nnz); dcsr_cp(A, &mgl[0].A);
    mgl[0].b = dvec_create(n); mgl[0].x = dvec_create(n);

    // setup preconditioner
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl, amgparam);
        break;

        case SA_AMG: // Smoothed Aggregation AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = amg_setup_sa(mgl, amgparam);
        break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl, amgparam);
        break;

    }

    if (status < 0) {
        printf("\n In function create_precond: AMG data failed to set up \n");
        return NULL;
    }

    // setup preconditioner
    param_amg_to_prec(pcdata, amgparam);
    pcdata->max_levels = mgl[0].num_levels;
    pcdata->mgl_data = mgl;

    pc->data = pcdata;
    switch (amgparam->cycle_type) {

        case V_CYCLE:
            pc->fct = precond_amg; break;
        case W_CYCLE:
            pc->fct = precond_amg; break;
        case AMLI_CYCLE:
            pc->fct = precond_amli; break;
        case NL_AMLI_CYCLE:
            pc->fct = precond_nl_amli; break;
        case ADD_CYCLE:
            pc->fct = precond_amg_add; break;
        default:
            pc->fct = precond_amg; break;

    }

    return pc;
}

precond* create_precond_amg(dCSRmat *A, AMG_param *amgparam)
{

    precond *pc = (precond*)calloc(1, sizeof(precond));

    precond_data *pcdata = (precond_data*)calloc(1, sizeof(precond_data));

    const INT nnz = A->nnz, m = A->row, n = A->col;
    short prtlvl = amgparam->print_level;
    INT      status = SUCCESS;

    AMG_data *mgl = amg_data_create(amgparam->max_levels);
    mgl[0].A = dcsr_create(m,n,nnz); dcsr_cp(A, &mgl[0].A);
    mgl[0].b = dvec_create(n); mgl[0].x = dvec_create(n);

    // setup AMG
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl, amgparam);
        break;

        case SA_AMG: // Smoothed Aggregation AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = amg_setup_sa(mgl, amgparam);
        break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl, amgparam);
        break;

    }

    if (status < 0) {
        printf("\n In function create_precond: AMG data failed to set up \n");
        return NULL;
    }

    /*// custom smoother
    smoother_matvec *smmv = smoother_matvec_alloc(1);
    smmv->type = amgparam->smoother; // there is no smoother_type, its smoother only; 
    if(amgparam->smoother_function) {
        smmv->fct = amgparam->smoother_function;
    }
    mgl->wdata = smmv; */

    // setup preconditioner
    param_amg_to_prec(pcdata, amgparam);
    pcdata->max_levels = mgl[0].num_levels;
    pcdata->mgl_data = mgl;

    pc->data = pcdata;
    switch (amgparam->cycle_type) {

        case V_CYCLE:
            pc->fct = precond_amg; break;
        case W_CYCLE:
            pc->fct = precond_amg; break;
        case AMLI_CYCLE:
            pc->fct = precond_amli; break;
        case NL_AMLI_CYCLE:
            pc->fct = precond_nl_amli; break;
        case ADD_CYCLE:
            pc->fct = precond_amg_add; break;
        default:
            pc->fct = precond_amg; break;

    }

    return pc;
}


precond* create_precond_famg(dCSRmat *A, dCSRmat *M, AMG_param *amgparam)
{

    precond *pc = (precond*)calloc(1, sizeof(precond));

    precond_data *pcdata = (precond_data*)calloc(1, sizeof(precond_data));

    const INT nnz = A->nnz, m = A->row, n = A->col, nnz_M = M->nnz;
    short prtlvl = amgparam->print_level;
    INT status = SUCCESS;

    AMG_data *mgl = amg_data_create(amgparam->max_levels);
    mgl[0].A = dcsr_create(m, n, nnz); dcsr_cp(A, &mgl[0].A);
    mgl[0].M = dcsr_create(m, n, nnz_M); dcsr_cp(M, &mgl[0].M);
    mgl[0].b = dvec_create(n); mgl[0].x = dvec_create(n);

    // setup AMG
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA FAMG ...\n");
            status = famg_setup_ua(mgl, amgparam);
        break;

        case SA_AMG: // Smoothed Aggregation AMG setup
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA FAMG ...\n");
            status = famg_setup_sa(mgl, amgparam);
        break;

        default: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA FAMG ...\n");
            status = famg_setup_ua(mgl, amgparam);
        break;

    }

    if (status < 0) {
        printf("\n In function create_precond: FAMG data failed to set up \n");
        return NULL;
    }

    // setup preconditioner
    param_amg_to_prec(pcdata, amgparam);
    pcdata->max_levels = mgl[0].num_levels;
    pcdata->mgl_data = mgl;

    pc->data = pcdata;
    switch (amgparam->cycle_type) {

        case V_CYCLE:
            pc->fct = precond_famg; break;
        case W_CYCLE:
            pc->fct = precond_famg; break;
        case AMLI_CYCLE:
            pc->fct = precond_famli; break;
        case ADD_CYCLE:
            pc->fct = precond_famg_add; break;
        default:
            pc->fct = precond_famg; break;

    }

    return pc;
}


precond* create_precond_ra(dCSRmat *A,
                           dCSRmat *M,
                           REAL s_frac_power,
                           REAL t_frac_power,
                           REAL alpha,
                           REAL beta,
                           REAL scaling_a,
                           REAL scaling_m,
                           AMG_param *amgparam)
{
    precond *pc = (precond*)calloc(1, sizeof(precond));

    precond_ra_data *pcdata = (precond_ra_data*)calloc(1, sizeof(precond_ra_data));

    const SHORT prtlvl = amgparam->print_level;
    const SHORT max_levels = amgparam->max_levels;
    const INT m = A->row, n = A->col, nnz = A->nnz, nnz_M = M->nnz;
    INT status = SUCCESS;
    INT i;
    //INT ii, jj;

    /*fprintf(stdout, "\n A: \n");
    fprintf(stdout,"\nrow, col, nnz: %d, %d, %d \n", A->row, A->col, A->nnz);
    for (ii=0; ii < A->row; ++ii) {
        fprintf(stdout, "\n");
        for (jj = A->IA[ii]; jj < A->IA[ii+1]; ++jj) {
            fprintf(stdout, "%.16e\t", A->val[jj]);
        }
    } // end for js
    fprintf(stdout, "\n");*/

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
    REAL16 AAA_tol=powl(2e0,-40e0);  // tolerance of the AAA algorithm
    INT k=-22; // k is the number of nodes in the final interpolation after tolerance is achieved or mmax is reached.
    INT print_level=0; // print level for AAA

    // output of the AAA algorithm.  It contains residues, poles, nodes, weights, function values
    REAL **rpnwf=malloc(5*sizeof(REAL *));

    // compute the rational approximation using AAA algorithms
    REAL err_max=get_cpzwf(frac_inv, (void *)func_param, rpnwf, &mbig, &mmax_in, &k, xmin_in, xmax_in, AAA_tol, print_level);
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
    if(prtlvl>5){
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
    //fprintf(stdout,"\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

    /* --------------------------------------------- */
    // scaling stiffness matrix
    pcdata->scaled_A = dcsr_create_p(m, n, nnz);
    dcsr_cp(A, pcdata->scaled_A);
    dcsr_axm(pcdata->scaled_A, scaling_a);

    /*fprintf(stdout, "\n A: \n");
    fprintf(stdout,"\nrow, col, nnz: %d, %d, %d \n", A->row, A->col, A->nnz);
    for (ii=0; ii < A->row; ++ii) {
        fprintf(stdout, "\n");
        for (jj = A->IA[i]; jj < A->IA[i+1]; ++jj) {
            fprintf(stdout, "%.16e\t", A->val[jj]);
        }
    } // end for js
    fprintf(stdout, "\n");*/

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
    /*fprintf(stdout, "\n A: \n");
    fprintf(stdout,"\nrow, col, nnz: %d, %d, %d \n", A->row, A->col, A->nnz);
    for (ii=0; ii < A->row; ++ii) {
        fprintf(stdout, "\n");
        for (jj = A->IA[i]; jj < A->IA[i+1]; ++jj) {
            fprintf(stdout, "%.16e\t", A->val[jj]);
        }
    } // end for js
    fprintf(stdout, "\n");*/
    /*fprintf(stdout, "\n M: \n");
    fprintf(stdout,"\nrow, col, nnz: %d, %d, %d \n", M->row, M->col, M->nnz);
    for (ii=0; ii < M->row; ++ii) {
        fprintf(stdout, "\n");
        for (jj = M->IA[i]; jj < M->IA[i+1]; ++jj) {
            fprintf(stdout, "%.16e\t", M->val[jj]);
        }
    } // end for js
    fprintf(stdout, "\n");*/


    for(i = 0; i < npoles; ++i) {
        //fprintf(stdout,"\nAMG for pole %d\n", i);
        pcdata->mgl[i] = amg_data_create(max_levels);
        dcsr_alloc(n, n, 0, &(pcdata->mgl[i][0].A));
        dcsr_add(A, scaling_a, M, -pcdata->poles->val[i]*scaling_m, &(pcdata->mgl[i][0].A));
        pcdata->mgl[i][0].b = dvec_create(n);
        pcdata->mgl[i][0].x = dvec_create(n);

        /*fprintf(stdout, "\n Matrix pole %d: \n", i);
        fprintf(stdout,"\nrow, col, nnz: %d, %d, %d \n", pcdata->mgl[i][0].A.row, pcdata->mgl[i][0].A.col, pcdata->mgl[i][0].A.nnz);
        for (ii=0; ii < pcdata->mgl[i][0].A.row; ++ii) {
            fprintf(stdout, "\n");
            for (jj = pcdata->mgl[i][0].A.IA[ii]; jj < pcdata->mgl[i][0].A.IA[ii+1]; ++jj) {
                fprintf(stdout, "%.16e\t", pcdata->mgl[i][0].A.val[jj]);
            }
        } // end for js
        fprintf(stdout, "\n");*/

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
    pcdata->amgparam = (AMG_param*)malloc(sizeof(AMG_param));
    param_amg_init(pcdata->amgparam);
    param_amg_cp(amgparam, pcdata->amgparam);
    //pcdata->amgparam = amgparam;

    // save scaled alpha and beta
    pcdata->scaled_alpha = scaled_alpha;
    pcdata->scaled_beta = scaled_beta;
    pcdata->s_power = s_frac_power;
    pcdata->t_power = t_frac_power;

    pc->data = pcdata;
    pc->fct = precond_ra_fenics;

    // clean
    if (rpnwf) free(rpnwf);

    return pc;
}


INT get_poles_no(precond* pc)
{
    precond_ra_data* data = (precond_ra_data*)(pc->data);

    return data->poles->row;
}


dvector* compute_ra_aaa(REAL s_frac_power,
                        REAL t_frac_power,
                        REAL alpha,
                        REAL beta,
                        REAL scaling_a,
                        REAL scaling_m)
{
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
    REAL xmin_in = 0.e0, xmax_in = 1.e0;  // interval for x
    INT mbig = (1<<14) + 1;  // initial number of points on the interval [x_min, x_max]
    INT mmax_in = (INT )(mbig/2);  // maximal final number of pole + 1
    REAL16 AAA_tol = powl(2e0,-40e0);  // tolerance of the AAA algorithm
    INT k = -22; // k is the number of nodes in the final interpolation after tolerance is achieved or mmax is reached.
    INT print_level = 0; // print level for AAA

    // output of the AAA algorithm.  It contains residues, poles, nodes, weights, function values
    REAL **rpnwf = malloc(5 * sizeof(REAL *));

    // compute the rational approximation using AAA algorithms
    REAL err_max = get_cpzwf(frac_inv, (void *)func_param, rpnwf, &mbig, &mmax_in, &k, xmin_in, xmax_in, AAA_tol, print_level);
    if(rpnwf == NULL) {
      fprintf(stderr,"\nUnsuccessful AAA computation of rational approximation\n");
      fflush(stderr);
      return 0;
    }
    printf("Approximation error: %.16e\n", err_max);
    printf("Number of poles: %d\n", k-1);

    // assign poles and residuals
    dvector *res = dvec_create_p(2*k - 1);
    array_cp(k-1, rpnwf[1], res->val);
    array_cp(k, rpnwf[0], &(res->val[k]));

    // check if poles are non negative
    REAL polez;
    for(i = 0; i < k-1; ++i) {
      polez = res->val[i];
      if(polez > 0e0) {
	    fprintf(stderr,"\n%%%%%% *** HAZMATH WARNING*** Positive pole in function=%s", \
	        __FUNCTION__);
	    fprintf(stdout,"\n%%%%%%  0 < pole(%d)=%.16e\n", i, polez);
	    break;
      }
    }

    // clean
    if(rpnwf) free(rpnwf);

    return res;
}


precond* create_precond_hxcurl(dCSRmat *Acurl,
                               dCSRmat *Pcurl,
                               dCSRmat *Grad,
                               SHORT prectype,
                               AMG_param *amgparam)
{
    precond *pc = (precond*)calloc(1, sizeof(precond));

    HX_curl_data *pcdata = (HX_curl_data*)calloc(1, sizeof(HX_curl_data));

    const SHORT prtlvl = amgparam->print_level;
    const SHORT max_levels = amgparam->max_levels;

    /*------------------------*/
    /* Local Variables */
    /*------------------------*/
    INT status = SUCCESS;

    /*------------------------*/
    /* setup vector Laplacian */
    /*------------------------*/

    // get transpose of P_curl
    dCSRmat *Pt_curl = dcsr_create_p(Pcurl->col, Pcurl->row, Pcurl->nnz);
    dcsr_trans(Pcurl, Pt_curl);

    // get A_vgrad
    dCSRmat *A_vgrad = (dCSRmat*)calloc(1, sizeof(dCSRmat));
    dcsr_rap(Pt_curl, Acurl, Pcurl, A_vgrad);

    // initialize A, b, x for mgl_vgrad[0]
    AMG_data *mgl_vgrad = amg_data_create(max_levels);
    mgl_vgrad[0].A = dcsr_create(A_vgrad->row, A_vgrad->col, A_vgrad->nnz);
    dcsr_cp(A_vgrad, &mgl_vgrad[0].A);
    mgl_vgrad[0].b = dvec_create(A_vgrad->row);
    mgl_vgrad[0].x = dvec_create(A_vgrad->col);

    // setup AMG for vector Laplacian
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_vgrad, amgparam); break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = amg_setup_sa(mgl_vgrad, amgparam); break;

        default: // Unsmoothed Aggregation AMG
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
    dCSRmat *Gradt = dcsr_create_p(Grad->col, Grad->row, Grad->nnz);
    dcsr_trans(Grad, Gradt);

    // get A_grad
    dCSRmat *A_grad = (dCSRmat*)calloc(1, sizeof(dCSRmat));
    dcsr_rap(Gradt, Acurl, Grad, A_grad);

    // initialize A, b, x for mgl_grad[0]
    AMG_data *mgl_grad = amg_data_create(max_levels);
    mgl_grad[0].A = dcsr_create(A_grad->row, A_grad->col, A_grad->nnz);
    dcsr_cp(A_grad, &mgl_grad[0].A);
    mgl_grad[0].b = dvec_create(A_grad->row);
    mgl_grad[0].x = dvec_create(A_grad->col);

    // setup AMG for scalar Laplacian
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_grad, amgparam);
        break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
            status = amg_setup_sa(mgl_grad, amgparam);
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

    // copy amgparam (this is for swig)
    AMG_param *amgparam2 = (AMG_param*)malloc(sizeof(AMG_param));
    param_amg_init(amgparam2);
    param_amg_cp(amgparam, amgparam2);

    /*------------------------*/
    pcdata->A = dcsr_create_p(Acurl->row, Acurl->col, Acurl->nnz); dcsr_cp(Acurl, pcdata->A);

    pcdata->smooth_type = 1;
    pcdata->smooth_iter = 1;

    pcdata->P_curl = dcsr_create_p(Pcurl->row, Pcurl->col, Pcurl->nnz); dcsr_cp(Pcurl, pcdata->P_curl);
    pcdata->Pt_curl = Pt_curl;
    pcdata->A_vgrad = A_vgrad;
    pcdata->amgparam_vgrad = amgparam2;
    pcdata->mgl_vgrad = mgl_vgrad;

    pcdata->Grad = dcsr_create_p(Grad->row, Grad->col, Grad->nnz); dcsr_cp(Grad, pcdata->Grad);
    pcdata->Gradt = Gradt;
    pcdata->A_grad = A_grad;
    pcdata->amgparam_grad = amgparam2;
    pcdata->mgl_grad = mgl_grad;

    pcdata->backup_r = (REAL*)calloc(Acurl->row, sizeof(REAL));
    pcdata->w = (REAL*)calloc(Acurl->row, sizeof(REAL));
    
    pc->data = pcdata;
    switch (prectype) {

        case PREC_HX_CURL_A:
            pc->fct = precond_hx_curl_additive; break;
        case PREC_HX_CURL_M:
            pc->fct = precond_hx_curl_multiplicative; break;
        default:
            pc->fct = precond_hx_curl_additive; break;

    }

    return pc;
}


precond* create_precond_hxdiv_3D(dCSRmat *Adiv, 
                                 dCSRmat *P_div, 
                                 dCSRmat *Curl,
                                 dCSRmat *P_curl,
                                 SHORT prectype,
                                 AMG_param *amgparam)
{   
    precond *pc = (precond*)calloc(1, sizeof(precond));

    HX_div_data *pcdata = (HX_div_data*)calloc(1, sizeof(HX_div_data));
    
    const SHORT prtlvl = amgparam->print_level;
    const SHORT max_levels = amgparam->max_levels;

    /*------------------------*/
    /* Local Variables */
    /*------------------------*/
    INT  status = SUCCESS;

    /*------------------------*/
    /* setup vector Laplacian */
    /*------------------------*/
    // get transpose of P_curl
    dCSRmat *Pt_curl = dcsr_create_p(P_curl->col, P_curl->row, P_curl->nnz);
    dcsr_trans(P_curl, Pt_curl);
    // get transpose of P_div
    dCSRmat *Pt_div = dcsr_create_p(P_div->col, P_div->row, P_div->nnz);
    dcsr_trans(P_div, Pt_div);
    // get transpose of Curl
    dCSRmat *Curlt = dcsr_create_p(Curl->col, Curl->row, Curl->nnz);
    dcsr_trans(Curl, Curlt);

    // get A_curl
    dCSRmat *A_curl = (dCSRmat*)calloc(1, sizeof(dCSRmat));
    dcsr_rap(Curlt, Adiv, Curl, A_curl);
    // get A_curlgrad
    dCSRmat *A_curlgrad = (dCSRmat*)calloc(1, sizeof(dCSRmat));
    dcsr_rap(Pt_curl, A_curl, P_curl, A_curlgrad);
    // get A_divgrad
    dCSRmat *A_divgrad = (dCSRmat*)calloc(1, sizeof(dCSRmat));
    dcsr_rap(Pt_div, Adiv, P_div, A_divgrad);

    // initialize A, b, x for mgl_curlgrad[0]
    AMG_data *mgl_curlgrad = amg_data_create(max_levels);
    mgl_curlgrad[0].A = dcsr_create(A_curlgrad->row, A_curlgrad->col, A_curlgrad->nnz);
    dcsr_cp(A_curlgrad, &mgl_curlgrad[0].A);
    mgl_curlgrad[0].b = dvec_create(A_curlgrad->row);
    mgl_curlgrad[0].x = dvec_create(A_curlgrad->col);

    // setup AMG for vector Laplacian
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_curlgrad, amgparam); break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
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
    mgl_divgrad[0].A = dcsr_create(A_divgrad->row, A_divgrad->col, A_divgrad->nnz);
    dcsr_cp(A_divgrad, &mgl_divgrad[0].A);
    mgl_divgrad[0].b = dvec_create(A_divgrad->row);
    mgl_divgrad[0].x = dvec_create(A_divgrad->col);

    // setup AMG for vector Laplacian
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_divgrad, amgparam); break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
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

    // copy amgparam (this is for swig)
    AMG_param *amgparam2 = (AMG_param*)malloc(sizeof(AMG_param));
    param_amg_init(amgparam2);
    param_amg_cp(amgparam, amgparam2);

    /*------------------------*/
    // setup preconditioner
    pcdata->A = dcsr_create_p(Adiv->row, Adiv->col, Adiv->nnz); dcsr_cp(Adiv, pcdata->A);

    pcdata->smooth_type = 1;
    pcdata->smooth_iter = 1;

    pcdata->P_curl = dcsr_create_p(P_curl->row, P_curl->col, P_curl->nnz); dcsr_cp(P_curl, pcdata->P_curl);
    pcdata->Pt_curl = Pt_curl;
    pcdata->P_div = dcsr_create_p(P_div->row, P_div->col, P_div->nnz); dcsr_cp(P_div, pcdata->P_div);
    pcdata->Pt_div = Pt_div;
    pcdata->Curl = dcsr_create_p(Curl->row, Curl->col, Curl->nnz); dcsr_cp(Curl, pcdata->Curl);
    pcdata->Curlt = Curlt;
    pcdata->A_curlgrad = A_curlgrad;
    pcdata->A_divgrad = A_divgrad;
    pcdata->amgparam_curlgrad = amgparam2;
    pcdata->mgl_curlgrad = mgl_curlgrad;
    pcdata->amgparam_divgrad = amgparam2;
    pcdata->mgl_divgrad = mgl_divgrad;

    pcdata->A_curl = A_curl;
    pcdata->A_grad = NULL;
    pcdata->amgparam_grad = amgparam2;
    pcdata->mgl_grad = NULL;

    pcdata->backup_r = (REAL*)calloc(Adiv->row, sizeof(REAL));
    pcdata->w = (REAL*)calloc(2*(A_curl->row), sizeof(REAL));
    
    pc->data = pcdata;
    switch (prectype) {

        case PREC_HX_DIV_A:
            pc->fct = precond_hx_div_additive; break;
        case PREC_HX_DIV_M:
            pc->fct = precond_hx_div_multiplicative; break;
        default:
            pc->fct = precond_hx_div_additive; break;

    }

    return pc;
}


precond* create_precond_hxdiv_2D(dCSRmat *Adiv,
                                 dCSRmat *P_div,
                                 dCSRmat *Curl,
                                 SHORT prectype,
                                 AMG_param *amgparam)
{
    precond *pc = (precond*)calloc(1, sizeof(precond));

    HX_div_data *pcdata = (HX_div_data*)calloc(1, sizeof(HX_div_data));
    
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
    dCSRmat *Pt_div = dcsr_create_p(P_div->col, P_div->row, P_div->nnz);
    dcsr_trans(P_div, Pt_div);
    // get transpose of Curl
    dCSRmat *Curlt = dcsr_create_p(Curl->col, Curl->row, Curl->nnz);
    dcsr_trans(Curl, Curlt);

    // get A_grad
    dCSRmat *A_grad = (dCSRmat*)calloc(1, sizeof(dCSRmat));
    dcsr_rap(Curlt, Adiv, Curl, A_grad);

    // get A_divgrad
    dCSRmat *A_divgrad = (dCSRmat*)calloc(1, sizeof(dCSRmat));;
    dcsr_rap(Pt_div, Adiv, P_div, A_divgrad);

    // initialize A, b, x for mgl_grad[0]
    AMG_data *mgl_grad = amg_data_create(max_levels);
    mgl_grad[0].A = dcsr_create(A_grad->row, A_grad->col, A_grad->nnz);
    dcsr_cp(A_grad, &mgl_grad[0].A);
    mgl_grad[0].b = dvec_create(A_grad->row);
    mgl_grad[0].x = dvec_create(A_grad->col);

    // setup AMG for scalar Laplacian
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_grad, amgparam); break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
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
    mgl_divgrad[0].A = dcsr_create(A_divgrad->row, A_divgrad->col, A_divgrad->nnz);
    dcsr_cp(A_divgrad, &mgl_divgrad[0].A);
    mgl_divgrad[0].b = dvec_create(A_divgrad->row);
    mgl_divgrad[0].x = dvec_create(A_divgrad->col);

    // setup AMG for vector Laplacian
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unsmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
            status = amg_setup_ua(mgl_divgrad, amgparam); break;

        case SA_AMG: // Smoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
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

    // copy amgparam (this is for swig)
    AMG_param *amgparam2 = (AMG_param*)malloc(sizeof(AMG_param));
    param_amg_init(amgparam2);
    param_amg_cp(amgparam, amgparam2);

    /*------------------------*/
    // setup preconditioner
    pcdata->A = dcsr_create_p(Adiv->row, Adiv->col, Adiv->nnz); dcsr_cp(Adiv, pcdata->A);

    pcdata->smooth_type = 1;
    pcdata->smooth_iter = 1;

    pcdata->P_curl = NULL;
    pcdata->Pt_curl = NULL;
    pcdata->P_div = dcsr_create_p(P_div->row, P_div->col, P_div->nnz); dcsr_cp(P_div, pcdata->P_div);
    pcdata->Pt_div = Pt_div;
    pcdata->Curl = dcsr_create_p(Curl->row, Curl->col, Curl->nnz); dcsr_cp(Curl, pcdata->Curl);
    pcdata->Curlt = Curlt;
    pcdata->A_curlgrad = NULL;
    pcdata->A_divgrad = A_divgrad;
    pcdata->amgparam_curlgrad = NULL;
    pcdata->mgl_curlgrad = NULL;
    pcdata->amgparam_divgrad = amgparam2;
    pcdata->mgl_divgrad = mgl_divgrad;
    pcdata->A_curl = NULL;
    pcdata->A_grad = A_grad;
    pcdata->amgparam_grad = amgparam2;
    pcdata->mgl_grad = mgl_grad;

    pcdata->backup_r = (REAL*)calloc(Adiv->row, sizeof(REAL));
    pcdata->w = (REAL*)calloc(Adiv->row, sizeof(REAL));

    pc->data = pcdata;
    switch (prectype) {

        case PREC_HX_DIV_A:
            pc->fct = precond_hx_div_additive_2D; break;
        case PREC_HX_DIV_M:
            pc->fct = precond_hx_div_multiplicative_2D; break;
        default:
            pc->fct = precond_hx_div_additive_2D; break;

    }

    return pc;
}

void apply_precond(REAL *r, REAL *z, precond *pc)
{
    //printf("calling pc->fct \n");
    pc->fct(r, z, pc->data);
    //printf("done calling pc->fct \n");

}


PyObject* py_callback_setup(PyObject* pyfunc, AMG_param *amgparam)
{
    printf("Here I am inside the callback setup \n");

    if(!PyCallable_Check(pyfunc)) {
        PyErr_SetString(PyExc_TypeError, "Parameter must be callable!");
        return NULL;
    }
    Py_XINCREF(pyfunc);
// LUDMIL: changed below because this makes no sense. 
//    (PyObject*)amgparam->smoother_function = pyfunc;
    amgparam->smoother_function = (void *)pyfunc;
    Py_INCREF(Py_None);

    return Py_None;
}


PyObject* py_callback_eval(REAL *r, REAL *x, smoother_matvec *smmv)
{
    printf("Here I am inside the callback eval\n");

    //import_array();
    // get data
    //LUDMIL    npy_intp D[1]; D[0] = smmv->data->A->col; t// 
    // LUDMIL Not exactly sure what this (ABOVE) is supposed to do (find number of columns?), so I rewrote it
    /* INT D[1];//LUDMIL */
    /* dCSRmat *A=(dCSRmat *)smmv->data; //LUDMIL */
    /* D[0] = A->col; //LUDMIL */
    
    
    // create new Python arrays
    //LUDMIL: this NPY_DOUBLE is not the right thing, so commenting out
    PyObject *rr=NULL; //LUDMIL ?   
    PyObject *xx=NULL;  //LUDMIL? 
    //LUDMIL     PyObject *rr = PyArray_SimpleNewFromData(1, D, NPY_DOUBLE, (void*)r);
    //LUDMIL     PyObject *xx = PyArray_SimpleNewFromData(1, D, NPY_DOUBLE, (void*)x);

    // memory management
    Py_INCREF(rr);
    Py_INCREF(xx);

    // get the function
    PyObject *pyfunc = (PyObject*)smmv->fct;
    // check if callable
    if(!PyCallable_Check(pyfunc)) {
        PyErr_SetString(PyExc_TypeError, "Smoother matvec function is not (python) callable!");
        return NULL;
    }
    Py_INCREF(pyfunc);

    // Build up the argument list...
    PyObject *arglist = Py_BuildValue("(OO)", rr, xx);

    // ...for calling the Python function
    PyObject *result = PyEval_CallObject(pyfunc, arglist);

    Py_DECREF(arglist);

    PyObject *return_obj = PyUnicode_FromString("everything is gone be ok");
    return return_obj;

}
