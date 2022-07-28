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
inline static REAL16 frac_inv(REAL16 x, REAL16 s1, REAL16 s2, REAL16 alpha, REAL16 beta)
{
  // fprintf(stdout,"\nf-param: s=%Le,t=%Le,alpha=%Le,beta=%Le\n",s1,s2,alpha,beta);
  // fflush(stdout);
  return 1./(alpha*powl(x,s1)+beta*powl(x,s2));
}

/**/

/*---------------------------------*/
/*--      Public Functions      --*/
/*---------------------------------*/



dCSRmat* create_matrix(double *A, int nnz, int *ja, int nnz2, int *ia, int n, int ncol)
{
  dCSRmat *mat = (dCSRmat *)calloc(1, sizeof(dCSRmat));
  mat->row = n-1;
  mat->col = ncol;
  mat->nnz = nnz;
  mat->IA = (INT *)calloc(n, sizeof(INT)); // ia is of length nrow + 1
  mat->JA = (INT *)calloc(nnz, sizeof(INT)); // ja is of length nnz
  mat->val = (REAL *)calloc(nnz, sizeof(REAL)); // val is of length nnz
  iarray_cp(n, ia, mat->IA);
  iarray_cp(nnz, ja, mat->JA);
  array_cp(nnz, A, mat->val);

  return mat;
}

dCOOmat* create_matrix_coo(double *A, int nnz, int *ja, int nnz2, int *ia, int n, int ncol)
{
  dCOOmat *mat = (dCOOmat *)calloc(1, sizeof(dCOOmat));
  int i, ij;
  mat->row = n-1;
  mat->col = ncol;
  mat->nnz = nnz;
  mat->rowind = (INT *)calloc(nnz, sizeof(INT)); // rowind is of length nnz
  mat->colind = (INT *)calloc(nnz, sizeof(INT)); // colind is of length nnz
  mat->val = (REAL *)calloc(nnz, sizeof(REAL)); // val is of length nnz
  for (i=0;i<mat->row;++i){
    for(ij=ia[i];ij<ia[i+1];++ij){
      // let us drop small entries
      //      if(fabs(A[ij])<1e-15) continue;
      mat->rowind[ij]=i;
      mat->colind[ij]=ja[ij];
      mat->val[ij]=A[ij];
    }
  }
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

ivector* create_ivector(int *x, int n)
{
  /* for now not copy arrays, make test that produce seg. fault  */
  ivector* vec;
  vec = (ivector* ) calloc(1, sizeof(ivector));
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
    REAL     setup_start, setup_end;
    pc->setup_time = 0.;
    get_time(&setup_start);

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
    get_time(&setup_end);
    pc->setup_time = setup_end - setup_start;

    return pc;
}

precond* create_precond_amg(dCSRmat *A, AMG_param *amgparam)
{

    precond *pc = (precond*)calloc(1, sizeof(precond));

    precond_data *pcdata = (precond_data*)calloc(1, sizeof(precond_data));

    const INT nnz = A->nnz, m = A->row, n = A->col;
    short prtlvl = amgparam->print_level;
    INT      status = SUCCESS;
    REAL     setup_start, setup_end;
    pc->setup_time = 0.;
    get_time(&setup_start);

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
    get_time(&setup_end);
    pc->setup_time = setup_end - setup_start;

    return pc;
}

precond* create_precond_famg(dCSRmat *A, dCSRmat *M, AMG_param *amgparam)
{

    precond *pc = (precond*)calloc(1, sizeof(precond));

    precond_data *pcdata = (precond_data*)calloc(1, sizeof(precond_data));

    const INT nnz = A->nnz, m = A->row, n = A->col, nnz_M = M->nnz;
    short prtlvl = amgparam->print_level;
    INT status = SUCCESS;
    REAL setup_start, setup_end;
    pc->setup_time = 0.;
    get_time(&setup_start);

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
    get_time(&setup_end);
    pc->setup_time = setup_end - setup_start;

    return pc;
}


precond* create_precond_amg_bsr(dBSRmat *A, AMG_param *amgparam)
{

    precond *pc = (precond*)calloc(1, sizeof(precond));

    precond_data_bsr *pcdata = (precond_data_bsr*)calloc(1, sizeof(precond_data_bsr));

    short prtlvl = amgparam->print_level;
    INT status = SUCCESS;
    REAL setup_start, setup_end;
    pc->setup_time = 0.;
    get_time(&setup_start);

    AMG_data_bsr *mgl = amg_data_bsr_create(amgparam->max_levels);
    mgl[0].A = dbsr_create(A->ROW, A->COL, A->NNZ, A->nb, A->storage_manner); dbsr_cp(A, &mgl[0].A);
    mgl[0].b = dvec_create(mgl[0].A.ROW*mgl[0].A.nb);
    mgl[0].x = dvec_create(mgl[0].A.COL*mgl[0].A.nb);

    // bsr amg setup
    switch (amgparam->AMG_type) {

        case UA_AMG: // Unmoothed Aggregation AMG
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG BSR...\n");
            status = amg_setup_ua_bsr(mgl, amgparam);
        break;

        default:
            if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG BSR...\n");
            status = amg_setup_ua_bsr(mgl, amgparam);
        break;

    }

    if (status < 0) {
        printf("\n In function create_precond: AMG BSR data failed to set up \n");
        return NULL;
    }

    // setup preconditioner
    // todo: make a function for this
    pcdata->print_level = amgparam->print_level;
    pcdata->maxit = amgparam->maxit;
    pcdata->tol = amgparam->tol;
    pcdata->cycle_type = amgparam->cycle_type;
    pcdata->smoother = amgparam->smoother;
    pcdata->presmooth_iter = amgparam->presmooth_iter;
    pcdata->postsmooth_iter = amgparam->postsmooth_iter;
    pcdata->relaxation = amgparam->relaxation;
    pcdata->coarse_scaling = amgparam->coarse_scaling;
    pcdata->amli_degree = amgparam->amli_degree;
    pcdata->amli_coef = amgparam->amli_coef;
    pcdata->tentative_smooth = amgparam->tentative_smooth;
    pcdata->max_levels = mgl[0].num_levels;
    pcdata->mgl_data = mgl;
    pcdata->A = A;

    pc->data = pcdata;

    switch (amgparam->cycle_type) {
        //case NL_AMLI_CYCLE: // Nonlinear AMLI AMG
        //    prec.fct = precond_dbsr_namli; break;
        default: // V,W-Cycle AMG
            pc->fct = precond_dbsr_amg;
        break;
    }
    get_time(&setup_end);
    pc->setup_time = setup_end - setup_start;

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
    INT i, j;
    REAL setup_start, setup_end;
    pc->setup_time = 0.;
    get_time(&setup_start);

    //------------------------------------------------
    // compute the rational approximation
    //------------------------------------------------
    // scaling parameters for rational approximation
    REAL scaled_alpha = alpha, scaled_beta = beta;
    // scale alpha = alpha*sa^(-s)*sm^(s-1)
    scaled_alpha = alpha*pow(scaling_a, -s_frac_power)*pow(scaling_m, s_frac_power-1.);
    // scale beta = beta*sa^(-t)*sm^(t-1)
    scaled_beta  = beta*pow(scaling_a, -t_frac_power)*pow(scaling_m, t_frac_power-1.);

    /* Get interpolation points and function values */
    // parameters used in the function
    REAL16 func_param[4];
    func_param[0] = (REAL16)s_frac_power;
    func_param[1] = (REAL16)t_frac_power;
    if (scaled_alpha > scaled_beta)
    {
      func_param[2] = 1.;
      func_param[3] = (REAL16)scaled_beta/scaled_alpha;
    }
    else
    {
      func_param[2] = (REAL16)scaled_alpha/scaled_beta;
      func_param[3] = 1.;
    }

    // get points and function values first
    INT numval = (1<<14)+1;  // initial number of points on the interval [x_min, x_max]
    REAL xmin_in = 0.e0, xmax_in = 1.e0;  // interval for x
    REAL16 **zf = set_f_values(frac_inv, func_param[0], func_param[1], func_param[2], func_param[3], &numval, \
                               xmin_in, xmax_in, 0);

    /* AAA algorithm for the rational approximation */
    // parameters used in the AAA algorithm
    INT mmax_in = 30;  // maximal final number of pole + 1
    REAL16 AAA_tol = powl(2e0,-40e0);  // tolerance of the AAA algorithm
    INT k = -22; // k is the number of nodes in the final interpolation after tolerance is achieved or mmax is reached.
    INT print_level = 0; // print level for AAA
    // output of the AAA algorithm.  It contains residues (Re + Im), poles (Re + Im), nodes, weights, function values
    REAL **rpnwf = malloc(7 * sizeof(REAL *));

    // compute the rational approximation using AAA algorithms
    REAL err_max=get_rpzwf(numval, zf[0], zf[1], rpnwf, &mmax_in, &k, AAA_tol, print_level);
    if(rpnwf == NULL) {
      fprintf(stderr,"\nUnsuccessful AAA computation of rational approximation\n");
      fflush(stderr);
      return 0;
    }
    printf(" HAZ ---- Rational approx error in interp points: %.16e\n", err_max);

    // assign poles and residues
    REAL drop_tol = AAA_tol;
    INT ii = 1; // skip first residue

    REAL *polesr = malloc((k-1) * sizeof(REAL));
    REAL *polesi = malloc((k-1) * sizeof(REAL));
    REAL *resr = malloc(k * sizeof(REAL));
    REAL *resi = malloc(k * sizeof(REAL));

    /* filter poles and residues smaller than some tolerance */
    // Note: free residual is always only real! also, it's always saved to preserve the numbering (N poles, N+1 residues)
    resi[0] = 0.;
    if(fabs(rpnwf[0][0]) < drop_tol) resr[0] = 0.; else resr[0] = rpnwf[0][0];

    for(i = 1; i < k; ++i) {
        if((fabs(rpnwf[0][i]) < drop_tol) && (fabs(rpnwf[1][i]) < drop_tol)) {
            // Case 1: remove poles and residues where abs(res) < tol
            fprintf(stderr,"\n%%%%%% *** HAZMATH WARNING*** Pole number reduced in function=%s \n", \
                    __FUNCTION__);
            fprintf(stdout,"%%%%%%  Removing pole[%d] = %.8e + %.8e i \t residue[%d] = %.8e + %.8e i\n", \
	                i-1, rpnwf[2][i-1], rpnwf[3][i-1], i, rpnwf[0][i], rpnwf[1][i]);
	        fflush(stdout);
        }
        else if((fabs(rpnwf[0][i]) > drop_tol) && (fabs(rpnwf[3][i-1]) < drop_tol)) {
            // Case 2: only real residues and poles (Note: if Im(pole) = 0, then Im(res) = 0.)
            resr[ii] = rpnwf[0][i]; resi[ii] = 0.; polesi[ii-1] = 0.;
            if(fabs(rpnwf[2][i-1]) < drop_tol) polesr[ii-1] = 0.; else polesr[ii-1] = rpnwf[2][i-1];
            ii++;
        }
        else {
            // Case 3: there is at least one pair of complex conjugate poles
            // Note: only save one pole per pair -- check if it is already saved: (not the best search ever but ehh)
            for(j = 0; j < ii-1; ++j) {
                if((fabs(polesr[j] - rpnwf[2][i-1]) < drop_tol) && (fabs(polesi[j] - rpnwf[3][i-1]) < drop_tol)) break;
            }
            // if we found it, skip it
            if(j < ii-1) {
                fprintf(stderr,"\n%%%%%% *** HAZMATH WARNING*** Pole number reduced in function=%s \n", \
                    __FUNCTION__);
                fprintf(stdout,"%%%%%%  Removing pole[%d] = %.8e + %.8e i \t residue[%d] = %.8e + %.8e i\n", \
                        i-1, rpnwf[2][i-1], rpnwf[3][i-1], i, rpnwf[0][i], rpnwf[1][i]);
                fflush(stdout);
            }
            else {
                polesi[ii-1] = rpnwf[3][i-1]; // this should be always > drop_tol in Case 3
                if(fabs(rpnwf[0][i]) > drop_tol) resr[ii] = rpnwf[0][i]; else resr[ii] = 0.;
                if(fabs(rpnwf[1][i]) > drop_tol) resi[ii] = rpnwf[1][i]; else resi[ii] = 0.;
                if(fabs(rpnwf[2][i-1]) > drop_tol) polesr[ii-1] = rpnwf[2][i-1]; else polesr[ii-1] = 0.;
                ii++;
            }
        }
    }
    // new number of poles+1
    k = ii;

    // the easiest way seems to first save real part and then imag part
    pcdata->residues = dvec_create_p(2*k);
    pcdata->poles = dvec_create_p(2*(k-1));

    // copy and append in pcdata
    array_cp(k, resr, pcdata->residues->val);
    array_cp(k, resi, &(pcdata->residues->val[k]));
    array_cp(k-1, polesr, pcdata->poles->val);
    array_cp(k-1, polesi, &(pcdata->poles->val[k-1]));

     // print poles, residues
    if(prtlvl > 1){
        printf("Poles:\n");
        for(i = 0; i < k-1; ++i) {
            printf("pole[%d] = %.10e + %.10e i\n", i, pcdata->poles->val[i], pcdata->poles->val[k-1+i]);
        }
        printf("\n");
        printf("Residues:\n");
        for(i = 0; i < k; ++i) {
            printf("res[%d] = %.10e + %.10e i\n", i, pcdata->residues->val[i], pcdata->residues->val[k+i]);
        }
        printf("\n");
    }
    // FIXME: is this necessary? maybe merge with previous or next loop
    /*REAL polez;
    for(i = 0; i < k-1; ++i) {
      polez = pcdata->poles->val[i];
      if(polez > 0e0){
        fprintf(stderr,"\n%%%%%% *** HAZMATH WARNING*** Positive pole in function=%s", \
            __FUNCTION__);
	    fprintf(stdout,"\n%%%%%%  0 < pole(%d)=%.16e\n", i, polez);
	    break;
      }
    }*/

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
    INT npoles = k - 1;
    pcdata->mgl = (AMG_data **)calloc(npoles, sizeof(AMG_data *));

    // assemble amg data for all shifted laplacians:
    // (scaling_a*A - poles[i] * scaling_m*M)
    // NOTE: ONLY REAL PART USED HERE.
    // Also, I didn't distinguish between zero and nonzero poles because
    // some pole can have only imag part nonzero (so it's still needed in the solve).
    for(i = 0; i < npoles; ++i) {
        //fprintf(stdout,"\nAMG for pole %d\n", i);
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
    // NOTE: all poles have the same amgparams
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
    get_time(&setup_end);
    pc->setup_time = setup_end - setup_start;


    // clean
    // if (rpnwf[0]) free(rpnwf); // FIXME
    free(resr); free(resi); free(polesr); free(polesi);
    if (rpnwf) free(rpnwf);

    return pc;
}

INT get_poles_no(precond* pc)
{
    precond_ra_data* data = (precond_ra_data*)(pc->data);

    return data->poles->row;
}

dvector* ra_aaa(int numval,
                double *z,
                double *f,
                double AAA_tol)
{
    INT i;
    //------------------------------------------------
    // compute the rational approximation
    //------------------------------------------------

    /*fprintf(stdout, "npoints: %d\n", numval); fflush(stdout);
    for(i = 0; i < numval-1; ++i) {
        fprintf(stdout, "%.10e \t %.10e\n", z[i], f[i]); fflush(stdout);
    }*/

    /* AAA algorithm for the rational approximation */
    // parameters used in the AAA algorithm
    INT mmax_in = 30;  // maximal final number of pole + 1
    INT k = -22; // k is the number of nodes in the final interpolation after tolerance is achieved or mmax is reached.
    INT print_level = 1; // print level for AAA

    // output of the AAA algorithm.  It contains residues (Re + Im), poles (Re + Im), nodes, weights, function values
    REAL **rpnwf = malloc(7 * sizeof(REAL *));

    // cast points and function values to long double
    REAL16 *zz = malloc(numval * sizeof(REAL16));
    REAL16 *ff = malloc(numval * sizeof(REAL16));

    for(i = 0; i < numval; ++i) {
        zz[i] = (REAL16)z[i];
        ff[i] = (REAL16)f[i];
    }

    /*fprintf(stdout, "npoints: %d\n", numval); fflush(stdout);
    for(i = 0; i < numval; ++i) {
        fprintf(stdout, "%Le \t %Le\n", zz[i], ff[i]); fflush(stdout);
    }*/

    // compute the rational approximation using AAA algorithms
    REAL err_max = get_rpzwf(numval, zz, ff, rpnwf, &mmax_in, &k, AAA_tol, print_level);
    if(rpnwf == NULL) {
      fprintf(stderr,"\nUnsuccessful AAA computation of rational approximation\n");
      fflush(stderr);
      return 0;
    }
    printf(" HAZ ---- Rational approx error in interp points: %.16e\n", err_max);
    // printf("Number of poles: %d\n", k-1);

    // print poles, residues
    /*printf("Poles:\n");
    for(i = 0; i < k-1; ++i) {
        printf("pole[%d] = %.10e + %.10e i\n", i, rpnwf[2][i], rpnwf[3][i]);
    }
    printf("\n");
    printf("Residues:\n");
    for(i = 0; i < k; ++i) {
        printf("res[%d] = %.10e + %.10e i\n", i, rpnwf[0][i], rpnwf[1][i]);
    }
    printf("\n");*/

    // assign poles and residuals
    dvector *res = dvec_create_p(4*k - 2);
    array_cp(k-1, rpnwf[2], res->val);
    array_cp(k-1, rpnwf[3], &(res->val[k-1]));
    array_cp(k, rpnwf[0], &(res->val[2*k-2]));
    array_cp(k, rpnwf[1], &(res->val[3*k-2]));

    // check if poles are non negative
    /*
    REAL polez;
    for(i = 0; i < k-1; ++i) {
      polez = res->val[i];
      if(polez > 0e0) {
	    fprintf(stderr,"\n%%%%%% *** HAZMATH WARNING*** Positive pole in function=%s", \
	        __FUNCTION__);
	    fprintf(stdout,"\n%%%%%%  0 < pole(%d)=%.16e\n", i, polez);
	    break;
      }
    }*/

    // clean
    if(rpnwf[0]) free(rpnwf[0]);
    if(rpnwf) free(rpnwf);
    free(zz); free(ff);

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
    REAL setup_start, setup_end;
    pc->setup_time = 0.;
    get_time(&setup_start);

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
    get_time(&setup_end);
    pc->setup_time = setup_end - setup_start;


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
    REAL setup_start, setup_end;
    pc->setup_time = 0.;
    get_time(&setup_start);

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
    get_time(&setup_end);
    pc->setup_time = setup_end - setup_start;

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
    REAL setup_start, setup_end;
    pc->setup_time = 0.;
    get_time(&setup_start);

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
    get_time(&setup_end);
    pc->setup_time = setup_end - setup_start;

    return pc;
}

void apply_precond(REAL *r, REAL *z, precond *pc)
{
    //printf("calling pc->fct \n");
    pc->fct(r, z, pc->data);
    //printf("done calling pc->fct \n");

}

/*
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
*/
/*
PyObject* py_callback_eval(REAL *r, REAL *x, smoother_matvec *smmv)
{
    printf("Here I am inside the callback eval\n");

    //import_array();
    // get data
    //LUDMIL    npy_intp D[1]; D[0] = smmv->data->A->col; t//
    // LUDMIL Not exactly sure what this (ABOVE) is supposed to do (find number of columns?), so I rewrote it
    // INT D[1];//LUDMIL
    // dCSRmat *A=(dCSRmat *)smmv->data; //LUDMIL
    // D[0] = A->col; //LUDMIL


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
*/

INT wrapper_krylov_amg(dCSRmat *mat, dvector *rhs, dvector *sol)
{
    INT niters = 0;

    /* set Parameters from Reading in Input File */
    input_param inparam;
    param_input_init(&inparam);
    param_input("./input.dat", &inparam);

    // Set parameters for linear iterative methods
    linear_itsolver_param itparam;  // parameters for linear itsolver
    param_linear_solver_init(&itparam);
    param_linear_solver_set(&itparam, &inparam);
    if (itparam.linear_print_level > PRINT_MIN) param_linear_solver_print(&itparam);

    // Set parameters for algebriac multigrid methods
    AMG_param amgparam; // parameters for AMG
    param_amg_init(&amgparam);
    param_amg_set(&amgparam, &inparam);
    if (amgparam.print_level > PRINT_MIN) param_amg_print(&amgparam);

    niters = linear_solver_dcsr_krylov_amg(mat, rhs, sol, &itparam, &amgparam);

    return niters;
}


INT wrapper_krylov_amg_schwarz(dCSRmat *mat, dvector *rhs, dvector *sol)
{
    INT niters = 0;

    /* set Parameters from Reading in Input File */
    input_param inparam;
    param_input_init(&inparam);
    param_input("./input_schwarz.dat", &inparam);

    /* Set parameters for linear iterative methods */
    linear_itsolver_param itparam;
    param_linear_solver_init(&itparam);
    param_linear_solver_set(&itparam, &inparam);
    if (itparam.linear_print_level > PRINT_MIN) param_linear_solver_print(&itparam);

    // Set parameters for algebraic multigrid methods
    AMG_param  amgparam; // parameters for AMG
    param_amg_init(&amgparam);
    param_amg_set(&amgparam, &inparam);
    if (amgparam.print_level > PRINT_MIN) param_amg_print(&amgparam);

    niters = linear_solver_dcsr_krylov_amg(mat, rhs, sol, &itparam, &amgparam);

    return niters;
}


INT fenics_bsr_solver(INT block_size, dCSRmat *A, dvector *b, dvector *sol)
{
    INT i;
    INT *perm = (INT*)calloc(2*block_size, sizeof(INT));
    for(i = 0; i < block_size; ++i)
    {
        perm[2*i] = i;
        perm[2*i+1] = i + block_size;
    }
    dCSRmat AT;
    dcsr_alloc(A->col, A->row, A->nnz, &AT);

    dcsr_transz(A, perm, &AT);
    dcsr_transz(&AT, perm, A);

    dBSRmat Absr = dcsr_2_dbsr(A, 2);
    dcsr_free(&AT);

    for(i = 0; i < 2*block_size; ++i)
    {
        sol->val[i] = b->val[perm[i]];
    }
    dvec_cp(sol, b);
    /* set Parameters from Reading in Input File */
    input_param inparam;
    param_input_init(&inparam);
    param_input("./input.dat", &inparam);

    /* set initial guess */
    dvec_set(b->row, sol, 0.0);

    /* Set Solver Parameters */
    INT solver_flag=-20;
    /* Set parameters for linear iterative methods */
    linear_itsolver_param linear_itparam;
    param_linear_solver_set(&linear_itparam, &inparam);
    if (linear_itparam.linear_print_level > PRINT_MIN) param_linear_solver_print(&linear_itparam);

    /* Set parameters for algebraic multigrid methods */
    AMG_param amgparam;
    param_amg_init(&amgparam);
    param_amg_set(&amgparam, &inparam);
    if (amgparam.print_level > PRINT_MIN) param_amg_print(&amgparam);

    printf("\n===========================================================================\n");
    printf("Solving the linear system \n");
    printf("===========================================================================\n");

    // Use Krylov Iterative Solver
    solver_flag = linear_solver_dbsr_krylov_amg(&Absr, b, sol, &linear_itparam, &amgparam);

    for(i = 0; i < 2*block_size; ++i)
    {
        b->val[perm[i]] = sol->val[i];
    }
    dvec_cp(b, sol);
    free(perm);

    return solver_flag;
}


INT fenics_metric_amg_solver(block_dCSRmat *A,
                             dvector *b,
                             dvector *x,
                             block_dCSRmat *AD,
                             block_dCSRmat *M,
                             dCSRmat *interface_dof)
{
    /* set Parameters from Reading in Input File */
    input_param inparam;
    param_input_init(&inparam);
    param_input("./input_metric.dat", &inparam);

    /* set initial guess */
    dvec_set(b->row, x, 0.0);

    // drop small entries
    dcsr_compress_inplace(AD->blocks[0], 1e-12);
    dcsr_compress_inplace(AD->blocks[3], 1e-12);
    dcsr_compress_inplace(M->blocks[0], 1e-12);
    //dcsr_compress_inplace(M->blocks[1], 1e-12);
    //dcsr_compress_inplace(M->blocks[2], 1e-12);
    dcsr_compress_inplace(M->blocks[3], 1e-12);

    /* Set Solver Parameters */
    INT solver_flag = -20;
    /* Set parameters for linear iterative methods */
    linear_itsolver_param linear_itparam;
    param_linear_solver_set(&linear_itparam, &inparam);
    if (linear_itparam.linear_print_level > PRINT_MIN) param_linear_solver_print(&linear_itparam);

    /* Set parameters for algebriac multigrid methods */
    AMG_param amgparam;
    param_amg_init(&amgparam);
    param_amg_set(&amgparam, &inparam);
    if(amgparam.print_level > PRINT_MIN) param_amg_print(&amgparam);

    printf("\n===========================================================================\n");
    printf("Solving the linear system \n");
    printf("===========================================================================\n");

    // Use Krylov Iterative Solver
    if ( (linear_itparam.linear_precond_type >= 10) && (linear_itparam.linear_precond_type < 15) ){
        // todo: interface_dof is not used in the function below - remove as requirement
        solver_flag = linear_solver_bdcsr_krylov_metric_amg(A, b, x, &linear_itparam, &amgparam, AD, M, interface_dof);
    }
    // No preconditioner
    else{
        solver_flag = linear_solver_bdcsr_krylov(A, b, x, &linear_itparam);
    }

    return solver_flag;
}


void print_bdcsr_matrix(block_dCSRmat *A)
{
    fprintf(stdout,"\n------------ A ---------- \n"); fflush(stdout);
    bdcsr_print_matlab(stdout, A);
    fflush(stdout);
}
