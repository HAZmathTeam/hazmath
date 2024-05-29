/*! \file src/solver/amg_setup_ua.c
 *
 *  Unsmoothed Aggregation AMG: SETUP phase
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 12/24/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note   Done cleanup for releasing -- Xiaozhe Hu 03/11/2017 & 08/27/2021
 *
 *  \todo   Add safe guard for the overall computatinal complexity -- Xiaozhe Hu
 *  \todo   Add maximal weighted matching coarsning -- Xiaozhe Hu
 *  \todo   Add maximal independent set aggregation -- Xiaozhe Hu
 *
 */

#include "hazmath.h"

static void form_tentative_p(ivector *vertices, dCSRmat *tentp, REAL **basis, INT levelNum, INT num_aggregations);
static void construct_strongly_coupled(dCSRmat *A, AMG_param *param, dCSRmat *Neigh);
static SHORT aggregation_hec(dCSRmat *A, ivector *vertices, AMG_param *param, dCSRmat *Neigh, INT *num_aggregations, INT lvl);
static INT heavy_edg(const REAL *wei,const INT *numb, const INT n0,const INT n1);
static SHORT aggregation_hem(dCSRmat *A, ivector *vertices, AMG_param *param, dCSRmat *Neigh, INT *num_aggregations, INT lvl);
static SHORT aggregation_vmb(dCSRmat *A, ivector *vertices, AMG_param *param, dCSRmat *Neigh, INT *num_aggregations, INT lvl);
static void smooth_aggregation_p(dCSRmat *A, dCSRmat *tentp, dCSRmat *P, AMG_param *param, INT levelNum, dCSRmat *N);
static SHORT amg_setup_unsmoothP_unsmoothR(AMG_data *, AMG_param *);
static SHORT amg_setup_smoothP_smoothR(AMG_data *, AMG_param *);
static SHORT famg_setup_unsmoothP_unsmoothR(AMG_data *, AMG_param *);
static SHORT famg_setup_smoothP_smoothR(AMG_data *, AMG_param *);
static void form_boolean_p_bsr(const ivector *vertices,dBSRmat *tentp,const AMG_data_bsr *mgl,const INT NumAggregates);
static void form_tentative_p_bsr(const ivector *vertices,dBSRmat *tentp, const AMG_data_bsr *mgl,const INT NumAggregates,const INT dim,REAL **basis);
static SHORT amg_setup_unsmoothP_unsmoothR_bsr(AMG_data_bsr *mgl, AMG_param *param);
static SHORT amg_setup_general_bdcsr(AMG_data_bdcsr *mgl, AMG_param *param);
static SHORT amg_setup_bdcsr_metric(AMG_data_bdcsr *mgl, AMG_param *param);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/***********************************************************************************************/
/**
 * \fn SHORT amg_setup_ua (AMG_data *mgl, AMG_param *param)
 *
 * \brief Set up phase of unsmoothed aggregation AMG
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       SUCCESS if successed; otherwise, error information.
 *
 * \author Xiaozhe Hu
 * \date   12/28/2011
 */
SHORT amg_setup_ua (AMG_data *mgl,
                    AMG_param *param)
{

    SHORT status = amg_setup_unsmoothP_unsmoothR(mgl, param);

    return status;
}

/***********************************************************************************************/
/**
 * \fn SHORT amg_setup_sa (AMG_data *mgl, AMG_param *param)
 *
 * \brief Set up phase of smoothed aggregation AMG
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       SUCCESS if successed; otherwise, error information.
 *
 * \author Xiaozhe Hu
 * \date   07/13/2020
 */
SHORT amg_setup_sa (AMG_data *mgl,
                    AMG_param *param)
{

    SHORT status = amg_setup_smoothP_smoothR(mgl, param);

    return status;
}

/***********************************************************************************************/
/**
 * \fn SHORT famg_setup_ua (AMG_data *mgl, AMG_param *param)
 *
 * \brief Set up phase of unsmoothed aggregation AMG with fractional smoothers
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       SUCCESS if successed; otherwise, error information.
 *
 * \author Ana Budisa
 * \date   2020-05-20
 */
SHORT famg_setup_ua (AMG_data *mgl,
                     AMG_param *param)
{

    SHORT status = famg_setup_unsmoothP_unsmoothR(mgl, param);

    return status;
}

/***********************************************************************************************/
/**
 * \fn SHORT famg_setup_sa (AMG_data *mgl, AMG_param *param)
 *
 * \brief Set up phase of smoothed aggregation AMG
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       SUCCESS if successed; otherwise, error information.
 *
 * \author Xiaozhe Hu
 * \date   07/13/2020
 */
SHORT famg_setup_sa (AMG_data *mgl,
                     AMG_param *param)
{

    SHORT status = famg_setup_smoothP_smoothR(mgl, param);

    return status;
}

/***********************************************************************************************/
/**
 * \fn INT amg_setup_ua_bsr (AMG_data_bdcsr *mgl, AMG_param *param)
 *
 * \brief Set up phase of unsmoothed aggregation AMG (BSR format)
 *
 * \param mgl    Pointer to AMG data: AMG_data_bcr
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       SUCCESS if successed; otherwise, error information.
 *
 * \author Xiaozhe Hu
 * \date   03/16/2012
 */
SHORT amg_setup_ua_bsr(AMG_data_bsr  *mgl,
                      AMG_param     *param)
{
#if DEBUG_MODE > 0
    printf("### DEBUG: [-Begin-] %s ...\n", __FUNCTION__);
#endif

    SHORT status = amg_setup_unsmoothP_unsmoothR_bsr(mgl, param);

#if DEBUG_MODE > 0
    printf("### DEBUG: [--End--] %s ...\n", __FUNCTION__);
#endif

    return status;
}


/***********************************************************************************************/
/**
 * \fn INT amg_setup_bdcsr (AMG_data_bdcsr *mgl, AMG_param *param)
 *
 * \brief Set up phase of AMG (block_dCSRmat format)
 *
 * \param mgl    Pointer to AMG data: AMG_data_bdcsr
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       SUCCESS if successed; otherwise, error information.
 *
 * \author Xiaozhe Hu
 * \date   04/16/2012
 */
SHORT amg_setup_bdcsr(AMG_data_bdcsr  *mgl,
                      AMG_param     *param)
{
#if DEBUG_MODE > 0
    printf("### DEBUG: [-Begin-] %s ...\n", __FUNCTION__);
#endif

    SHORT status = amg_setup_general_bdcsr(mgl, param);

#if DEBUG_MODE > 0
    printf("### DEBUG: [--End--] %s ...\n", __FUNCTION__);
#endif

    return status;
}

/***********************************************************************************************/
/**
 * \fn INT metric_amg_setup_bdcsr (AMG_data_bdcsr *mgl, AMG_param *param)
 *
 * \brief Set up phase of metric AMG (block_dCSRmat format)
 *
 * \param mgl    Pointer to AMG data: AMG_data_bdcsr
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       SUCCESS if successed; otherwise, error information.
 *
 * \author Xiaozhe Hu
 * \date   05/15/2012
 *
 * \note  special AMG setup phase for interface prolem only
 *
 */
SHORT metric_amg_setup_bdcsr(AMG_data_bdcsr  *mgl,
                            AMG_param     *param)
{
#if DEBUG_MODE > 0
    printf("### DEBUG: [-Begin-] %s ...\n", __FUNCTION__);
#endif

    SHORT status = amg_setup_bdcsr_metric(mgl, param);

#if DEBUG_MODE > 0
    printf("### DEBUG: [--End--] %s ...\n", __FUNCTION__);
#endif

    return status;
}


/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/
/***********************************************************************************************/
/**
 * \fn static void form_tentative_p (ivector *vertices, dCSRmat *tentp,
 *                                   REAL **basis, INT levelNum, INT num_aggregations)
 *
 * \brief Form aggregation based on strong coupled neighbors
 *
 * \param vertices           Pointer to the aggregation of vertices
 * \param tentp              Pointer to the prolongation operators
 * \param basis              Pointer to the near kernel
 * \param levelNum           Level number
 * \param num_aggregations   Number of aggregations
 *
 * \author Xiaozhe Hu
 * \date   09/29/2009
 *
 * \note Modified by Xiaozhe Hu on 05/25/2014
 *
 * \note this subroutine uses given near-null basis to form tentative prolongation
 *       and the basis is not necessary constant
 *
 */
static void form_tentative_p(ivector *vertices,
                             dCSRmat *tentp,
                             REAL **basis,
                             INT levelNum,
                             INT num_aggregations)
{
    INT i, j;

    /* Form tentative prolongation */
    tentp->row = vertices->row;
    tentp->col = num_aggregations;
    tentp->nnz = vertices->row;

    tentp->IA  = (INT *)calloc(tentp->row+1,sizeof(INT));

    // local variables
    INT  *IA = tentp->IA;
    INT  *vval = vertices->val;
    const INT row = tentp->row;

    // first run
    for ( i = 0, j = 0; i < row; i++ ) {
        IA[i] = j;
        if (vval[i] > UNPT) j++;
    }
    IA[row] = j;

    // allocate memory for P
    tentp->nnz = j;
    tentp->JA  = (INT *)calloc(tentp->nnz, sizeof(INT));
    tentp->val = (REAL *)calloc(tentp->nnz, sizeof(REAL));

    INT  *JA = tentp->JA;
    REAL *val = tentp->val;

    // second run
    for (i = 0, j = 0; i < row; i ++) {
        IA[i] = j;
        if (vval[i] > UNPT) {
            JA[j] = vval[i];
            val[j] = basis[0][i];
            j ++;
        }
    }
}
/*COMMENTED AS IT IS NOT USED*/
/* /\***********************************************************************************************\/ */
/* /\** */
/*  * \fn static void form_boolean_p (ivector *vertices, dCSRmat *tentp, INT levelNum, */
/*  *                                 INT num_aggregations) */
/*  * */
/*  * \brief Form aggregation based on strong coupled neighbors */
/*  * */
/*  * \param vertices           Pointer to the aggregation of vertices */
/*  * \param tentp              Pointer to the prolongation operators */
/*  * \param levelNum           Level number */
/*  * \param num_aggregations   Number of aggregations */
/*  * */
/*  * \author Xiaozhe Hu */
/*  * \date   09/29/2009 */
/*  * */
/*  * Modified by Xiaozhe Hu on 05/25/2014 */
/*  *\/ */
/* static void form_boolean_p(ivector *vertices, */
/*                            dCSRmat *tentp, */
/*                            INT levelNum, */
/*                            INT num_aggregations) */
/* { */
/*     INT i, j; */

/*     /\* Form tentative prolongation *\/ */
/*     tentp->row = vertices->row; */
/*     tentp->col = num_aggregations; */
/*     tentp->nnz = vertices->row; */

/*     tentp->IA  = (INT *)calloc(tentp->row+1,sizeof(INT)); */

/*     // local variables */
/*     INT  *IA = tentp->IA; */
/*     INT  *vval = vertices->val; */
/*     const INT row = tentp->row; */

/*     // first run */
/*     for ( i = 0, j = 0; i < row; i++ ) { */
/*         IA[i] = j; */
/*         if (vval[i] > UNPT) j++; */
/*     } */
/*     IA[row] = j; */

/*     // allocate memory for P */
/*     tentp->nnz = j; */
/*     tentp->JA  = (INT *)calloc(tentp->nnz, sizeof(INT)); */
/*     tentp->val = (REAL *)calloc(tentp->nnz, sizeof(REAL)); */

/*     INT  *JA = tentp->JA; */
/*     REAL *val = tentp->val; */

/*     // second run */
/*     for (i = 0, j = 0; i < row; i ++) { */
/*         IA[i] = j; */
/*         if (vval[i] > UNPT) { */
/*             JA[j] = vval[i]; */
/*             val[j] = 1.0; */
/*             j ++; */
/*         } */
/*     } */
/* } */


/***********************************************************************************************/
/**
 * \fn static construct_strongly_coupled (dCSRmat *A, AMG_param *param, dCSRmat *Neigh)
 *
 * \brief get strongly coupled matrix by dropping relative small entries
 *
 * \param A           Pointer to the matrix A
 * \param param       Pointer to AMG parameters (contain the parameter determines strong coulping)
 * \param Neigh       Pointer to the matrix which only contains the strongly coupled neighbors
 *
 * \author Xiaozhe Hu
 * \date   09/29/2009
 *
 * Modified by Xiaozhe Hu on 08/27/2021
 */
static void construct_strongly_coupled(dCSRmat *A,
                                       AMG_param *param,
                                       dCSRmat *Neigh)
{

  // local variables
  const INT  row = A->row, col = A->col, nnz = A->IA[row]-A->IA[0];
  const INT  *AIA = A->IA, *AJA = A->JA;
  const REAL *Aval = A->val;

  INT  i, j, index, row_start, row_end;
  REAL strongly_coupled = param->strong_coupled;
  if(0){
    if(fabs(strongly_coupled) < 1e-6) strongly_coupled=1e-6;
  }
  // REAL strongly_coupled2 = pow(strongly_coupled,2);
  REAL strongly_coupled2 = strongly_coupled*strongly_coupled;
  //
  INT  *NIA, *NJA;
  REAL *Nval;

  // get the diagonal entries
  dvector diag;
  dcsr_getdiag(0, A, &diag);

  // allocate Neigh
  dcsr_alloc(row, col, nnz, Neigh);

  NIA  = Neigh->IA; NJA  = Neigh->JA;
  Nval = Neigh->val;

  // set IA for Neigh
  for ( i = row; i >= 0; i-- ) NIA[i] = AIA[i];

  // main loop of finding strongly coupled neighbors
  for ( index = i = 0; i < row; ++i ) {
    NIA[i] = index;
    row_start = AIA[i]; row_end = AIA[i+1];
    //fprintf(stdout,"\nHere coupled index %lld \n", (long long )index);
    for ( j = row_start; j < row_end; ++j ) {
      //fprintf(stdout,"\nHere coupled row start %lld end %lld \n", (long long )row_start, (long long )row_end);
      if ( (AJA[j] == i)
	   || ( ((Aval[j]*Aval[j]) >= strongly_coupled2*ABS(diag.val[i]*diag.val[AJA[j]])) && (Aval[j] < 0e0) )
	   )
	{
	  //fprintf(stdout,"\nHere coupled \n");
	  NJA[index] = AJA[j];
	  Nval[index] = Aval[j];
	  index++;
	  //fprintf(stdout,"\nHere coupled \n");
	}

    } // end for ( j = row_start; j < row_end; ++j )
  } // end for ( index = i = 0; i < row; ++i )

  dvec_free(&diag); // free it here;
  NIA[row] = index;

  Neigh->nnz = index;
  Neigh->JA  = (INT*) realloc(Neigh->JA,  (Neigh->IA[row])*sizeof(INT));
  Neigh->val = (REAL*)realloc(Neigh->val, (Neigh->IA[row])*sizeof(REAL));
  //
  if(0){
    //begin finding connected components (ltz):
    iCSRmat *blk_dfs=run_dfs(Neigh->row,Neigh->IA, Neigh->JA);
    index=0;
    for(i=0;i<blk_dfs->row;++i){
      j=blk_dfs->IA[i+1]-blk_dfs->IA[i];
      if(j>1){
	/* fprintf(stdout,"\nnontrivial block:size(%lld)=%d",(long long )i,(long long )j); */
	index++;
      }
    }
    fprintf(stdout,"\n blocks(total)=%lld ; blocks(non-trivial:size>1)=%lld; strongly_coupled=%.5e\n",(long long )blk_dfs->row,(long long )index,strongly_coupled);
    icsr_free(blk_dfs);free(blk_dfs);
    //end finding connected components (ltz):
  } //end if(0);
}

/***********************************************************************************************/
/**
 * \fn static void smooth_aggregation_p(dCSRmat *A, dCSRmat *tentp, dCSRmat *P,
 *                                      AMG_param *param, INT levelNum, dCSRmat *N)
 *
 * \brief Smooth the tentative prolongation
 *
 * \param A         Pointer to the coefficient matrices
 * \param tentp     Pointer to the tentative prolongation operators
 * \param P         Pointer to the prolongation operators
 * \param param     Pointer to AMG parameters
 * \param levelNum  Current level number
 * \param N         Pointer to strongly coupled neighbors
 *
 * \author Xiaozhe Hu
 * \date   09/29/2009
 *
 */
static void smooth_aggregation_p(dCSRmat *A,
                                dCSRmat *tentp,
                                dCSRmat *P,
                                AMG_param *param,
                                INT levelNum,
                                dCSRmat *N)
{
    const SHORT filter = param->smooth_filter;
    const INT   row = A->row, col= A->col;
    const REAL  smooth_factor = param->tentative_smooth;

    dCSRmat S;
    dvector diag;  // diagonal entries

    REAL row_sum_A, row_sum_N;
    INT i,j;

    /* Step 1. Form smoother */
    /* Without filter: Using A for damped Jacobian smoother */
    if ( filter != ON ) {

        // copy structure from A
        S = dcsr_create(row, col, A->IA[row]);

        for ( i=0; i<=row; ++i ) S.IA[i] = A->IA[i];
        for ( i=0; i<S.IA[S.row]; ++i ) S.JA[i] = A->JA[i];

        dcsr_getdiag(0, A, &diag);  // get the diagonal entries of A

        // check the diagonal entries.
        // if it is too small, use Richardson smoother for the corresponding row
        for (i=0; i<row; ++i) {
            if (ABS(diag.val[i]) < 1e-6) diag.val[i] = 1.0;
        }

        for (i=0; i<row; ++i) {
            for (j=S.IA[i]; j<S.IA[i+1]; ++j) {
                if (S.JA[j] == i) {
                    S.val[j] = 1 - smooth_factor * A->val[j] / diag.val[i];
                }
                else {
                    S.val[j] = - smooth_factor * A->val[j] / diag.val[i];
                }
            }
        }
    }

    /* Using filtered A for damped Jacobian smoother */
    else {
        /* Form filtered A and store in N */
        for (i=0; i<row; ++i) {
          for (row_sum_A = 0.0, j=A->IA[i]; j<A->IA[i+1]; ++j) {
            if (A->JA[j] != i) row_sum_A += A->val[j];
          }

          for (row_sum_N = 0.0, j=N->IA[i]; j<N->IA[i+1]; ++j) {
            if (N->JA[j] != i) row_sum_N += N->val[j];
          }

          for (j=N->IA[i]; j<N->IA[i+1]; ++j) {
            if (N->JA[j] == i) {
              N->val[j] += row_sum_A - row_sum_N;
            }
          }
        }
        // copy structure from N (filtered A)
        S = dcsr_create(row, col, N->IA[row]);

        for (i=0; i<=row; ++i) S.IA[i] = N->IA[i];

        for (i=0; i<S.IA[S.row]; ++i) S.JA[i] = N->JA[i];

        dcsr_getdiag(0, N, &diag);  // get the diagonal entries of N (filtered A)

        // check the diagonal entries.
        // if it is too small, use Richardson smoother for the corresponding row
        for (i=0;i<row;++i) {
            if (ABS(diag.val[i]) < 1e-6) diag.val[i] = 1.0;
        }

        for (i=0;i<row;++i) {
            for (j=S.IA[i]; j<S.IA[i+1]; ++j) {
                if (S.JA[j] == i) {
                    S.val[j] = 1 - smooth_factor * N->val[j] / diag.val[i];
                }
                else {
                    S.val[j] = - smooth_factor * N->val[j] / diag.val[i];
                }
            }
        }

    }

    dvec_free(&diag);

    /* Step 2. Smooth the tentative prolongation P = S*tenp */
    dcsr_mxm(&S, tentp, P);
    P->nnz = P->IA[P->row];
    dcsr_free(&S);

}

/***********************************************************************************************/
/**
 * \fn static SHORT aggregation_vmb (dCSRmat *A, ivector *vertices, AMG_param *param,
 *                                   dCSRmat *Neigh, INT *num_aggregations, INT lvl)
 *
 * \brief Form aggregation based on strongly coupled neighbors using greedy method
 *
 * \param A                 Pointer to the coefficient matrices
 * \param vertices          Pointer to the aggregation of vertices
 * \param param             Pointer to AMG parameters
 * \param Neigh             Pointer to strongly coupled neighbors
 * \param num_aggregations  Pointer to number of aggregations
 * \param lvl               Level number
 *
 * \author Xiaozhe Hu
 * \date   09/29/2009
 *
 * \note Setup A, P, PT and levels using the unsmoothed aggregation algorithm;
 *       Reference: P. Vanek, J. Mandel and M. Brezina
 *       "Algebraic Multigrid on Unstructured Meshes", 1994
 *
 */
static SHORT aggregation_vmb(dCSRmat *A,
                             ivector *vertices,
                             AMG_param *param,
                             dCSRmat *Neigh,
                             INT *num_aggregations, INT lvl)
{
    // local variables
    const INT    row = A->row;
    const INT    max_aggregation = param->max_aggregation;

    // return status
    SHORT  status = SUCCESS;

    // local variables
    INT    num_left = row;
    INT    subset, count;
    INT    *num_each_agg;

    INT    i, j, row_start, row_end;
    INT    *NIA = NULL, *NJA = NULL;

    // find strongly coupled neighbors
    construct_strongly_coupled(A, param, Neigh);

    NIA  = Neigh->IA; NJA  = Neigh->JA;

    /*------------------------------------------*/
    /*             Initialization               */
    /*------------------------------------------*/
    ivec_alloc(row, vertices);
    iarray_set(row, vertices->val, -2);
    *num_aggregations = 0;

    /*-------------*/
    /*   Step 1.   */
    /*-------------*/
    //    for ( i = 0; i < row; ++i ) {
    for ( i=row-1; i>=0;i-- ) {
        if ( (NIA[i+1] - NIA[i]) == 1 ) {
            vertices->val[i] = UNPT;
            num_left--;
        }
        else {
            subset = TRUE;
            row_start = NIA[i]; row_end = NIA[i+1];
            for ( j = row_start; j < row_end; ++j ) {
                if ( vertices->val[NJA[j]] >= UNPT ) {
                    subset = FALSE;
                    break;
                }
            }
            if ( subset ) {
                count = 0;
                vertices->val[i] = *num_aggregations;
                num_left--;
                count++;
                row_start = NIA[i]; row_end = NIA[i+1];
                for ( j = row_start; j < row_end; ++j ) {
                    if ( (NJA[j]!=i) && (count < max_aggregation) ) {
                        vertices->val[NJA[j]] = *num_aggregations;
                        num_left--;
                        count ++;
                    }
                }
                (*num_aggregations)++;
            }
        }
    }

    /*-------------*/
    /*   Step 2.   */
    /*-------------*/
    INT *temp_C = (INT*)calloc(row,sizeof(INT));

    if ( *num_aggregations < MIN_CDOF ) {
        status = ERROR_AMG_COARSEING; goto END;
    }

    num_each_agg = (INT*)calloc(*num_aggregations,sizeof(INT));

    for ( i = row; i--; ) {
        temp_C[i] = vertices->val[i];
        if ( vertices->val[i] >= 0 ) num_each_agg[vertices->val[i]] ++;
    }

    //    for ( i = 0; i < row; ++i ) {
    for ( i = row-1; i >= 0; i-- ) {
        if ( vertices->val[i] < UNPT ) {
            row_start = NIA[i]; row_end = NIA[i+1];
            for ( j = row_start; j < row_end; ++j ) {
                if ( temp_C[NJA[j]] > UNPT
                    && num_each_agg[temp_C[NJA[j]]] < max_aggregation ) {
                    vertices->val[i] = temp_C[NJA[j]];
                    num_left--;
                    num_each_agg[temp_C[NJA[j]]] ++ ;
                    break;
                }
            }
        }
    }

    /*-------------*/
    /*   Step 3.   */
    /*-------------*/
    while ( num_left > 0 ) {
        for ( i = 0; i < row; ++i ) {
            if ( vertices->val[i] < UNPT ) {
                count = 0;
                vertices->val[i] = *num_aggregations;
                num_left--;
                count++;
                row_start = NIA[i]; row_end = NIA[i+1];
                for ( j = row_start; j < row_end; ++j ) {
                    if ( (NJA[j]!=i) && (vertices->val[NJA[j]] < UNPT)
						             && (count<max_aggregation) ) {
                        vertices->val[NJA[j]] = *num_aggregations;
                        num_left--;
                        count++;
                    }
                }
                (*num_aggregations)++;
            }
        }
    }

    free(num_each_agg);

END:
    free(temp_C);

    return status;
}

/***********************************************************************************************/
/**
 * \fn static SHORT aggregation_hec (dCSRmat *A, ivector *vertices, AMG_param *param,
 *                                   dCSRmat *Neigh, INT *num_aggregations,INT lvl)
 *
 * \brief Heavy edge coarsening aggregation based on strongly coupled neighbors
 *
 * \param A                 Pointer to the coefficient matrices
 * \param vertices          Pointer to the aggregation of vertices
 * \param param             Pointer to AMG parameters
 * \param Neigh             Pointer to strongly coupled neighbors
 * \param num_aggregations  Pointer to number of aggregations
 * \param lvl               Level number
 *
 * \author Xiaozhe Hu
 * \date   03/12/2017
 *
 * \TODO Add control of maximal size of each aggregates
 *
 * \note Refer to J. Urschel, X. Hu, J. Xu and L. Zikatanov
 *       "A Cascadic Multigrid Algorithm for Computing the Fiedler Vector of Graph Laplacians", 2015
 *
 */
static SHORT aggregation_hec(dCSRmat *A,
                             ivector *vertices,
                             AMG_param *param,
                             dCSRmat *Neigh,
                             INT *num_aggregations, INT lvl)
{
    // local variables
    const INT    row = A->row;
    //const INT    max_aggregation = param->max_aggregation;

    // return status
    SHORT  status = SUCCESS;

    INT  i, j, k, ii, jj, row_start, row_end;;
    INT  *NIA = NULL, *NJA = NULL;
    REAL *Nval = NULL;
    REAL maxval = 0.0;

    INT *perm =(INT *)calloc(row, sizeof(INT));

    // find strongly coupled neighbors
    construct_strongly_coupled(A, param, Neigh);

    NIA  = Neigh->IA; NJA  = Neigh->JA;
    Nval = Neigh->val;

    /*------------------------------------------*/
    /*             Initialization               */
    /*------------------------------------------*/
    ivec_alloc(row, vertices);
    iarray_set(row, vertices->val, -2);
    *num_aggregations = 0;

    // get random permutation
    for(i=0; i<row; i++) perm[i] = i;
    iarray_shuffle(row, perm);

    // main loop
    for ( ii = 0; ii < row; ii++ ) {
      i = perm[ii];
      if ( (NIA[i+1] - NIA[i]) == 1 ) {
	vertices->val[i] = UNPT;
      } else {
	// find the most strongly connected neighbor
	row_start = NIA[i]; row_end = NIA[i+1];
	maxval = 0.0;
	k = -1;
	for (jj = row_start; jj < row_end; jj++) {
	  if (NJA[jj] != i) {
	    if ( ABS(Nval[jj]) > maxval ) {
	      k = jj;
	      maxval = ABS(Nval[jj]);
	    }
	  }
	} // end for (jj = row_start+1; jj < row_end; jj++)
	if(k<0) {
	  // run again to find anything in this row:
	  maxval = -1e0;
	  for (jj = row_start; jj < row_end; jj++){
	    if (NJA[jj] != i) {
	      if ( ABS(Nval[jj]) > maxval ) {
		k = jj;
		maxval = ABS(Nval[jj]);
	      }
	    }
	      } // end for (jj = row_start+1; jj < row_end; jj++)
	}
	if(k<0) {
	  // well we should never come here, but if we do:
	  fprintf(stderr,"\n\n%%%% ERROR Stop: negative index encountered in %s()\n\n",__FUNCTION__);
	  exit(16);
	}
	j = NJA[k];
	// create a new aggregates if the most strongly conncected neighbor
	// is still avaliable
	if (vertices->val[j] < UNPT) {
	    vertices->val[j] = *num_aggregations;
	    vertices->val[i] = *num_aggregations;
                (*num_aggregations)++;
            }
            else
            {
                vertices->val[i] = vertices->val[j];
            } // end if (vertices->val[j] < UNPT)
        } // end if ( (NIA[i+1] - NIA[i]) == 1 )

    } // end for ( ii = 0; ii < row; ii++ )

    // free spaces
    free(perm);

    return status;

}
/***********************************************************************/
/* \fn static void heavy_edg(REAL *wei,INT *numb,INT *iheavy,INT n)
 *
 * \brief [iheavy]=argmax(wei(k),k in numb(1:n))
 *
 * \param wei               real array(n) with weights
 * \param numb              integer array(n) with indices
 * \param n                 size of wei and n;
 *
 * \return *iheavy           numb[k] where wei[k] is maximal, k=1:n.
 *
 * \author Ludmil Zikatanov
 * \date   20230328
 *
 * \note Refer to Kim, Xu, Zikatanov 2003: "A multigrid method based on graph matching for
 *                                          convection–di␁usion equations"
 *
 */
static INT heavy_edg(const REAL *wei,const INT *numb, const INT n0, const INT n1)
{
  /*====================================================================*/
  /*--------------------------------------------------------------------
  ...  Pick the heaviest WEI.
  --------------------------------------------------------------------*/
  INT j,nend,step=1;
  if(n0<0) return -1;
  if(n1<0) return -1;
  if(n0>n1) {step=-1;}
  nend=n1+step;
  /* while (numb[j]<0){ */
  /*   iheavy = numb[j]; */
  /*   temp = wei[j]; */
  /*   j+=step; */
  /* } */
  INT iheavy = numb[n0];
  REAL temp = wei[n0];
  j=n0;
  while(j!=nend){
    if(wei[j] > temp && numb[j] >=0){
      temp = wei[j];
      iheavy = numb[j];
    }
    j+=step;
  }
  return iheavy;
}
/**************************************************************************************/
/* \fn static SHORT aggregation_hem (dCSRmat *A, ivector *vertices, AMG_param *param,
 *                                   dCSRmat *Neigh, INT *num_aggregations,INT lvl)
 *
 * \brief heavy/light edge matching based on strongly
 *        coupled neighbors
 *
 * \param A                 Pointer to the coefficient matrices
 * \param vertices          Pointer to the aggregation of vertices
 * \param param             Pointer to AMG parameters
 * \param Neigh             Pointer to strongly coupled neighbors
 * \param num_aggregations  Pointer to number of aggregations
 * \param lvl               Level number
 *
 * \author Ludmil Zikatanov
 * \date   20230328
 *
 * \note Refer to Kim, Xu, Zikatanov 2003: "A multigrid method based on graph matching for
 *                                          convection–di␁usion equations"
 *
 */
static SHORT aggregation_hem(dCSRmat *A,
			     ivector *vertices,
			     AMG_param *param,
			     dCSRmat *Neigh,
			     INT *num_aggregations, INT lvl)
{
    // local variables
  const INT n= A->row;
  SHORT  status = SUCCESS;
  SHORT print_level=(SHORT )param->print_level;
  //
  INT  j, k, jk,iz,pick,row_start, row_end;
  REAL ajk;
  //ORDER the column indices in Neigh in increasing order; just in case (so bidomain
  //examples behave in a certain way); these work without such ordering too, but...
  ///////////////////////////////////////////////////
  // dCSRmat AT=dcsr_create(0,0,0); // ORDERING
  //construct_strongly_coupled(A, param, &AT);// ORDERING
  //  dcsr_alloc(AT.col,AT.row,AT.nnz,Neigh); // ORDERING
  //  dcsr_transz(&AT,NULL,Neigh); // ORDERING
  //  dcsr_free(&AT); // ORDERING
  ///////////////////////////////////////////////////
  construct_strongly_coupled(A, param, Neigh);// //No ORDERING of column indices
  /******************************/
  //
  INT *ia  = Neigh->IA,*ja  = Neigh->JA;
  REAL *a = Neigh->val;
  //INT *ia  =A->IA,*ja  =A->JA;
  //REAL *a = A->val;
  /*------------------------------------------*/
  /*             Initialization               */
  /*------------------------------------------*/
  ivec_alloc(n, vertices);
  iarray_set(n, vertices->val, -1);
  INT *mask=vertices->val;
  ivector num_els=ivec_create(n);
  iarray_set(n, num_els.val, 0);
  /*find first the "diagonal part" of A */
  INT maxdeg=0,l=0; //INT kdir = 0;
  for(k=0;k<n;++k){
    l=ia[k+1]-ia[k];
    if(l>maxdeg) maxdeg=l;
    if(l>1) continue;
    mask[k] = -2;
    //    kdir++;
  }
  REAL *work=calloc(maxdeg,sizeof(REAL));
  INT *iwork=calloc(maxdeg,sizeof(INT));
  //
  INT kmatch=-1;
  kmatch=kmatch+1;  // stupid way to make the complier happy (and James happy)... --Xiaozhe
  //  INT kiso=0;
  INT nc=0;
  for(k=0;k<n;++k){
    if((mask[k]<0) && (mask[k]>-2)){
      // this is interior and unmatched; count its unmatched neighbors:
      iz = 0;
      row_start=ia[k]; row_end=ia[k+1];
      for(jk = row_start;jk<row_end;++jk){
	j = ja[jk];
	ajk = a[jk];
	if(mask[j]<0 && j!=k){
	  // found an unmatched neighbor;
	  work[iz]=fabs(ajk);
	  iwork[iz]=j;
	  iz++;
	}
      }
      if(iz){
	//     Matched edges.
	pick=heavy_edg(work,iwork,iz-1,0);
	mask[k] = nc;
	mask[pick] = nc;
	num_els.val[nc]+=2;
	kmatch+=2;// two are matched
	// these are the only nc, so increment only here!
	nc++;	
      }
      /* else{ */
      /* 	//Isolated points; no change in mask (vertices->val)! */
      /* 	// num_els.val[nc]++; */
      /* 	INT kiso++;// one isolated; */
      /* } */
    }
  }
  if(print_level>10){
    fprintf(stdout,"\n%%%%num(aggregates(pass1))=%lld\n", (long long )nc); fflush(stdout);
  }
  num_els.row=nc;
  num_els.val=realloc(num_els.val,num_els.row*sizeof(INT));
  INT kc;
  /**/
  /*second run to remove the isolated*/
  for(k=0;k<n;++k){
    if(mask[k]!=(-1)) continue;
      // this is isolated
    iz=0;
    row_start=ia[k];
    row_end=ia[k+1];
    for(jk = row_start;jk<row_end;++jk){
      j = ja[jk];
      //Not needed ...      ajk = a[jk];
      kc=mask[j];
      if(kc>=0){// here by default we cannot have k=j
	//found an aggregate nearby
	work[iz]=(REAL )(-num_els.val[kc]);
	iwork[iz]=kc;
	iz++;
      }
    }
    if(iz){
      //add to the aggregate kc which has least number of elements
      kc=heavy_edg(work,iwork,0,iz-1);
      mask[k] =kc;
      kmatch++;
      num_els.val[kc]++;
    }else{
      //Isolated points again, this cannot happen, so it must be an error
      fprintf(stderr,"%%%%ERROR in %s: isolated point (=%d) on the second matching pass",__FUNCTION__,k);
      exit(17);
    }
  }
  //    fprintf(stdout,"\n%%%%After Pass2: nc=%d\n\n",nc); fflush(stdout);
  free(work);
  free(iwork);
  ivec_free(&num_els);
  //
  if(print_level>10){
    fprintf(stdout,"\n%%%%num(aggregates(pass2))=%lld (should be the same as pass1)\n", (long long )nc); fflush(stdout);
  }
  *num_aggregations = nc;
  for(k=0;k<n;++k){
    if(mask[k]<0){
      mask[k]=UNPT;
      //      fprintf(stdout,"\nvec[%d]=%d",k,mask[k]);
    }
  }
  return status;
}
/**
 * \fn static void form_boolean_p_bsr(const ivector *vertices, dBSRmat *tentp,
 *                                    const AMG_data_bsr *mgl,
 *                                    const INT NumAggregates)
 *
 * \brief Form boolean prolongations in dBSRmat (assume constant vector is in
 *        the null space)
 *
 * \param vertices           Pointer to the aggregation of vertices
 * \param tentp              Pointer to the prolongation operators
 * \param mgl                Pointer to AMG levels
 * \param NumAggregates      Number of aggregations
 *
 * \author Xiaozhe Hu
 * \date   05/27/2014
 */
static void form_boolean_p_bsr(const ivector       *vertices,
                               dBSRmat             *tentp,
                               const AMG_data_bsr  *mgl,
                               const INT            NumAggregates)
{
    INT i, j;

    /* Form tentative prolongation */
    tentp->ROW = vertices->row;
    tentp->COL = NumAggregates;
    tentp->nb  = mgl->A.nb;
    INT nb2    = tentp->nb * tentp->nb;

    tentp->IA  = (INT*)calloc(tentp->ROW+1, sizeof(INT));

    // local variables
    INT * IA = tentp->IA;
    INT *JA;
    REAL *val;
    INT *vval = vertices->val;

    const INT row = tentp->ROW;

    // first run
    for (i = 0, j = 0; i < row; i ++) {
        IA[i] = j;
        if (vval[i] > -1) {
            j ++;
        }
    }
    IA[row] = j;

    // allocate
    tentp->NNZ = j;

    tentp->JA = (INT*)calloc(tentp->NNZ, sizeof(INT));

    tentp->val = (REAL*)calloc(tentp->NNZ*nb2, sizeof(REAL));

    JA = tentp->JA;
    val = tentp->val;

    // second run
    for (i = 0, j = 0; i < row; i ++) {
        IA[i] = j;
        if (vval[i] > -1) {
            JA[j] = vval[i];
            ddense_identity (&(val[j*nb2]), tentp->nb, nb2);
            j ++;
        }
    }
}

/**
 * \fn static void form_tentative_p_bsr(const ivector *vertices, dBSRmat *tentp,
 *                                      const AMG_data_bsr *mgl, const INT NumAggregates,
 *                                      const const INT dim, REAL **basis)
 *
 * \brief Form tentative prolongation for BSR format matrix (use general basis for
 *        the null space)
 *
 * \param vertices           Pointer to the aggregation of vertices
 * \param tentp              Pointer to the prolongation operators
 * \param mgl                Pointer to AMG levels
 * \param NumAggregates      Number of aggregations
 * \param dim                Dimension of the near kernel space
 * \param basis              Pointer to the basis of the near kernel space
 *
 * \author Xiaozhe Hu
 * \date   05/27/2014
 */
static void form_tentative_p_bsr(const ivector       *vertices,
                                 dBSRmat             *tentp,
                                 const AMG_data_bsr  *mgl,
                                 const INT            NumAggregates,
                                 const INT            dim,
                                 REAL               **basis)
{
    INT i, j, k;

    INT p, q;

    const INT nnz_row = dim/mgl->A.nb; // nonzeros per row

    /* Form tentative prolongation */
    tentp->ROW = vertices->row;
    tentp->COL = NumAggregates*nnz_row;
    tentp->nb = mgl->A.nb;
    const INT nb = tentp->nb;
    const INT nb2 = nb * nb;

    tentp->IA  = (INT*)calloc(tentp->ROW+1, sizeof(INT));

    // local variables
    INT  *IA = tentp->IA;
    INT  *JA;
    REAL *val;

    const INT *vval = vertices->val;
    const INT  row = tentp->ROW;

    // first run
    for (i = 0, j = 0; i < row; i ++) {
        IA[i] = j;
        if (vval[i] > -1) {
            j = j + nnz_row;
        }
    }
    IA[row] = j;

    // allocate
    tentp->NNZ = j;
    tentp->JA = (INT*)calloc(tentp->NNZ, sizeof(INT));
    tentp->val = (REAL*)calloc(tentp->NNZ*nb2, sizeof(REAL));

    JA  = tentp->JA;
    val = tentp->val;

    // second run
    for (i = 0, j = 0; i < row; i ++) {
        IA[i] = j;
        if (vval[i] > -1) {

            for (k=0; k<nnz_row; k++) {

                JA[j] = vval[i]*nnz_row + k;

                for (p=0; p<nb; p++) {

                    for (q=0; q<nb; q++) {

                        val[j*nb2 + p*nb + q] = basis[k*nb+p][i*nb+q];

                    }

                }

                j++;

            }
        }
    }
}


/***********************************************************************************************/
/**
 * \fn static SHORT amg_setup_unsmoothP_unsmoothR (AMG_data *mgl, AMG_param *param)
 *
 * \brief Setup phase of plain aggregation AMG, using unsmoothed P and unsmoothed A
 *
 * \param mgl    Pointer to AMG_data
 * \param param  Pointer to AMG_param
 *
 * \return       SUCCESS if succeed, error otherwise
 *
 * \author Xiaozhe Hu
 * \date   02/21/2011
 *
 */
static SHORT amg_setup_unsmoothP_unsmoothR(AMG_data *mgl,
                                           AMG_param *param)
{
    // local variables
    const SHORT prtlvl     = param->print_level;
    const SHORT cycle_type = param->cycle_type;
    const SHORT csolver    = param->coarse_solver;
    const SHORT min_cdof   = MAX(param->coarse_dof,1);
    const INT   m          = mgl[0].A.row;

    // local variables
    SHORT         max_levels = param->max_levels, lvl = 0, status = SUCCESS;
    INT           i;
    REAL          setup_start, setup_end;
    Schwarz_param swzparam;
    //dCSRmat As;

    get_time(&setup_start);

    // level info (fine: 0; coarse: 1)
    ivector *vertices = (ivector *)calloc(max_levels,sizeof(ivector));

    // each level stores the information of the number of aggregations
    INT *num_aggs = (INT *)calloc(max_levels,sizeof(INT));

    // each level stores the information of the strongly coupled neighborhoods
    dCSRmat *Neighbor = (dCSRmat *)calloc(max_levels,sizeof(dCSRmat));

    // Initialize level information
    for ( i = 0; i < max_levels; ++i ) num_aggs[i] = 0;

    mgl[0].near_kernel_dim   = 1;
    mgl[0].near_kernel_basis = (REAL **)calloc(mgl->near_kernel_dim,sizeof(REAL*));

    for ( i = 0; i < mgl->near_kernel_dim; ++i ) {
        mgl[0].near_kernel_basis[i] = (REAL *)calloc(m,sizeof(REAL));
        array_set(m, mgl[0].near_kernel_basis[i], 1.0);
    }

    // Initialize Schwarz parameters
    mgl->Schwarz_levels = param->Schwarz_levels;
    if ( param->Schwarz_levels > 0 ) {
        swzparam.Schwarz_mmsize = param->Schwarz_mmsize;
        swzparam.Schwarz_maxlvl = param->Schwarz_maxlvl;
        swzparam.Schwarz_type   = param->Schwarz_type;
        swzparam.Schwarz_blksolver = param->Schwarz_blksolver;
    }

    // Initialize AMLI coefficients
    if ( cycle_type == AMLI_CYCLE ) {
        const INT amlideg = param->amli_degree;
        param->amli_coef = (REAL *)calloc(amlideg+1,sizeof(REAL));
        REAL lambda_max = 2e0, lambda_min = 0.25*lambda_max;
        amg_amli_coef(lambda_max, lambda_min, amlideg, param->amli_coef);
    }

/* #if DIAGONAL_PREF */
/*     dcsr_diagpref(&mgl[0].A); // reorder each row to make diagonal appear first */
/* #endif */

    /*----------------------------*/
    /*--- checking aggregation ---*/
    /*----------------------------*/
    // Main AMG setup loop
    while ( (mgl[lvl].A.row > min_cdof) && (lvl < max_levels-1) ) {

      /*-- Setup Schwarz smoother if necessary */
      if ( lvl < param->Schwarz_levels ) {
          mgl[lvl].Schwarz.A = dcsr_sympat(&mgl[lvl].A);
          Schwarz_setup(&mgl[lvl].Schwarz, &swzparam,NULL);
      }

        /*-- Aggregation --*/
        switch ( param->aggregation_type ) {

            case VMB: // VMB aggregation
                status = aggregation_vmb(&mgl[lvl].A, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl], lvl);
                break;

            case HEC: // Heavy edge coarsening aggregation
                status = aggregation_hec(&mgl[lvl].A, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl], lvl);
                break;

            case HEM: // Heavy edge matching
                status = aggregation_hem(&mgl[lvl].A, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl], lvl);
                break;

            default: // wrong aggregation type
                status = ERROR_AMG_AGG_TYPE;
                check_error(status, __FUNCTION__);
                break;
        }

        /*-- Choose strength threshold adaptively --*/
        if ( num_aggs[lvl]*4 > mgl[lvl].A.row )
            param->strong_coupled /= 2;
        else if ( num_aggs[lvl]*1.25 < mgl[lvl].A.row )
            param->strong_coupled *= 2;

        // Check 1: Did coarsening step succeed?
        if ( status < 0 ) {
            // When error happens, stop at the current multigrid level!
            if ( prtlvl > PRINT_MIN ) {
                printf("### HAZMATH WARNING: Forming aggregates on level-%lld failed!\n", (long long )lvl);
            }
            status = SUCCESS; break;
        }

        /*-- Form Prolongation --*/
        form_tentative_p(&vertices[lvl], &mgl[lvl].P, mgl[0].near_kernel_basis,
                         lvl+1, num_aggs[lvl]);

        // Check 2: Is coarse sparse too small?
        if ( mgl[lvl].P.col < MIN_CDOF ) break;

#if 0
        // Check 3: Does this coarsening step too aggressive?
        if ( mgl[lvl].P.row > mgl[lvl].P.col * MAX_CRATE ) {
            if ( prtlvl > PRINT_MIN ) {
                printf("### HAZMATH WARNING: Coarsening might be too aggressive!\n");
                printf("### HAZMATH: Fine level = %lld, coarse level = %(long long )d. Discard!\n",
                       mgl[lvl].P.row, mgl[lvl].P.col);
            }
            break;
        }
#endif

#if 0
        // Check 4: Is this coarsening ratio too small?
        if ( (REAL)mgl[lvl].P.col > mgl[lvl].P.row * MIN_CRATE ) {
            if ( prtlvl > PRINT_MIN ) {
                printf("### WARNING: Coarsening rate is too small!\n");
                printf("### WARNING: Fine level = %lld, coarse level = %lld. Discard!\n",
                       (long long )mgl[lvl].P.row, (long long )mgl[lvl].P.col);
            }
            break;
        }
#endif

        /*-- Form restriction --*/
        dcsr_trans(&mgl[lvl].P, &mgl[lvl].R);

        /*-- Form coarse level stiffness matrix --*/
        dcsr_rap_agg(&mgl[lvl].R, &mgl[lvl].A, &mgl[lvl].P,
                               &mgl[lvl+1].A);

        dcsr_free(&Neighbor[lvl]);
        ivec_free(&vertices[lvl]);

        ++lvl;

    }

/* #if DIAGONAL_PREF */
/*     dcsr_diagpref(&mgl[lvl].A); // reorder each row to make diagonal appear first */
/* #endif */

    // Setup coarse level systems for direct solvers
    switch (csolver) {

      //#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
            // Need to sort the matrix A for UMFPACK to work
	  dCSRmat A_tran=dcsr_create(mgl[lvl].A.col,	\
				      mgl[lvl].A.row,	\
				      mgl[lvl].A.nnz);
            dcsr_transz(&mgl[lvl].A, NULL, &A_tran);
            dcsr_cp(&A_tran, &mgl[lvl].A);
            dcsr_free(&A_tran);
            mgl[lvl].Numeric = hazmath_factorize(&mgl[lvl].A, 0);
            break;
        }
	  //#endif
        default:
            // Do nothing!
            break;
    }

    // setup total level number and current level
    mgl[0].num_levels = max_levels = lvl+1;
    mgl[0].w          = dvec_create(m);

    for ( lvl = 1; lvl < max_levels; ++lvl) {
        INT mm = mgl[lvl].A.row;
        mgl[lvl].num_levels = max_levels;
        mgl[lvl].b          = dvec_create(mm);
        mgl[lvl].x          = dvec_create(mm);

        mgl[lvl].cycle_type     = cycle_type; // initialize cycle type!

        if ( cycle_type == NL_AMLI_CYCLE )
            mgl[lvl].w = dvec_create(3*mm);
        else
            mgl[lvl].w = dvec_create(2*mm);
    }

    if ( prtlvl > PRINT_NONE ) {
        get_time(&setup_end);
        print_amg_complexity(mgl,prtlvl);
        print_cputime("Unsmoothed aggregation setup", setup_end - setup_start);
    }

    for (i=0; i<max_levels; i++) {
      dcsr_free(&Neighbor[i]);
      ivec_free(&vertices[i]);
    }

    free(Neighbor);
    free(vertices);
    free(num_aggs);

    return status;
}


/***********************************************************************************************/
/**
 * \fn static SHORT amg_setup_smoothP_smoothR (AMG_data *mgl, AMG_param *param)
 *
 * \brief Setup phase of smoothed aggregation AMG, using smoothed P and smoothed A
 *
 * \param mgl    Pointer to AMG_data
 * \param param  Pointer to AMG_param
 *
 * \return       SUCCESS if succeed, error otherwise
 *
 * \author Xiaozhe Hu
 * \date   07/13/2020
 *
 */
static SHORT amg_setup_smoothP_smoothR(AMG_data *mgl,
                                       AMG_param *param)
{
    // local variables
    const SHORT prtlvl     = param->print_level;
    const SHORT cycle_type = param->cycle_type;
    const SHORT csolver    = param->coarse_solver;
    const SHORT min_cdof   = MAX(param->coarse_dof,1);
    const INT   m          = mgl[0].A.row;
    //INT counter = 0;

    // local variables
    SHORT         max_levels = param->max_levels, lvl = 0, status = SUCCESS;
    INT           i;
    REAL          setup_start, setup_end;
    Schwarz_param swzparam;

    get_time(&setup_start);

    // level info (fine: 0; coarse: 1)
    ivector *vertices = (ivector *)calloc(max_levels,sizeof(ivector));

    // each level stores the information of the number of aggregations
    INT *num_aggs = (INT *)calloc(max_levels,sizeof(INT));

    // each level stores the information of the strongly coupled neighborhoods
    dCSRmat *Neighbor = (dCSRmat *)calloc(max_levels,sizeof(dCSRmat));

    // each level stores the information of the tentative prolongations
    dCSRmat *tentative_p = (dCSRmat *)calloc(max_levels,sizeof(dCSRmat));

    // Initialize level information
    for ( i = 0; i < max_levels; ++i ) num_aggs[i] = 0;

    mgl[0].near_kernel_dim   = 1;
    mgl[0].near_kernel_basis = (REAL **)calloc(mgl->near_kernel_dim,sizeof(REAL*));

    //fprintf(stdout,"Here %d\n", counter); ++counter;

    for ( i = 0; i < mgl->near_kernel_dim; ++i ) {
        mgl[0].near_kernel_basis[i] = (REAL *)calloc(m,sizeof(REAL));
        array_set(m, mgl[0].near_kernel_basis[i], 1.0);
    }

    // Initialize Schwarz parameters
    mgl->Schwarz_levels = param->Schwarz_levels;
    if ( param->Schwarz_levels > 0 ) {
        swzparam.Schwarz_mmsize = param->Schwarz_mmsize;
        swzparam.Schwarz_maxlvl = param->Schwarz_maxlvl;
        swzparam.Schwarz_type   = param->Schwarz_type;
        swzparam.Schwarz_blksolver = param->Schwarz_blksolver;
    }

    // Initialize AMLI coefficients
    if ( cycle_type == AMLI_CYCLE ) {
        const INT amlideg = param->amli_degree;
        param->amli_coef = (REAL *)calloc(amlideg+1,sizeof(REAL));
        REAL lambda_max = 2.0, lambda_min = lambda_max/4;
        amg_amli_coef(lambda_max, lambda_min, amlideg, param->amli_coef);
    }

/* #if DIAGONAL_PREF */
/*     dcsr_diagpref(&mgl[0].A); // reorder each row to make diagonal appear first */
/* #endif */

    /*----------------------------*/
    /*--- checking aggregation ---*/
    /*----------------------------*/
    // Main AMG setup loop
    while ( (mgl[lvl].A.row > min_cdof) && (lvl < max_levels-1) ) {

        /*-- Setup Schwarz smoother if necessary */
        if ( lvl < param->Schwarz_levels ) {
          mgl[lvl].Schwarz.A=dcsr_sympat(&mgl[lvl].A);
          Schwarz_setup(&mgl[lvl].Schwarz, &swzparam,NULL);
          //printf("Schwarz setup done!\n");
        }

        /*-- Aggregation --*/
        switch ( param->aggregation_type ) {

            case VMB: // VMB aggregation
                status = aggregation_vmb(&mgl[lvl].A, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl],lvl);
                break;

            case HEC: // Heavy edge coarsening aggregation
                status = aggregation_hec(&mgl[lvl].A, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl],lvl);
                break;

            case HEM: // Heavy edge matching
                status = aggregation_hem(&mgl[lvl].A, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl],lvl);
                break;

            default: // wrong aggregation type
                status = ERROR_AMG_AGG_TYPE;
                check_error(status, __FUNCTION__);
                break;
        }

        /*-- Choose strength threshold adaptively --*/
        if ( num_aggs[lvl]*4 > mgl[lvl].A.row )
            param->strong_coupled /= 2;
        else if ( num_aggs[lvl]*1.25 < mgl[lvl].A.row )
            param->strong_coupled *= 2;

        // Check 1: Did coarsening step succeed?
        if ( status < 0 ) {
            // When error happens, stop at the current multigrid level!
            if ( prtlvl > PRINT_MIN ) {
                printf("### HAZMATH WARNING: Forming aggregates on level-%lld failed!\n", (long long )lvl);
            }
            status = SUCCESS; break;
        }

        /*-- Form Prolongation --*/
        form_tentative_p(&vertices[lvl], &tentative_p[lvl], mgl[0].near_kernel_basis, lvl+1, num_aggs[lvl]);

        /* -- Form smoothed prolongation -- */
        smooth_aggregation_p(&mgl[lvl].A, &tentative_p[lvl], &mgl[lvl].P, param, lvl+1, &Neighbor[lvl]);

        // Check 2: Is coarse sparse too small?
        if ( mgl[lvl].P.col < MIN_CDOF ) break;

#if 0
        // Check 3: Does this coarsening step too aggressive?
        if ( mgl[lvl].P.row > mgl[lvl].P.col * MAX_CRATE ) {
            if ( prtlvl > PRINT_MIN ) {
                printf("### HAZMATH WARNING: Coarsening might be too aggressive!\n");
                printf("### HAZMATH: Fine level = %lld, coarse level = %lld. Discard!\n",
                       (long long )mgl[lvl].P.row, (long long )mgl[lvl].P.col);
            }
            break;
        }
#endif

#if 0
        // Check 4: Is this coarsening ratio too small?
        if ( (REAL)mgl[lvl].P.col > mgl[lvl].P.row * MIN_CRATE ) {
            if ( prtlvl > PRINT_MIN ) {
                printf("### WARNING: Coarsening rate is too small!\n");
                printf("### WARNING: Fine level = %lld, coarse level = %lld. Discard!\n",
                       (long long )mgl[lvl].P.row, (long long )mgl[lvl].P.col);
            }
            break;
        }
#endif

        /*-- Form restriction --*/
        dcsr_trans(&mgl[lvl].P, &mgl[lvl].R);

        /*-- Form coarse level stiffness matrix --*/
        dcsr_rap(&mgl[lvl].R, &mgl[lvl].A, &mgl[lvl].P,
                               &mgl[lvl+1].A);

        dcsr_free(&Neighbor[lvl]);
        ivec_free(&vertices[lvl]);

        ++lvl;

    }

/* #if DIAGONAL_PREF */
/*     dcsr_diagpref(&mgl[lvl].A); // reorder each row to make diagonal appear first */
/* #endif */

    // Setup coarse level systems for direct solvers
    switch (csolver) {

      //#if WITH_SUITESPARSE
    case SOLVER_UMFPACK: {
      // Need to sort the matrix A for UMFPACK to work
      dCSRmat Ac_tran;
      dcsr_trans(&mgl[lvl].A, &Ac_tran);
      // It is equivalent to do transpose and then sort
      //     fasp_dcsr_trans(&mgl[lvl].A, &Ac_tran);
      //     fasp_dcsr_sort(&Ac_tran);
      dcsr_cp(&Ac_tran, &mgl[lvl].A);
      dcsr_free(&Ac_tran);
      mgl[lvl].Numeric = hazmath_factorize(&mgl[lvl].A, 0);
      break;
    }
	    //#endif
        default:
            // Do nothing!
            break;
    }

    // setup total level number and current level
    mgl[0].num_levels = max_levels = lvl+1;
    mgl[0].w          = dvec_create(m);

    for ( lvl = 1; lvl < max_levels; ++lvl) {
        INT mm = mgl[lvl].A.row;
        mgl[lvl].num_levels = max_levels;
        mgl[lvl].b          = dvec_create(mm);
        mgl[lvl].x          = dvec_create(mm);

        mgl[lvl].cycle_type     = cycle_type; // initialize cycle type!

        if ( cycle_type == NL_AMLI_CYCLE )
            mgl[lvl].w = dvec_create(3*mm);
        else
            mgl[lvl].w = dvec_create(2*mm);
    }

    if ( prtlvl > PRINT_NONE ) {
        get_time(&setup_end);
        print_amg_complexity(mgl,prtlvl);
        print_cputime("Smoothed aggregation setup", setup_end - setup_start);
    }

    for (i=0; i<max_levels; i++) {
      dcsr_free(&Neighbor[i]);
      ivec_free(&vertices[i]);
      dcsr_free(&tentative_p[i]);
    }

    free(Neighbor);
    free(vertices);
    free(num_aggs);
    free(tentative_p);

    return status;
}


/***********************************************************************************************/
/**
 * \fn static SHORT famg_setup_unsmoothP_unsmoothR (AMG_data *mgl, AMG_param *param)
 *
 * \brief Setup phase of plain aggregation AMG with fractional smoothers,
 *        using unsmoothed P and unsmoothed R
 *
 * \param mgl    Pointer to AMG_data
 * \param param  Pointer to AMG_param
 *
 * \return       SUCCESS if succeed, error otherwise
 *
 * \author Ana Budisa
 * \date   2020-05-20
 *
 */
static SHORT famg_setup_unsmoothP_unsmoothR(AMG_data *mgl,
                                            AMG_param *param)
{
    // local variables
    const SHORT prtlvl     = param->print_level;
    const SHORT cycle_type = param->cycle_type;
    const SHORT csolver    = param->coarse_solver;
    const SHORT min_cdof   = MAX(param->coarse_dof,50);
    const INT   m          = mgl[0].A.row;

    // local variables
    SHORT         max_levels = param->max_levels, lvl = 0, status = SUCCESS;
    INT           i;
    REAL          setup_start, setup_end;
    Schwarz_param swzparam;

    get_time(&setup_start);

    // level info (fine: 0; coarse: 1)
    ivector *vertices = (ivector *)calloc(max_levels,sizeof(ivector));

    // each level stores the information of the number of aggregations
    INT *num_aggs = (INT *)calloc(max_levels,sizeof(INT));

    // each level stores the information of the strongly coupled neighborhoods
    dCSRmat *Neighbor = (dCSRmat *)calloc(max_levels,sizeof(dCSRmat));

    // Initialize level information
    for ( i = 0; i < max_levels; ++i ) num_aggs[i] = 0;

    mgl[0].near_kernel_dim   = 1;
    mgl[0].near_kernel_basis = (REAL **)calloc(mgl->near_kernel_dim,sizeof(REAL*));

    for ( i = 0; i < mgl->near_kernel_dim; ++i ) {
        mgl[0].near_kernel_basis[i] = (REAL *)calloc(m,sizeof(REAL));
        array_set(m, mgl[0].near_kernel_basis[i], 1.0);
    }

    // Initialize Schwarz parameters
    mgl->Schwarz_levels = param->Schwarz_levels;
    if ( param->Schwarz_levels > 0 ) {
        swzparam.Schwarz_mmsize = param->Schwarz_mmsize;
        swzparam.Schwarz_maxlvl = param->Schwarz_maxlvl;
        swzparam.Schwarz_type   = param->Schwarz_type;
        swzparam.Schwarz_blksolver = param->Schwarz_blksolver;
    }

    // Initialize AMLI coefficients
    if ( cycle_type == AMLI_CYCLE ) {
        const INT amlideg = param->amli_degree;
        param->amli_coef = (REAL *)calloc(amlideg+1,sizeof(REAL));
        REAL lambda_max = 2.0, lambda_min = lambda_max/4;
        amg_amli_coef(lambda_max, lambda_min, amlideg, param->amli_coef);
    }

/* #if DIAGONAL_PREF */
/*     dcsr_diagpref(&mgl[0].A); // reorder each row to make diagonal appear first */
/* #endif */

    /*----------------------------*/
    /*--- checking aggregation ---*/
    /*----------------------------*/
    // Main AMG setup loop
    while ( (mgl[lvl].A.row > min_cdof) && (lvl < max_levels-1) ) {

      /*-- Setup Schwarz smoother if necessary */
      if ( lvl < param->Schwarz_levels ) {
          mgl[lvl].Schwarz.A=dcsr_sympat(&mgl[lvl].A);
          dcsr_shift(&(mgl[lvl].Schwarz.A), 1);
          Schwarz_setup(&mgl[lvl].Schwarz, &swzparam,NULL);
      }

        /*-- Aggregation --*/
        switch ( param->aggregation_type ) {

            case VMB: // VMB aggregation
                status = aggregation_vmb(&mgl[lvl].A, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl],lvl);
                break;

            case HEC: // Heavy edge coarsening aggregation
                status = aggregation_hec(&mgl[lvl].A, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl],lvl);
                break;

            case HEM: // Heavy edge matching
                status = aggregation_hem(&mgl[lvl].A, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl],lvl);
                break;

            default: // wrong aggregation type
                status = ERROR_AMG_AGG_TYPE;
                check_error(status, __FUNCTION__);
                break;
        }

        /*-- Choose strength threshold adaptively --*/
        if ( num_aggs[lvl]*4 > mgl[lvl].A.row )
            param->strong_coupled /= 2;
        else if ( num_aggs[lvl]*1.25 < mgl[lvl].A.row )
            param->strong_coupled *= 2;

        // Check 1: Did coarsening step succeed?
        if ( status < 0 ) {
            // When error happens, stop at the current multigrid level!
            if ( prtlvl > PRINT_MIN ) {
                printf("### HAZMATH WARNING: Forming aggregates on level-%lld failed!\n", (long long )lvl);
            }
            status = SUCCESS; break;
        }

        /*-- Form Prolongation --*/
        form_tentative_p(&vertices[lvl], &mgl[lvl].P, mgl[0].near_kernel_basis,
                         lvl+1, num_aggs[lvl]);

        // Check 2: Is coarse sparse too small?
        if ( mgl[lvl].P.col < MIN_CDOF ) break;

#if 0
        // Check 3: Does this coarsening step too aggressive?
        if ( mgl[lvl].P.row > mgl[lvl].P.col * MAX_CRATE ) {
            if ( prtlvl > PRINT_MIN ) {
                printf("### HAZMATH WARNING: Coarsening might be too aggressive!\n");
                printf("### HAZMATH: Fine level = %lld, coarse level = %lld. Discard!\n",
                       (long long )mgl[lvl].P.row, (long long )mgl[lvl].P.col);
            }
            break;
        }
#endif

#if 0
        // Check 4: Is this coarsening ratio too small?
        if ( (REAL)mgl[lvl].P.col > mgl[lvl].P.row * MIN_CRATE ) {
            if ( prtlvl > PRINT_MIN ) {
                printf("### WARNING: Coarsening rate is too small!\n");
                printf("### WARNING: Fine level = %lld, coarse level = %lld. Discard!\n",
                       (long long )mgl[lvl].P.row, (long long )mgl[lvl].P.col);
            }
            break;
        }
#endif

        /*-- Form restriction --*/
        dcsr_trans(&mgl[lvl].P, &mgl[lvl].R);

        /*-- Form coarse level stiffness matrix --*/
        dcsr_rap_agg(&mgl[lvl].R, &mgl[lvl].A, &mgl[lvl].P,
                               &mgl[lvl+1].A);

        // added 2020-05-13 by Ana Budisa
        /*-- Form coarse level mass matrix --*/
        dcsr_rap_agg(&mgl[lvl].R, &mgl[lvl].M, &mgl[lvl].P, &mgl[lvl+1].M);

        dcsr_free(&Neighbor[lvl]);
        ivec_free(&vertices[lvl]);

        ++lvl;

    }

/* #if DIAGONAL_PREF */
/*     dcsr_diagpref(&mgl[lvl].A); // reorder each row to make diagonal appear first */
/* #endif */

    // Setup coarse level systems for direct solvers
    switch (csolver) {

      //#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
            // Need to sort the matrix A for UMFPACK to work
            dCSRmat Ac_tran;
            dcsr_trans(&mgl[lvl].A, &Ac_tran);
            // It is equivalent to do transpose and then sort
            //     fasp_dcsr_trans(&mgl[lvl].A, &Ac_tran);
            //     fasp_dcsr_sort(&Ac_tran);
            dcsr_cp(&Ac_tran, &mgl[lvl].A);
            dcsr_free(&Ac_tran);
            mgl[lvl].Numeric = hazmath_factorize(&mgl[lvl].A, 0);
            break;
        }
	  //#endif
        default:
            // Do nothing!
            break;
    }

    // setup total level number and current level
    mgl[0].num_levels = max_levels = lvl+1;
    mgl[0].w          = dvec_create(m);

    for ( lvl = 1; lvl < max_levels; ++lvl) {
        INT mm = mgl[lvl].A.row;
        mgl[lvl].num_levels = max_levels;
        mgl[lvl].b          = dvec_create(mm);
        mgl[lvl].x          = dvec_create(mm);

        mgl[lvl].cycle_type     = cycle_type; // initialize cycle type!

        if ( cycle_type == NL_AMLI_CYCLE )
            mgl[lvl].w = dvec_create(3*mm);
        else
            mgl[lvl].w = dvec_create(2*mm);
    }

    if ( prtlvl > PRINT_NONE ) {
        get_time(&setup_end);
        print_amg_complexity(mgl,prtlvl);
        print_cputime("Unsmoothed aggregation setup", setup_end - setup_start);
    }

    free(Neighbor);
    free(vertices);
    free(num_aggs);

    return status;
}


/***********************************************************************************************/
/**
 * \fn static SHORT famg_setup_smoothP_smoothR (AMG_data *mgl, AMG_param *param)
 *
 * \brief Setup phase of smoothed aggregation AMG, using smoothed P and smoothed A
 *
 * \param mgl    Pointer to AMG_data
 * \param param  Pointer to AMG_param
 *
 * \return       SUCCESS if succeed, error otherwise
 *
 * \author Xiaozhe Hu, Ana Budisa
 * \date   07/16/2020
 *
 */
static SHORT famg_setup_smoothP_smoothR(AMG_data *mgl,
                                       AMG_param *param)
{
    // local variables
    const SHORT prtlvl     = param->print_level;
    const SHORT cycle_type = param->cycle_type;
    const SHORT csolver    = param->coarse_solver;
    const SHORT min_cdof   = MAX(param->coarse_dof,1);
    const INT   m          = mgl[0].A.row;

    // local variables
    SHORT         max_levels = param->max_levels, lvl = 0, status = SUCCESS;
    INT           i;
    REAL          setup_start, setup_end;
    Schwarz_param swzparam;

    get_time(&setup_start);

    // level info (fine: 0; coarse: 1)
    ivector *vertices = (ivector *)calloc(max_levels,sizeof(ivector));

    // each level stores the information of the number of aggregations
    INT *num_aggs = (INT *)calloc(max_levels,sizeof(INT));

    // each level stores the information of the strongly coupled neighborhoods
    dCSRmat *Neighbor = (dCSRmat *)calloc(max_levels,sizeof(dCSRmat));

    // each level stores the information of the tentative prolongations
    dCSRmat *tentative_p = (dCSRmat *)calloc(max_levels,sizeof(dCSRmat));

    // Initialize level information
    for ( i = 0; i < max_levels; ++i ) num_aggs[i] = 0;

    mgl[0].near_kernel_dim   = 1;
    mgl[0].near_kernel_basis = (REAL **)calloc(mgl->near_kernel_dim,sizeof(REAL*));

    for ( i = 0; i < mgl->near_kernel_dim; ++i ) {
        mgl[0].near_kernel_basis[i] = (REAL *)calloc(m,sizeof(REAL));
        array_set(m, mgl[0].near_kernel_basis[i], 1.0);
    }

    // Initialize Schwarz parameters
    mgl->Schwarz_levels = param->Schwarz_levels;
    if ( param->Schwarz_levels > 0 ) {
        swzparam.Schwarz_mmsize = param->Schwarz_mmsize;
        swzparam.Schwarz_maxlvl = param->Schwarz_maxlvl;
        swzparam.Schwarz_type   = param->Schwarz_type;
        swzparam.Schwarz_blksolver = param->Schwarz_blksolver;
    }

    // Initialize AMLI coefficients
    if ( cycle_type == AMLI_CYCLE ) {
        const INT amlideg = param->amli_degree;
        param->amli_coef = (REAL *)calloc(amlideg+1,sizeof(REAL));
        REAL lambda_max = 2.0, lambda_min = lambda_max/4;
        amg_amli_coef(lambda_max, lambda_min, amlideg, param->amli_coef);
    }

/* #if DIAGONAL_PREF */
/*     dcsr_diagpref(&mgl[0].A); // reorder each row to make diagonal appear first */
/* #endif */

    /*----------------------------*/
    /*--- checking aggregation ---*/
    /*----------------------------*/
    // Main AMG setup loop
    while ( (mgl[lvl].A.row > min_cdof) && (lvl < max_levels-1) ) {

      /*-- Setup Schwarz smoother if necessary */
      if ( lvl < param->Schwarz_levels ) {
          mgl[lvl].Schwarz.A=dcsr_sympat(&mgl[lvl].A);
          Schwarz_setup(&mgl[lvl].Schwarz, &swzparam,NULL);
      }

        /*-- Aggregation --*/
        switch ( param->aggregation_type ) {

            case VMB: // VMB aggregation
                status = aggregation_vmb(&mgl[lvl].A, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl],lvl);
                break;

            case HEC: // Heavy edge coarsening aggregation
                status = aggregation_hec(&mgl[lvl].A, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl],lvl);
                break;

            case HEM: // Heavy edge matching
                status = aggregation_hem(&mgl[lvl].A, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl],lvl);
                break;

            default: // wrong aggregation type
                status = ERROR_AMG_AGG_TYPE;
                check_error(status, __FUNCTION__);
                break;
        }

        /*-- Choose strength threshold adaptively --*/
        if ( num_aggs[lvl]*4 > mgl[lvl].A.row )
            param->strong_coupled /= 2;
        else if ( num_aggs[lvl]*1.25 < mgl[lvl].A.row )
            param->strong_coupled *= 2;

        // Check 1: Did coarsening step succeed?
        if ( status < 0 ) {
            // When error happens, stop at the current multigrid level!
            if ( prtlvl > PRINT_MIN ) {
                printf("### HAZMATH WARNING: Forming aggregates on level-%lld failed!\n", (long long )lvl);
            }
            status = SUCCESS; break;
        }

        /*-- Form Prolongation --*/
        form_tentative_p(&vertices[lvl], &tentative_p[lvl], mgl[0].near_kernel_basis, lvl+1, num_aggs[lvl]);

        /* -- Form smoothed prolongation -- */
        smooth_aggregation_p(&mgl[lvl].A, &tentative_p[lvl], &mgl[lvl].P, param, lvl+1, &Neighbor[lvl]);

        // Check 2: Is coarse sparse too small?
        if ( mgl[lvl].P.col < MIN_CDOF ) break;

#if 0
        // Check 3: Does this coarsening step too aggressive?
        if ( mgl[lvl].P.row > mgl[lvl].P.col * MAX_CRATE ) {
            if ( prtlvl > PRINT_MIN ) {
                printf("### HAZMATH WARNING: Coarsening might be too aggressive!\n");
                printf("### HAZMATH: Fine level = %lld, coarse level = %lld. Discard!\n",
                       (long long )mgl[lvl].P.row, (long long )mgl[lvl].P.col);
            }
            break;
        }
#endif

#if 0
        // Check 4: Is this coarsening ratio too small?
        if ( (REAL)mgl[lvl].P.col > mgl[lvl].P.row * MIN_CRATE ) {
            if ( prtlvl > PRINT_MIN ) {
                printf("### WARNING: Coarsening rate is too small!\n");
                printf("### WARNING: Fine level = %lld, coarse level = %lld. Discard!\n",
                       (long long )mgl[lvl].P.row, (long long )mgl[lvl].P.col);
            }
            break;
        }
#endif

        /*-- Form restriction --*/
        dcsr_trans(&mgl[lvl].P, &mgl[lvl].R);

        /*-- Form coarse level stiffness matrix --*/
        dcsr_rap(&mgl[lvl].R, &mgl[lvl].A, &mgl[lvl].P,
                               &mgl[lvl+1].A);

        /*-- Form coarse level mass matrix --*/
        dcsr_rap(&mgl[lvl].R, &mgl[lvl].M, &mgl[lvl].P, &mgl[lvl+1].M);

        dcsr_free(&Neighbor[lvl]);
        ivec_free(&vertices[lvl]);

        ++lvl;

    }

/* #if DIAGONAL_PREF */
/*     dcsr_diagpref(&mgl[lvl].A); // reorder each row to make diagonal appear first */
/* #endif */

    // Setup coarse level systems for direct solvers
    switch (csolver) {

      //#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
            // Need to sort the matrix A for UMFPACK to work
            dCSRmat Ac_tran;
            dcsr_trans(&mgl[lvl].A, &Ac_tran);
            // It is equivalent to do transpose and then sort
            //     fasp_dcsr_trans(&mgl[lvl].A, &Ac_tran);
            //     fasp_dcsr_sort(&Ac_tran);
            dcsr_cp(&Ac_tran, &mgl[lvl].A);
            dcsr_free(&Ac_tran);
            mgl[lvl].Numeric = hazmath_factorize(&mgl[lvl].A, 0);
            break;
        }
	  //#endif
        default:
            // Do nothing!
            break;
    }

    // setup total level number and current level
    mgl[0].num_levels = max_levels = lvl+1;
    mgl[0].w          = dvec_create(m);

    for ( lvl = 1; lvl < max_levels; ++lvl) {
        INT mm = mgl[lvl].A.row;
        mgl[lvl].num_levels = max_levels;
        mgl[lvl].b          = dvec_create(mm);
        mgl[lvl].x          = dvec_create(mm);

        mgl[lvl].cycle_type     = cycle_type; // initialize cycle type!

        if ( cycle_type == NL_AMLI_CYCLE )
            mgl[lvl].w = dvec_create(3*mm);
        else
            mgl[lvl].w = dvec_create(2*mm);
    }

    if ( prtlvl > PRINT_NONE ) {
        get_time(&setup_end);
        print_amg_complexity(mgl,prtlvl);
        print_cputime("Smoothed aggregation setup", setup_end - setup_start);
    }

    for (i=0; i<max_levels; i++) {
      dcsr_free(&Neighbor[i]);
      ivec_free(&vertices[i]);
      dcsr_free(&tentative_p[i]);
    }

    free(Neighbor);
    free(vertices);
    free(num_aggs);
    free(tentative_p);

    return status;
}

/***********************************************************************************************/
/**
 * \fn static SHORT amg_setup_unsmoothP_unsmoothR_bsr (AMG_data_bsr *mgl,
 *                                                     AMG_param *param)
 *
 * \brief Set up phase of plain aggregation AMG, using unsmoothed P and unsmoothed A
 *        in BSR format
 *
 * \param mgl    Pointer to AMG data: AMG_data_bsr
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       SUCCESS if succeed, error otherwise
 *
 * \author Xiaozhe Hu
 * \date   03/16/2012
 *
 */
static SHORT amg_setup_unsmoothP_unsmoothR_bsr(AMG_data_bsr   *mgl,
                                               AMG_param      *param)
{
    const SHORT CondType = 1; // Condensation method used for AMG

    const SHORT prtlvl   = param->print_level;
    const SHORT csolver  = param->coarse_solver;
    const SHORT min_cdof = MAX(param->coarse_dof,50);
    const INT   m        = mgl[0].A.ROW;
    const INT   nb       = mgl[0].A.nb;

    SHORT     max_levels = param->max_levels;
    SHORT     i, lvl = 0, status = SUCCESS;
    REAL      setup_start, setup_end;

    //AMG_data *mgl_csr = amg_data_create(max_levels);

    dCSRmat   temp1, temp2;

#if DEBUG_MODE > 0
    printf("### DEBUG: [-Begin-] %s ...\n", __FUNCTION__);
    printf("### DEBUG: nr=%lld, nc=%lld, nnz=%lld\n",
           (long long )mgl[0].A.ROW, (long long )mgl[0].A.COL, (long long )mgl[0].A.NNZ);
#endif

    get_time(&setup_start);

    /*-----------------------*/
    /*--local working array--*/
    /*-----------------------*/
    // level info (fine: 0; coarse: 1)
    ivector *vertices = (ivector *)calloc(max_levels, sizeof(ivector));

    //each elvel stores the information of the number of aggregations
    INT *num_aggs = (INT *)calloc(max_levels, sizeof(INT));

    // each level stores the information of the strongly coupled neighborhoods
    dCSRmat *Neighbor = (dCSRmat *)calloc(max_levels, sizeof(dCSRmat));

    for ( i=0; i<max_levels; ++i ) num_aggs[i] = 0;

    /*------------------------------------------*/
    /*-- setup null spaces for whole Jacobian --*/
    /*------------------------------------------*/
    /*
     mgl[0].near_kernel_dim   = 1;
     mgl[0].near_kernel_basis = (REAL **)fasp_mem_calloc(mgl->near_kernel_dim, sizeof(REAL*));

     for ( i=0; i < mgl->near_kernel_dim; ++i ) mgl[0].near_kernel_basis[i] = NULL;
     */

    /*----------------------------*/
    /*--- checking aggregation ---*/
    /*----------------------------*/
    // Main AMG setup loop
    while ( (mgl[lvl].A.ROW > min_cdof) && (lvl < max_levels-1) ) {

        /*-- get the diagonal inverse --*/
        mgl[lvl].diaginv = dbsr_getdiaginv(&mgl[lvl].A);

        switch ( CondType ) {
            case 2:
                mgl[lvl].PP = condenseBSR(&mgl[lvl].A); break;
            default:
                mgl[lvl].PP = condenseBSRLinf(&mgl[lvl].A); break;
        }

        /*-- Aggregation --*/
        switch ( param->aggregation_type ) {

            case VMB: // VMB aggregation

                status = aggregation_vmb(&mgl[lvl].PP, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl], lvl);
                break;

            case HEC: // Heavy edge coarsening aggregation

                status = aggregation_hec(&mgl[lvl].PP, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl], lvl);
                break;

            case HEM: // Heavy edge matching
                status = aggregation_hem(&mgl[lvl].PP, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl], lvl);
                break;

            default: // wrong aggregation type
                status = ERROR_AMG_AGG_TYPE;
                check_error(status, __FUNCTION__);
                break;
        }

        if ( status < 0 ) {
            // When error happens, force solver to use the current multigrid levels!
            if ( prtlvl > PRINT_MIN ) {
                printf("### WARNING: Forming aggregates on level-%lld failed!\n", (long long )lvl);
            }
            status = SUCCESS; break;
        }

        /* -- Form Prolongation --*/
        if ( lvl == 0 && mgl[0].near_kernel_dim >0 ) {
            form_tentative_p_bsr(&vertices[lvl], &mgl[lvl].P, &mgl[0],
                                 num_aggs[lvl], mgl[0].near_kernel_dim,
                                 mgl[0].near_kernel_basis);
        }
        else {
            form_boolean_p_bsr(&vertices[lvl], &mgl[lvl].P, &mgl[0], num_aggs[lvl]);
        }

        /*-- Form resitriction --*/
        dbsr_trans(&mgl[lvl].P, &mgl[lvl].R);

        /*-- Form coarse level stiffness matrix --*/
        dbsr_rap(&mgl[lvl].R, &mgl[lvl].A, &mgl[lvl].P, &mgl[lvl+1].A);

        /* -- Form extra near kernal space if needed --*/
        if (mgl[lvl].A_nk != NULL){

            mgl[lvl+1].A_nk = (dCSRmat *)calloc(1, sizeof(dCSRmat));
            mgl[lvl+1].P_nk = (dCSRmat *)calloc(1, sizeof(dCSRmat));
            mgl[lvl+1].R_nk = (dCSRmat *)calloc(1, sizeof(dCSRmat));

            temp1 = dbsr_2_dcsr(&mgl[lvl].R);
            dcsr_mxm(&temp1, mgl[lvl].P_nk, mgl[lvl+1].P_nk);
            dcsr_trans(mgl[lvl+1].P_nk, mgl[lvl+1].R_nk);
            temp2 = dbsr_2_dcsr(&mgl[lvl+1].A);
            dcsr_rap(mgl[lvl+1].R_nk, &temp2, mgl[lvl+1].P_nk, mgl[lvl+1].A_nk);
            dcsr_free(&temp1);
            dcsr_free(&temp2);

        }

        dcsr_free(&Neighbor[lvl]);
        ivec_free(&vertices[lvl]);

        ++lvl;
    }

    // Setup coarse level systems for direct solvers (BSR version)
    switch (csolver) {

      //#if WITH_UMFPACK
        case SOLVER_UMFPACK: {
            // Need to sort the matrix A for UMFPACK to work
            mgl[lvl].Ac = dbsr_2_dcsr(&mgl[lvl].A);
	    dCSRmat Ac_tran=dcsr_create(mgl[lvl].Ac.col,	\
				      mgl[lvl].Ac.row,	\
				      mgl[lvl].Ac.nnz);
            dcsr_transz(&mgl[lvl].Ac, NULL, &Ac_tran);
            dcsr_cp(&Ac_tran, &mgl[lvl].Ac);
            dcsr_free(&Ac_tran);
            mgl[lvl].Numeric = hazmath_factorize(&mgl[lvl].Ac, 0);
            break;
        }
	  //#endif

        default:
            // Do nothing!
            break;
    }


    // setup total level number and current level
    mgl[0].num_levels = max_levels = lvl+1;
    mgl[0].w = dvec_create(3*m*nb);

    if (mgl[0].A_nk != NULL){

#if WITH_UMFPACK
        // Need to sort the matrix A_nk for UMFPACK
      dcsr_free(&temp1); // just in case:::
      temp1=dcsr_create(mgl[lvl].A_nk.col,	\
			mgl[lvl].A_nk.row,			\
			mgl[lvl].A_nk.nnz);
      dcsr_transz(mgl[0].A_nk, NULL, &temp1);
      dcsr_cp(&temp1, mgl[0].A_nk);
      dcsr_free(&temp1);
#endif

    }

    for ( lvl = 1; lvl < max_levels; lvl++ ) {
        const INT mm = mgl[lvl].A.ROW*nb;
        mgl[lvl].num_levels = max_levels;
        mgl[lvl].b          = dvec_create(mm);
        mgl[lvl].x          = dvec_create(mm);
        mgl[lvl].w          = dvec_create(3*mm);

        if (mgl[lvl].A_nk != NULL){

#if WITH_UMFPACK
            // Need to sort the matrix A_nk for UMFPACK
	  temp1=dcsr_create(mgl[lvl].A_nk.col,			\
			    mgl[lvl].A_nk.row,			\
			    mgl[lvl].A_nk.nnz);
	  dcsr_free(&temp1); // just in case:::
	  dcsr_transz(mgl[lvl].A_nk, NULL, &temp1);
	  dcsr_cp(&temp1, mgl[lvl].A_nk);
	  dcsr_free(&temp1);
#endif

        }

    }

    if ( prtlvl > PRINT_NONE ) {
        get_time(&setup_end);
        print_amgcomplexity_bsr(mgl,prtlvl);
        print_cputime("Unsmoothed aggregation (BSR) setup", setup_end - setup_start);
    }

    free(vertices); vertices = NULL;
    free(num_aggs); num_aggs = NULL;
    free(Neighbor); Neighbor = NULL;

#if DEBUG_MODE > 0
    printf("### DEBUG: [--End--] %s ...\n", __FUNCTION__);
#endif

    return status;
}

/***********************************************************************************************/
/**
 * \fn static SHORT amg_setup_general_bdcsr (AMG_data_bdcsr *mgl,
 *                                                     AMG_param *param)
 *
 * \brief Set up phase of plain aggregation AMG, using unsmoothed P and unsmoothed A
 *        in block_dCSR format
 *
 * \param mgl    Pointer to AMG data: AMG_data_bdcsr
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       SUCCESS if succeed, error otherwise
 *
 * \author Xiaozhe Hu
 * \date   04/16/2012
 *
 * \note    Assume block_row = block_col !! -- Xiaozhe Hu
 *
 */
static SHORT amg_setup_general_bdcsr(AMG_data_bdcsr *mgl,
                                     AMG_param      *param)
{
    const SHORT prtlvl   = param->print_level;
    const SHORT csolver  = param->coarse_solver;
    //const SHORT min_cdof = MAX(param->coarse_dof,50);
    const INT   brow     = mgl[0].A.brow;
    const INT   bcol     = mgl[0].A.bcol;

    SHORT max_levels = param->max_levels;
    SHORT i, j,lvl = 0, status = SUCCESS;
    REAL  setup_start, setup_end;
    INT   total_row, total_col, total_nnz;

    // get total size
    bdcsr_get_total_size(&(mgl[0].A), &total_row, &total_col, &total_nnz);

#if DEBUG_MODE > 0
    printf("### DEBUG: [-Begin-] %s ...\n", __FUNCTION__);
    printf("### DEBUG: nblock_row=%lld, nblock_col=%lld, total_row=%lld, total_col=%lld, total_nnz=%lld\n",
            (long long )brow, (long long )bcol, (long long )total_row, (long long )total_col, (long long )total_nnz);
#endif

    // diagonal matrices for coarsening
    dCSRmat *A_diag = mgl[0].A_diag;

    // AMG data for each diaongal block
    AMG_data **mgl_diag = (AMG_data **)calloc(brow, sizeof(AMG_data *));

    // local variable
    dCSRmat temp_mat;

    /*---------------------------*/
    /*--Main step for AMG setup--*/
    /*---------------------------*/
    get_time(&setup_start);

    // setup AMG for each diagonal block (given by A_diag)
    for (i=0; i<brow; i++){

        if ( prtlvl > PRINT_NONE ) printf("\n Diagonal block %lld ...\n", (long long )i);

        /* set AMG for diagonal blocks */
        mgl_diag[i] = amg_data_create(max_levels);
        dcsr_alloc(A_diag[i].row, A_diag[i].row, A_diag[i].nnz, &mgl_diag[i][0].A);
        dcsr_cp(&(A_diag[i]), &mgl_diag[i][0].A);
        mgl_diag[i][0].b=dvec_create(A_diag[i].row);
        mgl_diag[i][0].x=dvec_create(A_diag[i].row);

        switch (param->AMG_type) {

            case UA_AMG: // Unsmoothed Aggregation AMG
                if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                status = amg_setup_ua(mgl_diag[i], param);
                break;

            case SA_AMG: // Smoothed Aggregation AMG
                if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
                status = amg_setup_sa(mgl_diag[i], param);
                break;

            default: // UA AMG
                if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                status = amg_setup_ua(mgl_diag[i], param);
                break;

        }

    }

    // save mgl_diag data
    mgl[0].mgl_diag = mgl_diag;

    // set the total levels  (TODO: need better way to do this --Xiaozhe Hu)
    mgl[0].num_levels = mgl_diag[0]->num_levels;
    for (i=0; i<brow; i++){
        //mgl[0].num_levels = MAX(mgl[0].num_levels, mgl_diag[i]->num_levels);
        mgl[0].num_levels = MIN(mgl[0].num_levels, mgl_diag[i]->num_levels);
    }
    max_levels = mgl[0].num_levels;

    // construct P and R in the block_dCSR format
    for (lvl=0; lvl<max_levels-1; lvl++){

        //printf("level = %d\n", lvl);

        // allocate P
        bdcsr_alloc(brow, bcol, &(mgl[lvl].P));
        // allocate R
        bdcsr_alloc(brow, bcol, &(mgl[lvl].R));

        for (i=0; i<brow; i++){

            // copy P
            dcsr_alloc(mgl_diag[i][lvl].P.row, mgl_diag[i][lvl].P.col, mgl_diag[i][lvl].P.nnz, mgl[lvl].P.blocks[i*brow+i]);
            dcsr_cp(&mgl_diag[i][lvl].P, mgl[lvl].P.blocks[i*brow+i]);

            // copy R
            dcsr_alloc(mgl_diag[i][lvl].R.row, mgl_diag[i][lvl].R.col, mgl_diag[i][lvl].R.nnz, mgl[lvl].R.blocks[i*brow+i]);
            dcsr_cp(&mgl_diag[i][lvl].R, mgl[lvl].R.blocks[i*brow+i]);

            // set other blocks to be Null
            for (j=0; j<brow; j++){
                if (i != j){
                    mgl[lvl].P.blocks[i*brow+j] = NULL;
                    mgl[lvl].R.blocks[i*brow+j] = NULL;
                }
            }

        }

    }

    // form coarse level matrices
    for (lvl=0; lvl<max_levels-1; lvl++){

        // allocate
        bdcsr_alloc(brow, bcol, &(mgl[lvl+1].A));

        // form coarse level matrices
        for (i=0; i<brow; i++){
            for (j=0; j<brow; j++){

                if (i==j)  // diagonal block
                {
                    dcsr_rap(mgl[lvl].R.blocks[i*brow+i], mgl[lvl].A.blocks[i*brow+j], mgl[lvl].P.blocks[j*brow+j], mgl[lvl+1].A.blocks[i*brow+j]);
                }
                else // off-diagonal blocks
                {
                    // temp = R*A
                    dcsr_mxm(mgl[lvl].R.blocks[i*brow+i], mgl[lvl].A.blocks[i*brow+j], &temp_mat);
                    // Ac = temp*P
                    dcsr_mxm (&temp_mat, mgl[lvl].P.blocks[j*brow+j], mgl[lvl+1].A.blocks[i*brow+j]);
                    // cleam temp mat
                    dcsr_free(&temp_mat);
                }

            }
        }
    }

    // form coarse level A_diag
    for (lvl=1; lvl<max_levels; lvl++){
        // allocate
        mgl[lvl].A_diag = (dCSRmat *)calloc(brow, sizeof(dCSRmat));

        // form A_diag
        for (i=0; i<brow; i++){
            dcsr_alloc(mgl_diag[i][lvl].A.row, mgl_diag[i][lvl].A.col, mgl_diag[i][lvl].A.nnz, &mgl[lvl].A_diag[i]);
            dcsr_cp(&(mgl_diag[i][lvl].A), &mgl[lvl].A_diag[i]);
        }
    }

    // Setup coarse level systems for direct solvers (block_dCSRmat version)
    lvl = max_levels-1;
    switch (csolver) {

      //#if WITH_SUITESPARSE
    case SOLVER_UMFPACK: {
  // Need to sort the matrix A for UMFPACK to work
      mgl[lvl].Ac = bdcsr_2_dcsr(&mgl[lvl].A);
      dCSRmat Ac_tran=dcsr_create(mgl[lvl].Ac.col,mgl[lvl].Ac.row,mgl[lvl].Ac.nnz);
      dcsr_transz(&mgl[lvl].Ac, NULL, &Ac_tran);
      dcsr_cp(&Ac_tran, &mgl[lvl].Ac);
      dcsr_free(&Ac_tran);
      mgl[lvl].Numeric = hazmath_factorize(&mgl[lvl].Ac, 0);
      break;
    }
	  //#endif

        default:
            // Do nothing!
            break;
    }

    // allocate workspace on the fine level
    mgl[0].w = dvec_create(3*(mgl[0].b.row));

    // allocation on coarse levels
    for ( lvl = 1; lvl < max_levels; lvl++ ) {
        bdcsr_get_total_size(&(mgl[lvl].A), &total_row, &total_col, &total_nnz);
        mgl[lvl].num_levels = max_levels;
        mgl[lvl].b          = dvec_create(total_col);
        mgl[lvl].x          = dvec_create(total_col);
        mgl[lvl].w          = dvec_create(3*total_col);
    }

    if ( prtlvl > PRINT_NONE ) {
        get_time(&setup_end);
        print_cputime("Block dCSRmat AMG setup", setup_end - setup_start);
    }

#if DEBUG_MODE > 0
    printf("### DEBUG: [--End--] %s ...\n", __FUNCTION__);
#endif

    return status;
}


/***********************************************************************************************/
/**
 * \fn static SHORT amg_setup_metric_bdcsr (AMG_data_bdcsr *mgl,
 *                                          AMG_param *param)
 *
 * \brief Set up phase of metric AMG, using unsmoothed P and unsmoothed A
 *        in block_dCSR format
 *
 * \param mgl    Pointer to AMG data: AMG_data_bdcsr
 * \param param  Pointer to AMG parameters: AMG_param

 *
 * \return       SUCCESS if succeed, error otherwise
 *
 * \author Xiaozhe Hu
 * \date   05/15/2012
 *
 * \note    Assume block_row = block_col !! -- Xiaozhe Hu
 *
 */
static SHORT amg_setup_bdcsr_metric(AMG_data_bdcsr *mgl,
                                     AMG_param      *param)
{
    const SHORT prtlvl   = param->print_level;
    const SHORT csolver  = param->coarse_solver;
    //const SHORT min_cdof = MAX(param->coarse_dof,50);
    const INT   brow     = mgl[0].A.brow;
    const INT   bcol     = mgl[0].A.bcol;

    SHORT max_levels = param->max_levels;
    SHORT i, j,lvl = 0, status = SUCCESS;
    REAL  setup_start, setup_end;
    INT   total_row, total_col, total_nnz;

    // get total size
    bdcsr_get_total_size(&(mgl[0].A), &total_row, &total_col, &total_nnz);

#if DEBUG_MODE > 0
    printf("### DEBUG: [-Begin-] %s ...\n", __FUNCTION__);
    printf("### DEBUG: nblock_row=%lld, nblock_col=%lld, total_row=%lld, total_col=%lld, total_nnz=%lld\n",
            (long long )brow, (long long )bcol, (long long )total_row, (long long )total_col, (long long )total_nnz);
#endif

    // diagonal matrices for coarsening
    dCSRmat *A_diag = mgl[0].A_diag;

    // AMG data for each diaongal block
    AMG_data **mgl_diag = (AMG_data **)calloc(brow, sizeof(AMG_data *));

    // local variable
    dCSRmat temp_mat;

    /*---------------------------*/
    /*--Main step for AMG setup--*/
    /*---------------------------*/
    get_time(&setup_start);

    // setup AMG for each diagonal block (given by A_diag)
    for (i=0; i<brow; i++){

        if ( prtlvl > PRINT_NONE ) printf("\n Diagonal block %lld ...\n", (long long )i);

        /* set AMG for diagonal blocks */
        mgl_diag[i] = amg_data_create(max_levels);
        dcsr_alloc(A_diag[i].row, A_diag[i].row, A_diag[i].nnz, &mgl_diag[i][0].A);
        dcsr_cp(&(A_diag[i]), &mgl_diag[i][0].A);
        mgl_diag[i][0].b=dvec_create(A_diag[i].row);
        mgl_diag[i][0].x=dvec_create(A_diag[i].row);
        switch (param->AMG_type) {

            case MUA_AMG: // Unsmoothed Aggregation AMG
                if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                status = amg_setup_ua(mgl_diag[i], param);
                break;

            case MSA_AMG: // Smoothed Aggregation AMG
                if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
                status = amg_setup_sa(mgl_diag[i], param);
                break;

            default: // UA AMG
                if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
                status = amg_setup_ua(mgl_diag[i], param);
                break;

        }

    }

    // save mgl_diag data
    mgl[0].mgl_diag = mgl_diag;

    // set the total levels  (TODO: need better way to do this --Xiaozhe Hu)
    mgl[0].num_levels = mgl_diag[0]->num_levels;
    for (i=0; i<brow; i++){
        //mgl[0].num_levels = MAX(mgl[0].num_levels, mgl_diag[i]->num_levels);
        mgl[0].num_levels = MIN(mgl[0].num_levels, mgl_diag[i]->num_levels);
    }
    max_levels = mgl[0].num_levels;

    // construct P and R in the block_dCSR format
    for (lvl=0; lvl<max_levels-1; lvl++){

        //printf("level = %d\n", lvl);

        // allocate P
        bdcsr_alloc(brow, bcol, &(mgl[lvl].P));
        // allocate R
        bdcsr_alloc(brow, bcol, &(mgl[lvl].R));

        for (i=0; i<brow; i++){

            // copy P
            dcsr_alloc(mgl_diag[i][lvl].P.row, mgl_diag[i][lvl].P.col, mgl_diag[i][lvl].P.nnz, mgl[lvl].P.blocks[i*brow+i]);
            dcsr_cp(&mgl_diag[i][lvl].P, mgl[lvl].P.blocks[i*brow+i]);

            // copy R
            dcsr_alloc(mgl_diag[i][lvl].R.row, mgl_diag[i][lvl].R.col, mgl_diag[i][lvl].R.nnz, mgl[lvl].R.blocks[i*brow+i]);
            dcsr_cp(&mgl_diag[i][lvl].R, mgl[lvl].R.blocks[i*brow+i]);

        }

    }

    // form coarse level matrices
    for (lvl=0; lvl<max_levels-1; lvl++){
        //printf("level = %d\n", lvl);

        // allocate
        bdcsr_alloc(brow, bcol, &(mgl[lvl+1].A));

        // form coarse level matrices
        for (i=0; i<brow; i++){
            for (j=0; j<brow; j++){

                if (i==j)  // diagonal block
                {
                    dcsr_rap(mgl[lvl].R.blocks[i*brow+i], mgl[lvl].A.blocks[i*brow+j], mgl[lvl].P.blocks[j*brow+j], mgl[lvl+1].A.blocks[i*brow+j]);
                }
                else // off-diagonal blocks
                {
                    // temp = R*A
                    dcsr_mxm(mgl[lvl].R.blocks[i*brow+i], mgl[lvl].A.blocks[i*brow+j], &temp_mat);
                    // Ac = temp*P
                    dcsr_mxm (&temp_mat, mgl[lvl].P.blocks[j*brow+j], mgl[lvl+1].A.blocks[i*brow+j]);
                    // cleam temp mat
                    dcsr_free(&temp_mat);
                }
                //

                //printf("done!\n");
                //getchar();
            }
        }
    }

    // form coarse level A_diag
    for (lvl=1; lvl<max_levels; lvl++){
        // allocate
        mgl[lvl].A_diag = (dCSRmat *)calloc(brow, sizeof(dCSRmat));

        // form A_diag
        for (i=0; i<brow; i++){
            dcsr_alloc(mgl_diag[i][lvl].A.row, mgl_diag[i][lvl].A.col, mgl_diag[i][lvl].A.nnz, &mgl[lvl].A_diag[i]);
            dcsr_cp(&(mgl_diag[i][lvl].A), &mgl[lvl].A_diag[i]);
        }
    }

    // Setup coarse level systems for direct solvers (block_dCSRmat version)
    lvl = max_levels;
    switch (csolver) {

      //#if WITH_UMFPACK
        case SOLVER_UMFPACK: {
            // Need to sort the matrix A for UMFPACK to work
            mgl[lvl].Ac = bdcsr_2_dcsr(&mgl[lvl].A);
	    dCSRmat Ac_tran=dcsr_create(mgl[lvl].Ac.col,mgl[lvl].Ac.row,mgl[lvl].Ac.nnz);
            dcsr_transz(&mgl[lvl].Ac, NULL, &Ac_tran);
            dcsr_cp(&Ac_tran, &mgl[lvl].Ac);
            dcsr_free(&Ac_tran);
            mgl[lvl].Numeric = hazmath_factorize(&mgl[lvl].Ac, 0);
            break;
        }
	  //#endif

        default:
            // Do nothing!
            break;
    }

    // allocate workspace on the fine level
    mgl[0].w = dvec_create(6*(mgl[0].b.row));

    // get the interface submatrix
    //block_dCSRmat A_gamma;
    mgl[0].A_gamma = (block_dCSRmat *)calloc(1, sizeof(block_dCSRmat));
    bdcsr_alloc(brow, bcol, mgl[0].A_gamma);

    // assume 2 x 2 block structure
    INT size0   = mgl[0].interface_dof->row;
    INT size1   = mgl[0].interface_dof->row;
    //
    //Grab indices and values
    INT *gamma0 = calloc(mgl[0].interface_dof->nnz,sizeof(INT));
    INT *gamma1 = calloc(mgl[0].interface_dof->nnz,sizeof(INT));
    INT ij,nnz_g=0;
    for (i=0;i<mgl[0].interface_dof->row;++i){
      for(ij=mgl[0].interface_dof->IA[i];		\
	  ij<mgl[0].interface_dof->IA[i+1];++ij){
	// let us drop small entries, so we have 1-1 when r=0:
	if(fabs(mgl[0].interface_dof->val[ij])<1e-10) continue;
	gamma1[nnz_g]=i;
	gamma0[nnz_g]=mgl[0].interface_dof->JA[ij];
	nnz_g++;
      }
    }
    //dcsr_write_dcoo("A00.dat", mgl[0].A.blocks[0]);
    //for (i=0; i<mgl[0].interface_dof->nnz; i++) printf("idx[%d]=%d\n", i, mgl[0].interface_dof->JA[i]);
    //getchar();

    dcsr_getblk(mgl[0].A.blocks[0], gamma0, gamma0, size0, size0, mgl[0].A_gamma->blocks[0]);
    dcsr_getblk(mgl[0].A.blocks[1], gamma0, gamma1, size0, size1, mgl[0].A_gamma->blocks[1]);
    dcsr_getblk(mgl[0].A.blocks[2], gamma1, gamma0, size1, size0, mgl[0].A_gamma->blocks[2]);
    dcsr_getblk(mgl[0].A.blocks[3], gamma1, gamma1, size1, size1, mgl[0].A_gamma->blocks[3]);
    // free something that should not have been allocated
    free(gamma0);
    free(gamma1);
    /*
    dcsr_write_dcoo("AT00.dat", mgl[0].A_gamma->blocks[0]);
    dcsr_write_dcoo("AT01.dat", mgl[0].A_gamma->blocks[1]);
    dcsr_write_dcoo("AT10.dat", mgl[0].A_gamma->blocks[2]);
    dcsr_write_dcoo("AT11.dat", mgl[0].A_gamma->blocks[3]);
    getchar();
    */

    // convert to dBSRmat format
    mgl[0].A_gamma_bsr = bdcsr_2_dbsr(mgl[0].A_gamma);
    mgl[0].A_gamma_diaginv = dbsr_getdiaginv(&mgl[0].A_gamma_bsr);

    // allocation on coarse level
    for ( lvl = 1; lvl < max_levels; lvl++ ) {
        bdcsr_get_total_size(&(mgl[lvl].A), &total_row, &total_col, &total_nnz);
        mgl[lvl].num_levels = max_levels;
        mgl[lvl].b          = dvec_create(total_col);
        mgl[lvl].x          = dvec_create(total_col);
        mgl[lvl].w          = dvec_create(6*total_col);
    }

    if ( prtlvl > PRINT_NONE ) {
        get_time(&setup_end);
        print_cputime("Block dCSRmat AMG setup", setup_end - setup_start);
    }

#if DEBUG_MODE > 0
    printf("### DEBUG: [--End--] %s ...\n", __FUNCTION__);
#endif

    return status;
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
