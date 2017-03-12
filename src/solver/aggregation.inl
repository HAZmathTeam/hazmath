/*! \file amg_setup_ua.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 12/24/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *  ATTENTION: Do NOT use auto-indentation in this file!!!
 *
 * \todo    Add maximal weighted matching coarsning -- Xiaozhe Hu
 * \todo    Add maximal independent set aggregation -- Xiaozhe Hu
 * \todo    Add heavy edge coarsening               -- Xiaozhe Hu
 *
 * \note   Done cleanup for releasing -- Xiaozhe Hu 03/11/2017
 *
 */

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/
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

/**
 * \fn static void form_boolean_p (ivector *vertices, dCSRmat *tentp, INT levelNum,
 *                                 INT num_aggregations)
 *
 * \brief Form aggregation based on strong coupled neighbors
 *
 * \param vertices           Pointer to the aggregation of vertices
 * \param tentp              Pointer to the prolongation operators
 * \param levelNum           Level number
 * \param num_aggregations   Number of aggregations
 *
 * \author Xiaozhe Hu
 * \date   09/29/2009
 *
 * Modified by Xiaozhe Hu on 05/25/2014
 */
static void form_boolean_p(ivector *vertices,
                           dCSRmat *tentp,
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
            val[j] = 1.0;
            j ++;
        }
    }
}

/**
 * \fn static SHORT aggregation_vmb (dCSRmat *A, ivector *vertices, AMG_param *param,
 *                                   INT levelNum, dCSRmat *Neigh, INT *num_aggregations)
 *
 * \brief Form aggregation based on strong coupled neighbors
 *
 * \param A                 Pointer to the coefficient matrices
 * \param vertices          Pointer to the aggregation of vertices
 * \param param             Pointer to AMG parameters
 * \param levelNum          Level number
 * \param Neigh             Pointer to strongly coupled neighbors
 * \param num_aggregations  Pointer to number of aggregations
 *
 * \author Xiaozhe Hu
 * \date   09/29/2009
 *
 * \note Setup A, P, PT and levels using the unsmoothed aggregation algorithm;
 *       Refer to P. Vanek, J. Madel and M. Brezina
 *       "Algebraic Multigrid on Unstructured Meshes", 1994
 *
 */
static SHORT aggregation_vmb (dCSRmat *A,
                              ivector *vertices,
                              AMG_param *param,
                              INT levelNum,
                              dCSRmat *Neigh,
                              INT *num_aggregations)
{   
    const INT    row = A->row, col = A->col, nnz = A->IA[row]-A->IA[0];
    const INT  * AIA = A->IA, * AJA = A->JA;
    const REAL * Aval = A->val;
    const INT    max_aggregation = param->max_aggregation;
    
    // return status
    SHORT  status = SUCCESS;
    
    // local variables
    INT    num_left = row;
    INT    subset, count;
    INT  * num_each_agg;
    
    REAL   strongly_coupled, strongly_coupled2;
    INT    i, j, index, row_start, row_end;
    INT  * NIA = NULL, * NJA = NULL;
    REAL * Nval = NULL;
    
    dvector diag;
    dcsr_getdiag(0, A, &diag);  // get the diagonal entries
    
    //if ( GE(param->tentative_smooth, SMALLREAL) ) {
    //    strongly_coupled = param->strong_coupled * pow(0.5, levelNum-1);
    //}
    //else {
    strongly_coupled = param->strong_coupled;
    //}
    strongly_coupled2 = pow(strongly_coupled,2);
    
    /*------------------------------------------*/
    /*    Form strongly coupled neighborhood    */
    /*------------------------------------------*/
    dcsr_alloc(row, col, nnz, Neigh);
    
    NIA  = Neigh->IA;
    NJA  = Neigh->JA;
    Nval = Neigh->val;
    
    for ( i = row; i >= 0; i-- ) NIA[i] = AIA[i];
    
    for ( index = i = 0; i < row; ++i ) {
        NIA[i] = index;
        row_start = AIA[i]; row_end = AIA[i+1];
        for ( j = row_start; j < row_end; ++j ) {

            if ( (AJA[j] == i) 
              || ( (pow(Aval[j],2) >= strongly_coupled2*ABS(diag.val[i]*diag.val[AJA[j]])) && (Aval[j] < 0) )
               )
            {
                NJA[index] = AJA[j];
                Nval[index] = Aval[j];
                index++;
            }

        } // end for ( j = row_start; j < row_end; ++j )
    } // end for ( index = i = 0; i < row; ++i )
    NIA[row] = index;
    
    Neigh->nnz = index;
    Neigh->JA  = (INT*) realloc(Neigh->JA,  (Neigh->IA[row])*sizeof(INT));
    Neigh->val = (REAL*)realloc(Neigh->val, (Neigh->IA[row])*sizeof(REAL));
    
    NIA  = Neigh->IA;
    NJA  = Neigh->JA;
    Nval = Neigh->val;
    
    dvec_free(&diag);
    
    /*------------------------------------------*/
    /*             Initialization               */
    /*------------------------------------------*/
    ivec_alloc(row, vertices);
    iarray_set(row, vertices->val, -2);
    *num_aggregations = 0;
    
    /*-------------*/
    /*   Step 1.   */
    /*-------------*/
    for ( i = 0; i < row; ++i ) {
        if ( (AIA[i+1] - AIA[i]) == 1 ) {
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
    
    //for ( i = 0; i < *num_aggregations; i++ ) num_each_agg[i] = 0; // initialize
    
    for ( i = row; i--; ) {
        temp_C[i] = vertices->val[i];
        if ( vertices->val[i] >= 0 ) num_each_agg[vertices->val[i]] ++;
    }
    
    for ( i = 0; i < row; ++i ) {
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

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
