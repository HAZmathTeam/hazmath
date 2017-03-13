/*! \file src/solver/amg_setup_ua.c
 *
 *  Unsmoothed Aggregation AMG: SETUP phase
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 12/24/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note   Done cleanup for releasing -- Xiaozhe Hu 03/11/2017
 *
 *  \todo   Add safe guard for the overall computatinal complexity -- Xiaozhe Hu
 *
 */

#include "hazmath.h"
#include "aggregation.inl"

static SHORT amg_setup_unsmoothP_unsmoothR(AMG_data *, AMG_param *);

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

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

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
    const SHORT min_cdof   = MAX(param->coarse_dof,50);
    const INT   m          = mgl[0].A.row;
    
    // local variables
    SHORT         max_levels = param->max_levels, lvl = 0, status = SUCCESS;
    INT           i;
    REAL          setup_start, setup_end;
    
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
    
    // Initialize AMLI coefficients
    if ( cycle_type == AMLI_CYCLE ) {
        const INT amlideg = param->amli_degree;
        param->amli_coef = (REAL *)calloc(amlideg+1,sizeof(REAL));
        REAL lambda_max = 2.0, lambda_min = lambda_max/4;
        amg_amli_coef(lambda_max, lambda_min, amlideg, param->amli_coef);
    }
    
#if DIAGONAL_PREF
    dcsr_diagpref(&mgl[0].A); // reorder each row to make diagonal appear first
#endif
    
    /*----------------------------*/
    /*--- checking aggregation ---*/
    /*----------------------------*/
    // Main AMG setup loop
    while ( (mgl[lvl].A.row > min_cdof) && (lvl < max_levels-1) ) {
        
        //printf("level = %d\n", lvl);
        //dcsr_write_dcoo("A.dat", &mgl[lvl].A);
        //getchar();

        /*-- Aggregation --*/
        switch ( param->aggregation_type ) {
                
            case VMB: // VMB aggregation     
                status = aggregation_vmb(&mgl[lvl].A, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl]);
                break;

            case HEC: // Heavy edge coarsening aggregation
                status = aggregation_hec(&mgl[lvl].A, &vertices[lvl], param,
                                         &Neighbor[lvl], &num_aggs[lvl]);
                //printf("level = %d\n", lvl);
                //ivector_write("vt.dat", &vertices[lvl]);
                //dcsr_write_dcoo("N.dat", &Neighbor[lvl]);
                //getchar();
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
                printf("### HAZMATH WARNING: Forming aggregates on level-%d failed!\n", lvl);
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
                printf("### HAZMATH: Fine level = %d, coarse level = %d. Discard!\n",
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
                printf("### WARNING: Fine level = %d, coarse level = %d. Discard!\n",
                       mgl[lvl].P.row, mgl[lvl].P.col);
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
    
#if DIAGONAL_PREF
    dcsr_diagpref(&mgl[lvl].A); // reorder each row to make diagonal appear first
#endif
    
    // Setup coarse level systems for direct solvers
    switch (csolver) {
 
#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
            // Need to sort the matrix A for UMFPACK to work
            dCSRmat Ac_tran;
            Ac_tran = dcsr_create(mgl[lvl].A.row, mgl[lvl].A.col, mgl[lvl].A.nnz);
            dcsr_trans(&mgl[lvl].A, &Ac_tran);
            // It is equivalent to do transpose and then sort
            //     fasp_dcsr_trans(&mgl[lvl].A, &Ac_tran);
            //     fasp_dcsr_sort(&Ac_tran);
            dcsr_cp(&Ac_tran, &mgl[lvl].A);
            dcsr_free(&Ac_tran);
            mgl[lvl].Numeric = umfpack_factorize(&mgl[lvl].A, 0);
            break;
        }
#endif
            
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

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
