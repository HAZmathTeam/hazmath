/*! \file src/solver/amg_setup_c.c
 *
 *  Classical AMG: SETUP phase
 *
 *  Created by James Adler and Xiaozhe Hu on 12/26/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  Ref Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller
 *        Appendix P475 A.7 (by A. Brandt, P. Oswald and K. Stuben)
 *        Academic Press Inc., San Diego, CA, 2001.
 */

#include "hazmath.h"
#include "coarsening.inl"
#include "interpolation.inl"


/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/***********************************************************************************************/
SHORT amg_setup_c (AMG_data *mgl,
                   AMG_param *param)
{
    /**
     * \fn SHORT amg_setup_c(AMG_data *mgl, AMG_param *param)
     *
     * \brief Setup phase of classic AMG
     *
     * \param mgl    Pointer to AMG data: AMG_data
     * \param param  Pointer to AMG parameters: AMG_param
     *
     * \return       SUCCESS if successed; otherwise, error information.
     *
     */
    
    const SHORT prtlvl     = param->print_level;
    const SHORT cycle_type = param->cycle_type;
    const SHORT csolver    = param->coarse_solver;
    const SHORT min_cdof   = MAX(param->coarse_dof,MIN_CDOF);
    const INT   m          = mgl[0].A.row;

    // local variables
    SHORT         status = SUCCESS;
    INT           lvl = 0, max_lvls = param->max_levels;
    REAL          setup_start, setup_end;
    iCSRmat       Scouple; // strong n-couplings
    
    // level info (fine: 0; coarse: 1)
    ivector       vertices = ivec_create(m);
    
    get_time(&setup_start);
    
    // Make sure classical AMG will not call dcsr_mxv_agg!
    //param->tentative_smooth = 1.0;
    
    // If user want to use aggressive coarsening but did not specify number of
    // levels use aggressive coarsening, make sure apply aggressive coarsening
    // on the finest level only !!!
    if ( param->coarsening_type == COARSE_AC ) {
        param->aggressive_level = MAX(param->aggressive_level, 1);
    }
    
    // Initialize AMLI coefficients
    if ( cycle_type == AMLI_CYCLE ) {
        const INT amlideg = param->amli_degree;
        param->amli_coef = (REAL *)calloc(amlideg+1,sizeof(REAL));
        amg_amli_coef(2.0, 0.5, amlideg, param->amli_coef);
    }

#if DIAGONAL_PREF
    // Reorder each row to keep the diagonal entries appear first !!!
    dcsr_diagpref(&mgl[0].A);
#endif
    
    // Main AMG setup loop
    while ( (mgl[lvl].A.row > min_cdof) && (lvl < max_lvls-1) ) {
        
        
        /*-- Coarsening and form the structure of interpolation --*/
        status = amg_coarsening_c(&mgl[lvl].A, &vertices, &mgl[lvl].P,
		                                &Scouple, param);

           // Check 1: Did coarsening step succeeded?
        if ( status < 0 ) {
            // When error happens, stop at the current multigrid level!
            if ( prtlvl > 0 ) {
                printf("### WARNING: Could not find any C-variables!\n");
                printf("### WARNING: RS coarsening on level-%d failed!\n", lvl);
            }
	    free(Scouple.IA);
	    free(Scouple.JA);
            status = SUCCESS; break;
        }

        // Check 2: Is coarse sparse too small?
        if ( mgl[lvl].P.col < MIN_CDOF ) {
	  free(Scouple.IA);
	  free(Scouple.JA);
	  break;
        }
	
        // Check 3: Does this coarsening step too aggressive?
        if ( mgl[lvl].P.row > mgl[lvl].P.col * 10.0 ) {
            if ( prtlvl > 0 ) {
                printf("### WARNING: Coarsening might be too aggressive!\n");
                printf("### WARNING: Fine level = %d, coarse level = %d. Discard!\n",
                       mgl[lvl].P.row, mgl[lvl].P.col);
            }
	    free(Scouple.IA);
	    free(Scouple.JA);
            break;
        }
        
        /*-- Perform aggressive coarsening only up to the specified level --*/
        if ( mgl[lvl].P.col*1.5 > mgl[lvl].A.row ) param->coarsening_type = COARSE_C;
        if ( lvl == param->aggressive_level ) param->coarsening_type = COARSE_C;
        
        /*-- Store the C/F marker --*/
        {
            INT size = mgl[lvl].A.row;
            mgl[lvl].cfmark = ivec_create(size);
            memcpy(mgl[lvl].cfmark.val, vertices.val, size*sizeof(INT));
        }
        
        /*-- Form interpolation --*/
        amg_interp(&mgl[lvl].A, &vertices, &mgl[lvl].P, &Scouple, param);

        /*-- Form coarse level matrix: two RAP routines available! --*/
        dcsr_trans(&mgl[lvl].P, &mgl[lvl].R);

        dcsr_rap(&mgl[lvl].R, &mgl[lvl].A, &mgl[lvl].P, &mgl[lvl+1].A);
        
        /*-- Clean up Scouple generated in coarsening --*/
        free(Scouple.IA);
        free(Scouple.JA);
        
        // Check 4: Is the coarse matrix too dense?
        if ( mgl[lvl].A.nnz / mgl[lvl].A.row > mgl[lvl].A.col * 0.2 ) {
            if ( prtlvl > 0 ) {
                printf("### WARNING: Coarse matrix is too dense!\n");
                printf("### WARNING: m = n = %d, nnz = %d!\n",
                       mgl[lvl].A.col, mgl[lvl].A.nnz);
            }
            dcsr_free(&mgl[lvl+1].A);
            break;
        }
        
        ++lvl;
        
#if DIAGONAL_PREF
        // reorder each row to make diagonal appear first
        dcsr_diagpref(&mgl[lvl].A);
#endif
        
    }
    
    // Setup coarse level systems for direct solvers
    switch (csolver) {

            
#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
            // Need to sort the matrix A for UMFPACK to work
            dCSRmat Ac_tran;
            //Ac_tran = dcsr_create(mgl[lvl].A.row, mgl[lvl].A.col, mgl[lvl].A.nnz);
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
    mgl[0].num_levels = max_lvls = lvl+1;
    mgl[0].w          = dvec_create(m);
        
    for ( lvl = 1; lvl < max_lvls; ++lvl ) {
        INT mm = mgl[lvl].A.row;
        mgl[lvl].num_levels = max_lvls;
        mgl[lvl].b          = dvec_create(mm);
        mgl[lvl].x          = dvec_create(mm);
        
        mgl[lvl].cycle_type     = cycle_type; // initialize cycle type!
        
        // allocate work arrays for the solve phase
        if ( cycle_type == NL_AMLI_CYCLE )
            mgl[lvl].w = dvec_create(3*mm);
        else
            mgl[lvl].w = dvec_create(2*mm);
    }
    
    ivec_free(&vertices);

    if ( prtlvl > PRINT_NONE ) {
        get_time(&setup_end);
        print_amg_complexity(mgl, prtlvl);
        print_cputime("Classical AMG setup", setup_end - setup_start);
    }
    
    return status;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
