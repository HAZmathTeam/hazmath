/*
 *  mgcycle.c
 *
 *  Abstract multigrid cycle -- non-recursive version
 *
 *  Created by James Adler and Xiaozhe Hu on 12/25/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

//#include <math.h>
//#include <time.h>

#include "hazmat.h"
#include "mg_util.inl"


/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

void mgcycle (AMG_data *mgl,
              AMG_param *param)
{
    /**
     * \fn void mgcycle (AMG_data *mgl, AMG_param *param)
     *
     * \brief Solve Ax=b with non-recursive multigrid cycle
     *
     * \param mgl    Pointer to AMG data: AMG_data
     * \param param  Pointer to AMG parameters: AMG_param
     *
     * \author Xiaozhe Hu
     * \date   12/25/2015
     *
     */
    
    const SHORT  prtlvl = param->print_level;
    const SHORT  amg_type = param->AMG_type;
    const SHORT  smoother = param->smoother;
    const SHORT  smooth_order = param->smooth_order;
    const SHORT  cycle_type = param->cycle_type;
    const SHORT  coarse_solver = param->coarse_solver;
    const SHORT  nl = mgl[0].num_levels;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol * 1e-4;
    const SHORT  ndeg = param->polynomial_degree;
    
    /*
    // Schwarz parameters
    Schwarz_param swzparam;
    if ( param->Schwarz_levels > 0 ) {
        swzparam.Schwarz_blksolver = param->Schwarz_blksolver;
    }
     */
    
    // local variables
    REAL alpha = 1.0;
    INT  num_lvl[MAX_AMG_LVL] = {0}, l = 0;
    
    
ForwardSweep:
    while ( l < nl-1 ) {
        
        num_lvl[l]++;
        
        // pre-smoothing with ILU method
        if ( l < mgl->ILU_levels ) {
            //fasp_smoother_dcsr_ilu(&mgl[l].A, &mgl[l].b, &mgl[l].x, &mgl[l].LU);
        }
        
        // or pre-smoothing with Schwarz method
        else if ( l < mgl->Schwarz_levels ) {
            /*
            switch (mgl[l].Schwarz.Schwarz_type) {
                case SCHWARZ_SYMMETRIC:
                    fasp_dcsr_Schwarz_forward_smoother(&mgl[l].Schwarz, &swzparam, 
                                                       &mgl[l].x, &mgl[l].b);
                    fasp_dcsr_Schwarz_backward_smoother(&mgl[l].Schwarz, &swzparam, 
                                                        &mgl[l].x, &mgl[l].b);
                    break;
                default:
                    fasp_dcsr_Schwarz_forward_smoother(&mgl[l].Schwarz, &swzparam, 
                                                       &mgl[l].x, &mgl[l].b);
                    break;
            }
             */
        }
        
        // or pre-smoothing with standard smoothers
        else {
            dcsr_presmoothing(smoother, &mgl[l].A, &mgl[l].b, &mgl[l].x,
                                   param->presmooth_iter, 0, mgl[l].A.row-1, 1,
                                   relax, ndeg, smooth_order, mgl[l].cfmark.val);
        }
        
        // form residual r = b - A x
        array_cp(mgl[l].A.row, mgl[l].b.val, mgl[l].w.val);
        dcsr_aAxpy(-1.0,&mgl[l].A, mgl[l].x.val, mgl[l].w.val);
        
        // restriction r1 = R*r0
        switch ( amg_type ) {
            case UA_AMG:
                dcsr_mxv_agg(&mgl[l].R, mgl[l].w.val, mgl[l+1].b.val);
                break;
            default:
                dcsr_mxv(&mgl[l].R, mgl[l].w.val, mgl[l+1].b.val);
                break;
        }
        
        // prepare for the next level
        ++l; dvec_set(mgl[l].A.row, &mgl[l].x, 0.0);
        
    }
    
    // If AMG only has one level or we have arrived at the coarsest level,
    // call the coarse space solver:
    switch ( coarse_solver ) {

/*
#if WITH_MUMPS
        case SOLVER_MUMPS: {
            // use MUMPS direct solver on the coarsest level
            mgl[nl-1].mumps.job = 2;
            fasp_solver_mumps_steps(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, &mgl[nl-1].mumps);
            break;
        }
#endif
            
#if WITH_SuperLU
        case SOLVER_SUPERLU: {
            // use SuperLU direct solver on the coarsest level
            fasp_solver_superlu(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, 0);
            break;
        }
#endif
            
#if WITH_UMFPACK
        case SOLVER_UMFPACK: {
            // use UMFPACK direct solver on the coarsest level
            fasp_umfpack_solve(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, mgl[nl-1].Numeric, 0);
            break;
        }
#endif
 */
            
        default:
            // use iterative solver on the coarsest level
            coarse_itsolver(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, tol, prtlvl);
            
    }
    
    // BackwardSweep:
    while ( l > 0 ) {
        
        --l;
        
        // find the optimal scaling factor alpha
        if ( param->coarse_scaling == ON ) {
            alpha = array_dotprod(mgl[l+1].A.row, mgl[l+1].x.val, mgl[l+1].b.val)
                  / dcsr_vmv(&mgl[l+1].A, mgl[l+1].x.val, mgl[l+1].x.val);
            alpha = MIN(alpha, 1.0); // Add this for safety! --Chensong on 10/04/2014
        }
        
        // prolongation u = u + alpha*P*e1
        switch ( amg_type ) {
            case UA_AMG:
                dcsr_aAxpy_agg(alpha, &mgl[l].P, mgl[l+1].x.val, mgl[l].x.val);
                break;
            default:
                dcsr_aAxpy(alpha, &mgl[l].P, mgl[l+1].x.val, mgl[l].x.val);
                break;
        }
        
        // post-smoothing with ILU method
        if ( l < mgl->ILU_levels ) {
            //fasp_smoother_dcsr_ilu(&mgl[l].A, &mgl[l].b, &mgl[l].x, &mgl[l].LU);
        }
        
        // post-smoothing with Schwarz method
        else if ( l < mgl->Schwarz_levels ) {
            /*
            switch (mgl[l].Schwarz.Schwarz_type) {
                case SCHWARZ_SYMMETRIC:
                    fasp_dcsr_Schwarz_backward_smoother(&mgl[l].Schwarz, &swzparam, 
                                                        &mgl[l].x, &mgl[l].b);
                    fasp_dcsr_Schwarz_forward_smoother(&mgl[l].Schwarz, &swzparam, 
                                                       &mgl[l].x, &mgl[l].b);
                    break;
                default:
                    fasp_dcsr_Schwarz_backward_smoother(&mgl[l].Schwarz, &swzparam, 
                                                        &mgl[l].x, &mgl[l].b);
                    break;
            }
             */
        }
        
        // post-smoothing with standard methods
        else {
            dcsr_postsmoothing(smoother, &mgl[l].A, &mgl[l].b, &mgl[l].x,
                                    param->postsmooth_iter, 0, mgl[l].A.row-1, -1,
                                    relax, ndeg, smooth_order, mgl[l].cfmark.val);
        }
        
        if ( num_lvl[l] < cycle_type ) break;
        else num_lvl[l] = 0;
    }
    
    if ( l > 0 ) goto ForwardSweep;
    
    
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
