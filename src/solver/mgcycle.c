/*! \file src/solver/mgcycle.c
 *
 *  Abstract multigrid cycle
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 12/25/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 * \note   Done cleanup for releasing -- Xiaozhe Hu 03/12/2017
 *
 */

#include "hazmath.h"
#include "mg_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/**
 * \fn void mgcycle (AMG_data *mgl, AMG_param *param)
 *
 * \brief Solve Ax=b with non-recursive multigrid cycle (V- and W-cycle)
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \author Xiaozhe Hu
 * \date   12/25/2015
 *
 */
void mgcycle(AMG_data *mgl,
             AMG_param *param)
{    
    const SHORT  prtlvl = param->print_level;
    const SHORT  amg_type = param->AMG_type;
    const SHORT  smoother = param->smoother;
    const SHORT  cycle_type = param->cycle_type;
    const SHORT  coarse_solver = param->coarse_solver;
    const SHORT  nl = mgl[0].num_levels;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol * 1e-4;
    
    // local variables
    REAL alpha = 1.0;
    INT  num_lvl[MAX_AMG_LVL] = {0}, l = 0;
    
ForwardSweep:
    while ( l < nl-1 ) {
        
        num_lvl[l]++;
        
        // pre-smoothing with standard smoothers
        dcsr_presmoothing(smoother, &mgl[l].A, &mgl[l].b, &mgl[l].x,
                          param->presmooth_iter, 0, mgl[l].A.row-1, 1,
                          relax);
        
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
            
#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
            // use UMFPACK direct solver on the coarsest level
            umfpack_solve(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, mgl[nl-1].Numeric, 0);
            break;
        }
#endif
            
        default:
            // use iterative solver on the coarsest level
            coarse_itsolver(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, tol, prtlvl);
            break;
            
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

        // post-smoothing with standard methods
        dcsr_postsmoothing(smoother, &mgl[l].A, &mgl[l].b, &mgl[l].x,
                           param->postsmooth_iter, 0, mgl[l].A.row-1, -1,
                           relax);
        
        if ( num_lvl[l] < cycle_type ) break;
        else num_lvl[l] = 0;
    }
    
    if ( l > 0 ) goto ForwardSweep;
    
    
}


/**
 * \fn void amli (AMG_data *mgl, AMG_param *param, INT level)
 *
 * \brief Solve Ax=b with recursive AMLI-cycle
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 * \param level  Current level
 *
 * \author Xiaozhe Hu
 * \date   01/23/2011
 *
 * \note AMLI polynomial computed by the best approximation of 1/x.
 *       Refer to Johannes K. Kraus, Panayot S. Vassilevski, Ludmil T. Zikatanov,
 *       "Polynomial of best uniform approximation to $x^{-1}$ and smoothing in
 *        two-level methods", 2013.
 *
 */
void amli(AMG_data *mgl,
          AMG_param *param,
          INT level)
{    
    const SHORT  amg_type=param->AMG_type;
    const SHORT  prtlvl = param->print_level;
    const SHORT  smoother = param->smoother;
    const SHORT  coarse_solver = param->coarse_solver;
    const SHORT  degree= param->amli_degree;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol*1e-4;
    
    // local variables
    REAL   alpha  = 1.0;
    REAL * coef   = param->amli_coef;
    
    dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x;   // fine level b and x
    dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
    
    dCSRmat *A0 = &mgl[level].A;   // fine level matrix
    dCSRmat *A1 = &mgl[level+1].A; // coarse level matrix
    
    const INT m0 = A0->row, m1 = A1->row;
    
    REAL     *r        = mgl[level].w.val;      // work array for residual
    REAL     *r1       = mgl[level+1].w.val+m1; // work array for residual  
    
    if ( prtlvl >= PRINT_MOST )
        printf("AMLI level %d, smoother %d.\n", level, smoother);
    
    if ( level < mgl[level].num_levels-1 ) {
                
        // presmoothing
        dcsr_presmoothing(smoother,A0,b0,e0,param->presmooth_iter,
                          0,m0-1,1,relax);
        
        // form residual r = b - A x
        array_cp(m0,b0->val,r);
        dcsr_aAxpy(-1.0,A0,e0->val,r);
        
        // restriction r1 = R*r0
        switch (amg_type) {
            case UA_AMG:
                dcsr_mxv_agg(&mgl[level].R, r, b1->val); break;
            default:
                dcsr_mxv(&mgl[level].R, r, b1->val); break;
        }
        
        // coarse grid correction
        {
            INT i;
            
            array_cp(m1,b1->val,r1);
            
            for ( i=1; i<=degree; i++ ) {
                dvec_set(m1,e1,0.0);
                amli(mgl, param, level+1);
                
                // b1 = (coef[degree-i]/coef[degree])*r1 + A1*e1;
                // First, compute b1 = A1*e1
                dcsr_mxv(A1, e1->val, b1->val);
                // Then, compute b1 = b1 + (coef[degree-i]/coef[degree])*r1
                array_axpy(m1, coef[degree-i]/coef[degree], r1, b1->val);
            }
            
            dvec_set(m1,e1,0.0);
            amli(mgl, param, level+1);
        }
        
        // find the optimal scaling factor alpha
        array_ax(m1, coef[degree], e1->val);
        if ( param->coarse_scaling == ON ) {
            alpha = array_dotprod(m1, e1->val, r1)
            / dcsr_vmv(A1, e1->val, e1->val);
            alpha = MIN(alpha, 1.0);
        }
        
        // prolongation e0 = e0 + alpha * P * e1
        switch (amg_type) {
            case UA_AMG:
                dcsr_aAxpy_agg(alpha, &mgl[level].P, e1->val, e0->val);
                break;
            default:
                dcsr_aAxpy(alpha, &mgl[level].P, e1->val, e0->val);
                break;
        }
        
        // postsmoothing
        dcsr_postsmoothing(smoother,A0,b0,e0,param->postsmooth_iter,
                           0,m0-1,-1,relax);

    }
    
    else { // coarsest level solver
        
        switch (coarse_solver) {
                
#if WITH_SUITESPARSE
            case SOLVER_UMFPACK:
                // use UMFPACK direct solver on the coarsest level //
                umfpack_solve(A0, b0, e0, mgl[level].Numeric, 0);
                break;
#endif
                
            default:
                /* use iterative solver on the coarsest level */
                coarse_itsolver(A0, b0, e0, tol, prtlvl);
                
        }
        
    }
    
}

/**
 * \fn void nl_amli(AMG_data *mgl, AMG_param *param, INT level, INT num_levels)
 *
 * \brief Solve Ax=b with recursive nonlinear AMLI-cycle
 *
 * \param mgl         Pointer to AMG_data data
 * \param param       Pointer to AMG parameters
 * \param level       Current level
 * \param num_levels  Total number of levels
 *
 * \author Xiaozhe Hu
 * \date   04/06/2010
 *
 * \note Refer to Xiazhe Hu, Panayot S. Vassilevski, Jinchao Xu
 *       "Comparative Convergence Analysis of Nonlinear AMLI-cycle Multigrid", 2013.
 *
 */
void nl_amli (AMG_data *mgl,
              AMG_param *param,
              INT level,
              INT num_levels)
{    
    const SHORT  amg_type=param->AMG_type;
    const SHORT  prtlvl = param->print_level;
    const SHORT  smoother = param->smoother;
    const SHORT  coarse_solver = param->coarse_solver;
    const REAL   relax = param->relaxation;
    const INT    maxit = param->amli_degree+1;
    const REAL   tol = param->tol*1e-4;
    
    dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x;   // fine level b and x
    dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
    
    dCSRmat *A0 = &mgl[level].A;   // fine level matrix
    dCSRmat *A1 = &mgl[level+1].A; // coarse level matrix
    
    const INT m0 = A0->row, m1 = A1->row;
    
    REAL    *r = mgl[level].w.val;      // work array for residual
    
    dvector uH, bH;  // for coarse level correction
    uH.row = m1; uH.val = mgl[level+1].w.val + m1;
    bH.row = m1; bH.val = mgl[level+1].w.val + 2*m1;

    if ( prtlvl >= PRINT_MOST )
        printf("Nonlinear AMLI level %d, smoother %d.\n", num_levels, smoother);
    
    if ( level < num_levels-1 ) {
        
        // presmoothing
        dcsr_presmoothing(smoother,A0,b0,e0,param->presmooth_iter,
                          0,m0-1,1,relax);
        
        // form residual r = b - A x
        array_cp(m0,b0->val,r);
        dcsr_aAxpy(-1.0,A0,e0->val,r);
        
        // restriction r1 = R*r0
        switch (amg_type) {
            case UA_AMG:
                dcsr_mxv_agg(&mgl[level].R, r, b1->val);
                break;
            default:
                dcsr_mxv(&mgl[level].R, r, b1->val);
                break;
        }
        
        // call nonlinear AMLI-cycle recursively
        {
            dvec_set(m1,e1,0.0);
            
            // The coarsest problem is solved exactly.
            // No need to call krylov method on second coarest level
            if ( level == num_levels-2 ) {
                nl_amli(&mgl[level+1], param, 0, num_levels-1);
            }
            else { // recursively call preconditioned Krylov method on coarse grid
                
                precond_data pcdata;
                
                param_amg_to_prec(&pcdata, param);
                pcdata.maxit = 1;
                pcdata.max_levels = num_levels-1;
                pcdata.mgl_data = &mgl[level+1];
                
                precond pc;
                pc.data = &pcdata;
                pc.fct = precond_nl_amli;
                
                array_cp(m1, b1->val, bH.val);
                array_cp(m1, e1->val, uH.val);
                
                switch (param->nl_amli_krylov_type) {
                    case SOLVER_GCG: // Use GCG
                        dcsr_pgcg(A1,&bH,&uH,&pc,tol,maxit,1,PRINT_NONE);
                        break;                      
                    default: // Use FGMRES
                        dcsr_pvfgmres(A1,&bH,&uH,&pc,tol,maxit,30,1,PRINT_NONE);
                        break;
                }
                
                array_cp(m1, bH.val, b1->val);
                array_cp(m1, uH.val, e1->val);
            }
            
        }
        
        // prolongation e0 = e0 + P*e1
        switch (amg_type) {
            case UA_AMG:
                dcsr_aAxpy_agg(1.0, &mgl[level].P, e1->val, e0->val);
                break;
            default:
                dcsr_aAxpy(1.0, &mgl[level].P, e1->val, e0->val);
                break;
        }
        
        // postsmoothing
        dcsr_postsmoothing(smoother,A0,b0,e0,param->postsmooth_iter,
                           0,m0-1,-1,relax);
        
    }
    
    else { // coarsest level solver
        
        switch (coarse_solver) {
                
#if WITH_SUITESPARSE
            case SOLVER_UMFPACK:
                // use UMFPACK direct solver on the coarsest level //
                umfpack_solve(A0, b0, e0, mgl[level].Numeric, 0);
                break;
#endif
                
            default:
                /* use iterative solver on the coarsest level */
                coarse_itsolver(A0, b0, e0, tol, prtlvl);
                
        }
        
    }

}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
