/*! \file src/solver/mgcycle.c
 *
 *  Abstract multigrid cycle
 *
 *  Created by James Adler and Xiaozhe Hu on 12/25/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 */

#include "hazmath.h"
#include "mg_util.inl"

static SHORT krylov_cycle_dcsr_pgcg(dCSRmat *, dvector *, dvector *, precond *);
static SHORT krylov_cycle_dcsr_pgcr(dCSRmat *, dvector *, dvector *, precond *);

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
        
        // pre-smoothing with standard smoothers
        dcsr_presmoothing(smoother, &mgl[l].A, &mgl[l].b, &mgl[l].x,
                          param->presmooth_iter, 0, mgl[l].A.row-1, 1,
                          relax, ndeg, smooth_order, mgl[l].cfmark.val);
        
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
                           relax, ndeg, smooth_order, mgl[l].cfmark.val);
        
        if ( num_lvl[l] < cycle_type ) break;
        else num_lvl[l] = 0;
    }
    
    if ( l > 0 ) goto ForwardSweep;
    
    
}


void amli (AMG_data *mgl,
                       AMG_param *param,
                       INT level)
{
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
    
    const SHORT  amg_type=param->AMG_type;
    const SHORT  prtlvl = param->print_level;
    const SHORT  smoother = param->smoother;
    const SHORT  smooth_order = param->smooth_order;
    const SHORT  coarse_solver = param->coarse_solver;
    const SHORT  degree= param->amli_degree;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol*1e-4;
    const SHORT  ndeg = param->polynomial_degree;
    
    // local variables
    REAL   alpha  = 1.0;
    REAL * coef   = param->amli_coef;
    
    dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x;   // fine level b and x
    dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
    
    dCSRmat *A0 = &mgl[level].A;   // fine level matrix
    dCSRmat *A1 = &mgl[level+1].A; // coarse level matrix
    
    const INT m0 = A0->row, m1 = A1->row;
    
    INT      *ordering = mgl[level].cfmark.val; // smoother ordering
    REAL     *r        = mgl[level].w.val;      // work array for residual
    REAL     *r1       = mgl[level+1].w.val+m1; // work array for residual  
    
    if ( prtlvl >= PRINT_MOST )
        printf("AMLI level %d, smoother %d.\n", level, smoother);
    
    if ( level < mgl[level].num_levels-1 ) {
                
        // presmoothing
        dcsr_presmoothing(smoother,A0,b0,e0,param->presmooth_iter,
                          0,m0-1,1,relax,ndeg,smooth_order,ordering);
        
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
                           0,m0-1,-1,relax,ndeg,smooth_order,ordering);

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


void nl_amli (AMG_data *mgl,
              AMG_param *param,
              INT level,
              INT num_levels)
{
    /**
     * \fn void nl_amli (AMG_data *mgl, AMG_param *param,
     *                               INT level, INT num_levels)
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
    
    const SHORT  amg_type=param->AMG_type;
    const SHORT  prtlvl = param->print_level;
    const SHORT  smoother = param->smoother;
    const SHORT  smooth_order = param->smooth_order;
    const SHORT  coarse_solver = param->coarse_solver;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol*1e-4;
    const SHORT  ndeg = param->polynomial_degree;
    
    dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x;   // fine level b and x
    dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
    
    dCSRmat *A0 = &mgl[level].A;   // fine level matrix
    dCSRmat *A1 = &mgl[level+1].A; // coarse level matrix
    
    const INT m0 = A0->row, m1 = A1->row;
    
    INT      *ordering = mgl[level].cfmark.val; // smoother ordering
    REAL     *r        = mgl[level].w.val;      // work array for residual
    
    dvector uH;  // for coarse level correction
    uH.row = m1; uH.val = mgl[level+1].w.val + m1;

    if ( prtlvl >= PRINT_MOST )
        printf("Nonlinear AMLI level %d, smoother %d.\n", num_levels, smoother);
    
    if ( level < num_levels-1 ) {
        
        // presmoothing
        dcsr_presmoothing(smoother,A0,b0,e0,param->presmooth_iter,
                          0,m0-1,1,relax,ndeg,smooth_order,ordering);
        
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
            
            // V-cycle will be enforced when needed !!!
            if ( mgl[level+1].cycle_type <= 1 ) {
                
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
                
                array_cp (m1, e1->val, uH.val);
                
                switch (param->nl_amli_krylov_type) {
                        
                    case SOLVER_GCG: // Use GCG
                        krylov_cycle_dcsr_pgcg(A1,b1,&uH,&pc);
                        break;
                        
                    default: // Use GCR
                        krylov_cycle_dcsr_pgcr(A1,b1,&uH,&pc);
                        break;
                }
                
                array_cp (m1, uH.val, e1->val);
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
                           0,m0-1,-1,relax,ndeg,smooth_order,ordering);
        
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
/*--     Private Functions       --*/
/*---------------------------------*/

static SHORT krylov_cycle_dcsr_pgcg (dCSRmat *A,
                                     dvector *b,
                                     dvector *u,
                                     precond *pc)
{
    /**
     * \fn static SHORT krylov_cycle_dcsr_pgcg (dCSRmat *A, dvector *b,
     *                                               dvector *u, precond *pc)
     *
     * \brief A preconditioned GCR method for solving Au=b
     *
     * \param *A    Pointer to the coefficient matrix
     * \param *b    Pointer to the dvector of right hand side
     * \param *u    Pointer to the dvector of DOFs
     * \param *pre  Pointer to the structure of precondition (precond)
     *
     * \author Zheng Li, Chensong Zhang
     * \date   11/09/2014
     *
     * \note   Specified for unsmoothed aggregation cycle
     */
    
    REAL   absres, relres, normb;
    REAL   alpha1, alpha2, gamma1, gamma2, rho1, rho2, beta1, beta2, beta3, beta4;
    REAL   *work, *r, *x1, *v1, *v2;
    
    INT    m=A->row;
    REAL   *x = u->val;
    
    // allocate temp memory
    work = (REAL *)calloc(4*m,sizeof(REAL));
    r = work; x1 = r + m; v1 = r + 2*m; v2 = r + 3*m;
    
    normb = array_norm2(m, b->val);
    
    array_cp(m, b->val, r);
    
    // Preconditioning
    if (pc != NULL)
        pc->fct(r, x, pc->data);
    else
        array_cp(m, r, x);
    
    // v1 = A*p
    dcsr_mxv(A, x, v1);
    
    // rho1 = (p,v1)
    rho1 = array_dotprod (m, x, v1);
    
    // alpha1 = (p, r)
    alpha1 = array_dotprod (m, x, r);
    
    beta1 = alpha1/rho1;
    
    // r = r - beta1 *v1
    array_axpy(m, -beta1, v1, r);
    
    // norm(r)
    absres = array_norm2(m, r);
    
    // compute relative residual
    relres = absres/normb;
    
    // if relres reaches tol(0.2), pgcg will stop,
    // otherwise, another one pgcg iteration will do.
    if (relres < 0.2) {
        array_ax(m, beta1, x);
        free(work);
        return SUCCESS;
    }
    
    // Preconditioning
    if (pc != NULL)
        pc->fct(r, x1, pc->data);
    else
        array_cp(m, r, x1);
    
    //v2 = A*p
    dcsr_mxv(A, x1, v2);
    
    //gamma0 = (x1,v1)
    gamma1 = array_dotprod (m, x1, v1);
    
    //alpha2 = (x1,r)
    alpha2  = array_dotprod(m, x1, r);
    
    //rho2 = (x1,v2)
    rho2 = array_dotprod(m, x1, v2);
    
    gamma2 = gamma1;
    
    beta2 = rho2 - gamma1*gamma2/rho1;
    beta3 = (alpha1 - gamma2*alpha2/beta2)/rho1;
    beta4 = alpha2/beta2;
    
    array_ax(m, beta3, x);
    
    array_axpy(m, beta4, x1, x);
    
    // free
    free(work);
    return SUCCESS;
}

static SHORT krylov_cycle_dcsr_pgcr (dCSRmat *A,
                                     dvector *b,
                                     dvector *u,
                                     precond *pc)
{
    /**
     * \fn static SHORT krylov_cycle_dcsr_pgcr (dCSRmat *A, dvector *b,
     *                                               dvector *u, precond *pc)
     *
     * \brief A preconditioned GCR method for solving Au=b
     *
     * \param *A    Pointer to the coefficient matrix
     * \param *b    Pointer to the dvector of right hand side
     * \param *u    Pointer to the dvector of DOFs
     * \param *pre  Pointer to the structure of precondition (precond)
     *
     * \author zheng Li, Chensong Zhang
     * \date   11/09/2014
     *
     * \note   Specified for unsmoothed aggregation cycle.
     */

    
    REAL   absres = BIGREAL;
    REAL   relres  = BIGREAL, normb  = BIGREAL;
    REAL   alpha, alpha1, alpha2, alpha3, alpha4, beta, gamma, rho1, rho2;
    
    INT    m=A->row;
    REAL   *x = u->val;
    
    // allocate temp memory
    REAL *work, *r, *x1, *v1, *v2;
    work = (REAL *)calloc(4*m,sizeof(REAL));
    r = work; x1 = r + m; v1 = r + 2*m; v2 = r + 3*m;
    
    normb=array_norm2(m, b->val);
    array_cp(m, b->val, r);
    
    // Preconditioning
    if (pc != NULL)
        pc->fct(r, x, pc->data);
    else
        array_cp(m, r, x);
    
    // v1 = A*x
    dcsr_mxv(A, x, v1);
    // rho1 = (v1,v1)
    rho1 = array_dotprod (m, v1, v1);
    // alpha1 = (r, v1)
    alpha1 = array_dotprod (m, v1, r);
    
    alpha = alpha1/rho1;
    
    // r = r - alpha *v1
    array_axpy(m, -alpha, v1, r);
    
    // norm(r)
    absres = array_norm2(m, r);
    
    // compute relative residual
    relres = absres/normb;
    
    // if relres reaches tol(0.2), pgcr will stop,
    // otherwise, another one pgcr iteration will do.
    if (relres < 0.2) {
        array_ax(m, alpha, x);
        free(work);
        return SUCCESS;
    }
    
    // Preconditioning
    if (pc != NULL)
        pc->fct(r, x1, pc->data);
    else
        array_cp(m, r, x1);
    
    //v2 = A*x1
    dcsr_mxv(A, x1, v2);
    
    //gamma = (v1,v2)
    gamma = array_dotprod (m, v1, v2);
    //beta = (v2,v2)
    beta  = array_dotprod(m, v2, v2);
    //alpha2 = (r,v2)
    alpha2 = array_dotprod(m, r, v2);
    
    rho2 = beta - gamma*gamma/rho1;
    
    alpha3 = alpha1/rho1 - gamma*alpha2/(rho1*rho2);
    
    alpha4 = alpha2/rho2;
    
    // x = alpha3*x + alpha4*x1
    array_ax(m, alpha3, x);
    array_axpy(m, alpha4, x1, x);
    
    // free
    free(work);
    return SUCCESS;
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
