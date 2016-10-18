/*! \file src/nonlinear/newton.c
 *
 *  Created by James Adler and Xiaozhe Hu on 10/18/16.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 */

/* This code will contain all the tools needed to perform Newton stepping */

#include "hazmat.h"

/******************************************************************************************************/
void initialize_newton(newton *n_it,input_param *inparam)
{

    // Number of Newton steps
    n_it->tsteps = inparam->nonlinear_itsolver_maxit;

    // Current Step
    n_it->current_step = 0;

    // Tolerances 1 - ||nonlinear residual||<tol OR 2 - ||update|| < tol OR 0 - BOTH
    n_it->tol_type = inparam->nonlinear_itsolver_toltype;
    n_it->tol      = inparam->nonlinear_itsolver_tol;

    // Step Length: sol = sol_prev + step_length*update
    n_it->step_length = 1.0; // Default is full Newton

    // Matrices and Vectors
    // Assume the form A(sol) = f gets linearized to
    // Jacobian(sol_prev)[update] = f - A(sol_prev)
    n_it->Jac=NULL;
    n_it->sol=NULL;
    n_it->sol_prev=malloc(sizeof(struct dvector));
    n_it->update=malloc(sizeof(struct dvector));
    n_it->rhs=NULL;     /* f - A(sol_prev) */

    return;
}
/******************************************************************************************************/

/****************************************************************************************/
void free_newton(newton* n_it)
{
    /* frees memory of arrays of newton struct */

    if(n_it->Jac) {
        dcsr_free(n_it->Jac);
        n_it->Jac=NULL;
    }

    if(n_it->sol) {
        dvec_free(n_it->sol);
        n_it->sol=NULL;
    }

    if(n_it->sol_prev) {
        dvec_free(n_it->sol_prev);
        free(n_it->sol_prev);
        n_it->sol_prev=NULL;
    }

    if(n_it->update) {
        dvec_free(n_it->update);
        free(n_it->update);
        n_it->rhs_update=NULL;
    }

    if(n_it->rhs) {
        dvec_free(n_it->rhs);
        n_it->rhs=NULL;
    }

    return;
}
/****************************************************************************************/

/******************************************************************************************************/
void update_newtonstep(newton* n_it)
{
    /* Updates the newton data at each step */

    // Counters
    n_it->current_step++;

    // Solution
    if(n_it->current_step==1) {
        dvec_alloc(n_it->sol->row,n_it->sol_prev);
    }
    dvec_cp(n_it->sol,n_it->sol_prev);

    // Update
    if(n_it->current_step==1) {
        dvec_alloc(n_it->sol->row,n_it->update);
    }
    dvec_set(n_it->update,0.0);

    return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void update_sol_newton(newton *n_it)
{
    /********* Updates solution to nonlinear problem *********************
   *
   *     sol = sol_prev + step_length * update
   *
   */

    dvec_axpyz(n_it->step_length,n_it->update,n_it->sol_prev,n_it->sol);

    return;
}
/******************************************************************************************************/

/******************************************************************************************************/
bool check_newton_convergence(newton *n_it,fespace* FE,trimesh* mesh, qcoordinates* cq)
{
    /********* Checks if Newton has converged *********************
   *
   *     If tol_type = 1: Check if || nonlinear residual (rhs)|| < tol
   *                   2: Check if || update || < tol
   *                   0: Check both
   *
   */

    bool newton_stop = false;
    REAL res_norm = L2norm(n_it->rhs->val,FE,mesh,cq);
    REAL update_norm = L2norm(n_it->rhs->val,FE,mesh,cq);

    if(n_it->current_step>=n_it->max_steps) {
        newton_stop=true;
        printf("The Newton iterations have reached the max number of iterations (%d Newton Steps) \n",n_it->current_step);
        printf("Convergence may not be reached.\n");
        printf("Final Nonlinear Residual = %25.16e\tLast Update Norm = %25.16e\n",res_norm,update_norm);
        return newton_stop;
    }

    switch (n_it->toltype)
    {
    default:
        if(res_norm<tol || update_norm<tol) {
            newton_stop=true;
            printf("Convergence met after %d Newton Steps.\n",n_it->current_step);
            printf("Final Nonlinear Residual = %25.16e\tLast Update Norm = %25.16e\n",res_norm,update_norm);
        }
        break;
    case 1:
        if(res_norm<tol) {
            newton_stop=true;
            printf("Convergence met after %d Newton Steps.\n",n_it->current_step);
            printf("Final Nonlinear Residual = %25.16e\tLast Update Norm = %25.16e\n",res_norm,update_norm);
        }
        break;
    case 2:
        if(update_norm<tol) {
            newton_stop=true;
            printf("Convergence met after %d Newton Steps.\n",n_it->current_step);
            printf("Final Nonlinear Residual = %25.16e\tLast Update Norm = %25.16e\n",res_norm,update_norm);
        }
        break;
    }

    return newton_stop;
}
/******************************************************************************************************/



