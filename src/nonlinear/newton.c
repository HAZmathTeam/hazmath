/*! \file src/nonlinear/newton.c
 *
 * \brief This code will contain all the tools needed to perform Newton stepping
 *
 *  Created by James Adler and Xiaozhe Hu on 10/18/16.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 * \note modified by James Adler 11/11/2016
 */

#include "hazmath.h"

/******************************************************************************************************/
void initialize_newton(newton *n_it,input_param *inparam)
{
  /*!
   * \fn void initialize_newton(newton *n_it,input_param *inparam)
   *
   * \brief Initialize the Newton struct for nonlinear iterations.
   *
   * \param inparam       Input from input parameter list
   *
   * \return n_it         Struct for Newton Iterations
   *
   */

    // Number of Newton steps
    n_it->max_steps = inparam->nonlinear_itsolver_maxit;

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
    n_it->Jac_block=NULL;
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
  /*!
   * \fn void free_newton(newton* n_it)
   *
   * \brief Frees memory of arrays of newton struct
   *
   * \return n_it         Freed struct for Newton Iterations
   *
   */

    if(n_it->Jac) {
        dcsr_free(n_it->Jac);
        n_it->Jac=NULL;
    }

    if(n_it->Jac_block) {
        bdcsr_free(n_it->Jac_block);
        n_it->Jac_block=NULL;
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
        n_it->update=NULL;
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
  /*!
   * \fn void update_newtonstep(newton* n_it)
   *
   * \brief Updates the Newton data at each step.
   *
   * \return n_it     Updated Newton struct
   *
   */

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
    dvec_set(n_it->update->row,n_it->update,0.0);

    return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void update_sol_newton(newton *n_it)
{
  /*!
   * \fn void update_sol_newton(newton *n_it)
   *
   * \brief Updates the solution to the nonlinear problem.
   *        sol = sol_prev + step_length * update
   *
   * \return n_it.sol     Updated Newton solution
   *
   */

    dvec_axpyz(n_it->step_length,n_it->update,n_it->sol_prev,n_it->sol);

    return;
}
/******************************************************************************************************/

/******************************************************************************************************/
int check_newton_convergence(newton *n_it,fespace* FE,trimesh* mesh, qcoordinates* cq)
{
  /*!
   * \fn int check_newton_convergence(newton *n_it,fespace* FE,trimesh* mesh, qcoordinates* cq)
   *
   * \brief Checks if Newton has converged:
   *        If tol_type = 1: Check if ||nonlinear residual (rhs)|| < tol
   *                      2: Check if ||update|| < tol
   *                      0: Check both
   *
   * \param n_it     Newton struct
   * \param FE       FE space
   * \param mesh     Mesh struct
   * \param cq       Quadrature for computing norms
   *
   */

    int newton_stop = 0;
    REAL tol = n_it->tol;
    REAL res_norm = L2norm(n_it->rhs->val,FE,mesh,cq);
    REAL update_norm = L2norm(n_it->update->val,FE,mesh,cq);

    if(n_it->current_step>=n_it->max_steps) {
        newton_stop=1;
        printf("The Newton iterations have reached the max number of iterations (%d Newton Steps) \n",n_it->current_step);
        printf("Convergence may not be reached.\n");
        printf("Final Nonlinear Residual = %25.16e\tLast Update Norm = %25.16e\n",res_norm,update_norm);
        return newton_stop;
    }

    switch (n_it->tol_type)
    {
    default:
        if(res_norm<tol || update_norm<tol) {
            newton_stop=1;
            printf("Convergence met after %d Newton Steps.\n",n_it->current_step);
            printf("Final Nonlinear Residual = %25.16e\tLast Update Norm = %25.16e\n",res_norm,update_norm);
        }
        break;
    case 1:
        if(res_norm<tol) {
            newton_stop=1;
            printf("Convergence met after %d Newton Steps.\n",n_it->current_step);
            printf("Final Nonlinear Residual = %25.16e\tLast Update Norm = %25.16e\n",res_norm,update_norm);
        }
        break;
    case 2:
        if(update_norm<tol) {
            newton_stop=1;
            printf("Convergence met after %d Newton Steps.\n",n_it->current_step);
            printf("Final Nonlinear Residual = %25.16e\tLast Update Norm = %25.16e\n",res_norm,update_norm);
        }
        break;
    }

    return newton_stop;
}
/******************************************************************************************************/



