/*! \file src/nonlinear/newton.c
 *
 * \brief This code will contain all the tools needed to perform Newton stepping
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 10/18/16.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 * \note modified by James Adler 11/11/2016
 */

#include "hazmath.h"

/******************************************************************************************************/
/*!
 * \fn void initialize_newton(newton *n_it,input_param *inparam,INT ndof,INT blksize)
 *
 * \brief Initialize the Newton struct for nonlinear iterations.
 *
 * \param inparam       Input from input parameter list
 * \param ndof          Number of DOF in system
 * \param blksize       Block size of if block matrices (assumes square)
 *
 * \return n_it         Struct for Newton Iterations
 *
 */
void initialize_newton(newton *n_it,input_param *inparam,INT ndof,INT blksize)
{

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
  if(n_it->isblock) {// In block form the Jacobian is a block_dCSRmat
    n_it->Jac=NULL;
    n_it->Jac_block=malloc(sizeof(struct block_dCSRmat));
    bdcsr_alloc(blksize,blksize,n_it->Jac_block);
  } else { //The Jacobian is a regular dCSRmat
    n_it->Jac=malloc(sizeof(struct dCSRmat));
    n_it->Jac_block=NULL;
  }
  n_it->sol=malloc(sizeof(struct dvector));
  // dvector sol = dvec_create(ndof);
  // n_it->sol = &sol;
  n_it->sol_prev=malloc(sizeof(struct dvector));
  n_it->update=malloc(sizeof(struct dvector));
  n_it->rhs=malloc(sizeof(struct dvector));     /* f - A(sol_prev) */
  n_it->res_norm=0.0;
  n_it->update_norm=0.0;

  dvec_alloc(ndof,n_it->sol);
  dvec_alloc(ndof,n_it->rhs);
  dvec_alloc(ndof,n_it->sol_prev);
  dvec_alloc(ndof,n_it->update);

  return;
}
/******************************************************************************************************/

/****************************************************************************************/
/*!
 * \fn void free_newton(newton* n_it)
 *
 * \brief Frees memory of arrays of newton struct
 *
 * \return n_it         Freed struct for Newton Iterations
 *
 */
void free_newton(newton* n_it)
{
  if(n_it->Jac) {
    dcsr_free(n_it->Jac);
    free(n_it->Jac);
    n_it->Jac=NULL;
  }

  if(n_it->Jac_block) {
    bdcsr_free(n_it->Jac_block);
    free(n_it->Jac_block);
    n_it->Jac_block=NULL;
  }

  if(n_it->sol) {
    dvec_free(n_it->sol);
    free(n_it->sol);
    n_it->sol=NULL;
  }

  if(n_it->rhs) {
    dvec_free(n_it->rhs);
    free(n_it->rhs);
    n_it->rhs=NULL;
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

  // if(n_it->current_step>0) {
  //   if(n_it->sol_prev) {
  //     dvec_free(n_it->sol_prev);
  //     free(n_it->sol_prev);
  //     n_it->sol_prev=NULL;
  //   }

  //   if(n_it->update) {
  //     dvec_free(n_it->update);
  //     free(n_it->update);
  //     n_it->update=NULL;
  //   }
  // } else {
  //   if(n_it->sol_prev) {
  //     free(n_it->sol_prev);
  //     n_it->sol_prev=NULL;
  //   }
  //   if(n_it->update) {
  //     free(n_it->update);
  //     n_it->update=NULL;
  //   }
  // }

  return;
}
/****************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void update_newtonstep(newton* n_it)
 *
 * \brief Updates the Newton data at each step.
 *
 * \return n_it     Updated Newton struct
 *
 */
void update_newtonstep(newton* n_it)
{
  // Counters
  n_it->current_step++;

  // Solution
  // if(n_it->current_step==1) {
  //   dvec_alloc(n_it->sol->row,n_it->sol_prev);
  // }
  dvec_cp(n_it->sol,n_it->sol_prev);

  // Update
  // if(n_it->current_step==1) {
  //   dvec_alloc(n_it->sol->row,n_it->update);
  // }
  dvec_set(n_it->update->row,n_it->update,0.0);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void update_sol_newton(newton *n_it)
 *
 * \brief Updates the solution to the nonlinear problem.
 *        sol = sol_prev + step_length * update
 *
 * \return n_it.sol     Updated Newton solution
 *
 */
void update_sol_newton(newton *n_it)
{
  dvec_axpyz(n_it->step_length,n_it->update,n_it->sol_prev,n_it->sol);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn INT check_newton_convergence(newton *n_it)
 *
 * \brief Checks if Newton has converged:
 *        If tol_type = 1: Check if ||update||_L2 < tol
 *                      2: Check if (||nonlinear residual (rhs)||_l2) / sqrt(length(rhs)) < tol
 *                      3: Check both 1 and 2 (either)
 *
 * \param n_it     Newton struct
 *
 * \note If using the residual check, note that we compute the little l2 norms
 *       of the residual but then scale by the square root of the length of the
 *       residual vector, so that we are equivalent to the l_infinity norm of rhs
 *       and to scale for larger problems.
 */
INT check_newton_convergence(newton *n_it)
{

  INT newton_stop = 0;
  REAL tol = n_it->tol;

  // Get norms
  REAL res_norm = n_it->res_norm; // l2 norm of residual
  INT res_length = n_it->rhs->row;
  REAL res_norm_scaled = res_norm / sqrt(res_length); // scaled norm of residual
  REAL update_norm = n_it->update_norm; // L2 norm of update

  if(n_it->current_step>=n_it->max_steps) {
    newton_stop=1;
    printf("**********************************************************************************\n");
    printf("The Newton iterations have reached the max number of iterations (%lld Newton Steps) \n",(long long )n_it->current_step);
    printf("Convergence may not be reached.\n");
    printf("\nl2-norm of Nonlinear Residual = %25.16e\n",res_norm);
    printf("               Scaled Version = %25.16e\n\n",res_norm_scaled);
    printf("L2 Norm of Update             = %25.16e\n\n",update_norm);
    return newton_stop;
  }

  switch (n_it->tol_type)
  {
  default:
    if(update_norm<tol) {
      newton_stop=1;
      printf("**********************************************************************************\n");
      printf("Convergence met after %lld Newton Steps.\n",(long long )n_it->current_step);
      printf("\nl2-norm of Final Nonlinear Residual = %25.16e\n",res_norm);
      printf("                     Scaled Version = %25.16e\n\n",res_norm_scaled);
      printf("Final L2 Norm of Update             = %25.16e\n\n",update_norm);
    }
    break;
  case 2:
    if(res_norm_scaled<tol) {
      newton_stop=1;
      printf("**********************************************************************************\n");
      printf("Convergence met after %lld Newton Steps.\n",(long long )n_it->current_step);
      printf("\nl2-norm of Final Nonlinear Residual = %25.16e\n",res_norm);
      printf("                     Scaled Version = %25.16e\n\n",res_norm_scaled);
      printf("Final L2 Norm of Update             = %25.16e\n\n",update_norm);
    }
    break;
  case 3:
    if(res_norm_scaled<tol || update_norm<tol) {
      newton_stop=1;
      printf("**********************************************************************************\n");
      printf("Convergence met after %lld Newton Steps.\n",(long long )n_it->current_step);
      printf("\nl2-norm of Final Nonlinear Residual = %25.16e\n",res_norm);
      printf("                     Scaled Version = %25.16e\n\n",res_norm_scaled);
      printf("Final L2 Norm of Update             = %25.16e\n\n",update_norm);
    }
    break;
  }

  return newton_stop;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void get_residual_norm(newton *n_it)
 *
 * \brief Computes the (little) l2 norm of the nonlinear residual (rhs).
 *
 * \note This uses the little l2 norm, though it probably should be a dual space
 *       norm for the entire block system instead.  One should avoid this as a
 *       stopping tolerance check regardless.
 */
void get_residual_norm(newton *n_it)
{
  INT res_length = n_it->rhs->row;
  INT i=0;
  REAL resnorm = 0.0;
  for(i=0;i<res_length;i++) {
    resnorm += n_it->rhs->val[i]*n_it->rhs->val[i];
  }
  n_it->res_norm = sqrt(resnorm);
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void get_update_norm(newton *n_it,fespace* FE,mesh_struct* mesh, qcoordinates* cq)
 *
 * \brief Computes the L2 norm of the update.
 *
 * \param n_it     Newton struct
 * \param FE       FE space
 * \param mesh     Mesh struct
 * \param cq       Quadrature for computing norms
 *
 * \note Assumes non-block form
 */
void get_update_norm(newton *n_it,fespace* FE,mesh_struct* mesh, qcoordinates* cq)
{
  n_it->update_norm = L2norm(n_it->update->val,FE,mesh,cq);
  return;
}
/******************************************************************************************************/
// OUTDATED
/******************************************************************************************************/
// /*!
//  * \fn void get_blockresidual_norm(newton *n_it,block_fespace* FE,mesh_struct* mesh, qcoordinates* cq)
//  *
//  * \brief Computes the norms of the nonlinear residual (rhs) and the update.
//  *
//  * \param n_it     Newton struct
//  * \param FE       Block FE space
//  * \param mesh     Mesh struct
//  * \param cq       Quadrature for computing norms
//  *
//  * \note Assumes block form (and we store the total (combined) norm)
//  */
// void get_blockresidual_norm(newton *n_it,block_fespace* FE,mesh_struct* mesh, qcoordinates* cq)
// {
//   REAL* res_norm = (REAL *) calloc(FE->nspaces,sizeof(REAL));
//   L2norm_block(res_norm,n_it->rhs->val,FE,mesh,cq);
//
//   REAL total_res_norm = 0;
//   INT i;
//   for(i=0;i<FE->nspaces;i++) {
//     total_res_norm += res_norm[i]*res_norm[i];
//   }
//   n_it->res_norm = sqrt(total_res_norm);
//
//   if(res_norm) free(res_norm);
//   return;
// }
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void get_blockupdate_norm(newton *n_it,block_fespace* FE,mesh_struct* mesh, qcoordinates* cq)
 *
 * \brief Computes the L2 norm the update.
 *
 * \param n_it     Newton struct
 * \param FE       Block FE space
 * \param mesh     Mesh struct
 * \param cq       Quadrature for computing norms
 *
 * \note Assumes block form (and we store the total (combined) norm)
 */
void get_blockupdate_norm(newton *n_it,block_fespace* FE,mesh_struct* mesh, qcoordinates* cq)
{
  REAL* update_norm = (REAL *) calloc(FE->nspaces,sizeof(REAL));
  L2norm_block(update_norm,n_it->update->val,FE,mesh,cq);

  REAL total_update_norm = 0;
  INT i;
  for(i=0;i<FE->nspaces;i++) {
    total_update_norm += update_norm[i]*update_norm[i];
  }
  n_it->update_norm = sqrt(total_update_norm);

  if(update_norm) free(update_norm);
  return;
}
/******************************************************************************************************/
