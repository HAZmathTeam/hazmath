/*! \file src/nonlinear/fas.c
 *
 * \brief This code will contain all the tools needed to perform FAS cycles and its variants
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 02/23/21.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 */

#include "hazmath.h"

/******************************************************************************************************/
/*!
 * \fn void initialize_fas(fas_struct *fas_cyc,input_param *inparam,INT ndof)
 *
 * \brief Initialize the FAS struct for nonlinear iterations.
 *
 * \param inparam       Input from input parameter list
 * \param ndof          Number of DOF in system
 *
 * \return fas_cyc      Struct for FAS Cycles
 *
 */
void initialize_fas(fas_struct *fas_cyc,input_param *inparam,INT ndof)
{

  // Number of FAS steps max
  fas_cyc->max_steps = inparam->nonlinear_itsolver_maxit;

  // Current Step
  fas_cyc->current_step = 0;

  // Tolerances 1 - ||nonlinear residual||<tol
  fas_cyc->tol_type = inparam->nonlinear_itsolver_toltype;
  fas_cyc->tol      = inparam->nonlinear_itsolver_tol;

  // Smoothers
  fas_cyc->smooth_preits = inparam->fas_presmoothers;
  fas_cyc->smooth_postits = inparam->fas_postsmoothers;
  fas_cyc->smooth_tol = inparam->fas_smooth_tol;

  // Step Length: sol = sol_prev + step_length*update
  fas_cyc->step_length = 1.0; // Default is a full step

  // Matrices and Vectors
  fas_cyc->P=malloc(sizeof(struct dCSRmat));
  fas_cyc->R=malloc(sizeof(struct dCSRmat));
  fas_cyc->sol=malloc(sizeof(struct dvector));
  fas_cyc->nonlinear_res=malloc(sizeof(struct dvector));
  fas_cyc->res_norm=0.0;

  dvec_alloc(ndof,fas_cyc->sol);
  dvec_alloc(ndof,fas_cyc->nonlinear_res);

  return;
}
/******************************************************************************************************/

/****************************************************************************************/
/*!
 * \fn void free_fas(fas_struct* fas_cyc)
 *
 * \brief Frees memory of arrays of newton struct
 *
 * \return fas_cyc  Freed struct for FAS cycles
 *
 */
void free_fas(fas_struct* fas_cyc)
{
  if(fas_cyc->P) {
    dcsr_free(fas_cyc->P);
    free(fas_cyc->P);
    fas_cyc->P=NULL;
  }

  if(fas_cyc->R) {
    dcsr_free(fas_cyc->R);
    free(fas_cyc->R);
    fas_cyc->R=NULL;
  }

  if(fas_cyc->sol) {
    dvec_free(fas_cyc->sol);
    free(fas_cyc->sol);
    fas_cyc->sol=NULL;
  }

  if(fas_cyc->nonlinear_res) {
    dvec_free(fas_cyc->nonlinear_res);
    free(fas_cyc->nonlinear_res);
    fas_cyc->nonlinear_res=NULL;
  }

  return;
}
/****************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn INT check_fas_convergence(fas_struct *fas_cyc)
 *
 * \brief Checks if FAS has converged:
 *       For now just check:
 *       (||nonlinear residual (rhs)||_l2) / sqrt(length(rhs)) < tol
 *
 * \param fas_cyc     FAS struct
 *
 * \note If using the residual check, note that we compute the little l2 norms
 *       of the residual but then scale by the square root of the length of the
 *       residual vector, so that we are equivalent to the l_infinity norm of rhs
 *       and to scale for larger problems.
 */
INT check_fas_convergence(fas_struct *fas_cyc)
{

  INT fas_stop = 0;
  REAL tol = fas_cyc->tol;

  // Get norms
  REAL res_norm = fas_cyc->res_norm; // l2 norm of residual
  INT res_length = fas_cyc->nonlinear_res->row;
  REAL res_norm_scaled = res_norm / sqrt(res_length); // scaled norm of residual

  if(fas_cyc->current_step>=fas_cyc->max_steps) {
    fas_stop=1;
    printf("**********************************************************************************\n");
    printf("The FAS cycles have reached the max number of iterations (%d FAS Cycles) \n",fas_cyc->current_step);
    printf("Convergence may not be reached.\n");
    printf("\nl2-norm of Nonlinear Residual = %25.16e\n",res_norm);
    printf("               Scaled Version = %25.16e\n\n",res_norm_scaled);
    return fas_stop;
  }

  switch (fas_cyc->tol_type)
  {
  default:
    if(res_norm_scaled<tol) {
      fas_stop=1;
      printf("**********************************************************************************\n");
      printf("Convergence met after %d FAS Cycles.\n",fas_cyc->current_step);
      printf("\nl2-norm of Final Nonlinear Residual = %25.16e\n",res_norm);
      printf("                     Scaled Version = %25.16e\n\n",res_norm_scaled);
    }
    break;
  }

  return fas_stop;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void get_residual_norm(fas_struct* fas_cyc)
 *
 * \brief Computes the (little) l2 norm of the nonlinear residual (rhs).
 *
 * \param fas_cyc  FAS struct
 *
 * \note This uses the little l2 norm, though it probably should be a dual space
 *       norm for the entire block system instead.  One should avoid this as a
 *       stopping tolerance check regardless.
 */
void get_fasresidual_norm(fas_struct *fas_cyc)
{
  INT res_length = fas_cyc->nonlinear_res->row;
  INT i=0;
  REAL resnorm = 0.0;
  for(i=0;i<res_length;i++) {
    resnorm += fas_cyc->nonlinear_res->val[i]*fas_cyc->nonlinear_res->val[i];
  }
  fas_cyc->res_norm = sqrt(resnorm);
  return;
}
/******************************************************************************************************/
