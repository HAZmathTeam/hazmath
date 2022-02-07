/*! \file examples/solver/solver_bsr.c
 *
 *  Created by Xiaozhe Hu on 02/06/22.
 *  Copyright 2022_HAZMATH__. All rights reserved.
 *
 * \brief This program read in a matrix (in BSR format) and a right hand side and solve it by certain linear solvers
 *
 * \note
 *
 */

/************* HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
/***********************************************************************/

/****** MAIN DRIVER **************************************************/
int main (int argc, char* argv[])
{
  printf("\n===========================================================================\n");
  printf("Beginning Program to solve a linear system in BSR format.\n");
  printf("===========================================================================\n");

  /* matrix and right hand side */
  dBSRmat Absr;
  dvector b, x;

  printf("\n===========================================================================\n");
  printf("Reading the matrix, right hand side, and parameters\n");
  printf("===========================================================================\n");
  /* set Parameters from Reading in Input File */
  input_param inparam;
  param_input_init(&inparam);
  param_input("./input.dat", &inparam);

  /* read the matrix and right hand side */
  // Read the stiffness matrix from bsrmat_SPE01.dat
  dbsr_read("bsrmat_SPE01.dat", &Absr);

  // Read the RHS from rhs_SPE01.dat
  dvector_read("rhs_SPE01.dat", &b);
  /************************************************************/

  /*************** ACTION *************************************/
  /* set initial guess */
  dvec_alloc(b.row, &x);
  dvec_set(b.row, &x, 0.0);

  /* Set Solver Parameters */
  INT solver_flag=-20;
  /* Set parameters for linear iterative methods */
  linear_itsolver_param linear_itparam;
  param_linear_solver_set(&linear_itparam, &inparam);
  /* Set parameters for algebriac multigrid methods */
  AMG_param amgparam;
  param_amg_init(&amgparam);
  param_amg_set(&amgparam, &inparam);
  param_amg_print(&amgparam);

  printf("\n===========================================================================\n");
  printf("Solving the linear system \n");
  printf("===========================================================================\n");

  // Use Krylov Iterative Solver
  // Diagonal preconditioner
  if (linear_itparam.linear_precond_type == PREC_DIAG) {
      solver_flag = linear_solver_dbsr_krylov_diag(&Absr, &b, &x, &linear_itparam);
  }
  // AMG preconditioner
  else if (linear_itparam.linear_precond_type == PREC_AMG){
      solver_flag = linear_solver_dbsr_krylov_amg(&Absr, &b, &x, &linear_itparam, &amgparam);
  }
  // No preconditoner
  else{
      solver_flag = linear_solver_dbsr_krylov(&Absr, &b, &x, &linear_itparam);
  }

  // Clean up memory
  dbsr_free(&Absr);
  dvec_free(&b);
  dvec_free(&x);

}	/* End of Program */
/*******************************************************************/
