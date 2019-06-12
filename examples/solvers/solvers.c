/*! \file examples/Solver/Solver.c
 *
 *  Created by Xiaozhe Hu on 01/01/19.
 *  Copyright 2019_HAZMATH__. All rights reserved.
 *
 * \brief This program read in a matrix and a right hand side and solve it by certain linear solvers
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
  printf("Beginning Program to solve a linear system.\n");
  printf("===========================================================================\n");

  /* matrix and right hand side */
  dCSRmat A;
  dvector b;
  dvector x;

  printf("\n===========================================================================\n");
  printf("Reading the matrix, right hand side, and parameters\n");
  printf("===========================================================================\n");

  /* read the matrix and right hand side */
  dcoo_read_dcsr("A.dat", &A);
  dvector_read("b.dat", &b);

  /* set Parameters from Reading in Input File */
  input_param inparam;
  param_input_init(&inparam);
  param_input("./input.dat", &inparam);

  /* set initial guess */
  dvec_alloc(A.row, &x);
  dvec_set(x.row, &x, 0.0);

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

  // Use AMG as iterative solver
  if (linear_itparam.linear_itsolver_type == SOLVER_AMG){
    solver_flag = linear_solver_amg(&A, &b, &x, &amgparam);
  } else { // Use Krylov Iterative Solver
    // Diagonal preconditioner
    if (linear_itparam.linear_precond_type == PREC_DIAG) {
        solver_flag = linear_solver_dcsr_krylov_diag(&A, &b, &x, &linear_itparam);
    }
    // AMG preconditioner
    else if (linear_itparam.linear_precond_type == PREC_AMG){
        solver_flag = linear_solver_dcsr_krylov_amg(&A, &b, &x, &linear_itparam, &amgparam);
    }
    // No preconditoner
    else{
        solver_flag = linear_solver_dcsr_krylov(&A, &b, &x, &linear_itparam);
    }
  }

  // Clean up memory
  dcsr_free(&A);
  dvec_free(&b);

}	/* End of Program */
/*******************************************************************/
