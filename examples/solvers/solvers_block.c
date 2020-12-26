/*! \file examples/solver/solver_block.c
 *
 *  Created by Xiaozhe Hu on 01/01/19.
 *  Copyright 2019_HAZMATH__. All rights reserved.
 *
 * \brief This program read in a (block) matrix and a right hand side and solve it by certain linear solvers
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
  block_dCSRmat A;
  dvector b;
  dvector x;

  dCSRmat L;
  dCSRmat M;
  dCSRmat Lf;

  dCSRmat *A_diag;

  printf("\n===========================================================================\n");
  printf("Reading the matrix, right hand side, and parameters\n");
  printf("===========================================================================\n");
  /* read in matrix and right hand side from Simula test problem */
  bdcsr_alloc(2, 2, &A);
  dcoo_read_dcsr("simula/n32/A11.dat", A.blocks[0]);
  dcoo_read_dcsr("simula/n32/A12.dat", A.blocks[1]);
  dcoo_read_dcsr("simula/n32/A21.dat", A.blocks[2]);
  dcoo_read_dcsr("simula/n32/A22.dat", A.blocks[3]);
  dvector_read("simula/n32/b.dat", &b);

  dcoo_read_dcsr("simula/n32/As.dat", &L);  // Laplacian matrix
  dcoo_read_dcsr("simula/n32/Ms.dat", &M);  // Mass matrix
  //dcoo_read_dcsr("simula/n16/Af.dat", &Lf); // Fractional Laplacian matrix

  /* set diagonal blocks for block diagonal preconditioner */
  //A_diag = (dCSRmat *)calloc(2, sizeof(dCSRmat));
  // first block
  //dcsr_alloc(A.blocks[0]->row, A.blocks[0]->col, A.blocks[0]->nnz, &A_diag[0]);
  //dcsr_cp(A.blocks[0], &(A_diag[0]));
  // second block
  //dcsr_alloc(Lf.row, Lf.col, Lf.nnz, &A_diag[1]);
  //dcsr_cp(&Lf, &(A_diag[1]));

  /* set parameters for fractiaon problem */
  REAL s_frac_power = -0.5;
  REAL t_frac_power = 0.5;
  REAL alpha = 1.;
  REAL beta = 0.;
  REAL scaling_a = 1./8.00; // 1./8.0052
  REAL scaling_m = 512.;     // 32.

  //-----------------------------------------------------
  /* set Parameters from Reading in Input File */
  input_param inparam;
  param_input_init(&inparam);
  param_input("./input_block.dat", &inparam);

  /* set Solver Parameters */
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
  /************************************************************/
  /* set initial guess */
  dvec_alloc(b.row,&x); //  same as x=dvec_create_p(A->row);
  dvec_set(x.row, &x, 0.0);


  // block preconditioner
  if (linear_itparam.linear_precond_type >0 && linear_itparam.linear_precond_type <40) {
      //solver_flag = linear_solver_bdcsr_krylov_block(&A, &b, &x, &linear_itparam, &amgparam, A_diag);
      solver_flag = linear_solver_bdcsr_babuska_block_2(&A, &b, &x, &L, &M, &linear_itparam, &amgparam, &amgparam, s_frac_power, t_frac_power, alpha, beta, scaling_a, scaling_m);
  }
  else
  // no preconditioner
  {
      solver_flag = linear_solver_bdcsr_krylov(&A, &b, &x, &linear_itparam);
  }


  //dvector_write("x.dat", &x);

  // Clean up memory
  bdcsr_free(&A);
  dvec_free(&b);
  dvec_free(&x);

  //dcsr_free(&A_diag[0]);
  //dcsr_free(&A_diag[1]);
  //if(A_diag) free(A_diag);

  dcsr_free(&L);
  dcsr_free(&M);
  //dcsr_free(&Lf);


}	/* End of Program */
/*******************************************************************/
