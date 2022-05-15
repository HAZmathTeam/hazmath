/*! \file examples/solver/solver_bdcsr.c
 *
 *  Created by Xiaozhe Hu on 04/16/22.
 *  Copyright 2022_HAZMATH__. All rights reserved.
 *
 * \brief This program read in a matrix (in block_dCSRmat format) and a right hand side and solve it by certain linear solvers
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
  block_dCSRmat A;
  dvector b, x;

  dCSRmat *A_diag;

  INT i;

  printf("\n===========================================================================\n");
  printf("Reading the matrix, right hand side, and parameters\n");
  printf("===========================================================================\n");
  /* set Parameters from Reading in Input File */
  input_param inparam;
  param_input_init(&inparam);
  param_input("./input.dat", &inparam);

  /* read the matrix and right hand side (2 by 2) */
  INT brow = 2;
  INT bcol = 2;
  //bdcsr_alloc_minimal(brow, bcol, &A);
  bdcsr_alloc(brow, bcol, &A);
  A_diag = (dCSRmat *)calloc(brow, sizeof(dCSRmat));

  // Read the 00 block of the stiffness matrix
  dcoo_read_dcsr("A11.dat", A.blocks[0]);
  dcoo_read_dcsr("A12.dat", A.blocks[1]);
  dcoo_read_dcsr("A21.dat", A.blocks[2]);
  dcoo_read_dcsr("A22.dat", A.blocks[3]);

  // Read the RHS from rhs_SPE01.dat
  dvector_read("Arhs.dat", &b);
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

  // --------------------------------------------------------------------------------------------
  // Set diagonal blocks for AMG solver.  Coarsening is based on the blocks in A_diag.
  // They can be diagonal blocks of the block matrix A or approximations to the Schur complements
  // --------------------------------------------------------------------------------------------
  // Use first diagonal block directly in A_diag
  dcsr_alloc(A.blocks[0]->row, A.blocks[0]->col, A.blocks[0]->nnz, &A_diag[0]);
  dcsr_cp(A.blocks[0], &A_diag[0]);

  // Use approximated Schur complement of the second diagonal block in A_diag
  dvector diag_M;
  dCSRmat invM = dcsr_create(A.blocks[0]->row,A.blocks[0]->row,A.blocks[0]->row);
  dcsr_getdiag(A.blocks[0]->row, A.blocks[0], &diag_M);
  for (i=0;i<A.blocks[0]->row;i++)
  {
      invM.IA[i] = i;
      invM.JA[i] = i;
      if (diag_M.val[i] > SMALLREAL) invM.val[i]   = 1.0/diag_M.val[i];
      else invM.val[i] = 1.0;
  }
  invM.IA[A.blocks[0]->row] = A.blocks[0]->row;
  dCSRmat BTB;
  dcsr_rap(A.blocks[2], &invM, A.blocks[1], &BTB);
  dcsr_add(&BTB, 1.0, A.blocks[3], -1.0, &A_diag[1]);

  // Use Krylov Iterative Solver
  if (linear_itparam.linear_precond_type == PREC_AMG){
      solver_flag = linear_solver_bdcsr_krylov_amg(&A, &b, &x, &linear_itparam, &amgparam, A_diag);
  }
  // No preconditoner
  else{
       solver_flag = linear_solver_bdcsr_krylov(&A, &b, &x, &linear_itparam);
  }

  // Clean up memory
  bdcsr_free(&A);
  dvec_free(&b);
  dvec_free(&x);

  for (i=0; i<brow; i++) dcsr_free(&A_diag[i]);
  if (A_diag) free(A_diag);

  dvec_free(&diag_M);
  dcsr_free(&invM);
  dcsr_free(&BTB);

}	/* End of Program */
/*******************************************************************/
