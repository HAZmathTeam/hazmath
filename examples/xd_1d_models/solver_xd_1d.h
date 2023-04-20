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
INT solver_xd_1d(const char *finput_solver,const char *dir_matrices)
{
  /*************** ACTION *************************************/
  //  char *dir_matrices=strdup("./input/1d_matrices_2d/");
  //  char *finput_solver=strdup("./input/solver.input");
  /* Set Solver Parameters */
  input_param inparam;
  dCSRmat A; //
  dvector b,x;
  ivector idofs;
  //
  read_and_setup(finput_solver,dir_matrices,&inparam,&A,&b,&x,&idofs);
  //
  //  INT num_iters=-20;
  /* Set parameters for linear iterative methods */
  linear_itsolver_param linear_itparam;
  param_linear_solver_set(&linear_itparam, &inparam);
  /* Set parameters for algebriac multigrid methods */
  AMG_param amgparam;
  param_amg_init(&amgparam);
  param_amg_set(&amgparam, &inparam);
  param_amg_print(&amgparam);

  fprintf(stdout,"\n===========================================================================\n");
  fprintf(stdout,"Solving the linear system \n");
  fprintf(stdout,"===========================================================================\n");
  // --------------------------------------------------------------------------------------------
  // Set diagonal blocks for AMG solver.  Coarsening is based on the blocks in AD.
  // They can be diagonal blocks of the block matrix A or approximations to the Schur complements
  // --------------------------------------------------------------------------------------------
  if (linear_itparam.linear_precond_type == 16 ){
    linear_solver_dcsr_krylov_metric_amg(&A, &b, &x, &idofs, &linear_itparam, &amgparam);
  }
  else if (linear_itparam.linear_precond_type == PREC_AMG){
    linear_solver_dcsr_krylov_amg(&A, &b, &x, &linear_itparam, &amgparam);
  }
  // No preconditoner
  else{
    linear_itparam.linear_precond_type=0;
    linear_solver_dcsr_krylov(&A, &b, &x, &linear_itparam);
  }
  //  fprintf(stdout,"\nIters=%lld;Preconditioner=%lld\n\n",(long long )num_iters,(long long )linear_itparam.linear_precond_type);
  
  char *fsolution  = fname_set("output/","sol.txt");
  dvec_write(fsolution,&x);
  free(fsolution);
  dvec_free(&b);
  dvec_free(&x);
  dcsr_free(&A);
  return 0;
}
/*******************************************************************/
