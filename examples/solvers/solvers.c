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
// local include (temp)
// if READ_TO_EOF is 0 the first record in the input files is m,n, nnz and n for dvector
//otherwise the number of elements in (i,j,v) or dvector is found by reading until EOF.
#ifndef READ_TO_EOF
#define READ_TO_EOF 1
#endif
/****** MAIN DRIVER **************************************************/
int main (int argc, char* argv[])
{
  printf("\n===========================================================================\n");
  printf("Beginning Program to solve a linear system.\n");
  printf("===========================================================================\n");

  /* matrix and right hand side */
  dCSRmat *A=NULL;
  dvector *b=NULL;
  dvector *x=NULL;

  printf("\n===========================================================================\n");
  printf("Reading the matrix, right hand side, and parameters\n");
  printf("===========================================================================\n");

  /* set Parameters from Reading in Input File */
  input_param inparam;
  param_input_init(&inparam);
  param_input("./input.dat", &inparam);
  /* read the matrix and right hand side */
  SHORT read_to_eof=READ_TO_EOF;
  FILE *fp=NULL;
  char *fnamea=NULL,*fnameb=NULL;
  if(argc<3){
    fprintf(stderr,"\n\n=========================================================\n\n");
    fprintf(stderr,"***ERROR: %s called with wrong number of arguments!!!\n",argv[0]);
    fprintf(stderr,"Usage: %s filename_with_MATRIX(I,J,V) filename_with_RHS\n",argv[0]);
    fprintf(stderr,"\n***USING THE DEFAULTS:\n\t\t\t%s A.dat b.dat",argv[0]);
    fprintf(stderr,  "\n=========================================================\n\n");
    fnamea=strdup("A.dat");
    fnameb=strdup("b.dat");
    read_to_eof=0;
//    exit(129);
  } else {
    fnamea=strdup(argv[1]);
    fnameb=strdup(argv[2]);
  }
  if(read_to_eof){
    fprintf(stdout,"\n%s: reading file \"%s\" unitl EOF\n", __FUNCTION__,fnamea);
    fp = fopen(fnamea,"r");
    if (!fp) check_error(ERROR_OPEN_FILE, __FUNCTION__);
    if(fnamea) free(fnamea);
    A=dcoo_read_eof_dcsr_p(fp,NULL);
    fclose(fp);
    fprintf(stdout,"\n%s: reading file \"%s\" unitl EOF\n", __FUNCTION__,fnameb);
    fp = fopen(fnameb,"r");
    if (!fp) check_error(ERROR_OPEN_FILE, __FUNCTION__);
    if(fnameb) free(fnameb);
    b=dvector_read_eof_p(fp);
  } else {
    fp = fopen(fnamea,"r");
    if (!fp) check_error(ERROR_OPEN_FILE, __FUNCTION__);
    if(fnamea) free(fnamea);
    A=dcoo_read_dcsr_p(fp);
    fclose(fp);
    fp = fopen(fnameb,"r");
    if (!fp) check_error(ERROR_OPEN_FILE, __FUNCTION__);
    if(fnameb) free(fnameb);
    b=dvector_read_p(fp);
  }
  /************************************************************/
  /*************** ACTION *************************************/
  /* set initial guess */
  dvec_alloc_p(A->row,&x); //  same as x=dvec_create_p(A->row);
  dvec_set(x->row, x, 0.0);

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
    solver_flag = linear_solver_amg(A, b, x, &amgparam);
  } else { // Use Krylov Iterative Solver
    // Diagonal preconditioner
    if (linear_itparam.linear_precond_type == PREC_DIAG) {
        solver_flag = linear_solver_dcsr_krylov_diag(A, b, x, &linear_itparam);
    }
    // AMG preconditioner
    else if (linear_itparam.linear_precond_type == PREC_AMG){
        solver_flag = linear_solver_dcsr_krylov_amg(A, b, x, &linear_itparam, &amgparam);
    }
    // No preconditoner
    else{
        solver_flag = linear_solver_dcsr_krylov(A, b, x, &linear_itparam);
    }
  }

  // Clean up memory
  free(A);
  free(b);
  free(x);
}	/* End of Program */
/*******************************************************************/

