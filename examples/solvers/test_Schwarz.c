/*! \file examples/solver/test_Schwarz.c
 *
 *  Created by Ludmil Zikatanov on 20230116.
 *  Copyright 2022_HAZMATH__. All rights reserved.
 *
 * \brief This program reads in a matrix and solves it by the Schwarz method. 
 *
 * \note
 *
 */

/************* HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
/***********************************************************************/
// performs maxiter schwarz iterations.
// none of the arguments should be changed here except the approximation x!
/****** MAIN DRIVER **************************************************/
int main (int argc, char* argv[])
{
  /* matrix and right hand side */
  fprintf(stdout,"\n===========================================================================\n");
  fprintf(stdout,"Program to test multiplicative Schwarz method for a linear system in CSR format.\n");
  fprintf(stdout,"===========================================================================\n");
  fprintf(stdout,"\n===========================================================================\n");
  fprintf(stdout,"Reading the matrix, right hand side, and parameters\n");
  fprintf(stdout,"===========================================================================\n");
  /* set Parameters from Reading in Input File */
  input_param inparam;
  param_input_init(&inparam);
  param_input("./input_Schwarz.dat", &inparam);
  /* Set Solver Parameters */
  INT solver_flag;
  /* Set parameters for linear iterative methods */
  linear_itsolver_param linear_itparam;
  param_linear_solver_set(&linear_itparam, &inparam);
  REAL tol=linear_itparam.linear_tol;
  /* Set parameters for algebriac multigrid methods */
  AMG_param amgparam;
  param_amg_init(&amgparam);
  param_amg_set(&amgparam, &inparam);
  param_amg_print(&amgparam);
  Schwarz_param swzparam;
  Schwarz_data  Schwarz;
  //
  fprintf(stdout,"\n===========================================================================\n");
  fprintf(stdout,"Solving the linear system \n");
  fprintf(stdout,"===========================================================================\n");

  // --------------------------------------------------------------------------------------------
    // initialize the schwarz parameters:
    param_Schwarz_init(&swzparam);
    // Use Iterative Solver   
    param_linear_solver_print(&linear_itparam);
    //
    param_amg_print(&amgparam);
    //    if (amgparam.Schwarz_levels > 0 ) {
    swzparam.Schwarz_mmsize = amgparam.Schwarz_mmsize;
    swzparam.Schwarz_maxlvl = amgparam.Schwarz_maxlvl;
    swzparam.Schwarz_type   = amgparam.Schwarz_type;
    swzparam.Schwarz_blksolver = amgparam.Schwarz_blksolver;
      //    }
    Schwarz.Schwarz_type = swzparam.Schwarz_type;
    Schwarz.swzparam=&swzparam;
    //
    param_Schwarz_print(&swzparam);
  /************************************************************/
    dCSRmat *A=calloc(1,sizeof(dCSRmat));
    dcoo_read_dcsr("A.dat",A); //'A' is for ascii.
    fprintf(stdout,"\n%%%% Matrix stat: rows(A)=%d,cols(A)=%d,nnz(A)=%d\n",A->row,A->col,A->nnz); 
    //
    dvector b, x;
    b=dvec_create(A->row);
    x=dvec_create(A->row);    
    // Set the same matrix for Schwarz:
    Schwarz.A=A[0];
    // random RHS
    dvec_rand(b.row,&b);
    // x0=0
    memset(x.val,0,x.row*sizeof(REAL));
    INT iter=0,niter=1000;
    REAL relres=1e10,rnrm0=1e20,rnrm=1e20;
    rnrm0=residual_norm(A,&x,&b);
    ivector *seeds=calloc(1,sizeof(ivector));
    seeds[0]=ivec_create(0);
    seeds=NULL;
    /* seeds[0]=ivec_create(Schwarz.A.row); */
    /* for(i=0;i<seeds->row;++i){ */
    /*   seeds->val[i]=i; */
    /* } */
    Schwarz_setup(&Schwarz,&swzparam,seeds);
    while(relres>tol && iter <niter){
      rnrm=residual_norm(A,&x,&b);
      relres=rnrm/rnrm0;
      fprintf(stdout,"\nres(%d)=%.5e; relres(%d)=%.5e;",iter+1,rnrm,iter+1,relres);
      ++iter;
      smoother_dcsr_Schwarz(&Schwarz,&x,&b,1);
    }
    relres=rnrm/rnrm0;
    fprintf(stdout,"\nres(%d)=%.5e; relres(%d)=%.5e;",iter+1,rnrm,iter+1,rnrm/rnrm0);
    // end iterations
    fprintf(stdout,"\n%%%% end iters\n");
  /* Clean up memory */
  dvec_free(&b);
  dvec_free(&x);
  ivec_free(seeds);
  Schwarz.A=dcsr_create(0,0,0);
  schwarz_data_free(&Schwarz);  
  dcsr_free(A);// A is a void pointer so we can free it.
  free(A);
  free(seeds);
}
/*******************************************************************/
/* End of Program */

