/*! \file examples/Solver/eigen.c
 *
 *  Created by Xiaozhe Hu on 01/01/19.
 *  Copyright 2019_HAZMATH__. All rights reserved.
 *
 * \brief This program read a symmetric matrix A and a SPD B and
 *         solves the eigenvalue problem A*x=lambda*B*x using LAPACK.
 *
 * \note
 *
 */

/************* HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
/****** MAIN DRIVER **************************************************/
int main (int argc, char* argv[])
{
  dCSRmat *A=NULL;
  dCSRmat *B=NULL;
  printf("\n%%===========================================================================\n");
  printf("%%READING A and B matrix\n");
  printf("%%===========================================================================\n");
  FILE *fp=NULL;
  char *fnamea=NULL,*fnameb=NULL;
  fnamea=strdup("A20.dat");
  fnameb=strdup("B20.dat");
  fp = fopen(fnamea,"r");
  if (!fp) check_error(ERROR_OPEN_FILE, __FUNCTION__);
  if(fnamea) free(fnamea);
  A=dcoo_read_dcsr_p(fp);
  fclose(fp);
  fp = fopen(fnameb,"r");
  if (!fp) check_error(ERROR_OPEN_FILE, __FUNCTION__);
  if(fnameb) free(fnameb);
  B=dcoo_read_dcsr_p(fp);
  /*************** ACTION *************************************/
  REAL *evalues=calloc(A->row,sizeof(REAL));
  //  only eigenvalues:
  //  eigsymm(A,B,evalues,NULL);
  //  print_full_mat(1,A->row,evalues,"ev");
  REAL *evectors=calloc(A->row*A->row,sizeof(REAL));
  // both eigenvalues and eigenvectors
  eigsymm(A,B,evalues,evectors);
  print_full_mat(A->row,1,evalues,"evals");
  print_full_mat(A->row,A->col,evectors,"evecs");
  // Clean up memory
  free(A);
  free(B);
  free(evalues);
  free(evectors);
}	/* End of Program */
/*******************************************************************/

