/*! \file examples/solver/solver_frac.c
 *
 *  Created by Xiaozhe Hu on 09/27/20.
 *  Copyright 2020_HAZMATH__. All rights reserved.
 *
 * \brief This program solve fractional problem using rational approximation
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


/************************************************************************/
/*function to be approximated*/
/************************************************************************/
REAL16 f_to_approx(REAL16 x,void *param)
{
  REAL16 *s,s1,s2,alpha,beta;
  if(param!=NULL){
    s=(REAL16 *)param;
    s1=s[0];
    s2=s[1];
    alpha=s[2];
    beta=s[3];
  } else {
    s1=0.5e0;
    s2=-0.5e0;
    alpha=1e0;
    beta=2e0;
  }
  //  fprintf(stdout,"\ns1=%Lf; s2=%Lf; alpha=%Lf; beta=%Lf;",s1,s2,alpha,beta);
  return 1./(alpha*powl(x,s1)+beta*powl(x,s2));
}
/**/



/****** MAIN DRIVER **************************************************/
int main (int argc, char* argv[])
{
  printf("\n===========================================================================\n");
  printf("Beginning Program to solve a linear system.\n");
  printf("===========================================================================\n");

  /* matrix and right hand side */
  dCSRmat *A=NULL;
  dCSRmat *M=NULL;
  dvector *b=NULL;
  dvector *x=NULL;

  /* poles and residues of rational approximation */
  dvector poles;
  dvector residues;

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
  char *fnamea=NULL,*fnameb=NULL, *fnamem; //*fpoles=NULL,*fresidues=NULL;
  if(argc<3){
    fprintf(stderr,"\n\n=========================================================\n\n");
    fprintf(stderr,"***ERROR: %s called with wrong number of arguments!!!\n",argv[0]);
    fprintf(stderr,"Usage: %s filename_with_MATRIX(I,J,V) filename_with_RHS\n",argv[0]);
    fprintf(stderr,"\n***USING THE DEFAULTS:\n\t\t\t%s As.dat rd.dat Ms.dat",argv[0]);
    fprintf(stderr,  "\n=========================================================\n\n");
    fnamea=strdup("simula/As.dat");
    fnameb=strdup("simula/rd.dat");
    fnamem=strdup("simula/Ms.dat");
    //fpoles=strdup("poles.dat");
    //fresidues=strdup("residues.dat");
    read_to_eof=0;
//    exit(129);
  } else {
    fnamea=strdup(argv[1]);
    fnameb=strdup(argv[2]);
    fnamem=strdup(argv[3]);
    //fpoles=strdup(argv[3]);
    //fresidues=strdup(argv[4]);
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
    fclose(fp);
    fprintf(stdout,"\n%s: reading file \"%s\" unitl EOF\n", __FUNCTION__,fnamea);
    fp = fopen(fnamem,"r");
    if (!fp) check_error(ERROR_OPEN_FILE, __FUNCTION__);
    if(fnamem) free(fnamem);
    M=dcoo_read_eof_dcsr_p(fp,NULL);
    fclose(fp);
    //fprintf(stdout,"\n%s: reading file \"%s\" unitl EOF\n", __FUNCTION__,fpoles);
    //fp = fopen(fpoles,"r");
    //if (!fp) check_error(ERROR_OPEN_FILE, __FUNCTION__);
    //if(fpoles) free(fpoles);
    //poles=dvector_read_eof_p(fp);
    //fclose(fp);
    //fprintf(stdout,"\n%s: reading file \"%s\" unitl EOF\n", __FUNCTION__,fresidues);
    //fp = fopen(fresidues,"r");
    //if (!fp) check_error(ERROR_OPEN_FILE, __FUNCTION__);
    //if(fresidues) free(fresidues);
    //residues=dvector_read_eof_p(fp);
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
    fclose(fp);
    fp = fopen(fnamem,"r");
    if (!fp) check_error(ERROR_OPEN_FILE, __FUNCTION__);
    if(fnamem) free(fnamem);
    M=dcoo_read_dcsr_p(fp);
    fclose(fp);
    //fp = fopen(fpoles,"r");
    //if (!fp) check_error(ERROR_OPEN_FILE, __FUNCTION__);
    //if(fpoles) free(fpoles);
    //poles=dvector_read_p(fp);
    //fclose(fp);
    //fp = fopen(fresidues,"r");
    //if (!fp) check_error(ERROR_OPEN_FILE, __FUNCTION__);
    //if(fresidues) free(fresidues);
    //residues=dvector_read_p(fp);
  }
  /************************************************************/

  //dvector_print(stdout, poles);
  //dvector_print(stdout, residues);

  /************************************************************/
  /* compute the rational approximation
  /************************************************************/
  INT i;
  // parameters used in the function
  REAL16 func_param[4]={-0.5e0,0.5e0,1e0,0e0};  //original s, t, alpha,beta
  REAL s = -0.5, t = 0.5, alpha = 0.0, beta = 1.0;
  REAL alphas = alpha, betas = beta;

  // parameters used in the AAA algorithm
  REAL xmin_in=0e0, xmax_in=1e0;  // interval for x
  INT mbig=(1<<14)+1;  // initial number of points on the interval [x_min, x_max]
  INT mmax_in=(INT )(mbig/2);  // maximal number of pole + 1
  REAL16 AAA_tol=powl(2e0,-52e0);  // tolerance of the AAA algorithm
  INT k=-22; // k is the number of nodes in the final interpolation after tolerance is achieved or mmax is reached.
  INT print_level=0; // print level

  // output of the AAA algorithm
  REAL **rpnwf=malloc(5*sizeof(REAL *));  // output of the AAA algorithm.  It contains residues, poles, nodes, weights, function values

  // scaling so that the interval for x is [0,1]
  REAL sa = 1./8.0052;
  REAL sm = 32.;

  // scale alpha = alpha*sa^(-s)*sm^(s-1)
  alphas = alpha*pow(sa, -s)*pow(sm, s-1.);
  // scale beta = beta*sa^(-t)*sm^(t-1)
  betas  = beta*pow(sa, -t)*pow(sm, t-1.);

  // assign scaled parameters for the function
  func_param[0] = s;
  func_param[1] = t;
  if (alphas > betas)
  {
    func_param[2] = 1.;
    func_param[3] = betas/alphas;
  }
  else
  {
    func_param[2] = alphas/betas;
    func_param[3] = 1.;
  }

  // compute the rational approximation using AAA algorithms
  REAL err_max=get_cpzwf(f_to_approx, (void *)func_param,	rpnwf, &mbig, &mmax_in, &k, xmin_in, xmax_in, AAA_tol, print_level);

  // assign poles and residules
  dvec_alloc(k,  &residues);
  dvec_alloc(k-1, &poles);
  array_cp(k, rpnwf[0], residues.val);
  array_cp(k-1, rpnwf[1], poles.val);

  /*
  // print
  dvector_print(stdout, &poles);
  dvector_print(stdout, &residues);
  for(i=0;i<k;i++)printf("\n res(%d)=%.16e;",i+1,*(rpnwf[0]+i));
  for(i=0;i<k-1;i++)printf("\n pol(%d)=%.16e;",i+1,*(rpnwf[1]+i));
  */

  /************************************************************/

  /************************************************************/
  /* scale the whole problem
  /************************************************************/
  dCSRmat As = dcsr_create(A->row, A->col, A->nnz);
  dCSRmat Ms = dcsr_create(M->row, M->col, M->nnz);
  dvector bs = dvec_create(b->row);

  dcsr_cp(A, &As);
  dcsr_axm(&As, sa);

  dcsr_cp(M, &Ms);
  dcsr_axm(&Ms, sm);

  dvec_cp(b, &bs);
  if (alphas > betas)
  {
    dvec_ax(1./alphas, &bs);
  }
  else
  {
    dvec_ax(1./betas, &bs);
  }

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
  printf("Solving the fractional system (alpha*A^s + beta*A^t)x=b \n");
  printf("===========================================================================\n");

  // use rational approximation
  //solver_flag = linear_solver_dcsr_krylov_amg(A, b, x, &linear_itparam, &amgparam);
  solver_flag = linear_solver_frac_rational_approx(&As, &bs, x, &Ms, &amgparam, &linear_itparam, &poles, &residues);

  dvector_write("x.dat", x);

  // Clean up memory
  free(A);
  free(b);
  free(x);
  dcsr_free(&As);
  dcsr_free(&Ms);
  dvec_free(&bs);
  dvec_free(&poles);
  dvec_free(&residues);
}	/* End of Program */
/*******************************************************************/
