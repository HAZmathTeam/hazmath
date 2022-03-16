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
#include "solvers_frac_help.h"
/***********************************************************************/
/*MAIN DRIVER*/
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
  dvector poles_r,poles_i;
  dvector residues_r,residues_i;

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
    fprintf(stderr,"\n***USING THE DEFAULTS:\n\t\t\t%s A3.dat b3.dat M3.dat",argv[0]);
    fprintf(stderr,  "\n=========================================================\n\n");
    fnamea=strdup("./A3.dat");
    fnameb=strdup("./b3.dat");
    fnamem=strdup("./M3.dat");
    read_to_eof=0;
  } else {
    fnamea=strdup(argv[1]);
    fnameb=strdup(argv[2]);
    fnamem=strdup(argv[3]);
  }
  if(read_to_eof){
    fprintf(stdout,"\n%s: reading file \"%s\" until EOF\n", __FUNCTION__,fnamea);
    fp = fopen(fnamea,"r");
    if (!fp) check_error(ERROR_OPEN_FILE, __FUNCTION__);
    if(fnamea) free(fnamea);
    A=dcoo_read_eof_dcsr_p(fp,NULL);
    fclose(fp);
    fprintf(stdout,"\n%s: reading file \"%s\" until EOF\n", __FUNCTION__,fnameb);
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
  }
  /************************************************************/
  /* compute the rational approximation
  /************************************************************/
  // parameters used in the function
  /* 
     approximates inverse(alpha*x^s + beta*x^t), that is it
     preconditions (alpha*x^s + beta*x^t).
  */
  REAL s = FUNC_PARAM_S, t = FUNC_PARAM_T, \
    alpha = FUNC_PARAM_ALPHA,		   \
    beta = FUNC_PARAM_BETA;
  REAL16 func_param[4]={-1e20,-1e20,-1e20,-1e20};// not needed
  INT i;
  REAL xmin_in=0e0, xmax_in=1e0;  // interval for x
  INT numval=(1<<14)+1; // how many init points for the AAA.
  INT mmax_in=(INT )(numval/2);
  mmax_in=MAX_NUMBER_POLES;//
  // parameters used in the AAA algorithm
  REAL16 AAA_tol=powl(2e0,TOL_POWER);  // tolerance of the AAA algorithm
  INT k=-22; // k is the number of poles in the final interpolation after tolerance is achieved or mmax is reached.
  INT print_level=0; // print level
  ////////////////////////////////////////////////////////////////////////
  // scaling so that the interval for M^{-1}A is [0,1]
  //
  //  REAL sa = 1./8.0052; //(1/l1-norm of A)
  //  REAL sm = 32.; // should be inverse(min_diag_M).
  REAL sum0,sa = -1e0;
  REAL sm = -1e0; // should be inverse(min_diag_M).
  INT ij,j;
  sm=1e20;
  for (i=0;i<M->row;++i){    
    for (ij=M->IA[i];ij<M->IA[i+1];++ij){
      j=M->JA[ij];
      if(i!=j)continue;
      if(sm > M->val[ij]){
	sm=M->val[ij];
      }
    }
  }
  sa=-1e0;
  for (i=0;i<A->row;++i){
    sum0=0e0;
    for (ij=A->IA[i];ij<A->IA[i+1];++ij){
      j=A->JA[ij];
      sum0+=fabs(A->val[ij]);
    }
    if(sa<sum0) sa=sum0;
}
  REAL ddd=2;
  sm=1/sm/ddd/(ddd+1);
  sa=1e0/sa;
  fprintf(stdout,"\nscaling params:sm=%.16e,sa=%.16e\n",sm,sa);
  REAL alphas = alpha, betas = beta;  
  // scale alpha = alpha*sa^(-s)*sm^(s-1)
  alphas = alpha*pow(sa, -s)*pow(sm, s-1.);
  // scale beta = beta*sa^(-t)*sm^(t-1)
  betas  = beta*pow(sa, -t)*pow(sm, t-1.);
  // assign scaled parameters for the function
  func_param[0] = s;
  func_param[1] = t;
  if (alphas > betas) {
    func_param[2] = 1.;
    func_param[3] = betas/alphas;
  }else{
    func_param[2] = alphas/betas;
    func_param[3] = 1.;
  }
  // rational approximation arrays:
  REAL **rpzwf=malloc(7*sizeof(REAL *));
  // evaluate the function at many many points on the interval [0,1]
  REAL16 **zf=set_f_values(f_to_approx_local,		\
			   func_param[0],func_param[1],	\
			   func_param[2],func_param[3],	\
			   &numval,xmin_in,xmax_in,	\
			   print_level);
  // compute the rational approximation using AAA algorithm
  // rpzwf[0]-> real part of residues (resr[]); (also as last entry
  //            contains the free term resr[m]);  
  // rpzwf[1]-> imaginary part of residues (resi[]); (also as last
  //            entry contains the free term resi[m]);
  // rpzwf[2]-> real(poles) (polr[]);
  // rpzwf[3]-> imag(poles) (poli[]);
  // the rational approximation is:
  // r(z)=res[k-1] + \sum_{i=0}^{k-2} res[i]/(z-pol[i]);
  //
  REAL err_max=get_rpzwf(numval,zf[0],zf[1],	\
			 rpzwf,			\
			 &mmax_in, &k,		\
			 AAA_tol,print_level);
  dvec_alloc(k,  &residues_r);
  dvec_alloc(k,  &residues_i);
  dvec_alloc(k-1, &poles_r);
  dvec_alloc(k-1, &poles_i);
  array_cp(k, rpzwf[0], residues_r.val);
  array_cp(k, rpzwf[1], residues_i.val);
  array_cp(k-1, rpzwf[2], poles_r.val);
  array_cp(k-1, rpzwf[3], poles_i.val);
  // print
  /* for(i=0;i<k;i++)printf("\n res(%d)=%.16e + %.16e * i;",i+1,*(rpzwf[0]+i),*(rpzwf[1]+i)); */
  /* for(i=0;i<k-1;i++)printf("\n pol(%d)=%.16e  + %.16e * i;",i+1,*(rpzwf[2]+i),*(rpzwf[3]+i)); */
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
  if (alphas > betas) {
    dvec_ax(1./alphas, &bs);
  } else {
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
  //  solver_flag = linear_solver_frac_rational_approx(&As, &bs, x, &Ms, &amgparam, &linear_itparam, &poles, &residues);
  solver_flag = linear_solver_frac_rational_approx(&As,&bs,x,&Ms,	\
						   &amgparam,&linear_itparam, \
						   &poles_r, &poles_i,	\
						   &residues_r,&residues_i);
  fprintf(stdout,"\nim=sqrt(-1);\n");
  for (i=0;i<residues_r.row;++i){
    fprintf(stdout,"\nresidues(%d,1:1)=%.18e + im*(%.18e);",i+1,residues_r.val[i],residues_i.val[i]);
  }
  for (i=0;i<poles_r.row;++i){
    fprintf(stdout,"\npoles(%d,1:1)=%.18e+im*(%.18e);",i+1,poles_r.val[i],poles_i.val[i]);
  }
  fprintf(stdout,"\ns=%.16Le, t=%.16Le, alpha=%.16Le, beta=%.16Le\n",func_param[0],func_param[1],func_param[2],func_param[3]);
  dvector_write("x.dat", x);
  // Clean up memory
  free(A);
  free(b);
  free(x);
  dcsr_free(&As);
  dcsr_free(&Ms);
  dvec_free(&bs);
  dvec_free(&poles_r);
  dvec_free(&residues_r);
  dvec_free(&poles_i);
  dvec_free(&residues_i);
}	/* End of Program */
/*******************************************************************/
