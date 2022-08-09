#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "hazmath.h"
/*=====================================================================*/
/**
 * \fn static void dcsr_trilu_diag(const SHORT is_sym, dCSRmat *a, 
 *                          dCSRmat *al, dvector *adiag)
 *
 * \brief extracting the lower triangle by columns; upper triangle by
 *        rows and the diagonal of of a dcsr matrix (done in-place,
 *        a is overwritten by its upper triangle; also removes all
 *        zero elements).
 *
 * \param is_sym if the matrix is symmetric, then
 *               AL->row=A->col=A->nnz=0 and all pointers in al are
 *               free and NULL;
 *
 * \param *a pointer to the dCSRmat matrix (on output is the STRICT
 *           upper trianngle by rows).
 *
 * \param *adiag pointer to a dvector with the diagonal elements.
 *
 * \param *al pointer to the dCSRmat matrix: strict lower triangle
 *            stored by columns.
 *
 *
 * \return void; a and al are modified; a has its strict upper
 *               triangle; adiag contains the diagonal of a
 *
 * \author Ludmil
 * \date 20220802
 */
static void dcsr_trilu_diag(const SHORT is_sym,			\
		     dCSRmat *a,dCSRmat *al, dvector *adiag)
{
  INT k,j,kj,kj0,kj1;
  // transpose a: al should be allocated upon entry here. no permutation
  if(is_sym){
    al->row=0; al->col=0; al->nnz=0;
    if(al->IA) free(al->IA);
    if(al->JA) free(al->JA);
    if(al->val) free(al->val);
    al->IA=NULL; al->JA=NULL; al->val=NULL;
  } else {
    //AL here should not be null!
    dcsr_transz(a,NULL,al);
    // now grab the upper triangle in al which is the lower triangle of
    // a by columns.
    al->nnz=al->IA[0];
    for(k=0;k<al->row;k++){
      kj0=al->IA[k];
      kj1=al->IA[k+1];
      al->IA[k]=al->nnz;
      for(kj=kj0;kj<kj1;kj++){
	j=al->JA[kj];
	// skip even the diagonal here:
	if( ((k-j)>=0) || (al->val[kj]==0.) ) continue; // lower is rowind>=colind
	//fprintf(stdout,"\nlu=%d;(k,j)=(%d,%d):%.5e",lu,k,j,al->val[kj]);fflush(stdout);
	//      fprintf(stdout,"\nlu=%d;(k,j)=(%d,%d):%.5e",lu,k,j,al->val[kj]);fflush(stdout);
	al->JA[al->nnz]=j;
	al->val[al->nnz]=al->val[kj];
	al->nnz++;
      }
    }
    al->IA[al->row]=al->nnz;
    //fprintf(stdout,"nnnza=%d\n",al->nnz);fflush(stdout);
    if(al->nnz>0){
      al->JA=realloc(al->JA,al->nnz*sizeof(INT));
      al->val=realloc(al->val,al->nnz*sizeof(REAL));
    }else{
      free(al->JA);
      free(al->val);
      al->JA=NULL;
      al->val=NULL;
    }
  }
  /****************************************************************/
  // now upper triangle and the diagonal (if is_sym<>0 we come here too):
  a->nnz=a->IA[0];
  for(k=0;k<a->row;k++){
    kj0=a->IA[k];
    kj1=a->IA[k+1];
    a->IA[k]=a->nnz;
    for(kj=kj0;kj<kj1;kj++){
      j=a->JA[kj];
      if(k==j) {
	adiag->val[k]=a->val[kj];
	continue;
      }
      //      if(k<=j) for upper; (k>=j) for lower;
      if( ((k-j)>0) || (a->val[kj]==0.) ) continue; // lower is rowind>=colind
      //      fprintf(stdout,"\nlu=%d;(k,j)=(%d,%d):%.5e",lu,k,j,a->val[kj]);fflush(stdout);
      //      fprintf(stdout,"\nlu=%d;(k,j)=(%d,%d):%.5e",lu,k,j,a->val[kj]);fflush(stdout);
      a->JA[a->nnz]=j;
      a->val[a->nnz]=a->val[kj];
      a->nnz++;
    }
  }
  a->IA[a->row]=a->nnz;
  //  fprintf(stdout,"nnnza_upper=%d\n",a->nnz);fflush(stdout);
  if(a->nnz>0){
    a->JA=realloc(a->JA,a->nnz*sizeof(INT));
    a->val=realloc(a->val,a->nnz*sizeof(REAL));
  }else{
    free(a->JA);
    free(a->val);
    a->JA=NULL;
    a->val=NULL;
  }
  return;
}
/**/
/*=====================================================================*/
/**
 * \fn void symbolic_factor_symm(dCSRmat *AU, dCSRmat *U)
 *
 * \brief Symbolic factorization of a sparse matrix with a symmetric
 *        pattern. Gives the structure of the upper triagle (U) of the
 *        LU factorization; as input we have only the upper triangle
 *        of A and the diagonal of A.
 *
 * \param *AU pointer to the upper triangle of the dCSR matrix
 *
 * \param *U pointer to a dCSRmat: upper triangle of dvector containing the nonzero structure of the factor. U->val is NULL here (not needed). 
 *
 * \return void; *U is modified. it should be allocated as a dCSRmat
 *               struct upon entry here. The members of the struct are
 *               allocated here except U->val.
 *
 * \note Reference: Chapter 7 in Sergio Pissanetzky. Sparse matrix
 *                  technology. Academic Press Inc. [Harcourt Brace
 *                  Jovanovich Publishers], London, 1984.
 *
 * \note fortan77 (19960101); C99 (20220801)
 *
 * \author Ludmil 
 * \date 20220802
 *
 */
void symbolic_factor_symm(dCSRmat *AU, dCSRmat *U)
{
  INT k,kp1,i,j,jk,last_u;
  INT n=AU->row;
  INT *ichn=(INT *)calloc(n,sizeof(INT));
  //  setup U
  U->row=AU->row;
  U->col=AU->col;
  U->nnz=AU->nnz;
  U->IA=calloc(U->row+1,sizeof(INT));
  if(U->row<U->nnz){
    U->JA=calloc(U->nnz,sizeof(INT));
    last_u=U->nnz-1;
  }else {
    U->JA=calloc(U->row,sizeof(INT));
    last_u=U->row-1;
  }
  //
  INT nm1 = n - 1;
  for(i=0;i<n;i++){
    U->IA[i] = -1;
    ichn[i] = -1;
  }
  INT nzp = 0,nzpi,nzpp,minx,qend;
  INT ia_strt,ia_end,it_strt,it_end;
  SHORT flag;
  for(i=0;i<nm1;++i){
    nzpi = nzp;
    nzpp = n + nzp - i;
    minx = n; //
    ia_strt = AU->IA[i];
    ia_end = AU->IA[i+1];
    if(ia_end > ia_strt) {
      for(jk = ia_strt;jk<ia_end;++jk){
	j = AU->JA[jk];
	if(nzp>last_u){ 
	  U->JA=realloc(U->JA,(nzp+1)*sizeof(INT));
	  last_u=nzp; 
	} 
	U->JA[nzp] = j;
	nzp++;// = nzp + 1;
	if(j < minx) minx = j;
	U->IA[j] = i;// here j < n always
      }
    }
    qend = ichn[i];
    flag=0;
    if(qend >= 0) {
      k = qend;
      while(1){
	flag=0;
	k = ichn[k];
	kp1 = k + 1;
	it_strt = U->IA[k];
	it_end = U->IA[kp1]; 
	if(kp1==i) it_end = nzpi;
	U->IA[i] = i;
	for(jk = it_strt;jk<it_end;++jk){
	  j = U->JA[jk];
	  if(U->IA[j]!=i) {
	    if(nzp>last_u){
	      U->JA=realloc(U->JA,(nzp+1)*sizeof(INT));
	      last_u=nzp;
	    }
	    U->JA[nzp] = j;
	    nzp++;
	    U->IA[j] = i;
	    if(j < minx) minx = j;
	  }
	}
	if(nzp==nzpp) {flag=1;break;}
	if(k==qend) {flag=0;break;}
      }
    }
    if(flag || minx!=n) {
      k = ichn[minx];
      if(k>=0) { 
	ichn[i] = ichn[k];
	ichn[k] = i;
      }  else {
	ichn[minx] = i;
	ichn[i] = i;
      }
    }
    U->IA[i] = nzpi;
  }
  U->IA[n-1] = nzp;
  U->IA[n] = nzp;
  U->nnz = U->IA[n];  
  free(ichn);
  U->JA=realloc(U->JA,U->nnz*sizeof(INT));
  return;
}
/*=====================================================================*/
/**
 * \fn static void numeric_factor_symm(dCSRmat *AU,dCSRmat *U, 
 *                                     dvector *adiag,dvector *dinv)
 *
 * \brief Numerical factorization of a symmetric sparse matrix. JU
 *        needs to be ordered (column indices in every row to be
 *        ordered) for this to work.
 *
 * \param *AU pointer to the upper triangle of the dCSR matrix
 *
 * \param *U pointer to a dCSRmat: upper triangle of dvector
 *           containing the nonzero structure of the factor. U->val is
 *           created here. U->JA needs to have ordered indices (two
 *           transpositions before call to this.
 *
 * \param *dinv pointer to a dvector with the inverse of the diagonal
 *              elements of U^T*D*U
 *
 * \return void; *U is modified and the numerical values of the factor
 *               are calculated here.
 *
 * \note Reference: Chapter 7 in Sergio Pissanetzky. Sparse matrix
 *                  technology. Academic Press Inc. [Harcourt Brace
 *                  Jovanovich Publishers], London, 1984.
 *
 * \note fortan77 (19960101); C99 (20220801)
 * 
 * \author Ludmil 
 * \date 20220802
 *
 */
static void numeric_factor_symm(dCSRmat *AU,dCSRmat *U,		\
				dvector *adiag,dvector *dinv)
{
  /*
  */
  INT i,j,k,l,ip1,lp1,jk,qend;
  INT ia_strt,ia_end,it_strt,it_end,it1_strt,it1_end,line0;
  INT n=U->row;
  REAL um;
  /*Action:*/
  if(U->val){
    free(U->val);
    U->val=NULL;
  }
  U->val=(REAL *)calloc(U->nnz,sizeof(REAL));
  INT *ichn=NULL, *next=NULL;
  ichn=(INT *)calloc(n,sizeof(INT));
  next=(INT *)calloc(n,sizeof(INT));
  for(k=0;k<n;++k){
    ichn[k] = -1;
    next[k]=-1;
  }
  //
  for(i=0;i<n;++i){
    ip1 = i + 1;
    it_strt = U->IA[i];
    it_end = U->IA[ip1];
    if(it_end>it_strt) {
      for(jk = it_strt;jk<it_end;++jk){
	dinv->val[U->JA[jk]] = 0e0;
      }
      ia_strt = AU->IA[i];
      ia_end = AU->IA[ip1];
      if(ia_end > ia_strt) {
	for(jk=ia_strt;jk<ia_end;++jk){
	  dinv->val[AU->JA[jk]] = AU->val[jk];
	}
      }
    }
    dinv->val[i] = adiag->val[i];
    qend = ichn[i];
    if(qend >= 0) {
      line0 = ichn[qend];
      while(1){
	l = line0;
	lp1=l+1;
	line0 = ichn[l];
	it1_strt = next[l];
	it1_end = U->IA[lp1];
	um = U->val[it1_strt] * dinv->val[l];
	for(jk = it1_strt;jk<it1_end;++jk){
	  j = U->JA[jk];
	  dinv->val[j] = dinv->val[j] - U->val[jk] * um;
	}
	U->val[it1_strt] = um;
	next[l] = it1_strt + 1;
	if(it1_strt!=(it1_end-1)) {
	  jk = U->JA[it1_strt+1];
	  j = ichn[jk];
	  if(j>=0) {
	    ichn[l] = ichn[j];
	    ichn[j] = l;
	  } else{
	    ichn[jk] = l;
	    ichn[l] = l;
	  }
	}
	if(l==qend)
	  break;
      }
    }
    dinv->val[i] = 1e0/dinv->val[i];
    if(it_end>it_strt){
      for(jk = it_strt;jk<it_end;++jk){
	U->val[jk] = dinv->val[U->JA[jk]];
      }
      jk = U->JA[it_strt];
      j = ichn[jk];
      if(j>=0) {
	ichn[i] = ichn[j];
	ichn[j] = i;
      } else{
	ichn[jk] = i;
	ichn[i] = i;
      } 
    }
    next[i] = it_strt;
  }
  free(ichn);
  free(next);
  return;
}
/**************************************************************************/
/**
 * \fn void* run_hazmath_factorize(dCSRmat *A,INT print_level)
 *
 * \brief Performs factorization of A using HAZMATH (assumes A is
 *        symmetric and in CSR format)
 *
 * \note This routine does factorization only.
 *
 * \param A	       	        Matrix A to be factored
 * \param print_level       Flag to print messages to screen
 *
 * \return Numeric	        Stores U^T D U decomposition of A
 *
 * \author    Ludmil Zikatanov
 *
 */
void *run_hazmath_factorize(dCSRmat *A,INT print_level)
{
  const SHORT is_sym=1;
  // as of now, A is assumed to have symmetric pattern for the symbolic factorization and symmetric for the numeric factorization. 
  INT n,nnz;
  // size of the output structure. 
  size_t total=2*sizeof(dCSRmat) + 1*sizeof(dvector) + 3*sizeof(SHORT);
  void *Num=(void *)calloc(total/sizeof(char),sizeof(char));
  dCSRmat *U=NULL,*L=NULL;
  dvector *dinv=NULL;
  SHORT *extra=NULL;
  hazmath_get_numeric(Num, &U,&dinv,&extra,&L);
  extra[0]=is_sym;
  extra[1]=0;extra[2]=0;// should not be used.
  // end alloc for the Numeric structure. 
  /**/
  if(print_level>6){
    fprintf(stdout,"\nUsing HAZMATH factorize: A=U^T*D*U\n");
    fflush(stdout);
  }
  dCSRmat AU,AL;
  dvector adiag;
  n=A->row; nnz=A->nnz;
  adiag=dvec_create(n);
  dinv->row=n;
  dinv->val=calloc(n,sizeof(REAL));
  /***************************/
  if(is_sym){
    AL=dcsr_create(0,0,0);
  } else {
    AL=dcsr_create(n,n,nnz);
  }
  AU=dcsr_create(n,n,nnz);
  dcsr_cp(A,&AU);
  dcsr_trilu_diag(is_sym,&AU,&AL,&adiag);
  /***************************/
  U->row=AU.row;
  U->col=AU.col;
  symbolic_factor_symm(&AU,U); //symbolic factorization
  U->nnz=U->IA[U->row];
  // AL is used as working:  
  AL.row=U->row; AL.col=U->col;  AL.nnz=U->nnz;
  AL.IA=realloc(AL.IA,(AL.row+1)*sizeof(INT));
  AL.JA=realloc(AL.JA,(AL.nnz)*sizeof(INT));
  //SHOULD BE NULL:  AL.val=NULL;//realloc(AL.val,(AL.nnz)*sizeof(REAL));
  // order U:
  U->val=NULL;
  dcsr_transz(U,NULL,&AL);
  dcsr_transz(&AL,NULL,U);
  dcsr_free(&AL);
  numeric_factor_symm(&AU,U,&adiag,dinv);
  dcsr_free(&AU);
  dvec_free(&adiag);
  return (void *)Num; 
}
/********************************************************************/
/**
 * \fn INT run_hazmath_solve(dCSRmat *A,dvector *f,dvector *x, void* Numeric, INT print_level)
 *
 * \brief Performs Gaussian Elmination (Reduction and back
 *        substitution for solving the symmetric system with
 *        decomposed factorized matrix): Ax = f, A=U^T*dinv*U using
 *        HAZMATH (Assumes A is symmetric and in CSR format and that
 *        the factorization has been done) 
 *
 * \note This routine does a solve only
 *
 * \param A	       	      Matrix A of the system to be solved
 * \param f	       	      Right hand side vector
 * \param Numeric           contains the U^T*D*U decomposition of A
 * \param print_level       Flag to print messages to screen
 *
 * \return x	         	    Solution
 *
 * \param *U pointer to the dCSR matrix with the U factor in A=U^T*dinv*U
 *
 * \param *dinv pointer to a dvector containing the inverse of the
 *              diagonal of U.
 *
 * \param *x pointer to a dvector with the right hand side, on output
 *           it contains the solution.
 *
 * \return void; *x is modified
 *
 * \note fortan77 (19960301); C99 (20220801)
 * \author Ludmil Zikatanov
 *
 */
INT run_hazmath_solve(dCSRmat *A,dvector *f,dvector *x, \
		  void* Numeric, INT print_level)
{
  if(print_level>6){
    fprintf(stdout,"\nUsing HAZMATH solve for U^T*D*U*x=b\n");
    fflush(stdout);
  }
  // arrays
  dCSRmat *U,*L=NULL;  dvector *dinv;  SHORT *extra;
  // get them from *Numeric
  hazmath_get_numeric(Numeric, &U, &dinv,&extra, &L);
  //
  x->row=f->row;
  memcpy(x->val,f->val,x->row*sizeof(REAL));
  /*
    NOW action:
    --------------------------------------------------------------------
    ...  REDUCTION AND BACK SUBSTITUTION FOR SOLVING THE
    ...  SYMMETRIC SYSTEM WITH FACTORIZED MATRIX
    --------------------------------------------------------------------
  */
  INT k,i,it_strt,it_end;
  REAL xtmp;
  INT nm = U->row-1;
  for(k=0;k<nm;++k){
    it_strt = U->IA[k];
    it_end = U->IA[k+1];
    xtmp = x->val[k];
    if(it_end > it_strt){
      for(i= it_strt;i<it_end;++i){
	x->val[U->JA[i]] -= U->val[i]*xtmp;
      }
    }
    x->val[k] = xtmp*dinv->val[k];
  }
  x->val[nm] *= dinv->val[nm];
  for(k = nm;k>=0;--k){
    it_strt = U->IA[k];
    it_end = U->IA[k+1];
    if(it_end > it_strt){
      xtmp = x->val[k];
      for(i = it_strt;i<it_end;++i){
	xtmp -= U->val[i] *x->val[U->JA[i]];
      }
      x->val[k] = xtmp;
    }
  }
  return (INT )SUCCESS; 
}
/********************************************************************/
/**
 * \fn void hazmath_get_numeric(void *Numeric, dCSRmat **U, dvector **dinv, SHORT **extra,dCSRmat *L)
 *
 * \brief from the hazmath structure Numeric gets U and D
 *
 * \param Numeric   Pointer to the structure with numerical factorization
 *
 * \param **U      Pointer to the (*dCSRmat) which will be extracted from Numeric. 
 *
 * \param **dinv   pointer to *dinv which will be extracted from Numeric
 *
 * \author Ludmil Zikatanov
 * \date   20220802
 */
void hazmath_get_numeric(void *Numeric, dCSRmat **U, dvector **dinv, SHORT **extra,dCSRmat **L)
{
  void *wrk=(void *)Numeric;
  extra[0]=(SHORT *)wrk;
  *extra[0]=1; 
  wrk+=3*sizeof(SHORT);
  //  SHORT is_sym=extra[0];
  U[0]=(dCSRmat *)wrk;
  wrk+=sizeof(dCSRmat);
  dinv[0]=(dvector *)wrk;
  wrk+=sizeof(dvector);
  L[0]=(dCSRmat *)wrk;
  wrk+=sizeof(dCSRmat);
  //
  return; 
}
/*=====================================================================*/
/**
 * \fn static void symbolic_factor_symm(dCSRmat *AU, dCSRmat *U)
 *
 * \brief checks if a CSR  matrix is symmetric
 *
 * \param *A pointer to the dCSR matrix
 *
 * \return: 0 if the matrix is symmetric;
 *          1 pattern is symmetric but the matrix is not;
 *          2 pattern is not symmetric;
 *          3 matrix is not square
 *
 * \author Ludmil 
 * \date 20220802
 *
 */
INT chk_symmetry(dCSRmat *A)
{ 
  if(A->row != A->col){
    fprintf(stderr,"\n\n*** WARNING(%s): the matrix is NOT square: (%d x %d)\n",__FUNCTION__,A->row,A->col);
    return 3;
  }
  INT flag=0,i,j,jk,jt,irow;
  INT istrt,itstrt,iend,itend;
  REAL diff0=0e0,d=-1e20,tol=1e-14;
  dCSRmat Aord,AT;  
  AT=dcsr_create(A->col,A->row,A->nnz);
  dcsr_transz(A,NULL,&AT);
  Aord=dcsr_create(A->row,A->col,A->nnz);
  dcsr_transz(&AT,NULL,&Aord);
  for (i=0;i<Aord.row;++i){
    //
    istrt=Aord.IA[i];  iend=Aord.IA[i+1];
    //
    itstrt=AT.IA[i];    itend=AT.IA[i+1];
    //
    if((itend-itstrt) != (iend-istrt)){
      flag=2;
      irow=i;
      break;
    }
    flag=0;
    for(jk=istrt;jk<iend;++jk){
      j=Aord.JA[jk];
      jt=AT.JA[jk];
      if((j-jt)==0){
	d=fabs(Aord.val[jk]-AT.val[jk]);
	if(d>diff0)
	  diff0=d;
      } else {
	flag=4;// non-symmetric pattern
	break;
      }
    }
    if(flag !=0) {
      flag=2;// pattern is not symmetric
      irow=i;// record the row
      break;
    }
  }
  if(diff0>tol)
    flag=1;  
  switch(flag){
  case 1: 
    fprintf(stderr,"\n\n********************************************************");
    fprintf(stderr,"\n* WARNING(%s): the matrix seems nonsymmetric \n* %.3e=|A-A^T| > tol_symm=%.3e\n",__FUNCTION__,diff0,tol);
    fprintf(stderr,"\n********************************************************\n\n");
    break;
  case 2:
    fprintf(stderr,"\n\n********************************************************");
        fprintf(stderr,"\n* WARNING(%s): nonsymmetric pattern at row=%d\n",__FUNCTION__,irow);
    fprintf(stderr,"\n********************************************************\n\n");
  default:
    break;
  }
  dcsr_free(&Aord);
  dcsr_free(&AT);  
  return flag;
}
