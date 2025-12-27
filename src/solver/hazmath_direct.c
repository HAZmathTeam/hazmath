#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "hazmath.h"
// /**************************************************************************/
// /**
//  * \fn void do_permutation(dCSRmat* a, INT* perm, INT* iperm)
//  *
//  * \brief  gives a permutation using "algorithm"
//  *
//  * \param perm:      ivector which on output contains the permutation.
//  *
//  * \param ia,ja:     the structure of a icsr(dcsr) matrix a
//  *
//  * \param algorithm: currently ONLY a min degree algorithm in two
//  *                   seteps: 1. construct (n by ?)  sparse matrix
//  *                   D(i,degree)=1, where the degree gives the degree
//  *                   of the i-th node. 2. The CSR transpose of this
//  *                   automatically gives the ordering by degree:
//  *                   D^T(degree,:) degree=1,2,\ldots max_degree is in
//  *                   DT->JA * \author Ludmil Zikatanov (20220810)
//  *
//  */
// /**********************************************************************************/
void do_permutation(dCSRmat* a, INT* perm, INT* iperm) 
{
  INT n = a->row;
  for (INT i = 0; i < n; ++i) perm[i] = i;
  if (0) {
    INT *ia = a->IA, *ja = a->JA;
    INT j, iblk, istrt, iend, nb, nblk;
    //  INT nnz=ia[n];
    iCSRmat* idfs = run_dfs(n, ia, ja);
    //  icsr_write_icoo("DFS",idfs);

    memcpy(perm, idfs->JA, n * sizeof(INT));
    ///
    dCSRmat a_t, awrk;
    ///
    awrk = dcsr_create(0, 0, 0);
    awrk.IA = calloc(n + 1, sizeof(INT));
    awrk.JA = calloc(n, sizeof(INT));
    awrk.val = NULL;
    //
    a_t = dcsr_create(0, 0, 0);
    a_t.IA = calloc(n + 1, sizeof(INT));
    a_t.JA = calloc(n, sizeof(INT));
    a_t.val = NULL;
    ///
    INT num, ideg, nnzp;
    for (iblk = 0; iblk < idfs->row; ++iblk) {
      istrt = idfs->IA[iblk];
      iend = idfs->IA[iblk + 1];
      nblk = iend - istrt;
      if (nblk < 2) continue;
      nb = 0;
      nnzp = 0;
      awrk.IA[0] = nnzp;
      awrk.row = nblk;
      awrk.col = 0;
      for (j = istrt; j < iend; ++j) {
        num = idfs->JA[j];
        perm[nb] = num;
        ideg = ia[num + 1] - ia[num] - 1;
        if (ideg < 0) ideg = 0;
        if ((ideg + 1) > awrk.col) awrk.col = ideg + 1;
        awrk.JA[nnzp] = ideg;
        nnzp++;
        awrk.IA[nb + 1] = nnzp;
        nb++;
      }
      if (nb != nblk) {
        fprintf(stderr, "\n\nERROR: %lld=nb != nblk=%lld in %s", (long long)nb,
                (long long)nblk, __FUNCTION__);
        exit(15);
      }
      awrk.nnz = awrk.IA[nblk];
      dcsr_transz(&awrk, NULL, &a_t);
      for (j = istrt; j < iend; ++j) {
        idfs->JA[j] = perm[a_t.JA[j - istrt]];
      }
    }
    memcpy(perm, idfs->JA, n * sizeof(INT));
    dcsr_free(&a_t);
    dcsr_free(&awrk);
    icsr_free(idfs);
    free(idfs);
  }
  if (iperm != NULL)
    for (INT i = 0; i < n; ++i) iperm[perm[i]] = i;
  return;
}
/********************************************************************* */
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
static void dcsr_trilu_diag(const SHORT is_sym,				\
			    dCSRmat *a,dCSRmat *al, dvector *adiag)
{
  INT k,j,kj,kj0,kj1;
  // transpose a: al should be allocated upon entry here.
  // permuta a;
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
	if( ((k-j)>=0) || (al->val[kj]==0e0) ) continue; // lower is rowind>=colind
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
      if( ((k-j)>0) || (a->val[kj]==0e0) ) continue; // lower is rowind>=colind
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
static void symbolic_factor_symm(dCSRmat *AU, dCSRmat *U)
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
void numeric_factor_symm(dCSRmat *AU,dCSRmat *U,		\
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
/*=====================================================================*/
/**
 * \fn static void numeric_factor_gen(dCSRmat *AU,dCSRmat *AL,
 *                                    dCSRmat *U, dCSRmat *L,
 *                                    dvector *adiag,dvector *dinv)
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
static void numeric_factor_gen(const SHORT is_sym,		\
			       dCSRmat *AU,dCSRmat *AL,		\
			       dCSRmat *U,dCSRmat *L,		\
			       dvector *adiag,dvector *dinv)
{
  /*
  */
  INT i,j,k,l,ip1,lp1,jk,qend;
  INT ia_strt,ia_end,it_strt,it_end,it1_strt,it1_end,line0;
  INT n=U->row;
  REAL lm,um, *di=dinv->val;
  /*Action:*/
  if(U->val){
    free(U->val);
    U->val=NULL;
  }
  if(L->val){
    free(L->val);
    L->val=NULL;
  }
  L->row=U->row; L->col=U->col;  L->nnz=U->nnz;
  L->IA=calloc((L->row+1),sizeof(INT));
  L->JA=calloc(L->nnz,sizeof(INT));
  // order U:
  dcsr_transz(U,NULL,L);
  dcsr_transz(L,NULL,U);
  //end order U:
  REAL *xl=NULL;
  U->val=(REAL *)calloc(U->nnz,sizeof(REAL));
  if(is_sym){
    dcsr_free(L);
    L=U;
    numeric_factor_symm(AU,U,adiag,dinv);
    return;
  } else {
    // get the lower triangle by columns:
    memcpy(L->IA,U->IA,(U->row+1)*sizeof(INT));
    memcpy(L->JA,U->JA,(U->nnz)*sizeof(INT));
    L->val=(REAL *)calloc(L->nnz,sizeof(REAL));
    xl=(REAL *)calloc(n,sizeof(REAL));
  }
  INT *ichn=NULL, *next=NULL;
  ichn=(INT *)calloc(n,sizeof(INT));
  next=(INT *)calloc(n,sizeof(INT));
  for(k=0;k<n;++k){
    ichn[k] = -1;
    next[k]=-1;
  }
  //
  for(i=0;i<n;++i){
    //    fprintf(stdout,"\n\ni=%d:",i+1);fflush(stdout);
    ip1 = i + 1;
    it_strt = U->IA[i];
    it_end = U->IA[ip1];
    if(it_end>it_strt) {
      for(jk = it_strt;jk<it_end;++jk){
	di[U->JA[jk]] = 0e0;
	xl[L->JA[jk]] = 0e0;
      }
      ia_strt = AU->IA[i];
      ia_end = AU->IA[ip1];
      if(ia_end > ia_strt) {
	//	fprintf(stdout,"\n");
	for(jk=ia_strt;jk<ia_end;++jk){
	  di[AU->JA[jk]] = AU->val[jk];
	  xl[AL->JA[jk]] = AL->val[jk];
	}
      }
    }
    di[i] = adiag->val[i];
    qend = ichn[i];
    if(qend >= 0) {
      line0 = ichn[qend];
      while(1){
	l = line0;
	lp1=l+1;
	line0 = ichn[l];
	it1_strt = next[l];
	it1_end = U->IA[lp1];
	lm = L->val[it1_strt] * di[l];
	um = U->val[it1_strt] * di[l];
	for(jk = it1_strt;jk<it1_end;++jk){
	  j = U->JA[jk];
	  di[j] = di[j] - U->val[jk] * lm;
	  xl[j] = xl[j] - L->val[jk] * um;
	}
	L->val[it1_strt] = lm;
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
    di[i] = 1e0/di[i];
    xl[i] = 0e0;
    if(it_end>it_strt){
      for(jk = it_strt;jk<it_end;++jk){       
	U->val[jk] = di[U->JA[jk]];
	L->val[jk] = xl[L->JA[jk]];
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
  free(xl);
  return;
}
/**************************************************************************/
/**
 * \fn void* run_hazmath_factorize(dCSRmat *A,ivector *p,INT print_level,
 *                                  void *more_params)
 *
 * \brief Performs factorization of A using HAZMATH (assumes A is
 *        symmetric and in CSR format)
 *
 * \note This routine does factorization only.
 *
 * \param A	       	        Matrix A to be factored
 * \param A	       	        Matrix A to be factored
 * \param is_sym                [0/1] if the matrix is 
 *                              non-symmetric/symmetric;
 *
 * \param p                     if not NULL then it is a permutation 
 *                              of the vertices
 *
 *                              non-symmetric/symmetric;
 * \param print_level       Flag to print messages to screen
 *
 * \return Numeric	        Stores U^T D U decomposition of A
 *
 * \author    Ludmil Zikatanov
 *
 */
void *run_hazmath_factorize(dCSRmat *A, INT print_level,	\
			     void *more_params)
{
  // as of now, A is assumed to have symmetric pattern for the symbolic factorization numeric factorization.

  SHORT is_sym,use_perm,permute_algorithm;
  SHORT *mpar=(SHORT *)more_params;
  if(more_params!=NULL){
    is_sym=mpar[0];
    use_perm=mpar[1];
    permute_algorithm=mpar[2];
  } else {
    is_sym=0;
    use_perm=1;
    permute_algorithm=0;
  }
  INT n,nnz;
  // size of the output structure. 
  size_t total=2*sizeof(dCSRmat) + 1*sizeof(dvector) + 3*sizeof(SHORT) + 1*sizeof(ivector);
  void *Num=(void *)calloc(total/sizeof(char),sizeof(char));
  memset(Num,0,((size_t )(total/sizeof(char)))*sizeof(char));
  dCSRmat *U=NULL,*L=NULL;
  dvector *dinv=NULL;
  SHORT *extra=NULL;
  ivector *perm=NULL;
  hazmath_get_numeric(Num,&U,&dinv,&extra,&L,&perm);
  extra[0]=is_sym;
  extra[1]=use_perm;// this should always be 1, i.e. always use permutation;
  extra[2]=permute_algorithm;// not used for now.

  if(print_level>10) fprintf(stdout,"\nUsing HAZMATH factorize (on the coarsest grid): ");

    /*  if(print_level>10){ */
    /* if(extra[0] && extra[1]) */
    /*   fprintf(stdout,"A(p,p)=U^T*D*U "); */
    /* else if((!extra[0]) && extra[1]) */
    /*   fprintf(stdout,"A(p,p)=L*D*U "); */
    /* else if(extra[0] && (!extra[1])) */
    /*   fprintf(stdout,"A=U^T*D*U "); */
    /* else */
    /*   fprintf(stdout,"A=L*D*U "); */
    /* } */
  //    fprintf(stdout,"\n                                               is_symmetric[0/1]=%d; use_permutation[0/1]=%d\n",extra[0],extra[1]);
  //    fflush(stdout);
  //
  /* this is all set, if called again will set this again*/
  if(use_perm){
    perm->row=A->col;
    perm->val=calloc(perm->row,sizeof(INT));
//OLD    do_permutation(perm,A->col,A->IA,A->JA,(const SHORT )extra[2]);
    ivector iperm; iperm.row=0;iperm.val=NULL;
    do_permutation(A,perm->val,iperm.val);
  } else {
    perm->row=0;
    free(perm->val);perm->val=NULL;
  }
  /*end permutation*/
  /* if(print_level>10){ */
  /*   fprintf(stdout,"\nUsing HAZMATH factorize (is_sym=%d): A=L*D*U\n",is_sym); */
  /*   fflush(stdout); */
  /* } */
  dCSRmat AU,AL;
  dvector adiag;
  n=A->row; nnz=A->nnz;
  adiag=dvec_create(n);
  dinv->row=n;
  dinv->val=calloc(n,sizeof(REAL));// this is stored in Numeric. 
  /***************************/
  AU=dcsr_create(n,n,nnz);
  AL=dcsr_create(n,n,nnz);
  if((perm!=NULL) && (perm->val!=NULL) && perm->row){
    dcsr_transz(A,perm->val,&AL);
    dcsr_transz(&AL,perm->val,&AU);
  } else {
    dcsr_cp(A,&AU);
    if(!is_sym)
      dcsr_cp(A,&AL);
    else {
      dcsr_free(&AL);
      AL=dcsr_create(0,0,0);
    }
  }
  dcsr_trilu_diag(is_sym,&AU,&AL,&adiag);
  /* if(print_level>9){  */
  /*   dcsr_write_dcoo("AU.dat",&AU);  */
  /*   dcsr_write_dcoo("AL.dat",&AL); */
  /*   dvec_write("adiag.dat",&adiag);  */
  /* }  */
  /***************************/
  U->row=AU.row;
  U->col=AU.col;
  symbolic_factor_symm(&AU,U); //symbolic factorization, assuming symmetric pattern. 
  U->nnz=U->IA[U->row];
  //
  numeric_factor_gen(is_sym,&AU,&AL,U,L,&adiag,dinv);
  dcsr_free(&AU);
  dcsr_free(&AL);
  dvec_free(&adiag);
  return (void *)Num; 
}
/********************************************************************/
/**
 * \fn INT run_hazmath_solve(dCSRmat *A,dvector *f,dvector *x,     
 *                           void* Numeric, INT print_level)
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
 * \param Numeric           contains the L*D*U decomposition of A
 * \param print_level       Flag to print messages to screen
 *
 * \param *f pointer to a dvector with the right hand side.
 *           
 * \param *x pointer to a dvector the solution.
 *
 * \note fortan77 (19960301); C99 (20220801)
 * \author Ludmil Zikatanov
 *
 */
INT run_hazmath_solve(dCSRmat *A,dvector *f,dvector *x, \
		      void* Numeric, INT print_level)
{
  // arrays
  dCSRmat *U,*L=NULL;  dvector *dinv;  SHORT *extra; ivector *perm;
  // get them from *Numeric
  hazmath_get_numeric(Numeric, &U, &dinv,&extra, &L, &perm);
  if(print_level>10){
    fprintf(stdout,"\nUsing HAZMATH solve: ");
    if(extra[0] && extra[1])
      fprintf(stdout,"U^T*D*U*x(p)=f(p) ");
    else if((!extra[0]) && extra[1])
      fprintf(stdout,"L*D*U*x(p)=f(p) ");
    else if(extra[0] && (!extra[1]))
      fprintf(stdout,"A=U^T*D*U*x=f ");
    else
      fprintf(stdout,"A=L*D*U*x=f ");
    fprintf(stdout,"\n                    (is_symmetric[0/1]=%lld; use_permutation[0/1]=%lld)\n",(long long )extra[0],(long long )extra[1]);
    fflush(stdout);
  }
  //
  INT k,i,it_strt,it_end;
  REAL xtmp;
  INT nm = U->row-1;
  x->row=f->row;
  memcpy(x->val,f->val,x->row*sizeof(REAL));
  /*
    NOW action:
    --------------------------------------------------------------------
    ...  REDUCTION AND BACK SUBSTITUTION FOR SOLVING
    ...  A LINEAR SYSTEM WITH FACTORIZED MATRIX
    --------------------------------------------------------------------
  */
  if(extra[0]) // if symmetric, then L just points to U;
      L=U;
  if(print_level>6){
    fprintf(stdout,"\nSolve phase (is_symmetric(0/1)=%lld; use_permutation(0/1)=%lld)\n", \
	    (long long )extra[0],(long long )extra[1]);
  }
  // if we have a permutation:
  if(extra[1] && (perm!=NULL) && (perm->val !=NULL)&& (perm->row)){
    INT k0;
    for(k0=0;k0<nm;++k0){
      k=perm->val[k0];
      it_strt = L->IA[k0];
      it_end = L->IA[k0+1];
      xtmp = x->val[k];
      if(it_end > it_strt){
	for(i= it_strt;i<it_end;++i){
	  x->val[perm->val[L->JA[i]]] -= L->val[i]*xtmp;
	}
      }
      x->val[k] = xtmp*dinv->val[k0];
    }
    x->val[perm->val[nm]] *= dinv->val[nm];
    for(k0 = nm;k0>=0;--k0){
      k=perm->val[k0];
      it_strt = U->IA[k0];
      it_end = U->IA[k0+1];
      if(it_end > it_strt){
	xtmp = x->val[k];
	for(i = it_strt;i<it_end;++i){
	  xtmp -= U->val[i] *x->val[perm->val[U->JA[i]]];
	}
	x->val[k] = xtmp;
      }
    }
  } else {
    // no permutation:
    for(k=0;k<nm;++k){
      it_strt = L->IA[k];
      it_end = L->IA[k+1];
      xtmp = x->val[k];
      if(it_end > it_strt){
	for(i= it_strt;i<it_end;++i){
	  x->val[L->JA[i]] -= L->val[i]*xtmp;
	}
      }
      x->val[k] = xtmp*dinv->val[k];
    }
    //  x->val[nm] *= dinv->val[nm];
    x->val[nm] *= dinv->val[nm];
    for(k=nm;k>=0;--k){
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
  }
  return (INT )SUCCESS; 
}
/********************************************************************/
/**
 * \fn void hazmath_get_numeric(void *Numeric, dCSRmat **U, dvector **dinv, SHORT **extra,dCSRmat *L,ivector *perm)
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
void hazmath_get_numeric(void *Numeric, dCSRmat **U, dvector **dinv, SHORT **extra,dCSRmat **L,ivector **perm)
{
  void *wrk=(void *)Numeric;
  extra[0]=(SHORT *)wrk;
  //  *extra[0]=1; 
  wrk+=3*sizeof(SHORT);
  //  SHORT is_sym=extra[0];
  U[0]=(dCSRmat *)wrk;
  wrk+=sizeof(dCSRmat);
  dinv[0]=(dvector *)wrk;
  wrk+=sizeof(dvector);
  L[0]=(dCSRmat *)wrk;
  wrk+=sizeof(dCSRmat);
  perm[0]=(ivector *) wrk;
  wrk+=sizeof(ivector);
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
    fprintf(stderr,"\n\n*** WARNING(%s): the matrix is NOT square: (%lld x %lld)\n",__FUNCTION__,(long long )A->row,(long long )A->col);
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
        fprintf(stderr,"\n* WARNING(%s): nonsymmetric pattern at row=%lld\n",__FUNCTION__,(long long )irow);
    fprintf(stderr,"\n********************************************************\n\n");
  default:
    break;
  }
  dcsr_free(&Aord);
  dcsr_free(&AT);  
  return flag;
}
