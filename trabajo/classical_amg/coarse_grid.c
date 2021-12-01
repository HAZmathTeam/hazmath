#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "hazmath.h"
/***********************************************************************************************/
/*!
 * \fn static SHORT icoo_2_icsr (dCOOmat *A, dCSRmat *B)
 *
 * \brief Transform an integer COOmat matrix to a iCSRmat format.
 *
 * \param A   Pointer to iCOOmat matrix
 * \param B   Pointer to iCSRmat matrix (must be allocated with same row col and nnz
 *
 * \return    SUCCESS if successed; otherwise, error information.
 *
 */
static SHORT icoo_2_icsr (iCOOmat *A,
			  iCSRmat *B)
{
    const INT m=A->row, n=A->col, nnz=A->nnz;
    B->row=A->row;
    B->col=A->col;
    B->nnz=A->nnz;
    B->IA=realloc(B->IA,(B->row+1)*sizeof(INT));
    B->JA=realloc(B->JA,B->nnz*sizeof(INT));
    if(A->val!=NULL){
      B->val=realloc(B->val,B->nnz*sizeof(INT));
    } else {
      free(B->val);
      B->val=NULL;
    }
    // local variables
    INT *ia = B->IA;
    INT *ja = B->JA;
    INT *row_idx = A->rowind;
    INT *col_idx = A->colind;
    INT i, iind, jind;

    INT *ind = (INT *) calloc(m+1,sizeof(INT));

    // initialize
    memset(ind, 0, sizeof(INT)*(m+1));

    // count number of nonzeros in each row
    for (i=0; i<nnz; ++i) ind[row_idx[i]+1]++;

    // set row pointer
    ia[0] = 0;
    for (i=1; i<=m; ++i) {
        ia[i] = ia[i-1]+ind[i];
        ind[i] = ia[i];
    }

    if(A->val!=NULL){
      // set column index and values.
      for (i=0; i<nnz; ++i) {
        iind = row_idx[i];
        jind = ind[iind];
        ja[jind] = col_idx[i];
        B->val[jind] = A->val[i];
        ind[iind] = ++jind;
      }
    }else{
      // set column index.
      for (i=0; i<nnz; ++i) {
        iind = row_idx[i];
        jind = ind[iind];
        ja[jind] = col_idx[i];
        ind[iind] = ++jind;
      }
    }
    if (ind) free(ind);
    return SUCCESS;
}
/***********************************************************************************************/
/*!
 * \fn static void level_ordering(iCSRmat *A, iCSRmat *lbfs, INT *p, INT *level)
 *
 * \brief store in A the level structure contained in LBFS.
 *
 * \param A   Pointer to iCSRmat matrix (num rows = num levels and num nnz is n)
 * \param lbfs   Pointer to iCSRmat matrix (must be allocated with same row col and nnz
 * \param p     bfs permutation
 * \param level for each vertex, its level. 
 *
 * \return    SUCCESS if successed; otherwise, error information.
 *
 */

static void level_ordering(iCSRmat *A, iCSRmat *lbfs)
{
  INT k,n=lbfs->nnz;
  iCOOmat tmp_mat;  
  tmp_mat.col=n;
  tmp_mat.row=-1;
  for(k=0;k<n;k++){
    if(tmp_mat.row<lbfs->val[lbfs->JA[k]]){
      tmp_mat.row = lbfs->val[lbfs->JA[k]];
    }
  }
  // uncomment this if we want to go level by level;
  tmp_mat.row++;
  tmp_mat.nnz=n;
  tmp_mat.rowind=lbfs->val;
  tmp_mat.colind=lbfs->JA;
  tmp_mat.val=NULL;
  icsr_free(A);
  *A=icsr_create(0,0,0);
  icoo_2_icsr(&tmp_mat,A);
  return;
}
/**
 * \fn static ivector sparse_MIS0(dCSRmat *a, iCSRmat **dfs_bfs)
 *
 * \brief get the maximal independet set of a CSR matrix
 *
 * \param a pointer to the matrix
 * \param bfs_blk is an iCSR matrix with bfs ordering: 
 *        bfs_blk->IA: should give the pointers to the beginning of
 *                     every connected component.
 *        bfs_blk->JA: is the ordering: from root(s) if we have
 *                     several connected components.
 *        bfs_blk->val: contains the level number of every point. 
 *
 * \note: 20210622 (ltz);
 */
static ivector sparse_MIS0(dCSRmat *a, iCSRmat *lbfs, ivector *anc)
{
    //!  A 
    INT n = a->row;
    INT *ia = a->IA;
    INT *ja = a->JA;
    INT *iord = lbfs->JA;
    // local variables
    INT i,j,ii,jk;
    INT row_begin, row_end;
    INT count=0;
    INT *flag;
    flag = (INT *)calloc(n, sizeof(INT));
    memset(flag, 0, sizeof(INT)*n);

    //! return
    ivector mis;
    mis.row = 0; // this is temporary, not used till the end
    mis.val = (INT*)calloc(n,sizeof(INT));
    if(iord != NULL){
      /*main loop*/
      for (ii=0;ii<n;ii++) {
	i=iord[ii];
        if (!flag[i]) {
	  flag[i] = 1;
	  row_begin = ia[i]; row_end = ia[i+1];
	  if(row_end-row_begin<2) {
	    /* as before: isolated point, just skip it; */
	    flag[i] = -1;
	    fprintf(stdout,"\nflag[%d]=%d, row_size=%d",i,flag[i],row_end-row_begin);
	    continue;
	  }
	  for (jk = row_begin; jk<row_end; jk++) {
	    j=ja[jk];
	    if ((flag[j] > 0) && (j != i)) {
	      flag[i] = -1;
	      break;
	    }
	  }
	  if (flag[i]>0) {
	    mis.val[count] = i; count++;
	    for (jk = row_begin; jk<row_end; jk++) {
	      flag[ja[jk]] = -1;
	    }
	  }
        } // end if
      }// end for
    } else {
      // no ordering
      for (i=0;i<n;i++) {
        if (flag[i] == 0) {
	  flag[i] = 1;
	  row_begin = ia[i]; row_end = ia[i+1];
	  if(row_end-row_begin<2) {
	    /* as before: isolated point, just skip it; */
	    flag[i] = -1;
	    continue;
	  }
	  for (jk = row_begin; jk<row_end; jk++) {
	    j=ja[jk];
	    if ((flag[j] > 0) && (j != i)) {
	      flag[i] = -1;
	      break;
	    }
	  }
	  if (flag[i]>0) {
	    mis.val[count] = i; count++;
	    for (j = row_begin; j<row_end; j++) {
	      flag[ja[j]] = -1;
	    }
	  }
        } // end if
      }// end for
    }
    // form Maximal Independent Set
    fprintf(stdout,"\n\ncount=%d; coarsening ratio = %7.2f\n",count, ((REAL )n)/((REAL )count));
    mis.row = count;
    mis.val=(INT *)realloc(mis.val, count*sizeof(INT));
    //    ivector_print(stdout,&mis);
    // clean
    if (flag) free(flag);
    //return
    return mis;
}
/*********************************************************************/
INT main(INT argc, char **argv)
{
  INT dummy,ish,np1,n,m,nnz,k;
  FILE *fp   = NULL;
  INT  *ia   = NULL, *ja    = NULL;
  dCSRmat a;
  ivector anc,inv_ord;
  iCSRmat **ldfsbfs=NULL,*lbfs=NULL, *ldfs=NULL;
  //
  //  fprintf(stdout,"\n%%%%argc=%3d\n",argc);
  //  for(k=0;k<argc;k++)
  //    fprintf(stdout,"\n%%%%arg(%3d)=\'%s\'",k,argv[k]);
  //  fprintf(stdout,"\n%%%%==============\n");
  ish=0;
  //    fp=fopen(argv[1],"r");
  dcoo_read_dcsr(argv[1],&a);
  n=a.row;
  nnz=a.nnz;
  ia=a.IA;
  ja=a.JA;
  fprintf(stdout,"\nn=%d,m=%d,nnz=%d",n,a.col,nnz);fflush(stdout);
  fprintf(stdout,"\n%%%%n=%3d,nnz=%3d,shift=%3d\n",n,nnz,ish);fflush(stdout);
  if(ish){
    for(k=0;k<=n;k++) ia[k]+=ish;
    for(k=0;k<nnz;k++) ja[k]+=ish;
  }
  //
  anc=ivec_create(n);
  inv_ord=ivec_create(n);
  ldfsbfs=lex_bfs(n,ia,ja,&inv_ord,&anc);
  lbfs=ldfsbfs[1];
  // dfs is not needed, use it as work space if we would like to choose CG level by level. 
  //  level_ordering(ldfsbfs[0],ldfsbfs[1]);
  //  icsr_print_rows(stdout,0,ldfsbfs[0],"levels");
  // Now choose MIS
  ivector mis=sparse_MIS0(&a, ldfsbfs[1],&anc);
  icsr_free(ldfsbfs[0]);
  icsr_free(ldfsbfs[1]);
  free(ldfsbfs[0]);
  free(ldfsbfs[1]);
  free(ldfsbfs);
  ivec_free(&anc);
  ivec_free(&inv_ord);
  if(argc>1){
    dcsr_free(&a);
    ivec_free(&inv_ord);
    ivec_free(&anc);
  }else{
    free(ia);
    free(ja);
  }
    return 0;
}
  //
  //
  /* fprintf(stdout,"\np=["); */
  /* for(k=0;k<n;k++){ */
  /*   fprintf(stdout,"%3d ",lbfs->JA[k]+1); */
  /* } */
  /* fprintf(stdout,"];\n"); */
  /* fprintf(stdout,"\npinv=["); */
  /* for(k=0;k<n;k++){ */
  /*   fprintf(stdout,"%3d ",inv_ord.val[k]+1); */
  /* } */
  /* fprintf(stdout,"];\n"); */
  /* fprintf(stdout,"\nlevel=["); */
  /* for(k=0;k<n;k++){ */
  /*   fprintf(stdout,"%3d ",lbfs->val[lbfs->JA[k]]+1); */
  /* } */
  /* fprintf(stdout,"];\n"); */
  /* fprintf(stdout,"\ntree=["); */
  /* for(k=0;k<n;k++){ */
  /*   fprintf(stdout,"%3d ",anc.val[lbfs->JA[k]]+1); */
  /* } */
  /* fprintf(stdout,"];\n"); */
