#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "hazmath.h"
/**
 * \fn ivector sparse_MIS0(dCSRmat *a, iCSRmat **dfs_bfs)
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
static ivector *sparse_MIS0(dCSRmat *a, iCSRmat **dfs_bfs, ivector *anc)
{
    //!  A 
    INT n = a->row;
    INT *ia = a->IA;
    INT *ja = a->JA;
    INT *iord = dfs_bfs[1]->JA;

    // local variables
    INT i,j,ii;
    INT row_begin, row_end;
    INT count=0;
    INT *flag;
    flag = (INT *)calloc(n, sizeof(INT));
    memset(flag, 0, sizeof(INT)*n);

    //! return
    ivector *mis=malloc(1*sizeof(ivector));
    mis->row = 0; // this is temporary, not used till the end
    mis->val = (INT*)calloc(n,sizeof(INT));

    //main loop
    if(iord != NULL){
      for (ii=0;ii<n;ii++) {
	i=iord[ii];
        if (flag[i] == 0) {
	  flag[i] = 1;
	  row_begin = ia[i]; row_end = ia[i+1];
	  for (j = row_begin; j<row_end; j++) {
	    if (flag[ja[j]] > 0) {
	      flag[i] = -1;
	      break;
	    }
	  }
	  if (flag[i]>0) {
	    mis->val[count] = i; count++;
	    for (j = row_begin; j<row_end; j++) {
	      flag[ja[j]] = -1;
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
	  for (j = row_begin; j<row_end; j++) {
	    if (flag[ja[j]] > 0) {
	      flag[i] = -1;
	      break;
	    }
	  }
	  if (flag[i]>0) {
	    mis->val[count] = i; count++;
	    for (j = row_begin; j<row_end; j++) {
	      flag[ja[j]] = -1;
	    }
	  }
        } // end if
      }// end for
    }
    // form Maximal Independent Set
    mis->row = count;
    mis->val=(INT *)realloc(mis->val, count*sizeof(INT));
    // clean
    if (flag) free(flag);
    //return
    return mis;
}
/*********************************************************************/
INT main(INT argc, char **argv)
{
  INT dummy,ish,np1,n,m,nnz,k,iread;
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
  if(argc>1){
    //    fp=fopen(argv[1],"r");
    dcoo_read_dcsr(argv[1],&a);
    n=a.row;
    nnz=a.nnz;
    ia=a.IA;
    ja=a.JA;
  } else {
    fp=fopen("abc","r");
    iread=fscanf(fp,"%d %d %d",&n,&m,&nnz);
    ia=calloc(n+1,sizeof(INT));
    for(k=0;k<=n;k++) iread=fscanf(fp,"%d",ia+k);
    // shift if ia[0] is not 0.
    ish=0;
    if(ia[0]>0) ish=-ia[0];  
    nnz = ia[n]+ish;
    ja=calloc(nnz,sizeof(INT));
    for(k=0;k<nnz;k++) iread=fscanf(fp,"%d",ja+k);
    fclose(fp);
  }
  fprintf(stdout,"\n%%%%n=%3d,nnz=%3d,shift=%3d\n",n,nnz,ish);fflush(stdout);
  if(ish){
    for(k=0;k<=n;k++) ia[k]+=ish;
    for(k=0;k<nnz;k++) ja[k]+=ish;
  }
  //
  anc=ivec_create(n);
  inv_ord=ivec_create(n);
  ldfsbfs=lex_bfs(n,ia,ja,&inv_ord,&anc);
  ldfs=ldfsbfs[0];  
  lbfs=ldfsbfs[1];  
  //
  fprintf(stdout,"\np=[");
  for(k=0;k<n;k++){
    fprintf(stdout,"%3d ",lbfs->JA[k]+1);
  }
  fprintf(stdout,"];\n");
  fprintf(stdout,"\npinv=[");
  for(k=0;k<n;k++){
    fprintf(stdout,"%3d ",inv_ord.val[k]+1);
  }
  fprintf(stdout,"];\n");
  fprintf(stdout,"\nlevel=[");
  for(k=0;k<n;k++){
    fprintf(stdout,"%3d ",lbfs->val[lbfs->JA[k]]+1);
  }
  fprintf(stdout,"];\n");
  fprintf(stdout,"\ntree=[");
  for(k=0;k<n;k++){
    fprintf(stdout,"%3d ",anc.val[lbfs->JA[k]]+1);
  }
  fprintf(stdout,"];\n");
  icsr_free(ldfsbfs[0]);
  icsr_free(ldfsbfs[1]);
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
