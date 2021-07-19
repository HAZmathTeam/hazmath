/*! \file src/amr/gsearch.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note containing all essentials routines for mesh refinement
 *
 */
#include "hazmath.h"
/*************************************************************************/
// main
void main()
{
  dCSRmat *a=malloc(1*sizeof(dCSRmat));
  // read the matrix in coo format; convert it to csr and this is a:
  dcoo_read_dcsr("graphij0",a);
  fprintf(stdout,"\nDepth First Search permutation:\n");
  INT i,k;
  // find the connected components
  iCSRmat *blk_dfs=run_dfs(a);
  fprintf(stdout,"\nend of list: blk_dfs->IA[%d]=%d (.eq. %d \?)\n",blk_dfs->row,blk_dfs->IA[blk_dfs->row],blk_dfs->row);
  /* for(i=0;i<blk_dfs->row;i++){ */
  /*   fprintf(stdout,"\nblock[%d]; size=%d:  ",i,blk_dfs->IA[i+1]-blk_dfs->IA[i]); */
  /*   for(k=blk_dfs->IA[i];k<blk_dfs->IA[i+1];k++){  */
  /*     fprintf(stdout,"%d ",blk_dfs->JA[k]); */
  /*   }  */
  /* }   */
  fprintf(stdout,"\nblocks=%d; vertices=%d; vertices(orig)=%d\n",blk_dfs->row,blk_dfs->col,a->col);

  /*******************************************************************/
  ivector *roots=malloc(1*sizeof(ivector));
  // grab a point in every connected component:
  roots->row=blk_dfs->row;
  roots->val=(INT *)calloc(roots->row,sizeof(INT));
  for(i=0;i<blk_dfs->row;i++){
    k=blk_dfs->IA[i];
    roots->val[i]=blk_dfs->JA[k];
  }
  ivector *anc=malloc(1*sizeof(ivector));
  INT level_max=a->row+1;
  fprintf(stdout,"BFS:");
  iCSRmat *blk_bfs=run_bfs(a,roots,anc,level_max);
  fprintf(stdout,"\nAnother ord::");
  for(i=0;i<blk_bfs->row;i++){
    fprintf(stdout,"\nrow=%d: ",i);
    for(k=blk_bfs->IA[i];k<blk_bfs->IA[i+1];++k){
      fprintf(stdout,"%7d(%5d)",blk_bfs->JA[k],anc->val[blk_bfs->JA[k]]);
      //blk_bfs->val[blk_bfs->JA[k]]);
    }
  }
  fprintf(stdout,"\nBFS: levels=%d;vertices=%d; vertices(orig)=%d\n",blk_bfs->row,blk_bfs->IA[blk_bfs->row],a->col);
  fprintf(stdout,"\n");
  ivec_free(anc);
  ivec_free(roots);
  icsr_free(blk_dfs);
  icsr_free(blk_bfs);
  dcsr_free(a);
  /*******************************************************************/
  return;
}
/*END********************************************************************/
