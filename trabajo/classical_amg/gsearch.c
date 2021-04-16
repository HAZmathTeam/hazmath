/*! \file src/amr/gsearch.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note containing all essentials routines for mesh refinement
 *
 */
#include "hazmath.h"
/*********************************************************************/
// main
void main()
{
  dCSRmat *a=malloc(1*sizeof(dCSRmat));
  // read the matrix in coo format; convert it to csr and this is a:
  dcoo_read_dcsr("graphij0",a);
  fprintf(stdout,"\nDepth First Search permutation:\n");
  INT i,k;
  // find the connected components
  iCSRmat *blk=haz_dfs(a);
  fprintf(stdout,"\nend of list: blk->IA[%d]=%d (.eq. %d \?)\n",blk->row,blk->IA[blk->row],blk->row);
  for(i=0;i<blk->row;i++){
    fprintf(stdout,"\nblock[%d]; size=%d:  ",i,blk->IA[i+1]-blk->IA[i]);
    for(k=blk->IA[i];k<blk->IA[i+1];k++){ 
      fprintf(stdout,"%d ",blk->JA[k]);
    } 
  }  
  fprintf(stdout,"\nblocks=%d; vertices=%d; vertices(orig)=%d\n",blk->row,blk->col,a->col);
  return;
}
/*END********************************************************************/
