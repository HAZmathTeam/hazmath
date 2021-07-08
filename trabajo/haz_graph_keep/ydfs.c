/*! \file src/amr/dfs1.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note containing all essentials routines for mesh refinement
 *
 */
#include "hazmath.h"
void dfs_recurrence(int v, int *ia, int *ja, int *mask, int *jblk, int *pos)
{
  int i, vi;
  //  fprintf(stdout,"%d ",v);fflush(stdout);
  mask[v] = 1;
  jblk[pos[0]++]=v;
  for (vi = ia[v]; vi<ia[v+1]; ++vi){
    i=ja[vi];
    if (!mask[i]){
      dfs_recurrence(i,ia,ja,mask,jblk,pos);
    }
  }
  return;
}
iCSRmat *haz_dfs(dCSRmat *a)
{
  // depth first search: find all conncted components
  //  blk->val is the connected component number for the vertex blk->JA[i]
  INT i,j,k,pos;
  // short hand;
  INT *ia=a->IA,*ja=a->JA,*blk->IA=blk->IA,*jblk=blk->JA;
  INT *mask=blk->val;
  for (i=0;i<nr;++i) blk->val[i]=0;  
  fprintf(stdout,"\nconnected\?\?\?\n");
  blk->nrow=0;
  blk->IA[0]=0;
  pos=blk->IA[0];
  for (i=0;i<a->nrow;++i){
    if(!blk->val[i]){
      dfs_recurrence(i,a->IA,a->JA,blk->val,blk->JA,&pos);
      /* fprintf(stdout,"i=%d; pos=%d,blk->nrow=%d; blk->IA[%d]=%d\n",		\ */
      /* 	      i,pos,blk->nrow[0]+1,blk->nrow[0],blk->IA[blk->nrow[0]]);fflush(stdout); */
      blk->IA[++blk->nrow]=pos;
    }
  }
  blk->IA=realloc(blk->IA,(blk->nrow+1)*sizeof(INT));
  return blk;
}

