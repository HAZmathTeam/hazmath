/*! \file src/graphs/dfs.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note depth first search and related routines. 
 *
 */
#include "hazmath.h"
/*****************************************************************************/
static void dfs_recurrence(int v, int *ia, int *ja,		\
			   int *mask, int *jblk, int *pos)
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
/*****************************************************************************/
iCSRmat *haz_dfs(dCSRmat *a)
{
  // depth first search: find all conncted components
  //  blk->val is the connected component number for the vertex blk->JA[i]
  INT i,j,k,pos;
  // short hand;
  iCSRmat *blk=malloc(1*sizeof(iCSRmat));
  blk[0]=icsr_create(a->row,a->col,a->row);
  INT *mask=blk->val;
  for (i=0;i<a->row;++i) blk->val[i]=0;  
  blk->row=0;
  blk->IA[0]=0;
  pos=blk->IA[0];
  for (i=0;i<a->row;++i){
    if(!blk->val[i]){
      dfs_recurrence(i,a->IA,a->JA,blk->val,blk->JA,&pos);
      blk->IA[++blk->row]=pos;
    }
  }
  blk->IA=realloc(blk->IA,(blk->row+1)*sizeof(INT));
  return blk;
}
/*************************************************************************/
void dfs00_(INT *nin, INT *ia, INT *ja, INT *nblko,INT *iblk, INT *jblk)
{
  INT n,n1,nb,nblk,count,k1,sp,vp,v=0,wp=0,w=0,v1=0,i,ii,myflag;
  INT *numb=NULL,*iedge=NULL,*lowlink=NULL;
  /*---------------------------------------------------------------------
   This function is Tarjan's DFS algorithm for finding the strongly
   connected components of a diggraph.  Based on the FRED GUSTAVSON'S
   implementation of the algorithm.  The new numbering is in CSR
   format the pointers are in C...  IBLK, and the actual numbering is
   in JBLK. all diagonal elements must be nonzero in ia,ja.
    ---------------------------------------------------------------------*/
  /*C...  Initialization.*/
  n=*nin;
  n1 = n+1;
  iedge = (INT *) calloc(n1,sizeof(INT));
  numb = (INT *) calloc(n1,sizeof(INT));
  lowlink= (INT *) calloc(n1,sizeof(INT));
  /* Check */
  if(!(numb && iedge && lowlink)) {
    fprintf(stderr,"\nCould not allocate local mem in dfs\n");
    exit(19);
  }

  nblk = -1;
  nb = -1;
  count  = nb+1;
  k1 = count;
  vp = n;//
  sp = n;// may be n???
  for (i=0; i<n;i++) {
    iedge[i] = ia[i];
    numb[i]=-1;
    lowlink[i]=-1;
  }
  numb[n] = -1;
  lowlink[n] = -1;
  myflag=10;
  while (1) { //10 continue
    //    fprintf(stdout," myflag=%i %i, %i %i\n",myflag,nblk,count,n);fflush(stdout);
    if(myflag==10){
  /*
    C...  Get out when the  renumbering is done;
    C...  Otherwise begin at vertex K1 
  */
      if(count==n) {
	nblk++; //not sure why....
	iblk[nblk] = n;
	*nblko=nblk;
	if(iedge) free(iedge);
	if(numb) free(numb);
	if(lowlink) free(lowlink);
	return;
      }
      for (ii=k1;ii<n;ii++) {
	i=ii;
	if(numb[ii]<0) {
	  myflag=30;
	  break; //go to 30
	}
      }
      if(myflag !=30){
	fprintf (stderr,"There is an error in DEPTH FIRST SEARCH:  %i == %i\n",i,n);
	exit(254);
      } else {
	v = i;
	k1 = v + 1;
	myflag=50;
      }
    }
    /*C...  :::*/
    if(myflag==40) {
      vp--;
      iblk[vp] = v;
      v=w;
      myflag=50;
    }
    if(myflag==50) {
      nb++;
      numb[v] = nb;
      lowlink[v] = numb[v];
      sp--;
      jblk[sp] = v;
      v1 = v+1;
      myflag=60;
    }
    /*...  */
    while(60) { // 60   continue;
      if(myflag == 60) {
	wp = iedge[v];
	w = ja[wp];
	//	fprintf(stdout,"\nv=%d myflag=%i wp=%i, count=%i n=%i\n",v,myflag,wp,count,n);fflush(stdout);
	iedge[v] = wp+1;
	if(numb[w] >= numb[v])  {
	  myflag=70;
	} else if (numb[w]<0) {
	  myflag=40; 
	  break;
	} else {
	  if(lowlink[v]>numb[w])
	    lowlink[v]=numb[w];
	  myflag=70;
	}
      }
      if(iedge[v] < ia[v1]) {
	myflag=60;
	continue;
      }
      /*...*/  
      if(lowlink[v] >= numb[v]) {//don't {go to 90}
	nblk++;
	iblk[nblk] = count;
	/**/
	while(80) {     //    80 continue;
	  w = jblk[sp];
	  numb[w] = n;
	  sp++;
	  jblk[count] = w;
	  count++;
	  if(v == w) break; //{don't go to 80}
	}
	/*C...  */
	if(sp==n) {
	  myflag=10; 
	  break;
	}
      } 
      myflag=70;
      w = v;
      v = iblk[vp];
      vp++;
      v1 = v + 1;
      if(lowlink[v]>lowlink[w]) 
	lowlink[v]=lowlink[w];
    }
  }
}


