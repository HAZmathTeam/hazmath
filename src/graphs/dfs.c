/*! \file src/graphs/dfs.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note depth first search and related routines. 
 *
 */
#include "hazmath.h"
/***********************************************************************************************/
/*!
 * \fn static void dfs_recurrence(int v, int *ia, int *ja, int *mask,int *jdfs, int *pos)
 *
 * \brief recurrsive call of the dfs search: number a vertex if not
 *         numbered and move to one of its nonnumbered neighbors
 *
 * \param v                current vertex
 * \param (ia,ja):         adjacency structure of the graph;
 * \param mask:            array indicating if a vertex is numbered
 * \param jdfs:            array containing the permutation;
 * \param pos:             the position of the v in the permutation array;
 *
 *
 * \author Ludmil Zikatanov
 * \date   20210516
 */
static void dfs_recurrence(int v, int *ia, int *ja,		\
			   int *mask, int *jdfs, int *pos)
{
  int i, vi;
  //  fprintf(stdout,"%d ",v);fflush(stdout);
  mask[v] = 1;
  jdfs[pos[0]++]=v;
  for (vi = ia[v]; vi<ia[v+1]; ++vi){
    i=ja[vi];
    if (!mask[i]){
      dfs_recurrence(i,ia,ja,mask,jdfs,pos);
    }
  }
  return;
}
/***********************************************************************************************/
/*!
 * \fn iCSRmat *run_dfs(INT n, INT *ia, INT *ja)
 *
 * \brief recurrsive dfs search of a graph
 *
 * \param n                number of vertices
 * \param (ia,ja):         adjacency structure of the graph;
 *
 * \return pointer to iCSRmat of size n,n with n nonzeroes, which
 *         contains the connected components and the permutation of
 *         the vertices (component after component). The values in the
 *         matrix are connected component number for the vertex
 *         dfs->val[i]=connected component number(dfs->JA[i])
 *
 * \author Ludmil Zikatanov
 * \date   20210516
 */
iCSRmat *run_dfs(INT n, INT *ia, INT *ja)
{
  INT i,j,k,pos;
  // short hand;
  iCSRmat *dfs=malloc(1*sizeof(iCSRmat));
  dfs[0]=icsr_create(n,n,n);
  INT *mask=dfs->val;
  for (i=0;i<n;++i) dfs->val[i]=0;  
  dfs->row=0;
  dfs->IA[0]=0;
  pos=dfs->IA[0];
  for (i=0;i<n;++i){
    if(!dfs->val[i]){
      dfs_recurrence(i,ia,ja,dfs->val,dfs->JA,&pos);
      dfs->IA[++dfs->row]=pos;
    }
  }
  dfs->IA=realloc(dfs->IA,(dfs->row+1)*sizeof(INT));
  return dfs;
}
/*****************************************************************************/
/*!
 * \fn iCSRmat *dfs_di(void *a, const char c)
 *
 * \brief dfs on graphs given by INT or REAL CSR matrix.
 *
 * \param a:                  The CSR matrix 
 * \param c:                  a char ('r' or 'i' or 'R' or 'I')
 *                            indicating whether this is a REAL or INT matrix;
 *
 * \return returns the output of run_dfs(a->row,a->IA,a->JA)
 *
 * \author Ludmil Zikatanov
 * \date   20210516
 */
iCSRmat *dfs_di(void *a, const char c)
{
  /* 
     do a depth first search for INT or REAL matrix. It does not use
     a->val, so upon entry this can be null 
  */
  dCSRmat *ad=NULL;
  iCSRmat *ai=NULL;
  if(c=='R' || c=='r'){
    ad=(dCSRmat *)a;
    return run_dfs(ad->row, ad->IA, ad->JA);  
  } else if(c=='I' || c=='i'){
    ai=(iCSRmat *)a;
    return run_dfs(ai->row, ai->IA, ai->JA);  
  } else {
    fprintf(stderr,"### ERROR: Wrong value of c in %s: c=%c\n",__FUNCTION__,c);
    exit(ERROR_INPUT_PAR);
  }
}
/***********************************************************************************************/
/*!
 * \fn dfs00_(INT *nin, INT *ia, INT *ja, INT *nblko,INT *iblk, INT *jblk)
 *
 * \brief non-recurrsive dfs search of a graph.  This function
 *        implements Tarjan's DFS algorithm for finding the strongly
 *        connected components of a di-graph.  It is the FRED
 *        GUSTAVSON'S non-recurrsive version of the algorithm.  
 *
 * \param n                number of vertices
 * \param (ia,ja):         adjacency structure of the graph; all diagonal
 *                         elements must be nonzero in ia,ja.
 * \param nblko:           number of connected components
 * \param (iblk,jblk):     CSR matrix containing the permutation: vertices
 *                         are ordered ccomponent by compoent.
 *
 * \return The new numbering is in CSR format: the pointers are in
 *         IBLK[], and the actual numbering is in JBLK[].
 *
 * \author Ludmil Zikatanov
 * \date   20210516
 */
void dfs00_(INT *nin, INT *ia, INT *ja, INT *nblko,INT *iblk, INT *jblk)
{
  INT n,n1,nb,nblk,count,k1,sp,vp,v=0,wp=0,w=0,v1=0,i,ii,myflag;
  INT *numb=NULL,*iedge=NULL,*lowlink=NULL;
  /*---------------------------------------------------------------------
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


