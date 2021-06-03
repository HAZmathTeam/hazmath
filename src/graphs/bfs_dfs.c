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
 * \fn iCSRmat *run_bfs(INT n,INT *ia, INT *ja,ivector *roots,ivector *anc,const INT lvlmax)
 *
 * \brief recurrsive bfs search of a graph
 *
 * \param n                number of vertices
 * \param (ia,ja):         adjacency structure of the graph;
 * \param roots[]:         a vector containing starting points for BFS. Aiming
 *                         at having one starting point in every
 *                         connected component.
 * \param anc[]:           a vector containing the ancestors of the vertices:
 *                         v<->anc[v] is an edge in the bfs tree. 
 *
 * \return pointer to iCSRmat of size n,n with n nonzeroes, which
 *         contains the bfs ordering permutation of the vertices
 *         (component after component). The values in the matrix are
 *         the distance from the root vertex in every connected
 *         component: for the vertex v bfs->val[v]=distance to the
 *         root[j], j=1:ncomponents
 *
 * \author Ludmil Zikatanov
 * \date   20210516
 */
iCSRmat *run_bfs(INT n,INT *ia, INT *ja,	\
		 ivector *roots,		\
		 ivector *anc,			\
		 const INT lvlmax)
{
  /* n, ia, ja are the CSR matrix elements. */
  /* anc[v] is the ancestor of v; */
  /* roots[] is input; */
  INT i,k,q,v,vi,qbeg,qend,lvl;  
  iCSRmat *blk=malloc(1*sizeof(iCSRmat));
  blk[0]=icsr_create(n,n,n);
  anc->row=n;
  anc->val=(INT *)calloc(anc->row,sizeof(INT));
  for(i=0;i<blk->nnz;++i) blk->val[i]=0;
  for(i=0;i<anc->row;++i) anc->val[i]=-1;
  lvl=0;
  blk->IA[lvl]=0;
  k=blk->IA[lvl];
  if(roots->row<=0){
    /* take the first vertex as root if none are given as input */
    roots->row=1;
    roots->val=(INT *)realloc(roots->val,roots->row*sizeof(INT));
    roots->val[0]=0;
  }
  /* Now roots are set as they are either input or root[0]=0 */
  for(i=0;i<roots->row;++i){
    /* fprintf(stdout,"\nroots[%d]=%d",i,roots->val[i]); */
    blk->val[roots->val[i]]=lvl+1;
    blk->JA[k]=roots->val[i];
    k++;
  }
  blk->IA[lvl+1]=k;
  /* we need to repeat this */
  while(1){
    qbeg=blk->IA[lvl];
    qend=blk->IA[lvl+1];
    for(q=qbeg;q<qend;++q){
      v = blk->JA[q];
      for(vi=ia[v];vi<ia[v+1];++vi){
	i=ja[vi];
	if (!blk->val[i]){
	  blk->val[i]=lvl+1;// mark as visited;
	  blk->JA[k]=i;
	  k++;
	  anc->val[i]=v; // ancestor;
	}
	fprintf(stdout,"\nlvl=%d,v=%d; nbr=%d,blk->val=%d",lvl,v,i,blk->val[i]);fflush(stdout);
      }
    }
    lvl++;
    blk->IA[lvl+1]=k;    
    if(k<=qend) break;
  }
  /* fprintf(stdout,"\nord (k=%d):",k); */
  /* for(i=0;i<k;i++){ */
  /*   v=blk->JA[i]; */
  /*   fprintf(stdout,"\nblk->val[%d]=%d",v,blk->val[v]);fflush(stdout); */
  /* } */
  /* for(i=0;i<(lvl+1);i++){ */
  /*   fprintf(stdout,"\nblk->IA[%d]=%d",i,blk->IA[i]); */
  /* } */
  blk->row=lvl;
  blk->IA=(INT *)realloc(blk->IA,(blk->row+1)*sizeof(INT));
  //end
  return blk;
}
/*****************************************************************************/
/*****************************************************************************/
/*!
 * \fn iCSRmat *bfs_di(void *a, const char c,ivector *roots,
 *                     ivector *anc,const INT lvlmax)
 *
 * \brief dfs on graphs given by INT or REAL CSR matrix.
 *
 * \param a:                  The CSR matrix 
 * \param c:                  a char ('r' or 'i' or 'R' or 'I')
 *                            indicating whether this is a REAL or INT matrix;
 *
 * \return returns the output of bfs_di(a->row,a->IA,a->JA)
 *
 * \author Ludmil Zikatanov
 * \date   20210516
 */
iCSRmat *bfs_di(void *a, const char c,		\
		ivector *roots,			\
		ivector *anc,			\
		const INT lvlmax)
{
  /* 
     do a breadth first search for INT or REAL matrix. It does not use
     a->val, so upon entry these can be null 
  */
  dCSRmat *ad=NULL;
  iCSRmat *ai=NULL;
  if(c=='R' || c=='r'){
    ad=(dCSRmat *)a;
    return run_bfs(ad->row, ad->IA, ad->JA, roots, anc, lvlmax);  
  } else if(c=='I' || c=='i'){
    ai=(iCSRmat *)a;
    return run_bfs(ai->row, ai->IA, ai->JA, roots, anc, lvlmax);  
  } else {
    fprintf(stderr,"### ERROR: Wrong value of c in %s (c=%c)\n",__FUNCTION__,c);
    exit(ERROR_INPUT_PAR);
  }
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx*/
INT check0(weights *elem1, weights *elem2)
{
  if ( elem1->val < elem2->val)
    return -1;
   else if (elem1->val > elem2->val)
      return 1;
   else
      return 0;
}

INT check1(iweights *elem1, iweights *elem2)
{
  //descending order. 
  if  (elem1->mask > elem2->mask)
    return -1;
  else if(( elem1->mask == elem2->mask) && (elem1->val > elem2->val))
    return -1;
  else if(elem1->mask < elem2->mask)
    return 1;
  else if(( elem1->mask == elem2->mask) && (elem1->val < elem2->val))
    return 1;
  else
    return 0;
}

void getp(INT *ie, INT *je, REAL *w, INT ne, INT *p)
{
  INT k;
  weights *tosort=NULL;  
  tosort = (weights *)calloc(ne,sizeof(weights));
  for (k=0; k<ne;k++)
    {
      tosort[k].val=w[k];
      tosort[k].id=k;
    }
  //Sort now.
  qsort((void *) tosort,ne,sizeof(weights),(testit )check0 );                  
  for (k=0; k<ne;k++)
    p[k]=tosort[k].id;

  if(tosort) free(tosort);  
  return;
}
void getpz(REAL *z, INT nv, INT *p)
{
  INT k;
  weights *tosort=NULL;  
  tosort = (weights *)calloc(nv,sizeof(weights));
  for (k=0; k<nv;k++)
    {
      tosort[k].val=-z[k];
      tosort[k].id=k;
    }
  qsort((void *) tosort,nv,sizeof(weights),(testit )check0 );                  
  // after things are sorted,  we get the permutation
  for (k=0; k<nv;k++)
    p[k]=tosort[k].id;
  if(tosort) free(tosort);  
  return;
}
void getpi(INT *iz, INT *maskv, INT nv, INT *p)
{
  INT k;
  iweights *tosort=NULL;  
  tosort = (iweights *)calloc(nv,sizeof(iweights));
  for (k=0; k<nv;k++)
    {
      tosort[k].mask=maskv[k];
      tosort[k].val=iz[k];
      tosort[k].id=k;
    }
  qsort((void *) tosort,nv,sizeof(iweights),(testit )check1 );                  
  // after things are sorted,  we have the permutation
  for (k=0; k<nv;k++)
    p[k]=tosort[k].id;
  if(tosort) free(tosort);  
  return;
}
/*******************************************************************/
void bfsx(INT nv, INT *ia, INT *ja, INT *ibfs, INT *jbfs,	\
	INT *maske, INT *p, INT *et, INT *lev, REAL *w,REAL *z)
{
  INT i,j,k;
  INT i1,i2,k0,mj,kbeg,kend,ipoint,iai,iai1,klev;
  // intialize  
  memset(maske,0,nv*sizeof(INT));
  klev=1; //level number ; for indexing this should be klev-1;
  kbeg=ibfs[klev-1];
  kend=ibfs[klev];
  for(i1=kbeg;i1<kend;++i1){
    i=jbfs[i1];
    maske[i]=klev;//;
    et[i]=-1;
    p[i1-kbeg]=i1-kbeg;//trivial permutation.
  }
  //  fprintf(stdout,"level=%i; number of vertices: %i\n",klev,kend-kbeg);
  ipoint=ibfs[1];
  while(1) {
    k0=0;
    for(i2=kbeg;i2<kend;++i2){
      i1=p[i2-kbeg]+kbeg;
      i=jbfs[i1];
      iai = ia[i];
      iai1=ia[i+1];     
      for(k=iai;k<iai1;++k){
	j=ja[k];
	if(i==j) continue;
	mj=maske[j];
	if(!mj){
	  jbfs[ipoint]=j;
	  maske[j]=klev;
	  // bfs tree edge
	  et[j]=i; //father or mother of j is i.
	  w[k0]=z[j];
	  k0++;
	  ipoint++;
	}
      }	   
    }
    kbeg=kend;
    klev++;
    ibfs[klev]=ipoint;
    kend=ipoint;
    //    fprintf(stdout,"level=%i; number of vertices: %i\n",klev,kend-kbeg);
    //    fprintf(stdout,"level=%i; number of vertices: %i : %i\n",klev,kend-kbeg,k0);
    getpz(w, k0, p);
    if(ipoint >=nv)  break;
  }
  //  fprintf(stdout,"level=%i; TOTAL number of vertices: %i and %i\n",klev,ibfs[klev],nv);
  *lev=klev;
return;
}

void bfstree(INT root, INT nv, INT ne, INT *ia, INT *ja, INT *ibfs, INT *jbfs, \
	    INT *mask, INT *et, INT *lev, INT *ledge, REAL *w,REAL *wo) 
{
  /* bfs tree rooted at root */
  INT i,j,k,maxdeg,ih;
  INT i1,kbeg,kend,ipoint,iai,iai1,klev;
  if(root < 0) {
    maxdeg=-1;
    for(i=0;i<nv;++i){
      ih=ia[i+1]-ia[i];
      if(maxdeg<ih){
	maxdeg=ih;
	root=i;
      }
    }
    if(root<0){
      fprintf(stderr,"ERROR: could not find root for bfs\n");
      exit(129);
    }
  }
  // initialization
  /*************************************************************/
  memset(mask,0,nv*sizeof(INT));
  for(i=0;i<nv;++i) et[i]=-1;
  klev=1; //level number ; for indexing this should be klev-1;
  ibfs[0]=0;
  ibfs[1]=1;
  jbfs[ibfs[0]]=root;
  mask[root]=klev;//;
  wo[jbfs[ibfs[0]]]=-1.;// there is no edge associated with the root of the tree. 
  fprintf(stdout,"ia,ja %i:  %i, root=%i\n",ia[root],ia[root+1],root+1);
  ipoint=ibfs[1];
  kbeg=ibfs[0];
  kend=ibfs[1];
  while(1) {
    for(i1=kbeg;i1<kend;++i1){
      i=jbfs[i1];
      iai = ia[i];
      iai1=ia[i+1];
      ih=iai1-iai-1;
      if(!ih) continue;
      for(k=iai;k<iai1;++k){
	j=ja[k];
	//	fprintf(stdout,"%i : %i\n",j, ia[j+1]-ia[j]);
	if(i==j) continue;
	if(!mask[j]){
	  jbfs[ipoint]=j;
	  mask[j]=klev;
	  // bfs tree edge
	  et[j]=i; //parent of j is i.
	  wo[ipoint]=w[ledge[k]];
	  ipoint++;
	}
      }	   
    }
    //    fprintf(stdout,"level=%i; number of vertices: %i;
    //    test=%i\n", klev,kend-kbeg,ipoint-kend);
    if(ipoint==ibfs[klev]) break;
    kbeg=kend;
    klev++;
    ibfs[klev]=ipoint;
    kend=ipoint;
    if(ipoint>=nv)  break;
  }
  //  fprintf(stdout, "\nA1=[");
  //  for(i=0;i<ipoint;i++){
  //    j=jbfs[i];
  //    if(et[j]<0) continue;
    //        fprintf(stdout,
    //    	    "leaf: %i, node=%i with mask %i;  parent=%i;, weight=%12.4e\n",
    //    	    i+1,j,mask[j],et[j],wo[i]);
    //    fprintf(stdout, "%7i %7i %16.8e\n",j+1,et[j]+1,wo[i]);
    //fprintf(stdout, "%7i %7i %16.8e\n",et[j]+1,j+1,wo[i]);
    //  }
  //   fprintf(stdout, "];\n");
  *lev=klev;
  return;
}
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
  INT i,pos;
  //  INT j,k;
  // short hand;
  iCSRmat *dfs=malloc(1*sizeof(iCSRmat));
  dfs[0]=icsr_create(n,n,n);
  //  INT *mask=dfs->val;
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
/*********************************************************************************************************/
/*ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ*/
/***********************************************************************************************/
/*!
 * \fn void bfs00(const INT croot,iCSRmat *a, iCSRmat *bfs,INT *et, INT *mask)
 *
 * \brief bfs search of a graph corresponding to
 *
 * \param n                number of vertices
 * \param a:               adjacency matrix of the graph;
 * \param croots:          the starting vertex
 *
 * \param mask[]:
 *
 * \return On return iCSRmat *bfs is a iCSRmat of size n,n with n nonzeroes, which
 *         contains the bfs ordering permutation of the vertices
 *
 * \author Ludmil Zikatanov
 * \date   20210516
 */
void bfs00(const INT croot,			\
	   iCSRmat *a, iCSRmat *bfs,		\
	   INT *et, INT *mask)
{
  /* bfs for the connected component of the graph "a" containing
   * root. If root > number of vertices,
   * then root is chosen to be the max degree vertex.  on exit: the
   * number of rows in the returned icsrmat equals the number of
   * levels in bfs; the number of columns equals the number of rows in
   * a.
*/
  INT *ia=a->IA, *ja=a->JA;
  INT root=croot,nv=a->row,i,ih;
  if(root<0){
    fprintf(stderr,"ERROR: could not use %d as a root for bfs\n",root);
    exit(129);
  }
  et[root]=-1;
  INT *ibfs=bfs->IA, *jbfs=bfs->JA;
  INT i1,j,k,iai,iai1,ipoint,klev,kbeg,kend;
  /*************************************************************/
  // initialization
  klev=1; //level number ; for indexing this should be klev-1;
  ibfs[0]=0;
  ibfs[1]=1;
  jbfs[ibfs[0]]=root;
  mask[root]=klev;//;
  //  fprintf(stdout,"ia,ja %i:  %i, root=%i\n",ia[root],ia[root+1],root);
  ipoint=ibfs[1];
  kbeg=ibfs[0];
  kend=ibfs[1];
  while(1) {
    for(i1=kbeg;i1<kend;++i1){
      i=jbfs[i1];
      iai = ia[i];
      iai1=ia[i+1];
      ih=iai1-iai-1;
      //      fprintf(stdout,"ih=%d,iai=%d,iai1=%d;i=%d\n",ih,iai,iai1,i);
      if(ih<0) continue; //diagonals only or empty rows are ignored;
      for(k=iai;k<iai1;++k){
	j=ja[k];
	if(i==j) continue; // no self edges are counted;
	//	fprintf(stdout,"\ni=%d,j=%d,mask[%d]=%d; ipoint=%d",i,j,j,mask[j],ipoint);
	if(!mask[j]){
	  jbfs[ipoint]=j;
	  mask[j]=klev+1;
	  et[j]=i; // tree back edge pointing to ancestor
	  ipoint++;
	}
      }
    }
    if(ipoint==ibfs[klev]){
      /* fprintf(stdout,"\nexiting before the end: ipoint=%d\n;",ipoint);fflush(stdout); */
      /* print_full_mat_int(1,klev,ibfs,"ibfs1"); */
      /* print_full_mat_int(1,ibfs[klev],jbfs,"jbfs1"); */
      break;
    }
    kbeg=kend;
    klev++;
    ibfs[klev]=ipoint;
    kend=ipoint;
    if(ipoint>=nv){
      /* fprintf(stdout,"\nexiting at the end: ipoint=%d\n;",ipoint); */
      /* fflush(stdout); */
      /* print_full_mat_int(1,klev,ibfs,"ibfs2"); */
      /* print_full_mat_int(1,ibfs[klev],jbfs,"jbfs2"); */
      break;
    }
  }
  bfs->row=klev;
  bfs->col=nv;
  bfs->nnz=ibfs[klev];
  return;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx*/
/***********************************************************************/
iCSRmat *bfscc(INT nblk,INT *iblk, INT *jblk, iCSRmat *a, INT *et)
{
  /*
     bfs for the graph with possibly multiple connected components
     (cc); uses bfs00 to traverse every connected component. Also can
     be run for one connected component.
  */
  INT *ia=a->IA;
  INT nv=a->row;
  INT root,lvlend,bfsend,nvblk,i0,ib,ijb,i,ih,maxdeg,iablki,iablki1;
  iCSRmat *bfs=malloc(1*sizeof(iCSRmat));
  bfs[0]=icsr_create(nv+nblk,nv,nv);
  iCSRmat *bfstmp=malloc(1*sizeof(iCSRmat));
  bfstmp[0]=icsr_create(0,0,0);
  bfstmp->IA=bfs->IA;
  bfstmp->JA=bfs->JA;
  INT *mask=bfs->val;// this will be the val in the bfs.
  bfsend=0;
  lvlend=0;
  //  INT nblk1=nblk-1,iablki11=-22;
  for(ib=0;ib<nblk;ib++){
    iablki=iblk[ib];
    iablki1=iblk[ib+1];
    //    fprintf(stdout,"\nnodes=%d",iablki1-iablki);
    for(ijb=iablki;ijb<iablki;ijb++){
      //      fprintf(stdout,"\nzzzzz=jblk[%d]=%d",ijb,jblk[ijb]);
      mask[jblk[ijb]]=0;
      et[jblk[ijb]]=-1;
    }
    maxdeg=-1;root=-1;
    for(ijb=iablki;ijb<iablki1;ijb++){
      i=jblk[ijb];
      ih=ia[i+1]-ia[i];
      if(maxdeg<ih){
	maxdeg=ih;
	root=i;
      }
      //      root=i;
      //      break;
    }
    fprintf(stdout,"\nblock=%d;root=%d",ib,root);
    //    root=jblk[iablki];
    nvblk=iablki1-iablki;
    bfstmp->row=nvblk;
    bfstmp->col=nvblk;
    bfstmp->nnz=nvblk;
    bfs00(root,a,bfstmp,et,mask);
    //    print_full_mat_int(1,a->row,mask,"mask");
    //    print_full_mat_int(1,a->row,et,"etree");
    //    print_full_mat_int(1,bfstmp->IA[bfs->row],bfstmp->JA,"bfstmpja");
    bfstmp->IA+=(bfstmp->row+1);
    bfstmp->JA+=(bfstmp->nnz);
    bfsend+=bfstmp->nnz;
    lvlend+=(bfstmp->row+1);
  }
  //  print_full_mat_int(1,lvlend,bfs->IA,"bfsia");
  ib=1;
  for(i=1;i<lvlend;i++){
    if(bfs->IA[i]) {
      bfs->IA[ib]=bfs->IA[ib-1]+(bfs->IA[i]-bfs->IA[i-1]);
      ib++;
    }
  }
  ib--;
  //  fprintf(stdout,"\nlvlend=%d; ib=%d",lvlend,ib);
  bfs->row=ib;
  bfs->col=nv;
  bfs->nnz=bfs->IA[bfs->row];
  /**/
  //  print_full_mat_int(1,(bfs->row+1),bfs->IA,"bfsia");
  //  print_full_mat_int(1,bfs->nnz,bfs->JA,"bfsja");
  /**/
  bfs->IA=realloc(bfs->IA,(bfs->row+1)*sizeof(INT));
  bfs->JA=realloc(bfs->JA,(bfs->nnz)*sizeof(INT));
  bfs->val=realloc(bfs->val,(bfs->nnz)*sizeof(INT));
  /**/
  i0=0;
  for(i=0;i<bfs->row;i++){
    ih=i0;
    for(ijb=bfs->IA[i];ijb<bfs->IA[i+1];ijb++){
      ib=bfs->JA[ijb];
      if(et[ib]<0) {i0=i;}
      //      fprintf(stdout,"\netree[%d]=%d; mask[%d]=%d; mmm=%d ; i0=%d, ih=%d", ib,et[ib],ib,mask[ib],mask[ib]+i0,ih);
      mask[ib]+=i0;
    }
  }
  /* fprintf(stdout,"\nbfs0=["); */
  /* icsr_print_matlab_val(stdout,bfs); */
  /* fprintf(stdout,"];bfs=sparse(bfs0(:,1),bfs0(:,2),bfs0(:,3),%d,%d);\n\n",bfs->row,bfs->col); */
  return bfs;
}
/* /\***********************************************************************************************\/ */
/* /\*! */
/*  * \fn dfs00_(INT *nin, INT *ia, INT *ja, INT *nblko,INT *iblk, INT *jblk) */
/*  * */
/*  * \brief non-recurrsive dfs search of a graph.  This function */
/*  *        implements Tarjan's DFS algorithm for finding the strongly */
/*  *        connected components of a di-graph.  It is the FRED */
/*  *        GUSTAVSON'S non-recurrsive version of the algorithm.   */
/*  * */
/*  * \param n                number of vertices */
/*  * \param (ia,ja):         adjacency structure of the graph; all diagonal */
/*  *                         elements must be nonzero in ia,ja. */
/*  * \param nblko:           number of connected components */
/*  * \param (iblk,jblk):     CSR matrix containing the permutation: vertices */
/*  *                         are ordered ccomponent by compoent. */
/*  * */
/*  * \return The new numbering is in CSR format: the pointers are in */
/*  *         IBLK[], and the actual numbering is in JBLK[]. */
/*  * */
/*  * \author Ludmil Zikatanov */
/*  * \date   20210516 */
/*  *\/ */
/* void dfs00_(INT *nin, INT *ia, INT *ja, INT *nblko,INT *iblk, INT *jblk) */
/* { */
/*   INT n,n1,nb,nblk,count,k1,sp,vp,v=0,wp=0,w=0,v1=0,i,ii,myflag; */
/*   INT *numb=NULL,*iedge=NULL,*lowlink=NULL; */
/*   /\*--------------------------------------------------------------------- */
/*     ---------------------------------------------------------------------*\/ */
/*   /\*C...  Initialization.*\/ */
/*   n=*nin; */
/*   n1 = n+1; */
/*   iedge = (INT *) calloc(n1,sizeof(INT)); */
/*   numb = (INT *) calloc(n1,sizeof(INT)); */
/*   lowlink= (INT *) calloc(n1,sizeof(INT)); */
/*   /\* Check *\/ */
/*   if(!(numb && iedge && lowlink)) { */
/*     fprintf(stderr,"\nCould not allocate local mem in dfs\n"); */
/*     exit(19); */
/*   } */

/*   nblk = -1; */
/*   nb = -1; */
/*   count  = nb+1; */
/*   k1 = count; */
/*   vp = n;// */
/*   sp = n;// may be n??? */
/*   for (i=0; i<n;i++) { */
/*     iedge[i] = ia[i]; */
/*     numb[i]=-1; */
/*     lowlink[i]=-1; */
/*   } */
/*   numb[n] = -1; */
/*   lowlink[n] = -1; */
/*   myflag=10; */
/*   while (1) { //10 continue */
/*     //    fprintf(stdout," myflag=%i %i, %i %i\n",myflag,nblk,count,n);fflush(stdout); */
/*     if(myflag==10){ */
/*   /\* */
/*     C...  Get out when the  renumbering is done; */
/*     C...  Otherwise begin at vertex K1  */
/*   *\/ */
/*       if(count==n) { */
/* 	nblk++; //not sure why.... */
/* 	iblk[nblk] = n; */
/* 	*nblko=nblk; */
/* 	if(iedge) free(iedge); */
/* 	if(numb) free(numb); */
/* 	if(lowlink) free(lowlink); */
/* 	return; */
/*       } */
/*       for (ii=k1;ii<n;ii++) { */
/* 	i=ii; */
/* 	if(numb[ii]<0) { */
/* 	  myflag=30; */
/* 	  break; //go to 30 */
/* 	} */
/*       } */
/*       if(myflag !=30){ */
/* 	fprintf (stderr,"There is an error in DEPTH FIRST SEARCH:  %i == %i\n",i,n); */
/* 	exit(254); */
/*       } else { */
/* 	v = i; */
/* 	k1 = v + 1; */
/* 	myflag=50; */
/*       } */
/*     } */
/*     /\*C...  :::*\/ */
/*     if(myflag==40) { */
/*       vp--; */
/*       iblk[vp] = v; */
/*       v=w; */
/*       myflag=50; */
/*     } */
/*     if(myflag==50) { */
/*       nb++; */
/*       numb[v] = nb; */
/*       lowlink[v] = numb[v]; */
/*       sp--; */
/*       jblk[sp] = v; */
/*       v1 = v+1; */
/*       myflag=60; */
/*     } */
/*     /\*...  *\/ */
/*     while(60) { // 60   continue; */
/*       if(myflag == 60) { */
/* 	wp = iedge[v]; */
/* 	w = ja[wp]; */
/* 	//	fprintf(stdout,"\nv=%d myflag=%i wp=%i, count=%i n=%i\n",v,myflag,wp,count,n);fflush(stdout); */
/* 	iedge[v] = wp+1; */
/* 	if(numb[w] >= numb[v])  { */
/* 	  myflag=70; */
/* 	} else if (numb[w]<0) { */
/* 	  myflag=40;  */
/* 	  break; */
/* 	} else { */
/* 	  if(lowlink[v]>numb[w]) */
/* 	    lowlink[v]=numb[w]; */
/* 	  myflag=70; */
/* 	} */
/*       } */
/*       if(iedge[v] < ia[v1]) { */
/* 	myflag=60; */
/* 	continue; */
/*       } */
/*       /\*...*\/   */
/*       if(lowlink[v] >= numb[v]) {//don't {go to 90} */
/* 	nblk++; */
/* 	iblk[nblk] = count; */
/* 	/\**\/ */
/* 	while(80) {     //    80 continue; */
/* 	  w = jblk[sp]; */
/* 	  numb[w] = n; */
/* 	  sp++; */
/* 	  jblk[count] = w; */
/* 	  count++; */
/* 	  if(v == w) break; //{don't go to 80} */
/* 	} */
/* 	/\*C...  *\/ */
/* 	if(sp==n) { */
/* 	  myflag=10;  */
/* 	  break; */
/* 	} */
/*       }  */
/*       myflag=70; */
/*       w = v; */
/*       v = iblk[vp]; */
/*       vp++; */
/*       v1 = v + 1; */
/*       if(lowlink[v]>lowlink[w])  */
/* 	lowlink[v]=lowlink[w]; */
/*     } */
/*   } */
/* } */
/* /\************************************************************\/ */
