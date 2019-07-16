/*! \file src/graphs/bfs.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note breadth first search and related routines. 
 *
 */
#include "hazmath.h"
/***********************************************************************/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx*/
/***********************************************************************/
void bfs00(const INT croot,			\
	   iCSRmat *a, iCSRmat *bfs,		\
	    INT *et, INT *mask)
{
  /* bfs for the connected component of the graph "a" containing
   * root. If root is negative or bigger than the number of vertices,
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
    }
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
  /* for(i=0;i<nv;i++){ */
  /*   iai=ia[i];iai1=ia[i+1]; */
  /*   for(k=iai;k<iai1;k++){ */
  /*     j=ja[k]; */
  /*     if(i==j)continue; */
  /*     fprintf(stdout, "%7i %7i %16.8e\n",i+1,j+1,w[ledge[k]]); */
  /*   } */
  /* } */
  /* fprintf(stdout, "];\n"); */
