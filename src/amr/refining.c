/*! \file src/amr/refining.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 * \note: routines used to refine a simplicial grid grid ref_level times.

 */
#include "hazmath.h"
unsigned int reflect2(INT n, INT is, INT it,				\
		      INT* sv1, INT *sv2, INT* stos1, INT* stos2,	\
		      INT visited, INT *wrk)
/*************************************************************/
/* it works for all meshes that can be consistently ordered. */
/********************************************************************/
{
  /* 
     sv1 are the n+1 vertices of is; sv2 are the (n+1) vertices of
     (it).  
     
     stos1 are the n+1 neighbors of is; stos2 are the (n+1) neighbors of
     (it).
     
     This routine checks whether is is reflected neighbor of it
     and reorders is if it was not visited before One main assumption is
     that (is) and (it) intersect in (n-1) dimensional simplex with (n)
     vertices, so that only one element of sv1 (say, k1) is not present
     in sv2 and only one element of sv2 (say, k2) is not present in
     sv1. Then we make an important assumption here, that stos1[k1] = it
     and stos[k2]= is. This is always achievable when creating the stos
     (simplex to simplex) relation.

     wrk is working space of size n+2, n is the spatial dimension 
  */
  INT n1=n+1,n2=n+2;
  INT *wrk1=NULL,*invp=NULL,*p=NULL,*pw=NULL,*wrk2=NULL;
  INT i,j;
  /**/
  /* we also check if we have reflected neighbors */
  if(visited){
    for (i=0; i<n1;i++){
      if(stos1[i]!=it){
	if(sv1[i]-sv2[i]) {
	  /* not reflected neighbors */
	  fprintf(stderr,"\n***ERROR in %s ; (is)=%d(vis=%d) and (it) = %d are both visited but are not reflected neighbors.\n\n",__FUNCTION__,is,visited,it);
	  fprintf(stderr,"\n***The problem is at node %d, (sv1=%d,sv2=%d)\n",i,sv1[i],sv2[i]);
	  return 2;
	}
      }
    }
    /* we have reflected neighbors, so we just return */
    return 0;
  }
  INT kv1=-1,kv2=-1;
  for (i=0; i<n1;i++){
    if(stos1[i] == it){
      kv1=sv1[i];
      break;
    }      
  }
  for (i=0; i<n1;i++){
    if(stos2[i] == is){
      kv2=sv2[i];
      break;
    }      
  }
  if (kv1<0 || kv2 < 0) {
    fprintf(stderr,"\n***ERROR in %s ; kv1=%d, kv2=%d must be positive.\n\n",__FUNCTION__,kv1,kv2);
    return 3;
  }
  wrk1=wrk; wrk2=wrk1+n2; p=wrk2+n2;invp=p+n2;  pw=invp+n2;
  memcpy(wrk1,sv1,n1*sizeof(INT));wrk1[n1] = kv2; 
  isi_sortp(n2,wrk1,p,pw);
  /*  returrn wrk1 to the initial state */
  memcpy(wrk1,sv1,n1*sizeof(INT));wrk1[n1] = kv2; 
  /* second array*/
  memcpy(wrk2,sv2,n1*sizeof(INT)); wrk2[n1] = kv1;
  isi_sortp(n2,wrk2,pw,invp);
  /* return wrk2 to init state */
  memcpy(wrk2,sv2,n1*sizeof(INT)); wrk2[n1] = kv1;
  /*
    We use here identity1: sv1[p1[k]] = sv2[p2[k]] for all k. Hence we have:
    
    s2[j] = s2[p2[invp2[j]]] = s1[p1[invp2[j]]]
    
    where we used the identity1 with k=invp2[j]
    
  */
  /* now we can use p2 and wrk1 to move around what we need */
  for (i=0; i<n1;i++){
    j=p[invp[i]]; /* this is the new index of sv1[i] */
    if(wrk2[i] != kv2){
      sv1[i] = wrk1[j];
    } else {
      sv1[i]=kv1;
    }
  }
  for (i=0; i<n1;i++){
    j=p[invp[i]]; /* this is the new index of sv1[i] */
    if(wrk2[i] != kv2){
      wrk1[i] = stos1[j];
    } else {
      wrk1[i] = it;
    }
  }
  for (i=0; i<n1;i++){
    stos1[i] = wrk1[i];
  }
  return 0;
}
/******************************************************************/
/*using bfs to get the reflected mesh*/
void abfstree(INT it, scomplex *sc,INT *wrk) 
{
  /* from the scomplex the mask are the marked and the */
  /* bfs tree rooted at 0 */
  INT n=sc->n,n1=n+1,ns=sc->ns;
  INT i,j,k,iii,is,isn1,itn1;
  INT i1,in1,kbeg,kend,nums,iai,iai1,klev;
  INT *mask = sc->marked;
  INT *jbfs = sc->child0;
  INT ireflect=-10;
  // initialization
  /*************************************************************/
  INT *isnbr,*itnbr,*isv,*itv;
  /* caution: uses sc->marked as working array of size ns */
  memset(mask,0,ns*sizeof(INT));
  //  is=0;//(INT )(ns/2);
  nums=0;
  klev=1; //level number ; for indexing this should be klev-1;
  jbfs[nums]=it; // thit it an input simplex where to begin. 
  mask[it]=klev;
  nums++;
  kbeg=0; kend=1;
  while(1) {
    for(i1=kbeg;i1<kend;++i1){
      it=jbfs[i1];
      iai  = it*n1;
      iai1 = iai+n1;
      itn1 = it*n1;
      itnbr=(sc->nbr+itn1); 
      itv=(sc->nodes+itn1);
      for(k=iai;k<iai1;++k){
	is=sc->nbr[k];
	//	fprintf(stdout,"%i(nbr=%i) ",i,j);fflush(stdout);	  
	if(is<0) continue;
	isn1=is*n1;
	isnbr=(sc->nbr+isn1);
	isv=(sc->nodes+isn1);
	ireflect = reflect2(n,is,it,isv,itv,isnbr,itnbr,mask[is],wrk);
	switch(ireflect) {
	case 2 :
	case 3 :
	  exit(ireflect);
	  break;
	case 0 :
	  /* if it was visited and all is OK, we just return */
	  /* if(mask[it]) { */
	  /*   return 0; */
	  /* } else { */
	  /*   break; */
	  /* } */
	  break;
	default :
	  fprintf(stderr,						\
		  "Invalid return from reflect2 in %s; return value = %d\n", \
		  __FUNCTION__,ireflect);
	  exit(4);
	}
	if(!mask[is]){
	  jbfs[nums]=is;
	  mask[is]=klev;
	  //	  fprintf(stdout,"%i(%i,%i)",i,j,mask[j]);fflush(stdout);      
	  nums++;
	}
	
      }
    }
    kbeg=kend; kend=nums;klev++;
    if(nums >= ns)  break;
  }
  //  fprintf(stdout,"%%BFS levels for reflect: %d; ",klev-1);
  for(i=0;i<ns;i++){
    jbfs[i]=-1;
  }
  return;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
void refining(INT ref_levels, scomplex *sc, INT nstar, REAL *xstar)
{
/* ATTN: it refines ref_levels times and then sets ref_levels=0 so
 *  when called again it will not refine automatically unless
 *  ref_levels is reset to new value; 
 *
 *xstar is an array nstar by dim
 *  which contains coordinates of points where we want to refine.
 *  here n is used instead of dim.  bail out if no refinement is
 *  required
*/
  if(ref_levels<=0) return;
  INT n=sc->n,j=-1,k=-1;
  INT ns,nv,n1,nsold,nvold,level;
  /* working space used in reflect_mesh */
  INT *wrk=(INT *)calloc(5*(sc->n+2),sizeof(INT));  
  // longest edge?
  abfstree(0,sc,wrk);
  n=sc->n; n1=n+1; level=0;
  /* we first mark everything */
  for (j=0;j<sc->ns;j++) {sc->marked[j]=1;}
  fprintf(stdout,"refine: ");
  while(level < ref_levels){
    nsold=sc->ns;
    nvold=sc->nv;
    /* 
       mark some simplices. if nstar=0 nothing is marked additionally,
       i.e. whatever was marked and its children stays marked: the
       refinement is automatic up to level=ref_levels.
    */
    if(nstar)
      markstar(level,sc,nstar,xstar);
    else
      marks(level,sc);
    for(j = 0;j < nsold;j++) {
      if(sc->marked[j] && (sc->childn[j]<0)){
	haz_refine_simplex(sc, j, -1);
      }
    }
    /* new mesh */
    ns=sc->ns; nv=sc->nv;
    level++;
    fprintf(stdout,".%d.",level);
  }
  fprintf(stdout,"\n");
  ns=0;
  /*  
      remove hierarchy, only last sc is needed for computations...
  */
  for (j=0;j<sc->ns;j++){
    if(sc->child0[j]<0 || sc->childn[j]<0){
      for (k=0;k<n1;k++) {
	sc->nodes[ns*n1+k]=sc->nodes[j*n1+k];
	sc->nbr[ns*n1+k]=sc->nbr[j*n1+k];
      }
      sc->gen[ns]=sc->gen[j];
      sc->flags[ns]=sc->flags[j];
      ns++;
    }
  }
  sc->ns=ns;
  sc->nv=nv;
  fprintf(stdout,"%%After %d levels of refinement:\tsimplices=%d ; vertices=%d\n",ref_levels,sc->ns,sc->nv); fflush(stdout);
  if(wrk) free(wrk);
  return;
}
