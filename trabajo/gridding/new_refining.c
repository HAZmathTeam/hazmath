/*! \file src/amr/refining.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 * \note: routines used to refine a simplicial grid grid ref_level times.

 */
#include "hazmath.h"
//#include "grid_defs.h"
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
void abfstree(INT it, scomplex *sc,INT *wrk);
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
void n_refine(INT ref_type, INT ref_levels, scomplex *sc,	\
	      dvector *errors, 
	      void (*solving)(INT , scomplex *, void *),		\
	      void (*estimating)(INT , scomplex *, void *),	\
	      void (*marking)(INT , scomplex *, void *))
{
/* 
 *
 * refine ref_levels. Follows the algorithm:
 * SOLVE-->ESTIMATE-->MARK-->REFINE. 
 *
 */
  if(ref_levels<=0) return;
  /* 
     "anything" is a void pointer to structure or array or bunch of
     arrays of mixed type which are used in solve, estimate and mark
     phases.
  */
  void *anything = (void *)errors;
  /**/
  INT n=sc->n,i=-1,j=-1,k=-1,i123=-10;
  INT ns,nv,n1,nsold,nvold,level;
  if(!sc->level){
    /*form neighboring list; */
    find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
    //    haz_scomplex_print(sc,0,__FUNCTION__);  fflush(stdout);
    INT *wrk=calloc(5*(sc->n+2),sizeof(INT));
    /* construct bfs tree for the dual graph */
    abfstree(0,sc,wrk);
    //    haz_scomplex_print(sc,0,__FUNCTION__);fflush(stdout);
    //    exit(100);
    if(wrk) free(wrk);
  }
  n=sc->n; n1=n+1; level=0;
  fprintf(stdout,"refine: ");
  while(sc->level < ref_levels && TRUE){
    //    if(dxstar->row && dxstar->val) markstar(sc,dxstar);
    /* SOLVE */
    (*solving)(ref_type, sc, anything);
    /* ESTIMATE */
    (*estimating)(ref_type,sc, anything);
    /* MARK */
    (*marking)(ref_type, sc, anything);
      /* for(i = 0;i<sc->ns;i++){ */
      /* //    if(sc->gen[i] < level) continue;     */
      /* fprintf(stdout,"\n%s: ZZZZZZZZZZZzsimplex %d (mark=%d, gen=%d, c0=%d)", \ */
      /* 	      __FUNCTION__,i,sc->marked[i],sc->gen[i],sc->child0[i]); */
      /* } */
    /* REFINE FOLLOWS : */
    nsold=sc->ns;
    nvold=sc->nv;
    /* 
     * refine everything that is marked on the finest level and is
     * not yet refined: (marked>0 and child<0)
     */
    for(j = 0;j < nsold;j++)
      if(sc->marked[j] && (sc->child0[j]<0||sc->childn[j]<0))
	haz_refine_simplex(sc, j, -1);
    /* new mesh */
    ns=sc->ns; nv=sc->nv;
    sc->level++;
    fprintf(stdout,".%d.",sc->level);//,nsold,ns,nv);
  }
  fprintf(stdout,"\n");
  return;
}
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
  /* bfs tree rooted at simplex it */
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
void scfinalize(scomplex *sc)
{
  /* copy the final grid at position 1*/
  INT ns,j=-10,k=-10,n=sc->n,n1=sc->n+1;
  /*  
      store the finest mesh in sc structure. 
      on input sc has all the hierarchy, on return sc only has the final mesh. 
  */
  ns=0;
  for (j=0;j<sc->ns;j++){
    /*
      On the last grid are all simplices that were not refined, so
      these are the ones for which child0 and childn are not set. 
    */
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
  //  sc->nv=nv;
  fprintf(stdout,"\n%%After %d levels of refinement:\tsimplices=%d ; vertices=%d\n",sc->level,sc->ns,sc->nv); fflush(stdout);
  return;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
void sc2mesh(scomplex *sc,trimesh *mesh)
{
  /* copy the final grid at position 1*/
  INT ns,n=sc->n,n1=sc->n+1,jk=-10,k=-10,j=-10;
  /*  
      store the finest mesh in trimesh structure. 
      sc has all the hierarchy, trimesh will have only the last mesh. 
  */
  ns=0;
  for (j=0;j<sc->ns;j++){
    /*  On the last grid are all simplices that were not refined, so
      these are the ones for which child0 and childn are not set.  */
    if(sc->child0[j]<0 || sc->childn[j]<0){
      mesh->el_flag[ns]=sc->flags[j];
      mesh->el_vol[ns]=sc->vols[j];
      ns++;
    }
  }
  mesh->nv=sc->nv;
  for (j=0;j<mesh->nv;j++){
    mesh->v_flag[j]=sc->bndry[j];
  }
  mesh->el_v->IA = (INT *) calloc(ns+1,sizeof(INT)); 
  mesh->el_v->JA = (INT *) calloc(ns*n1,sizeof(INT));
  mesh->nelm=ns;
  mesh->el_v->IA[0]=0;   
  for (j=0;j<mesh->nelm;j++){
    mesh->el_v->IA[j+1]=mesh->el_v->IA[j]+n1;       
  }
  jk=0;
  for (j=0;j<sc->ns;j++){
    /*  copy el_v map using only the top grid;    */
    if(sc->child0[j]<0 || sc->childn[j]<0){
      for (k=0;k<n1;k++)
	memcpy(mesh->el_v->JA+jk*n1,sc->nbr+j*n1,n1*sizeof(INT));
      jk++;
  }
  fprintf(stdout,"\n%%After %d levels of refinement:\tsimplices=%d ; vertices=%d\n",sc->level,sc->nv,ns); fflush(stdout);
  }
  return;
}
