/*! \file src/utilities/amr_utils.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190115.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note contains some utility functions for mesh refinement.
 *
 */
#include "hazmath.h"
/**********************************************************************/
/*!
 * \fn INT aresame(INT *a, INT *b, INT n)
 *
 * \brief checks if two arrays have same elements up to a permutation.
 *
 * \param a:   input array
 * \param b:   input array to compare with a.
 * \param n:   the size of a and b;
 *
 * \return
 *
 * \note
 *
 */
INT aresame(INT *a, INT *b, INT n)
{
  /*
     checks (n^2 algorithm) if two have the same elements (up to a
     permutation), if they are returns 1, otherwise 0
  */
  INT i,j,flag,ai,bj;
  for (i=0;i<n;i++){
    ai=a[i];
    flag=0;
    for(j=0;j<n;j++){
      bj=b[j];
      if(ai==bj){
	flag=1; break;
      }
    }
    if(!flag) return 0;
  }
  return 1;
}
/**********************************************************************/
/*!
 * \fn INT aresamep(INT *a, INT *b, INT n, INT *p)
 *
 * \brief checks (n^2 algorithm) if two have the same elements (up to a
     permutation);
 *
 * \param a:   input array
 * \param b:   input array to compare with a.
 * \param n:   the size of a and b;
 * \param p:   the permutation which takes a into b if they are the same.
 *
 *
 * \return returns 0 if the arrays are not the same; returns 2 if they
 *                    are the same and p is the permutation
 *                    so that a[i]=b[p[i]]. If there is no permutation,
 *                    i.e. p[i]=i, then returns 1.
 *
 *
 * \note
 *
 */
INT aresamep(INT *a, INT *b, INT n, INT *p)
{
  INT i,j,ai,bj;
  INT flag=-1,iret=1;
  for (i=0;i<n;i++)p[i]=-1;
  for (i=0;i<n;i++){
    ai=a[i];
    flag=0;
    for(j=0;j<n;j++){
      bj=b[j];
      if(ai==bj){
	p[i]=j;
	flag=1;
	if(p[i]!=i) iret=2;
	break;
      }
    }
    if(!flag) return 0;
  }
  return iret;
}
/**********************************************************************/
/*!
 * \fn INT xins(INT n, INT *nodes, REAL *xs, REAL *xstar)
 *
 * \brief     In dimension "n" constructs the map from reference simplex to
 *    simplex with coordinates xs[0..n].  Then solves a linear system
 *    with partial pivoting to determine if a point given with
 *    coordinates xstar[0..n-1] is in the (closed) simplex defined by
 *    "nodes[0..n] and xs[0..n]"
 *
 * \param n:         dimension of the simplex
 * \param nodes:     global numbering of the vertices of the simplex
 *
 * \param xs[]:      coordinates of all vertices in the simplicial
 *                   complex. This corresponds to the global
 *                   numbering. So the xs[nodes[0...(n+1)]]=coords of
 *                   the simplex vertices
 *
 * \param xstar[]:   node for which we would like to check its incidence with the simplex.
 *
 * \param nbig:       dimension of the space where the simplex is embedded
 * \return  0: the point is in the simplex; nonzero: not in the simplex.
 *
 * \note
 *
 */
INT xins(INT n, INT *nodes, REAL *xs, REAL *xstar)
{
  // does not support different dimensions of the simplex and the
  // space containing it.
  INT n1=n+1,i,j,l0n,ln,j1;
  INT *p=NULL;
  REAL *A=NULL,*xhat=NULL, *piv=NULL;
  A=(REAL *)calloc(n*n,sizeof(REAL));
  xhat=(REAL *)calloc(n,sizeof(REAL));
  piv=(REAL *)calloc(n,sizeof(REAL));
  p=(INT *)calloc(n,sizeof(INT));
  /* fprintf(stdout,"\nj=%d; vertex=%d\n",0+1,nodes[0]+1); fflush(stdout); */
  /* fprintf(stdout,"\nvertex=%d\n",nodes[0]+1); */
  /* for(i=0;i<n;i++){ */
  /*   fprintf(stdout,"xyz=%f ",xs[l0n+i]); */
  /* } */
  /* fprintf(stdout,"\n"); */
  l0n=nodes[0]*n;
  for (j = 1; j<n1;j++){
    /* grab the vertex */
    ln=nodes[j]*n;
    j1=j-1;
    for(i=0;i<n;i++){
      /*A_{ij} = x_{i,nodes(j)}-x0_{i,nodes(j)},j=1:n (skip 0)*/
      A[i*n+j1] = xs[ln+i]-xs[l0n+i];
    }
  }
  //  fflush(stdout);
  for(i=0;i<n;i++)
    xhat[i] = xstar[i]-xs[l0n+i];
  ddense_solve_pivot(1, n, A, xhat, p, piv);
  REAL xhatn=1e0,eps0=1e-10,xmax=1e0+eps0;
  /* check the solution if within bounds */
  INT flag = 0;
  for(j=0;j<n;j++){
    if((xhat[j] < -eps0) || (xhat[j] > xmax)){
      flag=(j+1);
      //      fprintf(stdout,"\nNOT FOUND: xhat(%d)=%e\n\n",flag,xhat[j]);
      break;
    }
    xhatn -= xhat[j];
    if((xhatn<-eps0) || (xhatn>xmax)) {
      flag=n+1;
      //          fprintf(stdout,"\nNOT FOUND: xhat(%d)=%e\n\n",flag,xhatn);
      break;
    }
  }
  if(A) free(A);
  if(xhat) free(xhat);
  if(p) free(p);
  if(piv) free(piv);
  return flag;
}



/**********************************************************************/
/*!
 * \fn void marks(scomplex *sc,dvector *errors)
 *
 * \brief marks simplices based on input vector with errors.
 *
 * \param sc:            simplicial complex
 * \param errors[]:      dvector with errors.
 *
 * \return               changes sc->marked[]: simplices marked
 *                       for refinemend have sc->marked[s]=true;
 *
 * \note
 *
 */
void marks(scomplex *sc,dvector *errors)
{
  /* mark simplices depending on the value of an estimator */
  /* the estimator here is the aspect ratio of the simplex */
  INT n=sc->n,nbig=sc->nbig,n1=n+1,ns=sc->ns,level=sc->level;
  //INT kbadel;
  INT ke,i,j,j1,k,p,ni,mj,mk;
  INT ne=(INT )((n*n1)/2);
  REAL slmin,slmax,asp,aspmax=-10.;;
  REAL *sl=(REAL *)calloc(ne,sizeof(REAL));
  INT kbad=0;
  for(i = 0;i<ns;i++){
    if(sc->gen[i] < level) continue;
    ni=n1*i;
    ke=0;
    for (j = 0; j<n;j++){
      mj=sc->nodes[ni+j]*n;
      j1=j+1;
      for (k = j1; k<n1;k++){
	mk=sc->nodes[ni+k]*n;
	sl[ke]=0e0;
	for(p=0;p<nbig;p++){
	  sl[ke]+=(sc->x[mj+p]-sc->x[mk+p])*(sc->x[mj+p]-sc->x[mk+p]);
	}
	sl[ke]=sqrt(sl[ke]);
	ke++;
      }
    }
    slmin=sl[0];
    slmax=sl[0];
    for(j=1;j<ne;j++){
      if(sl[j]<slmin) slmin=sl[j];
      if(sl[j]>slmax) slmax=sl[j];
    }
    asp=slmax/slmin;
    if(asp>1e1){
      sc->marked[i]=1;
      kbad++;
      //      fprintf(stdout,"\nlev=%d, gen=%d, asp= %e(%e/%e)\n",level,sc->gen[i],asp,slmax,slmin);
      if(asp>aspmax){
        //kbadel=i;
        aspmax=asp;
      }
    }
  }
  //  fprintf(stdout,"\nbad:%d, aspectmax=%e (at el=%d)\n",kbad,aspmax,kbadel);
  //  exit(33);
  if(sl)free(sl);
  return;
}
/******************************************************************************/
/*!
 * \fn unsigned int reflect2(INT n, INT is, INT it, INT* sv1, INT
 *		      *sv2, INT* stos1, INT* stos2, INT visited, INT *wrk)
 *
 * \brief         This routine checks is and it are reflected neighbors
 *                reorders is if it was not visited before. One main
 *                assumption is that (is) and (it) intersect in (n-1)
 *                dimensional simplex with (n) vertices, so that only
 *                one element of sv1 (say, k1) is not present in sv2
 *                and only one element of sv2 (say, k2) is not present in
 *                sv1. Then we make an important assumption here, that
 *                stos1[k1] = it and stos[k2]= is. This is always
 *                achievable when creating the simplex-to-simplex map.
 *
 * \param sv1[]:       the n+1 vertices of (is);
 * \param sv2[]        the (n+1) vertices of (it).
 *
 * \param stos1[]:     are the n+1 neighbors of (is); stos2 are the
 *                     (n+1) neighbors of (it).
 *
 *  \param             wrk[] is working space of size n+2,
 *                     n is the spatial dimension
 *
 * \return             returns 0 if is and it are reflected neighbors
 *
 * \note
 *
 */
unsigned int reflect2(INT n, INT is, INT it,				\
		      INT* sv1, INT *sv2, INT* stos1, INT* stos2,	\
		      INT visited, INT *wrk)
/*************************************************************/
/* it works for all meshes that can be consistently ordered. */
/********************************************************************/
{
  /*
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
/**********************************************************************/
/*!
 * \fn void abfstree(const INT it0, scomplex *sc,INT *wrk,const INT
 *                   print_level)
 *
 * \brief  uses the simplex-to-simplex map to create a bfs tree. Then
 *         the edges of the tree are followed to try to consistently
 *         order the locally (see reflect2). bfs tree: constructs all
 *         bfs trees for each connected component of the element
 *         neighboring list in the connected component containing it;
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
/*using bfs to get the reflected mesh*/
void abfstree(const INT it0, scomplex *sc,INT *wrk,const INT print_level)
{
  //  haz_scomplex_print(sc,0,__FUNCTION__);  fflush(stdout);
  INT it=it0, n=sc->n,n1=n+1,ns=sc->ns,cc=sc->cc;
  INT i,j,k,iii,is,isn1,itn1;
  INT i1,kcc;
  //INT in1;
  INT kbeg,kend,nums,iai,iai1,klev;
  iCSRmat neib=icsr_create(ns,ns,(n+1)*ns);
  iii=0;
  neib.IA[0]=iii;
  for(i=0;i<ns;i++){
    isn1=i*n1;
    for(j=0;j<n1;j++){
      is=sc->nbr[isn1+j];
      if(is>=0){
	///	fprintf(stdout,"\niii=%d,ns=%d,elemnts=%d(=%d?)",iii,ns,neib.IA[ns],neib.nnz); fflush(stdout);
	neib.JA[iii]=is;
	iii++;
      }
    }
    neib.IA[i+1]=iii;
  }
  neib.nnz=neib.IA[neib.row];
  neib.JA=realloc(neib.JA,neib.nnz*sizeof(INT));
  neib.val=realloc(neib.val,neib.nnz*sizeof(INT));
  //  fprintf(stdout,"\nns=%d, nnz_ia=%d,nnz=%d\n",ns,neib.IA[neib.row],neib.nnz); fflush(stdout);
  // assuming neib.val has more than 2*num_simplices
  INT *mask,*jbfs;
  if(neib.nnz<(2*ns+1)) {
    neib.val=realloc(neib.val,(2*ns+1)*sizeof(INT));
  }
  mask=neib.val;
  jbfs = mask+ns;
  // find the connected components.
  //zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
  //  fprintf(stdout,"\nNNNNNNNNNNNNNNNN(%d)=%s",ns,__FUNCTION__); fflush(stdout);
  iCSRmat *blk_dfs=run_dfs(ns,neib.IA, neib.JA);
  /********************************************************************************************************/
  /* fprintf(stdout,"\nns=%d,nnz_ia=%d,nnz=%d",ns,neib.IA[ns],neib.nnz); fflush(stdout); */
  /* fprintf(stdout,"\nneib1=["); */
  /* icsr_print_matlab(stdout,neib); */
  /* fprintf(stdout,"];"); */
  /* fprintf(stdout,"\nneib=sparse(neib1(:,1),neib1(:,2),neib1(:,3));\n"); */
  /********************************************************************************************************/
  cc=blk_dfs->row;
  //zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
  INT ireflect;
  // initialization
  /*************************************************************/
  INT *isnbr,*itnbr,*isv,*itv;
  memset(mask,0,ns*sizeof(INT));
  //  is=0;//(INT )(ns/2);
  for(kcc=0;kcc<cc;kcc++){
    ireflect=-10;
    it=blk_dfs->JA[blk_dfs->IA[kcc]];
    if(print_level>5){
      fprintf(stdout,
	      "\n%s: Component=%d; root=%d;\n",__FUNCTION__,kcc,it);fflush(stdout);
    }
    nums=0;
    klev=1; //level number ; for indexing this should be klev-1;
    jbfs[nums]=it; // this is an input simplex where to begin.
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
	    fprintf(stderr,						\
		    "Invalid return from reflect2 in %s; return value = %d\n", \
		    __FUNCTION__,ireflect);
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
      //    fprintf(stdout,"\nkbeg=%i,kend=%i,nums=%i",kbeg,kend,nums);fflush(stdout);
      if(kend>=nums) break;/* exit here even if we have not visited all simplices as they could be disconnected */
      kbeg=kend; kend=nums;klev++;
      // this below only works if the domain is connected;
      if(nums >= ns)  break;
    }
  }
  //  fprintf(stdout,"%%BFS levels for reflect: %d; ",klev-1);
  //  for(i=0;i<ns;i++){
  //    jbfs[i]=-1;
  //  }
  icsr_free(&neib);
  icsr_free(blk_dfs);free(blk_dfs);
  return;
}
/******************************************************************/
/*!
 * \fn scomplex *scfinest(scomplex *sc)
 *
 * \brief
 *
 * \param sc: scomplex containing the whole hierarchy of refinements
 *
 * \return the simplicial complex corresponding to all simplices which
 *         were not refined.
 *
 */
scomplex *scfinest(scomplex *sc)
{
  INT ns,i=0,j=-10,k=-10,n=sc->n,nbig=sc->nbig,n1=sc->n+1,nv=sc->nv;
  scomplex *sctop=NULL;
  /*
      store the finest mesh in and return the sc structure. save the
      correspondence between elements in an sc->child0[] as a negative
      number.  sc has all the hierarchy, on return sctop only has
      only the final mesh.
  */
  /*first step: compute the number of simplices on the final level */
  ns=0;
  for (j=0;j<sc->ns;j++){
    /* On the last grid are all simplices that were not refined,
       so these are the ones for which child0 and childn are not
       set. */
    if(sc->child0[j]<0 || sc->childn[j]<0)ns++;
  }
  /*allocate the scomplex for the finest level*/
  sctop=haz_scomplex_init(n,ns,nv,nbig);
  /* we dont need bunch of these, at least for assembly, so we free them*/
  free(sctop->parent);sctop->parent=NULL;
  free(sctop->childn);sctop->childn=NULL;
  free(sctop->child0);sctop->child0=NULL;
  ns=0;
  for (j=0;j<sc->ns;j++){
    if(sc->child0[j]<0 || sc->childn[j]<0){
      sc->child0[j]=-(ns+1); //save this for future reference;
      for (k=0;k<n1;k++) {
	sctop->nodes[ns*n1+k]=sc->nodes[j*n1+k];
	sctop->nbr[ns*n1+k]=sc->nbr[j*n1+k];
      }
      sctop->flags[ns]=sc->flags[j];
      sctop->marked[ns]=sc->marked[j];// making sure nothing is marked on the top for refinement;
      sctop->vols[ns]=sc->vols[j];// making sure nothing is marked on the top for refinement;
      ns++;
    }
  }
  //  if(sctop->ns!=ns) {
    /*    issue an error here and stop */
  //  }
  /* connected components, these should not change */
  sctop->cc=sc->cc;
  sctop->bndry_cc=sc->bndry_cc;
  /* copy the boudary codes and the coordinates*/
  for(i=0;i<nv;i++){
    sctop->bndry[i]=sc->bndry[i];
    for(j=0;j<nbig;j++)
      sctop->x[i*n+j]=sc->x[i*nbig+j];
  }
  free(sctop->csys);sctop->csys=NULL;
  return sctop;
}
/**********************************************************************/
/*!
 * \fn void scfinalize(scomplex *sc,const INT set_bndry_codes)
 *
 * \brief Removes all hierachy and make sc to represent only the final
 *        grid. computes connected components and connected components
 *        on the boundary.
 *
 * \param sc: simplicial complex
 * \param set_bndry_codes: if set to 1, all boundary vertices get a code of 128+(connected component number);
 *
 * \return
 *
 * \note
 *
 */
void scfinalize(scomplex *sc,const INT set_bndry_codes)
{
  // INT n=sc->n;
  INT ns,j=-10,k=-10;
  INT n1=sc->n+1;
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
      }
      sc->child0[ns]=-1;
      sc->childn[ns]=-1;
      sc->gen[ns]=sc->gen[j];
      sc->flags[ns]=sc->flags[j];
      ns++;
    }
  }
  sc->ns=ns;
  sc->nodes=realloc(sc->nodes,n1*sc->ns*sizeof(INT));
  sc->nbr=realloc(sc->nbr,n1*sc->ns*sizeof(INT));
  sc->vols=realloc(sc->vols,sc->ns*sizeof(REAL));
  sc->child0=realloc(sc->child0,sc->ns*sizeof(INT));
  sc->childn=realloc(sc->childn,sc->ns*sizeof(INT));
  sc->gen=realloc(sc->gen,sc->ns*sizeof(INT));
  sc->flags=realloc(sc->flags,sc->ns*sizeof(INT));
  find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
  // this also can be called separately
  find_cc_bndry_cc(sc,set_bndry_codes);
  //
  if(set_bndry_codes){
    for(j=0;j<sc->nv;++j){
      if(sc->bndry[j]>128) sc->bndry[j]-=128;
    }
  }
  //fprintf(stdout,"\n%%After %d levels of refinement:\tsimplices=%d ; vertices=%d\n",sc->level,sc->ns,sc->nv); fflush(stdout);
  return;
}
/**********************************************************************/
/*!
 * \fn void cube2simp_free(cube2simp *c2s)
 *
 * \brief free all arrays in the cube-to-simplex structure
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void cube2simp_free(cube2simp *c2s)
{
  if(c2s->bits)free(c2s->bits);
  if(c2s->nodes)free(c2s->nodes);
  if(c2s->edges)free(c2s->edges);
  if(c2s->faces)free(c2s->faces);
  if(c2s->perms)free(c2s->perms);
  if(c2s)free(c2s);
  return;
}
/**********************************************************************/
/*!
 * \fn static void binary0(cube2simp *c2s)
 *
 * \brief stores in an array the coordinates of the vertices of the
 *        unit cube in dimension (dim). Lexicographical ordering from
 *        0,0,...,0 to 1,1,...,1
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
static void binary0(cube2simp *c2s)
{
  INT nvcube=c2s->nvcube;
  INT shift,i,j,k,kn,nbits=c2s->n-1;
  for(k = 0;k<nvcube;k++){
    kn=k*c2s->n;
    for (i=nbits ; i >=0; --i){
      c2s->bits[kn+i] = (unsigned INT )((k >> i & 1));
    }
  }
  shift=(1<<(c2s->n-1));
  INT nperm,jp=-22,jpo=-22,mid=(int)(c2s->nvcube/2);
  for(k=0;k<nvcube;k++) c2s->perms[k]=k;
  /* form all n+1 permutations in reverse order! why in reverse order?...*/
  nperm=1;
  for(j=c2s->n-1;j>=0;j--){
    jp=nperm*nvcube; jpo=jp+mid;
    //    fprintf(stdout,"\nnperm=%d,jp=%d,jpo=%d,face=%d; shift=%d",nperm,jp,jpo,c2s->n-j-1,shift);
    for(k = 0;k<nvcube;k++){
      kn=k*c2s->n;
      if((int)c2s->bits[kn+j]){
	c2s->perms[jp]=k;
	c2s->perms[jpo]=k-shift;
	jp++;jpo++;
      }
    }
    shift>>=1;
    nperm++;
  }
  /* fprintf(stdout,"\nNumber of permutations=%d",nperm); */
  /* for(j=0;j<nperm;j++){ */
  /*   jp=j*nvcube; */
  /*   fprintf(stdout,"\nPermutation=%d:",j+1); */
  /*   for(i=0;i<nvcube;i++){ */
  /*     fprintf(stdout," %d",c2s->perms[jp+i]+1); */
  /*   } */
  /* } */
  /* fprintf(stdout,"\n"); */
  return;
}
/***************************************************************************/
/*!
 * \fn static unsigned INT bitdiff(const INT dim, unsigned INT *bits1,\
 *                                 unsigned INT *bits2)
 *
 * \brief returns the l1-norm of the difference between two arrays of
 *        unsigned integers.  this should be changed to have a void
 *        array as input.
 *
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
static unsigned INT bitdiff(const INT dim, unsigned INT *bits1,unsigned INT *bits2)
{
  INT j;
  unsigned INT numbits=0;
  for(j=0;j<dim;j++){
    numbits+=abs(bits1[j]-bits2[j]);
  }
  return numbits;
}
/**********************************************************************/
/*!
 * \fn void reverse(void *arr,INT length, size_t elsize)
 *
 * \brief permutes a void array whose elements are of size elsize
 *        a[0],...a_[length-1]-->a[length-1],...a[0].
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void reverse(void *arr,INT length, size_t elsize)
{
  INT i,nnn=(INT)(length/2);
  void *swap=(void *)malloc(elsize);
  //  reverses ordering in an INT array;
  void *arrk=arr+elsize*(length-1);
  void *arri=arr;
  for(i=0;i<nnn;i++){
    memcpy(swap,arri,elsize);
    memcpy(arri,arrk,elsize);
    memcpy(arrk,swap,elsize);
    arri+=elsize;
    arrk-=elsize;
  }
  if(swap)free(swap);
  return;
}
/**********************************************************************/
/*!
 * \fn cube2simp *cube2simplex(INT dim)
 *
 * \brief in dimension dim splits the cube in dim factorial
 *        dim-dimensional simplices. stores everything in a structure
 *        cube2simp. It also outputs all local permutations of
 *        vertices which can be used to create a criss-cross mesh.
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
cube2simp *cube2simplex(INT dim)
{
  /* */
  INT i;
  /* allocation */
  cube2simp *c2s=malloc(sizeof(cube2simp));
  c2s->n=dim;
  c2s->ns=1;
  c2s->nvcube = (1 << c2s->n);
  c2s->nvface = (1 << (c2s->n-1));
  i=1; for(i=1;i<=c2s->n;i++)c2s->ns*=i;
  c2s->ne=c2s->n*(1<<(c2s->n-1)); /* number of edges in the n-cube.*/
  c2s->nf=2*c2s->n; /* number of n-1 dimensional faces in the n-cube */
  /////////////////////////////////////////////////////
  c2s->edges=(INT *)calloc(2*c2s->ne,sizeof(INT));
  c2s->bits=(unsigned INT *)calloc(c2s->n*c2s->nvcube,sizeof(unsigned INT));
  c2s->faces=(INT *)calloc(2*c2s->n*c2s->nvface,sizeof(INT));
  c2s->nodes=(INT *)calloc(c2s->ns*(c2s->n+1),sizeof(unsigned INT));
  c2s->perms=(INT *)calloc(c2s->nvcube*(c2s->n+1),sizeof(unsigned INT));
  memset(c2s->nodes,0,(c2s->n+1)*c2s->ns*sizeof(INT));
  /*end of allocation*/
  INT k1,kn1,k2,kn2,dim1=c2s->n+1,nvcube=c2s->nvcube;
  /***********************************************/
  binary0(c2s);
  /***********************************************/
  /****EDGES**********/
  INT *edges=c2s->edges;
  memset(edges,0,c2s->ne*sizeof(INT));
  unsigned INT numbits=22;
  unsigned INT *b1,*b2;
  //  fprintf(stdout,"Memory: edges=%d,ns=%d\n",ne,ns);
  INT nedge=0,nvcubem1=nvcube-1;
  for(k1 = 0;k1<nvcubem1;k1++){
    kn1=k1*dim;
    b1=c2s->bits+kn1;
    for(k2 = k1+1;k2<nvcube;k2++){
      kn2=k2*dim;
      b2=c2s->bits+kn2;
      numbits=bitdiff(dim,b1,b2)-1;
      if(numbits) continue;
      /* we found an edge, store it. */
      edges[nedge*2]=k1;
      edges[nedge*2+1]=k2;
      nedge++;
    }
  }
  /****SIMPLICES**********/
  INT root=0,j,m,node,nq0,nq;
  m=2; for(i=2;i<dim1;i++) m=1+i*m;
  INT *queue=(INT *)calloc(m,sizeof(INT));
  memset(queue,0,m*sizeof(INT));
  INT *parent=(INT *)calloc(m,sizeof(INT));
  memset(parent,0,m*sizeof(INT));
  // form a tree. every path in the tree is a simplex. the tree has
  // dim_factorial leaves.
  nq0=0;nq=1;parent[0]=-1;
  queue[nq0]=root;
  while(1){
    m=nq;
    for(j=nq0;j<nq;j++){
      node=queue[j];
      if(node==(nvcubem1)) continue;
      for(i=0;i<c2s->ne;i++){
	/*
	  for a given first end of an edge, collect all the second
	  ends in the queue;
	*/
	if(edges[2*i]==node){
	  queue[m]=edges[2*i+1];
	  parent[m]=j;
	  //	  fprintf(stdout,"(m=%d;%d)",m,edges[2*i+1]);
	  m++;
	}
      }
    }
    if(nq>=m) break;
    nq0=nq;
    nq=m;
  }
  //  fprintf(stdout,"\nlast:=%d\n",nq-nq0);
  k1=0;// simplex number;
  for(j=nq0;j<nq;j++){
    i=c2s->n;
    node=queue[j];
    m=j;
    while(parent[m]>=0){
      c2s->nodes[k1*dim1+i]=queue[m];
      m=parent[m];
      i--;
    }
    k1++;// increment simplex number;
  }
  if(queue)free(queue);
  if(parent)free(parent);
  // finally reverse all bits to have the last coordinate ordered first
  // in the local numbering.
  for(j=0;j<c2s->nvcube;j++){
    reverse((c2s->bits+dim*j),dim,sizeof(INT));
    /* fprintf(stdout,"\nj:%d ",j); */
    /* for(i=0;i<dim;i++){ */
    /*   fprintf(stdout,"%d",c2s->bits[dim*j+i]); */
    /* } */
  }
  /****FACES**********/
  INT *faces=c2s->faces;
  INT j0,j1;
  memset(faces,0,c2s->nf*sizeof(INT));
  //  fprintf(stdout,"\n");
  //    fprintf(stdout,"\nk2=%d bits=(",k2);
  //	fprintf(stdout,"%d,",k1);
  //    fprintf(stdout,")");
  //  fprintf(stdout,"\n");fflush(stdout);
  for(k2 = 0;k2<c2s->n;k2++){
    kn2=k2*c2s->nvface;
    j0=0;j1=0;
    for(k1 = 0;k1<nvcube;k1++){
      kn1=k1*c2s->n;
      if(!c2s->bits[kn1+k2]){
	c2s->faces[kn2+j0]=k1;
	j0++;
      } else {
	c2s->faces[kn2+c2s->n*c2s->nvface+j1]=k1;
	j1++;
      }
    }
  }
  //  print_full_mat_int(c2s->nf,c2s->nvface,c2s->faces,"UCubef");fflush(stdout);
  return c2s;
}
/******************************************************************/
/*!
 * \fn INT dvec_set_amr(const REAL value, scomplex *sc, INT npts, REAL *pts,
 * dvector *toset)
 *
 * \brief given a dvector of size sc->ns sets the values of a vector
 *        to equal value at every simplex that is on the last level of
 *        refinement (not refined simplex) and contains a point from
 *        dvector pts (note that the size of pts->val should be
 *        sc->n*pts->row)
 *
 * \param dvector toset;
 *
 * \return number of simplices where the value was assigned
 *
 */
INT dvec_set_amr(const REAL value, scomplex *sc, INT npts, REAL *pts, REAL *toset)
{
  INT j,k,jpts,n=sc->n,n1=sc->n+1,ns=sc->ns;
  REAL *pval0=NULL; /*place holder*/
  INT *scnjn=NULL; /*place holder*/
  k=0;
  for(j=0;j<ns;j++) {
    scnjn = sc->nodes+j*n1; /* beginning of local vertex numbering for
			       simplex j.*/
    for(jpts=0;jpts<npts;jpts++){
      pval0=pts + jpts*n;
      if(!xins(n,scnjn,sc->x,pval0)){
	//	fprintf(stdout,"\nel=%d, found: %d",j,jpts);
	toset[j]=value;
	k++;
	break;
      }
    }
  }
  return k;
}
/*COMPUTING CONNECTED COMPONENTS*/
/**********************************************************************/
/*!
 * \fn void find_cc_bndry_cc(scomplex *sc, const INT set_bndry_codes)
 *
 * \brief in a simplicial complex, finds all connected components and
 *        the connected components on the boundary.  the arrays
 *        nodes[] and nbr[] should be set. The neighbor of simplex i
 *        sharing the face formed by nodes[1,...,j-1,j+1,...] should
 *        be at place j in the neighbors of i.
 *
 * \param sc: a simplicial complex; sc->bndry, sc->neib, sc->nbr must
 *            be allocated and filled in on entry here.
 *
 * \param set_bndry_codes : if true, then sets all boundary codes on
 *                          every connected component to be different.
 *
 *
 * \note
 *
 */
void find_cc_bndry_cc(scomplex *sc,INT set_bndry_codes)
{
  //
  INT ns = sc->ns, dim=sc->n;
  INT dim1=dim+1,iii,i,j,k,m,isn1,is,nbf,nnzbf;
  iCSRmat s2s=icsr_create(ns,ns,dim1*ns+ns);
  nbf=0;
  nnzbf=0;
  iii=0;
  s2s.IA[0]=iii;
  for(i=0;i<ns;i++){
    isn1=i*dim1;
    for(j=0;j<dim1;j++){
      is=sc->nbr[isn1+j];
      if(is>=0){
	s2s.JA[iii]=is;
	s2s.val[iii]=1;
	iii++;
      } else {
	nbf++;
	nnzbf+=dim;
      }
    }
    s2s.IA[i+1]=iii;
  }
  //  fprintf(stdout,"\n");fflush(stdout);
  sc->cc=-10;
  // find the connected components in the domain:
  iCSRmat *blk_dfs=run_dfs(ns,s2s.IA, s2s.JA);
  sc->cc=blk_dfs->row;
  for(i=0;i<sc->cc;++i){
    for(k=blk_dfs->IA[i];k<blk_dfs->IA[i+1];++k){
      j=blk_dfs->JA[k];
      sc->flags[j]=i+1;
    }
  }
  // all this was to comopute the connected components. Let us free the memory now:
  icsr_free(&s2s);// no need of this anymore.
  //
  // now working on the boundary:
  //
  iCSRmat f2v=icsr_create(nbf,sc->nv,nbf*dim);
  INT nbfnew=0;
  INT nnzf2v=0;
  f2v.IA[0]=nnzf2v;
  for(i=0;i<ns;i++){
    for(j=0;j<dim1;j++){
      if(sc->nbr[i*dim1+j]<0) {
	for(m=0;m<dim1;m++){
	  if(m==j) continue;
	  f2v.JA[nnzf2v]=sc->nodes[i*dim1+m];
	  f2v.val[nnzf2v]=1;
	  nnzf2v++;
	}
	nbfnew++;
	f2v.IA[nbfnew]=nnzf2v;
      }
    }
  }
  f2v.nnz=nnzf2v;
  if(nbf!=nbfnew){
    fprintf(stderr,"\n%%***ERROR(1): num. bndry faces mismatch (nbf=%d .ne. nbfnew=%d) in %s",nbf,nbfnew,__FUNCTION__);
    exit(65);
  }
  // FIX numbering (ignoring all interior vertices):
  INT *indx    = calloc(sc->nv,sizeof(INT));
  INT *indxinv = calloc(sc->nv,sizeof(INT));
  for(i=0;i<sc->nv;++i) indx[i]=-1;
  for(i=0;i<nnzf2v;++i) indx[f2v.JA[i]]++;
  f2v.col=0;
  for(i=0;i<sc->nv;++i){
    if(indx[i]<0) continue;
    indx[i]=f2v.col;
    indxinv[f2v.col]=i;
    f2v.col++;
  }
  if(f2v.col<sc->nv)
    indxinv=realloc(indxinv,f2v.col*sizeof(INT));
  for(i=0;i<f2v.nnz;++i)
    f2v.JA[i]=indx[f2v.JA[i]];
  // end fix numbering. f2v is constructed
  iCSRmat v2f,f2f;
  INT j0,j1,ke,je,found;
  icsr_trans(&f2v,&v2f);
  /*******************************************************************/    
  icsr_mxm(&f2v,&v2f,&f2f);
  /*******************************************************************/    
  icsr_free(&v2f);// we do not need v2f now
  /*******************************************************************/    
  /* 
     now remove all rows in f2f that correspond to interior faces and
     all entries that are not dim, i.e. the number of vertices in
     a (n-2)-simplex;
  */
  /* fprintf(stdout,"\nnrf2f=%d ; ncf2f=%d;\n f2f=[",f2f.row,f2f.col); */
  /* icsr_print_matlab_val(stdout,&f2f); */
  /* fprintf(stdout,"];\n"); */
  /* fprintf(stdout,"f2f=sparse(f2f(:,1),f2f(:,2),f2f(:,3),nrf2f,ncf2f);\n\n"); */
  f2f.nnz=f2f.IA[0];
  for(i=0;i<f2f.row;i++){    
    j0=f2f.IA[i];
    j1=f2f.IA[i+1];
    f2f.IA[i]=f2f.nnz;
    for(ke=j0;ke<j1;ke++){
      je=f2f.JA[ke];
      if(je==i){
	f2f.JA[f2f.nnz]=i;
	f2f.val[f2f.nnz]=1;
	f2f.nnz++;
      	continue;
      }
      if(f2f.val[ke]!=(dim-1)) continue;
      //      if((je==i)||(isbface[je]==0)||(f2f.val[ke]!=(nvface/2))) continue;
      f2f.JA[f2f.nnz]=je;
      f2f.val[f2f.nnz]=f2f.val[ke];
      f2f.nnz++;
    }
  }
  //icsr_nodiag(f2f);
  f2f.IA[f2f.row]=f2f.nnz;
  f2f.JA=realloc(f2f.JA,f2f.nnz*sizeof(INT));
  f2f.val=realloc(f2f.val,f2f.nnz*sizeof(INT));
  /* fprintf(stdout,"\nnrf2f1=%d ; ncf2f1=%d;\n f2f1=[",f2f.row,f2f.col); */
  /* icsr_print_matlab_val(stdout,&f2f); */
  /* fprintf(stdout,"];\n"); */
  /* fprintf(stdout,"f2f1=sparse(f2f1(:,1),f2f1(:,2),f2f1(:,3),nrf2f1,ncf2f1);\n\n"); */
  /*******************************************************************/    
  icsr_free(blk_dfs);free(blk_dfs);
  blk_dfs=run_dfs(f2f.row,f2f.IA, f2f.JA);
  sc->bndry_cc=0;
  for(i=0;i<blk_dfs->row;++i){
    found=blk_dfs->IA[i+1]-blk_dfs->IA[i];
    if(found>1) sc->bndry_cc++;
  }
  icsr_free(&f2f);
  /*******************************************************************/    
  if(set_bndry_codes) {
    for(i=0;i<sc->bndry_cc;++i){
      for(k=blk_dfs->IA[i];k<blk_dfs->IA[i+1];++k){
	j=blk_dfs->JA[k];
	for(m=0;m<dim;m++){
	  sc->bndry[indxinv[f2v.JA[dim*j+m]]]=i+1+128;
	}
      }
    }
  } else {
    // here we find the boundary codes of every added point.
    // now do nothing
  }
  icsr_free(blk_dfs);free(blk_dfs);
  icsr_free(&f2v);
  fprintf(stdout,"%%%%--> number of connected components in the bulk=%d\n",sc->cc);
  //  fprintf(stdout,"%%number of boundary faces=%d (nnzbf=%d)\n",nbf,nnzbf);
  fprintf(stdout,"%%%%--> number of connected components on the boundary=%d\n",sc->bndry_cc);
  free(indx);
  free(indxinv);
  return;
}
/*EOF*/
