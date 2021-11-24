/*! \file src/amr/scomplex.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note containing all essentials routines for mesh refinement
 *
 */
#include "hazmath.h"
/**********************************************************************/
/*!
 * \fn REAL void haz_scomplex_realloc(scomplex *sc)
 *
 * \brief  reallocates memory if not allocated with haz_scomplex_init;
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void haz_scomplex_realloc(scomplex *sc)
{
  INT i,j,ns=sc->ns,nv=sc->nv,n=sc->n;
  INT n1=n+1;
  sc->print_level=0;
  sc->factorial=1.;
  for (j=2;j<n1;j++) sc->factorial *= ((REAL )j);
  sc->nbr=realloc(sc->nbr,n1*ns*sizeof(INT));
  sc->marked=realloc(sc->marked,ns*sizeof(INT));
  sc->gen=realloc(sc->gen,ns*sizeof(INT));
  sc->parent=realloc(sc->parent,ns*sizeof(INT));
  sc->child0=realloc(sc->child0,ns*sizeof(INT));
  sc->childn=realloc(sc->childn,ns*sizeof(INT));
  sc->bndry=realloc(sc->bndry,nv*sizeof(INT));
  sc->csys=realloc(sc->csys,nv*sizeof(INT));
  sc->flags=realloc(sc->flags,ns*sizeof(INT)); // element flags
  sc->vols=realloc(sc->vols,ns*sizeof(REAL)); // element volumes
  for (i = 0;i<sc->ns;i++) {
    sc->marked[i] = FALSE; // because first this array is used as working array.
    sc->gen[i] = 0;
    sc->parent[i]=-1;
    sc->child0[i]=-1;
    sc->childn[i]=-1;
    sc->flags[i]=-1;
  }
  for (i = 0;i<nv;i++) {
    sc->bndry[i]=0;
    sc->csys[i]=0;
  }
  return;
}
/**********************************************************************/
/*!
 * \fn REAL chk_sign(const int it, const int nbrit)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
REAL chk_sign(const int it, const int nbrit)
{
/*
  nbrit is a neighbor of it and this picks the sign of the normal
  vector for the face shared by both. Used to construct dcsrmat atf
  and local matrices to build aff. If the face is on the boundary,
  i.e. nbrit <0, we always return 1 (i.e. the outward normal).
*/
  if(it>nbrit) return 1e0;
  return -1e0;
}
/*********************************************************************/
/*!
 * \fn REAL volume_compute(INT dim, REAL factorial, REAL *xs,void *wrk);
 *
 * \brief Computes the volume of the simplex in Rn, which is det(B)/d_factorial.
 *
 * \param dim        I: The dimension of the problem.
 * \param factorial  I: dim! (dim factorial)
 * \param xs         I: coordinates of the simplex
 * \param wrk        W: working array of dimension
 *                      (dim+1)*(dim*sizeof(REAL) + sizeof(INT))
 *
 * \return the volume of the simplex.
 *
 * \note
 */
REAL volume_compute(INT dim, REAL factorial, REAL *xs,void *wrk)
{
  INT dim1 = dim+1,i,j,ln,ln1;
  REAL *bt=(REAL *)wrk;
  REAL *piv=bt+dim*dim;
  INT *p = (INT *)(wrk+(dim*dim + dim)*sizeof(REAL));
  REAL vol=-1e20;
  // construct bt using xs;
  for (j = 1;j<dim1;j++){
    ln=j*dim; ln1=ln-dim;
    for(i=0;i<dim;i++){
      bt[ln1+i] = xs[ln+i]-xs[i];  // k-th row of bt is [x(k)-x(0)]^T. x(k) are coords of vertex k.
    }
  }
  //  SHORT flag=ddense_lu(1, dim, &vol, bt,p,piv);
  //  if(flag){
  if(ddense_lu(1, dim, &vol, bt,p,piv)) {
    return 0e0; // degeneraate simplex;
  } else
    return fabs(vol)/factorial;
}
/*********************************************************************/
/*!
 * \fn void sc_vols(scomplex *sc);
 *
 * \brief Fills in the array with volumes of the simpleices in
 *        n-dimansional simplicial complex.
 *
 * \param sc: pointer to the simplicial complex sc.
 *
 * \note
 */
void sc_vols(scomplex *sc)
{
  INT dim = sc->n, ns=sc->ns;
  INT dim1 = dim+1,i,j,node,idim1;
  void *wrk=malloc(dim1*(dim*sizeof(REAL) + sizeof(INT)));
  REAL *xs=calloc(dim1*dim,sizeof(REAL));
  for(i=0;i<ns;++i){
    idim1=i*dim1;
    for (j = 0;j<dim1;++j){
      node=sc->nodes[idim1+j];
      memcpy((xs+j*dim),(sc->x+node*dim),dim*sizeof(REAL));
    }
    sc->vols[i]=volume_compute(dim,sc->factorial,xs,wrk);
  }
  free(wrk);
  free(xs);
  return;
}
/**********************************************************************/
/*!
 * \fn scomplex *haz_scomplex_init(const INT n,INT ns, INT nv, const INT nbig)
 *
 * \brief Initialize simplicial complex in dimension n with ns
 *        simplices and nv vertices.
 *
 * \param n is the dimension of the simplicial complex;
 * \param ns is the number of simplices
 * \param nv is the number of vertices
 * \param nbig is the dimension of the space where the simplicial
 *        complex is embedded;
 *
 * \return initialized structure of type scomplex
 *
 * \note
 *
 */
scomplex *haz_scomplex_init(const INT n,INT ns, INT nv,const INT nbig)
{
  /*
     n = dimension of the simplicial complex;
     nbig = dimension of the space where the complex is embedded;
     Future work: we should think of making this for different
     dimensions, e.g. 2-homogenous complex in 3d
  */
  scomplex *sc=(scomplex *) malloc(sizeof(scomplex));
  sc->nbig=nbig; sc->n=n;
  if(sc->nbig<=0)sc->nbig=n;
  INT n1=sc->n+1,i,j,in1;
  sc->level=0;
  sc->ref_type=0;
  sc->print_level=0;
  sc->factorial=1.;
  for (j=2;j<n1;j++) sc->factorial *= ((REAL )j);
  //
  sc->marked=(INT *) calloc(ns,sizeof(INT));
  sc->gen=(INT *) calloc(ns,sizeof(INT));
  sc->nbr=(INT *) calloc(ns*n1,sizeof(INT));
  sc->parent=(INT *)calloc(ns,sizeof(INT));
  sc->child0=(INT *)calloc(ns,sizeof(INT));
  sc->childn=(INT *)calloc(ns,sizeof(INT));
  sc->nodes=(INT *)calloc(ns*n1,sizeof(INT));
  sc->bndry=(INT *)calloc(nv,sizeof(INT));
  sc->csys=(INT *)calloc(nv,sizeof(INT));/* coord sys: 1 is polar, 2
					    is cyl and so on */
  sc->bndry_v=malloc(sizeof(iCSRmat));
  sc->bndry_v[0]=icsr_create(0,0,0); // only alloc pointers; two are null. 
  sc->parent_v=malloc(sizeof(iCSRmat));
  sc->parent_v[0]=icsr_create(nv,nv,nv);
  sc->flags=(INT *)calloc(ns,sizeof(INT));
  sc->x=(REAL *)calloc(nv*nbig,sizeof(REAL));
  sc->vols=(REAL *)calloc(ns,sizeof(REAL));
  for (i = 0;i<ns;i++) {
    sc->marked[i] = FALSE; // because first is used for something else.
    sc->gen[i] = 0;
    sc->parent[i]=-1;
    sc->child0[i]=-1;
    sc->childn[i]=-1;
    sc->flags[i]=-1;
    sc->vols[i]=-1e20;
    in1=i*n1;
    for(j=0;j<n1;j++){
      sc->nodes[in1+j]=-1;
      sc->nbr[in1+j]=-1;
    }
  }
  INT nnz_pv=0;
  sc->parent_v->IA[0]=nnz_pv;
  for (i = 0;i<nv;i++) {
    sc->bndry[i]=0;
    sc->csys[i]=0;
    sc->parent_v->JA[nnz_pv]=i;
    nnz_pv++;
    sc->parent_v->IA[i+1]=nnz_pv;
  }
  // not needed for now, it will be freed later.
  if(nnz_pv) memset(sc->parent_v->val,0,nnz_pv*sizeof(INT));
  //////////////////////////////////////
  sc->nv=nv;
  sc->ns=ns;
  sc->bndry_cc=1; // one connected component on the boundary for now.
  sc->cc=1; // one connected component in the bulk for now.
  // NULL pointers for the rest
  sc->etree=NULL;
  sc->bfs=malloc(sizeof(iCSRmat));
  sc->bfs[0]=icsr_create(0,0,0);
  // the parent_v->val is not needed for now
  if(sc->parent_v->val) {
    free(sc->parent_v->val);
    sc->parent_v->val=NULL;
  }
  sc->bndry_v=malloc(sizeof(iCSRmat));
  sc->bndry_v[0]=icsr_create(0,0,0);
  return sc;
}
/**********************************************************************/
/*!
 * \fn scomplex haz_scomplex_null(const INT n,INT ns, INT nv, const INT nbig)
 *
 * \brief Initialize all pointers in simplicial complex to NULL;
 * \param n is the dimension of the simplicial complex;
 * \param nbig is the dimension of the space where the simplicial
 *        complex is embedded;
 *
 * \return initialized structure of type scomplex
 *
 * \note
 *
 */
scomplex haz_scomplex_null(const INT n,const INT nbig)
{
  /*
     n = dimension of the simplicial complex;
     nbig = dimension of the space where the complex is embedded;
     Future work: we should think of making this for different
     dimensions, e.g. 2-homogenous complex in 3d
  */
  scomplex sc;
  sc.nbig=nbig; sc.n=n;
  if(sc.nbig<=0)sc.nbig=n;
  INT n1=sc.n+1,j;
  sc.level=0;
  sc.ref_type=0;
  sc.print_level=0;
  sc.factorial=1.;
  for (j=2;j<n1;j++) sc.factorial *= ((REAL )j);
  //
  sc.marked=NULL;
  sc.gen=NULL;
  sc.nbr=NULL;
  sc.parent=NULL;
  sc.child0=NULL;
  sc.childn=NULL;
  sc.nodes=NULL;
  sc.bndry=NULL;
  sc.csys=NULL;
  sc.bndry_v=malloc(sizeof(iCSRmat));
  sc.bndry_v[0]=icsr_create(0,0,0); //
  sc.parent_v=malloc(sizeof(iCSRmat));
  sc.parent_v[0]=icsr_create(0,0,0);
  sc.flags=NULL;
  sc.x=NULL;
  sc.vols=NULL;
  sc.bndry_cc=1; // one connected component on the boundary for now.
  sc.cc=1; // one connected component in the bulk for now.
  // NULL pointers for the rest
  sc.etree=NULL;
  sc.bfs=malloc(sizeof(iCSRmat));
  sc.bfs[0]=icsr_create(0,0,0);
  return sc;
}
/**********************************************************************/
/*!
 * \fn scomplex *haz_scomplex_init_part(scomplex *sc)
 *
 * \brief Initialize simplicial complex which has already been
 *        allocated partially.
 * Assumption is: nv,ns,nodes,nbr, cc,
 *        bndry_cc are known. in dimension n with ns simplices and nv
 *        vertices.
 *
 * \param sc pointer to partially allocated simplician complex;
 *
 * \return full structure of type scomplex with all allocations done
 *
 * \note
 *
 */
void haz_scomplex_init_part(scomplex *sc)
{
  INT nv=sc->nv,ns=sc->ns,n1=sc->n+1,i,j;
  sc->level=0;
  sc->ref_type=0;
  sc->print_level=0;
  sc->factorial=1.;
  for (j=2;j<n1;j++) sc->factorial *= ((REAL )j);
  //
  sc->marked=(INT *) calloc(ns,sizeof(INT));
  sc->gen=(INT *) calloc(ns,sizeof(INT));
  sc->parent=(INT *)calloc(ns,sizeof(INT));
  sc->child0=(INT *)calloc(ns,sizeof(INT));
  sc->childn=(INT *)calloc(ns,sizeof(INT));
  sc->bndry=(INT *)calloc(nv,sizeof(INT));
  sc->csys=(INT *)calloc(nv,sizeof(INT));/* coord sys: 1 is polar, 2
					    is cyl and so on */
  sc->flags=(INT *)calloc(ns,sizeof(INT)); // element flags
  //  sc->vols=(REAL *)calloc(ns,sizeof(REAL));// simplex volumes
  for (i = 0;i<ns;i++) {
    sc->marked[i] = FALSE; // because first is used for something else.
    sc->gen[i] = 0;
    sc->parent[i]=-1;
    sc->child0[i]=-1;
    sc->childn[i]=-1;
    sc->flags[i]=-1;
  }
  for (i = 0;i<nv;i++) {
    sc->bndry[i]=0;
    sc->csys[i]=0;
  }
  return;
}
/**********************************************************************/
/*!
 * \fn void vol_simplex(INT dim, REAL fact, REAL *xf, REAL *volt, void *wrk)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void vol_simplex(INT dim, REAL fact, REAL *xf, REAL *volt, void *wrk)
{
  /*
     computes the volume of a simplex; wrk should be at least
     dim*(dim+1) REALS and dim integers.
  */
  INT dim1 = dim+1,i,j,ln,ln1;
  REAL *bt=(REAL *)wrk;
  REAL *piv=bt+dim*dim;
  INT *p = (INT *)(wrk+(dim*dim + dim)*sizeof(REAL));
  // construct bt using xf;
  for (j = 1;j<dim1;j++){
    ln=j*dim; ln1=ln-dim;
    for(i=0;i<dim;i++){
      bt[ln1+i] = xf[ln+i]-xf[i];
    }
  }
  //  print_full_mat(dim,dim,bt,"bt");
  if(ddense_lu(1, dim, volt, bt,p,piv))
    *volt=0e0;
  else
    *volt=fabs(*volt)/fact;
  return;
}
/**********************************************************************/
/*!
 * \fn scomplex *haz_scomplex_read(FILE *fp)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
scomplex *haz_scomplex_read(FILE *fp,INT print_level)
{
  INT i,ns,nv,n,dummy;
  i=fscanf(fp,"%d %d %d %d",&ns,&nv,&n,&dummy);
  INT n1=n+1,j,k,n1kj=-10,nbig=n;// we can only read same dimension complexes now.
  scomplex *sc=(scomplex *)haz_scomplex_init(n,ns,nv,n);
  for (j=0;j<n1;j++) {
    for (k=0;k<ns;k++){
      n1kj=n1*k+j;
      dummy=fscanf(fp," %d ", sc->nodes+n1kj);
      /* shift if needed ; this should not be here: later CHANGE */
      sc->nodes[n1kj]=sc->nodes[n1kj]-1;
    }
  }
  for (k=0;k<ns;k++){
    dummy=fscanf(fp," %d ", sc->flags+k);
  }
  for(j=0;j<nbig;j++){
    for(i=0;i<nv;i++){
      dummy=fscanf(fp,"%lg",sc->x+i*nbig+j);
    }
  }
  for(i=0;i<nv;i++){
    dummy=fscanf(fp,"%i",sc->bndry+i);
  }
  sc->print_level=0;
  return sc;
}
/**********************************************************************/
/*!
 * \fn void haz_scomplex_print(scomplex *sc, const INT ns0,const char *infor)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void haz_scomplex_print(scomplex *sc, const INT ns0,const char *infor)
{
  // print simplicial complex, starting with ns0.
  INT i,j,in,in1;
  INT n=sc->n,n1=n+1,ns=sc->ns,nv=sc->nv,nbig=sc->nbig;
  if (ns0 < 0 || ns0>ns) return;
  fprintf(stdout,"\nN=%d,NBIG=%d, NS=%d, NV=%d\n",sc->n,sc->nbig,sc->ns,sc->nv);fflush(stdout);
  fprintf(stdout,"\n%s printout: %s\n",__FUNCTION__,infor);
  fprintf(stdout,"\nNODES list:\n");
  if(sc->parent){
    for(i=ns0;i<ns;i++){
      in1=i*n1;
      fprintf(stdout,"Element: %d ; vol=%e, Parent=%d; NODES=",i-ns0,sc->vols[i],sc->parent[i-ns0]);
      for(j=0;j<n1;j++)
	fprintf(stdout,"%d  ",sc->nodes[in1+j]);
      fprintf(stdout,"\n");
    }
  } else {
    for(i=ns0;i<ns;i++){
      in1=i*n1;
      fprintf(stdout,"Element: %d ; vol=%e, NODES=",i-ns0,sc->vols[i]);
      for(j=0;j<n1;j++)
	fprintf(stdout,"%d  ",sc->nodes[in1+j]);
      fprintf(stdout,"\n");
    }
  }
  fprintf(stdout,"\nNBR list:\n");
  if(sc->gen){
    for(i=ns0;i<ns;i++){
      in1=i*n1;
      fprintf(stdout,"Element: %d (%d) ; NBR=",i-ns0,sc->gen[i-ns0]);
      for(j=0;j<n1;j++)
	fprintf(stdout,"%d  ",sc->nbr[in1+j]-ns0);
      fprintf(stdout,"\n");
    }
  } else {
    for(i=ns0;i<ns;i++){
      in1=i*n1;
      fprintf(stdout,"Element: %d ; NBR=",i-ns0);
      for(j=0;j<n1;j++)
	fprintf(stdout,"%d  ",sc->nbr[in1+j]-ns0);
      fprintf(stdout,"\n");
    }
  }
  //
  if(sc->bndry){
    for(i=0;i<nv;i++){
      in=i*nbig;
      fprintf(stdout,"Node: %d ; Code: %d ; COORDS=",i,sc->bndry[i]);
      for(j=0;j<nbig;j++){
	fprintf(stdout,"%e  ",sc->x[in+j]);
      }
      fprintf(stdout,"\n");
    }
  } else {
    for(i=0;i<nv;i++){
      in=i*nbig;
      fprintf(stdout,"Node: %d ; COORDS=",i);
      for(j=0;j<nbig;j++){
	fprintf(stdout,"%e  ",sc->x[in+j]);
      }
      fprintf(stdout,"\n");
    }
  }
  return;
}
/**********************************************************************/
/*!
 * \fn void haz_scomplex_free(scomplex *sc)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void haz_scomplex_free(scomplex *sc)
{
  if(sc->marked) free(sc->marked);
  if(sc->gen) free(sc->gen);
  if(sc->parent) free(sc->parent);
  if(sc->child0) free(sc->child0);
  if(sc->childn) free(sc->childn);
  if(sc->bndry) free(sc->bndry);
  if(sc->flags) free(sc->flags);
  if(sc->nodes) free(sc->nodes);
  if(sc->x) free(sc->x);
  if(sc->vols) free(sc->vols);
  if(sc->nbr) free(sc->nbr);
  if(sc->csys) free(sc->csys);
  if(sc->etree) free(sc->etree);
  if(sc->bndry_v) {
    icsr_free(sc->bndry_v);free(sc->bndry_v);sc->bndry_v=NULL;
  }
  if(sc->parent_v) {
    icsr_free(sc->parent_v);free(sc->parent_v);sc->parent_v=NULL;
  }
  if(sc->bfs) {
    icsr_free(sc->bfs);free(sc->bfs);sc->bfs=NULL;
  }
  if(sc) free(sc);
  return;
}
/**********************************************************************/
/*!
 * \fn static unsigned int cmp_simplex(INT n, INT sim1, INT sim2, INT
 *			 *sv1, INT *sv2, INT *stos1, INT *stos2)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
static unsigned int cmp_simplex(INT n, INT sim1, INT sim2,	\
			 INT *sv1, INT *sv2, INT *stos1, INT *stos2)
{
  //sim1 and sim2 are pointers to the neighboring lists of two
  //simplices stos1 and stos2 are pointers to rows in
  //simplex-to-simplex incidence; sv1 and sv2 are pointers to rows in
  //simplex-vertex incidence for two neighboring simplices.
  //
  //rearrange to meet the structural condition.  compare two sets of
  //length n. they should overlap at exactly n-1 members.  If they do
  //not, 0 is returned, otherwise 1. no check if the members of the
  //sets are distinct (they should be for the purpose of comparing two
  //simplices.
  unsigned int fnd;
  unsigned int nf=0;
  INT i1,k1=-10,i2,k2=-10;
  INT i,j;
  for (i=0;i<n;i++){
    fnd=0;
    i1=sv1[i];
    for (j = 0;j < n;j++){
      if(i1 == sv2[j]){
	fnd=1;
	break;
      }
    }
    if(fnd) {
      continue;
    } else {
      // not found
      nf++;
      if(nf > 1) return 0;
      k1=i;
    }
  }
  // same with sv1 and sv2 swapped.
  nf=0;
  for (i=0;i<n;i++){
    fnd=0;
    i2=sv2[i];
    for (j=0;j<n;j++){
      if(i2 == sv1[j]){
	fnd=1;
	break;
      }
    }
    if(fnd) {
      continue;
    } else {
      // not found
      nf++;
      if(nf > 1) return 0;
      k2=i;
    }
  }
  /* NOW put the neightbors at the right places ******* */
  if(k1<0||k2<0){
    fprintf(stderr,"\n***ERROR in %s; k1=%d,k2=%d\n\n",__FUNCTION__,k1,k2);
    exit(65);
  } else {
    stos1[k1]=sim2;
    stos2[k2]=sim1;
    return 1;
  }
}
/**********************************************************************/
/*!
 * \fn void find_nbr(INT ns,INT nv,INT n,INT *sv,INT *stos)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void find_nbr(INT ns,INT nv,INT n,INT *sv,INT *stos)
{
  // find neighboring list
  INT* ivs=NULL, *jvs=NULL, *stosi=NULL, *stosk=NULL, *svi=NULL, *svk=NULL;
  INT i,j,k,jp,jia,jib,iabeg,iaend,ibbeg,ibend,kn1;
  INT n1 = n+1,nv1=nv+1,nsv=ns*n1;
  /* allocate */
  ivs=(INT *) calloc(nv1,sizeof(INT));
  jvs=(INT *) calloc(nsv,sizeof(INT));
  /* transpose to get the vertex--simplex incidence */
  for (i = 0; i < nv1; ++i)
    ivs[i] = 0;
  //
  for (k = 0; k < ns; ++k) {
    kn1=k*n1;
    for (i = 0; i < n1; i++) {
      j = sv[kn1+i] + 2;
      if (j < nv1) ivs[j]++;
    }
  }
  ivs[0] = 0;
  ivs[1] = 0;
  if (nv != 1) {
    for (i= 2; i < nv1; ++i) {
      ivs[i] += ivs[i-1];
    }
  }
  for (i = 0; i < ns; ++i) {
    iabeg = i*n1;
    iaend = i*n1+n1;
    for (jp = iabeg; jp < iaend; ++jp) {
      j = sv[jp]+1;
      k = ivs[j];
      jvs[k] = i;
      ivs[j] = k+1;
    }
  }
  /**/
  INT *icp=(INT *) calloc(ns,sizeof(INT));
  for (i = 0; i < ns; ++i) icp[i] = -1;
  for (i = 0; i < nsv; ++i) stos[i] = -1;
  for (i = 0; i < ns; ++i) {
      iabeg = i*n1;
      iaend = iabeg + n1;
      stosi=stos+iabeg; svi=sv+iabeg;
      for (jia = iabeg; jia < iaend; ++jia) {
	j = sv[jia]; // vertex number
	ibbeg = ivs[j];
	ibend = ivs[j+1];
	// loop over all simplices with this j as a vertex.
	for (jib = ibbeg; jib < ibend; ++jib) {
	  k = jvs[jib];
	  if(k<=i) continue;
	  if (icp[k] != i) {
	    icp[k] = i;
	    kn1=k*n1;
	    stosk=stos+kn1; svk=sv+kn1;
	    if(!cmp_simplex(n1,i,k,svi,svk,stosi,stosk)) continue;
	  }//if
	} //for
      } //for
    } //for (i...
    if (icp) free(icp);
    if (ivs) free(ivs);
    if (jvs) free(jvs);
    //  }
    return;
}
/**********************************************************************/
/*!
 * \fn INT haz_add_simplex(INT is, scomplex *sc,REAL *xnew, INT *pv,
 *                         INT ibnew,INT csysnew,INT nsnew, INT nvnew)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
INT haz_add_simplex(INT is, scomplex *sc,REAL *xnew,	\
		    INT *pv,INT ibnew,INT csysnew,	\
		    INT nsnew, INT nvnew)
{
  /* adds nodes and coords as well */
  INT n=sc->n, nbig=sc->nbig, n1 = n+1,nv=sc->nv;//ns=sc->ns;
  INT ks0=sc->child0[is], ksn=sc->childn[is];
  INT isc0=ks0*n1, iscn=ksn*n1;
  //  INT *dsti,*srci;
  INT j,j0,jn,nnz_pv;
  REAL *dstr;
  /* nodes  AND neighbors */
  sc->nbr=realloc(sc->nbr,(nsnew*n1)*sizeof(INT));
  sc->nodes=realloc(sc->nodes,(nsnew*n1)*sizeof(INT));
  for(j=0;j<n1;j++){
    j0=isc0+j;
    jn=iscn+j;
    //
    sc->nodes[j0]=-1;
    sc->nbr[j0]=-1;
    sc->nodes[jn]=-1;
    sc->nbr[jn]=-1;
  }
  //new vertex (if any!!!)
  if(nvnew != nv) {
    sc->x=realloc(sc->x,(nvnew*nbig)*sizeof(REAL));
    dstr=(sc->x+nv*nbig); memcpy(dstr,xnew,n*sizeof(REAL));
    if(xnew) free(xnew);
    sc->bndry=realloc(sc->bndry,(nvnew)*sizeof(INT));
    sc->bndry[nv]=ibnew;
    sc->csys=realloc(sc->csys,(nvnew)*sizeof(INT));
    sc->csys[nv]=csysnew;
    nnz_pv=sc->parent_v->nnz;
    sc->parent_v->row=nvnew;
    sc->parent_v->col=nv;
    sc->parent_v->JA=realloc(sc->parent_v->JA,(nnz_pv+2)*sizeof(REAL));
    sc->parent_v->JA[nnz_pv]=pv[0];
    sc->parent_v->JA[nnz_pv+1]=pv[1];
    sc->parent_v->nnz+=2;
    sc->parent_v->IA=realloc(sc->parent_v->IA,(nvnew+1)*sizeof(REAL));
    sc->parent_v->IA[nvnew]=sc->parent_v->nnz;
    /* fprintf(stdout,"\nnv=%d; nvnew=%d;nnz_pv=%d(pv[0]=%d,pv[1]=%d)",nv,nvnew,sc->parent_v->nnz,pv[0],pv[1]); */
  }
  //generation
  sc->gen=realloc(sc->gen,(nsnew)*sizeof(INT));
  sc->gen[ks0]=sc->gen[is]+1;
  sc->gen[ksn]=sc->gen[is]+1;
  //marked
  sc->marked=realloc(sc->marked,(nsnew)*sizeof(INT));
  sc->marked[ks0]=sc->marked[is];
  sc->marked[ksn]=sc->marked[is];
  //flags
  sc->flags=realloc(sc->flags,(nsnew)*sizeof(INT));
  sc->flags[ks0]=sc->flags[is];
  sc->flags[ksn]=sc->flags[is];
  //parents
  sc->parent=realloc(sc->parent,(nsnew)*sizeof(INT));
  sc->parent[ks0]=is;
  sc->parent[ksn]=is;
  //child0
  sc->child0=realloc(sc->child0,(nsnew)*sizeof(INT));
  sc->child0[ks0]=-1;
  sc->child0[ksn]=-1;
  //childn
  sc->childn=realloc(sc->childn,(nsnew)*sizeof(INT));
  sc->childn[ks0]=-1; sc->childn[ksn]=-1;
  //volumes
  sc->vols=realloc(sc->vols,(nsnew)*sizeof(REAL));// we can calculate all volumes at the end, so just allocate here.
  //scalars
  sc->ns=nsnew;
  sc->nv=nvnew;
  return 0;
}
/**********************************************************************/
/*!
 * \fn INT haz_refine_simplex(scomplex *sc, const INT is, const INT it)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note The algorithm is found in Traxler, C. T. An algorithm for
 *       adaptive mesh refinement in n-dimensions. Computing 59
 *       (1997), no. 2, 115–137 (MR1475530)
 *
 */
INT haz_refine_simplex(scomplex *sc, const INT is, const INT it)
{
  INT n=sc->n, nbig=sc->nbig,ns=sc->ns,nv=sc->nv;
  INT nsnew=ns,nvnew=nv;
  INT itype,nodnew=-10;
  INT n1=n+1,j,i,p,p0,pn,isn1,snbri,snbrp,snbrn,snbr0;
  INT jt,jv0,jvn,ks0,ksn,s0nbri,snnbri;//,isn;
  REAL *xnew;
  INT pv[2];
  INT csysnew,ibnew;
  if(is<0) return 0;
  if(sc->child0[is] >= 0) return 0;
  //  isn=is*n;
  isn1=is*n1;
  for (i=1;i<n;i++){
    snbri=sc->nbr[isn1+i] ; // the on-axis neighbor.
    if(snbri<0) continue; //this is a boundary
    if (sc->gen[snbri]<sc->gen[is]){//this was wrong in the code in the Traxler's paper
      haz_refine_simplex(sc,snbri,-1);
      nsnew=sc->ns;
      nvnew=sc->nv;
    }
  }
  if (it<0){
    xnew = (REAL *)calloc(nbig,sizeof(REAL));
    jv0=sc->nodes[isn1];
    jvn=sc->nodes[isn1+n];
    // here we can store the edge the new vertex comes from
    for(j=0;j<nbig;j++){
      xnew[j] = 0.5*(sc->x[jv0*nbig+j]+sc->x[jvn*nbig+j]);
    }
    // parents of the vertex:
    pv[0]=jv0;
    pv[1]=jvn;
    // boundary codes (these are also fixed later when connected components on the boundary are found.
    /* if(sc->bndry[jv0] > sc->bndry[jvn]) */
    /*   ibnew=sc->bndry[jv0]; */
    /* else */
    /*   ibnew=sc->bndry[jvn]; */
    ibnew=0; //added vertex is considered interior vertex by default. ;
    if(sc->csys[jv0] < sc->csys[jvn])
      csysnew=sc->csys[jv0];
    else
      csysnew=sc->csys[jvn];
    /* we have added a vertex, let us indicate this */
    nodnew = nvnew;
    nvnew++;
  } else {
    //    jx=    newvertex=t->child0->vertex[1];
    jt=sc->child0[it]; // child simplex number
    jv0 = sc->nodes[jt*n1+1]; // global vertex number of vertex 1.
    xnew = (sc->x + jv0*nv); /* xnew is the same pointer, the vertex
				has already been added. */
    nodnew=jv0;
    ibnew=sc->bndry[nodnew];
    csysnew=sc->csys[nodnew];
  }
  ks0=nsnew;
  sc->child0[is]=ks0; // child0 simplex number
  nsnew++;
  ksn=nsnew;
  sc->childn[is]=ksn; // childn simplex number
  nsnew++;
  /*
    Add two new simplices and initialize their parents, etc
  */
  haz_add_simplex(is,sc,xnew,pv,ibnew,csysnew,nsnew,nvnew);
  /*
    Initialize all vertex pointers of the children according to the
    scheme. Also initialize all pointers to bring simplices as long
    as they do not require a recursive call for subdivision.  Always
    remember: The mesh is supposed to meet the structural condition when
    this routine is called, and it will meet the structural condition
    again after this routine has terminated.
  */
  INT isc0=-100,iscn=-100;
  //ks0 = sc->child0[is];
  //ksn = sc->childn[is];
  isc0=ks0*n1;
  iscn=ksn*n1;
  sc->nodes[isc0+0]=sc->nodes[isn1];   // nodes[is][0]
  sc->nodes[iscn+0]=sc->nodes[isn1+n]; // nodes[is][n1] nodes is (ns x n1)
  /*backup:   sn->nodes[1]=s0->nodes[1]=nodnew;*/
  sc->nodes[iscn+1]=sc->nodes[isc0+1]=nodnew;
  sc->nbr[isc0]=sc->childn[is];
  sc->nbr[iscn]=sc->child0[is];
  /*
    We know the structure of the nbrs children, if existent already,
    thanks to the structural condition.
  */
  snbr0=sc->nbr[isn1];
  snbrn=sc->nbr[isn1+n];
  if(snbrn>=0){
    if (sc->child0[snbrn]>=0){
      //      s0->nbr[1]=sc->child0[snbrn];
      sc->nbr[isc0+1]=sc->child0[snbrn];
    } else {
      //      s0->nbr[1]=snbrn;
      sc->nbr[isc0+1]=snbrn;
    }
  }
  if(snbr0>=0){
    if(sc->childn[snbr0]>=0)
      //      sn->nbr[1]=sc->childn[snbr0];
      sc->nbr[iscn+1]=sc->childn[snbr0];
    else
      //      sn->nbr[1]=snbr0;
      sc->nbr[iscn+1]=snbr0;
  }
  /* Compute the simplex type. */
  itype=(sc->gen[is]) % n;
  // for each vertex p=1..n-1 of the parent simplex S
  for (p=1;p<n;p++)
    {
      /*
	 p0 is the index of S->vertex[p] in child0
	 pn is the index of S->vertex[p] in childn
      */
      pn = p+1;
      p0 = p+1;
      if (p > itype) pn=n-p+itype+1;
      /*       YYYYYYYYYYYYYYYYYYYYYY */
      /* initialize children % vertex according to structural cond */
      //      sn->nodes[pn]=s0->nodes[p0]=sc->nodes[isn1+p];
      sc->nodes[isc0+p0]=sc->nodes[isn1+p];
      sc->nodes[iscn+pn]=sc->nodes[isn1+p];
      snbrp=sc->nbr[isn1+p]; /* s->nbr[p] */
      if(snbrp<0) continue;
      if (sc->child0[snbrp]>=0) {
	/*
	   S-> nbr[p] is subdivided. The corresponding nbr pointers of the
	   children of S should then point to S nbr[p] -> children (structural
	   condition). It might however be that the vertices 0 and n of S have
	   exchanged local indices in the nbr. That has to be tested.
	*/
	if (sc->nodes[snbrp*n1+0]==sc->nodes[isn1+0]) {
	  //	  s0->nbr[p0]=sc->child0[snbrp];
	  //	  sn->nbr[pn]=sc->childn[snbrp];
	  sc->nbr[isc0+p0]=sc->child0[snbrp];
	  sc->nbr[iscn+pn]=sc->childn[snbrp];
	} else {
	  //	  s0->nbr[p0]=sc->childn[snbrp];
	  //	  sn->nbr[pn]=sc->child0[snbrp];
	  sc->nbr[isc0+p0]=sc->childn[snbrp];
	  sc->nbr[iscn+pn]=sc->child0[snbrp];
	}
      } else {
	/*
	  s->nbr[p] is not subdivided. The corresponding neighbor pointers
	  of the children of s are now set to s->nbr[p], which will be
	  corrected later, when this simplex is divided as well.
	*/
	//	sn->nbr[pn]=snbrp;
	//	s0->nbr[p0]=snbrp;
	sc->nbr[isc0+p0]=snbrp;
	sc->nbr[iscn+pn]=snbrp;
      }
    }
  for(i=0;i<n1;i++) {
     /*
	The children OF NEIGHBORS, if existent, still point to s as
	one of their nbrs. This contradicts the structural condition
	and is corrected here.
     */
    //    s0nbri=s0->nbr[i];
    //    snnbri=sn->nbr[i];
    s0nbri=sc->nbr[isc0+i];    /*s->child0->neighbor[i]*/
    snnbri=sc->nbr[iscn+i]; /*s->childn->neighbor[i]*/
    if(s0nbri>=0){
      if(s0nbri >=sc->ns) {
	fprintf(stderr,"\n\nSTOPPING: nsnew,s0nbri,snnbri,ns: %i %i %i %i\n\n",nsnew,snnbri,s0nbri,sc->ns); fflush(stdout);
	exit(32);
      }
      //      if(sc->gen[s0nbri]==s0->gen)
      if(sc->gen[s0nbri]==sc->gen[ks0])
	/* sc->nbr[s0nbri*n1+i] = s->child0->neighbor[i]->neighbor[i]*/
	sc->nbr[s0nbri*n1+i]=ks0;
    }
    if(snnbri>=0){
      if(snnbri >=sc->ns) {
	fprintf(stderr,"\n\nSTOPPING2: s0nbri,snnbri,ns: %i %i %i %i\n",nsnew,snnbri,s0nbri,sc->ns); fflush(stdout);
	exit(33);
      }
      //      if(sc->gen[snnbri]==sn->gen)
      if(sc->gen[snnbri]==sc->gen[ksn])
	/* sc->nbr[snnbri*n1+i] = s->childn->neighbor[i]->neighbor[i]*/
	sc->nbr[snnbri*n1+i]=ksn;
    }
  }
  /*
     NOW call the on-axis nbrs for refinement, passing to them a
     pointer to this simplex S so they find our new vertex.  Skip the
     neighbors opposite to x0 and xn, they do not have to be refined
     and refine the rest of the "on-axis" neighbors */
  for(i=1 ; i < n; i++){
    haz_refine_simplex(sc, sc->nbr[isn1+i], is);
  }
  return 0;
}
/******************************************************************/
/*!
 * \fn void refine(const INT ref_levels, scomplex *sc,ivector *marked)
 *
 * \brief Refines a simplicial complex.
 *
 * \param sc: scomplex containing the whole hierarchy of refinements
 *
 * \param marked: input ivector containing all simplices from the
 *               finest level which are marked for refinement; If
 *               marked is null then uniformly renfines the grid
 *               ref_levels.
 *
 *
 * \param ref_levels number of refinements. If marked is not null,
 * then this is ignored and only one refinement is done;
 *
 * \return void
 *
 */
void refine(const INT ref_levels, scomplex *sc,ivector *marked)
{
  if(ref_levels<=0) return;
  /*somethind to be done*/
  INT j=-1,i,nsold,print_level=0,nsfine=-1;
  if(!sc->level){
    /* sc->level this is set to 0 in haz_scomplex_init */
    /* form neighboring list on the coarsest level */
    find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
    INT *wrk=calloc(5*(sc->n+2),sizeof(INT));
    /* construct bfs tree for the dual graph */
    abfstree(0,sc,wrk,print_level);
    if(wrk) free(wrk);
  }
  if((!marked)){
    // just refine everything that was not refined:
    for (i=0;i<ref_levels;i++){
      nsold=sc->ns;
      for(j = 0;j < nsold;j++)
	if((sc->child0[j]<0||sc->childn[j]<0))
	  haz_refine_simplex(sc, j, -1);
      sc->level++;
    }
    for(j=0;j<sc->ns;j++) sc->marked[j]=TRUE; // not sure we need this.
    // we are done here;
    return;
  } else if(!sc->level){
    // we have not refined anything yet and marked is set, so
    if((marked->row)&&(marked->val))
      for(j=0;j<sc->ns;j++) sc->marked[j]=marked->val[j];
    else {
      //      issue a warning and mark everything for refinement;
      for(j=0;j<sc->ns;j++) sc->marked[j]=TRUE;
    }
  } else {
    /* we come here if we have refined few times and in such case we
       need to re-mark our simplices on the finest level using the
       values of marked at abs(sc->child0[j]+1) which were set from
       the previous visit here */
    for(j=0;j<sc->ns;j++) {
      if(sc->child0[j]<0||sc->childn[j]<0){
	nsfine=abs(sc->child0[j]+1);
	sc->marked[j]=marked->val[nsfine];
      }
    }
  }
  /*
   * refine everything that is marked on the finest level and is
   * not yet refined: (marked>0 and child<0)
   */
  nsold=sc->ns;
  for(j = 0;j < nsold;j++)
    if(sc->marked[j] && (sc->child0[j]<0||sc->childn[j]<0))
      haz_refine_simplex(sc, j, -1);
  /*
   *  compute volumes (the volumes on the coarsest grid should be set in
   * generate_initial_grid, but just in case we are coming here from
   * some other function we compute these here too.
   */
  INT node,in1;
  void *wrk1=malloc((sc->n+1)*(sc->n*sizeof(REAL) + sizeof(INT)));
  REAL *xs=calloc((sc->n+1)*sc->n,sizeof(REAL));
  for(i=0;i<sc->ns;++i){
    if(sc->child0[i]<0||sc->childn[i]<0){
      in1=i*(sc->n+1);
      for (j = 0;j<=sc->n;++j){
	node=sc->nodes[in1+j];
	memcpy((xs+j*sc->n),(sc->x+node*sc->n),sc->n*sizeof(REAL));
      }
      sc->vols[i]=volume_compute(sc->n,sc->factorial,xs,wrk1);
    }
  }
  free(wrk1);
  free(xs);
  sc->level++;
  return;
}
/******************************************************************/
/*!
 * \fn void sc2mesh(scomplex *sc,mesh_struct *mesh)
 *
 * \brief copies scomplex structure to mesh struct (not all but the
 *        bare minimum needed to define mesh_struct.
 *
 * \param scomplex sc;
 *
 * \param mesh_struct mesh;
 *
 *
 */
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
mesh_struct sc2mesh(scomplex *sc)
{
  /*********************************************************************/
  /* INPUT: pointer to a simplicial complex sc; returns mesh_struct
     mesh.

     Store the finest mesh in mesh_struct structure.  sc has all the
    hierarchy, mesh_struct will have only the last mesh.
  */
  /*********************************************************************/
  /* copy the final grid at position 1*/
  mesh_struct mesh;
  initialize_mesh(&mesh);
  INT ns=0,nv=sc->nv,n1=sc->n+1,dim=sc->n,dimbig=sc->nbig;
  if(dimbig!=dim)
    return mesh;
  INT jk=-10,k=-10,j=-10,i=-10;
  /*
    count the number of elements on the last level of refinement.
    On the last grid are all simplices that were not refined, so
    these are the ones for which child0 and childn are not set.
  */
  ns=0;
  for (j=0;j<sc->ns;j++)
    if(sc->child0[j]<0 || sc->childn[j]<0) ns++;
  /*Update mesh with known quantities*/
  mesh.dim = sc->n;
  mesh.nelm = ns; //do not ever put sc->ns here
  mesh.nv=nv;
  mesh.nconn_reg = sc->cc; //
  mesh.nconn_bdry = sc->bndry_cc;// the so called number of holes is this number minus 1.
  mesh.v_per_elm = n1;
  /*Allocate all pointers needed to init the mesh struct*/
  // these are initialized to NULL, so we can use realloc.
  mesh.el_flag = (INT *)realloc(mesh.el_flag,ns*sizeof(INT));
  //mesh.el_vol = (REAL *)realloc(mesh.el_vol,ns*sizeof(REAL));
  mesh.v_flag = (INT *)realloc(mesh.v_flag,nv*sizeof(INT));
  mesh.cv=allocatecoords(nv,dim);
  mesh.el_v=(iCSRmat *)malloc(sizeof(iCSRmat));
  mesh.el_v[0]=icsr_create(mesh.nelm,mesh.nv,mesh.nelm*mesh.v_per_elm);
  free(mesh.el_v->val); mesh.el_v->val=NULL;
  INT chk=(INT )(!(mesh.el_flag && mesh.v_flag && \
		 mesh.cv && mesh.cv->x && \
		 mesh.el_v && mesh.el_v->IA && mesh.el_v->JA
		   ));
  if(chk){
    fprintf(stderr,"\nCould not allocate memory for mesh in %s\n", \
	    __FUNCTION__);
    return mesh;
  }
  /********************************************************************/
  /*element flag and element volumes; volumes are recomputed later in
    build_mesh_all()*/
  INT nsfake=0;// this one must be ns at the end.
  for (j=0;j<sc->ns;j++){
    if(sc->child0[j]<0 || sc->childn[j]<0){
      mesh.el_flag[nsfake]=sc->flags[j];
      //mesh.el_vol[nsfake]=sc->vols[j];
      nsfake++;
    }
  }
  /*boundary flags*/
  mesh.nbv=0;
  for (j=0;j<mesh.nv;j++){
    if(sc->bndry[j]!=0) mesh.nbv++;
    mesh.v_flag[j]=sc->bndry[j];
  }
  /*Coordinates*/
  for(j=0;j<mesh.dim;j++){
    for(i=0;i<nv;i++){
      mesh.cv->x[j*nv+i]=sc->x[i*sc->n+j];
    }
  }
  // el_v
  mesh.el_v->IA[0]=0;
  for (j=0;j<mesh.nelm;j++)
    mesh.el_v->IA[j+1]=mesh.el_v->IA[j]+n1;
  jk=0;
  for (j=0;j<sc->ns;j++){
    /*  copy el_v map using only the top grid;    */
    if(sc->child0[j]<0 || sc->childn[j]<0){
      for (k=0;k<n1;k++)
        memcpy(mesh.el_v->JA+jk*n1,sc->nodes+j*n1,n1*sizeof(INT));
      jk++;
    }
  }
  return mesh;
}
/*********************************************************************/
/*!
 * \fn scomplex *sc_bndry(scomplex *sc)
 *
 * \brief creates a boundary simplicial complex from a given simplicial complex.
 *
 * \param sc         I: the simplicial complex whose boundary we want to find
 *
 * \return  the coundary simplicial complex
 * \note    future: should be combined with find_cc_bndry_cc() in amr_utils.c.
 */
scomplex sc_bndry(scomplex *sc)
{
  scomplex dsc;
  INT ns = sc->ns,nv=sc->nv,dim=sc->n;
  /*
     first find the neighbors so that we have a consistent ordering of
     the vertices and faces. These may already be found, but if not we
     do it again to be sure the ordering is consistent: neighbor
     sharing face [j] is stored at the same place as vertex [j] in
     nodes[:,j]
   */
  find_nbr(ns,nv,dim,sc->nodes,sc->nbr);
  /**/
  INT dim1=dim+1,i,j,k,m,ns_b1,in1,jn1;
  INT ns_b=0;
  for(i=0;i<ns;i++){
    for(j=0;j<dim1;j++){
      if(sc->nbr[i*dim1+j]<0)
	ns_b++;
    }
  }
  dsc.ns=ns_b;
  // end: computation of the nuumber of boundary faces. now store the vertices for every face.
  INT *fnodes=calloc(ns_b*dim,sizeof(INT));
  ns_b1=0;
  for(i=0;i<ns;i++){
    for(j=0;j<dim1;j++){
      if(sc->nbr[i*dim1+j]<0) {
	k=0;
	for(m=0;m<dim1;m++){
	  if(m==j) continue;
	  fnodes[ns_b1*dim+k]=sc->nodes[i*dim1+m];
	  k++;
	}
	ns_b1++;
      }
    }
  }
  if(ns_b!=ns_b1){
    fprintf(stderr,"\n%%***ERROR(65): num. bndry faces mismatch (ns_b=%d .ne. ns_b=%d) in %s",ns_b1,ns_b,__FUNCTION__);
    exit(65);
  }
  // FIX numbering, there is global numbering of nodes and local numbering of nodes:
  INT *indx=calloc(sc->nv,sizeof(INT));
  INT *indxinv=calloc(sc->nv,sizeof(INT));
  memset(indxinv,0,sc->nv*sizeof(INT));
  // memset(indx,0,sc->nv*sizeof(INT));
  for(i=0;i<sc->nv;++i) indx[i]=-1;
  for(i=0;i<ns_b*dim;++i) indx[fnodes[i]]++;// for interior pts this does not get incremented
  INT nv_b=0;
  for(i=0;i<sc->nv;++i){
    if(indx[i]<0) continue;// interior pt
    indx[i]=nv_b;
    indxinv[nv_b]=i;
    nv_b++;
  }
  fprintf(stdout,"\n%%number of boundary vertices=%d (total nv=%d)\n",nv_b,sc->nv);
  dsc=haz_scomplex_null((sc->n-1),sc->n);
  dsc.nv=nv_b;
  dsc.ns=ns_b;
  ////////////////
  if(dsc.nbig>dsc.n){
    fprintf(stdout,"\n%%%%In %s:Simplicial complex of dimension %d embedded in sc of dimension %d\n\n",__FUNCTION__,dsc.n,dsc.nbig);
  }
  // there are things we dont need:
  //  free(dsc.nodes);
  dsc.nodes=fnodes;
  if(dsc.nv<sc->nv)
    indxinv=realloc(indxinv,dsc.nv*sizeof(INT));
  // now we can init the complex and then remove the unnecessary stuff:
  for(i=0;i<dsc.ns*(dsc.n+1);++i)
    fnodes[i]=indx[fnodes[i]];
  // set x, sc->bndry and so on:
  dsc.x=(REAL *) calloc(dsc.nv*(dsc.nbig),sizeof(REAL));
  for(i=0;i<dsc.nv;i++){
    in1=i*dsc.nbig;
    j=indxinv[i];
    jn1=j*sc->n; // original coords:
    memcpy((dsc.x+in1),sc->x+jn1,sc->n*sizeof(REAL));
  }
  free(indx);
  free(indxinv);
  haz_scomplex_realloc(&dsc);
  find_nbr(dsc.ns,dsc.nv,dsc.n,dsc.nodes,dsc.nbr);
  return dsc;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*EOF*/
