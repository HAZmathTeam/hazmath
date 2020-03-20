/*! \file src/amr/scomplex.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note containing all essentials routines for mesh refinement
 *
 */
#include "hazmath.h"
void lexsort(REAL *x, INT nr, INT nc, INT *p, INT *invp);
/*===============================================*/
#ifndef USE_RANDOM
  #define USE_RANDOM 0
#endif
/*===============================================*/
#ifndef MAX_NV
  #define MAX_NV 1024
#endif
/*===============================================*/
REAL chk_sign0(const int it, const int nbrit)
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
/**********************************************************************/
/*!
 * \fn scomplex *haz_scomplex_init(INT n,INT ns, INT nv)
 *
 * \brief Initialize simplicial complex in dimension n with ns
 *        simplices and nv vertices.
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void vol_simplex0(INT dim, REAL fact, REAL *xf, REAL *volt, void *wrk)
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
  if(lufull(1, dim, volt, bt,p,piv))
    *volt=0e0;
  else
    *volt=fabs(*volt)/fact;
  return;
}
/**********************************************************************/
/*!
 * \fn void faces_cnt(subscomplex *subsc)
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
void faces_cnt0(subscomplex *subsc)
{
  /*
     dim1 here is dim + 1: number of vertices in a simplex This
     function finds face count and creates element to face map. It
     preserves the ordering in sc->nbr, which means that in element
     el, face[el*dim1+k] is oposize to the vertex
     sc->nodes[ele*dim1+k] and is shared with the element number
     sitting at sc->nbr[ele*dim1+k]
     Such ordering (for sc->nbr) is needed for the refinement of an
n-dimensional simplicial grid so it is also used here to construct
"face-face", etc.
*/
  INT dim=subsc->nbig, dim1=dim+1;
  INT is,it,ir,di,dj,j,k;
  scomplex *sc=subsc->parent;
  INT ns=sc->ns;
  INT *nbr=sc->nbr, *elf=subsc->elf;
  INT nf=0;
  /*
      fprintf(stdout,"\n------Elements: vertices, n=%d, %d, %d\n",subsc->parent->ns,subsc->parent->nv,sc->n);fflush(stdout);
  */
  for(it=0;it<ns;it++){
    di=dim1*it;
    for(j=0;j<dim1;j++){
      is=nbr[di+j];
      if(it>is){
	       elf[di+j]=nf;
	       /*
	         in the neighboring element, find the place of i in the
	          neighboring list of j and put the face number in that spot
	          in elf[].
	       */
	        if(is>=0){
	           dj=dim1*is;
	           ir=-1;
	           for(k = 0;k<dim1;k++){
	              ir=nbr[dj+k];
	              if(ir==it) { elf[dj+k]=nf; break; }
	            }
	            if(ir < 0) {
	               fprintf(stderr,						\
		                 "\n\n*** ERROR in %s: incompatible neighbors: %d is neighbor of %d but %d is not a neighbor of %d\n", \
		                 __FUNCTION__,is,it,it,is);
	               exit(127);
	            }
	        }
	      nf++;
      }
    }
  }
  // set the number of simplices in the subsc (faces)
  subsc->ns=nf;
  return;
}
/**********************************************************************/
/*!
 * \fn SHORT area_face0(INT dim, REAL fact, REAL *xf, REAL *sn,  REAL
 *	       *areas,REAL *volt,  void *wrk)
 *
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
 SHORT area_face0(INT dim, REAL fact, REAL *xf, REAL *sn,	\
 	       REAL *areas,REAL *volt,			\
 	       void *wrk)
{
  /*
     computes areas of all faces (n-1) dimensional simplices and their
     normal vectors for a given simplex. it also computes the volume
     of the simplex.

     work space: wrk should be at least dim*(dim+1) REALS and dim
     integers.
  */
  INT dim1 = dim+1,i,j,j1,ln,ln1;
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
//  print_full_mat(dim,1,piv,"piv");
  if(lufull(1, dim, volt, bt,p,piv)) {
    //print_full_mat(dim,dim,bt,"bt");
    volt[0]=0.;
    return 1; // this is a degenerate
  } else
    volt[0]=fabs(volt[0])/fact;
  memset(sn,0,dim1*dim*sizeof(REAL));
  for (j = 1; j<dim1;j++){
    j1=j-1;
    ln=j*dim;
    sn[ln+j1]=-1.;
    solve_pivot(0, dim, bt, (sn+ln), p,piv);
    //areas[] contains the norm of the normal vectors first (this is 1/altitude);
    areas[j]=0.;
    for(i=0;i<dim;i++){
      sn[i]-=sn[ln+i];
      areas[j]+=sn[ln+i]*sn[ln+i];
    }
    areas[j]=sqrt(areas[j]);
  }
  areas[0]=0.; for(i=0;i<dim;i++) areas[0]+=sn[i]*sn[i];
  areas[0]=sqrt(areas[0]);
  //  print_full_mat(dim1,dim,sn,"snsn22");
  //rescale all normal vectors:
  for(j=0;j<dim1;j++){
    ln=j*dim;
    for(i=0;i<dim;i++) { sn[ln+i]/=areas[j]; }
    areas[j]=fabs((*volt)*areas[j])*((REAL )dim);
  }
  //  print_full_mat(dim1,1,areas,"areas");
  //  print_full_mat(dim1,dim,sn,"snsn33");
  return 0;
}
static unsigned int cmp_simplex0(INT n, INT sim1, INT sim2,	\
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
 void find_nbr0(INT ns,INT nv,INT n,INT *sv,INT *stos)
{
  // find neighboring list stos
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
  /* for (i = 0; i < nv; ++i) { */
  /*   iabeg = ivs[i]; */
  /*   iaend = ivs[i+1]; */
  /*   fprintf(stdout,"row %d: ",i+1); */
  /*   for (jia = iabeg; jia < iaend; ++jia) { */
  /*     fprintf(stdout,"%d ",jvs[jia]+1); */
  /*   } */
  /*   fprintf(stdout,"\n"); */
  /* } */
  /*
   */
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
	    if(!cmp_simplex0(n1,i,k,svi,svk,stosi,stosk)) continue;
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
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*EOF*/
void init_s(INT dim,INT nv,REAL *x, INT *nodes, REAL *snsn, INT *perm)
{
  // chose the first dim nodes and then the n+1 so that t he resulting simplex has the largest volume
  INT dim1=dim+1,node,i,j,k,l;
  INT nlength=4*dim1*dim*sizeof(REAL)+2*dim*sizeof(REAL)+dim*sizeof(INT);
  REAL vols,vmax;
  for(i=0;i<dim1;i++) nodes[i]=perm[i];
  // choose initial simplex for constructing the convex hill
  REAL fact=1.;
  for (j=2;j<dim1;j++) fact *= ((REAL )j);
  REAL *xf=(REAL *)calloc(dim1*dim,sizeof(REAL)); // coordinates of vertices of the initial tet
  REAL *snloc=(REAL *)calloc(dim1*dim,sizeof(REAL)); // coordinates of vertices of the initial tet
  REAL *areas=(REAL *)calloc(dim1,sizeof(REAL));
  // length of the work space needed : 3 matrices(dim*dim), mass
  // center(dim); pivots(dim),integer permutation(dim) ;
  void *wrk=(void *)calloc(nlength,sizeof(char));
  for(j=0;j<dim1;j++){ //loop around nodes;
      i=nodes[j];
      memcpy((xf+j*dim),(x+i*dim),dim*sizeof(REAL));
  }
  SHORT iz=area_face0(dim,fact,xf,snsn,areas,&vmax,wrk);
//  fprintf(stdout,"\nvmax=%g",vmax);
//  print_full_mat_int(1,dim1,nodes,"nodes");
//  print_full_mat(1,dim1,areas,"areas");
  i=dim1+1;
  while(i<nv){
    k=perm[i];
    memcpy(xf,(x+k*dim),dim*sizeof(REAL));
//    print_full_mat(dim1,dim,xf,"xf");
    nodes[0]=k;
    area_face0(dim,fact,xf,snloc,areas,&vols,wrk);
    // fprintf(stdout,"\nvmax=%g",vols);
    // print_full_mat_int(1,dim1,nodes,"nodes");
    // print_full_mat(1,dim1,areas,"areas");
    if(fabs(vols)>fabs(vmax)){
      memcpy(snsn,snloc,dim*dim1*sizeof(REAL));
      // fprintf(stdout,"\n(k)=(%i); vmax=%e,vols=%e",k,vmax,vols);
      vmax=vols;
      node=k;
    }
//    fprintf(stdout,"\n%i; %e; %e",k,vols,vmax);
    i++;
  }
  nodes[0]=node;
  print_full_mat_int(dim1,1,nodes,"final_nodes");
  print_full_mat(dim1,dim,snsn,"normals");
  if(wrk) free(wrk);
  if(snloc) free(snloc);
  if(areas) free(areas);
  if(xf) free(xf);
  exit(129);
  return;
}
int main(void)
{
  INT i,j,k,l,m,nfaces,nv=8,ns=-1;
  INT use_random=USE_RANDOM;
  INT  dim=2,dim1=dim+1,dim_m1=dim-1;
  char *nameio=strndup("coord.mesh",10);
  FILE *fio;
  REAL *x=(REAL *)calloc(nv*dim,sizeof(REAL));
  INT *p=(INT *)calloc(2*nv,sizeof(INT));
  INT *invp=p+nv;
  ////////////////////////////////////////
  if(use_random){
  // generate random x and y
    for(j=0;j<dim;j++){
      srand(time(0)+j);
      for (i=0;i<nv;i++){
        x[i*dim+j]=((REAL )rand())/((REAL )RAND_MAX);
      }
    }
  } else {
    /*quadrilateral*/
    x[0*dim+0]=0.8;  x[0*dim+1]=0.75; //x[0*dim+2]=0.;
//    x[0*dim+0]=1.;  x[0*dim+1]=1.; //x[0*dim+2]=0.;
    x[1*dim+0]=0.;  x[1*dim+1]=1.; //x[0*dim+2]=0.;
    x[2*dim+0]=1.;  x[2*dim+1]=0.; //x[0*dim+2]=0.;
    x[3*dim+0]=0.;  x[3*dim+1]=0.; //x[0*dim+2]=0.;
/*L-shaped*/
    x[0*dim+0]= -1.; x[0*dim+1]= -1.; 	// vertex number, coord system, coords
    x[1*dim+0]= -1.; x[1*dim+1]=  0.;
    x[2*dim+0]=  0.; x[2*dim+1]= -1.;
    x[3*dim+0]=  0.; x[3*dim+1]=  0.;
    x[4*dim+0]=  1.; x[4*dim+1]= -1.;
    x[5*dim+0]=  1.; x[5*dim+1]=  0.;
    x[6*dim+0]= -1.; x[6*dim+1]=  1.;
    x[7*dim+0]=  0.; x[7*dim+1]=  1.;
  }
  for(i=0;i<nv;i++) {p[i]=i;invp[i]=i;}
  lexsort(x,nv,dim, p,invp);
  /*INITIALIZATION*/
  ns=1;
  nfaces=dim1*ns;
  INT *nodes=(INT *)calloc(ns*dim1,sizeof(INT));
  REAL *snsn=(REAL *)calloc(dim*nfaces,sizeof(REAL));
  init_s(dim,nv,x, nodes,snsn,p);
  // to keep which vertices are processed
  SHORT *mask=(SHORT *)calloc(nv,sizeof(SHORT));
  memset(mask,0,nv*sizeof(SHORT));
  // to keep face2vertex map
  iCSRmat *fv=(iCSRmat *)malloc(sizeof(iCSRmat));
  *fv=icsr_create(nfaces,nv,nfaces*dim);
  // initialize face vertex map.
  INT iptr=0;
  fv->IA[0]=iptr;
  for(j=0;j<dim1;j++){
      for (i=0;i<dim1;i++){
        if(i==j) continue;
        fv->JA[iptr]=nodes[i];
        mask[nodes[i]]=1;
        iptr++;
      }
      fv->IA[j+1]=iptr;
    }
    icsr_print_matlab(stdout,fv);
//  }
//  INT *stos=(INT *)calloc(nfaces*dim,sizeof(INT));
//  find_nbr0(nfaces,nv,dim_m1,faces,stos);
//  print_full_mat_int(dim1,dim,stos,"nbrs");
  // mass center of the faces
  /*get initial simplex*/
  REAL *xmass=(REAL *)calloc(dim*nfaces,sizeof(REAL));
  memset(xmass,0,dim*nfaces*sizeof(REAL));
  INT iaa,iab;
  for(j=0;j<nfaces;j++){
    iaa=fv->IA[j];
    iab=fv->IA[j+1];
    for(k=iaa;k<iab;k++) {
      i=fv->JA[k];
      for(l=0;l<dim;l++)
        xmass[j*dim+l]+=x[i*dim+l];
    }
    xmass[j]/=((REAL )dim);
  }
  print_full_mat(nfaces,dim,xmass,"mass");
  // matrix points,faces distances; every point is assigned a face.
  dCSRmat *f2p=malloc(1*sizeof(dCSRmat));
  dCSRmat *f2pt=malloc(1*sizeof(dCSRmat));
  *f2p=dcsr_create(nv,nfaces,nv);
  /*ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ*/
  // fio=fopen(nameio,"w");
  // if(fio){
    //    fprintf(stdout,"I/O on file \"%s\"",nameio);
    // for (i=0;i<nv;i++){
    //   fprintf(stdout,"\n%10i: (",i+1);
    //   for(j=0;j<dim-1;j++){
    //     fprintf(stdout,"%16.8e,",x[i*dim+j]);
    //   }
    //   fprintf(stdout,"%16.8e)",x[i*dim+dim-1]);
    // }
    for (i=0;i<nv;i++){
      fprintf(stdout,"\n%10i-->%10i: (",i,p[i]);
      for(j=0;j<dim-1;j++){
        fprintf(stdout,"%16.8e,",x[p[i]*dim+j]);
      }
      fprintf(stdout,"%16.8e)",x[p[i]*dim+dim-1]);
    }
    // fclose(fio);
  // } else {
  //   fprintf(stderr,"****ERROR:Could not open file \"%s\" for I/O\n",nameio);
  //   //return;
  //   exit(127);
  // }
  fprintf(stdout,"\n");
  return 0;
}
