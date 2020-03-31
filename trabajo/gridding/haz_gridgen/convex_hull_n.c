/*! \file src/amr/scomplex.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note containing all essentials routines for mesh refinement
 *
 */
#include "hazmath.h"
void visible(INT root, INT node, iCSRmat *stos, \
                iCSRmat *fv, SHORT *maskf, \
                INT dim, REAL fact,     \
                REAL *x, REAL **xmass, REAL **sn, void *wrk);
/************************************************************/
void lexsort(REAL *x, INT nr, INT nc, INT *p, INT *invp);
/*===============================================*/
#ifndef DIM
  #define DIM 2
#endif
#ifndef USE_RANDOM
  #define USE_RANDOM 0
#endif
#ifndef MAX_NV
  #define MAX_NV 1024
#endif
/******************************************************************/
static INT fmaxr(REAL *a,INT n)
{
  INT i,imax=0;
  REAL amax=a[0];
  for(i=1;i<n;i++){
    if(a[i]>amax){
      amax=a[i];
      imax=i;
    }
  }
  return imax;
}
/******************************************************************/
static INT fmaxi(INT *ia,INT n)
{
  INT i,imax=0,amax=ia[0];
  for(i=1;i<n;i++){
    if(ia[i]>amax){
      amax=ia[i];
      imax=i;
    }
  }
  return imax;
}
/******************************************************************/
static void lpr(FILE* fid,dCSRmat *A,const char *varname)
{
  // print a csr matrix in (i,j,val)
  INT i,j1,j2,j;
  // main loop
  fprintf(fid,"\n%s(rows=%d;cols=%d;nnz=%d)=\n",varname,A->row,A->col,A->nnz);
  for(i=0;i<A->row;i++) {
    j1 = A->IA[i];
    j2 = A->IA[i+1];
    for(j=j1;j<j2;j++) {
      fprintf(fid,"%d\t%d\t%25.16e\n",i,A->JA[j],A->val[j]);
    }
  }
  fprintf(fid,"\n");
  return;
}
/******************************************************************/
static void lpri(FILE* fid,iCSRmat *A,const char *varname)
{
  // print a integer csr matrix in (i,j)
  INT i,j1,j2,j;
  fprintf(fid,"\n%s(rows=%d;cols=%d;nnz=%d)=\n",varname,A->row,A->col,A->nnz);
  if(A->val){
    for(i=0;i<A->row;i++) {
      j1 = A->IA[i];
      j2 = A->IA[i+1];
      for(j=j1;j<j2;j++) {
        fprintf(fid,"%d\t%d\t%d\n",i,A->JA[j],A->val[j]);
      }
    }
  } else {
    for(i=0;i<A->row;i++) {
      j1 = A->IA[i];
      j2 = A->IA[i+1];
      for(j=j1;j<j2;j++) {
        fprintf(fid,"%d\t%d\t\n",i,A->JA[j]);
      }
    }
  }
  fprintf(fid,"\n");
  return;
}
/************************************************************/
REAL *massc(INT n, INT nv,REAL *x,INT *svi)
{
  // compute the mass center of nv points with indexes given in svi[]
  INT k,l,kmove;
  REAL invdim=1e0/((REAL )n);
  REAL *xmass=(REAL *)calloc(n,sizeof(REAL));
  memset(xmass,0,n*sizeof(REAL));
  for(k=0;k<nv;k++){
    kmove=svi[k]*n;
    for(l=0;l<n;l++){
      xmass[l]+=x[kmove+l]*invdim;
    }
  }
  return xmass;
}
/************************************************************/
static unsigned int ridges(INT n,INT pos,INT *sv1,INT *sv2,INT *stosloc)
{
  //from neighboring simplices finds the vertices opposite to the shared face
  // and places the neighboring simplices at this position.
  unsigned int fnd;
  unsigned int nf=0;
  INT i1,k1=-10;
  INT i,j,swp;
  for (i=0;i<n;i++){
    fnd=0;
    i1=sv1[i];
    for (j = 0;j < n;j++){
      if(i1 == sv2[j]){
	       fnd=1;
	       break;
      }
    }
    if(fnd) continue;
    else {
      // not found
      nf++;
      if(nf > 1) {k1=-11; break;}
      k1=i;
    }
  }
  /* NOW put the neighbors at the right places ******* */
  if(k1<0){
    fprintf(stderr,"\n***ERROR in %s; k1=%d\n\n",__FUNCTION__,k1);
    exit(65);
  } else {
    fprintf(stdout,"\nnewpos=%d ; oldpos=%d\nstosloc=",k1,pos);
    for(j=0;j<n;j++) fprintf(stdout," %d ",stosloc[j]);
    fprintf(stdout,"\n");
    swp=stosloc[pos];stosloc[pos]=stosloc[k1];stosloc[k1]=swp;
    for(j=0;j<n;j++) fprintf(stdout," %d ",stosloc[j]);
    fprintf(stdout,"\n222stos1=");
    fprintf(stdout,"\n");
    return 1;
  }
}
/************************************************************/
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
/******************************************************************/
void init_s(INT dim, REAL fact, INT nv,REAL *x, INT *nodes, REAL **snsn, \
            INT *perm, INT *invp, SHORT *mask, void *wrk)
{
  // chose the first dim nodes and then the n+1 so that t he resulting simplex has the largest volume
  INT dim1=dim+1,node,i,j,k,l;
  REAL vols,vmax;
  for(i=0;i<dim1;i++) nodes[i]=perm[i];
  // choose initial simplex for constructing the convex hill
  REAL *xf=(REAL *)calloc(dim1*dim,sizeof(REAL)); // coordinates of vertices of the initial tet
  REAL *snloc=(REAL *)calloc(dim1*dim,sizeof(REAL)); // coordinates of vertices of the initial tet
  REAL *areas=(REAL *)calloc(dim1,sizeof(REAL));
  // length of the work space needed : 3 matrices(dim*dim), mass
  // center(dim); pivots(dim),integer permutation(dim) ;
  for(j=0;j<dim1;j++){ //loop around nodes;
      i=nodes[j];
      memcpy((xf+j*dim),(x+i*dim),dim*sizeof(REAL));
  }
  SHORT iz=area_face0(dim,fact,xf,snloc,areas,&vmax,wrk);
  for(j=0;j<dim1;j++)
    memcpy(snsn[j],(snloc+j*dim),dim*sizeof(REAL));
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
      for(j=0;j<dim1;j++)
        memcpy(snsn[j],(snloc+j*dim),dim*sizeof(REAL));
      // fprintf(stdout,"\n(k)=(%i); vmax=%e,vols=%e",k,vmax,vols);
      vmax=vols;
      node=k;
    }
//    fprintf(stdout,"\n%i; %e; %e",k,vols,vmax);
    i++;
  }
  nodes[0]=node;
  INT sw;
  //we now swap p[0:dim] with what is needed
  for (j=0;j<dim1;j++){
    k=invp[nodes[j]];
    // swap p[k] and p[j]
    sw=perm[j]; perm[j]=perm[k]; perm[k]=sw;
    mask[perm[j]]=1;
  }
  for(j=0;j<nv;j++) invp[perm[j]]=j;
//  print_full_mat_int(dim1,1,nodes,"final_nodes");
//  print_full_mat(dim1,dim,snsn,"normals");
  if(snloc) free(snloc);
  if(areas) free(areas);
  if(xf) free(xf);
  return;
}
/******************************************************************/
iCSRmat *find_nbr0(INT ns,INT nv,INT n,iCSRmat *sv)
{
// find neighboring list stos
  INT i,j,k,jp,jia,jib,ibeg,iend,ibbeg,ibend,kn1;
  INT n1 = n+1,nv1=nv+1;
  INT *stosi=NULL,*stosk=NULL,*svi=NULL,*svk=NULL;
  /* allocate */
  iCSRmat *stos=(iCSRmat *)malloc(1*sizeof(iCSRmat));
  iCSRmat *vs=(iCSRmat *)malloc(1*sizeof(iCSRmat));
  /* transpose to get the vertex--simplex incidence */
  icsr_trans(sv,vs);
//  lpri(stdout,vs,"vs");
  icsr_mxm(sv,vs,stos);
  icsr_free(vs);
//  lpri(stdout,stos,"stos1");
// remove all entries not equal to n;
  ibeg=stos->IA[0];
  kn1=ibeg;
  for (i = 0; i < ns; i++){
    iend=stos->IA[i+1];
    for(j = ibeg; j < iend; j++){
      if(stos->val[j]==n) {
        stos->val[kn1]=stos->val[j];
        stos->JA[kn1]=stos->JA[j];
        kn1++;
      }
    }
    ibeg=stos->IA[i+1];
    stos->IA[i+1]=kn1;
  }
  stos->nnz=stos->IA[ns];
  if(vs) icsr_free(vs);
  stos->JA=realloc(stos->JA,stos->nnz*sizeof(INT));
  stos->val=realloc(stos->val,stos->nnz*sizeof(INT));
  // Now we order each entry to correspond to th sv array.
  for (i = 0; i < ns; i++){
    ibeg=stos->IA[i];
    iend=stos->IA[i+1];
    stosi=stos->JA+ibeg; svi=sv->JA+sv->IA[i];
    for(j = ibeg; j < iend; j++){
      k=stos->JA[j];
      svk=sv->JA+sv->IA[k];
      kn1=j-ibeg;
      ridges(n1,kn1,svi,svk,stosi);
    }
  }
  lpri(stdout,sv,"sv");
  lpri(stdout,stos,"stos2");
  return stos;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
int main(void)
{
  INT i,j,k,l,m,ln,nfaces,nv=8,ns=-1;
  INT use_random=USE_RANDOM;
  INT  dim=DIM;
  INT dim1=dim+1,dim_m1=dim-1;
  char *nameio=strndup("coord.mesh",10);
  FILE *fio;
  REAL scpro;
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
    memset(x,0,nv*dim*sizeof(REAL));
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
  // for(i=0;i<nv;i++){
  //   scpro=0e0;
  //   for(j=0;j<dim_m1;j++)
  //     scpro+=x[i*dim+j]*x[i*dim+j];
  //   x[i*dim+dim-1]=scpro;
  // }
  //  print_full_mat(nv,dim,x,"x");
  //  exit(4);
  //lexsort(x,nv,dim, p,invp);
  /*INITIALIZATION*/
  ns=1;
  nfaces=dim1*ns;
  INT *nodes=(INT *)calloc(ns*dim1,sizeof(INT));
  REAL **snsn=(REAL **)calloc(nfaces,sizeof(REAL *));
  for(j=0;j<nfaces;j++)
    snsn[j]=(REAL *)calloc(dim,sizeof(REAL));
  // to keep which vertices are processed
  SHORT *mask=(SHORT *)calloc(nv,sizeof(SHORT));
  memset(mask,0,nv*sizeof(SHORT));
  // initial simplex
  REAL fact=1.;
  for (j=2;j<dim1;j++) fact *= ((REAL )j);
  // work space
  //  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  INT nlength=4*dim1*dim*sizeof(REAL)+2*dim*sizeof(REAL)+dim*sizeof(INT);
  void *wrk=(void *)calloc(nlength,sizeof(char));
  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  init_s(dim,fact,nv,x, nodes,snsn,p,invp,mask,wrk);
  // to keep face2vertex map
  SHORT *maskf=(SHORT *)calloc(nfaces,sizeof(SHORT));
  memset(maskf,0,nv*sizeof(SHORT));
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
    for(j=0;j<fv->IA[nfaces];j++){
      fv->val[j]=1;
    }
    //    lpri(stdout,fv);
    // matrix to hold for each point one face to which this point is in the outside set.
  iCSRmat *stos=(iCSRmat *)find_nbr0(nfaces,nv,dim_m1,fv);
  // to store all mass centers.
  REAL **xmass=(REAL **)calloc(nfaces,sizeof(REAL *));
  INT ibeg,npf;
  for(j=0;j<nfaces;j++){
    ibeg=fv->IA[j];
    npf=fv->IA[j+1]-ibeg;
    xmass[j]=massc(dim,npf,x,(fv->JA+ibeg));
  }
  //  print_full_mat(nfaces,dim,xmass,"mass");
  // matrix: points,faces distances; every point is assigned a face.
  dCSRmat *p2f=malloc(1*sizeof(dCSRmat));
  INT nnzp2f=nv-dim1;
  *p2f=dcsr_create(nv,nfaces,nnzp2f);
  REAL disti,dist0=0e0;
  INT j0=-1;
  iptr=0;
  p2f->IA[0]=iptr;
  for (i=0;i<nv;i++){
    if(mask[i]) {
      p2f->IA[i+1]=iptr;
      continue;
    }
    j0=-1;
    dist0=0e0;
    for(j=0;j<nfaces;j++){
      // construct the tet:
      disti=0e0;
      for (l=0;l<dim;l++){
        // take the inner product of (x-x0).n for this face;
        disti+=snsn[j][l]*(x[i*dim+l]-xmass[j][l]);
      }
      fprintf(stdout,"\n%i::%i;(node,face)=(%i,%i); disti=%e; dist0=%e",iptr,nnzp2f,i,j,disti,dist0);
      if(disti<dist0) continue;
      j0=j;
      dist0=disti;
    }
    if(j0>=0){
      p2f->JA[iptr]=j0;
      p2f->val[iptr]=dist0;
      iptr++;
    }
    p2f->IA[i+1]=iptr;
  }
  INT kbegin,kend,q,iq,nbr,iqmax=2*dim1; // some number for neighbors...
  INT *vset=calloc(iqmax,sizeof(INT));
  lpr(stdout,p2f,"p2f");
  dCSRmat *f2p=malloc(1*sizeof(dCSRmat));
  j0=dcsr_trans(p2f,f2p);
  lpr(stdout,f2p,"f2p");
  // main loop around the faces.
  i=0;
  while(1){
    if(i>=nfaces) break;
    j0=f2p->IA[i];
    l=f2p->IA[i+1]-j0;// length;
    j=j0+fmaxr((f2p->val+j0),l);
    iq=f2p->JA[j];// point at a maximal distance from the outside list of face i.
    visible(i,iq, stos, \
            fv, maskf,  \
            dim, fact,  \
            x, xmass, snsn,wrk);
    break;
    fprintf(stdout,"\n(face,node)=(%d,%d)\t%e",i,m,f2p->val[j]);
    i++;// main loop around faces.
  }
  if(vset) free(vset);
  if(mask) free(mask);
  if(maskf) free(maskf);
  if(p2f) dcsr_free(p2f);
  if(f2p) dcsr_free(f2p);
  if(fv) icsr_free(fv);
  fprintf(stdout,"\n\n");
  return 0;
}
