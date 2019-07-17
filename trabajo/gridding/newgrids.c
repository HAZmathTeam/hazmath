#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <getopt.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>
#include "hazmath.h"
#include "grid_defs.h"
#ifndef INT
#define INT int
#endif
#ifndef REAL
#define REAL double
#endif
/********************************FINCTIONS:*********************/
void cube2simp_free(cube2simp *c2s);
INT reverse(void *arr,INT length, size_t elsize);
cube2simp *cube2simplex(INT dim);
scomplex *umesh(const INT dim, INT *nd, cube2simp *c2s, const INT intype);
void polar2cart(INT dim, REAL *px, REAL *cx);
REAL interp8(cube2simp *c2s, REAL *u, REAL *ue, REAL *xhat);
REAL interp4(cube2simp *c2s, REAL *u, REAL *xhat);
void unirefine(INT *nd,scomplex *sc);
/*********************************************************************/
void set_edges(input_grid *g0,cube2simp *c2s)
{
  /*adds missing edges to g0->seg*/
  INT *newseg=calloc(3*c2s->ne*g0->nel,sizeof(INT));
  INT *mnodes=calloc(c2s->nvcube,sizeof(INT));
  INT newne,i,j,k,kel,j0,j1,k0,k1,ke,swp,nvcube=c2s->nvcube;
  newne=0;
  INT found=0;
  for(kel=0;kel<g0->nel;kel++){
    memcpy(mnodes,(g0->mnodes+kel*(nvcube+1)),(nvcube+1)*sizeof(INT));
    for(i=0;i<c2s->ne;i++){
      j0=mnodes[c2s->edges[2*i]];
      j1=mnodes[c2s->edges[2*i+1]];
      if(j0>j1){swp=j0;j0=j1;j1=swp;}
      found=0;
      for(ke=0;ke<g0->ne;ke++){
	k0=g0->seg[3*ke];
	k1=g0->seg[3*ke+1];
	if((k0==j0)&&(k1==j1)){
  	  fprintf(stdout,"\n%%FOUND:**** e=(%d,%d),[%d,%d]",j0,j1,k0,k1); fflush(stdout);
	  found=1;break;
	}      
      }
      if(!found){
	newseg[3*newne]=j0;
	newseg[3*newne+1]=j1;
	newseg[3*newne+2]=1;
	newne++;
      }
    }
  }
  fprintf(stdout,"\n**** newne=%d",newne); fflush(stdout);
  if(!newne){
    free(newseg);
    free(mnodes);
    return;
  }
  newseg=realloc(newseg,3*newne*sizeof(INT));
  INT *p=calloc(newne,sizeof(INT));
  ilexsort(newne, 3,newseg,p);  
  //  print_full_mat_int(newne,3,newseg,"newseg");
  free(p);
  // remove dupps
  INT m,li1,ic[3]; //l1 distance
  k=newne-1;
  i=0;j=0;
  while (i<k){
    if(j==0)
      for(m=0;m<3;m++)
  	newseg[m]=newseg[3*i+m];
    for(m=0;m<3;m++)
      ic[m]=newseg[3*j+m];
    while(1 && i<k) {
      li1=0;
      for(m=0;m<3;m++){
  	li1+=abs(ic[m]-newseg[3*i+3+m]);
      }
      if(li1>0){
  	j++;i++;
  	for(m=0;m<3;m++){newseg[3*j+m]=newseg[3*i+m];}
  	break;
      }
      i++;
      //      fprintf(stdout,"i=%i\n",i);
    }
    //    fprintf(stdout,"i=%i, j=%i\n",i,j);
  }
  i++;j++; newne=j;
  //  print_full_mat_int(g0->ne,3,g0->seg,"newseg");
  g0->seg=realloc(g0->seg,(3*(g0->ne+newne))*sizeof(INT));
  memcpy((g0->seg+3*g0->ne),newseg,3*newne*sizeof(INT));
  g0->ne+=newne;
  free(newseg);
  free(mnodes);
  //  print_full_mat_int(g0->ne,3,g0->seg,"newseg");
  return;
}
void map2mac(scomplex *sc,cube2simp *c2s, input_grid *g)
{
  /* 
     maps a uniform grid from the n-dimensional cube to a hexagonal 
     macroelement from an initial grid of macroelements 
     given by its coordinates xmac[1:nvcube*dim]
     xmac[nvcube][dim]
  */
  INT i,j,k1,k2,k1c,kf,dim=sc->n;
  //  INT k2c;
  INT ksys;
  REAL *xmac=g->xv;  
  REAL *xhat = (REAL *)calloc(dim,sizeof(REAL));
  REAL *xemac=(REAL *)calloc(c2s->ne*dim,sizeof(REAL));
  // convert midpoints from polar to cartesian.
  for(i=0;i<c2s->ne;i++){
    k1=c2s->edges[2*i];
    k2=c2s->edges[2*i+1];
    k1c=g->systypes[g->csysv[k1]];
    //    k2c=g->systypes[g->csysv[k2]];
    //    fprintf(stdout,"\nverts=(%d,%d); coord_sys=(%d,%d)",k1,k2,k1c,k2c);fflush(stdout);
    if(g->csysv[k1]==g->csysv[k2] && k1c==1){
      //use xhat as a temp array:
      xhat[0]=0.5*(xmac[k1*dim]+xmac[k2*dim]);// this is rho
      // take half angles;
      for(j=1;j<dim;j++) {
	xhat[j]=0.5*(xmac[k1*dim+j]+xmac[k2*dim+j]);
      }
      polar2cart(dim,xhat,xemac+(i*dim));
      // translate by adding the origin. 
      ksys=g->csysv[k1];// k1c and k2c should be the same below. 
			      //      k2c=g->csysv[k2];
      for(j=0;j<dim;j++) {
	xemac[i*dim+j]+=g->ox[ksys*dim+j];
      }
    }
  }
  // end of mid points in polar;
  // now convert all vertices in cartesian as well. 
  for(i=0;i<c2s->nvcube;i++){
    k1c=g->systypes[g->csysv[i]];
    if(k1c==1){
      memcpy(xhat,xmac+i*dim,dim*sizeof(REAL));
      polar2cart(dim,xhat,xmac+(i*dim));
      //      translate
    }
    ksys=g->csysv[i];
    for(j=0;j<dim;j++) {
      xmac[i*dim+j]+=g->ox[ksys*dim+j];
    }
  }
  // now everything is in cartesian, so midpoints that are
  // not yet attended to are just averages.
  for(i=0;i<c2s->ne;i++){
    k1=c2s->edges[2*i];
    k2=c2s->edges[2*i+1];
    k1c=g->systypes[g->csysv[k1]];
    //    k2c=g->systypes[g->csysv[k2]];
    //skip all polar mid points
    if(g->csysv[k1]==g->csysv[k2] && k1c==1) continue;
    //    fprintf(stdout,"\ncart:verts=(%d,%d); coord_sys=(%d,%d)",k1,k2,k1c,k2c);fflush(stdout);
    for(j=0;j<dim;j++) {
      xemac[i*dim+j]=0.5*(xmac[k1*dim+j]+xmac[k2*dim+j]);
    }
  }
  //  print_full_mat(c2s->nvcube,dim,xmac,"X");
  r2c(c2s->nvcube,dim,sizeof(REAL),xmac); // we need xmac by rows here
  r2c(c2s->ne,dim,sizeof(REAL),xemac); // we need xemac (mid points of
				       // edges) also by rows
  for(kf=0;kf<sc->nv;kf++){
    for(i=0;i<dim;i++)xhat[i]=sc->x[kf*dim+i];
    for(i=0;i<dim;i++){
      //      fprintf(stdout,"coordinate:%d",i);
      sc->x[kf*dim+i]=interp8(c2s,xmac+i*c2s->nvcube,xemac+i*c2s->ne,xhat);
      //      sc->x[kf*dim+i]=interp4(c2s,xmac+i*c2s->nvcube,xhat);
    }
  }
  //  r2c(dim,c2s->nvcube,sizeof(REAL),xmac); // we need xmac by columns here
  //  r2c(dim,c2s->ne,sizeof(REAL),xemac); // we need xemac by rows agin
  if(xhat) free(xhat);
  if(xemac) free(xemac);
  return;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx*/
static INT *set_input_grid1(input_grid *g,cube2simp *c2s)
{
  /* 
     Every edge is put into a subset, i.e. two edges (v(i1),v(i2)) and
     (v(j1),v(j2)) are considered equivalent iff (i2-i1)=(j2-j1).  The
     number of divisions in an equivalent set of edges is taken to be
     the largest from the equivalence class.  OUTPUT array is a "dim"
     array and for each direction gives the number of partitions.
  */
  INT i,j,k,iri,ici,pmem;
  pmem=2*g->nv;
  if(pmem<2*g->ne) pmem=2*g->ne;
  if(pmem<2*g->nel) pmem=2*g->nel;
  if(pmem<2*g->nf) pmem=2*g->nf;
  INT *p=calloc(pmem,sizeof(INT));// permutation and inverse permutation;
  //
  memset(p,0,pmem);
  for (i=0;i<g->ne;i++){
    iri=g->seg[3*i];
    ici=g->seg[3*i+1];
    if(iri<ici){
      g->seg[3*i]=iri;
      g->seg[3*i+1]=ici;
    } else {
      g->seg[3*i]=ici;
      g->seg[3*i+1]=iri;
    }
    /* set up divisions */
    j=g->seg[3*i+1]-g->seg[3*i]; // should be always positive;
    //    fprintf(stdout,"\n%%z123=%d:(%d,%d);%d",i,3*i,3*i+1,g0->seg[3*efound[i]+2]);
    if(g->seg[3*i+2]>p[j])
      p[j]=g->seg[3*i+2];
  } 
  for (i=0;i<g->ne;i++){
    j=g->seg[3*i+1]-g->seg[3*i];
    g->seg[3*i+2]=p[j]; 
    //    fprintf(stdout,"\n[%d,%d]:div=%d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
  }
  for (i=0;i<g->ne;i++){
    j=g->seg[3*i+1]-g->seg[3*i];
    g->seg[3*i+2]=p[j]; 
    //    fprintf(stdout,"\n[%d,%d]:div=%d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
  }
  /*ORDER*/
  ilexsort(g->ne, 3,g->seg,p);
  k=0;
  for (i=0;i<g->ne;i++){
    if(g->seg[3*i]) continue;
    //    j=g->seg[3*i+1]-g->seg[3*i]-1;
    p[k]=g->seg[3*i+2];
    k++;
    //    fprintf(stdout,"\n[%d,%d]:div=%d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
  }
  p=realloc(p,g->dim*sizeof(INT)); // realloc to dimension g->dim
  //  for (i=0;i<g->dim;i++){
  //    fprintf(stdout,"\ndirection:%d; div=%d",i,p[i]);
  //  }
  //  input_grid_print(g);
  //  print_full_mat_int(g->ne,3,g->seg,"med");
  //  print_full_mat_int(g->nf,(c2s->nvface+1),g->mfaces,"mf");
  //  print_full_mat_int(g->nel,(c2s->nvcube+1),g->mnodes,"mel");
  return p;
}
/***********************************************************************/
static INT set_ndiv_edges(input_grid *g,	\
			  input_grid *g0,	\
			  cube2simp *c2s,	\
			  INT **nd,		\
			  const INT iter)
{
  /* 
     For a given global input grid g0 creates local input grids for
     every macroelement and computes the divisions for it. It is used
     iteratively in macro_split to se the correct divisions for every
     macroelement. efound should be of length c2s->ne +g0->ne. The
     input_grid g0 should be all set, and the input_grid g should have
     all its scalar values set.  nd is the array with the divisions,
     it must be g0->nel by c2s->n.
  */
  INT kel0,i,j0,j1,swp,kel,ke,k0,k1,ndiv;
  INT nel0=g0->nel,nvcube=c2s->nvcube;
  /*for easier reference*/
  INT *efound=calloc(c2s->ne*(g0->nel+1),sizeof(INT));
  INT *e0found=efound + c2s->ne;
  for(i=0;i<g0->ne;i++)e0found[i]=-1;
  // make all divisions > 0
  for(ke=0;ke<g0->ne;ke++){
    ndiv=abs(g0->seg[3*ke+2]);
    if(ndiv<=0)ndiv=1;
    g0->seg[3*ke+2]=ndiv;
  }
  for(ke=0;ke<g0->ne;ke++)      
    e0found[ke]=g0->seg[3*ke+2];
  //  print_full_mat_int(g0->ne,3,g0->seg,"seg0");
  for(kel0=0;kel0<nel0;kel0++){
    if((iter%2)) kel=nel0-kel0-1; else kel=kel0;
    // macroelement by macroelement try to find the edge divisions 
    for(i=0;i<c2s->ne;i++){
      g->seg[3*i]=c2s->edges[2*i];
      g->seg[3*i+1]=c2s->edges[2*i+1];
      g->seg[3*i+2]=-1;
      efound[i]=-1;
    }
    memcpy(g->mnodes,(g0->mnodes+kel*(nvcube+1)),(nvcube+1)*sizeof(INT));
    for(i=0;i<c2s->ne;i++){
      j0=g->mnodes[c2s->edges[2*i]];
      j1=g->mnodes[c2s->edges[2*i+1]];
      if(j0>j1){swp=j0;j0=j1;j1=swp;}
      for(ke=0;ke<g0->ne;ke++){
	k0=g0->seg[3*ke];
	k1=g0->seg[3*ke+1];
	if((k0==j0)&&(k1==j1)){
	  g->seg[3*i+2]=g0->seg[3*ke+2];
	  efound[i]=ke;
	}
      }
      //      if(iter==0)
      //fprintf(stdout,"\n%%iter=%d;Element:%d, edge=(%d,%d);",iter,kel,j0,j1);
    }
    //    if(iter==0)
    //      input_grid_print(g);
    nd[kel]=set_input_grid1(g,c2s);
    for(i=0;i<g->ne;i++){      
      if(efound[i]<0)
	continue;
      ke=efound[i];
      g0->seg[3*ke+2]=g->seg[3*i+2];
    }
  }
  INT chng=0;
  for(ke=0;ke<g0->ne;ke++){
    k1=abs(e0found[ke]-g0->seg[3*ke+2]);
    if(k1>chng)chng=k1;
  }
  //  print_full_mat_int(g0->ne,3,g0->seg,"seg1");
  //print_full_mat_int(1,c2s->n,nd[kel],"ndnd");
  //  fprintf(stderr,"\nchng=%d",chng);
  free(efound);
  return chng;
}
/************************************************************************/
iCSRmat *set_mmesh(input_grid *g0,				\
		   cube2simp *c2s,				\
		   ivector *etree0,				\
		   INT **elneib, INT **el2fnum,				\
		   ivector *isbface0, ivector *bcodesf0,		\
		   INT *wrk)
{
  /*prepare macro element mesh (mmesh) for passing to the mesh generator*/
  /*wrk is working integer array should have at least size 2*nvcube+2 */
  INT i,j0,j1,kel,jel,ke;
  INT nel0=g0->nel,nvcube=c2s->nvcube,nvface=c2s->nvface;
  INT je,kj,k2,iel2v,jel2v,k1,kface,kbnd,found;
  INT *p=wrk; 
  INT *mnodes=p+nvcube+1; 
  //  print_full_mat_int(g0->nel,(c2s->nvcube+1),g0->mnodes,"ex");
  ilexsort(g0->nf, (c2s->nvface+1),g0->mfaces,p);
  ilexsort(g0->nel,(c2s->nvcube+1),g0->mnodes,p);
  /*-------------------------------------------------------------------*/
  iCSRmat *el2v=malloc(1*sizeof(iCSRmat));
  el2v[0]=icsr_create(g0->nel,g0->nv,nvcube*g0->nel);
  el2v->IA[0]=0;
  for(kel=0;kel<nel0;kel++){
    memcpy((el2v->JA+el2v->IA[kel]),(g0->mnodes+kel*(nvcube+1)),nvcube*sizeof(INT));
    el2v->IA[kel+1]=el2v->IA[kel]+nvcube;
  }
  for(i=0;i<el2v->IA[nel0];i++)
    el2v->val[i]=1;
  iCSRmat *v2el=malloc(1*sizeof(iCSRmat));
  icsr_trans(el2v,v2el);	
  iCSRmat *el2el=malloc(1*sizeof(iCSRmat));
  icsr_mxm(el2v,v2el,el2el);
  icsr_free(v2el);
  /*   shrink the el2el matrix by removing any entry with value not
       equal to nvface; */
  el2el->nnz=el2el->IA[0];
  for(kel=0;kel<nel0;kel++){    
    j0=el2el->IA[kel];
    j1=el2el->IA[kel+1];
    el2el->IA[kel]=el2el->nnz;
    for(ke=j0;ke<j1;ke++){
      jel=el2el->JA[ke];
      if(jel==kel) {
	el2el->JA[el2el->nnz]=kel;
	el2el->val[el2el->nnz]=1;
	el2el->nnz++;
	continue;
      }
      if(el2el->val[ke]!=nvface) continue;
      el2el->JA[el2el->nnz]=jel;
      el2el->val[el2el->nnz]=el2el->val[ke];
      el2el->nnz++;
    }
  }  
  el2el->IA[nel0]=el2el->nnz;
  el2el->JA=realloc(el2el->JA,el2el->nnz*sizeof(INT));
  el2el->val=realloc(el2el->val,el2el->nnz*sizeof(INT));
  /* find the connected components:*/
  INT *iblk=calloc((nel0+1),sizeof(INT));
  INT *jblk=calloc((nel0),sizeof(INT));
  // number of connected domains and boundaries;
  INT nblkdom,nblkbnd;
  dfs00_(&nel0,el2el->IA, el2el->JA,&nblkdom,iblk,jblk);
  iblk=realloc(iblk,(nblkdom+1)*sizeof(INT));
  fprintf(stdout,"\n%%DFS(domains): %d connected components",nblkdom);
  /*  prepare for bfs: remove diagonal in el2el */
  icsr_nodiag(el2el);
  etree0->row=el2el->row+1;
  etree0->val=calloc(etree0->row,sizeof(INT));
  INT *etree=etree0->val;
  iCSRmat *bfs0=bfscc(nblkdom,iblk,jblk,el2el,etree);
  /* construct element to face to element maps for macroelements */
  /*FACES******************************************************/
  INT nfaceall=g0->nel*c2s->nf-(el2el->nnz/2);
  INT nfacei=(INT )(el2el->nnz/2);
  INT nfaceb=nfaceall-nfacei;
  // face to vertex:
  iCSRmat *f2v=malloc(sizeof(iCSRmat));
  f2v[0]=icsr_create(nfaceall,g0->nv,nvface*nfaceall);
  /**/
  bcodesf0->val=calloc(nfaceall,sizeof(INT));
  isbface0->val=calloc(nfaceall,sizeof(INT));
  INT *bcodesf=bcodesf0->val;
  INT *isbface=isbface0->val;
  isbface0->row=nfaceall;
  bcodesf0->row=nfaceall;
  /*init*/
  /*-------------------------------------------------------------------*/
  for(kel=0;kel<g0->nel;kel++){
    for(i=0;i<c2s->nf;i++){
      elneib[kel][i]=-1;
      el2fnum[kel][i]=-1;
    }
  }
  /*-------------------------------------------------------------------*/
  f2v->IA[0]=0;
  for(i=0;i<nfaceall;i++){
    f2v->IA[i+1]=f2v->IA[i]+nvface;
    isbface[i]=0;
  }
  INT *facei=calloc(2*nvface,sizeof(INT));
  INT *facej=facei+nvface;
  kface=0;
  for(kel=0;kel<nel0;kel++){
    iel2v=el2v->IA[kel];
    for (ke=0;ke<c2s->nf;ke++){
      found=0;
      k1=ke*nvface;
      for(i=0;i<nvface;i++){
	facei[i]=el2v->JA[iel2v+c2s->faces[k1+i]];
      }
      j0=el2el->IA[kel];
      j1=el2el->IA[kel+1];
      for(kj=j0;kj<j1;kj++){
	jel=el2el->JA[kj];
	jel2v=el2v->IA[jel];
	for (je=0;je<c2s->nf;je++){
	  k2=je*nvface;
	  for(i=0;i<nvface;i++){
	    facej[i]=el2v->JA[jel2v+c2s->faces[k2+i]];
	  }
	  if(aresame(facei,facej,nvface)){
	    if(kel<jel){
	      /* fprintf(stdout,"\nptr=%d;FOUND: fi(ei),fj(ej)=%d(%d),%d(%d);",kface,ke,kel,je,jel); */
	      /* print_full_mat_int(1,nvface,facei,"facei"); */
	      /* print_full_mat_int(1,nvface,facej,"facej"); */
	      /* fprintf(stdout,"\nadd nonzero at:%d ",f2v->IA[kface]); */
	      memcpy((f2v->JA+f2v->IA[kface]),facei,nvface*sizeof(INT));
	      elneib[kel][ke]=jel;
	      elneib[jel][je]=kel;
	      el2fnum[kel][ke]=kface;
	      el2fnum[jel][je]=kface;
	      kface++;// pointer for the face2 vertex matrix
	    }
	    found=1;
	  }
	}
      }
      if(!found){
	/* fprintf(stdout,"\nptr=%d;NOT FOUND: fi(ei)=%d(%d);",kface,ke,kel); */
	/* print_full_mat_int(1,nvface,facei,"facei"); */
	/* fprintf(stdout,"\nadd nonzero at:%d",f2v->IA[kface]); */
	isbface[kface]=1;
	memcpy((f2v->JA+f2v->IA[kface]),facei,nvface*sizeof(INT));
/* 
in this way bcodesf[1-elneib[kel][ke]] gives us the code of the corresponding face;; this can also be linked to f2v and so on. 
*/
	el2fnum[kel][ke]=kface;
	kface++;
      }     
    }
  }
  /*****************************************************/
  for(i=0;i<f2v->IA[nfaceall];i++)
    f2v->val[i]=1;
  /*******************************************************************/    
  //  fprintf(stdout,"\nkface=%d  ? = ? %d\n",kface,f2v->nnz);
  /*ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ*/
  for(i=0;i<nfaceall;i++){
    bcodesf[i]=-1;
    j1=f2v->IA[i+1]-f2v->IA[i];
    //    fprintf(stdout,"\nnnz(%d)=%d",i,j1); 
    if(j1>0){
      memcpy(facei,(f2v->JA+f2v->IA[i]),j1*sizeof(INT));
      for(je=0;je<g0->nf;je++){
	memcpy(facej,(g0->mfaces+je*(nvface+1)),nvface*sizeof(INT));
	kbnd=g0->mfaces[je*(nvface+1)+nvface];
	if(aresame(facei,facej,nvface)){
	  bcodesf[i]=kbnd;
	  //	  fprintf(stdout,"\ncode:%d ",bcodesf[i]);
	  //	  print_full_mat_int(1,nvface,facei,"facei");
	  //	  print_full_mat_int(1,nvface,facej,"facej");
	  break;
	}
      }
    }
  }
  /***********************************************************************/
  /*Connected components on the boundaries. First find face to face map*/
  iCSRmat *v2f=malloc(1*sizeof(iCSRmat));
  icsr_trans(f2v,v2f);
  iCSRmat *f2f=malloc(1*sizeof(iCSRmat));
  /*******************************************************************/    
  icsr_mxm(f2v,v2f,f2f);
  /* now remove all rows in f2f that correspond to interior faces and
     all entries that are not nvface/2, i.e. the number of vertices in
     a (n-2) cube; */
  f2f->nnz=f2f->IA[0];
  for(i=0;i<f2f->row;i++){    
    j0=f2f->IA[i];
    j1=f2f->IA[i+1];
    f2f->IA[i]=f2f->nnz;
    //    if(isbface[i]==0) continue;
    for(ke=j0;ke<j1;ke++){
      je=f2f->JA[ke];
      if(je==i){
	f2f->JA[f2f->nnz]=i;
	f2f->val[f2f->nnz]=0;
	f2f->nnz++;
	continue;
      }
      if((isbface[je]==0)||(f2f->val[ke]!=(nvface/2))) continue;
      //      if((je==i)||(isbface[je]==0)||(f2f->val[ke]!=(nvface/2))) continue;
      f2f->JA[f2f->nnz]=je;
      f2f->val[f2f->nnz]=f2f->val[ke];
      f2f->nnz++;
    }
  }    
  f2f->IA[f2f->row]=f2f->nnz;
  f2f->JA=realloc(f2f->JA,f2f->nnz*sizeof(INT));
  f2f->val=realloc(f2f->val,f2f->nnz*sizeof(INT));
  /*******************************************************************/    
  icsr_free(v2f);
  /* icsr_free(el2f); NOT USED */
  /*connected comps on the boundary*/
  iblk=realloc(iblk,(nfaceall+1)*sizeof(INT));
  jblk=realloc(jblk,(nfaceall)*sizeof(INT));
  dfs00_(&nfaceall,f2f->IA, f2f->JA,&nblkbnd,iblk,jblk);
  fprintf(stdout,"\n%%DFS(boundaries): %d connected components",nblkbnd-nfacei);
  //icsr_nodiag(f2f);
  icsr_free(f2f);
  free(iblk); 
  free(jblk);   

  /*now use the bfs*/
  INT lvl,keok,swp,keswp;
  for(lvl=0;lvl<bfs0->row;lvl++){
    j0=bfs0->IA[lvl];
    j1=bfs0->IA[lvl+1];
    for(kj=j0;kj<j1;kj++){
      jel=bfs0->JA[kj];
      kel=etree[jel];// ancestor, this stays unchanged
      if(kel>=0){
	je=locate0(jel,elneib[kel], c2s->nf);
	ke=locate0(kel,elneib[jel], c2s->nf);	  
	//	fprintf(stdout,"\nel=%d on face %d in el %d",jel,je,kel);
	keok=(je+c2s->n)%c2s->nf;
	if(keok!=ke){
	  //	  fprintf(stdout,"\nel=%d on face %d (should be %d) in el *** %d ***",kel,ke,keok,jel);
	  //	  swap in jel:
	  swp=elneib[jel][ke];
	  elneib[jel][ke]=elneib[jel][keok];
	  elneib[jel][keok]=swp;	  
	  swp=el2fnum[jel][ke];
	  el2fnum[jel][ke]=el2fnum[jel][keok];
	  el2fnum[jel][keok]=swp;	  
	  // we now need to swap vertices in g0->mnodes
	  for(i=0;i<nvface;i++){
	    facei[i]=c2s->faces[ke*nvface+i];
	    keswp=(ke+c2s->n)%c2s->nf;
	    facei[i+nvface]=c2s->faces[keswp*nvface+i];
	    // use mnodes as work space here;
	    mnodes[i]=c2s->faces[keok*nvface+i];
	    keswp=(keok+c2s->n)%c2s->nf;
	    mnodes[i+nvface]=c2s->faces[keswp*nvface+i];	    
	  }
	  for(i=0;i<nvcube;i++)
	    el2v->JA[i]=g0->mnodes[jel*(nvcube+1)+facei[i]];	  
	  for(i=0;i<nvcube;i++)
	    g0->mnodes[jel*(nvcube+1)+mnodes[i]]=el2v->JA[i];
	  /* fprintf(stdout,"\nEL=%d: ",jel); */
	  /* for(i=0;i<nvcube;i++){ */
	  /*   k1=g0->mnodes[jel*(nvcube+1)+i]; */
	  /*   fprintf(stdout,"%d ",k1); */
	  /* } */
	}
      }
    }
  }
  print_full_mat_int(g0->nel,c2s->nvcube+1,g0->mnodes,"mel0");
  /**** FINAL REORDER: make the vertices order in shared faces the
	same!!! ***/
  for(lvl=0;lvl<bfs0->row;lvl++){
    j0=bfs0->IA[lvl];
    j1=bfs0->IA[lvl+1];
    for(kj=j0;kj<j1;kj++){
      jel=bfs0->JA[kj];
      kel=etree[jel];// ancestor, this stays unchanged
      if(kel<0){
	fprintf(stdout,"\n%%splitting element=%d",jel);
      } else {
	je=locate0(jel,elneib[kel], c2s->nf);
	ke=locate0(kel,elneib[jel], c2s->nf);	  
	/* 
	   in kel, we have the face je; in jel we have the face ke. we
	   want to make ke in jel same as je in kel; also the opposite
	   face needs to be reordered. 
	*/
	for(i=0;i<nvface;i++){
	  facei[i]=c2s->faces[ke*nvface+i];
	  facei[i]=g0->mnodes[jel*(nvcube+1)+facei[i]];
	  mnodes[i]=c2s->faces[je*nvface+i];
	  mnodes[i]=g0->mnodes[kel*(nvcube+1)+mnodes[i]];
	}	
	k1=aresamep(mnodes,facei,nvface,p);
	if(!k1){
	  fprintf(stderr,"\nERROR: faces must have same vertices in: %s\n",__FUNCTION__);
	  exit(127);
	} else{
	  /* fprintf(stdout,"\nface(%d)in el %d",ke,jel); */
	  /* fprintf(stdout,"\nOLD face(%d)in el %d",je,kel); */
	  /* print_full_mat_int(1,nvface,facei,"faceswp"); */
	  /* print_full_mat_int(1,nvface,p,"p");	 */
	  for(i=0;i<nvface;i++){ 
	    facei[c2s->faces[ke*nvface+p[i]]]=c2s->faces[ke*nvface+i];
	    keswp=(ke+c2s->n)%c2s->nf;
	    facei[c2s->faces[keswp*nvface+p[i]]]=c2s->faces[keswp*nvface+i]; 
	  }
	  /* fprintf(stdout,"\nface(%d)in el %d",ke,jel); */
	  /* print_full_mat_int(1,nvcube,facei,"facei"); */
	  for(i=0;i<nvcube;i++)
	    el2v->JA[i]=g0->mnodes[jel*(nvcube+1)+facei[i]];
	  for(i=0;i<nvcube;i++)
	    g0->mnodes[jel*(nvcube+1)+i]=el2v->JA[i];
	}
      }
    }
  }
  print_full_mat_int(g0->nel,c2s->nvcube+1,g0->mnodes,"mel1");
  icsr_free(el2v);
  icsr_free(f2v);
  icsr_free(el2el);
  /*****************************************************/    
  fprintf(stdout,"\n"); 
  /*****************************************************/
  return bfs0;
}
/*******************************************************************/
scomplex **macro_split(input_grid *g0,cube2simp *c2s)
{
  /* 
     From an input grid loops over the macroelements and sets up the
     divisions in every dimension. First makes the array of edges with
     their division consistent (the input can be inconsistent) grids
     each having a single macroelement.
  */
  input_grid *g;
  INT i,j,j0,j1,kel,ke,pmem;
  INT nel0=g0->nel,nvcube=c2s->nvcube,nvface=c2s->nvface;
  scomplex **sc=malloc(g0->nel*sizeof(scomplex *));
  pmem=2*g0->nv;
  if(pmem<2*g0->ne) pmem=2*g0->ne;
  if(pmem<2*g0->nel) pmem=2*g0->nel;
  if(pmem<2*g0->nf) pmem=2*g0->nf;
  INT *p=calloc(pmem,sizeof(INT));
  //  print_full_mat_int(g0->nel,(c2s->nvcube+1),g0->mnodes,"ex");
  ilexsort(g0->nf, (c2s->nvface+1),g0->mfaces,p);
  ilexsort(g0->nel,(c2s->nvcube+1),g0->mnodes,p);
  /*-------------------------------------------------------------------*/
  //  INT *efound=calloc(c2s->ne*(g0->nel),sizeof(INT));
  INT **nd=calloc(g0->nel,sizeof(INT *));
  for(i=0;i<g0->nel;i++){
    nd[i]=calloc(c2s->n,sizeof(INT)); /* to hold the number of
 					 divisions in every coordinate
					 direction */
  }
  for(kel=0;kel<g0->nel;kel++)
    for(i=0;i<c2s->nf;i++){
      nd[kel][i]=-1;
    }
  g=malloc(1*sizeof(input_grid));
  /**/
  g->title=g0->title;
  g->dgrid=g0->dgrid;
  g->fgrid=g0->fgrid;
  g->dvtu=g0->dvtu;
  g->fvtu=g0->fvtu;
  g->print_level=g0->print_level;
  g->ref_type=g0->ref_type;
  g->nref=g0->nref;
  g->err_stop=g0->err_stop;
  /**/
  g->dim=c2s->n;
  g->ncsys=g0->ncsys;
  g->nv=c2s->nvcube;
  g->nf=c2s->nf;
  g->ne=c2s->ne;
  g->nel=1;
  input_grid_arrays(g);
  /* reassign this as these are the same as g0 */
  free(g->systypes);   g->systypes=g0->systypes;
  free(g->syslabels);   g->syslabels=g0->syslabels;
  free(g->ox); g->ox=g0->ox;
  INT chng=1,iter=0,maxiter=1024;  
  //  INT je,kj,k2,iel2v,jel2v,k1,kface,kbnd,found;
  set_edges(g0,c2s);
  while(chng && (iter<maxiter)){
    iter++;
    // make the divisions in g0->seg consistent;
    chng=set_ndiv_edges(g,g0,c2s,nd,iter);
  }
  /*-------------------------------------------------------------------*/
  // here we should end the function which sets all divisions
  /*-------------------------------------------------------------------*/  
  // form macroelement neighboring list. Use transposition; the sparse
  // matrices here are only used to compute the number of connected
  // components in the domain or on the boundary.  
  INT **elneib=calloc(g0->nel,sizeof(INT *));/* element neighboring
					       list where the
					       position of the
					       neighbor is the same
					       as the position of the
					       face in c2s->faces
					       shared by the two
					       el.
					    */
  INT **el2fnum=calloc(g0->nel,sizeof(INT *)); /* for every element
						 this gives the local
						 to global face number
						 map;
					      */
  for(i=0;i<g0->nel;i++){
    elneib[i]=calloc(c2s->nf,sizeof(INT)); /* to hold the neighbors */
    el2fnum[i]=calloc(c2s->nf,sizeof(INT)); /* to hold the face
					       numbers for every
					       element */
  }
  ivector *bcodesf=malloc(1*sizeof(ivector));
  ivector *isbface=malloc(1*sizeof(ivector));
  ivector *etree=malloc(1*sizeof(ivector));
  iCSRmat *bfs0=set_mmesh(g0,c2s,etree,		\
			  elneib, el2fnum,	\
			  isbface, bcodesf,p);
  /* fprintf(stdout,"\nbfs0=["); */
  /* icsr_print_matlab_val(stdout,bfs0); */
  /* fprintf(stdout,"];"); */
  /* fprintf(stdout,"\nbfs=sparse(bfs0(:,1),bfs0(:,2),bfs0(:,3));\n"); */
  /* print_full_mat_int(1,isbface->row,isbface->val,"isbface"); */
  /* fprintf(stdout,"\n*****   nf=%d (%d)******* \n",g0->nf,isbface->row); */
  /***********************************************************************/    
  /* set the divisions on every edge now; since they are consistent we
   have: */
  if(set_ndiv_edges(g,g0,c2s,nd,0)) {
    fprintf(stderr,"\n\n***ERR in %s: the divisions of the edges cannod be inconsistent during second call of set_ndiv_edges()\n\n",__FUNCTION__);
    exit(4);
  }
  /* for(kel=0;kel<g0->nel;kel++){ */
  /*   print_full_mat_int(1,c2s->nf,elneib[kel],"neib"); */
  /* } */
  for(kel=0;kel<g0->nel;kel++){
    for(i=0;i<c2s->n;i++){
      if(nd[kel][i]<=0) nd[kel][i]=1;
    }
    print_full_mat_int(1,c2s->n,nd[kel],"ndnd");
  }
  // use bfs to split elements:
  INT lvl,intype=1,kj,je,jel;  
  INT nsall,nvall;
  nsall=0;
  nvall=0;
  for(lvl=0;lvl<bfs0->row;lvl++){    
    j0=bfs0->IA[lvl];
    j1=bfs0->IA[lvl+1];
    for(kj=j0;kj<j1;kj++){
      jel=bfs0->JA[kj];
      //      fprintf(stdout,"\n **** %d",jel);
      //      print_full_matrix_int(1,c2s->nf,elneib[jel],"neib");
      //      print_full_mat_int(1,c2s->n,nd[jel],"ndnd");
      kel=etree->val[jel];// ancestor, this stays unchanged
      memcpy(g->mnodes,(g0->mnodes+jel*(nvcube+1)),(nvcube+1)*sizeof(INT));
      for(i=0;i<c2s->nvcube;i++){
	j=g->mnodes[i];// vertex number (global)
	g->csysv[i]=g0->csysv[j]; 
	g->labelsv[i]=g0->csysv[j];
	memcpy((g->xv+i*g->dim),(g0->xv+j*g0->dim),g->dim*sizeof(REAL));
      }
      //      print_full_mat(g->nv,g->dim,g->xv,"xv{1}"); fflush(stdout);
      if(kel<0){
	sc[jel]=umesh(g->dim,nd[jel],c2s,intype);
	//        print_full_mat_int(1,c2s->n,nd[jel],"ndnd");
      } else {
	je=locate0(jel,elneib[kel], c2s->nf);
	ke=locate0(kel,elneib[jel], c2s->nf);
	//	fprintf(stdout,"\nface(%d) in a divided el(%d) is also face(%d) in the current el(%d)",je,kel,ke,jel);
	sc[jel]=umesh(g->dim,nd[jel],c2s,intype);
      }
      /*copy vertices and coord systems*/
      nsall+=sc[jel]->ns;
      nvall+=sc[jel]->nv;
      fprintf(stdout,"\n%%Generated a uniform mesh in dim=%d; vertices: %d, simplexes %d",g->dim,sc[jel]->nv,sc[jel]->ns);      
      fprintf(stdout,"\n%%Mapping back to the macroelement...\n");
      map2mac(sc[jel],c2s,g);
      //      haz_scomplex_print(sc[jel],0,"HAHA");
    }
  }
  /**/
  //  for(i=0;i<g0->nel;i++){
  //    
  //  }
  free(p);
  icsr_free(bfs0);
  ivec_free(etree);
  ivec_free(bcodesf);
  ivec_free(isbface);
  for(i=0;i<g0->nel;i++){
    free(nd[i]);
    free(elneib[i]);
    free(el2fnum[i]);
  }
  free(nd);
  free(elneib);
  free(el2fnum);
  return sc;  
}
/****************************************************************/
/*   sc->nbig=n; sc->n=n; */
/*   INT n1=sc->n+1,i,j,in1; */
/*   sc->factorial=1.; */
/*   sc->level=0; */
/*   for (j=2;j<n1;j++) sc->factorial *= ((REAL )j); */
/*   // fprintf(stdout,"\nIMPORTANT: NS=%d (%d!)=%f",ns,n,sc->factorial); */
/*   sc->marked=(INT *) calloc(ns,sizeof(INT)); */
/*   sc->gen=(INT *) calloc(ns,sizeof(INT)); */
/*   sc->nbr=(INT *) calloc(ns*n1,sizeof(INT)); */
/*   sc->parent=(INT *)calloc(ns,sizeof(INT)); */
/*   sc->child0=(INT *)calloc(ns,sizeof(INT)); */
/*   sc->childn=(INT *)calloc(ns,sizeof(INT)); */
/*   sc->nodes=(INT *)calloc(ns*n1,sizeof(INT)); */
/*   sc->bndry=(INT *)calloc(nv,sizeof(INT)); */
/*   sc->csys=(INT *)calloc(nv,sizeof(INT));/\* coord sys: 1 is polar, 2 */
/* 					    is cyl and so on *\/ */
/*   sc->flags=(INT *)calloc(ns,sizeof(INT)); */
/*   sc->x=(REAL *)calloc(nv*n,sizeof(REAL)); */
/*   sc->vols=(REAL *)calloc(ns,sizeof(REAL)); */
/*   sc->fval=(REAL *)calloc(nv,sizeof(REAL)); // function values at every vertex; not used in general; */
/*   for (i = 0;i<ns;i++) { */
/*     sc->marked[i] = FALSE; // because first is used for something else.  */
/*     sc->gen[i] = 0; */
/*     sc->parent[i]=-1; */
/*     sc->child0[i]=-1; */
/*     sc->childn[i]=-1; */
/*     sc->flags[i]=-1; */
/*     in1=i*n1; */
/*     for(j=0;j<n1;j++){ */
/*       sc->nodes[in1+j]=-1; */
/*       sc->nbr[in1+j]=-1; */
/*     } */
/*   } */
/*   for (i = 0;i<nv;i++) { */
/*     sc->bndry[i]=0; */
/*     sc->csys[i]=0; */
/*     sc->fval[i]=0.; */
/*   } */
/*   sc->nv=nv; */
/*   sc->ns=ns; */
/*   return sc; */
/* } */
/**/
INT main(INT argc, char **argv)
{
  INT i=-1;
  input_grid *g=parse_input_grid("grid.input");
  INT dim=g->dim;
  cube2simp *c2s=cube2simplex(dim);
  scomplex **sc=macro_split(g,c2s);
  /*this can be used to generate grids in a different way, but not now:*/
  /* if(intype<-1){ */
  /*   for(i=0;i<dim;i++) ndd[i]=1; */
  /*   sc[0]=umesh(dim,ndd,c2s,intype); */
  /*   unirefine(nd,sc);     */
  /* }else{ */
  /*   sc[0]=umesh(dim,nd,c2s,intype); */
  /* } */
  /*GENERATE UNIFORM GRID: nodes in each direction: ND*/
  /*OUTPUT */
  if(dim==2||dim==3) {
    fprintf(stdout,"Writing vtk file...\n");
    vtkw("newmesh.vtu",sc[0],0,0,1.);
  }
  /*FREE*/
  /* for(i=0;i<g->nel;i++){ */
  /*       haz_scomplex_free(sc[i]); */
  /* } */
  free(sc);
  input_grid_free(g); 
  fprintf(stdout,"\nDone.\n");
  cube2simp_free(c2s);
  return 0;
}
/*********************EOF**********************************/
