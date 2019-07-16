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
iCSRmat *bfs01(INT nblk,INT *iblk, INT *jblk, iCSRmat *a, INT *et)
{
  /* 
     bfs for the graph with possibly multiple connected components;
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
  fprintf(stdout,"\nlvlend=%d; ib=%d",lvlend,ib);
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
  /* exit(77); */
  return bfs;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx*/
INT locate0(INT needle, INT *haystack, INT n)
{
  /* 
     finds an element in an array. on output gives the index in the
     array where the element is found
   */
  INT i;
  for (i=0;i<n;i++)
    if(needle==haystack[i])
      return i;
  return -1;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx*/
INT *set_input_grid1(input_grid *g,cube2simp *c2s)
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
    j=g->seg[3*i+1]-g->seg[3*i]-1;
    p[k]=g->seg[3*i+2];
    k++;
    //  fprintf(stdout,"\n[%d,%d]:div=%d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
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
INT set_ndiv_edges(input_grid *g,		\
		   input_grid *g0,		\
		   cube2simp *c2s,		\
		   INT *efound,			\
		   INT *ffound,			\
		   INT **nd,
		   const INT iter)
{
  /* 
     For a given global input grid g0 creates local input grids for
     every macroelement and computes the divisions for it. It is used
     iteratively in macro_split to se the correct divisions for every
     macroelement. efound should be of length c2s->ne +g0->ne and
     ffound should be of length c2s->nf +g0->nf.  The input_grid g0
     should be all set, and the input_grid g should have all its
     scalar values set.
     nd is the array with the divisions, it must be g0->nel by c2s->n. 
  */
  INT kel0,i,j0,j1,swp,kel,ke,k0,k1,ndiv;
  INT nel0=g0->nel,nvcube=c2s->nvcube;
  INT *e0found=efound + c2s->ne;
  INT *f0found=ffound + c2s->nf;
  /*foe easier reference*/
  for(i=0;i<g0->ne;i++)e0found[i]=-1;
  for(i=0;i<g0->nf;i++)f0found[i]=-1;
  // make all divisions > 0
  for(ke=0;ke<g0->ne;ke++){
    ndiv=abs(g0->seg[3*ke+2]);
    if(ndiv==0)ndiv=1;
    g0->seg[3*ke+2]=ndiv;
  }
  for(ke=0;ke<g0->ne;ke++)      
    e0found[ke]=g0->seg[3*ke+2];
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
      //	  fprintf(stdout,"\nElement:%d, edge=(%d,%d);",kel,j0,j1);
    }
      //    input_grid_print(g);
    nd[kel]=set_input_grid1(g,c2s);    
    for(i=0;i<g->ne;i++){      
      if(efound[i]<0) continue;
      ke=efound[i];
      g0->seg[3*ke+2]=g->seg[3*i+2];
    }
    //    print_full_mat_int(1,c2s->n,nd[kel],"ndnd");
  }
  INT chng=0;
  for(ke=0;ke<g0->ne;ke++){
    k1=abs(e0found[ke]-g0->seg[3*ke+2]);
    if(k1>chng)chng=k1;
  }
  //  fprintf(stderr,"\nchng=%d",chng);
  return chng;
}
/*******************************************************/
scomplex *macro_split(input_grid *g0,cube2simp *c2s)
{
  /* 
     From an input grid loops over the macroelements and sets up the
     divisions in every dimension. First makes the array of edges with
     their division consistent (the input can be inconsistent) grids
     each having a single macroelement.
  */
  input_grid *g;
  scomplex *sc;
  INT i,j0,j1,kel,jel,ke,pmem;
  INT nel0=g0->nel,nvcube=c2s->nvcube,nvface=c2s->nvface;
  pmem=2*g0->nv;
  if(pmem<2*g0->ne) pmem=2*g0->ne;
  if(pmem<2*g0->nel) pmem=2*g0->nel;
  if(pmem<2*g0->nf) pmem=2*g0->nf;
  INT *p=calloc(pmem,sizeof(INT));
  ilexsort(g0->nel,(c2s->nvcube+1),g0->mnodes,p);
  //  print_full_mat_int(g0->nel,(c2s->nvcube+1),g0->mnodes,"ex");
  ilexsort(g0->nf, (c2s->nvface+1),g0->mfaces,p);
  /*-------------------------------------------------------------------*/
  INT *efound=calloc(c2s->ne+g0->ne,sizeof(INT));
  INT *ffound=calloc(c2s->nf+g0->nf,sizeof(INT));
  /*-------------------------------------------------------------------*/
  INT **nd=calloc(g0->nel,sizeof(INT *));
  INT **elneib=calloc(g0->nel,sizeof(INT *));// element neighboring
					    // list where the position
					    // of the neighbor is the
					    // same as the position of
					    // the face in c2s->faces
					    // shared by the two el.
  INT **el2fnum=calloc(g0->nel,sizeof(INT *)); // for every element
					       // this gives the local
					       // to grobal face
					       // number map;
  for(i=0;i<g0->nel;i++){
    nd[i]=calloc(c2s->n,sizeof(INT)); /* to hold the number of
 					 divisions in every coordinate
					 direction */
    elneib[i]=calloc(c2s->nf,sizeof(INT)); /* to hold the neighbors */
    el2fnum[i]=calloc(c2s->nf,sizeof(INT)); /* to hold the neighbors */
 }
  for(kel=0;kel<g0->nel;kel++)
    for(i=0;i<c2s->nf;i++){
      elneib[kel][i]=-1;
      el2fnum[kel][i]=-1;
    }
  /*-------------------------------------------------------------------*/
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
  INT je,kj,k2,iel2v,jel2v,k1,kface,kbnd,found;
  while(chng&&(iter<maxiter)){
    iter++;
    // make the divisions in g0->seg consistent;
    chng=set_ndiv_edges(g,g0,c2s,efound,ffound,nd,iter);
  }
  /* set the divisions on every edge now; since they are consistent we
   have: */
  if(set_ndiv_edges(g,g0,c2s,efound,ffound,nd,0)) {
    fprintf(stderr,"\n\n***ERR in %s: the divisions of the edges cannod be inconsistent during second call of set_ndiv_edges()\n\n",__FUNCTION__);
    exit(4);
  }
  for(kel=0;kel<g0->nel;kel++)
    for(i=0;i<c2s->n;i++)
      if(nd[kel][i]<=0) nd[kel][i]=1;  
  // form macroelement neighboring list. Use transposition; the sparse
  // matrices here are only used to compute the number of connected
  // components in the domain or on the boundary.
  iCSRmat *el2v=malloc(1*sizeof(iCSRmat));
  el2v[0]=icsr_create(g0->nel,g0->nv,nvcube*g0->nel);
  el2v->IA[0]=0;
  for(kel=0;kel<nel0;kel++){
    memcpy((el2v->JA+el2v->IA[kel]),(g0->mnodes+kel*(nvcube+1)),nvcube*sizeof(INT));
    el2v->IA[kel+1]=el2v->IA[kel]+nvcube;
  }
  fprintf(stdout,"\n\n YYY************ %d %d *************",el2v->nnz,el2v->IA[nel0]);fflush(stdout);
  for(i=0;i<el2v->IA[nel0];i++)
    el2v->val[i]=1;
  iCSRmat *v2el=malloc(1*sizeof(iCSRmat));
  icsr_trans(el2v,v2el);	
  /* fprintf(stdout,"\nv2el=["); */
  /* icsr_print_matlab_val(stdout,v2el); */
  /* fprintf(stdout,"];\n"); */
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
  //  print_full_mat_int(1,el2el->nnz,el2el->JA,"JA");
  /* fprintf(stdout,"\nel2el0=["); */
  /* icsr_print_matlab_val(stdout,el2el); */
  /* fprintf(stdout,"];el2el=sparse(el2el0(:,1),el2el0(:,2),el2el0(:,3),%d,%d);\n\n",g0->nel,g0->nel); */
  /* fprintf(stdout,"\n\n XXX************ %d %d *************",el2el->nnz,nel0);fflush(stdout); */
  // find the connected components:
  INT *iblk=calloc((nel0+1),sizeof(INT));
  INT *jblk=calloc((nel0),sizeof(INT));
  INT nblkdom,nblkbnd; // number of connected domains and this boundaries;
  // boundaries
  dfs00_(&nel0,el2el->IA, el2el->JA,&nblkdom,iblk,jblk);
  iblk=realloc(iblk,(nblkdom+1)*sizeof(INT));
  fprintf(stdout,"\nDFS(domains): %d connected components",nblkdom);  
  for(ke=0;ke<nblkdom;ke++){
    fprintf(stdout,"\ncc=%d:  v=",ke);
    j0=iblk[ke];
    j1=iblk[ke+1];
    for(kj=j0;kj<j1;kj++){
      fprintf(stdout,"%d ",jblk[kj]);
    }
  }
  fprintf(stdout,"\n");
  /* fprintf(stdout,"\nel2elx=["); */
  /* icsr_print_matlab_val(stdout,el2el); */
  /* fprintf(stdout,"];el2el=sparse(el2elx(:,1),el2elx(:,2),el2elx(:,3),%d,%d);\n\n",g0->nel,g0->nel); */
  /* fprintf(stdout,"\n\n ************ nonzeroes in el2el=%d (%d) *************",el2el->nnz,el2el->IA[nel0]);fflush(stdout); */
  icsr_nodiag(el2el);
  //  print_full_mat_int(1,el2el->row+1,el2el->IA,"iaeee");
  //  print_full_mat_int(1,el2el->nnz,el2el->JA,"jaeee");
  INT *etree=calloc((el2el->row+1),sizeof(INT));
  iCSRmat *bfs0=bfs01(nblkdom,iblk,jblk,el2el,etree);  
  fprintf(stdout,"\nbfs0=[");
  icsr_print_matlab_val(stdout,bfs0);
  fprintf(stdout,"];bfs=sparse(bfs0(:,1),bfs0(:,2),bfs0(:,3),%d,%d);\n\n",bfs0->row,bfs0->col);
  fprintf(stdout,"\n\n ************ nonzeroes in el2el=%d (%d) *************",el2el->nnz,el2el->IA[nel0]);fflush(stdout);
  /*FACES******************************************************/
  INT nfaceall=g0->nel*c2s->nf-(el2el->nnz/2);
  INT nfacei=(INT )(el2el->nnz/2);
  INT nfaceb=nfaceall-nfacei;
  // face to vertex:
  iCSRmat *f2v=malloc(sizeof(iCSRmat));
  f2v[0]=icsr_create(nfaceall,g0->nv,nvface*nfaceall);
  INT *bcodesf=calloc(2*nfaceall,sizeof(INT));
  INT *isbface=bcodesf+nfaceall;
  f2v->IA[0]=0;
  for(i=0;i<nfaceall;i++){
    f2v->IA[i+1]=f2v->IA[i]+nvface;
    isbface[i]=0;
  }
  /*  print_full_mat_int(1,nfaceall+1,f2v->IA,"ia");   */
  /* fprintf(stdout,"\n *** int_faces=%d; bnd_faces=%d\n",nfacei,nfaceb);   */
  /* fprintf(stdout,"\n *** all_faces: %d ?= %d; (nnznew?=nnz):%d?=%d\n",nfaceall,f2v->row,f2v->nnz,f2v->IA[f2v->row]);   */
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
  fprintf(stdout,"\nkface=%d  ? = ? %d\n",kface,f2v->nnz);
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
  for (i=0;i<nel0;i++){
    print_full_mat_int(1,c2s->nf,elneib[i],"elneib");
    //    print_full_mat_int(1,c2s->n,nd[i],"nd");fflush(stdout);
  }
  /***********************************************************************/    
  iCSRmat *v2f=malloc(1*sizeof(iCSRmat));
  icsr_trans(f2v,v2f);
  /*******************************************************************/    
  /* fprintf(stdout,"\nf2v0=["); */
  /* icsr_print_matlab_val(stdout,f2v); */
  /* fprintf(stdout,"];"); */
  /* fprintf(stdout,"\nf2v=sparse(f2v0(:,1),f2v0(:,2),f2v0(:,3));\n"); */
  /* print_full_mat_int(1,nfaceall,isbface,"isbface"); */
  /* fprintf(stdout,"\n*****   nf=%d (%d)******* \n",g0->nf,nfaceall);  */
  /***********************************************************************/    
  iCSRmat *el2f=malloc(1*sizeof(iCSRmat));
  icsr_mxm(el2v,v2f,el2f);
   /* now remove all entries that are not 2 to obtain the el2f
      matrix; */
  el2f->nnz=el2f->IA[0];
  for(kel=0;kel<el2f->row;kel++){    
    j0=el2f->IA[kel];
    j1=el2f->IA[kel+1];
    el2f->IA[kel]=el2f->nnz;
    for(ke=j0;ke<j1;ke++){
      if(el2f->val[ke]!=(nvface)) continue;
      el2f->JA[el2f->nnz]=el2f->JA[ke];
      el2f->val[el2f->nnz]=el2f->val[ke];
      el2f->nnz++;
    }
  }    
  el2f->IA[el2f->row]=el2f->nnz;
  el2f->JA=realloc(el2f->JA,el2f->nnz*sizeof(INT));
  el2f->val=realloc(el2f->val,el2f->nnz*sizeof(INT));
  /* fprintf(stdout,"\nel2f=["); */
  /* icsr_print_matlab_val(stdout,el2f); */
  /* fprintf(stdout,"];"); */
  /* fprintf(stdout,"\nel2f=sparse(el2f(:,1),el2f(:,2),el2f(:,3));\n"); */
  /***********************************************************************/    
  iCSRmat *f2f=malloc(1*sizeof(iCSRmat));
  /*******************************************************************/    
  icsr_mxm(f2v,v2f,f2f);
  // now remove all entries in f2f that are not 1;
  f2f->nnz=el2f->IA[0];
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
  icsr_free(el2f);
  /*connected comps on the boundary*/
  iblk=realloc(iblk,(nfaceall+1)*sizeof(INT));
  jblk=realloc(jblk,(nfaceall)*sizeof(INT));
  dfs00_(&nfaceall,f2f->IA, f2f->JA,&nblkbnd,iblk,jblk);
  fprintf(stdout,"\nDFS(boundaries): %d connected components",nblkbnd-nfacei);
  //icsr_nodiag(f2f);
  icsr_free(f2f);
  free(iblk); 
  free(jblk);   
  /* fprintf(stdout,"\nf2fb=["); */
  /* icsr_print_matlab_val(stdout,f2f); */
  /* fprintf(stdout,"];"); */
  /* fprintf(stdout,"\nf2f=sparse(f2fb(:,1),f2fb(:,2),f2fb(:,3));\n"); */
  //
  /*now use the bfs*/
  // we first split the root element;
  INT lvl,keok,swp,keswp;
  for(lvl=0;lvl<bfs0->row;lvl++){
    j0=bfs0->IA[lvl];
    j1=bfs0->IA[lvl+1];
    for(kj=j0;kj<j1;kj++){
      jel=bfs0->JA[kj];
      kel=etree[jel];// ancestor, this stays unchanged
      if(kel<0){
	fprintf(stdout,"\nsplitting element=%d",jel);
      } else {
	je=locate0(jel,elneib[kel], c2s->nf);
	ke=locate0(kel,elneib[jel], c2s->nf);	  
	fprintf(stdout,"\nel=%d on face %d in el %d",jel,je,kel);
	keok=(je+c2s->n)%c2s->nf;
	if(keok!=ke){
	  fprintf(stdout,"\nel=%d on face %d (should be %d) in el *** %d ***",kel,ke,keok,jel);
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
	    // use g->mnodes as work space here;
	    g->mnodes[i]=c2s->faces[keok*nvface+i];
	    keswp=(keok+c2s->n)%c2s->nf;
	    g->mnodes[i+nvface]=c2s->faces[keswp*nvface+i];	    
	  }
	  for(i=0;i<nvcube;i++)
	    el2v->JA[i]=g0->mnodes[jel*(nvcube+1)+facei[i]];	  
	  for(i=0;i<nvcube;i++)
	    g0->mnodes[jel*(nvcube+1)+g->mnodes[i]]=el2v->JA[i];
	  fprintf(stdout,"\nEL=%d: ",jel);
	  for(i=0;i<nvcube;i++){
	    k1=g0->mnodes[jel*(nvcube+1)+i];
	    fprintf(stdout,"%d ",k1);
	  }
	}
      }
    }
  }
  /**** FINAL REORDER ***/
  for(lvl=0;lvl<bfs0->row;lvl++){
    j0=bfs0->IA[lvl];
    j1=bfs0->IA[lvl+1];
    for(kj=j0;kj<j1;kj++){
      jel=bfs0->JA[kj];
      kel=etree[jel];// ancestor, this stays unchanged
      if(kel<0){
	fprintf(stdout,"\nsplitting element=%d",jel);
      } else {
	je=locate0(jel,elneib[kel], c2s->nf);
	ke=locate0(kel,elneib[jel], c2s->nf);	  
	fprintf(stdout,"\nel=%d on face %d in el %d",jel,je,kel);
	keok=(je+c2s->n)%c2s->nf;
	if(keok!=ke){
	  fprintf(stdout,"\nel=%d on face %d (should be %d) in el *** %d ***",kel,ke,keok,jel);
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
	    // use g->mnodes as work space here;
	    g->mnodes[i]=c2s->faces[keok*nvface+i];
	    keswp=(keok+c2s->n)%c2s->nf;
	    g->mnodes[i+nvface]=c2s->faces[keswp*nvface+i];	    
	  }
	  for(i=0;i<nvcube;i++)
	    el2v->JA[i]=g0->mnodes[jel*(nvcube+1)+facei[i]];	  
	  for(i=0;i<nvcube;i++)
	    g0->mnodes[jel*(nvcube+1)+g->mnodes[i]]=el2v->JA[i];
	  fprintf(stdout,"\nEL=%d: ",jel);
	  for(i=0;i<nvcube;i++){
	    k1=g0->mnodes[jel*(nvcube+1)+i];
	    fprintf(stdout,"%d ",k1);
	  }
	}
      }
    }
  }
  print_full_mat_int(g0->nel,c2s->nvcube+1,g0->mnodes,"mel");
  icsr_free(bfs0);
  free(etree);
  icsr_free(el2v);
  icsr_free(f2v);
  icsr_free(el2el);
  /*****************************************************/    
  fprintf(stdout,"\n"); 
  /*****************************************************/    
  for(kel=0;kel<nel0;kel++){
    memcpy(g->mnodes,(g0->mnodes+kel*(nvcube+1)),(nvcube+1)*sizeof(INT));
    for(i=0;i<g->nf;i++){
      for(ke=0;ke<nvface;ke++){
	g->mfaces[i*(nvface+1)+ke]=c2s->faces[i*nvface+ke];
      }
      g->mfaces[i*(nvface+1)+nvface]=-1;
    }
  }
  /*copy vertices and coord systems*/
  for(kel=0;kel<nel0;kel++){
    memcpy(g->mnodes,(g0->mnodes+kel*(nvcube+1)),(nvcube+1)*sizeof(INT));
    for(i=0;i<g->nv;i++){
      j0=g->mnodes[i];// vertex number (global)
      g->csysv[i]=g0->csysv[j0]; 
      g->labelsv[i]=g0->csysv[j0];
      memcpy((g->xv+i*g->dim),(g0->xv+j0*g0->dim),g->dim*sizeof(REAL));
    }
    //    print_full_mat(g->nv,g->dim,g->xv,"xv{1}"); fflush(stdout);
  }
  free(efound);
  free(ffound);
  for(i=0;i<g0->nel;i++)
    free(nd[i]);
  free(nd);
  exit(33);
  return sc;  
}
/*********************************************************************/
void map2mac(scomplex *sc,cube2simp *c2s, input_grid *g)
{
  /* 
     maps a uniform grid from the n-dimensional cube to a hexagonal 
     macroelement given by its coordinates xmacro[1:nvcube*dim]
     xmac[nvcube][dim]
  */
  INT i,j,k1,k2,k1c,k2c,kf,dim=sc->n;
  INT ksys;
  REAL *xmac=g->xv;  
  REAL *xhat = (REAL *)calloc(dim,sizeof(REAL));
  REAL *xemac=(REAL *)calloc(c2s->ne*dim,sizeof(REAL));
  // convert midpoints from polar to cartesian.
  for(i=0;i<c2s->ne;i++){
    k1=c2s->edges[2*i];
    k2=c2s->edges[2*i+1];
    k1c=g->systypes[g->csysv[k1]];
    k2c=g->systypes[g->csysv[k2]];
    fprintf(stdout,"\nverts=(%d,%d); coord_sys=(%d,%d)",k1,k2,k1c,k2c);fflush(stdout);
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
    k2c=g->systypes[g->csysv[k2]];
    //skip all polar mid points
    if(g->csysv[k1]==g->csysv[k2] && k1c==1) continue;
    fprintf(stdout,"\ncart:verts=(%d,%d); coord_sys=(%d,%d)",k1,k2,k1c,k2c);fflush(stdout);
    for(j=0;j<dim;j++) {
      xemac[i*dim+j]=0.5*(xmac[k1*dim+j]+xmac[k2*dim+j]);
    }
  }
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
INT main(INT argc, char **argv)
{
  //  INT j=-1;
  input_grid *g=parse_input_grid("grid.input");
  INT dim=g->dim;
  cube2simp *c2s=cube2simplex(dim);
  /*------------------------------------------------------*/
  scomplex *sc;
  INT intype=0;
  /*------------------------------------------------------*/
  macro_split(g,c2s);
  INT *nd=set_input_grid1(g,c2s);  
  /*this can be used to generate grids in a different way, but not now:*/
  /* if(intype<-1){ */
  /*   for(i=0;i<dim;i++) ndd[i]=1; */
  /*   sc=umesh(dim,ndd,c2s,intype); */
  /*   unirefine(nd,sc);     */
  /* }else{ */
  /*   sc=umesh(dim,nd,c2s,intype); */
  /* } */
  /*GENERATE UNIFORM GRID: nodes in each direction: ND*/
  sc=umesh(dim,nd,c2s,intype);
  fprintf(stdout,"\nGenerated a uniform mesh in dim=%d; vertices: %d, simplexes %d",dim,sc->nv,sc->ns);
  if(nd) free(nd);
  fprintf(stdout,"\nedges=%d",c2s->ne);
  //  INT k1,k2,j1,j2,l1,l2;
  //KEEP KEEP  gcomplex *macs=form_macro(g);
  fprintf(stdout,"\nMapping back to the macroelement...\n");
  map2mac(sc,c2s,g);
  input_grid_free(g); 
  fprintf(stdout,"\nDone.\n");
  //  haz_scomplex_print(sc,0,"HAHA");
  if(dim==2||dim==3) {
    fprintf(stdout,"Writing vtk file...\n");
    vtkw("newmesh.vtu",sc,0,0,1.);
  }
  cube2simp_free(c2s);
  haz_scomplex_free(sc);
  return 0;
}
/*********************EOF**********************************/
