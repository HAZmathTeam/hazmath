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
/************************************************************************/
//void abfstree(INT it, scomplex *sc,INT *wrk);
//void scfinalize(scomplex *sc);
/***********************************************************************/
REAL interp4(cube2simp *c2s, REAL *u, REAL *xhat)
{
  /*INTerpolate d-linearly in d dimensions on the UNIT cube */
  INT dim=c2s->n,k,i,kdim;
  REAL phik;
  REAL s=0.;
  for(k = 0;k<c2s->nvcube;k++){
    kdim=k*dim;
    phik=1e0;
    for(i = 0;i<dim;++i){      
      if(!c2s->bits[kdim+i])
	phik*=(1e0-xhat[i]);
      else	
	phik*=xhat[i];
    }
    s+=u[k]*phik;
  }
  return s;
}
/************************************************************************/
REAL interp8(cube2simp *c2s, REAL *u, REAL *ue, REAL *xhat)
{
  /*INTerpolate quadratically in d dimensions on the UNIT cube */
  INT dim=c2s->n,k,i,k1,k2,kdim;  
  REAL phik,zmid,psimid;
  REAL s=0.;
  for(k = 0;k<c2s->nvcube;k++){
    kdim=k*dim;
    phik=1e0;
    for(i = 0;i<dim;++i){      
      if(!c2s->bits[kdim+i]){
	phik*=(1e0-xhat[i]);
	//	fprintf(stdout,"\nold:%d:vert:(1-x[%d])",k,i);
      } else {
	phik*=xhat[i];
	//	fprintf(stdout,"\nold:%d:vert:x[%d]",k,i);
      }
    }
    psimid=0.;zmid=0.;
    //    fprintf(stdout,"\nvertex=%d; normal=",k);
    for(i=0;i<dim;++i){
      if(!c2s->bits[kdim+i]){
	psimid-=xhat[i];
	zmid+=0.25e00;
	//	fprintf(stdout,"\n%d:vert:(-x[%d])",k,i);
      } else {
	psimid+=xhat[i];
	zmid-=0.75e00;
	//	fprintf(stdout,"\n%d:vert:x[%d]",k,i);
      }
    }
    //    fprintf(stdout,"\n%d:zmid=%5.1f",k,zmid);
    phik*=2e00*(zmid+psimid);
    s+=u[k]*phik;
  }
  REAL se,phie;
  se=0.;
  for(k=0;k<c2s->ne;k++){
    k2=c2s->edges[2*k];    
    k1=c2s->edges[2*k+1];
    //    fprintf(stdout,"\nmid=(%d,%d);",k1,k2);
    phie=1e0;
    for(i = 0;i<dim;++i){      
      if(!c2s->bits[k1*dim+i]){
	phie*=(1e0-xhat[i]);
	//	fprintf(stdout,"\n1:(1-x[%d])",i);
      } else {
	phie*=xhat[i];
	//	fprintf(stdout,"\n1:x[%d]",i);
      }
    }
    for(i = 0;i<dim;++i){      
      if(c2s->bits[k2*dim+i]==c2s->bits[k1*dim+i]) continue;
      if(!c2s->bits[k2*dim+i]){
	phie*=(1e0-xhat[i]);
	//	fprintf(stdout,"\n2:(1-x[%d])=%e",i,1.-xhat[i]);
      } else {
	phie*=xhat[i];
	//	fprintf(stdout,"\n2:x[%d]=%e",i,xhat[i]);
      }
    }
    phie*=4e00;
    //    fprintf(stdout,"\nedge=(%d,%d),coord=%e,phie=%e",k1,k2,ue[k],phie);
    se+=ue[k]*phie;
  }
  //  fprintf(stdout,"\n");
  //  fprintf(stdout,"\n***************** xhat=(%e,%e):%e",xhat[0],xhat[1],se);
  //  fprintf(stdout,"\ns=%e",s);
  //  fprintf(stdout,"\n");
  return (s+se);
}
/**********************************************************************/
static void coord_perm(SHORT type, INT n,void *x, size_t elsize)
{
  /*
    permutes coordinate arrays (or any array whose elements are of size): 
    if type=1: x[0],x[1]...x[n-1]--> x[1]...x[n-1],x[0] 
    if type=0(inverse): y[0],y[1]...y[n-1]--> y[n-1],y[0]...y[n-2] 
  */
  INT i;  
  void *xx0n=(void *)calloc(1,elsize*sizeof(void));
  if(type){
    memcpy(xx0n,x,elsize);
    for(i=0;i<n;i++)
      memcpy(x+i*elsize,x+(i+1)*elsize,elsize);
    memcpy(x+(n-1)*elsize,xx0n,elsize);
  } else {
    memcpy(xx0n,x+(n-1)*elsize,elsize);
    for(i=(n-1);i>=1;i--)
      memcpy(x+i*elsize,x+(i-1)*elsize,elsize);
    memcpy(x,xx0n,elsize);
  }
  if(xx0n)free(xx0n);
  return;
}
/************************************************************************/
void polar2cart(INT dim, REAL *px, REAL *cx)
{
  // polar is r, theta1,...theta[n-1]; cart is x[0]...x[n-1] px are
  // n-coordinates in polar coord system converts polar coordnates to
  // cartesian in d dimensions.  origin is set at 0,0,0 so if it is
  // different, translation needs to be done after return from here.
  INT i,j;
  if(dim==1) return;
  REAL rho = px[0];
  REAL cend=rho; 
  memset(cx,0,dim*sizeof(REAL));
  for(i=0;i<(dim-1);i++){
    cx[i]=rho*cos(px[i+1]);
    for(j=0;j<i;j++){
      cx[i]*=sin(px[j+1]);
      //      fprintf(stdout,"\niiiii=%d,jjjjj=%d",i,j+1);
    }
    //    print_full_mat(1,dim,cx,"c1");
    cend*=sin(px[i+1]);    
  }
  cx[dim-1]=cend;
  /* cx[0]=rho*cos(px[1]); */
  /* cx[1]=rho*sin(px[1])*cos(px[2]); */
  /* cx[2]=rho*sin(px[1])*sin(px[2]); */
  // the conversion above puts cx[n-1] first, so put it back at the end.  
  coord_perm(1,dim,cx,sizeof(REAL));
  return;
}
/************************************************************************/
INT cart2polar(INT dim, REAL *c,REAL *p)
{
  INT i,dimm1=dim-1;
  REAL rl,r;
  // first put c[n-1] first to agree with the polar ordering;
  coord_perm(0,dim,c,sizeof(REAL));
  r=0.;
  for(i=0;i<dim;i++){
    r+=(c[i]*c[i]);
    p[i]=0e0;
  }  
  if(fabs(r)<1e-14){
    for(i=1;i<dim;i++) p[i]=-1e20;
    return 1;
  }
  r=sqrt(r);
  rl=r;
  INT flag=1;
  for(i=1;i<dim;i++){
    p[i]=acos(c[i-1]/rl);
    if(fabs(sin(p[i]))<1e-14){
      flag=0;
      break;
    }    
    rl/=sin(p[i]);
  }
  if(flag) p[dimm1]=atan2(c[dimm1],c[dimm1-1]);
  p[0]=r;
  // permute c back;
  coord_perm(1,dim,c,sizeof(REAL));
  return 0;
}
/************************************************************************/
void unirefine(INT *nd,scomplex *sc)  
{
/* 
 * refine uniformly l levels, where 2^l>max_m nd[m] using the generic
 * algorithm for refinement.
 Works in the following way: first construct a grid with refinements up to
 2^{l} such that 2^{l}>max_m nd[m]. then remove all x such that
 x[k]>nd[k]*2^{-l} and then remap to the unit square. 
*/
  INT ndmax=-1,i=-1,j=-1;
  for(i=0;i<sc->n;i++)
    if(ndmax<nd[i]) ndmax=nd[i];
  fprintf(stdout,"\nmax split=%d",ndmax);
  REAL sref=log2((REAL )ndmax);
  if(sref-floor(sref)<1e-3)
    sref=floor(sref);
  else
    sref=floor(sref)+1.;
  INT ref_levels= sc->n*((INT )sref);
  fprintf(stdout,"\nlog2 of the max=%e, l=%d",log2((REAL )ndmax)+1,ref_levels);
  ref_levels=0;
  if(ref_levels<=0) return;
  INT nsold;//ns,nvold,level;
  if(!sc->level){
    /* form neighboring list; */
    find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
    //    haz_scomplex_print(sc,0,__FUNCTION__);  fflush(stdout);
    INT *wrk=calloc(5*(sc->n+2),sizeof(INT));
    /* construct bfs tree for the dual graph */
    abfstree(0,sc,wrk);
    //    haz_scomplex_print(sc,0,__FUNCTION__);fflush(stdout);
    //    exit(100);
    if(wrk) free(wrk);
  }
  //INT n=sc->n, n1=n+1,level=0;
  fprintf(stdout,"refine: ");
  while(sc->level < ref_levels && TRUE){
    nsold=sc->ns;
    //    nvold=sc->nv;
    for(j = 0;j < nsold;j++)sc->marked[j]=TRUE;
    for(j = 0;j < nsold;j++)
      if(sc->marked[j] && (sc->child0[j]<0||sc->childn[j]<0))
	haz_refine_simplex(sc, j, -1);
    /* new mesh */
    //    ns=sc->ns; 
    sc->level++;
    fprintf(stdout,"u%du",sc->level);//,nsold,ns,nv);
  }
  fprintf(stdout,"\n");
  scfinalize(sc);
  return;
}
/********************************LIBLIB */
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

static void binary0(cube2simp *c2s)
{
  // stores in an array the coordinates of the vertices of the unit  
  // cube in dimension (dim). Lexicographical ordering from 0,0,...0
  // to 1,1,...,1
  
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
  /* form all n+1 permutations in reverse order! it is unclear why in reverse order...*/
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
static unsigned INT bitdiff(const INT dim, unsigned INT *bits1,unsigned INT *bits2){
  /*
    returns the l1-norm of the difference two arrays of unsigned 
    integers.  this should be changed to have a void array as input.
  */
  INT j;
  unsigned INT numbits=0;
  for(j=0;j<dim;j++){
    numbits+=abs(bits1[j]-bits2[j]);
  }
  return numbits;
}
/************************************************************************/
void reverse(void *arr,INT length, size_t elsize)
{
  /* 
     permutes a void array whose elements are of size elsize 
     a[0],...a_[length-1]-->a[length-1],...a[0]. 
  */
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
/************************************************************************/
cube2simp *cube2simplex(INT dim)
{
  /*
    in dimension dim splits the cube in dim factorial dim-dimensional
    simplices. stores everything in a structure cube2simp. It also
    outputs all local permutations of vertices which can be used to
    create a criss-cross mesh. 
  */
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
  INT k1,kn1,k2,kn2,dim1=c2s->n+1,		\
    nvcube=c2s->nvcube;
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
  //  print_full_mat_int(c2s->nf,c2s->nvface,c2s->faces,"UCubef");
  return c2s;
}
scomplex *umesh(const INT dim, INT *nd, cube2simp *c2s,	\
		const INT face, const INT face_parent,	\
		const scomplex *sc_parent,		\
		INT *nd_parent,			\
		const INT intype)
{
  /* face is the face that matches the face_parent in the neighboring
     element.  uniform simplicial mesh of the unit cube in dimension
     dim.  dim is the dimension, nd is the number of grid points in
     each dimension.  ordering is lexicographically by
     name=(x[0],...,x[n]).  more than 3D is not fully tested xmacro[]
     are the coordinates of a domain isomorphic to the cube via a
     bilinear or "Q2" change of coordinates. output is a simplicial
     complex sc.

     if(intype == -2) use unirefine() function (from unigrid.c in src/amr)

     if(intype == -1)construct grid using diagonals pointing
     0-7(0...0)-->(1...1).

     if (intype>0) starting with intype the mesh is constructed like
     criss-cross grid. This works in 2D and 3D, and is unclear whether
     it works in d>3.
  */
  INT iz1;
  INT jperm,i,j,flag,kf,type;
  INT dim1 = dim+1;
  // m is dim+1 so that we can handle even dimensions
  INT *m = (INT *)calloc(dim1,sizeof(INT));
  INT *mm = (INT *)calloc(dim1,sizeof(INT));
  INT *mp = (INT *)calloc(dim1,sizeof(INT));
  INT *cnodes = (INT *)calloc(c2s->nvcube,sizeof(INT));  
  //  INT *icycle = (INT *)calloc(dim+1,sizeof(INT));
  INT nv=1,ns=1;
  for(i=0;i<dim;i++){
    nv*=(nd[i]+1);
    ns*=nd[i];
  }
  ns*=c2s->ns; /*multiply by the number of simplices in the unit cube
		 (2 in 2D and 6 in 3d and 24 in 4d*/
  scomplex *sc = (scomplex *)haz_scomplex_init(dim,ns,nv);
  //  fprintf(stdout,"\nFaces=(%d,%d)=(face,face_parent)\n",face,face_parent);fflush(stdout);
  for(kf=0;kf<sc->nv;kf++){
    coord_lattice(m,dim,kf,sc->nv,nd);
    for(i=0;i<dim;i++){
      /* THIS HERE HAS THE X COORD FIRST WHICH IS NOT WHAT ONE HAS IF
	 USING THE BIJECTION BETWEEN BINARY NUMBERS AND THE
	 COORDINATES OF VERTICES IN THE UNIT CUBE> SO WE REVERSE THE
	 ORDERING OF DIVISIONS SO THAT WE PARTITION FIRST X and so
	 on. so basically we come here with the last coordinate
	 first. That is why we also have (dim-i-1) instaed of i*/
      sc->x[kf*dim+(dim-i-1)]=((REAL )m[i])/((REAL )nd[i]);
      /* OLD: sc->x[kf*dim+i]=((REAL )m[i])/((REAL )nd[i]); */
    }    
    //      print_full_mat_int(1,dim,m,"m1=");
    //      fprintf(stdout,"; iglobal=%d;",kf);
  }
  ns=0;
  for(kf=0;kf<sc->nv;kf++){
    coord_lattice(m,dim,kf,sc->nv,nd);
    flag=0;
    for(i=0;i<dim;i++){
      if(m[i]==nd[i]) {flag=1; break;}
    }
    if(flag) continue;
    if(intype==-1) {
      type=0;
    }else{
      /*criss-cross in any D*/
      /* determine type; this is not fully rigorously justified, but
       works in d=2,3*/
      type=(m[0]+intype)%2;
      for(i=1;i<dim-1;i++){
	type+=2*(m[i]%2);
      }
      // this is a hack here to work in 2D. unclear how to do in 2D yet or 4D. 
      if(dim==2){type=(abs(m[1]-m[0])+intype)%dim;}
      if((m[dim-1]%2)) type=dim-type;
      if((!(dim%2)) && (type>=(dim))) {type%=(dim);}
      if(dim%2 && (type>(dim+1))) {type%=(dim+1);}
    }
    //    type=0;
    /*depending on the type, split a cube in simplices*/
    //    fprintf(stdout,"\ntype=%d\n",type+1);
    for(j=0;j<c2s->nvcube;j++){
      //            fprintf(stdout,"j:%d; ",j+1);
      for(i=0;i<dim;i++){
    	mm[i]=m[i]+(c2s->bits[dim*j+i]);
	//		fprintf(stdout,"mm[%d]=%d; ",i+1,mm[i]+1);
      }
      cnodes[j]=num_lattice(mm,dim,nd);
      //            fprintf(stdout,"\n"); fflush(stdout);
    }   
    //    fprintf(stdout,"\n"); fflush(stdout);
    //  
    for(i=0;i<c2s->ns;i++){
      for(j=0;j<dim1;j++){
	iz1=c2s->nodes[i*dim1+j];
	jperm=c2s->perms[type*c2s->nvcube+iz1];
	//	fprintf(stdout,"\ntype=%d,ns=%d,jperm=%d,iz1=%d",type,ns,jperm,iz1);
	sc->nodes[ns*dim1+j]=cnodes[jperm];
      }
      ns++;
    }    
  }
  if(face>=0 && face_parent>=0 && (sc_parent!=NULL)){
    //    REAL shift[2]={0.,-1.};
    INT kfp,ijk,mi,mip,toskip,toadd;
    //    print_full_mat_int(1,c2s->nvface,(c2s->faces+face*c2s->nvface),"face");
    //print_full_mat_int(1,c2s->nvface,(c2s->faces+face_parent*c2s->nvface),"facep");
    if(face<dim){
      mi=dim-(face+1);
      toskip=0;
    } else {
      mi=dim-((face%dim)+1);
      toskip=nd[mi];
    }
    if(face_parent<dim){
      mip=dim-(face_parent+1);
      toadd=0;
    } else {
      mip=dim-((face_parent%dim)+1);
      toadd=nd_parent[mip];
    }
    //    fprintf(stdout,"\ntoskip =%d, toadd=%d,mi=%d,mip=%d",toskip,toadd,mi,mip);
    nv=0;  
    for(kf=0;kf<sc->nv;kf++){
      coord_lattice(m,dim,kf,sc->nv,nd);
      if(m[mi]==toskip){
	memcpy(mp,m,dim*sizeof(INT));
	mp[mip]=toadd;
	kfp=num_lattice(mp,dim,nd_parent);
	//	print_full_mat_int(1,c2s->n,nd,"ndnd");
	//	print_full_mat_int(1,c2s->n,nd_parent,"ndndp");
	//	print_full_mat_int(1,c2s->n,m,"m");
	//	print_full_mat_int(1,c2s->n,mp,"mp");
	/* fprintf(stdout,"\nskipping kf=%d, kfp=%d",kf,kfp); */
	/* fprintf(stdout,"\noldx=("); */
	/* for(ijk=0;ijk<dim;ijk++) */
	/*   fprintf(stdout,"%.5e ", sc_parent->x[kfp*dim+ijk]); */
	/* fprintf(stdout,"); newx=["); */
	/* for(ijk=0;ijk<dim;ijk++) */
	/*   fprintf(stdout,"%.5e ", sc->x[kf*dim+ijk]); */
	/* fprintf(stdout,"]\n"); */
      }
    }
  }
  if(m) free(m);
  if(mp) free(mp);
  if(mm) free(mm);
  return sc;
}

