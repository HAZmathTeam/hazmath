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
#ifndef INT
#define INT int
#endif
#ifndef REAL
#define REAL double
#endif
/********************************LIBLIB */
void cube2simp_free(cube2simp *c2s)
{
  if(c2s->bits)free(c2s->bits);
  if(c2s->nodes)free(c2s->nodes);
  if(c2s->perms)free(c2s->perms);
  if(c2s)free(c2s);
  return;
}

static void binary0(cube2simp *c2s)
{
  // stores in an array the coordinates of the vertices of the unit  
  // cube in dimension (dim). Lexicographical ordering from 0,0,...0
  // to 1,1,...,1
  
  INT dim=c2s->n,nvcube=c2s->nvcube;
  INT shift,oposite,i,j,k,kn,nbits=c2s->n-1;
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
INT reverse(void *arr,INT length, size_t elsize)
{
  /* 
     permutes a void array whose elements are of size elsize 
     a[0],...a_[length-1]-->a[length-1],...a[0]. 
  */
  INT i,j,k,nnn=(INT)(length/2);
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
  return 0;
}
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
  i=1; for(i=1;i<=c2s->n;i++)c2s->ns*=i;
  c2s->ne=c2s->n*(1<<(c2s->n-1)); /* number of edges in the cube.*/
  c2s->nf=2*c2s->n; /* number of n-1 dimensional faces in the cube */
  /////////////////////////////////////////////////////
  c2s->edges=(INT *)calloc(2*c2s->ne,sizeof(INT));
  c2s->bits=(unsigned INT *)calloc(c2s->n*c2s->nvcube,sizeof(unsigned INT));
  c2s->nodes=(INT *)calloc(c2s->ns*(c2s->n+1),sizeof(unsigned INT));
  c2s->perms=(INT *)calloc(c2s->nvcube*(c2s->n+1),sizeof(unsigned INT));
  memset(c2s->nodes,0,(c2s->n+1)*c2s->ns*sizeof(INT));
  /*end of allocation*/
  INT k1,kn1,k2,kn2,dim1=c2s->n+1,		\
    ns=c2s->ns,nvcube=c2s->nvcube;
  binary0(c2s);
  INT *edges=c2s->edges;
  memset(edges,0,c2s->ne*sizeof(INT));
  unsigned INT numbits=22;
  unsigned INT *b1,*b2;
  //  fprintf(stdout,"Memory: edges=%d,ns=%d\n",ne,ns);
  INT nedge=0,nvcube1=nvcube-1;
  for(k1 = 0;k1<nvcube1;k1++){
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
      if(node==(nvcube1)) continue;
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
  return c2s;
}
scomplex *umesh(const INT dim, INT *nd, cube2simp *c2s, const SHORT intype)
{
  /* 
     uniform simplicial mesh of the unit cube in dimension dim.
     dim is the dimension, nd is the number of grid points in each
     dimension.  ordering is lexicographically by
     name=(x[0],...,x[n]).  more than 3D is not fully tested xmacro[]
     are the coordinats of a domain isomorphic to the cube via a
     bilinear change of coordinates.
     output is a simplicial complex sc. 
  */
  INT iz1;
  INT jperm,k,i,j,flag,kf,type;
  INT dim1 = dim+1;
  // m is dim+1 so that we can handle even dimensions
  INT *m = (INT *)calloc(dim1,sizeof(INT));
  INT *mm = (INT *)calloc(dim1,sizeof(INT));
  INT *cnodes = (INT *)calloc(c2s->nvcube,sizeof(INT));  
  //  INT *icycle = (INT *)calloc(dim+1,sizeof(INT));
  INT nv=1,ns=1;
  for(i=0;i<dim;i++){
    nv*=(nd[i]+1);
    ns*=nd[i];
  }
  ns*=c2s->ns; /**/
  fprintf(stdout,"\nGenerating uniform mesh in dim=%d; vertices: %d, simplexes %d\n",dim,nv,ns);fflush(stdout);
  scomplex *sc = (scomplex *)haz_scomplex_init(dim,ns,nv);
  nv=0;
  //  icycle[dim]=0;
  for(kf=0;kf<sc->nv;kf++){
    coord_lattice(m,dim,kf,sc->nv,nd);
    for(i=0;i<dim;i++){
      sc->x[kf*dim+i]=((REAL )m[i])/((REAL )nd[i]);
    }    
  }
  ns=0;
  for(kf=0;kf<sc->nv;kf++){
    coord_lattice(m,dim,kf,sc->nv,nd);
    flag=0;
    for(i=0;i<dim;i++){
      if(m[i]==nd[i]) {flag=1; break;}
    }
    if(flag) continue;
    /* determine type; this is not fully rigorously justified, but
       works in d=2,3*/
    type=(m[0]%2+intype);
    for(i=1;i<dim-1;i++){
      type+=2*(m[i]%2);
    }
    if((m[dim-1]%2)) type=dim-type;
    if(dim%2 && (type>(dim+1))) {type%=(dim+1);}
    if(!(dim%2) && (type>=(dim+1))) {type%=(dim);}
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
	fprintf(stdout,"\ntype=%d,ns=%d,jperm=%d,iz1=%d",type,ns,jperm,iz1);
	sc->nodes[ns*dim1+j]=cnodes[jperm];
      }
      ns++;
    }
  }
  if(m) free(m);
  return sc;
}
void polar2cart(INT dim, coordsystem *polar,REAL *px, REAL *cx)
{
  // polar is r, theta1,...theta[n-1]; cart is x[0]...x[n-1]
  // px are n-coordinates in polar coord system
  // converts polar coordnates to cartesian in d dimensions. 
  INT i,j;
  if(dim==1) return;
  if(polar->type!=1) return;
  REAL rho = px[0];
  for(i=0;i<dim;i++){
    cx[i]=rho*cos(px[i]);
    for(j=0;j<i;j++){
      cx[i]*=sin(px[j+1]);
    }
    if(i<(dim-1)) {
      cx[i]*=cos(px[i+1]);
    }
    cx[i]+=polar->o[i];// add the coords of the origin.
  }  
  return;
}
  
