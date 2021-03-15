#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hazmath.h"
INT main()
{
  // construction from a list of points :
  INT level=3;
  INT ncir=1<<level,nsq=1<<2; // ncir=2^k,nsq=4
  INT np=-20,neold=-20,ne=-20,dim,ii,jj,i,j,k;
  REAL seedpt[2]={0.5,0.5}; // seed point which is inside a hole of the domain;
  REAL pi=M_PI, hmax=1e0, shape_bound=0.125e0;
  REAL r, a, h, twopi=2e0*pi;
  //////////////////////////////////////////////////////////
  dim=2; /// dimension
  np=ncir+nsq; //total number of points to describe the domain
  REAL *x_in=(REAL *)calloc(np*dim,sizeof(REAL));
  h=twopi/((REAL ) ncir);
  r=0.1e0;
  a=0.;
  // points(circle)
  for(i=0;i<ncir;++i){
    ii=dim*i;
    x_in[ii]  =seedpt[0] + r*cos(a);
    x_in[ii+1]=seedpt[1] + r*sin(a);
    a += h;
  }
  // (1,1);(0,1);(0,0);(1,0)
  REAL *xx=x_in+ncir*dim;
  xx[0] = 1.; xx[1]= 1.;
  xx[2] = 0.; xx[3]= 1.;
  xx[4] = 0.; xx[5]= 0.;
  xx[6] = 1.; xx[7]= 0.;
  for(i=0;i<np;++i){
    fprintf(stdout,"\nx[%d]=(",i);
    fprintf(stdout,"%.4f,%.4f)",x_in[dim*i],x_in[dim*i+1]);
  }
  fprintf(stdout,"\n");
  //edges(circle)
  neold=np; // as many edges as points (2 loops around circle and
  // around the square,
  // array for edge constraints:
  INT *edges=(INT *)calloc(2*neold,sizeof(INT));
  ne=0;
  for(i=0;i<(ncir-1);++i){
    edges[2*ne]=i;
    edges[2*ne+1]=i+1;
    ne++;
  }
  // ncir-1 edges are done, we add the last one:
  edges[2*ne]=0;
  edges[2*ne+1]=ncir-1;
  ne++;
  // edges(square):
  for(i=0;i<(nsq-1);++i){
    ii=ncir+i;
    edges[2*ne]=ii;
    edges[2*ne+1]=ii+1;
    ne++;
  }
  edges[2*ne]=ncir;
  edges[2*ne+1]=ncir+nsq-1;
  ne++;
  fprintf(stdout,"\nCheck:neold-ne=%d\nEdges:\n",neold-ne);
  /* for(i=0;i<ne;++i){ */
  /*   ii=edges[2*i]; */
  /*   jj=edges[2*i+1]; */
  /*   fprintf(stdout,"\ne(%d)=(%d,%d)=[(",i,ii,jj); */
  /*   fprintf(stdout,"%.3f,%.3f)--(%.3f,%.3f)]",	\ */
  /* 	      x_in[dim*ii],x_in[dim*ii+1],		\ */
  /* 	      x_in[dim*jj],x_in[dim*jj+1]);	       */
  /* } */
  /* fprintf(stdout,"\n");fflush(stdout); */
  shape_bound=0.125;
  REAL d,dxy;
  // first edge lenngth:
  ii=dim*edges[0];jj=dim*edges[1];
  d=(x_in[ii]-x_in[jj])*(x_in[ii]-x_in[jj])		\
    + (x_in[ii+1]-x_in[jj+1])*(x_in[ii+1]-x_in[jj+1]);
  hmax=sqrt(d);
  // compute the min edge length to put it as a constraint for the
  // edge length in the triangulation
  for(i=0;i<ne;++i){
    ii=dim*edges[2*i];
    jj=dim*edges[2*i+1];
    fprintf(stdout,"\ne(%d)=(%d,%d)=[(",i,ii/dim,jj/dim);
    fprintf(stdout,"%.3f,%.3f)--(%.3f,%.3f)]",	\
	    x_in[ii],x_in[ii+1],		\
	    x_in[jj],x_in[jj+1]);	      
    d=0e0;
    for(j=0;j<dim;++j){
      dxy=x_in[ii+j]-x_in[jj+j];
      d+=dxy*dxy;
    }
    if(hmax>sqrt(d)) hmax=sqrt(d);
  }
  fprintf(stdout,"\nhmax=%f\n\n",hmax);
  //scomplex *sc=(scomplex *)malloc(sizeof(scomplex));
  scomplex *sc=(scomplex *)haz_scomplex_init(dim,0,0);    
  //    exit(55);
  sc->nv=0;
  sc->ns=0;
  c_do_tri_2d(x_in,edges,					\
	      seedpt,						\
	      np, ne, dim,					\
	      shape_bound,					\
	      hmax,
	      &sc->nodes,&sc->nbr,&sc->x,	\
	      &sc->nv,&sc->ns);
  free(x_in);
  free(edges);
  return 0;
}
