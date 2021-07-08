/*! \file examples/basic_elliptic/triangulation2d.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 2019/01/09.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program generates simplicial grids in 2,3,4... dimension.
 *
 * \note This example highlights some of the features of the simple
 * mesh generator included with HAZmath. It is only to illustrate how
 * to use the mesh refinement.
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
/*********************************************************************/
#include "connected_components.h"
/*********************************************************************/
INT main()
{
  // construction from a list of points :
  INT ncir=-20,nsq=-20;
  INT np=-20,neold=-20,ne=-20,dim,ii,jj,i,j,k;
  REAL pi=M_PI, hmax=1e0, shape_bound=0.125e0;
  REAL r, a, h, sx,sy, twopi=2e0*pi;
  //////////////////////////////////////////////////////////
  dim=2; /// dimension
  REAL *c=calloc(dim,sizeof(REAL)); // seed point which is inside
					 // a hole of the domain; for
					 // a circular "hole" this is
  					 // the center;
  FILE *fp; 
  // fp=stdin; 
  fp=fopen("hole_grid.input","r");
  //read input:
  for(j=0;j<dim;j++){
    fscanf(fp,"%lg",(c+j));
  }
  fscanf(fp,"%lg",&r);// radius
  fscanf(fp,"%i",&ncir);//number of pts on circle
  fscanf(fp,"%lg",&sx); //length of square edge parallel to x
  fscanf(fp,"%lg",&sy); //length of square edge parallel to y
  fclose(fp);
  //end read input
  nsq=4;// points on the square
  np=ncir+nsq; //total number of points to describe the domain
  REAL *x_in=(REAL *)calloc(np*dim,sizeof(REAL));// points on the
						 // boundary of the
						 // domain.
  h=twopi/((REAL ) ncir);
  a=0.;
  // points(circle)
  for(i=0;i<ncir;++i){
    ii=dim*i;
    x_in[ii]  =c[0] + r*cos(a);
    x_in[ii+1]=c[1] + r*sin(a);
    a += h;
  }
  // we need half lengths of the sides of the rectangle
  sx=sx*0.5;
  sy=sy*0.5;
  // (1,1);(0,1);(0,0);(1,0)
  REAL *xx=x_in+ncir*dim;
  // Corners of the square centered at the center of the circle:
  // NE
  xx[0] = c[0]+sx; xx[1]= c[1]+sy;
  // NW
  xx[2] = c[0]-sx; xx[3]= c[1]+sy;
  // SW
  xx[4] = c[0]-sx; xx[5]= c[1]-sy;
  // SE
  xx[6] = c[0]+sx; xx[7]= c[1]-sy;
  //
  /* for(i=0;i<np;++i){ */
  /*   fprintf(stdout,"\nx[%d]=(",i); */
  /*   fprintf(stdout,"%.4f,%.4f)",x_in[dim*i],x_in[dim*i+1]); */
  /* } */
  /* fprintf(stdout,"\n"); */
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
  //  fprintf(stdout,"\nCheck:neold-ne=%d\nEdges:\n",neold-ne);
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
    /* fprintf(stdout,"\ne(%d)=(%d,%d)=[(",i,ii/dim,jj/dim); */
    /* fprintf(stdout,"%.3f,%.3f)--(%.3f,%.3f)]",	\ */
    /* 	    x_in[ii],x_in[ii+1],		\ */
    /* 	    x_in[jj],x_in[jj+1]);	       */
    d=0e0;
    for(j=0;j<dim;++j){
      dxy=x_in[ii+j]-x_in[jj+j];
      d+=dxy*dxy;
    }
    if(hmax>sqrt(d)) hmax=sqrt(d);
  }
  //  fprintf(stdout,"\nhmax=%f\n\n",hmax);
  //scomplex *sc=(scomplex *)malloc(sizeof(scomplex));
  scomplex *sc=(scomplex *)haz_scomplex_init(dim,0,0,dim);    
  //    exit(55);
  sc->n=dim;
  sc->nbig=dim;
  INT dim1=dim+1;
  sc->nv=0;
  sc->ns=0;
  c_do_tri_2d(x_in,edges,					\
	      c,						\
	      np, ne, dim,					\
	      shape_bound,					\
	      hmax,
	      &sc->nodes,&sc->x,	\
	      &sc->nv,&sc->ns);
  free(c);
  free(x_in);
  free(edges);
  //
  sc->nbr=(INT *)calloc(dim1*sc->ns,sizeof(INT));
  find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
  haz_scomplex_init_part(sc);
  find_cc_bndry_cc(sc);
  // add space at the end but 12 is one byte less so space is
  // overwritten with \0 as it should.
  char *fname=strndup("hole_XXX.yyy ",12*sizeof(char));
  size_t nbyt=(strlen(fname)+1)*sizeof(char);
  snprintf(fname,nbyt,"hole%03i.haz",ncir);
  //  fprintf(stdout,"\n\nFNAME=|%s|\n",fname);
  /// output-haz
  hazw(fname,sc,0);
  /// output-vtu
  snprintf(fname,nbyt,"hole%03i.vtu",ncir);
  vtkw(fname,sc,0,1.);
  /*FREE*/
  haz_scomplex_free(sc);
  free(fname);
  return 0;
}
