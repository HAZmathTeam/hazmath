#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "hazmath.h"
#include "grid_defs.h"
#include "grid_params.h"
/*****************************************************************/
/******* main_mesh.c ****************************************************
 *
 * Routines to generate meshes for HAZMATH.  Outputs in a "HAZMATH" format
 * and VTK format.  Can do 2D triangular mesh and 3D tetrahedral mesh.
 * Also capable of making domains with holes.
 *
 *************************************************************************/
/*******************************************************/
INT main(INT   argc,   char *argv[])
{
  INT dim=DIM;
  REAL *xcoord=NULL, *ycoord=NULL, *zcoord=NULL; 
  INT *je=NULL, *ib=NULL, *mask=NULL;             // structure of mesh
  INT nel,dim1;
  INT nx, ny, nz,nvert;
  INT i=-10,j=-10,k=-10;
  // Set Paramaters
  INT ref_levels=MAXREFLEVELS; // export also a vtk file with the mesh
  INT idovtk=VTKDO; // export also a vtk file with the mesh
  /* Filenames */
  char *nameout=NULL, *namevtk=NULL;

  char *prefix0=NULL; /*prefix for all filenames */
  size_t lenname, slen;
  /* deal with filenames */ 
  slen=strlen(OPREFIX);
  if(slen > (FILENAMESIZE-6))
    slen=FILENAMESIZE-6;
  prefix0=strndup( OPREFIX , slen);
  lenname=strlen(prefix0)+6;
  /* this is a haz out name, should always be there */
  nameout = (char *)calloc(lenname,sizeof(char));
  (void)strcpy(nameout,prefix0);
  if(idovtk){
    namevtk = (char *)calloc(lenname,sizeof(char));
    (void)strcpy(namevtk,prefix0);
  }
  /* max points in one direction */
  INT m0x = 1025, m0y = 1025, m0z = 83;
  /* haz file name */
  (void)strcat(nameout,".grd");
  /* print filenames */  
  fprintf(stdout,"%%Creating mesh: output goes on: mesh-out=%s",nameout);
  fprintf(stdout, "\n%%INPUT number of points in X-direction:\n");
  i=fscanf(stdin,"%d", &nx);  
  nx=chkn(nx, 2, m0x);
  if(!(nx%2)){
    nx++;
    fprintf(stdout, "\n%%Even number of points in X-direction. Changing to ODD, i.e. nx=%d\n",nx);
  }
  fprintf(stdout, "\n%%INPUT number of points in Y-direction:\n");
  i=fscanf(stdin,"%d", &ny);
  ny=chkn(ny, 2, m0y);
  if(!(ny%2)){
    ny++;
    fprintf(stdout, "\n%%Even number of points in Y-direction. Changing to ODD, i.e. ny=%d\n",ny);
  }
  nvert = nx*ny;             /* Total of DOF (Including the
				boundary) */
  nel = (nx-1)*(ny-1) * 2;   /* Total num. of elements 2 tris in a
				square */

  if(dim>2){
    fprintf(stdout, "\nINPUT number of points in Z-direction:\n");
    i=fscanf(stdin,"%d", &nz);
    nz=chkn(nz, 2, m0z);
    if(!(nz%2)){
      ny++;
      fprintf(stdout, "\n%%Even number of points in Z-direction. Changing to ODD., i.e. nz=%d\n",nz);
    }
    nvert = nvert*nz; /* Total of DOF (including the boundary)*/
    nel = nel*(nz-1)*3; /* Total # of elements: six tets in a cube */
  }  
  fprintf(stdout,"\n");
  dim1 = dim + 1; /* four vertecies in a tet; three in a triangle */    
  fprintf(stdout,"\n-------------------------------------------\n");
    
  fprintf(stdout,"%%NX=%d,NY=%d,NZ=%d\n",nx,ny,nz);
  fprintf(stdout,"%%TOTAL NUMBER OF VERTICES ON THE LATICE MESH=%d\n",nvert);
  
  fprintf(stdout,    "-------------------------------------------\n");
    
  // Simplicial complex: 
  scomplex *sc = (scomplex *)haz_scomplex_init(dim,nel,nvert);
  je = sc->nodes;
  ib = sc->bndry;
  /* je = (INT *)calloc(dim1*nel,sizeof(INT)); */
  /* ib = (INT *)calloc(nvert,sizeof(INT)); */
  mask = (INT *)calloc(nvert,sizeof(INT));    
  /* one array for all coords */
  /* REAL *x = (REAL *)malloc(dim*nvert*sizeof(REAL)); */
  REAL *x = sc->x;
  xcoord = x;
  ycoord = xcoord+nvert;
  if(dim > 2) 
    zcoord = ycoord+nvert;
  else
    zcoord=NULL;
  REAL *xo=(REAL *)calloc(2*dim,sizeof(REAL));
  REAL *xn=xo+dim;  
  if(je && ib && mask && (sc->flags) && xcoord && ycoord){
    if(dim>2 && zcoord) {
      getm3_(&nx,&ny,&nz,&nvert,&nel,xcoord,ycoord,zcoord,	\
	     je,sc->flags,ib,mask);
    } else {
      getm2_(&nx,&ny,&nvert,&nel,xcoord,ycoord,	\
	     je,sc->flags,ib,mask);
    }
    for(j=0;j<(sc->n+1)*sc->ns;j++){
      sc->nodes[j]-=1;
    }
    r2c(dim,nvert,sizeof(REAL),sc->x); /*
					 sc->x by cols is the same as
					 xcoord by rows;   
				       */
    //    swne(dim,nvert,xo,xn,sc->x);
    if(ref_levels){
      /*******************************************/
      features *feat=malloc(sizeof(feat));
      feat->nbig=dim;
      feat->n=dim;
      feat->nf=1;
      feat->x = (REAL *)calloc(feat->nbig*feat->nf,sizeof(REAL));
      for(i=0;i<feat->nf;i++){
	for(j=0;j<dim-1;j++)
	  feat->x[dim*i+j]=0.5;
	feat->x[dim*i+(dim-1)]=0.;
      }
      find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
      refining(ref_levels,sc,feat->nf,feat->x);
      free(feat->x);
    }
    hazw(nameout,sc,0,1);
    if(idovtk) {
      strcat(namevtk,".vtu");
      vtkw(namevtk,sc,0,0,1.);
    }
  }
  if(nameout) free(nameout);
  if(namevtk) free(namevtk);
  if(xo) free(xo); 
  if(mask) free(mask);
  return 0;
}

