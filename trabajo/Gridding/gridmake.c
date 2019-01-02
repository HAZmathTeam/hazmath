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
void scfinalize(scomplex *sc);
/*****************************************************************/
scomplex *unimesh(INT dim, INT *nd, REAL *xo, REAL *xn,	const INT ishift){
  INT i,j,dfactorial,nv=-10,ns=-10;
  /*CONSTANTS and BOUNDARY CODES:*/
  /* max points in every direction; max dimension = 10*/  
  INT m0[10]={1025,1025,1025,1,1,1,1,1,1,1};  
  INT  minneu=(INT )MARKER_NEUMANN;
  INT  maxneu=((INT )MARKER_ROBIN)-1;
  INT *ibcode=(INT *)calloc(dim*(dim-1),sizeof(INT));
  ibcode[0]=(INT )LEFTBC;
  ibcode[1]=(INT )RIGHTBC;
  ibcode[2]=(INT )FRONTBC;
  ibcode[3]=(INT )BACKBC;
  ibcode[4]=0;
  ibcode[5]=0;
  if(dim>2) {
    ibcode[4]=(INT )BOTTOMBC;
    ibcode[5]=(INT )TOPBC;
  }
  /*ONLY 3D above; ANY d below after the fortran is removed*/
  for(j=0;j<dim;j++){
    nd[j]=chkn(nd[j], 2, m0[j]);
    if(!(nd[j]%2)){
      nd[j]++;
      /*no warning*/
      /*      fprintf(stdout,"\n%%WARNING: Even number of points direction x[%d]. Changing to ODD, i.e. nd[%d]=%d\n",j,j,nd[j]);*/
    }
  }
  nv=nd[0];
  ns=(nd[0]-1);
  dfactorial=1;
  for(j=1;j<dim;j++){
    nv *= nd[j];
    ns *= (nd[j]-1);
    dfactorial *= (j+1);
  }
  /*  If the cube has coordinates x[0,...,n-1], where n is the spatial
      dimension, then for every permutation p[0,...,n-1] there is a
      simplex given by
      
      (0 .le. x[p[0]] .le. x[p[1]] .le. ... .le. x[p[n-1]] .le. 1)
  */
  ns*=dfactorial; /* n-factorial simplices in each cube. */
  // form the initial scomplex
  scomplex *sc = (scomplex *)haz_scomplex_init(dim,ns,nv);
  INT *mask = (INT *)calloc(sc->nv,sizeof(INT));    
  REAL *xcoord = sc->x;
  REAL *ycoord = xcoord+nv;
  REAL *zcoord = ycoord+nv;
  if(dim<=2) {
    getm2_(nd,&nv,&ns,xcoord,ycoord,		\
	   sc->nodes,sc->flags,sc->bndry,mask,
	   ibcode,&minneu,&maxneu);
    zcoord=NULL;
  } else {
    getm3_(nd,&nv,&ns,xcoord,ycoord,zcoord,			\
	   sc->nodes,sc->flags,sc->bndry,mask,			\
	   ibcode,&minneu,&maxneu);
  }
  if(mask) free(mask);
  if(ibcode) free(ibcode);
  r2c(sc->n,sc->nv,sizeof(REAL),sc->x); /*sc->x by cols is the same as
					  xcoord by rows;*/
  if(ishift){
    for(j=0;j<(sc->n+1)*sc->ns;j++){
      sc->nodes[j]-=ishift;
    }
  }
  return sc;
}
/******* main_mesh.c ****************************************************
 *
 * Routines to generate meshes for HAZMATH.  Outputs in a "HAZMATH" format
 * and VTK format.  Can do 2D triangular mesh and 3D tetrahedral mesh.
 *
 *************************************************************************/
/*******************************************************/
INT main(INT   argc,   char *argv[]){
  INT dim=DIM;
  // Set Paramaters
  INT ref_levels=MAXREFLEVELS; // export also a vtk file with the mesh
  INT use_features=(INT )USE_FEATURES;
  INT idovtk=VTKDO; // export also a vtk file with the mesh
  /* Filenames */
  char *nameout=NULL, *namevtk=NULL;
  char *prefix0=NULL; /* prefix for all filenames */
  size_t lenname, slen;
  /* deal with filenames */ 
  slen=strlen(OPREFIX);
  if(slen > (FILENAMESIZE-6))
    slen=FILENAMESIZE-6;
  prefix0=strndup( OPREFIX , slen);
  lenname=strlen(prefix0)+6;
  /* this is a haz out name, should always be assigned */
  nameout = (char *)calloc(lenname,sizeof(char));
  (void)strcpy(nameout,prefix0);
  if(idovtk){
    namevtk = (char *)calloc(lenname,sizeof(char));
    (void)strcpy(namevtk,prefix0);
  }
  /* haz file name */
  (void)strcat(nameout,".grd");
  /* print filenames */
  /**********************************************************************************/
  /* SHORT status; */

  /* // Overall CPU Timing */
  /* clock_t clk_overall_start = clock(); */

  /* // Set Parameters from Reading in Input File */
  /* input_param inparam; */
  /* param_input_init(&inparam); */
  /* param_input("./input.dat", &inparam); */

  /* // Open gridfile for reading */
  /* printf("\nCreating mesh and FEM spaces:\n"); */
  /* FILE* gfid = HAZ_fopen(inparam.gridfile,"r"); */

  /* // Create the mesh */
  /**********************************************************************************/
  INT i,j;
  INT *nd=(INT *)calloc(dim,sizeof(INT));
  fprintf(stdout,"%%Creating mesh: output goes on: mesh-out=%s",nameout);
  for(j=0;j<dim;j++){
    fprintf(stdout, "\n%%INPUT number of points in direction x[%d]:\n",j);
    i=fscanf(stdin,"%d", (nd+j));
  }
  // coordinates of the point where (0,0,...,0) needs to be mapped;
  // coordinates of the point where (1,1,...,1) is mapped;
  REAL xo[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  REAL xn[10]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};
  // boundary codes
  const INT ishift=1;
  scomplex *sc=unimesh(dim,nd,xo,xn,ishift);
  fprintf(stdout,"\n-------------------------------------------\n");
  fprintf(stdout,"%% UNIFORM MESH: ");
  for(j=0;j<(dim-1);j++){
    fprintf(stdout,"nd[%d]=%d, ",j,nd[j]);
  }
  fprintf(stdout,"nd[%d]=%d",dim-1,nd[dim-1]);
  fprintf(stdout,"\n-------------------------------------------\n");    
  if(ref_levels){
    dvector *errest=(dvector *)malloc(sizeof(dvector));
    long int lots = (1<<ref_levels)*sc->ns;
    //fprintf(stdout,"\n*** NNN=%li %d\n",lots,(1<<ref_levels));
    /* create a vector to contain the errors of all refinements*/
    /*LOTS=number of element if every element was refined ref_level times*/
    *errest=dvec_create((INT )lots);
    /*******************************************/
    /* if we have to refine only elements containing points from an
       array called "features array" */
    if(use_features){
      features *feat=malloc(sizeof(features));
      feat->nbig=dim;
      feat->n=dim;
      REAL vvv=0.;
      features_r(dim,use_features,feat,vvv);
      errest->row=feat->nf;
      errest->val=feat->x;
    } else {
      errest->row=0;
      errest->val=NULL;
    }
    n_refine(ref_levels,sc,errest);  
    // finalize grid (take only the top level and forget the hierarchy:
    scfinalize(sc);
  }
  //  haz_scomplex_print(sc,0,__FUNCTION__);  fflush(stdout);
  // write the output mesh file:    
  hazw(nameout,sc,0,1);
  if(idovtk) {
    strcat(namevtk,".vtu");
    //    fprintf(stdout,"\nvert=%d; simp=%d",sc->nv,sc->ns);
    vtkw(namevtk,sc,0,0,1.);
  }
  if(nameout) free(nameout);
  if(namevtk) free(namevtk);
  if(nd)free(nd);
  if(sc)free(sc);
  return 0;
}
