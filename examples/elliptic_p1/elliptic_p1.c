/*! \file examples/elliptic_p1/elliptic_p1.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note Solving Laplace equation with linear finite elements on sequence of grids. 
 *
 *  \note: modified by ltz on 20210531
 *
 */
/*********************************************************************/
#include "hazmath.h"
#include "supporting_elliptic_p1.h"
/****************************************************************************/
/****************************************************************************/

/* 
 * refinement type: .gt. 10 is uniform refinement and .le. 10
 *                  (typically 0) is the newest vertex bisection
*/
#ifndef REFINEMENT_TYPE
#define REFINEMENT_TYPE 0
#endif
/**/
#ifndef REFINEMENT_LEVELS
#define REFINEMENT_LEVELS 7
#endif
/**/
#ifndef SPATIAL_DIMENSION
#define SPATIAL_DIMENSION 3
#endif
/**/
#ifndef SET_BNDRY_CODES
#define SET_BNDRY_CODES 1
#endif
/****************************************************************************/
int main(int argc, char *argv[])
{
  INT ref_levels=REFINEMENT_LEVELS;
  INT dim=SPATIAL_DIMENSION;
  INT ref_type=REFINEMENT_TYPE; // >10 uniform refinement;
  //
  INT jlevel,k;
  scomplex **sc_all=NULL,*sc=NULL,*sctop=NULL;
  /**/
  ivector marked;  marked.row=0;  marked.val=NULL;
  dvector sol;  sol.row=0;  sol.val=NULL;
  sc_all=mesh_cube_init(dim,ref_type);
  sc=sc_all[0];
  /**/
  if(sc->ref_type>10){
    if(dim==3){
      for(jlevel=0;jlevel<ref_levels;++jlevel){
	uniformrefine3d(sc);
	sc_vols(sc);
      }
    } else if(dim ==2){
      for(jlevel=0;jlevel<ref_levels;++jlevel){
	uniformrefine2d(sc);
	sc_vols(sc);
      }
    } else {
      fprintf(stderr,"\n%%%% ERROR: WRONG spatial dimension for uniform refinement");
    }
    find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
    sc_vols(sc);
  } else {
    for(jlevel=0;jlevel<ref_levels;++jlevel){
      //      fprintf(stdout,".%d.",jlevel+1);fflush(stdout);
      /* choose the finest grid */
      sctop=scfinest(sc);
      // solve the FE
      //no need to solve here....      sol=fe_sol(sctop,1.0,1.0);
      /* mark everything; or use an estimator */
      marked.row=sctop->ns;
      marked.val=realloc(marked.val,marked.row*sizeof(INT));
      /*mark everything*/
      for(k=0;k<marked.row;k++) marked.val[k]=TRUE;
      // now we refine.
      refine(1,sc,&marked);
      /* free */
      dvec_free(&sol);
      haz_scomplex_free(sctop);
      /*  MAKE sc to be the finest grid only */
    }
  }
  scfinalize(sc,(INT )1);
  sc_vols(sc);
  ////////////////////////////////////////////// END MAP IT. 
  /* // now let us map it to a different domain: */
  /* // */
  /* REAL *vc=NULL; */
  /* if(dim==2){ */
  /*   REAL vc[]={-1.00, -2.00,			\ */
  /*   	       -2.00, -1.10,			\ */
  /*   	       0.00,  1.00,			\ */
  /*   	       0.50,  0.75};			\ */
  /*   mapit(sc,vc); */
  /* } else if(dim==3){ */
  /*   REAL vc[]={-1.00, -2.00, -1.00,		\ */
  /*   	       -2.00, -1.10, -1.00,		\ */
  /*   	       0.00,   1.00, -1.00,		\ */
  /*   	       0.50,   0.75, -1.00,		\ */
  /*   	       -1.00, -2.00,  1.00,		\ */
  /*   	       -2.00, -1.10,  1.00,		\ */
  /*   	       0.50,   0.75,  0.55,		\ */
  /*   	       0.00,   1.00,  0.75}; */
  /*   mapit(sc,vc); */
  /* } */
  ////////////////////////////////////////////// END MAP IT. 
  sol=fe_sol(sc,1.0,1.0);
  short todraw=1;
  draw_grids(todraw, sc,&sol);
  /* write the output mesh file:    */
  /* hazw("output/mesh.haz",sc,0); */
  dvec_free(&sol);
  ivec_free(&marked);
  haz_scomplex_free(sc_all[0]);  
  free(sc_all);  
  return 0;
}
