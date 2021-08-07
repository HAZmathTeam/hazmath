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
#include "meshes_inline.h"
#include "elliptic_p1_supporting.h"
/****************************************************************************/
// refinement type: 1 is uniform and 0 is newest vertex bisection
#ifndef UNIFORM_REFINEMENT
#define UNIFORM_REFINEMENT 1
#endif
/////////////////////////////////////////////////////////////////
#ifndef REFINEMENT_LEVELS
#define REFINEMENT_LEVELS 9
#endif

#ifndef SPATIAL_DIMENSION
#define SPATIAL_DIMENSION 4
#endif
/****************************************************************************/
int main(int argc, char *argv[])
{
  SHORT uniref=UNIFORM_REFINEMENT;
  INT ref_levels=REFINEMENT_LEVELS;
  INT dim=SPATIAL_DIMENSION;// 2d,3d,4d... example
  INT jlevel,k;
  scomplex *sc=NULL;
  switch(dim){
  case 5:
    sc=mesh5d();
    uniref=0;
    break;
  case 4:
    sc=mesh4d();    
    uniref=0;
    break;
  case 3:
    sc=mesh3d();    
    if( uniref ){
      fprintf(stdout,"\ndim=%dd(uniform):level=",dim);
      for(jlevel=0;jlevel<ref_levels;++jlevel){
	fprintf(stdout,".%d.",jlevel+1);fflush(stdout);
	uniformrefine3d(sc);
	sc_vols(sc);
      }
      find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
      sc_vols(sc);
    }
    break;
  default:
    sc=mesh2d();
    if( uniref ){
      fprintf(stdout,"\ndim=%dd(uniform):level=",dim);
      for(jlevel=0;jlevel<ref_levels;++jlevel){
	fprintf(stdout,".%d.",jlevel+1);fflush(stdout);
	uniformrefine2d(sc);
      }
      find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
      sc_vols(sc);
    }
  }
  fprintf(stdout,"\n");
  scomplex *sctop=NULL;
  ivector marked;
  marked.row=0;
  marked.val=NULL;
  dvector sol;
  // end intialization
  if(! (uniref) ){ 
    //    fprintf(stdout,"\nlevels=%d",ref_levels);fflush(stdout);
    find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
    sc_vols(sc);
    //    fprintf(stdout,"\ndim=%dd(newest):level=",dim);
    for(jlevel=0;jlevel<ref_levels;++jlevel){
      //      fprintf(stdout,".%d.",jlevel+1);fflush(stdout);
      /* choose the finest grid */
      sctop=scfinest(sc);
      // solve the FE
      sol=fe_sol(sctop,1.0,1.0);
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
    scfinalize(sc);
    sc_vols(sc);
  }
  //  icsr_print_matlab(stdout,sc->parent_v);
  //  haz_scomplex_print(sc,0,__FUNCTION__);
  sol=fe_sol(sc,1.0,1.0);
  scomplex *dsc=malloc(sizeof(scomplex));//boundary scomplex
  dsc[0]=sc_bndry(sc);
  draw_grids((SHORT )1, sc,dsc, &sol);
  /* write the output mesh file:    */
  //  hazw("output/mesh.haz",sc,0);
  //  hazw("output/dmesh.haz",dsc,0);
  dvec_free(&sol);
  ivec_free(&marked);
  haz_scomplex_free(sc);  
  haz_scomplex_free(dsc);
  return 0;
}
