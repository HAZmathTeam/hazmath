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
#include "elliptic_p1_supporting.h"
/****************************************************************************/
/* 
 * refinement type: .gt. 10 is uniform refinement and .le. 10
 *                  (typically 0) is the newest vertex bisection
*/
#ifndef REFINEMENT_TYPE
#define REFINEMENT_TYPE 11
#endif
/**/
#ifndef REFINEMENT_LEVELS
#define REFINEMENT_LEVELS 4
#endif

#ifndef SPATIAL_DIMENSION
#define SPATIAL_DIMENSION 3
#endif
/****************************************************************************/
int main(int argc, char *argv[])
{
  INT ref_levels=REFINEMENT_LEVELS;
  INT dim=SPATIAL_DIMENSION;
  INT ref_type=REFINEMENT_TYPE; // >10 uniform refinement;
  //
  INT jlevel,k;
  scomplex *sc=NULL,*sctop=NULL;
  /**/
  ivector marked;  marked.row=0;  marked.val=NULL;
  dvector sol;  sol.row=0;  sol.val=NULL;
  sc=mesh_cube_init(dim,ref_type);
  /**/
  if(sc->ref_type>10){
    //      fprintf(stdout,"\ndim=%dd(uniform):level=",dim);
    if(dim==3){
      for(jlevel=0;jlevel<ref_levels;++jlevel){
	//	fprintf(stdout,".%d.",jlevel+1);fflush(stdout);
	uniformrefine3d(sc);
	sc_vols(sc);
      }
    } else if(dim ==2){
      for(jlevel=0;jlevel<ref_levels;++jlevel){
	//	fprintf(stdout,".%d.",jlevel+1);fflush(stdout);
	uniformrefine2d(sc);
	sc_vols(sc);
      }
    } else {
      fprintf(stderr,"\n%%%% ERROR: WRONG spatial dimension for uniform refinement");
    }
    find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
    sc_vols(sc);
  } else {
    //  fprintf(stdout,"\n");
    // end intialization
    //    fprintf(stdout,"\nlevels=%d",ref_levels);fflush(stdout);
    //    find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
    // sc_vols(sc); these shouls all be found already
    //    fprintf(stdout,"\ndim=%dd(newest):level=",dim);
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
    scfinalize(sc);
    sc_vols(sc);
  }
  //  icsr_print_matlab(stdout,sc->parent_v);
  //  haz_scomplex_print(sc,0,__FUNCTION__);
  sol=fe_sol(sc,1.0,1.0);
  scomplex *dsc=malloc(sizeof(scomplex));//boundary scomplex
  dsc[0]=sc_bndry(sc);
  //  haz_scomplex_print(dsc,0,__FUNCTION__);
  SHORT todraw=1;
  draw_grids(todraw, sc,dsc, &sol);
  /* write the output mesh file:    */
  //  hazw("output/mesh.haz",sc,0);
  //  hazw("output/dmesh.haz",dsc,0);
  dvec_free(&sol);
  ivec_free(&marked);
  haz_scomplex_free(sc);  
  haz_scomplex_free(dsc);
  return 0;
}
