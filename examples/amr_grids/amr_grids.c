/*! \file examples/amr_grids/amr_grids.c
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
/* This macro definition below is amr_marking_type=44; and SHOULD BE
   MOVED TO MACROS or elseware later (ltz)*/
/*
  If set to something less than 2 it will refine non-stop until ref_levels is reached if there is at least one simplex containing at least one point. It is the maximum allowed number of features (nodes) per element. Any element containing more than this number of features is refined. 

*/
#ifndef MAX_NODES_PER_SIMPLEX
#define MAX_NODES_PER_SIMPLEX 5
#endif
/*********************************************************************/
#include "supporting_amr_grids.h"
/*********************************************************************/
//
INT main(INT   argc,   char *argv[])
{
  INT i;
  FILE *fp;
  fp=stdin;
  //no   fp=HAZ_fopen("inputs/2d_ann.input","r");
  //   fp=HAZ_fopen("inputs/2d_2L.input","r");
  // fp=HAZ_fopen("inputs/3d_fichera.input","r");
  // fp=HAZ_fopen("inputs/3d_2cubes_edge.input","r");
  // fp=HAZ_fopen("inputs/3d_2cubes_vertex.input","r");
  //  fp=HAZ_fopen("inputs/5d_cube.input","r");
  //  fp=HAZ_fopen("input/3d_cube.input","r");
  /*
    PARSE THE INPUT.
  */
  input_grid *g=parse_input_grid(fp);
  fclose(fp);
  scomplex **sc_all=generate_initial_grid(g);
  fprintf(stdout,"\nInitial mesh:\nElements = %12lld;\nVertices=%12lld\n",(long long )sc_all[0]->ns,(long long )sc_all[0]->nv); fflush(stdout);
  scomplex *sc=sc_all[0];
  INT ref_levels=g->nref,amr_marking_type=g->mark_type,j,k,kmarked;
  scomplex *sctop=NULL;
  dvector solfem,estimator;
  ivector marked;
  void *all=NULL;
  REAL *xstar=NULL;
  INT nstar,dim=sc->n;
  features feat;
  feat.n=0;feat.nbig=0;feat.x=NULL;feat.fill=-1e20;feat.fpf=NULL;
  iCSRmat node_ins;
  //NNNNNNNNNNNNNNNNNNNNNNNNNNNN
  fprintf(stdout,"\n****AMR_MARKING_TYPE=%d\n",amr_marking_type);
  if(amr_marking_type==44){
    char *data_file=strdup("./try_features.txt");
    feat.fpf=fopen(data_file,"r");
    feat.nbig=sc->n;
    feat.n=sc->n;
    feat.fill=-1e20;
    // last argument below is whether to map the simplicial complex to
    // a cube enclosing the data.
    j=features_r(&feat,sc,(INT )1,(REAL )1.1e0);
    free(data_file);
  }
  //NNNNNNNNNNNNNNNN
  if(amr_marking_type==0){
    // refine ref_levels;
    refine(ref_levels,sc,NULL);
  } else if (amr_marking_type==44){
    node_ins=icsr_create(0,0,0);
    nstar = feat.nf;
    xstar = feat.x;
    //
    INT max_nodes=(INT )MAX_NODES_PER_SIMPLEX ;
    if(max_nodes<=0) max_nodes=1;
    for(j=0;j<ref_levels;j++){
      sctop=scfinest(sc);
      /* MARK: marked is an ivector with num.rows=the number of
       *       simplices; its componenets are nonzero if the simplex
       *       is marked for refinement and 0 if it is not marked for
       *       refinement.
       */
      marked=mark_around_pts(sctop,sc,nstar,xstar,&node_ins,(const INT)max_nodes);
      kmarked=0;
      for(k=0;k<marked.row;++k)
	if(marked.val[k]) kmarked++;
      fprintf(stdout,"\n|lvl=%2lld|simplices[%2lld:%2lld]=%12lld|simplices[%2lld]=%12lld|",(long long int)j,0LL,(long long int)j,(long long int)sc->ns,(long long int)j,(long long int)sctop->ns);fflush(stdout);      
      refine(1,sc,&marked);
      if(!kmarked){
	fprintf(stdout,"\nthere were no simplices containing > %lld points. Exiting",(long long )MAX_NODES_PER_SIMPLEX);
		ivec_free(&marked);
		haz_scomplex_free(sctop);
		break;
      }
      ivec_free(&marked);
      haz_scomplex_free(sctop);
    }
    fprintf(stdout,"\n");
    icsr_free(&node_ins);
    ivec_free(&marked);
    free(feat.x);
  } else if(amr_marking_type==33){    
    REAL h = 0.05;  // step distance of points
    REAL threshold = h; // threshold for close to the points or not
    //
    nstar = g->num_refine_points;
    xstar = g->data_refine_points;
    //
    for(j=0;j<ref_levels;j++){
      /*
       * SELECT the finest grid:
       */
      sctop=scfinest(sc);
      /* MARK: marked is an ivector with num.rows=the number of
       *       simplices; its componenets are nonzero if the simplex
       *       is marked for refinement and 0 if it is not marked for
       *       refinement.
       */
      marked=mark_near_points(sctop,nstar,xstar, threshold);
      fprintf(stdout,"\n|lvl=%2lld|simplices[%2lld:%2lld]=%12lld|simplices[%2lld]=%12lld|",(long long int)j,0LL,(long long int)j,(long long int)sc->ns,(long long int)j,(long long int)sctop->ns);fflush(stdout);
      refine(1,sc,&marked);
      /* free */
      ivec_free(&marked);
      haz_scomplex_free(sctop);
    }
    ivec_free(&marked);
    //free(xstar);
  } else {
    /*
      Use "all" here can pass data around. Below we make 4 dvectors
      and one ivector, just as example. A good example for using the
      array will be to pass the whole hierarchy, not only the fine
      grid via all, e.g. all=(void *)sc
    */
    all=(void *)malloc(5*sizeof(dvector)+sizeof(ivector));
    /**/
    for(j=0;j<ref_levels;j++){
      /*
       * SELECT the finest grid:
       */
      sctop=scfinest(sc);
      /*
       * SOLVE on the finest grid
       */
      solfem=exmpl_solve(sctop,all);
      /*
       * ESTIMATE using the numerical solution and the data stored in *all
       */
      estimator=exmpl_estimate(sctop,&solfem,all);
      /*
       * MARK: marked is an ivector with num.rows=the number of
       *       simplices; its componenets are nonzero if the simplex
       *       is marked for refinement and 0 if it is not marked for
       *       refinement.
       */
      marked=exmpl_mark(sctop,&estimator,all);
      /*
       *  Refine the grid. this always refines 1 time, but since we
       *  are in a loop, it will refine ref_levels times;
       */
      //      haz_scomplex_print(sctop,0,__FUNCTION__);
      refine(1,sc,&marked);
      /* free */
      haz_scomplex_free(sctop);
      dvec_free(&solfem);
      dvec_free(&estimator);
      ivec_free(&marked);
    }
    free(all);
  }
  /*  MAKE sc to be the finest grid only */
  scfinalize(sc,(INT )1);
  /* WRITE THE OUTPUT MESH FILE:    */
  //  hazw(g->fgrid,sc,0); //last argument is the shift
  if(dim < 4)
    mshw(g->fgrid,sc,0); //  the choice of format should also be given
			 //  as part of a structure.
  /* WRITE THE OUTPUT vtu file for paraview:    */
  if(dim <4) {
    vtu_data vdata;
    vtu_data_init(sc,&vdata);
    vtkw(g->fvtu,&vdata);
    //    vtkw("output/1d_graph.vtu",&vdata); //0,(REAL )1);
    vtu_data_free(&vdata);
  }
  /*FREE: the input grid is freed here, because it has the filenames in it*/
  input_grid_free(g);
  haz_scomplex_free(sc);
  free(sc_all);
  return 0;
}
