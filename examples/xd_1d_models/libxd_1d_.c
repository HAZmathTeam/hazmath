/*! \file examples/amr_grids/xd_1d.c
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
  If set to something less than 2 it will refine non-stop until ref_levels is
  reached if there is at least one simplex containing at least one point. It is
  the maximum allowed number of features (nodes) per element. Any element
  containing more than this number of features is refined.

*/
/* 
 * refinement type: .gt. 10 is uniform refinement and .le. 10
 *                  (typically 0) is the newest vertex bisection
*/
#include "definitions_xd_1d.h"
#include "supporting_xd_1d.h"
#include "solver_xd_1d.h"
///////////////////////////////////////////////////////////////////////////////
void xd_1d_lib(const INT dimbig,const INT max_nodes_in,const INT ref_levels_in, \
	       const char *idir,const char *odir)
{
  INT j, k;
  data_1d g;
  //init 1d struct
  /* fprintf(stdout,"\n%%%%%%In %s:\n%%%%%%DIMBIG=%lld\nMAX_NODES=%lld\nREF_LEVELS=%lld\nINPUT=%s\nOUTPUT=%s", \ */
  /* 	  __FUNCTION__, */
  /* 	  (long long)dimbig,					\ */
  /* 	  (long long)max_nodes_in,				\ */
  /* 	  (long long)ref_levels_in,				\ */
  /* 	  idir,							\ */
  /* 	  odir); */
  data_1d_init(dimbig,idir,odir,&g);
  INT dim=g.dimbig;// 
  getdata_1d(&g);
  ////////////// read all 1d data.
  /* fprintf(stdout, "\nnv=%d,nvadd=%d,nseg=%d\n", g.nv, g.nvadd, g.nseg);  fflush(stdout); */
  //  form the complex:
  scomplex *sc_dim = haz_scomplex_init(g.dim, g.nvadd, g.nvadd, g.dimbig);
  dvector seg_r = dvec_create(sc_dim->ns);  // segment_radius;
  clock_t clk_1dmesh_start = clock();
  special_1d(sc_dim, &g, &seg_r);
  clock_t clk_1dmesh_end = clock();
  //
  /* fp=fopen("coords.txt","w"); */
  /* for(k=0;k<sc_dim->nv;k++){ */
  /*   for(j=0;j<sc_dim->nbig;j++){ */
  /*     fprintf(fp," %23.16e",sc_dim->x[k*sc_dim->nbig+j]); */
  /*   } */
  /*   fprintf(fp,"\n"); */
  /* } */
  /* fclose(fp); */
  // mshw("1d_graph.msh",sc_dim,0); //
  /* WRITE THE OUTPUT vtu file for paraview:    */
  vtu_data vdata;
  vtu_data_init(sc_dim, &vdata);
  vdata.dcell = realloc(vdata.dcell, (vdata.ndcell + 1) * sizeof(REAL *));
  vdata.names_dcell=realloc(vdata.names_dcell, (vdata.ndcell + 1) * sizeof(char *));
  vdata.dcell[vdata.ndcell] = seg_r.val;
  vdata.names_dcell[vdata.ndcell] = strdup("thickness");
  vdata.ndcell++;                         // increase with one;
  // write the 1d vtu:  
  vtkw(g.fvtu_1d, &vdata);
  //free the vdata with variable names etc
  vtu_data_free(&vdata);
  // free 1d data;
  dvec_free(&seg_r);
  //
  //  fprintf(stdout,"\n%%%%Meshing in dimension=%lld ...",(long long )dim);
  clock_t clk_3dmesh_start = clock();
  scomplex **sc_all=mesh_cube_init(dim,(INT )1,(INT )0); //number of divisions =1
  scomplex *sc_dimbig = sc_all[0];
  scomplex *sctop = NULL;
  ivector marked;
  void *all = NULL;
  INT kmarked;
  INT nstar = sc_dim->nv;
  REAL *xstar = calloc(nstar * dim, sizeof(REAL));
  memcpy(xstar, sc_dim->x, nstar * dim * sizeof(REAL));
  if (sc_dim) haz_scomplex_free(sc_dim);
  iCSRmat node_ins;
  init_pts(dim, nstar, xstar, sc_dimbig, (REAL)1.1);// scale by 1.1 to
						    // contain the
						    // network's
						    // interpolationm
						    // radius.
  node_ins = icsr_create(0, 0, 0);// sparse matrix giving the correspondence: (1d-node)->(3d-simplex)
  //
  INT max_nodes=max_nodes_in,ref_levels=ref_levels_in;
  if (max_nodes_in <= 0) max_nodes = 1;
  for (j = 0; j < ref_levels; j++) {
    sctop = scfinest(sc_dimbig);
    /* MARK: marked is an ivector with num.rows=the number of
     *       simplices; its componenets are nonzero if the simplex
     *       is marked for refinement and 0 if it is not marked for
     *       refinement.
     */
    marked = mark_around_pts(sctop, sc_dimbig, nstar, xstar, &node_ins, (const INT)max_nodes);
    kmarked = 0;
    for (k = 0; k < marked.row; ++k)
      if (marked.val[k]) kmarked++;
    fprintf(stdout,
            "\n|lvl=%2lld|simplices[%2lld:%2lld]=%12lld|simplices[%2lld]=%"
            "12lld|vertices[%2lld]=%12lld|marked=%lld",
            (long long int)j, 0LL, (long long int)j,
            (long long int)sc_dimbig->ns, (long long int)j,
            (long long int)sctop->ns, (long long int)j,
            (long long int)sc_dimbig->nv, (long long int)kmarked);
    fflush(stdout);
    refine(1, sc_dimbig, &marked);
    if (!kmarked) {
      fprintf(stdout,
              "\nDone:There were no simplices containing > %lld bif. points from the 1d network.\n",
              (long long)(max_nodes));
      ivec_free(&marked);
      haz_scomplex_free(sctop);
      break;
    }
  }
  icsr_free(&node_ins);
  free(xstar);
  free(all);
  /*  MAKE sc to be the finest grid only */
  //  icsr_write_icoo("parent_v.zcoo",sc_dimbig->parent_v);
  scfinalize(sc_dimbig, (INT)1);
  //
  clock_t clk_3dmesh_end = clock();
  //
  vtu_data_init(sc_dimbig, &vdata);
  vtkw(g.fvtu_3d,&vdata);
  vtu_data_free(&vdata);
  /*FREE: the input grid is freed here, because it has the filenames in it*/
  //  input_grid_free(g3d);
  //  haz_scomplex_print(sc_dimbig,0,__FUNCTION__);
  fprintf(stdout,"\n%%%%%%CPUtime(1d-mesh)     = %10.3f sec\n%%%%%%CPUtime(%dd-mesh) = %10.3f sec\n", \
	  (REAL ) (clk_1dmesh_end - clk_1dmesh_start)/CLOCKS_PER_SEC,	g.dimbig,
	  (REAL ) (clk_3dmesh_end - clk_3dmesh_start)/CLOCKS_PER_SEC);
  ////////////////////////////////////////////////////////////////////////////////
  data_1d_free(&g);
  haz_scomplex_free(sc_dimbig);
  free(sc_all);
  return;
}
