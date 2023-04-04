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
  If set to something less than 2 it will refine non-stop until ref_levels is
  reached if there is at least one simplex containing at least one point. It is
  the maximum allowed number of features (nodes) per element. Any element
  containing more than this number of features is refined.

*/
/* 
 * refinement type: .gt. 10 is uniform refinement and .le. 10
 *                  (typically 0) is the newest vertex bisection
*/
#ifndef MAX_NODES_PER_SIMPLEX
#define MAX_NODES_PER_SIMPLEX  1
#endif
/**/
#ifndef OUTER_SPATIAL_DIMENSION
#define OUTER_SPATIAL_DIMENSION 3
#endif
/**/
#ifndef REFINEMENT_LEVELS
#define REFINEMENT_LEVELS 150
#endif
/**/
/*********************************************************************/
typedef struct /* n-homogenous simplicial complex */
{
  INT dimbig;     /* the dimension of the space in which SC is embedded */
  INT dim;        /* the dimension of SC */
  INT nv;         /* number of 0-dimensional simplices */
  INT nvadd;      /* number of 0-dimensional simplices added on every segment */
  INT nseg;       /* number of segments */
  INT *seg;       /* nv boundary codes for vertices */
  INT *divisions; /*divisions per segment */
  REAL *xv; /*(nv times dimbig) array to hold the coordinates of vertices */
  REAL *pt_thickness; /* points attribute */
  REAL *seg_radius;   /* segments attribute */
  char *fv_coords;   /*input_file: coords of bifurcations*/
  char *fseg;       /*input_file: segments definition */
  char *fdivisions; /*input_file: divisions per segment*/
  char *fvtmp_coords; /*input_file: coordinates of all points (biffurcations or not */
  char *fpt_thickness;/*input_file: segment thickness*/
  char *fseg_radius;/*input_file: radius */
  char *fvtu_3d;  /*output_file: for the 3d grid in vtu*/
  char *fvtu_1d;  /*output_file: for the 1d grid on vtu */
} data_1d;
/*********************************************************************/
#include "supporting_xd_1d.h"
/*********************************************************************************************/
static INT init_pts(const INT dim, const INT npts, REAL *pts, scomplex *sc, const REAL scale)
{
  INT k = 0, i, j;
  cube2simp *c2s = cube2simplex(dim);  // now we have the vertices of the unit cube in bits
  REAL *vc = calloc(4 * dim * c2s->nvcube, sizeof(REAL));
  REAL *xmintmp = vc;                                 // maps to [0...0]
  REAL *xmaxtmp = xmintmp + dim * (c2s->nvcube - 1);  // last vertex
  REAL *xmin = xmaxtmp + dim * (c2s->nvcube - 1);     // last vertex
  REAL *xmax = xmin + dim * (c2s->nvcube - 1);        // last vertex
  INT kdimi;
  for (i = 0; i < dim; i++) {
    xmax[i] = pts[i];
    xmin[i] = xmax[i];
    kdimi = dim + i;
    for (k = 1; k < npts; k++) {
      if (pts[kdimi] > xmax[i]) {
        xmax[i] = pts[kdimi];
      }
      if (pts[kdimi] < xmin[i]) {
        xmin[i] = pts[kdimi];
      }
      kdimi += dim;
    }
  }
  fprintf(stdout,"\n\ndim=%d;scale=%.16e\n",dim,scale);fflush(stdout);
  print_full_mat(1,dim,xmin,"xmin");
  print_full_mat(1,dim,xmax,"xmax");
  for (i = 0; i < dim; i++) {
    xmaxtmp[i] = xmax[i] + (scale - 1e0) * (xmax[i] - xmin[i]);
    xmintmp[i] = xmin[i] - (scale - 1e0) * (xmax[i] - xmin[i]);
  }
  for (j = 1; j < c2s->nvcube - 1; j++) {
    for (i = 0; i < dim; i++) {
      vc[j * dim + i] =
          xmintmp[i] + (xmaxtmp[i] - xmintmp[i]) * (c2s->bits[dim * j + i]);
    }
  }
  mapit(sc, vc);
  free(vc);
  cube2simp_free(c2s);
  return 0;
}
/****************************************************************************************/
static void special_1d(scomplex *sc, data_1d *g, dvector *seg_r) {
  // construct 1-homogenous simplicial complex embedded in 2d or 3d
  INT i, j, k, l, m, dimbig;
  dimbig = g->dimbig;
  REAL *xvtmp = g->xv + g->nv * dimbig;
  REAL r = -1e20;
  INT i1, i2, ko, bego, endo, ptrn, ptre;
  bego = 0, ptrn = 0;
  ptre = 0;
  for (i = 0; i < g->nseg; ++i) {
    l = g->divisions[i] - 2;
    i1 = g->seg[2 * i];
    i2 = g->seg[2 * i + 1];
    r = g->seg_radius[i];
    if (!l) {
      endo = bego + 1;
      bego = endo + 1;
      sc->nodes[2 * ptre] = i1;
      sc->nodes[2 * ptre + 1] = i2;
      seg_r->val[ptre] = r;
      ptre++;
      continue;
    }
    endo = bego + l + 1;
    sc->nodes[2 * ptre] = i1;
    sc->nodes[2 * (ptre + l) + 1] = i2;
    seg_r->val[ptre + l] = r;
    k = ptrn;
    ko = bego + 1;
    for (j = ptre; j < (ptre + l); ++j) {
      sc->nodes[2 * j + 1] = k + g->nv;
      sc->nodes[2 * (j + 1)] = k + g->nv;
      // fprintf(stdout,"\nt123(%d,1:2)=[%d
      // %d];",j+1,sc->nodes[2*j]+1,sc->nodes[2*j+1]+1);
      seg_r->val[j] = r;
      for (m = 0; m < dimbig; ++m) {
        xvtmp[dimbig * k + m] = xvtmp[dimbig * ko + m];
      }
      k++;
      ko++;
    }
    // fprintf(stdout,"\nt123(%d,1:2)=[%d
    // %d];",ptre+l+1,sc->nodes[2*(ptre+l)]+1,sc->nodes[2*(ptre+l)+1]+1);
    bego = endo + 1;
    ptrn += l;
    ptre += (l + 1);
  }
  sc->nv = ptrn + g->nv;
  sc->ns = ptre;
  seg_r->val = realloc(seg_r->val, (sc->ns) * sizeof(REAL));
  sc->nodes = realloc(sc->nodes, (sc->n + 1) * (sc->ns) * sizeof(INT));
  free(sc->x);
  g->xv = realloc(g->xv, sc->nv * dimbig * sizeof(REAL));  //
  sc->x = g->xv;
  return;
}
INT main(INT argc, char *argv[]) {
  INT j, k;
  data_1d g;
  //init 1d struct
  data_1d_init(&g);
  INT dim=g.dimbig;// 
  getdata_1d(&g);
  ////////////// read all 1d data.
  //  form the complex:
  fprintf(stdout, "\nnv=%d,nvadd=%d,nseg=%d\n", g.nv, g.nvadd, g.nseg);
  fflush(stdout);
  scomplex *sc_dim = haz_scomplex_init(g.dim, g.nvadd, g.nvadd, g.dimbig);
  dvector seg_r = dvec_create(sc_dim->ns);  // seg_radius;
  //
  special_1d(sc_dim, &g, &seg_r);
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
  INT ref_levels=REFINEMENT_LEVELS;
  //
  fprintf(stdout,"\nMeshing in dimension=%lld ...",(long long )dim);
  //  clock_t clk_mesh_start = clock();
  scomplex **sc_all=mesh_cube_init(dim,(INT )1,(INT )0);
  scomplex *sc_dimbig = sc_all[0];
  //haz_scomplex_print(sc_dimbig,0,"ZZZ");
  fprintf(stdout,"\nElements = Simplexes = %12lld;\nDoF      = Vertices  = %12lld\n",(long long )sc_dimbig->ns,(long long )sc_dimbig->nv); fflush(stdout);
  /**/
  //  scfinalize(sc_dimbig,(INT )1);  
  //  sc_vols(sc_dimbig);

  /* input_grid *g3d = parse_input_grid(fp); */
  /* fclose(fp); */
  /* scomplex **sc_all = generate_initial_grid(g3d); */
  /* fprintf(stdout, "\nInitial mesh:\nElements = %12lld;\nVertices=%12lld\n", */
  /*         (long long)sc_all[0]->ns, (long long)sc_all[0]->nv); */
  /* fflush(stdout); */
  scomplex *sctop = NULL;
  ivector marked;
  void *all = NULL;
  INT kmarked;
  INT nstar = sc_dim->nv;
  REAL *xstar = calloc(nstar * dim, sizeof(REAL));
  memcpy(xstar, sc_dim->x, nstar * dim * sizeof(REAL));
  if (sc_dim) haz_scomplex_free(sc_dim);
  iCSRmat node_ins;
  init_pts(dim, nstar, xstar, sc_dimbig, (REAL)1.1);
  node_ins = icsr_create(0, 0, 0);
  //
  INT max_nodes = (INT)MAX_NODES_PER_SIMPLEX;
  if (max_nodes <= 0) max_nodes = 1;
  for (j = 0; j < ref_levels; j++) {
    sctop = scfinest(sc_dimbig);
    /* MARK: marked is an ivector with num.rows=the number of
     *       simplices; its componenets are nonzero if the simplex
     *       is marked for refinement and 0 if it is not marked for
     *       refinement.
     */
    marked = mark_around_pts(sctop, sc_dimbig, nstar, xstar, &node_ins,
                             (const INT)max_nodes);
    kmarked = 0;
    for (k = 0; k < marked.row; ++k)
      if (marked.val[k]) kmarked++;
    fprintf(stdout,
            "\n|lvl=%2lld|simplices[%2lld:%2lld]=%12lld|simplices[%2lld]=%"
            "12lld|vertices[%2lld]=%12lld|kmarked=%lld",
            (long long int)j, 0LL, (long long int)j,
            (long long int)sc_dimbig->ns, (long long int)j,
            (long long int)sctop->ns, (long long int)j,
            (long long int)sc_dimbig->nv, (long long int)kmarked);
    fflush(stdout);
    refine(1, sc_dimbig, &marked);
    if (!kmarked) {
      fprintf(stdout,
              "\nthere were no simplices containing > %lld points. Exiting",
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
  vtu_data_init(sc_dimbig, &vdata);
  vtkw(g.fvtu_3d,&vdata);
  vtu_data_free(&vdata);
  /*FREE: the input grid is freed here, because it has the filenames in it*/
  //  input_grid_free(g3d);
  //  haz_scomplex_print(sc_dimbig,0,__FUNCTION__);
  data_1d_free(&g);
  haz_scomplex_free(sc_dimbig);
  free(sc_all);
  //  clock_t clk_mesh_end = clock();
  return 0;
}
