/*! \file examples/amr_grids/amr_grids.c
 *
 *  Authors: James Adler, Xiaozhe Hu, and Ludmil Zikatanov
 *           HAZmath (https://hazmath.net)
 *           Created with the help of Claude (Anthropic)
 *
 *  Created 2019/01/09.  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief Generates simplicial grids in 2, 3, 4, ... dimensions.
 *        Supports DGS bisection, uniform Bey, and selective Bey
 *        refinement with conforming closure.
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
INT main(INT argc, char* argv[]) {
  INT i;
  const char* input_file = "input/3d_cube.input";
  if (argc > 1) input_file = argv[1];
  /*
    PARSE THE INPUT.
  */
  config_z config = get_input(input_file);
  input_grid* g = input_grid_alloc();
  config2vars_amr(&config, g);
  free_config(&config);
  input_grid_set_output(g, input_file);
  scomplex** sc_all = generate_initial_grid(g);
  fprintf(stdout, "\nInitial mesh:\nElements = %12lld;\nVertices=%12lld\n", (long long)sc_all[0]->ns, (long long)sc_all[0]->nv); fflush(stdout);
  scomplex* sc = sc_all[0];
  INT ref_levels = g->nref, amr_marking_type = g->mark_type, j, k, kmarked;
  scomplex* sctop = NULL;
  dvector solfem, estimator;
  ivector marked;
  void* all = NULL;
  REAL* xstar = NULL;
  INT nstar, dim = sc->dim, n1 = dim + 1;
  features feat;
  feat.n = 0; feat.nbig = 0; feat.x = NULL; feat.fill = -1e20; feat.fpf = NULL;
  iCSRmat node_ins;
  //NNNNNNNNNNNNNNNNNNNNNNNNNNNN
  fprintf(stdout, "\n****AMR_MARKING_TYPE=%d\n", amr_marking_type);
  if (amr_marking_type == 44) {
    char* data_file = strdup("./try_features.txt");
    feat.fpf = fopen(data_file, "r");
    feat.nbig = sc->dim;
    feat.n = sc->dim;
    feat.fill = -1e20;
    // last argument below is whether to map the simplicial complex to
    // a cube enclosing the data.
    j = features_r(&feat, sc, (INT)1, (REAL)1.1e0);
    free(data_file);
  }
  //NNNNNNNNNNNNNNNN
  if (amr_marking_type == 0) {
    if (sc->ref_type == 21 || sc->ref_type == 22) {
      // Marked Bey + face-Bey/bisection closure
      // ref_type 21: mark odd-indexed simplices (checkerboard test)
      // ref_type 22: mark simplices near the origin (re-entrant corner)
      for (j = 0; j < ref_levels; j++) {
        ivector mark = ivec_create(sc->ns);
        INT nmarked = 0;
        if (sc->ref_type == 21) {
          for (k = 0; k < sc->ns; k++) { mark.val[k] = (k % 2); if (mark.val[k]) nmarked++; }
        } else {
          /* mark simplices whose barycenter is within threshold of origin */
          /* threshold shrinks each level so refinement concentrates */
          REAL threshold = 1.0;
          for (INT jj = 0; jj < j; jj++) threshold *= 0.6;
          for (k = 0; k < sc->ns; k++) {
            REAL dist2 = 0.0;
            for (INT dd = 0; dd < dim; dd++) {
              REAL c = 0.0;
              for (INT vv = 0; vv < n1; vv++)
                c += sc->x[sc->nbig * sc->nodes[n1 * k + vv] + dd];
              c /= (REAL)n1;
              dist2 += c * c;
            }
            mark.val[k] = (dist2 < threshold * threshold) ? 1 : 0;
            if (mark.val[k]) nmarked++;
          }
        }
        fprintf(stdout, "\n%% Marked Bey (lvl %lld): ns=%lld, marked=%lld",
          (long long)j, (long long)sc->ns, (long long)nmarked);
        uniformrefine_marked(sc, &mark);
        ivec_free(&mark);
        sc_vols(sc);
        fprintf(stdout, " -> ns=%lld, nv=%lld", (long long)sc->ns, (long long)sc->nv);
      }
      fprintf(stdout, "\n");
      find_nbr(sc->ns, sc->nv, sc->dim, sc->nodes, sc->nbr);
      sc_vols(sc);
    } else if (sc->ref_type > 10) {
      // uniform (Freudenthal) refinement — works in any dimension
      for (j = 0; j < ref_levels; j++) {
        uniformrefine(sc);
        sc_vols(sc);
      }
      find_nbr(sc->ns, sc->nv, sc->dim, sc->nodes, sc->nbr);
      sc_vols(sc);
    } else {
      // DGS bisection refinement
      refine(ref_levels, sc, NULL);
    }
  } else if (amr_marking_type == 44) {
    node_ins = icsr_create(0, 0, 0);
    nstar = feat.nf;
    xstar = feat.x;
    //
    INT max_nodes = (INT)MAX_NODES_PER_SIMPLEX;
    if (max_nodes <= 0) max_nodes = 1;
    for (j = 0; j < ref_levels; j++) {
      sctop = scfinest(sc);
      /* MARK: marked is an ivector with num.rows=the number of
       *       simplices; its componenets are nonzero if the simplex
       *       is marked for refinement and 0 if it is not marked for
       *       refinement.
       */
      marked = mark_around_pts(sctop, sc, nstar, xstar, &node_ins, (const INT)max_nodes);
      kmarked = 0;
      for (k = 0; k < marked.row; ++k)
        if (marked.val[k]) kmarked++;
      fprintf(stdout, "\n|lvl=%2lld|simplices[%2lld:%2lld]=%12lld|simplices[%2lld]=%12lld|", (long long int)j, 0LL, (long long int)j, (long long int)sc->ns, (long long int)j, (long long int)sctop->ns); fflush(stdout);
      refine(1, sc, &marked);
      if (!kmarked) {
        fprintf(stdout, "\nthere were no simplices containing > %lld points. Exiting", (long long)MAX_NODES_PER_SIMPLEX);
        ivec_free(&marked);
        haz_scomplex_free(sctop);
        break;
      }
      ivec_free(&marked);
      haz_scomplex_free(sctop);
    }
    fprintf(stdout, "\n");
    icsr_free(&node_ins);
    ivec_free(&marked);
    free(feat.x);
  } else if (amr_marking_type == 33) {
    REAL h = 0.05;  // step distance of points
    REAL threshold = h; // threshold for close to the points or not
    //
    nstar = g->num_refine_points;
    xstar = g->data_refine_points;
    //
    for (j = 0; j < ref_levels; j++) {
      /*
       * SELECT the finest grid:
       */
      sctop = scfinest(sc);
      /* MARK: marked is an ivector with num.rows=the number of
       *       simplices; its componenets are nonzero if the simplex
       *       is marked for refinement and 0 if it is not marked for
       *       refinement.
       */
      marked = mark_near_points(sctop, nstar, xstar, threshold);
      fprintf(stdout, "\n|lvl=%2lld|simplices[%2lld:%2lld]=%12lld|simplices[%2lld]=%12lld|", (long long int)j, 0LL, (long long int)j, (long long int)sc->ns, (long long int)j, (long long int)sctop->ns); fflush(stdout);
      refine(1, sc, &marked);
      /* free */
      ivec_free(&marked);
      haz_scomplex_free(sctop);
    }
    ivec_free(&marked);
    //free(xstar);
  } else if (amr_marking_type == 34) {
    // Refine near specified points using Bey + face-Bey/bisection closure.
    // Mark simplices whose barycenter is within a shrinking threshold
    // of any specified point. Threshold shrinks by 0.6x per level.
    nstar = g->num_refine_points;
    xstar = g->data_refine_points;
    REAL threshold = 1.0;
    for (j = 0; j < ref_levels; j++) {
      marked = ivec_create(sc->ns);
      INT nmarked = 0;
      for (k = 0; k < sc->ns; k++) {
        REAL mindist2 = 1e30;
        for (INT s = 0; s < nstar; s++) {
          REAL dist2 = 0.0;
          for (INT dd = 0; dd < dim; dd++) {
            REAL c = 0.0;
            for (INT vv = 0; vv < n1; vv++)
              c += sc->x[sc->nbig * sc->nodes[n1 * k + vv] + dd];
            c /= (REAL)n1;
            REAL d = c - xstar[s * dim + dd];
            dist2 += d * d;
          }
          if (dist2 < mindist2) mindist2 = dist2;
        }
        marked.val[k] = (mindist2 < threshold * threshold) ? 1 : 0;
        if (marked.val[k]) nmarked++;
      }
      fprintf(stdout, "\n%% Marked Bey (lvl %lld): ns=%lld, marked=%lld",
        (long long)j, (long long)sc->ns, (long long)nmarked);
      uniformrefine_marked(sc, &marked);
      ivec_free(&marked);
      sc_vols(sc);
      fprintf(stdout, " -> ns=%lld, nv=%lld", (long long)sc->ns, (long long)sc->nv);
      threshold *= 0.6;
    }
    fprintf(stdout, "\n");
    find_nbr(sc->ns, sc->nv, sc->dim, sc->nodes, sc->nbr);
    sc_vols(sc);
  } else {
    /*
      Use "all" here can pass data around. Below we make 4 dvectors
      and one ivector, just as example. A good example for using the
      array will be to pass the whole hierarchy, not only the fine
      grid via all, e.g. all=(void *)sc
    */
    all = (void*)malloc(5 * sizeof(dvector) + sizeof(ivector));
    /**/
    for (j = 0; j < ref_levels; j++) {
      /*
       * SELECT the finest grid:
       */
      sctop = scfinest(sc);
      /*
       * SOLVE on the finest grid
       */
      solfem = exmpl_solve(sctop, all);
      /*
       * ESTIMATE using the numerical solution and the data stored in *all
       */
      estimator = exmpl_estimate(sctop, &solfem, all);
      /*
       * MARK: marked is an ivector with num.rows=the number of
       *       simplices; its componenets are nonzero if the simplex
       *       is marked for refinement and 0 if it is not marked for
       *       refinement.
       */
      marked = exmpl_mark(sctop, &estimator, all);
      /*
       *  Refine the grid. this always refines 1 time, but since we
       *  are in a loop, it will refine ref_levels times;
       */
      //      haz_scomplex_print(sctop,0,__FUNCTION__);
      refine(1, sc, &marked);
      /* free */
      haz_scomplex_free(sctop);
      dvec_free(&solfem);
      dvec_free(&estimator);
      ivec_free(&marked);
    }
    free(all);
  }
  /*  MAKE sc to be the finest grid only */
  scfinalize(sc, NULL, (INT)1);
  /* conformity check (expensive for large meshes — off by default) */
  INT do_conformity_check = 0;
  if (do_conformity_check) {
    INT nerr = sc_conformity_check(sc);
    if (nerr)
      fprintf(stderr, "\n%% FAIL: non-conforming mesh (%lld bad facets)\n", (long long)nerr);
    else
      fprintf(stdout, "\n%% conformity check PASSED (ns=%lld, nv=%lld)\n", (long long)sc->ns, (long long)sc->nv);
  }
  /* WRITE OUTPUT FILES: .haz, .msh, .vtu
   * g->fgrid is the base (e.g. "output/3d_fichera_rl3_rt20")
   * g->fvtu already has .vtu appended */
  {
    char fname[MAXFILENAMESIZE];
    /* hazw(fname, sc, 0); -- commented out, use sc_write_gmsh */
    snprintf(fname, sizeof(fname), "%s.msh", g->fgrid);
    sc_write_gmsh(fname, sc, 0);
    if (dim < 4) {
      vtu_data vdata;
      vtu_data_init(sc, &vdata);
      sc_write_vtk(g->fvtu, &vdata);
      vtu_data_free(&vdata);
    }
  }
  /*FREE: the input grid is freed here, because it has the filenames in it*/
  input_grid_free(g);
  haz_scomplex_free(sc);
  free(sc_all);
  return 0;
}
