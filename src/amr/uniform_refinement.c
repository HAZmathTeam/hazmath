/*! \file src/amr/uniform_refinement.c
 *
 *  Authors: James Adler, Xiaozhe Hu, and Ludmil Zikatanov
 *           HAZmath (https://hazmath.net)
 *           Created with the help of Claude (Anthropic)
 *
 *  Created 20170715.  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note Uniform and selective Bey (Freudenthal) refinement in any
 *        spatial dimension: get_edge_nd, uniformrefine,
 *        uniformrefine_marked (selective Bey + face-Bey/bisection closure).
 *        Original 2D/3D routines by Ludmil and Yuwen (20210604).
 *        N-dimensional generalization 20260322.
 */
#include "hazmath.h"
/********************************************************************************/
/*!
 * \fn void get_edge2d(iCSRmat *e2v, iCSRmat *el2e, scomplex *sc)
 *
 * \brief Returns the edge-to-vertex matrix e2v and the
 *        element-to-edge vertex el2e in 2d.  The input sc is the
 *        pointer to the scomplex with element-to-vertex matrix and
 *        coordinates of grid vertices. e2v is an incidence matrix of
 *        ne by nv, where ne is the number of edges and nv the number
 *        of vertices. el2e is an incidence matrix of ns by ne, where
 *        ns is the number of simplexes. The (i,j)-entry of e2v is 1
 *        if the i-th edge contains the j-th vertex. The (i,j)-entry
 *        of el2e is 1 if the i-th simplex contains the j-th edge. The
 *        edges are ascendingly ordered in a global lexicographical
 *        way.  The order of elements for el2e->JA corresponding to
 *        the i-th row is given in local lexicographic order, (0, 1),
 *        (0, 2), (1, 2).
 *
 * \param e2v  Pointer to iCSRmat edge-to-vertex matrix
 * \param el2e  Pointer to iCSRmat element-to-edge matrix
 *
 *  \note Yuwen 20210606.
 */
void get_edge2d(iCSRmat* e2v, iCSRmat* el2e, scomplex* sc) {
  if (sc->dim != 2) {
    fprintf(stderr, "%% *** ERROR in %s: the dimension = %lld and is not 3  !\n", __FUNCTION__, (long long)sc->dim);
    exit(2);
  }
  INT nv = sc->nv, ns = sc->ns;
  INT ns2 = 2 * ns, ns3 = 3 * ns, i, j, k, ih, n0, n1, n2, tmp;
  ivector ii = ivec_create(ns3);
  ivector jj = ivec_create(ns3);
  iCSRmat v2e, el2e0;
  iCSRmat U = icsr_create(0, 0, 0);
  iCSRmat el2v_csr = icsr_create(ns, nv, ns3);
  el2v_csr.IA[0] = 0;
  INT* el2v = sc->nodes;
  /* The pair  (ii,jj) encodes three edges in each element in ascending lexicographic order */
  for (k = 0; k < ns; k++) {
    /* sorting the local three vertices in ascending order */
    n0 = 3 * k; n1 = n0 + 1; n2 = n1 + 1;
    if (el2v[n0] > el2v[n1]) {tmp = n0; n0 = n1; n1 = tmp;}
    if (el2v[n0] > el2v[n2]) {tmp = n0; n0 = n2; n2 = tmp;}
    if (el2v[n1] > el2v[n2]) {tmp = n1; n1 = n2; n2 = tmp;}

    ii.val[k] = el2v[n0]; jj.val[k] = el2v[n1];
    ii.val[k + ns] = el2v[n0]; jj.val[k + ns] = el2v[n2];
    ii.val[k + ns2] = el2v[n1]; jj.val[k + ns2] = el2v[n2];

    el2v_csr.IA[k + 1] = el2v_csr.IA[k] + 3;
  }
  memcpy(el2v_csr.JA, el2v, (sc->dim + 1) * ns * sizeof(INT));
  for (i = 0; i < ns3; i++) {el2v_csr.val[i] = 1;}
  icsr_uniqueij(&U, &ii, &jj);
  ivec_free(&ii); ivec_free(&jj);
  /* Form the edge to vertex csr matrix */
  INT ne = U.nnz, counter = 0;
  e2v->row = ne;
  e2v->col = nv;
  e2v->nnz = 2 * ne;

  INT* ia = U.IA, *ja = U.JA;
  e2v->JA = (INT*)malloc(2 * ne * sizeof(INT));
  e2v->val = (INT*)malloc(2 * ne * sizeof(INT));
  for (i = 0; i < 2 * ne; i++) {e2v->val[i] = 1;}
  for (i = 0; i < nv; i++) {
    ih = ia[i + 1] - ia[i];
    for (j = 0; j < ih; j++) {
      e2v->JA[counter + 2 * j] = i;
      e2v->JA[counter + 2 * j + 1] = ja[ia[i] + j];
    }
    counter = counter + 2 * ih;
  }

  e2v->IA = (INT*)malloc((ne + 1) * sizeof(INT));
  e2v->IA[0] = 0;
  for (i = 0; i < ne; i++) {
    e2v->IA[i + 1] = e2v->IA[i] + 2;
  }
  icsr_free(&U);
  icsr_trans(e2v, &v2e);

  /* The i-th row of el2e0 enumerates all edges that shares at least one vertex with the i-th element*/
  icsr_mxm(&el2v_csr, &v2e, &el2e0);
  icsr_free(&v2e);
  icsr_free(&el2v_csr);

  ia = el2e0.IA, ja = el2e0.JA;
  INT* val = el2e0.val;

  /* The i-th row of el2e enumerates all edges that shares two vertices with the i-th element */
  el2e->row = ns; el2e->col = ne; el2e->nnz = 3 * ns;
  el2e->IA = (INT*)malloc((ns + 1) * sizeof(INT));
  el2e->JA = (INT*)malloc(3 * ns * sizeof(INT));
  el2e->val = NULL;
  el2e->IA[0] = 0;
  counter = 0;
  for (i = 0; i < ns; i++) {
    el2e->IA[i + 1] = el2e->IA[i] + 3;
    ih = ia[i + 1] - ia[i];
    for (j = 0; j < ih; j++) {
      if (val[ia[i] + j] == 2) {
        el2e->JA[counter] = ja[ia[i] + j];
        counter = counter + 1;
      }
    }
  }

  /* sorting each row of el2e such that the order of three edges is locally ascending lexicographic*/
  INT i3, ej1, ej2, *jtmp = (INT*)calloc(3, sizeof(INT));
  for (i = 0; i < ns; i++) {
    i3 = 3 * i;
    memcpy(jtmp, &el2e->JA[i3], 3 * sizeof(INT));
    for (j = 0; j < 3; j++) {
      ej1 = e2v->JA[2 * jtmp[j]];
      ej2 = e2v->JA[2 * jtmp[j] + 1];
      if (((sc->nodes[i3] == ej1) && (sc->nodes[i3 + 1] == ej2))
         || ((sc->nodes[i3] == ej2) && (sc->nodes[i3 + 1] == ej1)))
      {el2e->JA[i3] = jtmp[j]; continue;}

      if (((sc->nodes[i3] == ej1) && (sc->nodes[i3 + 2] == ej2))
         || ((sc->nodes[i3] == ej2) && (sc->nodes[i3 + 2] == ej1)))
      {el2e->JA[i3 + 1] = jtmp[j]; continue;}

      if (((sc->nodes[i3 + 1] == ej1) && (sc->nodes[i3 + 2] == ej2))
         || ((sc->nodes[i3 + 1] == ej2) && (sc->nodes[i3 + 2] == ej1)))
      {el2e->JA[i3 + 2] = jtmp[j]; continue;}
    }
  }
  free(jtmp);
  icsr_free(&el2e0);
}
/*******************************************************************************************/
/*!
 * \fn void get_edge3d(iCSRmat *e2v, iCSRmat *el2e, scomplex *sc)
 *
 * \brief Returns the edge-to-vertex matrix e2v and the element-to-edge vertex el2e in 3d.
 *        The input sc is the pointer to the scomplex with element-to-vertex matrix and
 *        coordinates of grid vertices. e2v is an incidence matrix of ne by nv, where ne
 *        is the number of edges and nv the number of vertices. el2e is an incidence matrix
 *        of ns by ne, where ns is the number of simplexes. The (i,j)-entry of e2v is 1 if
 *        the i-th edge contains the j-th vertex. The (i,j)-entry of el2e is 1 if the i-th simplex
 *        contains the j-th edge. The edges are ascendingly ordered in a global lexicographical way.
 *        The order of elements for el2e->JA corresponding to the i-th row
 *        is given in local lexicographic order, (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3).
 *
 * \param e2v  Pointer to iCSRmat edge-to-vertex matrix
 * \param el2e  Pointer to iCSRmat element-to-edge matrix
 *
 *  \note Yuwen 20210606.
 */
void get_edge3d(iCSRmat* e2v, iCSRmat* el2e, scomplex* sc) {
  if (sc->dim != 3) {
    fprintf(stderr, "%% *** ERROR in %s: the dimension = %lld and is not 3  !\n", __FUNCTION__, (long long)sc->dim);
    exit(3);
  }
  INT nv = sc->nv, ns = sc->ns;
  INT ns2 = 2 * ns, ns3 = 3 * ns, ns4 = 4 * ns, ns5 = 5 * ns, i, j, k, ih, n0, n1, n2, n3, tmp;
  ivector ii = ivec_create(6 * ns), jj = ivec_create(6 * ns);
  iCSRmat U = icsr_create(0, 0, 0), el2v_csr = icsr_create(ns, nv, 4 * ns), v2e, el2e0;
  INT* el2v = sc->nodes;
  el2v_csr.IA[0] = 0;
  /* The pair  (ii,jj) encodes 6 edges in each element in ascending lexicographic order */
  for (k = 0; k < ns; k++) {
    n0 = 4 * k; n1 = n0 + 1; n2 = n1 + 1; n3 = n2 + 1;
    if (el2v[n0] > el2v[n1]) {tmp = n0; n0 = n1; n1 = tmp;}
    if (el2v[n0] > el2v[n2]) {tmp = n0; n0 = n2; n2 = tmp;}
    if (el2v[n0] > el2v[n3]) {tmp = n0; n0 = n3; n3 = tmp;}
    if (el2v[n1] > el2v[n2]) {tmp = n1; n1 = n2; n2 = tmp;}
    if (el2v[n1] > el2v[n3]) {tmp = n1; n1 = n3; n3 = tmp;}
    if (el2v[n2] > el2v[n3]) {tmp = n2; n2 = n3; n3 = tmp;}
    /*printf("%d\t%d\t%d\t%d\n",n0,n1,n2,n3);*/
    ii.val[k] = el2v[n0]; jj.val[k] = el2v[n1];
    ii.val[k + ns] = el2v[n0]; jj.val[k + ns] = el2v[n2];
    ii.val[k + ns2] = el2v[n0]; jj.val[k + ns2] = el2v[n3];
    ii.val[k + ns3] = el2v[n1]; jj.val[k + ns3] = el2v[n2];
    ii.val[k + ns4] = el2v[n1]; jj.val[k + ns4] = el2v[n3];
    ii.val[k + ns5] = el2v[n2]; jj.val[k + ns5] = el2v[n3];

    el2v_csr.IA[k + 1] = el2v_csr.IA[k] + 4;
  }
  /*el2v_csr.JA = el2v;*/
  memcpy(el2v_csr.JA, el2v, 4 * ns * sizeof(INT));
  for (i = 0; i < 4 * ns; i++) {el2v_csr.val[i] = 1;}
  icsr_uniqueij(&U, &ii, &jj);
  ivec_free(&ii);
  ivec_free(&jj);

  INT ne = U.nnz, counter = 0;
  e2v->row = ne;
  e2v->col = nv;
  e2v->nnz = 2 * ne;
  /* form the edge to vertex csr matrix */
  INT* ia = U.IA, *ja = U.JA;
  e2v->JA = (INT*)calloc(2 * ne, sizeof(INT));
  e2v->val = (INT*)calloc(2 * ne, sizeof(INT)); /*memset(e2v->val,1,2*ne*sizeof(INT));*/
  for (i = 0; i < 2 * ne; i++) {e2v->val[i] = 1;}
  for (i = 0; i < nv; i++) {
    ih = ia[i + 1] - ia[i];
    for (j = 0; j < ih; j++) {
      e2v->JA[counter + 2 * j] = i;
      e2v->JA[counter + 2 * j + 1] = ja[ia[i] + j];
    }
    counter = counter + 2 * ih;
  }

  e2v->IA = (INT*)calloc((ne + 1), sizeof(INT));
  e2v->IA[0] = 0;
  for (i = 0; i < ne; i++) {
    e2v->IA[i + 1] = e2v->IA[i] + 2;
  }
  icsr_free(&U);
  icsr_trans(e2v, &v2e);
  /* The i-th row of el2e0 enumerates all edges that shares at least one vertex with the i-th element*/
  icsr_mxm(&el2v_csr, &v2e, &el2e0);
  icsr_free(&v2e);
  icsr_free(&el2v_csr);

  ia = el2e0.IA, ja = el2e0.JA;
  INT* val = el2e0.val;
  /* The i-th row of el2e enumerates all edges that shares two vertices with the i-th element*/
  el2e->row = ns; el2e->col = ne; el2e->nnz = 6 * ns;
  el2e->IA = (INT*)calloc((ns + 1), sizeof(INT));
  el2e->JA = (INT*)calloc(6 * ns, sizeof(INT));
  el2e->val = NULL;
  el2e->IA[0] = 0;
  counter = 0;
  for (i = 0; i < ns; i++) {
    el2e->IA[i + 1] = el2e->IA[i] + 6;
    ih = ia[i + 1] - ia[i];
    for (j = 0; j < ih; j++) {
      if (val[ia[i] + j] == 2) {
        el2e->JA[counter] = ja[ia[i] + j];
        counter = counter + 1;
      }
    }
  }

  /* rearrange each row of el2e such that the local edges are ordered in local lexicographic order */
  INT i6, i4, ej1, ej2, *jtmp = (INT*)calloc(6, sizeof(INT));
  for (i = 0; i < ns; i++) {
    i6 = 6 * i; i4 = 4 * i;
    memcpy(jtmp, &el2e->JA[i6], 6 * sizeof(INT));

    for (j = 0; j < 6; j++) {
      ej1 = e2v->JA[2 * jtmp[j]];
      ej2 = e2v->JA[2 * jtmp[j] + 1];
      if (((sc->nodes[i4] == ej1) && (sc->nodes[i4 + 1] == ej2))
         || ((sc->nodes[i4] == ej2) && (sc->nodes[i4 + 1] == ej1)))
      {el2e->JA[i6] = jtmp[j]; continue;}

      if (((sc->nodes[i4] == ej1) && (sc->nodes[i4 + 2] == ej2))
         || ((sc->nodes[i4] == ej2) && (sc->nodes[i4 + 2] == ej1)))
      {el2e->JA[i6 + 1] = jtmp[j]; continue;}

      if (((sc->nodes[i4] == ej1) && (sc->nodes[i4 + 3] == ej2))
         || ((sc->nodes[i4] == ej2) && (sc->nodes[i4 + 3] == ej1)))
      {el2e->JA[i6 + 2] = jtmp[j]; continue;}

      if (((sc->nodes[i4 + 1] == ej1) && (sc->nodes[i4 + 2] == ej2))
         || ((sc->nodes[i4 + 1] == ej2) && (sc->nodes[i4 + 2] == ej1)))
      {el2e->JA[i6 + 3] = jtmp[j]; continue;}

      if (((sc->nodes[i4 + 1] == ej1) && (sc->nodes[i4 + 3] == ej2))
         || ((sc->nodes[i4 + 1] == ej2) && (sc->nodes[i4 + 3] == ej1)))
      {el2e->JA[i6 + 4] = jtmp[j]; continue;}

      if (((sc->nodes[i4 + 2] == ej1) && (sc->nodes[i4 + 3] == ej2))
         || ((sc->nodes[i4 + 2] == ej2) && (sc->nodes[i4 + 3] == ej1)))
      {el2e->JA[i6 + 5] = jtmp[j]; continue;}
    }
  }
  icsr_free(&el2e0);
  free(jtmp);
  return;
}

/*!
 * \fn void uniformrefine2d(scomplex *sc)
 *
 * \brief Returns a uniform refinement of the input 2d grid sc and update its contents,
 *        in particular the element-to-vertex matrix and coordinates of grid vertices.
 *        This function simply addes new vertices at midpoints of each edge and divides
 *        each triangle into four subtriangles.
 *
 * \param sc  Pointer to scomplex grid structure.
 *
 *  \note Yuwen 20210606.
 */
void uniformrefine2d(scomplex* sc) {
  if (sc->dim != 2) {
    fprintf(stderr, "\n%% *** ERROR: %s called but the spatial dimension is not 2   !\n",
	    __FUNCTION__); exit(2);
  }
  iCSRmat e2v, el2e;
  get_edge2d(&e2v, &el2e, sc);
  INT nv = sc->nv, nv2 = nv * 2, ns = sc->ns, ne = e2v.row, i, i12, i3, nnz_p;
  INT* el2v = (INT*)calloc(3 * ns, sizeof(INT));
  memcpy(el2v, sc->nodes, 3 * ns * sizeof(INT));

  sc->ns = 4 * ns; sc->nv = nv + ne;
  sc->x = (REAL*)realloc(sc->x, 2 * (nv + ne) * sizeof(REAL));
  sc->nodes = (INT*)realloc(sc->nodes, 3 * 4 * ns * sizeof(INT));
  for (i = 0; i < ne; i++) {
    sc->x[nv2 + 2 * i] = (sc->x[2 * e2v.JA[2 * i]] + sc->x[2 * e2v.JA[2 * i + 1]]) * 0.5e0;
    sc->x[nv2 + 2 * i + 1] = (sc->x[2 * e2v.JA[2 * i] + 1] + sc->x[2 * e2v.JA[2 * i + 1] + 1]) * 0.5e0;
  }
  //
  sc->parent_v->row = sc->nv;// new num vertices
  sc->parent_v->col = nv;// old num vertices;
  nnz_p = sc->parent_v->IA[nv]; // where to start adding:
  sc->parent_v->nnz += 2 * ne;
  sc->parent_v->IA = (INT*)realloc(sc->parent_v->IA, (sc->nv + 1) * sizeof(INT));
  sc->parent_v->JA = (INT*)realloc(sc->parent_v->JA, (sc->parent_v->nnz) * sizeof(INT));
  sc->parent_v->val = (INT*)realloc(sc->parent_v->val, (sc->parent_v->nnz) * sizeof(INT));
  //
  for (i = 0; i < ne; i++) {
    // one added vertex between e2v.JA[2*i] and e2v.JA[2*i+1];
    sc->parent_v->JA[nnz_p] = e2v.JA[2 * i];
    sc->parent_v->val[nnz_p] = 0;
    nnz_p++;
    sc->parent_v->JA[nnz_p] = e2v.JA[2 * i + 1];
    sc->parent_v->val[nnz_p] = 0;
    nnz_p++;
    sc->parent_v->IA[nv + i + 1] = nnz_p;
  }
  //  fprintf(stdout,"\n%%%%nnz_ia=%d;nnz_p=%d\n\n",sc->parent_v->IA[sc->nv],nnz_p);fflush(stdout);
  //
  for (i = 0; i < ns; i++) {
    i12 = 12 * i; i3 = 3 * i;

    sc->nodes[i12] = el2v[i3];
    sc->nodes[i12 + 1] = nv + el2e.JA[i3]; sc->nodes[i12 + 2] = nv + el2e.JA[i3 + 1];

    sc->nodes[i12 + 3] = el2v[i3 + 1];
    sc->nodes[i12 + 4] = nv + el2e.JA[i3]; sc->nodes[i12 + 5] = nv + el2e.JA[i3 + 2];

    sc->nodes[i12 + 6] = el2v[i3 + 2];
    sc->nodes[i12 + 7] = nv + el2e.JA[i3 + 1]; sc->nodes[i12 + 8] = nv + el2e.JA[i3 + 2];

    sc->nodes[i12 + 9] = nv + el2e.JA[i3];
    sc->nodes[i12 + 10] = nv + el2e.JA[i3 + 1]; sc->nodes[i12 + 11] = nv + el2e.JA[i3 + 2];
  }
  icsr_free(&e2v);
  icsr_free(&el2e);
  free(el2v);
  // reallocate arrays:
  haz_scomplex_realloc(sc);
  // find neighbors
  //  find_nbr(sc->ns,sc->nv,sc->dim,sc->nodes,sc->nbr);
  return;
}

/***********************************************************************************************/
/*!
 * \fn void uniformrefine3d(scomplex *sc)
 *
 * \brief Returns a uniform refinement of the input 3d grid sc and update its contents,
 *        in particular the element-to-vertex matrix and coordinates of grid vertices.
 *        This function is based on the paper
 *        "J. Bey: Simplicial grid refinement: on Freudenthal's algorithm
 *         and the optimal number of congruence classes". It addes new vertices at midpoints
 *        of all edges and divides each tetrahedron into eight smaller tetrahedra of equal volumes.
 *
 * \param sc  Pointer to scomplex grid structure.
 *
 * \note Yuwen 20210606.
 */
void uniformrefine3d(scomplex* sc) {
  if (sc->dim != 3) {
    fprintf(stderr, "\n%% *** ERROR: %s called but the spatial dimension = %lld is not 3  !\n",
	    __FUNCTION__, (long long)sc->dim); exit(3);
  }
  iCSRmat e2v, el2e;
  get_edge3d(&e2v, &el2e, sc);
  INT nv = sc->nv, nv3 = nv * 3, ns = sc->ns, ne = e2v.row, i, i32, i2, i3, i4, i6, nnz_p;
  // INT tmp,j;
  sc->ns = 8 * ns; sc->nv = nv + ne;
  sc->vols = (REAL*)realloc(sc->vols, sc->ns * sizeof(REAL));
  sc->x = (REAL*)realloc(sc->x, 3 * sc->nv * sizeof(REAL));
  for (i = 0; i < ne; i++) {
    i3 = 3 * i; i2 = 2 * i;
    sc->x[nv3 + i3] = (sc->x[3 * e2v.JA[i2]] + sc->x[3 * e2v.JA[i2 + 1]]) * 0.5e0;
    sc->x[nv3 + i3 + 1] = (sc->x[3 * e2v.JA[i2] + 1] + sc->x[3 * e2v.JA[i2 + 1] + 1]) * 0.5e0;
    sc->x[nv3 + i3 + 2] = (sc->x[3 * e2v.JA[i2] + 2] + sc->x[3 * e2v.JA[i2 + 1] + 2]) * 0.5e0;
  }
  //
  sc->parent_v->row = sc->nv;// new num vertices
  sc->parent_v->col = nv;// old num vertices;
  nnz_p = sc->parent_v->IA[nv]; // where to start adding:
  sc->parent_v->nnz += 2 * ne;
  sc->parent_v->IA = (INT*)realloc(sc->parent_v->IA, (sc->nv + 1) * sizeof(INT));
  sc->parent_v->JA = (INT*)realloc(sc->parent_v->JA, (sc->parent_v->nnz) * sizeof(INT));
  sc->parent_v->val = (INT*)realloc(sc->parent_v->val, (sc->parent_v->nnz) * sizeof(INT));
  //
  for (i = 0; i < ne; i++) {
    i2 = 2 * i;
    // one added vertex between e2v.JA[2*i] and e2v.JA[2*i+1];
    sc->parent_v->JA[nnz_p] = e2v.JA[i2];
    sc->parent_v->val[nnz_p] = 0;
    nnz_p++;
    sc->parent_v->JA[nnz_p] = e2v.JA[i2 + 1];
    sc->parent_v->val[nnz_p] = 0;
    nnz_p++;
    sc->parent_v->IA[nv + i + 1] = nnz_p;
  }
  //  fprintf(stdout,"\n%%%%nnz_ia=%d;nnz_p=%d\n\n",sc->parent_v->IA[sc->nv],nnz_p);fflush(stdout);
  //
  icsr_free(&e2v);
  INT* el2v = (INT*)calloc(4 * ns, sizeof(INT));

  memcpy(el2v, sc->nodes, 4 * ns * sizeof(INT));
  free(sc->nodes);
  sc->nodes = (INT*)calloc(4 * sc->ns, sizeof(INT));

  for (i = 0; i < ns; i++) {
    i32 = 32 * i; i6 = 6 * i; i4 = 4 * i;
    /* four sub-tetrahedra at four corners, local ordering is important*/
    sc->nodes[i32] = el2v[i4]; sc->nodes[i32 + 1] = nv + el2e.JA[i6];
    sc->nodes[i32 + 2] = nv + el2e.JA[i6 + 1]; sc->nodes[i32 + 3] = nv + el2e.JA[i6 + 2];

    sc->nodes[i32 + 4] = nv + el2e.JA[i6]; sc->nodes[i32 + 5] = el2v[i4 + 1];
    sc->nodes[i32 + 6] = nv + el2e.JA[i6 + 3]; sc->nodes[i32 + 7] = nv + el2e.JA[i6 + 4];

    sc->nodes[i32 + 8] = nv + el2e.JA[i6 + 1]; sc->nodes[i32 + 9] = nv + el2e.JA[i6 + 3];
    sc->nodes[i32 + 10] = el2v[i4 + 2]; sc->nodes[i32 + 11] = nv + el2e.JA[i6 + 5];

    sc->nodes[i32 + 12] = nv + el2e.JA[i6 + 2]; sc->nodes[i32 + 13] = nv + el2e.JA[i6 + 4];
    sc->nodes[i32 + 14] = nv + el2e.JA[i6 + 5]; sc->nodes[i32 + 15] = el2v[i4 + 3];
    /* four sub-tetrahedra from the inner octahedron, local ordering is important */
    sc->nodes[i32 + 16] = nv + el2e.JA[i6]; sc->nodes[i32 + 17] = nv + el2e.JA[i6 + 1];
    sc->nodes[i32 + 18] = nv + el2e.JA[i6 + 2]; sc->nodes[i32 + 19] = nv + el2e.JA[i6 + 4];

    sc->nodes[i32 + 20] = nv + el2e.JA[i6]; sc->nodes[i32 + 21] = nv + el2e.JA[i6 + 1];
    sc->nodes[i32 + 22] = nv + el2e.JA[i6 + 3]; sc->nodes[i32 + 23] = nv + el2e.JA[i6 + 4];

    sc->nodes[i32 + 24] = nv + el2e.JA[i6 + 1]; sc->nodes[i32 + 25] = nv + el2e.JA[i6 + 2];
    sc->nodes[i32 + 26] = nv + el2e.JA[i6 + 4]; sc->nodes[i32 + 27] = nv + el2e.JA[i6 + 5];

    sc->nodes[i32 + 28] = nv + el2e.JA[i6 + 1]; sc->nodes[i32 + 29] = nv + el2e.JA[i6 + 3];
    sc->nodes[i32 + 30] = nv + el2e.JA[i6 + 4]; sc->nodes[i32 + 31] = nv + el2e.JA[i6 + 5];
  }
  icsr_free(&el2e);
  free(el2v);
  //  reallocate unallocated arrays
  haz_scomplex_realloc(sc);
  //  find_nbr(sc->ns,sc->nv,sc->dim,sc->nodes,sc->nbr);
  return;
}
/***********************************************************************************************/
/*!
 * \fn void get_edge_nd(iCSRmat *e2v, iCSRmat *el2e, scomplex *sc)
 *
 * \brief Returns the edge-to-vertex matrix e2v and the element-to-edge
 *        matrix el2e for an n-dimensional simplicial complex.
 *        Generalizes get_edge2d and get_edge3d to arbitrary dimension.
 *
 *        Edges are ordered in global ascending lexicographic order.
 *        Each row of el2e lists edges in local lexicographic order:
 *        (0,1),(0,2),...,(0,n),(1,2),...,(n-1,n).
 *
 * \param e2v   Pointer to iCSRmat edge-to-vertex matrix (ne x nv)
 * \param el2e  Pointer to iCSRmat element-to-edge matrix (ns x ne)
 * \param sc    Pointer to scomplex
 *
 * \note Ludmil 20260322.
 */
void get_edge_nd(iCSRmat *e2v, iCSRmat *el2e, scomplex *sc)
{
  INT dim = sc->dim;
  INT n1 = dim + 1;                     /* vertices per simplex */
  INT ne_local = dim * n1 / 2;          /* edges per simplex = C(n1,2) */
  INT nv = sc->nv, ns = sc->ns;
  INT total_pairs = ne_local * ns;
  INT i, j, k, p, q, ih, counter;
  INT *el2v = sc->nodes;
  /* --- build (ii,jj) pairs: for each simplex, sort vertices ascending,
         enumerate all C(n1,2) pairs (vi,vj) with i<j --- */
  ivector ii = ivec_create(total_pairs);
  ivector jj = ivec_create(total_pairs);
  iCSRmat U = icsr_create(0, 0, 0);
  iCSRmat el2v_csr = icsr_create(ns, nv, n1 * ns);
  iCSRmat v2e, el2e0;
  el2v_csr.IA[0] = 0;
  INT *sv = (INT *)malloc(n1 * sizeof(INT)); /* sorted local vertices */
  for (k = 0; k < ns; k++) {
    /* copy and sort the n1 vertex indices of simplex k */
    memcpy(sv, el2v + n1 * k, n1 * sizeof(INT));
    /* insertion sort (n1 is small) */
    for (i = 1; i < n1; i++) {
      INT tmp = sv[i];
      j = i - 1;
      while (j >= 0 && sv[j] > tmp) { sv[j + 1] = sv[j]; j--; }
      sv[j + 1] = tmp;
    }
    /* enumerate C(n1,2) pairs in lexicographic order */
    counter = 0;
    for (p = 0; p < n1; p++) {
      for (q = p + 1; q < n1; q++) {
        ii.val[k + counter * ns] = sv[p];
        jj.val[k + counter * ns] = sv[q];
        counter++;
      }
    }
    el2v_csr.IA[k + 1] = el2v_csr.IA[k] + n1;
  }
  free(sv);
  memcpy(el2v_csr.JA, el2v, n1 * ns * sizeof(INT));
  for (i = 0; i < n1 * ns; i++) el2v_csr.val[i] = 1;
  icsr_uniqueij(&U, &ii, &jj);
  ivec_free(&ii); ivec_free(&jj);
  /* --- form e2v (ne x nv, 2 entries per row) --- */
  INT ne = U.nnz;
  e2v->row = ne; e2v->col = nv; e2v->nnz = 2 * ne;
  INT *ia = U.IA, *ja = U.JA;
  INT nv_U = U.row; /* may be < nv if highest-numbered vertices never appear first in a sorted pair */
  e2v->JA  = (INT *)malloc(2 * ne * sizeof(INT));
  e2v->val = (INT *)malloc(2 * ne * sizeof(INT));
  for (i = 0; i < 2 * ne; i++) e2v->val[i] = 1;
  counter = 0;
  for (i = 0; i < nv_U; i++) {
    ih = ia[i + 1] - ia[i];
    for (j = 0; j < ih; j++) {
      e2v->JA[counter + 2 * j]     = i;
      e2v->JA[counter + 2 * j + 1] = ja[ia[i] + j];
    }
    counter += 2 * ih;
  }
  e2v->IA = (INT *)malloc((ne + 1) * sizeof(INT));
  e2v->IA[0] = 0;
  for (i = 0; i < ne; i++) e2v->IA[i + 1] = e2v->IA[i] + 2;
  icsr_free(&U);
  /* --- el2e via sparse multiply: el2v * v2e, keep entries == 2 --- */
  icsr_trans(e2v, &v2e);
  icsr_mxm(&el2v_csr, &v2e, &el2e0);
  icsr_free(&v2e);
  icsr_free(&el2v_csr);
  ia = el2e0.IA; ja = el2e0.JA;
  INT *val = el2e0.val;
  el2e->row = ns; el2e->col = ne; el2e->nnz = ne_local * ns;
  el2e->IA  = (INT *)malloc((ns + 1) * sizeof(INT));
  el2e->JA  = (INT *)malloc(ne_local * ns * sizeof(INT));
  el2e->val = NULL;
  el2e->IA[0] = 0;
  counter = 0;
  for (i = 0; i < ns; i++) {
    el2e->IA[i + 1] = el2e->IA[i] + ne_local;
    ih = ia[i + 1] - ia[i];
    for (j = 0; j < ih; j++) {
      if (val[ia[i] + j] == 2) {
        el2e->JA[counter] = ja[ia[i] + j];
        counter++;
      }
    }
  }
  /* --- reorder each row of el2e to local lexicographic order --- */
  INT *jtmp = (INT *)malloc(ne_local * sizeof(INT));
  for (i = 0; i < ns; i++) {
    INT row_off = ne_local * i;
    INT el_off  = n1 * i;
    memcpy(jtmp, el2e->JA + row_off, ne_local * sizeof(INT));
    for (j = 0; j < ne_local; j++) {
      INT ej1 = e2v->JA[2 * jtmp[j]];
      INT ej2 = e2v->JA[2 * jtmp[j] + 1];
      /* find which local pair (p,q) this edge matches */
      INT loc = 0;
      for (p = 0; p < n1; p++) {
        for (q = p + 1; q < n1; q++) {
          INT vp = el2v[el_off + p], vq = el2v[el_off + q];
          if ((vp == ej1 && vq == ej2) || (vp == ej2 && vq == ej1)) {
            el2e->JA[row_off + loc] = jtmp[j];
            goto next_edge;
          }
          loc++;
        }
      }
      next_edge: ;
    }
  }
  free(jtmp);
  icsr_free(&el2e0);
}
/***********************************************************************************************/
/*!
 * \fn void uniformrefine(scomplex *sc)
 *
 * \brief Uniform (Freudenthal) refinement in arbitrary dimension.
 *        Each n-simplex is subdivided into 2^n children by inserting
 *        midpoints on all edges. Uses binary-vector enumeration of
 *        children consistent with Bey's algorithm for n=3.
 *
 * \param sc  Pointer to scomplex grid structure.
 *
 * \note Ludmil 20260322.
 */
void uniformrefine(scomplex *sc)
{
  INT dim  = sc->dim;
  INT n1   = dim + 1;
  INT ne_local = dim * n1 / 2;
  INT nv   = sc->nv;
  INT ns   = sc->ns;
  INT nbig = sc->nbig;
  iCSRmat e2v, el2e;
  get_edge_nd(&e2v, &el2e, sc);
  INT ne = e2v.row;
  INT nchildren = (1 << dim);  /* 2^dim */
  INT ns_new = nchildren * ns;
  INT nv_new = nv + ne;
  INT i, d, nnz_p;
  /* --- compute midpoint coordinates --- */
  sc->x = (REAL *)realloc(sc->x, nbig * nv_new * sizeof(REAL));
  for (i = 0; i < ne; i++) {
    INT va = e2v.JA[2 * i], vb = e2v.JA[2 * i + 1];
    for (d = 0; d < nbig; d++) {
      sc->x[nbig * (nv + i) + d] =
        0.5 * (sc->x[nbig * va + d] + sc->x[nbig * vb + d]);
    }
  }
  /* --- update parent_v --- */
  sc->parent_v->row = nv_new;
  sc->parent_v->col = nv;
  nnz_p = sc->parent_v->IA[nv];
  sc->parent_v->nnz += 2 * ne;
  sc->parent_v->IA  = (INT *)realloc(sc->parent_v->IA,  (nv_new + 1) * sizeof(INT));
  sc->parent_v->JA  = (INT *)realloc(sc->parent_v->JA,  sc->parent_v->nnz * sizeof(INT));
  sc->parent_v->val = (INT *)realloc(sc->parent_v->val, sc->parent_v->nnz * sizeof(INT));
  for (i = 0; i < ne; i++) {
    sc->parent_v->JA[nnz_p]  = e2v.JA[2 * i];
    sc->parent_v->val[nnz_p] = 0;
    nnz_p++;
    sc->parent_v->JA[nnz_p]  = e2v.JA[2 * i + 1];
    sc->parent_v->val[nnz_p] = 0;
    nnz_p++;
    sc->parent_v->IA[nv + i + 1] = nnz_p;
  }
  icsr_free(&e2v);
  /* --- save old connectivity and allocate new --- */
  INT *el2v_old = (INT *)malloc(n1 * ns * sizeof(INT));
  memcpy(el2v_old, sc->nodes, n1 * ns * sizeof(INT));
  free(sc->nodes);
  sc->nodes = (INT *)calloc(n1 * ns_new, sizeof(INT));
  sc->ns = ns_new;
  sc->nv = nv_new;
  /* --- precompute local edge-index lookup: rank(p,q) for 0<=p<q<=dim --- */
  /* rank(p,q) = p*(2*dim - p - 1)/2 + (q - p - 1) */
  /* --- fill child connectivity using binary-vector formula --- */
  for (i = 0; i < ns; i++) {
    INT el_off  = ne_local * i;   /* offset into el2e.JA */
    INT *old_verts = el2v_old + n1 * i;
    INT *edges = el2e.JA + el_off;
    for (INT s = 0; s < nchildren; s++) {
      INT child_off = n1 * (nchildren * i + s);
      /* popcount of s */
      INT pc = 0;
      { INT tmp = s; while (tmp) { pc += tmp & 1; tmp >>= 1; } }
      /* compute L[0], U[0] */
      INT L = 0, U = pc;
      /* vertex w_0 */
      if (L == U)
        sc->nodes[child_off] = old_verts[L];
      else {
        INT rank_lu = L * (2 * dim - L + 1) / 2 + (U - L - 1);
        sc->nodes[child_off] = nv + edges[rank_lu];
      }
      /* vertices w_1 .. w_dim */
      for (INT kk = 1; kk <= dim; kk++) {
        INT bit_k = (s >> (dim - kk)) & 1;
        if (bit_k == 0) {
          /* L stays, U increments */
          U = U + 1;
        } else {
          /* L increments, U stays */
          L = L + 1;
        }
        if (L == U)
          sc->nodes[child_off + kk] = old_verts[L];
        else {
          INT rank_lu = L * (2 * dim - L + 1) / 2 + (U - L - 1);
          sc->nodes[child_off + kk] = nv + edges[rank_lu];
        }
      }
    }
  }
  icsr_free(&el2e);
  free(el2v_old);
  /* reallocate auxiliary arrays */
  haz_scomplex_realloc(sc);
}
/***********************************************************************************************/
/*!
 * \fn void uniformrefine_marked(scomplex *sc, ivector *marked)
 *
 * \brief Refinement dispatcher honoring the HAZmath AMR design rule
 *        (20260530): Bey (Freudenthal) refinement is used ONLY for
 *        UNIFORM refinement (all simplices marked); LOCAL/selective
 *        refinement is delegated to DGS newest-vertex bisection
 *        (refine()), which has a proven conforming completion.
 *
 *        The former hand-written "selective Bey + face-Bey/bisection
 *        closure" was REMOVED: it produced non-conforming meshes. This is
 *        the classic red/green incompatibility -- a pure-bisection closure
 *        cannot reproduce a Bey-refined face (whose central sub-simplex
 *        has all-midpoint vertices), and the face-Bey closure left cracks
 *        (facet-conforming but non-manifold boundary). Doing local
 *        refinement with DGS avoids this class of bug entirely.
 *
 * \param sc      Pointer to scomplex grid structure.
 * \param marked  ivector of length sc->ns; nonzero entries = marked.
 *                If NULL, refine all (uniform Bey).
 *
 * \note Ludmil 20260322; redesigned 20260530.
 *
 * \warning When marking is PARTIAL this calls refine(1,sc,marked), which
 *          builds a DGS refinement hierarchy in sc (sc->ns then includes
 *          parents). Callers needing a flat leaf mesh must follow the DGS
 *          usage pattern (scfinest()/scfinalize()), as the example
 *          amr_grids.c and the test programs do.
 */
void uniformrefine_marked(scomplex *sc, ivector *marked)
{
  if (!marked) { uniformrefine(sc); return; }
  INT i, nmarked = 0;
  for (i = 0; i < sc->ns; i++) if (marked->val[i]) nmarked++;
  if (nmarked == 0) return;                              /* nothing marked */
  if (nmarked == sc->ns) { uniformrefine(sc); return; }  /* all -> uniform Bey */
  /* partial -> local refinement via DGS bisection (conforming completion) */
  refine((INT)1, sc, marked);
}
