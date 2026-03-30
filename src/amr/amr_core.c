/*! \file src/amr/amr_core.c
 *
 *  Authors: James Adler, Xiaozhe Hu, and Ludmil Zikatanov
 *           HAZmath (https://hazmath.net)
 *           Created with the help of Claude (Anthropic)
 *
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note Core refinement routines: DGS bisection (Diening-Gehring-Storn
 *        2025), vertex coloring, simplex addition, and make_uniform_mesh.
 *        Extracted from scomplex.c 20260322.
 */
#include "hazmath.h"
/**********************************************************************/
/*!
 * \fn INT haz_add_simplex(INT is, scomplex *sc,REAL *xnew, INT *pv,
 *                         INT ibnew,INT csysnew,INT nsnew, INT nvnew)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
INT haz_add_simplex(INT is, scomplex* sc, REAL* xnew, INT* pv, INT ibnew,
                    INT csysnew, INT nsnew, INT nvnew) {
  /* adds nodes and coords as well */
  INT dim = sc->dim, nbig = sc->nbig, n1 = dim + 1, nv = sc->nv;  // ns=sc->ns;
  INT ks0 = sc->child0[is], ksn = sc->childn[is];
  INT isc0 = ks0 * n1, iscn = ksn * n1;
  //  INT *dsti,*srci;
  INT j, j0, jn, nnz_pv;
  REAL* dstr;
  /* nodes  AND neighbors */
  sc->nbr = realloc(sc->nbr, (nsnew * n1) * sizeof(INT));
  sc->nodes = realloc(sc->nodes, (nsnew * n1) * sizeof(INT));
  for (j = 0; j < n1; j++) {
    j0 = isc0 + j;
    jn = iscn + j;
    //
    sc->nodes[j0] = -1;
    sc->nbr[j0] = -1;
    sc->nodes[jn] = -1;
    sc->nbr[jn] = -1;
  }
  // new vertex (if any!!!)
  if (nvnew != nv) {
    sc->x = realloc(sc->x, (nvnew * nbig) * sizeof(REAL));
    dstr = (sc->x + nv * nbig);
    memcpy(dstr, xnew, dim * sizeof(REAL));
    if (xnew) free(xnew);
    sc->bndry = realloc(sc->bndry, (nvnew) * sizeof(INT));
    sc->bndry[nv] = ibnew;
    sc->csys = realloc(sc->csys, (nvnew) * sizeof(INT));
    sc->csys[nv] = csysnew;
    nnz_pv = sc->parent_v->nnz;
    sc->parent_v->row = nvnew;
    sc->parent_v->col = nv;
    sc->parent_v->JA = realloc(sc->parent_v->JA, (nnz_pv + 2) * sizeof(REAL));
    sc->parent_v->JA[nnz_pv] = pv[0];
    sc->parent_v->JA[nnz_pv + 1] = pv[1];
    sc->parent_v->val = realloc(sc->parent_v->val, (nnz_pv + 2) * sizeof(REAL));
    sc->parent_v->val[nnz_pv] = sc->level + 1;
    sc->parent_v->val[nnz_pv + 1] = sc->level + 1;
    sc->parent_v->nnz += 2;
    sc->parent_v->IA = realloc(sc->parent_v->IA, (nvnew + 1) * sizeof(REAL));
    sc->parent_v->IA[nvnew] = sc->parent_v->nnz;
    /* fprintf(stdout,"\nnv=%d;
     * nvnew=%d;nnz_pv=%d(pv[0]=%d,pv[1]=%d)",nv,nvnew,sc->parent_v->nnz,pv[0],pv[1]);
     */
  }
  // generation
  sc->gen = realloc(sc->gen, (nsnew) * sizeof(INT));
  sc->gen[ks0] = sc->gen[is] + 1;
  sc->gen[ksn] = sc->gen[is] + 1;
  // marked
  sc->marked = realloc(sc->marked, (nsnew) * sizeof(INT));
  sc->marked[ks0] = sc->marked[is];
  sc->marked[ksn] = sc->marked[is];
  // flags
  sc->flags = realloc(sc->flags, (nsnew) * sizeof(INT));
  sc->flags[ks0] = sc->flags[is];
  sc->flags[ksn] = sc->flags[is];
  // parents
  sc->parent = realloc(sc->parent, (nsnew) * sizeof(INT));
  sc->parent[ks0] = is;
  sc->parent[ksn] = is;
  // child0
  sc->child0 = realloc(sc->child0, (nsnew) * sizeof(INT));
  sc->child0[ks0] = -1;
  sc->child0[ksn] = -1;
  // childn
  sc->childn = realloc(sc->childn, (nsnew) * sizeof(INT));
  sc->childn[ks0] = -1;
  sc->childn[ksn] = -1;
  // volumes
  sc->vols = realloc(
      sc->vols, (nsnew) * sizeof(REAL));  // we can calculate all volumes at the
                                          // end, so just allocate here.
  // scalars
  sc->ns = nsnew;
  sc->nv = nvnew;
  return 0;
}
/**********************************************************************/
/*!
 * \fn static INT set_color(scomplex *sc, INT *color)
 *
 * \brief Generalized coloring with N+1 colors (Algorithm 2 from
 *        Diening, Gehring, Storn, "Adaptive Mesh Refinement for
 *        Arbitrary Initial Triangulations", Found. Comput. Math.,
 *        2025, DOI: 10.1007/s10208-024-09642-1).
 *
 *        Assigns to each vertex the smallest color not already
 *        attained by a neighboring vertex. The number of colors N+1
 *        is bounded by the maximal vertex degree + 1.
 *
 * \param sc     I: simplicial complex (uses nodes, ns, nv, n)
 * \param color  O: array of size nv; color[v] in {0,...,N}
 *
 * \return N (the largest color assigned)
 */
static INT set_color(scomplex* sc, INT* color) {
  INT nv = sc->nv, ns = sc->ns, dim = sc->dim;
  INT n1 = dim + 1, i, j, k, v, w, c;
  /* Build vertex-to-vertex adjacency from the element connectivity.
     We use a simple approach: for each simplex, all pairs of its
     vertices are neighbors. We store adjacency in CSR format. */
  /* Step 1: count edges per vertex (upper bound via element connectivity) */
  INT* deg = (INT*)calloc(nv, sizeof(INT));
  for (i = 0; i < ns; i++) {
    INT* el = sc->nodes + i * n1;
    for (j = 0; j < n1; j++) {
      deg[el[j]] += dim; /* at most dim neighbors per element */
    }
  }
  /* Allocate adjacency lists (with room for duplicates; we handle them) */
  INT* adjptr = (INT*)calloc(nv + 1, sizeof(INT));
  adjptr[0] = 0;
  for (i = 0; i < nv; i++) adjptr[i + 1] = adjptr[i] + deg[i];
  INT nnz_adj = adjptr[nv];
  INT* adjind = (INT*)malloc(nnz_adj * sizeof(INT));
  INT* pos = (INT*)calloc(nv, sizeof(INT)); /* current insert position */
  for (i = 0; i < ns; i++) {
    INT* el = sc->nodes + i * n1;
    for (j = 0; j < n1; j++) {
      v = el[j];
      for (k = 0; k < n1; k++) {
        if (k == j) continue;
        w = el[k];
        adjind[adjptr[v] + pos[v]] = w;
        pos[v]++;
      }
    }
  }
  free(pos);
  free(deg);
  /* Step 2: greedy coloring */
  INT max_color = 0;
  /* max possible colors is bounded by max vertex degree + 1;
     for safety allocate nv+1 entries */
  INT* used =
      (INT*)calloc(nv + 1, sizeof(INT)); /* flag array: used[c]=v+1 means color
                                            c is used by a neighbor of v */
  for (v = 0; v < nv; v++)
    color[v] = -1; /* uncolored = infinity in the paper */
  for (v = 0; v < nv; v++) {
    /* Mark colors used by neighbors of v */
    for (j = adjptr[v]; j < adjptr[v + 1]; j++) {
      w = adjind[j];
      if (color[w] >= 0) used[color[w]] = v + 1; /* mark as used for vertex v */
    }
    /* Find smallest color not used by any neighbor */
    for (c = 0;; c++) {
      if (used[c] != v + 1) {
        color[v] = c;
        break;
      }
    }
    if (c > max_color) max_color = c;
  }
  free(used);
  free(adjind);
  free(adjptr);
  return max_color; /* N */
}
/**********************************************************************/
/*!
 * \fn static void dgs_initialize(scomplex *sc, INT *color, INT N)
 *
 * \brief Maubach initialization using a generalized coloring
 *        (Definition 2 from Diening, Gehring, Storn, "Adaptive Mesh
 *        Refinement for Arbitrary Initial Triangulations",
 *        Found. Comput. Math., 2025,
 *        DOI: 10.1007/s10208-024-09642-1).
 *
 *        For each initial simplex T with vertices v0,...,vn, we sort
 *        the vertices so that c(v_j) < c(v_{j+1}) for all j.  Then
 *        the tagged simplex is:
 *          T = [vn, v0, v1,...,v_{n-1}]_n  if c(vn) == N
 *          T = [v0, v1,...,vn]_n            otherwise
 *        The tag gamma = n for all initial simplices.
 *
 * \param sc     I/O: simplicial complex (nodes are reordered in place)
 * \param color  I: vertex coloring array of size nv
 * \param N      I: the largest color
 */
static void dgs_initialize(scomplex* sc, INT* color, INT N) {
  INT ns = sc->ns, nv = sc->nv, dim = sc->dim, n1 = dim + 1;
  INT i, j;
  /* Allocate gen_N and set initial values: gen_N(v) = -color(v) */
  sc->gen_N = (INT*)calloc(nv, sizeof(INT));
  sc->dgs_N = N;
  for (i = 0; i < nv; i++) sc->gen_N[i] = -color[i];
  /* Sort each simplex by DECREASING gen_N */
  for (i = 0; i < ns; i++) {
    INT* el = sc->nodes + i * n1;
    for (j = 0; j < n1; j++)
      for (INT k = j + 1; k < n1; k++)
        if (sc->gen_N[el[j]] < sc->gen_N[el[k]]) {
          INT t = el[j]; el[j] = el[k]; el[k] = t;
        }
    sc->gen[i] = 0;
  }
  /* nbr is not reordered — it will be set to -1 for children
     and rebuilt by find_nbr in scfinalize */
}
/**********************************************************************/
/*!
 * \brief Build v2s index for current LEAF simplices only.
 *        Children created during refinement are NOT added.
 *        This prevents infinite recursion in the DGS closure.
 */
static void v2s_build(scomplex* sc) {
  INT ns = sc->ns, nv = sc->nv, n1 = sc->dim + 1;
  INT i, j, v;
  /* Count leaf simplex-vertex pairs */
  INT total = 0;
  for (i = 0; i < ns; i++)
    if (sc->child0[i] < 0) total += n1;
  sc->v2s_alloc = total > 0 ? total : 1;
  sc->v2s_head = (INT*)malloc(nv * sizeof(INT));
  sc->v2s_next = (INT*)malloc(sc->v2s_alloc * sizeof(INT));
  sc->v2s_simp = (INT*)malloc(sc->v2s_alloc * sizeof(INT));
  for (i = 0; i < nv; i++) sc->v2s_head[i] = -1;
  sc->v2s_count = 0;
  sc->v2s_nv = nv;
  for (i = 0; i < ns; i++) {
    if (sc->child0[i] >= 0) continue; /* skip non-leaf */
    INT* el = sc->nodes + i * n1;
    for (j = 0; j < n1; j++) {
      v = el[j];
      INT k = sc->v2s_count++;
      sc->v2s_simp[k] = i;
      sc->v2s_next[k] = sc->v2s_head[v];
      sc->v2s_head[v] = k;
    }
  }
}
/**********************************************************************/
static void v2s_free_data(scomplex* sc) {
  if (sc->v2s_head) { free(sc->v2s_head); sc->v2s_head = NULL; }
  if (sc->v2s_next) { free(sc->v2s_next); sc->v2s_next = NULL; }
  if (sc->v2s_simp) { free(sc->v2s_simp); sc->v2s_simp = NULL; }
  sc->v2s_count = 0; sc->v2s_alloc = 0; sc->v2s_nv = 0;
}
/**********************************************************************/
/* Algorithm 4 helpers: lvl_N and type_N from gen_N (eq 8 in paper)  */
/* gen_N = N*(lvl_N - 1) + type_N,  type_N in {1,...,N}              */
/**********************************************************************/
static inline INT dgs_type_N(INT gen, INT N) {
  return ((gen - 1) % N + N) % N + 1;
}
static inline INT dgs_lvl_N(INT gen, INT N) {
  INT t = dgs_type_N(gen, N);
  return (gen - t) / N + 1;
}
/**********************************************************************/
/*!
 * \brief Get the bisection edge of simplex is using Algorithm 4.
 *        Simplex vertices must be sorted by DECREASING gen_N.
 *        Returns the two vertex indices (v_a, v_b) of the bse.
 */
static void dgs_get_bse(scomplex* sc, INT is, INT *va_out, INT *vb_out) {
  INT dim = sc->dim, n1 = dim + 1, N = sc->dgs_N;
  INT *el = sc->nodes + is * n1;
  INT vm = el[dim], vm1 = el[dim - 1];
  INT lvl_m = dgs_lvl_N(sc->gen_N[vm], N);
  INT lvl_m1 = dgs_lvl_N(sc->gen_N[vm1], N);
  if (lvl_m != lvl_m1) {
    *va_out = vm1; *vb_out = vm;
  } else {
    INT j = dim;
    for (INT k = 0; k < n1; k++)
      if (dgs_lvl_N(sc->gen_N[el[k]], N) == lvl_m) { j = k; break; }
    *va_out = el[j]; *vb_out = vm;
  }
}
/**********************************************************************/
/*!
 * \brief Bisect simplex is using Algorithm 4, creating a NEW midpoint.
 *        Children are sorted by decreasing gen_N.
 * \return index of the new midpoint vertex
 */
static INT haz_bisect_new(scomplex* sc, INT is) {
  INT dim = sc->dim, nbig = sc->nbig, n1 = dim + 1, N = sc->dgs_N;
  INT *el = sc->nodes + is * n1;
  /* Algorithm 4: determine bse and gen_N of midpoint */
  INT vm = el[dim], vm1 = el[dim - 1];
  INT lvl_m = dgs_lvl_N(sc->gen_N[vm], N);
  INT lvl_m1 = dgs_lvl_N(sc->gen_N[vm1], N);
  INT va_idx, gen_new; /* va_idx = index into el[] */
  if (lvl_m != lvl_m1) {
    va_idx = dim - 1;
    gen_new = sc->gen_N[vm1] + N;
  } else {
    va_idx = dim;
    for (INT k = 0; k < n1; k++)
      if (dgs_lvl_N(sc->gen_N[el[k]], N) == lvl_m) { va_idx = k; break; }
    gen_new = sc->gen_N[vm] + 2*N + 1 - dgs_type_N(sc->gen_N[el[va_idx]], N);
  }
  INT v_a = el[va_idx], v_b = el[dim];
  /* Create midpoint vertex via haz_add_simplex */
  INT nsnew = sc->ns + 2, nvnew = sc->nv + 1, nodnew = sc->nv;
  REAL* xnew = (REAL*)calloc(nbig, sizeof(REAL));
  for (INT j = 0; j < nbig; j++)
    xnew[j] = 0.5 * (sc->x[v_a * nbig + j] + sc->x[v_b * nbig + j]);
  /* Project onto curve if both endpoints share a polar coordinate system */
  INT csysnew = (sc->csys[v_a] < sc->csys[v_b]) ? sc->csys[v_a] : sc->csys[v_b];
  if (sc->csys[v_a] == sc->csys[v_b] && sc->systypes &&
      sc->csys[v_a] < sc->ncsys && sc->systypes[sc->csys[v_a]] == 1) {
    /* Both endpoints in the same polar coordinate system.
       Compute midpoint on the arc: use Cartesian midpoint's direction
       but at the average radius from the origin. */
    INT sys = sc->csys[v_a];
    REAL *ox = sc->csys_ox + sys * dim;
    REAL ra = 0, rb = 0, rm_dir = 0;
    for (INT j = 0; j < dim; j++) {
      REAL da = sc->x[v_a * nbig + j] - ox[j];
      REAL db = sc->x[v_b * nbig + j] - ox[j];
      ra += da * da;
      rb += db * db;
    }
    ra = sqrt(ra); rb = sqrt(rb);
    REAL rho_avg = 0.5 * (ra + rb);
    /* Direction from origin to Cartesian midpoint */
    for (INT j = 0; j < dim; j++) {
      REAL d = xnew[j] - ox[j];
      rm_dir += d * d;
    }
    rm_dir = sqrt(rm_dir);
    if (rm_dir > 1e-14) {
      REAL scale = rho_avg / rm_dir;
      for (INT j = 0; j < dim; j++)
        xnew[j] = ox[j] + scale * (xnew[j] - ox[j]);
    }
  }
  INT pv[2] = {v_a, v_b};
  INT ks0 = sc->ns, ksn = sc->ns + 1;
  sc->child0[is] = ks0;
  sc->childn[is] = ksn;
  haz_add_simplex(is, sc, xnew, pv, 0, csysnew, nsnew, nvnew);
  /* Set gen_N for the new vertex (gen_N array was grown by haz_add_simplex
     via realloc of bndry/csys but gen_N needs separate growth) */
  sc->gen_N = realloc(sc->gen_N, nvnew * sizeof(INT));
  sc->gen_N[nodnew] = gen_new;
  /* Build children: replace va with v', replace vb with v' */
  el = sc->nodes + is * n1;
  INT *c0 = sc->nodes + ks0 * n1;
  INT *cn = sc->nodes + ksn * n1;
  for (INT j = 0; j < n1; j++) c0[j] = el[j];
  c0[va_idx] = nodnew;
  for (INT j = 0; j < n1; j++) cn[j] = el[j];
  cn[dim] = nodnew;
  /* Sort children by decreasing gen_N */
  for (INT i = 0; i < n1; i++)
    for (INT j = i+1; j < n1; j++) {
      if (sc->gen_N[c0[i]] < sc->gen_N[c0[j]]) { INT t=c0[i]; c0[i]=c0[j]; c0[j]=t; }
      if (sc->gen_N[cn[i]] < sc->gen_N[cn[j]]) { INT t=cn[i]; cn[i]=cn[j]; cn[j]=t; }
    }
  for (INT j = 0; j < n1; j++) {
    sc->nbr[ks0 * n1 + j] = -1;
    sc->nbr[ksn * n1 + j] = -1;
  }
  return nodnew;
}
/**********************************************************************/
/*!
 * \brief Bisect simplex is using Algorithm 4, REUSING existing midpoint.
 */
static void haz_bisect_reuse(scomplex* sc, INT is, INT nodnew,
                             INT edge_va, INT edge_vb) {
  INT dim = sc->dim, n1 = dim + 1;
  INT *el = sc->nodes + is * n1;
  /* Find the local indices of the shared edge endpoints in this simplex.
     We bisect along the shared edge {edge_va, edge_vb}, NOT along the
     simplex's own bse — the bse may differ for children created during
     the conforming closure of a different edge. */
  INT ia = -1, ib = -1;
  for (INT k = 0; k < n1; k++) {
    if (el[k] == edge_va) ia = k;
    if (el[k] == edge_vb) ib = k;
  }
  INT v_a = edge_va, v_b = edge_vb;
  INT nsnew = sc->ns + 2, nvnew = sc->nv;
  INT pv[2] = {v_a, v_b};
  INT ks0 = sc->ns, ksn = sc->ns + 1;
  sc->child0[is] = ks0;
  sc->childn[is] = ksn;
  haz_add_simplex(is, sc, NULL, pv, 0, 0, nsnew, nvnew);
  el = sc->nodes + is * n1;
  INT *c0 = sc->nodes + ks0 * n1;
  INT *cn = sc->nodes + ksn * n1;
  for (INT j = 0; j < n1; j++) c0[j] = el[j];
  c0[ia] = nodnew;
  for (INT j = 0; j < n1; j++) cn[j] = el[j];
  cn[ib] = nodnew;
  for (INT i = 0; i < n1; i++)
    for (INT j = i+1; j < n1; j++) {
      if (sc->gen_N[c0[i]] < sc->gen_N[c0[j]]) { INT t=c0[i]; c0[i]=c0[j]; c0[j]=t; }
      if (sc->gen_N[cn[i]] < sc->gen_N[cn[j]]) { INT t=cn[i]; cn[i]=cn[j]; cn[j]=t; }
    }
  for (INT j = 0; j < n1; j++) {
    sc->nbr[ks0 * n1 + j] = -1;
    sc->nbr[ksn * n1 + j] = -1;
  }
}
/**********************************************************************/
/*!
 * \fn INT haz_refine_simplex(scomplex *sc, const INT is, const INT it)
 *
 * \brief Bisects simplex is with DGS conforming closure (Algorithm 3)
 *        using Algorithm 4 bisection (handles N+1 > n+1 colors).
 *        Scans ALL simplices for the conforming closure (brute force).
 *
 * \param sc  I/O: simplicial complex with gen_N set
 * \param is  I: index of the simplex to bisect
 * \param it  I: unused (kept for API compatibility), pass -1
 *
 * \return 0 on success
 */
/**********************************************************************/
/* v2s helpers: build once per refinement pass, children NOT added.  */
/* The v2s index only tracks pre-existing leaves.                    */
/* For Algorithm 4, this is correct because children of a bisected   */
/* simplex do NOT share the parent's bisection edge (the midpoint    */
/* replaces one endpoint), so they don't need closure treatment.     */
/**********************************************************************/
static inline void v2s_add(scomplex* sc, INT v, INT s) {
  if (sc->v2s_count >= sc->v2s_alloc) {
    sc->v2s_alloc = sc->v2s_alloc > 0 ? sc->v2s_alloc * 2 : 256;
    sc->v2s_next = realloc(sc->v2s_next, sc->v2s_alloc * sizeof(INT));
    sc->v2s_simp = realloc(sc->v2s_simp, sc->v2s_alloc * sizeof(INT));
  }
  INT k = sc->v2s_count++;
  sc->v2s_simp[k] = s;
  sc->v2s_next[k] = sc->v2s_head[v];
  sc->v2s_head[v] = k;
}
/**********************************************************************/
INT haz_refine_simplex(scomplex* sc, const INT is, const INT it) {
  if (is < 0 || sc->child0[is] >= 0) return 0;
  INT dim = sc->dim, n1 = dim + 1;
  INT v0, vg;
  dgs_get_bse(sc, is, &v0, &vg);
  /* Conforming closure (Algorithm 3): scan ALL simplices (including
     children created during closure) for leaves sharing edge {v0,vg}
     whose bse differs — these must be refined before we bisect {v0,vg}.
     The static v2s index misses children created during closure, so we
     use a full scan to ensure correctness. */
  {
    INT restart = 1;
    while (restart) {
      restart = 0;
      for (INT s = 0; s < sc->ns; s++) {
        if (s == is || sc->child0[s] >= 0) continue;
        INT *el = sc->nodes + s * n1;
        INT h0 = 0, hg = 0;
        for (INT a = 0; a < n1; a++) {
          if (el[a] == v0) h0 = 1;
          if (el[a] == vg) hg = 1;
        }
        if (!h0 || !hg) continue;
        INT bv0, bvg;
        dgs_get_bse(sc, s, &bv0, &bvg);
        if ((bv0 == v0 && bvg == vg) || (bv0 == vg && bvg == v0)) continue;
        haz_refine_simplex(sc, s, -1);
        restart = 1;
        break;
      }
    }
  }
  if (sc->child0[is] >= 0) return 0;
  /* Bisect is, creating new midpoint */
  INT vnew = haz_bisect_new(sc, is);
  /* Bisect ALL remaining leaf simplices sharing edge {v0,vg}.
     This includes children created during the closure phase above,
     which are not in the static v2s index. */
  {
    INT nscan = sc->ns;
    for (INT s = 0; s < nscan; s++) {
      if (s == is || sc->child0[s] >= 0) continue;
      INT *el = sc->nodes + s * n1;
      INT h0 = 0, hg = 0;
      for (INT a = 0; a < n1; a++) {
        if (el[a] == v0) h0 = 1;
        if (el[a] == vg) hg = 1;
      }
      if (h0 && hg) {
        haz_bisect_reuse(sc, s, vnew, v0, vg);
      }
    }
  }
  return 0;
}
/******************************************************************/
/*!
 * \fn void refine(const INT ref_levels, scomplex *sc,ivector *marked)
 *
 * \brief Refines a simplicial complex.
 *
 * \param sc: scomplex containing the whole hierarchy of refinements
 *
 * \param marked: input ivector containing all simplices from the
 *               finest level which are marked for refinement; If
 *               marked is null then uniformly renfines the grid
 *               ref_levels.
 *
 *
 * \param ref_levels number of refinements. If marked is not null,
 * then this is ignored and only one refinement is done;
 *
 * \return void
 *
 * \note On first call (sc->level==0), performs DGS initialization
 *       (vertex coloring + simplex reordering) from:
 *       Diening, L., Gehring, L., Storn, J. Adaptive Mesh
 *       Refinement for Arbitrary Initial Triangulations.
 *       Found. Comput. Math. (2025).
 *       DOI: 10.1007/s10208-024-09642-1
 *
 */
void refine(const INT ref_levels, scomplex* sc, ivector* marked) {
  if (ref_levels <= 0) return;
  /*somethind to be done*/
  INT j = -1, i, nsold, nsfine = -1;
  if (!sc->level) {
    /* sc->level this is set to 0 in haz_scomplex_init */
    /* form neighboring list on the coarsest level */
    find_nbr(sc->ns, sc->nv, sc->dim, sc->nodes, sc->nbr);
    /*
     * DGS initialization (Diening-Gehring-Storn, FoCM 2025):
     * Color the vertices and reorder each simplex accordingly.
     * This works for ANY conforming initial triangulation and
     * does not require the reflected-neighbor (Traxler) ordering.
     */
    {
      INT* color = (INT*)calloc(sc->nv, sizeof(INT));
      INT N = set_color(sc, color);
      /* dgs_initialize permutes both nodes and nbr arrays consistently,
         so no need to recompute neighbors. */
      dgs_initialize(sc, color, N);
      free(color);
      if (sc->print_level > 0)
        fprintf(stdout,
                "\n%%%% DGS initialization: %lld colors used (N=%lld)\n",
                (long long)(N + 1), (long long)N);
    }
    /* Traxler BFS tree initialization commented out — using DGS only */
    /* { */
    /*   INT *wrk=calloc(5*(sc->dim+2),sizeof(INT)); */
    /*   /\* construct bfs tree for the dual graph *\/ */
    /*   abfstree(0,sc,wrk,print_level); */
    /*   if(wrk) free(wrk); */
    /* } */
  }
  if ((!marked)) {
    // just refine everything that was not refined:
    for (i = 0; i < ref_levels; i++) {
      v2s_build(sc);
      nsold = sc->ns;
      for (j = 0; j < nsold; j++)
        if (sc->child0[j] < 0 && sc->childn[j] < 0)
          haz_refine_simplex(sc, j, -1);
      v2s_free_data(sc);
      sc->level++;
    }
    for (j = 0; j < sc->ns; j++)
      sc->marked[j] = TRUE;  // not sure we need this.
    // we are done here;
    return;
  } else if (!sc->level) {
    // we have not refined anything yet and marked is set, so
    if ((marked->row) && (marked->val))
      for (j = 0; j < sc->ns; j++) sc->marked[j] = marked->val[j];
    else {
      //      issue a warning and mark everything for refinement;
      for (j = 0; j < sc->ns; j++) sc->marked[j] = TRUE;
    }
  } else {
    /* we come here if we have refined few times and in such case we
       need to re-mark our simplices on the finest level using the
       values of marked at abs(sc->child0[j]+1) which were set from
       the previous visit here */
    for (j = 0; j < sc->ns; j++) {
      if (sc->child0[j] < 0 || sc->childn[j] < 0) {
        nsfine = abs(sc->child0[j] + 1);
        sc->marked[j] = marked->val[nsfine];
      }
    }
  }
  /*
   * refine everything that is marked on the finest level and is
   * not yet refined: (marked>0 and child<0)
   */
  v2s_build(sc);
  nsold = sc->ns;
  for (j = 0; j < nsold; j++)
    if (sc->marked[j] && (sc->child0[j] < 0 && sc->childn[j] < 0))
      haz_refine_simplex(sc, j, -1);
  v2s_free_data(sc);
  /*
   *  compute volumes (the volumes on the coarsest grid should be set in
   * generate_initial_grid, but just in case we are coming here from
   * some other function we compute these here too.
   */
  INT node, in1;
  void* wrk1 = malloc((sc->dim + 1) * (sc->dim * sizeof(REAL) + sizeof(INT)));
  REAL* xs = calloc((sc->dim + 1) * sc->dim, sizeof(REAL));
  for (i = 0; i < sc->ns; ++i) {
    if (sc->child0[i] < 0 || sc->childn[i] < 0) {
      in1 = i * (sc->dim + 1);
      for (j = 0; j <= sc->dim; ++j) {
        node = sc->nodes[in1 + j];
        memcpy((xs + j * sc->dim), (sc->x + node * sc->dim), sc->dim * sizeof(REAL));
      }
      sc->vols[i] = volume_compute(sc->dim, sc->factorial, xs, wrk1);
    }
  }
  free(wrk1);
  free(xs);
  sc->level++;
  return;
}
/**********************************************************************/
/*!
* \fn scomplex* make_uniform_mesh(const INT dim,const INT mesh_ref_levels,const INT mesh_ref_type,const INT set_bndry_codes)
*
* \brief makes a mesh_ref_levels of refined mesh on the unit cube in dimension dim.
*
* \param dim               Dimension of computational domain
* \param mesh_ref_levels   Number of refinement levels
*
* \param mesh_ref_type     if > 10: uniform refinement ;
*                          if < 10: nearest vertex bisection ;
*
* \param set_bndry_codes   if .eq. 0: the boundary codes come from sc->bndry[];
*                          if .ne. 0  the boundary codes are set
*
* \return scomplex* with the mesh and FEM data.
*
* \note 20210815 (ltz)
* \note Works in 2D and 3D
*
*/
scomplex* make_uniform_mesh(const INT dim,const INT mesh_ref_levels,const INT mesh_ref_type,const INT set_bndry_codes)
{

  // Loop Counters
  INT jlevel;

  // Create a simplicial complex
  scomplex **sc_all=NULL,*sc=NULL;
  fprintf(stdout,"\n%%%%---------------------------------------------------------------------");
  fprintf(stdout,"\n%%%%Meshing...");
  // Get the coarsest mesh on the cube in dimension dim and set the refinement type.
  sc_all=mesh_cube_init(dim,1,mesh_ref_type);
  sc=sc_all[0];
  if(sc->ref_type>10){
    // Uniform (Freudenthal) refinement for any dimension
    for(jlevel=0;jlevel<mesh_ref_levels;++jlevel){
      uniformrefine(sc);
      sc_vols(sc);
    }
    // Get boundaries
    find_nbr(sc->ns,sc->nv,sc->dim,sc->nodes,sc->nbr);
    sc_vols(sc);
  } else {
    // DGS bisection refinement (works for any dimension)
    refine(mesh_ref_levels, sc, NULL);
  }
  fprintf(stdout,"Done.\n");
  scfinalize(sc, NULL, (INT)1);
  sc_vols(sc);
  if (FALSE) {// (do not export)/(export) the mesh to vtu: [FALSE/TRUE]
    vtu_data vdata;
    vtu_data_init(sc, &vdata);
    sc_write_vtk("mesh.vtu", &vdata);
    vtu_data_free(&vdata);
  }
  fprintf(stdout,"%%%%---------------------------------------------------------------------\n");
  // Build FEM data directly on the simplicial complex (dim <= 3 only)
  if (dim <= 3) sc_build_fem_data(sc);
  // Free the pointer array but NOT the scomplex itself (caller owns it)
  free(sc_all);
  return sc;
}
