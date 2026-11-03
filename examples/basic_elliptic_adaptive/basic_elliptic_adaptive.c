/*! \file examples/basic_elliptic_adaptive/basic_elliptic_adaptive.c
 *
 *  Created by Ludmil Zikatanov 20261103.
 *  Copyright 2026__HAZMATH__. All rights reserved.
 *  Created with the help of Claude (Anthropic).
 *
 * \brief Adaptive P1 finite elements for -Delta u = f on the L-shaped
 *        domain (2D) and the Fichera domain (3D), with the a posteriori
 *        error estimator of
 *
 *        Y. Li and L. Zikatanov, A posteriori error estimates of finite
 *        element methods by preconditioning, Comput. Math. Appl. 91
 *        (2021) 192-201.
 *
 *        Loop: SOLVE -> ESTIMATE -> MARK (Doerfler) -> REFINE (DGS bisection).
 *
 *        Two estimators are computed on every mesh (see doc/):
 *
 *          eta  : sum over vertex patches of the energy of the local
 *                 solutions eta_k in span{element bubbles, interior face
 *                 bubbles of the patch}, eq. (4.6) of the paper.
 *          zeta : the explicit residual estimator, sum of zeta_k^2.
 *
 *        Both are split into element indicators whose sum is the
 *        estimator. Either can drive the marking.
 *
 * \note  For constant f (problems 1 and 2) everything (P1 stiffness,
 *        load, bubble matrices, residuals) is integrated exactly with
 *          int_T lambda^alpha = |T| d! alpha! / (d+|alpha|)!,
 *        so no quadrature is used in the discretization or in the
 *        estimators. Quadrature (HAZmath local data, barycentric
 *        coordinates at the quadrature points filled once) is used for
 *        the moments of f in problem 3 and for the error against the
 *        exact solution in problems 1 and 3.
 *
 * \note  The local patch problems are small SPD systems (size bounded by
 *        shape regularity) and are solved by dense LU with scaled
 *        partial pivoting (ddense_solve_pivot).
 *
 * \note  FE_type in input.dat selects the element as in basic_elliptic:
 *          1  P1 (grad):              this file;
 *          20 Nedelec (curl), 30 Raviart-Thomas (div): (eps D u, D v) +
 *             (kappa u, v) = (f, v) with the robust residual estimator of
 *             Li-Zikatanov, Math. Comp. 2025, in
 *             basic_elliptic_adaptive_hcurl_hdiv.h.
 *
 *  Usage: ./basic_elliptic_adaptive.ex [input/2d_L.input]
 *         (the solver parameters and FE_type are read from ./input.dat)
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
#include <sys/resource.h>
#include "basic_elliptic_adaptive_data.h"
/*********************************************************************/

/* Leaf mesh data needed by the exact P1/bubble integration. */
typedef struct {
  INT dim, n1, ns, nv;
  INT* nodes;  /* [ns*n1] element -> vertices (points into the scomplex) */
  REAL* x;     /* [nv*dim] coordinates (points into the scomplex) */
  REAL* vol;   /* [ns] element volumes */
  REAL* glam;  /* [ns*n1*dim] gradients of the barycentric coordinates */
  REAL* hT;    /* [ns] element diameters */
  INT* v2e_ia; /* [nv+1] vertex -> elements (CSR) */
  INT* v2e_ja;
  INT* fnbr;   /* [ns*n1] element across the face opposite local vertex i, -1 on the boundary */
  INT* fid;    /* [ns*n1] interior face number of that face, -1 on the boundary */
  INT nfi;     /* number of interior faces */
} p1mesh;

/*********************************************************************/
static REAL factorial_r(INT n) {
  REAL f = 1.0;
  for (INT i = 2; i <= n; i++) f *= (REAL)i;
  return f;
}
/* int_T prod lambda_i^alpha_i = |T| d! alpha! / (d+|alpha|)! */
static REAL int_bary(INT dim, REAL vol, INT* alpha) {
  INT i, s = 0;
  REAL p = 1.0;
  for (i = 0; i <= dim; i++) { s += alpha[i]; p *= factorial_r(alpha[i]); }
  return vol * factorial_r(dim) * p / factorial_r(dim + s);
}
/*********************************************************************/
static void p1mesh_build(p1mesh* m, scomplex* sc) {
  INT dim = sc->dim, n1 = dim + 1, ns = sc->ns, nv = sc->nv;
  INT T, i, j, k, d;
  m->dim = dim; m->n1 = n1; m->ns = ns; m->nv = nv;
  m->nodes = sc->nodes;
  m->x = sc->x;
  m->vol = (REAL*)calloc(ns, sizeof(REAL));
  m->glam = (REAL*)calloc(ns * n1 * dim, sizeof(REAL));
  m->hT = (REAL*)calloc(ns, sizeof(REAL));
  REAL* xv = (REAL*)calloc(n1 * dim, sizeof(REAL));
  REAL* refmap = (REAL*)calloc(dim * dim, sizeof(REAL));
  for (T = 0; T < ns; T++) {
    INT* el = sc->nodes + T * n1;
    for (i = 0; i < n1; i++)
      for (d = 0; d < dim; d++) xv[i * dim + d] = sc->x[el[i] * dim + d];
    compute_refelm_mapping(refmap, m->glam + T * n1 * dim, xv, dim);
    m->vol[T] = fabs(sc->vols[T]);
    REAL h2 = 0.0;
    for (i = 0; i < n1; i++)
      for (j = i + 1; j < n1; j++) {
        REAL s = 0.0;
        for (d = 0; d < dim; d++) {
          REAL t = xv[i * dim + d] - xv[j * dim + d];
          s += t * t;
        }
        if (s > h2) h2 = s;
      }
    m->hT[T] = sqrt(h2);
  }
  free(xv); free(refmap);
  /* vertex -> element incidence */
  m->v2e_ia = (INT*)calloc(nv + 1, sizeof(INT));
  for (T = 0; T < ns; T++)
    for (i = 0; i < n1; i++) m->v2e_ia[sc->nodes[T * n1 + i] + 1]++;
  for (k = 0; k < nv; k++) m->v2e_ia[k + 1] += m->v2e_ia[k];
  m->v2e_ja = (INT*)calloc(m->v2e_ia[nv], sizeof(INT));
  INT* pos = (INT*)calloc(nv, sizeof(INT));
  for (k = 0; k < nv; k++) pos[k] = m->v2e_ia[k];
  for (T = 0; T < ns; T++)
    for (i = 0; i < n1; i++) {
      k = sc->nodes[T * n1 + i];
      m->v2e_ja[pos[k]++] = T;
    }
  free(pos);
  /* face neighbors: the face opposite local vertex i of T consists of the
     other dim vertices; the neighbor is the other element containing all
     of them */
  m->fnbr = (INT*)calloc(ns * n1, sizeof(INT));
  m->fid = (INT*)calloc(ns * n1, sizeof(INT));
  for (j = 0; j < ns * n1; j++) { m->fnbr[j] = -1; m->fid[j] = -1; }
  m->nfi = 0;
  for (T = 0; T < ns; T++) {
    INT* el = sc->nodes + T * n1;
    for (i = 0; i < n1; i++) {
      if (m->fid[T * n1 + i] >= 0) continue; /* already numbered from the other side */
      INT v0 = el[(i == 0) ? 1 : 0];
      INT nb = -1, inb = -1;
      for (j = m->v2e_ia[v0]; j < m->v2e_ia[v0 + 1] && nb < 0; j++) {
        INT E = m->v2e_ja[j];
        if (E == T) continue;
        INT* eE = sc->nodes + E * n1;
        INT common = 0, a, b, miss = -1;
        for (a = 0; a < n1; a++) {
          INT inface = 0;
          if (a == i) continue;
          for (b = 0; b < n1; b++)
            if (eE[b] == el[a]) { inface = 1; break; }
          common += inface;
        }
        if (common == dim) {
          nb = E;
          /* local index in E of the vertex not on the face */
          for (b = 0; b < n1; b++) {
            INT found = 0;
            for (a = 0; a < n1; a++)
              if (a != i && el[a] == eE[b]) { found = 1; break; }
            if (!found) { miss = b; break; }
          }
          inb = miss;
        }
      }
      if (nb >= 0) {
        m->fnbr[T * n1 + i] = nb;
        m->fnbr[nb * n1 + inb] = T;
        m->fid[T * n1 + i] = m->nfi;
        m->fid[nb * n1 + inb] = m->nfi;
        m->nfi++;
      }
    }
  }
  return;
}
static void p1mesh_free(p1mesh* m) {
  free(m->vol); free(m->glam); free(m->hT);
  free(m->v2e_ia); free(m->v2e_ja);
  free(m->fnbr); free(m->fid);
}
/*********************************************************************/
/* Element moments of the right-hand side f:
     fl[T*n1+i] = int_T f lambda_i      (P1 load)
     fb[T*np+p] = int_T f psi_p         (bubble loads; psi_0=b_T, psi_{1+m}=b_{F_m})
     f2[T] = int_T f^2,  f1[T] = int_T f.
   Constant f (problems 1, 2): exact formula, no quadrature.
   Problem 3: quadrature; the barycentric coordinates at the reference
   quadrature points are filled once (initialize_localdata_elm) and only the
   element map and the weights are updated per element (get_elmlocaldata). */
typedef struct {
  REAL *fl, *fb, *f2, *f1;
} rhs_moments;
static INT bubble_has(INT p, INT i) { return (p == 0) || (i != p - 1); }
static void rhs_moments_build(rhs_moments* r, p1mesh* m, scomplex* leaf, fespace* FE,
                              INT problem, REAL fconst, INT nq1d) {
  INT dim = m->dim, n1 = m->n1, ns = m->ns, np = dim + 2;
  INT T, i, p, q;
  r->fl = (REAL*)calloc(ns * n1, sizeof(REAL));
  r->fb = (REAL*)calloc(ns * np, sizeof(REAL));
  r->f2 = (REAL*)calloc(ns, sizeof(REAL));
  r->f1 = (REAL*)calloc(ns, sizeof(REAL));
  INT* alpha = (INT*)calloc(n1, sizeof(INT));
  if (problem != 3) {
    for (T = 0; T < ns; T++) {
      REAL vol = m->vol[T];
      for (i = 0; i < n1; i++) r->fl[T * n1 + i] = fconst * vol / (REAL)n1;
      for (p = 0; p < np; p++) {
        for (i = 0; i < n1; i++) alpha[i] = bubble_has(p, i);
        r->fb[T * np + p] = fconst * int_bary(dim, vol, alpha);
      }
      r->f2[T] = fconst * fconst * vol;
      r->f1[T] = fconst * vol;
    }
    free(alpha);
    return;
  }
  /* quadrature with the HAZmath local data structures */
  block_fespace bfe;
  memset(&bfe, 0, sizeof(block_fespace));
  bfe.nspaces = 1; bfe.nun = 1; bfe.ndof = FE->ndof; bfe.nelm = FE->nelm;
  bfe.var_spaces = (fespace**)calloc(1, sizeof(fespace*));
  bfe.var_spaces[0] = FE;
  simplex_local_data elm;
  fe_local_data fed;
  memset(&elm, 0, sizeof(simplex_local_data));
  memset(&fed, 0, sizeof(fe_local_data));
  initialize_localdata_elm(&elm, &fed, leaf, &bfe, nq1d);
  REAL wmax_err = 0.0;
  for (T = 0; T < ns; T++) {
    get_elmlocaldata(&elm, leaf, T);
    INT nq = elm.quad_local->nq;
    REAL wsum = 0.0;
    for (q = 0; q < nq; q++) {
      REAL* xq = elm.quad_local->x + q * dim;
      REAL w = elm.quad_local->w[q];
      REAL* lam = elm.lams + q * n1; /* reference values, filled once */
      REAL fq;
      rhs_sines(&fq, xq, 0.0, NULL);
      wsum += w;
      for (i = 0; i < n1; i++) r->fl[T * n1 + i] += w * fq * lam[i];
      for (p = 0; p < np; p++) {
        REAL b = 1.0;
        for (i = 0; i < n1; i++) if (bubble_has(p, i)) b *= lam[i];
        r->fb[T * np + p] += w * fq * b;
      }
      r->f2[T] += w * fq * fq;
      r->f1[T] += w * fq;
    }
    REAL rel = fabs(wsum - m->vol[T]) / m->vol[T];
    if (rel > wmax_err) wmax_err = rel;
  }
  if (wmax_err > 1e-12)
    fprintf(stderr, "%% WARNING: quadrature weights do not sum to |T| (max rel. err %.2e)\n", wmax_err);
  {
    INT nbad = 0, first = -1;
    for (T = 0; T < ns; T++) {
      INT ok = isfinite(r->f2[T]) && isfinite(r->f1[T]);
      for (i = 0; i < n1; i++) ok = ok && isfinite(r->fl[T * n1 + i]);
      for (p = 0; p < np; p++) ok = ok && isfinite(r->fb[T * np + p]);
      if (!ok) { nbad++; if (first < 0) first = T; }
    }
    if (nbad)
      fprintf(stderr, "%% WARNING: non-finite rhs moments on %lld elements (first %lld)\n",
              (long long)nbad, (long long)first);
  }
  free_simplexlocaldata(&elm);
  free_felocaldata(&fed);
  free(bfe.var_spaces);
  free(alpha);
}
static void rhs_moments_free(rhs_moments* r) {
  free(r->fl); free(r->fb); free(r->f2); free(r->f1);
}
/*********************************************************************/
/* P1 stiffness matrix (exact) and load vector from the moments:
   A_ij = sum_T |T| grad(lambda_i).grad(lambda_j),  b_i = sum_T int_T f lambda_i. */
static void assemble_p1_exact(dCSRmat* A, dvector* b, p1mesh* m, REAL* fl) {
  INT dim = m->dim, n1 = m->n1, nv = m->nv;
  INT k, j, a, c, d, T;
  INT* mark = (INT*)calloc(nv, sizeof(INT));
  for (k = 0; k < nv; k++) mark[k] = -1;
  /* sparsity pattern from the vertex -> element incidence */
  INT* ia = (INT*)calloc(nv + 1, sizeof(INT));
  for (k = 0; k < nv; k++) {
    INT cnt = 0;
    for (j = m->v2e_ia[k]; j < m->v2e_ia[k + 1]; j++) {
      INT* el = m->nodes + m->v2e_ja[j] * n1;
      for (a = 0; a < n1; a++)
        if (mark[el[a]] != k) { mark[el[a]] = k; cnt++; }
    }
    ia[k + 1] = ia[k] + cnt;
  }
  INT nnz = ia[nv];
  dcsr_alloc(nv, nv, nnz, A);
  memcpy(A->IA, ia, (nv + 1) * sizeof(INT));
  free(ia);
  memset(A->val, 0, nnz * sizeof(REAL));
  for (k = 0; k < nv; k++) mark[k] = -1;
  dvec_alloc(nv, b);
  dvec_set(nv, b, 0.0);
  for (k = 0; k < nv; k++) {
    INT p = A->IA[k];
    /* columns of row k; mark[] holds the position in the row */
    for (j = m->v2e_ia[k]; j < m->v2e_ia[k + 1]; j++) {
      INT* el = m->nodes + m->v2e_ja[j] * n1;
      for (a = 0; a < n1; a++)
        if (mark[el[a]] < A->IA[k]) { mark[el[a]] = p; A->JA[p++] = el[a]; }
    }
    for (j = m->v2e_ia[k]; j < m->v2e_ia[k + 1]; j++) {
      T = m->v2e_ja[j];
      INT* el = m->nodes + T * n1;
      REAL* gl = m->glam + T * n1 * dim;
      INT ik = -1;
      for (a = 0; a < n1; a++) if (el[a] == k) ik = a;
      for (c = 0; c < n1; c++) {
        REAL s = 0.0;
        for (d = 0; d < dim; d++) s += gl[ik * dim + d] * gl[c * dim + d];
        A->val[mark[el[c]]] += m->vol[T] * s;
      }
      b->val[k] += fl[T * n1 + ik];
    }
  }
  free(mark);
  return;
}
/*********************************************************************/
/* Computes the element indicators etaT2[T] (local bubble problems),
   zetaT2[T] (explicit residual) and the oscillation oscT2[T] =
   h_T^2 ||f - mean_T f||_T^2, from the right-hand-side moments r. */
static void estimate(p1mesh* m, REAL* u, rhs_moments* r, REAL* etaT2, REAL* zetaT2,
                    REAL* oscT2, REAL* eta2_patches) {
  INT dim = m->dim, n1 = m->n1, ns = m->ns, nv = m->nv;
  INT np = dim + 2; /* b_T and the dim+1 face bubbles of T */
  INT T, i, l, p, q, d, k, j, a;
  REAL* Ke = (REAL*)calloc(ns * np * np, sizeof(REAL));
  REAL* le = (REAL*)calloc(ns * np, sizeof(REAL));
  REAL* gu = (REAL*)calloc(ns * dim, sizeof(REAL));
  INT* S = (INT*)calloc(np * n1, sizeof(INT)); /* S[p*n1+i]=1 if lambda_i is a factor of psi_p */
  INT* alpha = (INT*)calloc(n1, sizeof(INT));
  for (p = 0; p < np; p++)
    for (i = 0; i < n1; i++) S[p * n1 + i] = (p == 0) ? 1 : (i != p - 1);
  /* element matrices, loads and gradients of u_h */
  for (T = 0; T < ns; T++) {
    INT* el = m->nodes + T * n1;
    REAL* gl = m->glam + T * n1 * dim;
    REAL vol = m->vol[T];
    REAL* g = gu + T * dim;
    for (d = 0; d < dim; d++) {
      g[d] = 0.0;
      for (i = 0; i < n1; i++) g[d] += u[el[i]] * gl[i * dim + d];
    }
    for (p = 0; p < np; p++) {
      /* load: int f psi_p - int grad u_h . grad psi_p */
      for (i = 0; i < n1; i++) alpha[i] = S[p * n1 + i];
      REAL lp = r->fb[T * np + p];
      for (i = 0; i < n1; i++) {
        if (!S[p * n1 + i]) continue;
        REAL gdl = 0.0;
        for (d = 0; d < dim; d++) gdl += g[d] * gl[i * dim + d];
        alpha[i] = 0;
        lp -= gdl * int_bary(dim, vol, alpha);
        alpha[i] = 1;
      }
      le[T * np + p] = lp;
      /* matrix */
      for (q = p; q < np; q++) {
        REAL s = 0.0;
        for (i = 0; i < n1; i++) {
          if (!S[p * n1 + i]) continue;
          for (l = 0; l < n1; l++) {
            if (!S[q * n1 + l]) continue;
            REAL Gil = 0.0;
            for (d = 0; d < dim; d++) Gil += gl[i * dim + d] * gl[l * dim + d];
            for (j = 0; j < n1; j++) alpha[j] = S[p * n1 + j] + S[q * n1 + j];
            alpha[i]--; alpha[l]--;
            s += Gil * int_bary(dim, vol, alpha);
          }
        }
        Ke[(T * np + p) * np + q] = s;
        Ke[(T * np + q) * np + p] = s;
      }
    }
  }
  /* explicit residual indicators:
     zeta_T^2 = (d+1) h_T^2 ||f||_T^2 + (d/2) sum_{interior F of T} h_F |r_F|^2 |F| */
  for (T = 0; T < ns; T++) {
    INT* el = m->nodes + T * n1;
    REAL* gl = m->glam + T * n1 * dim;
    REAL h2 = m->hT[T] * m->hT[T];
    REAL o = h2 * (r->f2[T] - r->f1[T] * r->f1[T] / m->vol[T]);
    oscT2[T] = (o > 0.0) ? o : 0.0;
    REAL z = (REAL)n1 * h2 * r->f2[T];
    for (i = 0; i < n1; i++) {
      INT E = m->fnbr[T * n1 + i];
      if (E < 0) continue;
      REAL ng = 0.0;
      for (d = 0; d < dim; d++) ng += gl[i * dim + d] * gl[i * dim + d];
      ng = sqrt(ng);
      REAL areaF = (REAL)dim * m->vol[T] * ng; /* |F_i| = d |T| |grad lambda_i| */
      REAL rF = 0.0;                         /* -(grad u_T - grad u_E).n_T, n_T=-grad lambda_i/|.| */
      for (d = 0; d < dim; d++)
        rF += (gu[T * dim + d] - gu[E * dim + d]) * gl[i * dim + d] / ng;
      REAL hF2 = 0.0;
      for (a = 0; a < n1; a++)
        for (j = a + 1; j < n1; j++) {
          if (a == i || j == i) continue;
          REAL s = 0.0;
          for (d = 0; d < dim; d++) {
            REAL t = m->x[el[a] * dim + d] - m->x[el[j] * dim + d];
            s += t * t;
          }
          if (s > hF2) hF2 = s;
        }
      if (dim == 1) hF2 = 0.0; /* a face is a point */
      z += 0.5 * (REAL)dim * sqrt(hF2) * rF * rF * areaF;
    }
    zetaT2[T] = z;
  }
  /* local problems on the vertex patches */
  INT maxne = 0;
  for (k = 0; k < nv; k++)
    if (m->v2e_ia[k + 1] - m->v2e_ia[k] > maxne) maxne = m->v2e_ia[k + 1] - m->v2e_ia[k];
  INT maxm = maxne * (1 + dim);
  REAL* Aloc = (REAL*)calloc(maxm * maxm, sizeof(REAL));
  REAL* bloc = (REAL*)calloc(maxm, sizeof(REAL));
  REAL* cloc = (REAL*)calloc(maxm, sizeof(REAL));
  REAL* piv = (REAL*)calloc(maxm, sizeof(REAL)); /* work space for the LU */
  INT* perm = (INT*)calloc(maxm, sizeof(INT));   /* row permutation of the LU */
  INT* flist = (INT*)calloc(maxm, sizeof(INT));
  INT* lmap = (INT*)calloc(maxne * np, sizeof(INT));
  REAL* cT = (REAL*)calloc(np, sizeof(REAL));
  *eta2_patches = 0.0;
  for (T = 0; T < ns; T++) etaT2[T] = 0.0;
  for (k = 0; k < nv; k++) {
    INT ja0 = m->v2e_ia[k], ne = m->v2e_ia[k + 1] - ja0, nf = 0;
    /* local numbering: elements 0..ne-1, then the interior faces through x_k */
    for (a = 0; a < ne; a++) {
      T = m->v2e_ja[ja0 + a];
      INT* el = m->nodes + T * n1;
      INT jk = -1;
      for (i = 0; i < n1; i++) if (el[i] == k) jk = i;
      lmap[a * np + 0] = a;
      for (i = 0; i < n1; i++) {
        INT f = m->fid[T * n1 + i];
        if (i == jk || f < 0) { lmap[a * np + 1 + i] = -1; continue; }
        INT ff = -1;
        for (j = 0; j < nf; j++) if (flist[j] == f) { ff = j; break; }
        if (ff < 0) { ff = nf; flist[nf++] = f; }
        lmap[a * np + 1 + i] = ne + ff;
      }
    }
    INT msz = ne + nf;
    memset(Aloc, 0, msz * msz * sizeof(REAL));
    memset(bloc, 0, msz * sizeof(REAL));
    for (a = 0; a < ne; a++) {
      T = m->v2e_ja[ja0 + a];
      for (p = 0; p < np; p++) {
        INT ip = lmap[a * np + p];
        if (ip < 0) continue;
        bloc[ip] += le[T * np + p];
        for (q = 0; q < np; q++) {
          INT iq = lmap[a * np + q];
          if (iq < 0) continue;
          Aloc[ip * msz + iq] += Ke[(T * np + p) * np + q];
        }
      }
    }
    /* dense LU with scaled partial pivoting; Aloc is overwritten by the
       factors and cloc (= bloc on entry) by the solution */
    memcpy(cloc, bloc, msz * sizeof(REAL));
    ddense_solve_pivot((INT)1, msz, Aloc, cloc, perm, piv);
    REAL ek = 0.0;
    for (j = 0; j < msz; j++) ek += cloc[j] * bloc[j];
    *eta2_patches += ek;
    /* split ||grad eta_k||^2 over the elements of the patch */
    for (a = 0; a < ne; a++) {
      T = m->v2e_ja[ja0 + a];
      for (p = 0; p < np; p++) cT[p] = (lmap[a * np + p] >= 0) ? cloc[lmap[a * np + p]] : 0.0;
      REAL s = 0.0;
      for (p = 0; p < np; p++)
        for (q = 0; q < np; q++) s += cT[p] * Ke[(T * np + p) * np + q] * cT[q];
      etaT2[T] += s;
    }
  }
  free(Ke); free(le); free(gu); free(S); free(alpha);
  free(Aloc); free(bloc); free(cloc); free(piv); free(perm); free(flist); free(lmap); free(cT);
}
/*********************************************************************/
/* Doerfler marking: shortest set of largest indicators with
   sum_{marked} mu_T^2 >= theta sum mu_T^2 */
static REAL* g_sortkey = NULL;
static int cmp_desc(const void* a, const void* b) {
  REAL x = g_sortkey[*(const INT*)a], y = g_sortkey[*(const INT*)b];
  return (x < y) - (x > y);
}
static ivector dorfler_mark(INT ns, REAL* mu2, REAL theta, INT* nmarked) {
  INT T;
  ivector marked = ivec_create(ns);
  INT* idx = (INT*)calloc(ns, sizeof(INT));
  REAL total = 0.0, acc = 0.0;
  for (T = 0; T < ns; T++) { idx[T] = T; total += mu2[T]; marked.val[T] = 0; }
  g_sortkey = mu2;
  qsort(idx, ns, sizeof(INT), cmp_desc);
  g_sortkey = NULL;
  *nmarked = 0;
  for (T = 0; T < ns && acc < theta * total; T++) {
    marked.val[idx[T]] = 1;
    acc += mu2[idx[T]];
    (*nmarked)++;
  }
  free(idx);
  return marked;
}
/*********************************************************************/
static REAL cfg_real(config_z* c, const char* name, REAL dflt) {
  for (INT i = 0; i < c->num_fields; i++) {
    field_z* f = &c->fields[i];
    if (strcmp(f->name, name)) continue;
    if (f->type == TYPE_INTEGER) return (REAL)f->value.int_val;
    if (f->type == TYPE_DOUBLE) return f->value.dbl_val;
  }
  return dflt;
}
/*********************************************************************/
/* FE_type 20 (Ned0, curl) and 30 (RT0, div): uses dorfler_mark above */
#include "basic_elliptic_adaptive_hcurl_hdiv.h"
/****** MAIN DRIVER **************************************************/
int main(int argc, char* argv[]) {
  const char* grid_file = "input/2d_L.input";
  if (argc > 1) grid_file = argv[1];
  clock_t clk_overall_start = clock();

  /* solver parameters and the finite element (FE_type 1, 20 or 30) */
  input_param inparam;
  param_input_init(&inparam);
  param_input("./input.dat", &inparam);
  INT fet = inparam.FE_type;
  if (fet != 1 && fet != 20 && fet != 30) {
    fprintf(stderr, "ERROR: FE_type must be 1 (P1, grad), 20 (Nedelec, curl) or 30 (Raviart-Thomas, div)\n");
    exit(1);
  }
  printf("\n===========================================================================\n");
  if (fet == 1)
    printf("Adaptive P1 FEM for -Delta u = f (estimator by preconditioning, Li-Zikatanov 2021)\n");
  else
    printf("Adaptive %s FEM for (eps %s u, %s v) + (kappa u, v) = (f, v)\n(robust residual estimator, Li-Zikatanov 2025)\n",
           (fet == 20) ? "Nedelec" : "Raviart-Thomas", (fet == 20) ? "curl" : "div", (fet == 20) ? "curl" : "div");
  printf("===========================================================================\n");

  /* grid and adaptivity parameters */
  config_z config = get_input(grid_file);
  input_grid* g = input_grid_alloc();
  config2vars_amr(&config, g);
  INT problem = (INT)cfg_real(&config, "ae_problem", 1.0);
  REAL theta = cfg_real(&config, "ae_theta", 0.5);
  INT mark_with = (INT)cfg_real(&config, "ae_estimator", 0.0); /* 0: eta (local problems); 1: zeta (residual) */
  REAL Eref = cfg_real(&config, "ae_reference_energy", -1.0);  /* problem 2: int_Omega u */
  INT max_dofs = (INT)cfg_real(&config, "ae_max_dofs", 1.0e6);
  INT dump_meshes = (INT)cfg_real(&config, "ae_dump_meshes", 0.0); /* 1: mesh_stage_XXX.vtu at every stage */
  REAL eps = cfg_real(&config, "ae_eps", 1.0);   /* FE_type 20/30 only */
  REAL kap = cfg_real(&config, "ae_kappa", 1.0); /* FE_type 20/30 only */
  free_config(&config);
  input_grid_set_output(g, grid_file);
  INT nstages = g->nref;
  REAL fconst = rhs_const_lshape(problem);
  /* problems 2 and 3 have u = 0 on the boundary */
  void (*bcfun)(REAL*, REAL*, REAL, void*) = (problem == 1) ? bc_lshape_p1 : bc_lshape_p2;

  if (fet != 1) {
    if (problem != 2 && problem != 3) {
      fprintf(stderr, "ERROR: FE_type %lld: use ae_problem{2} (f = (1,...,1)) or ae_problem{3} (smooth exact solution)\n",
              (long long)fet);
      exit(1);
    }
    printf(" grid file: %s\n problem = %lld, eps = %g, kappa = %g, Doerfler theta = %g, nq1d = %lld\n",
           grid_file, (long long)problem, eps, kap, theta, (long long)inparam.nquad);
    scomplex** sc_all = generate_initial_grid(g);
    run_hcurl_hdiv(sc_all[0], &inparam, fet, problem, eps, kap, theta, Eref, max_dofs, nstages, dump_meshes);
    input_grid_free(g);
    haz_scomplex_free(sc_all[0]);
    free(sc_all);
    printf("\nEnd of Program: Total CPU Time = %f seconds.\n",
           (REAL)(clock() - clk_overall_start) / CLOCKS_PER_SEC);
    return 0;
  }

  if (problem == 3)
    printf(" grid file: %s\n problem = 3 (u = prod sin(pi x_i), f = d pi^2 u, quadrature nq1d = %lld)\n",
           grid_file, (long long)inparam.nquad);
  else
    printf(" grid file: %s\n problem = %lld (f = %g)\n", grid_file, (long long)problem, fconst);
  printf(" Doerfler theta = %g, marking with %s\n", theta, mark_with ? "zeta (residual)" : "eta (local problems)");
  if (problem == 2 && Eref <= 0.0)
    printf(" ae_reference_energy not set: the error column is not computed; (f,u_h) is printed\n");

  scomplex** sc_all = generate_initial_grid(g);
  scomplex* sc = sc_all[0];
  INT dim = sc->dim;
  if (problem == 1 && dim != 2) {
    fprintf(stderr, "ERROR: ae_problem{1} (exact L-shape solution) is 2D only; use ae_problem{2} or {3}\n");
    exit(1);
  }
  g_dim_sines = dim;

  /* convergence table */
  char ftab[MAXFILENAMESIZE];
  snprintf(ftab, sizeof(ftab), "%sconvergence_problem%lld.txt", inparam.output_dir, (long long)problem);
  FILE* fptab = fopen(ftab, "w");
  const char* hdr0 = "%% stage k of the adaptive loop SOLVE -> ESTIMATE -> MARK -> REFINE: mesh after k refinements (stage 0 = initial mesh)\n";
  const char* hdr = "%% stage      ndof        ns      err_A           eta            zeta       eta/err  zeta/err   (f,u_h)                 marked   osc\n";
  printf("\n"); printf(hdr0, ""); printf(hdr, "");
  if (fptab) { fprintf(fptab, hdr0, ""); fprintf(fptab, hdr, ""); }

  linear_itsolver_param linear_itparam;
  param_linear_solver_set(&linear_itparam, &inparam);
  /* The UA AMG setup changes param->strong_coupled in place (it doubles or
     halves it per level). amgparam0 keeps the input values and is copied
     into amgparam at every stage of the adaptive loop; otherwise the
     threshold drifts over the stages until no coupling is strong and AMG
     stops at a single level. */
  AMG_param amgparam0, amgparam;
  param_amg_init(&amgparam0);
  param_amg_set(&amgparam0, &inparam);
  if (inparam.print_level < 3) { /* keep the output of each stage to one line */
    linear_itparam.linear_print_level = 0;
    amgparam0.print_level = 0;
  }

  /*==================================================================*/
  /* ADAPTIVE LOOP:  SOLVE -> ESTIMATE -> MARK -> REFINE               */
  /* stage k works on the mesh obtained after k refinements           */
  /*==================================================================*/
  INT stage;
  for (stage = 0; stage <= nstages; stage++) {
    amgparam = amgparam0;
    /* leaf mesh with FEM data; the hierarchy in sc is kept for refinement */
    /* CPU time per phase (printed for print_level > 1) */
    clock_t tph[10];
    INT nph = 0;
    tph[nph++] = clock();
    scomplex* leaf = (scomplex*)malloc(sizeof(scomplex));
    scfinalize(sc, leaf, (INT)1);
    sc_vols(leaf);
    if (dump_meshes) {
      char fvtu[MAXFILENAMESIZE];
      snprintf(fvtu, sizeof(fvtu), "%smesh_stage_%03lld.vtu", inparam.output_dir, (long long)stage);
      vtu_data vdata;
      vtu_data_init(leaf, &vdata);
      sc_write_vtk(fvtu, &vdata);
      vtu_data_free(&vdata);
    }
    tph[nph++] = clock();
    sc_build_fem_data(leaf);
    tph[nph++] = clock();

    fespace FE;
    create_fespace(&FE, leaf, 1);
    set_dirichlet_bdry(&FE, leaf, 1, 16);
    if (FE.ndof != leaf->nv) {
      fprintf(stderr, "ERROR: P1 dofs (%lld) != vertices (%lld)\n", (long long)FE.ndof, (long long)leaf->nv);
      exit(1);
    }
    tph[nph++] = clock();

    p1mesh m;
    p1mesh_build(&m, leaf);
    rhs_moments rm;
    rhs_moments_build(&rm, &m, leaf, &FE, problem, fconst, inparam.nquad);

    /*---------------- 1. SOLVE ---------------------------------------*/
    dCSRmat A;
    dvector b;
    assemble_p1_exact(&A, &b, &m, rm.fl);
    tph[nph++] = clock();
    eliminate_DirichletBC(bcfun, &FE, leaf, &b, &A, 0.0);
    tph[nph++] = clock();
    dvector u = dvec_create(FE.ndof);
    dvec_set(u.row, &u, 0.0);
    INT solver_flag;
    if (linear_itparam.linear_itsolver_type == 0) {
      solver_flag = directsolve_HAZ(&A, &b, &u, 0);
    } else if (linear_itparam.linear_itsolver_type == SOLVER_AMG) {
      solver_flag = linear_solver_amg(&A, &b, &u, &amgparam);
    } else if (linear_itparam.linear_precond_type == PREC_AMG) {
      solver_flag = linear_solver_dcsr_krylov_amg(&A, &b, &u, &linear_itparam, &amgparam);
    } else if (linear_itparam.linear_precond_type == PREC_DIAG) {
      solver_flag = linear_solver_dcsr_krylov_diag(&A, &b, &u, &linear_itparam);
    } else {
      solver_flag = linear_solver_dcsr_krylov(&A, &b, &u, &linear_itparam);
    }
    if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %lld!\n", (long long)solver_flag);
    {
      INT nnan = 0;
      for (INT jj = 0; jj < u.row; jj++) if (!isfinite(u.val[jj])) nnan++;
      if (nnan) {
        fprintf(stderr, "### ERROR: %lld non-finite entries in the discrete solution (stage %lld, ndof %lld); stopping.\n",
                (long long)nnan, (long long)stage, (long long)FE.ndof);
        exit(2);
      }
    }
    tph[nph++] = clock();

    /*---------------- 2. ESTIMATE ------------------------------------*/
    REAL* etaT2 = (REAL*)calloc(m.ns, sizeof(REAL));
    REAL* zetaT2 = (REAL*)calloc(m.ns, sizeof(REAL));
    REAL* oscT2 = (REAL*)calloc(m.ns, sizeof(REAL));
    REAL eta2_patches;
    estimate(&m, u.val, &rm, etaT2, zetaT2, oscT2, &eta2_patches);
    REAL eta2 = 0.0, zeta2 = 0.0, osc2 = 0.0;
    INT T, i;
    for (T = 0; T < m.ns; T++) { eta2 += etaT2[T]; zeta2 += zetaT2[T]; osc2 += oscT2[T]; }
    tph[nph++] = clock();
    if (inparam.print_level > 2) {
      /* checks: element indicators sum to the patch energies; mesh is conforming */
      printf("%%   sum_T eta_T^2 = %.15e, sum_k ||eta_k||^2 = %.15e, rel.diff = %.2e\n",
             eta2, eta2_patches, fabs(eta2 - eta2_patches) / (eta2 + 1e-300));
      INT nbad = sc_conformity_check(leaf);
      printf("%%   conformity check: %s (%lld bad facets)\n", nbad ? "FAILED" : "passed", (long long)nbad);
    }

    /* (f,u_h) = sum_T sum_i u_i int_T f lambda_i */
    REAL fuh = 0.0;
    for (T = 0; T < m.ns; T++)
      for (i = 0; i < m.n1; i++) fuh += u.val[m.nodes[T * m.n1 + i]] * rm.fl[T * m.n1 + i];

    /* ERROR (verification only) */
    REAL errA = -1.0;
    if (problem == 1 || problem == 3) {
      qcoordinates* cq = get_quadrature(leaf, inparam.nquad);
      errA = HDsemierror(u.val, (problem == 1) ? D_exactsol_lshape : D_exactsol_sines, &FE, leaf, cq, 0.0);
      free_qcoords(cq); free(cq);
    } else if (Eref > 0.0) {
      REAL e2 = Eref - fuh;
      errA = (e2 > 0.0) ? sqrt(e2) : -sqrt(-e2);
    }
    tph[nph++] = clock();

    /*---------------- 3. MARK (Doerfler) -----------------------------*/
    INT nmarked = 0;
    ivector marked;
    INT last = (stage == nstages) || (FE.ndof > max_dofs);
    if (!last) marked = dorfler_mark(m.ns, mark_with ? zetaT2 : etaT2, theta, &nmarked);

    const char* fmt = "%7lld %10lld %10lld  %.6e  %.6e  %.6e  %7.4f  %7.4f   %.15e %8lld  %.3e\n";
    REAL osc = sqrt(osc2);
    REAL eta = sqrt(eta2), zeta = sqrt(zeta2);
    REAL effe = (errA > 0.0) ? eta / errA : 0.0, effz = (errA > 0.0) ? zeta / errA : 0.0;
    printf(fmt, (long long)stage, (long long)FE.ndof, (long long)m.ns, errA, eta, zeta, effe, effz, fuh, (long long)nmarked, osc);
    if (fptab) { fprintf(fptab, fmt, (long long)stage, (long long)FE.ndof, (long long)m.ns, errA, eta, zeta, effe, effz, fuh, (long long)nmarked, osc); fflush(fptab); }

    if (last) {
      char fsol[MAXFILENAMESIZE];
      snprintf(fsol, sizeof(fsol), "%ssol_problem%lld.vtu", inparam.output_dir, (long long)problem);
      dump_sol_vtk(fsol, "u", leaf, &FE, u.val);
      char fmsh[MAXFILENAMESIZE];
      snprintf(fmsh, sizeof(fmsh), "%smesh_problem%lld.msh", inparam.output_dir, (long long)problem);
      sc_write_gmsh(fmsh, leaf, 0);
      printf("\n solution: %s\n mesh:     %s\n table:    %s\n", fsol, fmsh, ftab);
    } else {
      /*-------------- 4. REFINE (DGS bisection) -----------------------*/
      /* scfinest numbers the leaves (same order as scfinalize) for refine() */
      scomplex* sctop = scfinest(sc);
      haz_scomplex_free(sctop);
      refine(1, sc, &marked);
      ivec_free(&marked);
    }
    tph[nph++] = clock();
    if (inparam.print_level > 1) {
      static const char* phname[] = {"scfinalize+vols", "sc_build_fem_data", "fespace+dirichlet",
                                     "p1mesh+assembly", "eliminate_BC", "solve", "estimate",
                                     "error", "mark+refine/output"};
      printf("%%   CPU seconds:");
      for (INT ip = 1; ip < nph; ip++)
        printf(" %s=%.3f", phname[ip - 1], (REAL)(tph[ip] - tph[ip - 1]) / CLOCKS_PER_SEC);
      printf("\n");
    }

    free(etaT2); free(zetaT2); free(oscT2);
    rhs_moments_free(&rm);
    dvec_free(&u); dvec_free(&b); dcsr_free(&A);
    p1mesh_free(&m);
    free_fespace(&FE);
    haz_scomplex_free(leaf);
    if (last) break;
  }
  if (fptab) fclose(fptab);

  input_grid_free(g);
  haz_scomplex_free(sc);
  free(sc_all);

  clock_t clk_overall_end = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n",
         (REAL)(clk_overall_end - clk_overall_start) / CLOCKS_PER_SEC);
  {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
    double peak_mb = (double)ru.ru_maxrss / (1024.0 * 1024.0);
#else
    double peak_mb = (double)ru.ru_maxrss / 1024.0;
#endif
    printf("Peak RSS = %.2f MB\n\n", peak_mb);
  }
  (void)dim;
  return 0;
}
/*******************************************************************/
