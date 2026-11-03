/*! \file examples/basic_elliptic_adaptive/basic_elliptic_adaptive_hcurl_hdiv.h
 *
 *  Created by Ludmil Zikatanov 20261103.
 *  Copyright 2026__HAZMATH__. All rights reserved.
 *  Created with the help of Claude (Anthropic).
 *
 * \brief FE_type 20 (lowest order Nedelec) and 30 (lowest order
 *        Raviart-Thomas) in 2D and 3D, for
 *
 *   curl:  (eps curl u, curl v) + (kappa u, v) = (f, v),  u x n = 0,
 *   div:   (eps div u, div v)   + (kappa u, v) = (f, v),  u . n = 0,
 *
 *   with the robust residual estimator of Y. Li and L. Zikatanov, Nodal
 *   auxiliary a posteriori error estimates, Math. Comp. 96 (2027) 47-75,
 *   Thm 5.8 and (5.19)-(5.20):
 *
 *   eta_T^2 = kappa^{-1} h_T^2 ||R1||_T^2 + hbar_T^2 ||R2||_T^2
 *           + sum_{S in dT interior} kappa^{-1} h_S ||J1||_S^2 + eps^{-1/2} hbar_S ||J2||_S^2,
 *   hbar = min(h / sqrt(eps), 1 / sqrt(kappa)),  h_T = diam(T),  h_S = diam(S),
 *
 *   curl:  R1 = -div(f - kappa u_h),  J1 = [f - kappa u_h] . n,
 *          R2 = f - kappa u_h,        J2 = [eps curl u_h] x n  (2D: the jump of eps curl u_h),
 *   div:   R1 = curl(f - kappa u_h),  J1 = [f - kappa u_h] x n  (2D: tangential component),
 *          R2 = f - kappa u_h,        J2 = [eps div u_h].
 *
 *   Each interior face enters the indicators of both of its elements.
 *
 *   The mass and curl-curl / div-div matrices are computed exactly from the
 *   barycentric gradients; the right-hand side, the error and the element
 *   residuals use quadrature with the HAZmath local data (barycentric
 *   coordinates at the reference quadrature points filled once); the face
 *   residuals are exact (u_h is affine and D u_h constant on each element).
 *
 *   Bases (the HAZmath conventions, so that the HX operators apply):
 *     Ned0: phi_E = |E| (lambda_a grad lambda_b - lambda_b grad lambda_a), a < b
 *           (global vertex numbers);
 *     RT0:  phi_F = sigma_F |F| (x - x_opp), sigma_F = (d-1)! det[grad lambda
 *           at the vertices of F, in the order of f_v].
 *
 *   Global solve with 10000 or more dofs: CG + HX (curl) for Ned0, CG + HX
 *   (div) for RT0 in 3D; the sparse direct solver for smaller systems and for
 *   RT0 in 2D. linear_itsolver_type = 0 selects the direct solver always.
 */

typedef struct {
  INT dim, n1, fet, nl, nd; /* nl local dofs; nd components of D u */
  INT v[4];
  REAL xv[12];
  REAL g[12];               /* g[i*dim+j] = d lambda_i / d x_j */
  REAL vol, hT;
  INT dof[6];
  INT a[6], b[6];           /* Ned0: oriented edge; RT0: a = opposite vertex */
  REAL len[6];              /* Ned0: |E|; RT0: sigma_F |F| */
  REAL D[18];               /* D phi_k (nd components each) */
} hd_elm;

static REAL hd_dot(const REAL* x, const REAL* y, INT n) {
  REAL s = 0.0;
  for (INT i = 0; i < n; i++) s += x[i] * y[i];
  return s;
}
static void hd_cross3(const REAL* a, const REAL* b, REAL* c) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}
static REAL hd_det(const REAL* M, INT dim) { /* row major, dim = 2 or 3 */
  if (dim == 2) return M[0] * M[3] - M[1] * M[2];
  return M[0] * (M[4] * M[8] - M[5] * M[7]) - M[1] * (M[3] * M[8] - M[5] * M[6])
       + M[2] * (M[3] * M[7] - M[4] * M[6]);
}

static void hd_elm_get(hd_elm* e, scomplex* sc, INT T, INT fet) {
  sc_fem* fem = sc->fem;
  INT dim = sc->dim, n1 = dim + 1, i, j, k;
  e->dim = dim; e->n1 = n1; e->fet = fet;
  e->nl = (fet == 20) ? dim * n1 / 2 : n1;
  e->nd = (fet == 20 && dim == 3) ? 3 : 1;
  INT p = fem->el_v->IA[T];
  for (i = 0; i < n1; i++) {
    e->v[i] = fem->el_v->JA[p + i];
    for (j = 0; j < dim; j++) e->xv[i * dim + j] = sc->x[e->v[i] * dim + j];
  }
  /* barycentric gradients: rows of the inverse of [x_1-x_0, ..., x_d-x_0] */
  REAL col[3][3];
  for (i = 1; i < n1; i++)
    for (j = 0; j < dim; j++) col[i - 1][j] = e->xv[i * dim + j] - e->xv[j];
  if (dim == 2) {
    REAL det = col[0][0] * col[1][1] - col[0][1] * col[1][0];
    e->g[2] = col[1][1] / det;  e->g[3] = -col[1][0] / det;
    e->g[4] = -col[0][1] / det; e->g[5] = col[0][0] / det;
    e->vol = 0.5 * fabs(det);
  } else {
    REAL r[3];
    hd_cross3(col[1], col[2], r);
    REAL det = hd_dot(col[0], r, 3);
    for (j = 0; j < 3; j++) e->g[3 + j] = r[j] / det;
    hd_cross3(col[2], col[0], r);
    for (j = 0; j < 3; j++) e->g[6 + j] = r[j] / det;
    hd_cross3(col[0], col[1], r);
    for (j = 0; j < 3; j++) e->g[9 + j] = r[j] / det;
    e->vol = fabs(det) / 6.0;
  }
  for (j = 0; j < dim; j++) {
    REAL s = 0.0;
    for (i = 1; i < n1; i++) s += e->g[i * dim + j];
    e->g[j] = -s;
  }
  e->hT = 0.0;
  for (i = 0; i < n1; i++)
    for (k = i + 1; k < n1; k++) {
      REAL d2 = 0.0;
      for (j = 0; j < dim; j++) { REAL d = e->xv[i * dim + j] - e->xv[k * dim + j]; d2 += d * d; }
      if (sqrt(d2) > e->hT) e->hT = sqrt(d2);
    }
  if (fet == 20) {
    p = fem->el_ed->IA[T];
    for (k = 0; k < e->nl; k++) {
      INT E = fem->el_ed->JA[p + k];
      INT q = fem->ed_v->IA[E];
      INT va = fem->ed_v->JA[q], vb = fem->ed_v->JA[q + 1];
      if (va > vb) { INT t = va; va = vb; vb = t; }
      INT la = -1, lb = -1;
      for (i = 0; i < n1; i++) { if (e->v[i] == va) la = i; if (e->v[i] == vb) lb = i; }
      e->dof[k] = E; e->a[k] = la; e->b[k] = lb;
      REAL l2 = 0.0;
      for (j = 0; j < dim; j++) { REAL d = e->xv[lb * dim + j] - e->xv[la * dim + j]; l2 += d * d; }
      e->len[k] = sqrt(l2);
      REAL* ga = e->g + la * dim;
      REAL* gb = e->g + lb * dim;
      if (dim == 2) {
        e->D[k] = 2.0 * e->len[k] * (ga[0] * gb[1] - ga[1] * gb[0]);
      } else {
        hd_cross3(ga, gb, e->D + 3 * k);
        for (j = 0; j < 3; j++) e->D[3 * k + j] *= 2.0 * e->len[k];
      }
    }
  } else {
    REAL fac = (dim == 2) ? 1.0 : 2.0; /* (dim-1)! */
    p = fem->el_f->IA[T];
    for (k = 0; k < n1; k++) {
      INT F = fem->el_f->JA[p + k];
      INT q = fem->f_v->IA[F];
      INT ef[3], on[4] = {0, 0, 0, 0};
      for (i = 0; i < dim; i++) {
        INT vf = fem->f_v->JA[q + i];
        for (j = 0; j < n1; j++) if (e->v[j] == vf) { ef[i] = j; on[j] = 1; }
      }
      INT opp = -1;
      for (j = 0; j < n1; j++) if (!on[j]) opp = j;
      REAL M[9];
      for (INT c = 0; c < dim; c++)
        for (i = 0; i < dim; i++) M[c * dim + i] = e->g[ef[i] * dim + c];
      e->dof[k] = F; e->a[k] = opp; e->b[k] = -1;
      e->len[k] = fac * hd_det(M, dim) * fem->f_area[F];
      e->D[k] = (REAL)dim * e->len[k];
    }
  }
}

/* phi_k at barycentric coordinates lam */
static void hd_phi(const hd_elm* e, INT k, const REAL* lam, REAL* phi) {
  INT dim = e->dim, j;
  if (e->fet == 20) {
    INT a = e->a[k], b = e->b[k];
    for (j = 0; j < dim; j++) phi[j] = e->len[k] * (lam[a] * e->g[b * dim + j] - lam[b] * e->g[a * dim + j]);
  } else {
    INT o = e->a[k];
    for (j = 0; j < dim; j++) {
      REAL x = 0.0;
      for (INT i = 0; i < e->n1; i++) x += lam[i] * e->xv[i * dim + j];
      phi[j] = e->len[k] * (x - e->xv[o * dim + j]);
    }
  }
}

/* local matrix eps (D phi_k, D phi_l) + kappa (phi_k, phi_l), exact */
static void hd_elm_matrix(const hd_elm* e, REAL eps, REAL kap, REAL* Aloc) {
  INT dim = e->dim, n1 = e->n1, nl = e->nl, nd = e->nd, k, l, i, j;
  REAL c0 = e->vol / (REAL)((dim + 1) * (dim + 2));
#define MIJ(i, j) (c0 * (((i) == (j)) ? 2.0 : 1.0))
  for (k = 0; k < nl; k++)
    for (l = 0; l < nl; l++) {
      REAL m = 0.0;
      if (e->fet == 20) {
        const REAL* ga = e->g + e->a[k] * dim; const REAL* gb = e->g + e->b[k] * dim;
        INT c = e->a[l], d = e->b[l];
        const REAL* gc = e->g + c * dim; const REAL* gd = e->g + d * dim;
        m = hd_dot(gb, gd, dim) * MIJ(e->a[k], c) - hd_dot(gb, gc, dim) * MIJ(e->a[k], d)
          - hd_dot(ga, gd, dim) * MIJ(e->b[k], c) + hd_dot(ga, gc, dim) * MIJ(e->b[k], d);
      } else {
        INT ok = e->a[k], ol = e->a[l];
        for (i = 0; i < n1; i++)
          for (j = 0; j < n1; j++) {
            REAL di[3], dj[3];
            for (INT c = 0; c < dim; c++) {
              di[c] = e->xv[i * dim + c] - e->xv[ok * dim + c];
              dj[c] = e->xv[j * dim + c] - e->xv[ol * dim + c];
            }
            m += hd_dot(di, dj, dim) * MIJ(i, j);
          }
      }
      m *= e->len[k] * e->len[l];
      Aloc[k * nl + l] = eps * e->vol * hd_dot(e->D + k * nd, e->D + l * nd, nd) + kap * m;
    }
#undef MIJ
}

/* exact global assembly, row-wise over the dof -> element incidence */
static void hd_assemble(dCSRmat* A, scomplex* sc, INT fet, REAL eps, REAL kap) {
  sc_fem* fem = sc->fem;
  iCSRmat* el_d = (fet == 20) ? fem->el_ed : fem->el_f;
  INT ns = fem->ns_leaf, nd = (fet == 20) ? fem->nedge : fem->nface;
  INT i, k, l, T, E;
  INT* ia = (INT*)calloc(nd + 1, sizeof(INT));
  for (T = 0; T < ns; T++)
    for (k = el_d->IA[T]; k < el_d->IA[T + 1]; k++) ia[el_d->JA[k] + 1]++;
  for (E = 0; E < nd; E++) ia[E + 1] += ia[E];
  INT* ja = (INT*)calloc(ia[nd], sizeof(INT));
  INT* fill = (INT*)calloc(nd, sizeof(INT));
  for (T = 0; T < ns; T++)
    for (k = el_d->IA[T]; k < el_d->IA[T + 1]; k++) { E = el_d->JA[k]; ja[ia[E] + fill[E]++] = T; }
  INT* pos = (INT*)malloc(nd * sizeof(INT));
  for (E = 0; E < nd; E++) pos[E] = -1;
  INT* IA = (INT*)calloc(nd + 1, sizeof(INT));
  INT nnz = 0;
  for (E = 0; E < nd; E++) {
    for (i = ia[E]; i < ia[E + 1]; i++) {
      T = ja[i];
      for (l = el_d->IA[T]; l < el_d->IA[T + 1]; l++) {
        INT F = el_d->JA[l];
        if (pos[F] != E) { pos[F] = E; nnz++; }
      }
    }
    IA[E + 1] = nnz;
  }
  INT* JA = (INT*)calloc(nnz, sizeof(INT));
  REAL* val = (REAL*)calloc(nnz, sizeof(REAL));
  for (E = 0; E < nd; E++) pos[E] = -1;
  REAL Aloc[36];
  hd_elm e;
  for (E = 0; E < nd; E++) {
    INT p = IA[E];
    for (i = ia[E]; i < ia[E + 1]; i++) {
      T = ja[i];
      hd_elm_get(&e, sc, T, fet);
      hd_elm_matrix(&e, eps, kap, Aloc);
      INT kk = -1;
      for (k = 0; k < e.nl; k++) if (e.dof[k] == E) kk = k;
      for (l = 0; l < e.nl; l++) {
        INT F = e.dof[l];
        if (pos[F] < IA[E]) { pos[F] = p; JA[p] = F; p++; }
        val[pos[F]] += Aloc[kk * e.nl + l];
      }
    }
  }
  A->row = nd; A->col = nd; A->nnz = nnz; A->IA = IA; A->JA = JA; A->val = val;
  free(ia); free(ja); free(fill); free(pos);
}

/* homogeneous essential condition: identity rows/columns on boundary dofs */
static void hd_apply_bc(dCSRmat* A, dvector* b, INT* dir) {
  for (INT i = 0; i < A->row; i++) {
    for (INT k = A->IA[i]; k < A->IA[i + 1]; k++) {
      INT j = A->JA[k];
      if (dir[i] || dir[j]) A->val[k] = (i == j) ? 1.0 : 0.0;
    }
    if (dir[i]) b->val[i] = 0.0;
  }
}

/* library coefficients (only for the check of the exact assembly) */
static REAL g_hd_eps = 1.0, g_hd_kap = 1.0;
static void hd_coef_eps(REAL* v, REAL* x, REAL t, void* p) { *v = g_hd_eps; }
static void hd_coef_kap(REAL* v, REAL* x, REAL t, void* p) { *v = g_hd_kap; }

/*********************************************************************/
/* The adaptive loop for FE_type 20 / 30. */
static void run_hcurl_hdiv(scomplex* sc, input_param* inparam, INT fet, INT problem,
                           REAL eps, REAL kap, REAL theta, REAL Eref, INT max_dofs,
                           INT nstages, INT dump_meshes) {
  INT dim = sc->dim, n1 = dim + 1, nq1d = inparam->nquad;
  g_hd_eps = eps; g_hd_kap = kap;

  linear_itsolver_param itparam;
  param_linear_solver_set(&itparam, inparam);
  itparam.linear_precond_type = (fet == 20) ? PREC_HX_CURL_M : PREC_HX_DIV_M;
  INT direct = (itparam.linear_itsolver_type == 0);
  /* The UA AMG setup changes strong_coupled in place: copy amgparam0 at every stage. */
  AMG_param amgparam0, amgparam;
  param_amg_init(&amgparam0);
  param_amg_set(&amgparam0, inparam);
  if (inparam->print_level < 3) { itparam.linear_print_level = 0; amgparam0.print_level = 0; }

  char ftab[MAXFILENAMESIZE];
  snprintf(ftab, sizeof(ftab), "%sconvergence_fe%lld_problem%lld.txt", inparam->output_dir,
           (long long)fet, (long long)problem);
  FILE* fptab = fopen(ftab, "w");
  const char* hdr0 = "%% stage k of the adaptive loop SOLVE -> ESTIMATE -> MARK -> REFINE: mesh after k refinements (stage 0 = initial mesh)\n";
  const char* hdr = "%% stage      ndof         ns  iter      err_V           eta       eta/err   (f,u_h)                 marked\n";
  printf("\n"); printf(hdr0, ""); printf(hdr, "");
  if (fptab) { fprintf(fptab, hdr0, ""); fprintf(fptab, hdr, ""); }

  /*==================================================================*/
  /* ADAPTIVE LOOP:  SOLVE -> ESTIMATE -> MARK -> REFINE               */
  /* stage k works on the mesh obtained after k refinements           */
  /*==================================================================*/
  INT stage;
  for (stage = 0; stage <= nstages; stage++) {
    clock_t t0 = clock();
    scomplex* leaf = (scomplex*)malloc(sizeof(scomplex));
    scfinalize(sc, leaf, (INT)1);
    sc_vols(leaf);
    if (dump_meshes) {
      char fvtu[MAXFILENAMESIZE];
      snprintf(fvtu, sizeof(fvtu), "%smesh_stage_%03lld.vtu", inparam->output_dir, (long long)stage);
      vtu_data vdata;
      vtu_data_init(leaf, &vdata);
      sc_write_vtk(fvtu, &vdata);
      vtu_data_free(&vdata);
    }
    sc_build_fem_data(leaf);
    sc_fem* fem = leaf->fem;
    INT ns = fem->ns_leaf, nf = fem->nface;
    fespace FE;
    create_fespace(&FE, leaf, fet);
    set_dirichlet_bdry(&FE, leaf, 1, 16);
    INT ndof = FE.ndof;

    /*---------------- 1. SOLVE ---------------------------------------*/
    dCSRmat A;
    hd_assemble(&A, leaf, fet, eps, kap);
    if (inparam->print_level > 2 && stage < 2) { /* compare with the library assembly */
      qcoordinates* cq = get_quadrature(leaf, nq1d);
      dCSRmat K = dcsr_create(0, 0, 0), M = dcsr_create(0, 0, 0), L, Dm;
      assemble_global_single(&K, NULL, &FE, leaf, cq, local_assembly_DuDv, NULL, NULL, hd_coef_eps, 0.0);
      assemble_global_single(&M, NULL, &FE, leaf, cq, local_assembly_mass, NULL, NULL, hd_coef_kap, 0.0);
      dcsr_add(&K, 1.0, &M, 1.0, &L);
      dcsr_add(&L, 1.0, &A, -1.0, &Dm);
      REAL dmax = 0.0, amax = 0.0;
      for (INT k = 0; k < Dm.nnz; k++) if (fabs(Dm.val[k]) > dmax) dmax = fabs(Dm.val[k]);
      for (INT k = 0; k < A.nnz; k++) if (fabs(A.val[k]) > amax) amax = fabs(A.val[k]);
      printf("%%   check: max |A_lib - A_exact| / max |A_exact| = %.3e\n", dmax / amax);
      dcsr_free(&K); dcsr_free(&M); dcsr_free(&L); dcsr_free(&Dm);
      free_qcoords(cq); free(cq);
    }
    /* right-hand side by quadrature */
    block_fespace bfe;
    memset(&bfe, 0, sizeof(block_fespace));
    bfe.nspaces = 1; bfe.nun = 1; bfe.ndof = FE.ndof; bfe.nelm = FE.nelm;
    bfe.var_spaces = (fespace**)calloc(1, sizeof(fespace*));
    bfe.var_spaces[0] = &FE;
    simplex_local_data elm;
    fe_local_data fed;
    memset(&elm, 0, sizeof(simplex_local_data));
    memset(&fed, 0, sizeof(fe_local_data));
    initialize_localdata_elm(&elm, &fed, leaf, &bfe, nq1d);
    dvector b = dvec_create(ndof);
    dvec_set(ndof, &b, 0.0);
    hd_elm e;
    REAL fq[3], phi[3], xchk = 0.0;
    INT T, q, k, j, i;
    for (T = 0; T < ns; T++) {
      hd_elm_get(&e, leaf, T, fet);
      get_elmlocaldata(&elm, leaf, T);
      for (q = 0; q < elm.quad_local->nq; q++) {
        REAL* xq = elm.quad_local->x + q * dim;
        REAL* lam = elm.lams + q * n1;
        REAL w = elm.quad_local->w[q];
        for (j = 0; j < dim; j++) { /* lam and xq must describe the same point */
          REAL s = 0.0;
          for (i = 0; i < n1; i++) s += lam[i] * e.xv[i * dim + j];
          if (fabs(s - xq[j]) > xchk) xchk = fabs(s - xq[j]);
        }
        hd_rhs(fet, dim, problem, eps, kap, xq, fq);
        for (k = 0; k < e.nl; k++) { hd_phi(&e, k, lam, phi); b.val[e.dof[k]] += w * hd_dot(fq, phi, dim); }
      }
    }
    if (xchk > 1e-12) fprintf(stderr, "%% WARNING: barycentric/physical quadrature points differ by %.2e\n", xchk);
    hd_apply_bc(&A, &b, FE.dirichlet);
    dvector u = dvec_create(ndof);
    dvec_set(ndof, &u, 0.0);
    amgparam = amgparam0;
    INT iter;
    /* the HX (curl) preconditioner does not converge on very coarse 2D
       meshes (a few hundred dofs); small systems go to the direct solver.
       The HX (div) operators (get_curl_NedtoRT, get_Pigrad_H1toRT) are 3D
       only, so RT0 in 2D also uses the direct solver. */
    if (direct || ndof < 10000 || (fet == 30 && dim == 2)) {
      iter = directsolve_HAZ(&A, &b, &u, 0);
    } else if (fet == 20) {
      dCSRmat P_curl, Grad;
      get_Pigrad_H1toNed(&P_curl, leaf);
      get_grad_H1toNed(&Grad, leaf);
      iter = linear_solver_dcsr_krylov_hx_curl(&A, &b, &u, &itparam, &amgparam, &P_curl, &Grad);
      dcsr_free(&P_curl); dcsr_free(&Grad);
    } else {
      dCSRmat P_curl, Curl, P_div;
      get_Pigrad_H1toNed(&P_curl, leaf);
      get_curl_NedtoRT(&Curl, leaf);
      get_Pigrad_H1toRT(&P_div, &P_curl, &Curl, leaf);
      iter = linear_solver_dcsr_krylov_hx_div(&A, &b, &u, &itparam, &amgparam, &P_curl, &P_div, &Curl);
      dcsr_free(&P_curl); dcsr_free(&Curl); dcsr_free(&P_div);
    }
    if (iter < 0) printf("### ERROR: Solver does not converge with error code = %lld!\n", (long long)iter);
    for (i = 0; i < ndof; i++)
      if (!isfinite(u.val[i])) { fprintf(stderr, "### ERROR: non-finite discrete solution; stopping.\n"); exit(2); }
    REAL fuh = 0.0;
    for (i = 0; i < ndof; i++) fuh += b.val[i] * u.val[i];

    /*---------------- 2. ESTIMATE ------------------------------------*/
    INT nd = e.nd;
    REAL* Dh = (REAL*)calloc(ns * nd, sizeof(REAL));         /* D u_h on T */
    REAL* Uv = (REAL*)calloc(ns * n1 * dim, sizeof(REAL));   /* u_h at the vertices of T */
    REAL* eta2 = (REAL*)calloc(ns, sizeof(REAL));
    INT* fel = (INT*)malloc(2 * nf * sizeof(INT));
    for (i = 0; i < 2 * nf; i++) fel[i] = -1;
    REAL err2 = 0.0, uq[3], uh[3], du[3], r1v[3];
    INT nr1 = (fet == 30 && dim == 3) ? 3 : 1;
    for (T = 0; T < ns; T++) {
      hd_elm_get(&e, leaf, T, fet);
      get_elmlocaldata(&elm, leaf, T);
      REAL c[6];
      for (k = 0; k < e.nl; k++) c[k] = u.val[e.dof[k]];
      for (j = 0; j < nd; j++) {
        REAL s = 0.0;
        for (k = 0; k < e.nl; k++) s += c[k] * e.D[k * nd + j];
        Dh[T * nd + j] = s;
      }
      for (i = 0; i < n1; i++) {
        REAL lam[4] = {0.0, 0.0, 0.0, 0.0};
        lam[i] = 1.0;
        for (j = 0; j < dim; j++) Uv[(T * n1 + i) * dim + j] = 0.0;
        for (k = 0; k < e.nl; k++) {
          hd_phi(&e, k, lam, phi);
          for (j = 0; j < dim; j++) Uv[(T * n1 + i) * dim + j] += c[k] * phi[j];
        }
      }
      REAL r1 = 0.0, r2 = 0.0;
      for (q = 0; q < elm.quad_local->nq; q++) {
        REAL* xq = elm.quad_local->x + q * dim;
        REAL* lam = elm.lams + q * n1;
        REAL w = elm.quad_local->w[q];
        for (j = 0; j < dim; j++) uh[j] = 0.0;
        for (k = 0; k < e.nl; k++) {
          hd_phi(&e, k, lam, phi);
          for (j = 0; j < dim; j++) uh[j] += c[k] * phi[j];
        }
        hd_rhs(fet, dim, problem, eps, kap, xq, fq);
        REAL rr = 0.0;
        for (j = 0; j < dim; j++) rr += (fq[j] - kap * uh[j]) * (fq[j] - kap * uh[j]);
        r2 += w * rr;
        hd_R1(fet, dim, problem, eps, kap, xq, r1v);
        r1 += w * hd_dot(r1v, r1v, nr1);
        if (problem == 3) { /* error against the exact solution (verification) */
          hd_exact_u(fet, dim, xq, uq);
          hd_exact_du(fet, dim, xq, du);
          REAL ec = 0.0, eu = 0.0;
          for (j = 0; j < nd; j++) ec += (du[j] - Dh[T * nd + j]) * (du[j] - Dh[T * nd + j]);
          for (j = 0; j < dim; j++) eu += (uq[j] - uh[j]) * (uq[j] - uh[j]);
          err2 += w * (eps * ec + kap * eu);
        }
      }
      REAL hb = fmin(e.hT / sqrt(eps), 1.0 / sqrt(kap));
      eta2[T] = e.hT * e.hT * r1 / kap + hb * hb * r2;
      for (k = fem->el_f->IA[T]; k < fem->el_f->IA[T + 1]; k++) {
        INT f = fem->el_f->JA[k];
        if (fel[2 * f] < 0) fel[2 * f] = T; else fel[2 * f + 1] = T;
      }
    }
    /* face residuals on interior faces */
    for (INT f = 0; f < nf; f++) {
      INT T1 = fel[2 * f], T2 = fel[2 * f + 1];
      if (T1 < 0 || T2 < 0) continue;
      REAL* n = fem->f_norm + f * dim;
      REAL area = fem->f_area[f];
      INT p = fem->f_v->IA[f];
      INT fv[3];
      for (i = 0; i < dim; i++) fv[i] = fem->f_v->JA[p + i];
      REAL hS = 0.0;
      for (i = 0; i < dim; i++)
        for (INT i2 = i + 1; i2 < dim; i2++) {
          REAL d2 = 0.0;
          for (j = 0; j < dim; j++) { REAL d = leaf->x[fv[i] * dim + j] - leaf->x[fv[i2] * dim + j]; d2 += d * d; }
          if (sqrt(d2) > hS) hS = sqrt(d2);
        }
      /* J1 = kappa [u_h].n (curl) or kappa [u_h] x n (div): affine on S,
         int_S |w|^2 = |S| / (d (d+1)) (sum_v |w_v|^2 + |sum_v w_v|^2) */
      REAL s1 = 0.0, s2[3] = {0.0, 0.0, 0.0};
      INT p1 = fem->el_v->IA[T1], p2 = fem->el_v->IA[T2];
      for (i = 0; i < dim; i++) {
        INT l1 = -1, l2 = -1;
        for (j = 0; j < n1; j++) {
          if (fem->el_v->JA[p1 + j] == fv[i]) l1 = j;
          if (fem->el_v->JA[p2 + j] == fv[i]) l2 = j;
        }
        REAL jv[3], w[3] = {0.0, 0.0, 0.0};
        for (j = 0; j < dim; j++) jv[j] = kap * (Uv[(T1 * n1 + l1) * dim + j] - Uv[(T2 * n1 + l2) * dim + j]);
        if (fet == 20) w[0] = hd_dot(jv, n, dim);
        else if (dim == 2) w[0] = -jv[0] * n[1] + jv[1] * n[0];
        else hd_cross3(jv, n, w);
        s1 += hd_dot(w, w, 3);
        for (j = 0; j < 3; j++) s2[j] += w[j];
      }
      REAL J1 = area / (REAL)(dim * (dim + 1)) * (s1 + hd_dot(s2, s2, 3));
      /* J2 = [eps curl u_h] x n (curl) or [eps div u_h] (div), constant on S */
      REAL J2;
      if (nd == 3) {
        REAL w[3], r[3];
        for (j = 0; j < 3; j++) w[j] = eps * (Dh[T1 * 3 + j] - Dh[T2 * 3 + j]);
        hd_cross3(w, n, r);
        J2 = hd_dot(r, r, 3) * area;
      } else {
        REAL jc = eps * (Dh[T1] - Dh[T2]);
        J2 = jc * jc * area;
      }
      REAL hbS = fmin(hS / sqrt(eps), 1.0 / sqrt(kap));
      REAL a = hS * J1 / kap + hbS * J2 / sqrt(eps);
      eta2[T1] += a; eta2[T2] += a;
    }
    REAL eta = 0.0;
    for (T = 0; T < ns; T++) eta += eta2[T];
    eta = sqrt(eta);
    REAL err = -1.0;
    if (problem == 3) err = sqrt(err2);
    else if (Eref > 0.0) err = (Eref - fuh > 0.0) ? sqrt(Eref - fuh) : -sqrt(fuh - Eref);

    /*---------------- 3. MARK (Doerfler) -----------------------------*/
    INT last = (stage == nstages) || (ndof > max_dofs);
    INT nmarked = 0;
    ivector marked;
    if (!last) marked = dorfler_mark(ns, eta2, theta, &nmarked);
    const char* fmt = "%7lld %10lld %10lld %5lld  %.6e  %.6e  %7.4f   %.15e %8lld\n";
    printf(fmt, (long long)stage, (long long)(ndof - FE.nbdof), (long long)ns, (long long)iter,
           err, eta, (err > 0.0) ? eta / err : 0.0, fuh, (long long)nmarked);
    if (fptab) {
      fprintf(fptab, fmt, (long long)stage, (long long)(ndof - FE.nbdof), (long long)ns, (long long)iter,
              err, eta, (err > 0.0) ? eta / err : 0.0, fuh, (long long)nmarked);
      fflush(fptab);
    }
    if (inparam->print_level > 1)
      printf("%%   CPU seconds for this stage: %.2f\n", (REAL)(clock() - t0) / CLOCKS_PER_SEC);
    fflush(stdout);

    if (!last) {
      /*-------------- 4. REFINE (DGS bisection) -----------------------*/
      /* scfinest numbers the leaves (same order as scfinalize) for refine() */
      scomplex* sctop = scfinest(sc);
      haz_scomplex_free(sctop);
      refine(1, sc, &marked);
      ivec_free(&marked);
    } else {
      printf("\n table: %s\n", ftab);
    }
    free(Dh); free(Uv); free(eta2); free(fel);
    dvec_free(&u); dvec_free(&b); dcsr_free(&A);
    free_simplexlocaldata(&elm); free_felocaldata(&fed); free(bfe.var_spaces);
    free_fespace(&FE);
    haz_scomplex_free(leaf);
    if (last) break;
  }
  if (fptab) fclose(fptab);
}
