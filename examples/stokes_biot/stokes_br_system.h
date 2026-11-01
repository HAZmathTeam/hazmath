/*! \file examples/stokes_biot/stokes_br_system.h
 *
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief Local assembly and error computation for the Bernardi-Raugel
 *        (P1 + face bubbles) / P0 discretization of Stokes:
 *
 *          a(u,v) - (p, div v) - (div u, q) = (f,v),
 *          a(u,v) = 2 mu (eps(u), eps(v)),  eps(u) = (grad u + grad u^T)/2.
 *
 *        Static condensation of the bubbles (default): the bubble-bubble
 *        block A_bb is replaced by its diagonal D and the bubbles are
 *        eliminated, which gives a P1/P0 system with the pressure block
 *        -B_b D^{-1} B_b^T; the bubbles are recovered after the solve.
 *
 *        Block layout: [ bubble (FEtype 61, vector), u_1, ..., u_d (FEtype 1),
 *                        p (FEtype 0) ].
 *
 * \note Basis layout in fe_local_data:
 *         P1:      dphi[t*dim + j]           = d_j phi_t
 *         bubble:  phi[t*dim + i]            = (phi_t)_i
 *                  dphi[t*dim*dim + j*dim+i] = d_j (phi_t)_i
 *
 *        Pressure-robust right hand side (default, Li and Zikatanov,
 *        CMAME 393 (2022) 114815): the source is tested with Pi_h v,
 *        Pi_h the BDM1 interpolant.  Pi_h v = v on P1 and for a face bubble
 *        b_F = c lam..lam n_F (c = 4 in 2D, 8 in 3D), Pi_h b_F = c_F psi_F with
 *        psi_F the RT0 function with int_F psi_F.n_F = 1 and
 *        c_F = int_F b_F.n_F = 2|F|/3 (2D), 2|F|/15 (3D).
 *
 * \note Adapted from Stokes_elements/common/stokes_br_system.h and
 *       stokes_br_error.h.
 *
 */

// 1: test the source with Pi_h v (pressure robust); 0: with v
INT stokes_br_pr = 1;
// 1: bubbles condensed, only the diagonal of the bubble-bubble block is assembled
INT stokes_br_condense = 1;

/* Reference-element tables, built once (br_ref_init) and used on every
 * element.  Face t of T is identified by the local vertex o opposite to it;
 * its bubble is b_t = 2^d prod_{j != o} lam_j n_t and
 *   grad(prod_{j != o} lam_j) = sum_{k != o} (prod_{j != o,k} lam_j) grad lam_k.
 * Matrix rule: d points per direction (exact to degree 2d-1), weights sum to 1/d!:
 *   mw[q]          reference weights,
 *   mgr[q][o][k]   prod_{j != o,k} lam_j at point q (0 for k == o).
 * Face rule for the traction (Test 2): nf1 points on [-1,1] (2D) or the
 *   nf1 x nf1 rule on the reference triangle (3D).
 * The rhs uses the nquad rule of elm_data, whose barycentric coordinates
 * elm_data->lams are already reference values. */
typedef struct {
  INT dim;
  INT nq;
  REAL mw[27];
  REAL mgr[27][4][4];
  INT nf1;
  INT nfq;
  REAL fgp[2 * 49];
  REAL fgw[49];
} br_ref_tables;
static br_ref_tables br_ref = {0, 0};

static void br_ref_init(INT dim, INT nf1)
{
  INT q, o, k, j;
  if (br_ref.dim != dim) {
    INT nq1 = dim;
    INT nq = (dim == 2) ? nq1 * nq1 : nq1 * nq1 * nq1;
    REAL gp[3 * 27];
    if (dim == 2) triquad_(gp, br_ref.mw, nq1);
    else tetquad_(gp, br_ref.mw, nq1);
    for (q = 0; q < nq; q++) {
      REAL lam[4];
      lam[0] = 1.0;
      for (j = 0; j < dim; j++) {
        lam[j + 1] = gp[j * nq + q];
        lam[0] -= lam[j + 1];
      }
      for (o = 0; o < dim + 1; o++)
        for (k = 0; k < dim + 1; k++) {
          REAL p = 0.0;
          if (k != o) {
            p = 1.0;
            for (j = 0; j < dim + 1; j++)
              if (j != o && j != k) p *= lam[j];
          }
          br_ref.mgr[q][o][k] = p;
        }
    }
    br_ref.nq = nq;
    br_ref.dim = dim;
    br_ref.nf1 = 0;
  }
  if (nf1 > ((dim == 2) ? 5 : 7)) nf1 = (dim == 2) ? 5 : 7;
  if (br_ref.nf1 != nf1) {
    if (dim == 2) quad1d(br_ref.fgp, br_ref.fgw, nf1);
    else triquad_(br_ref.fgp, br_ref.fgw, nf1);
    br_ref.nfq = (dim == 2) ? nf1 : nf1 * nf1;
    br_ref.nf1 = nf1;
  }
}

/* Face data of T: ef[t][k] local index of the k-th vertex of face t,
 * opp[t] the local vertex opposite to face t, and
 * sF[t] = sign(n_t . (x_f - x_opp)), x_f a vertex of face t. */
static void br_face_map(simplex_local_data* elm_data, INT ef[4][3], INT* opp, REAL* sF)
{
  INT dim = elm_data->dim;
  INT t, k, j;
  for (t = 0; t < dim + 1; t++) {
    opp[t] = -1;
    for (k = 0; k < dim; k++)
      for (j = 0; j < dim + 1; j++)
        if (elm_data->v_on_f[t * dim + k] == elm_data->local_v[j]) ef[t][k] = j;
    for (j = 0; j < dim + 1 && opp[t] < 0; j++) {
      INT onf = 0;
      for (k = 0; k < dim; k++)
        if (elm_data->v_on_f[t * dim + k] == elm_data->local_v[j]) onf = 1;
      if (!onf) opp[t] = j;
    }
    INT jf = (opp[t] == 0) ? 1 : 0;
    REAL sgn = 0.0;
    for (k = 0; k < dim; k++)
      sgn += (elm_data->xv[jf * dim + k] - elm_data->xv[opp[t] * dim + k]) * elm_data->f_norm[t * dim + k];
    sF[t] = (sgn > 0.0) ? 1.0 : -1.0;
  }
}

/*!
 * \fn static void local_rhs_br(REAL* bLoc, simplex_local_data* elm_data, fe_local_data* fe_data, void (*rhs)(REAL*, REAL*, REAL, void*), REAL time)
 *
 * \brief Local rhs (f, v) + <g_N, v>_{Gamma_N} (or with Pi_h v if
 *        stokes_br_pr) for the velocity test functions and -(g, q) for the
 *        pressure.
 *        Pi_h b_F = c_F psi_F, psi_F|_T = s_F (x - x_opp) / (d |T|), with
 *        s_F = n_F.n_out = +-1 and x_opp the vertex of T opposite to F.
 *        The P1 basis functions are the barycentric coordinates
 *        elm_data->lams (reference values) and the pressure is P0.
 */
static void local_rhs_br(REAL* bLoc, simplex_local_data* elm_data, fe_local_data* fe_data,
                         void (*rhs)(REAL*, REAL*, REAL, void*), REAL time)
{
  INT dim = elm_data->dim;
  INT nq = elm_data->quad_local->nq;
  INT dpeB = fe_data->n_dof_per_space[0];
  INT dpeP = fe_data->n_dof_per_space[dim + 1];
  INT nvel = fe_data->n_dof - dpeP;
  INT q, k, t, j, l;
  REAL f[7];

  REAL sF[4];
  INT ef[4][3], opp[4];
  br_face_map(elm_data, ef, opp, sF);
  REAL cF = (dim == 2) ? 2.0 / 3.0 : 2.0 / 15.0;
  REAL cb = (dim == 2) ? 4.0 : 8.0;
  // Pi_h b_t = scl[t] (x - x_opp[t]) (rhs=pr), or b_t = cb prod lam n_t (rhs=std)
  REAL scl[4];
  for (t = 0; t < dpeB; t++) scl[t] = cF * elm_data->f_area[t] * sF[t] / (dim * elm_data->vol);

  for (q = 0; q < nq; q++) {
    REAL* x = elm_data->quad_local->x + q * dim;
    REAL w = elm_data->quad_local->w[q];
    REAL* lam = elm_data->lams + q * (dim + 1);
    (*rhs)(f, x, time, NULL);

    // bubbles
    for (t = 0; t < dpeB; t++) {
      REAL ft = 0.0;
      REAL* xo = elm_data->xv + opp[t] * dim;
      if (stokes_br_pr) {
        for (k = 0; k < dim; k++) ft += f[k] * (x[k] - xo[k]);
        ft *= scl[t];
      } else {
        REAL b = cb;
        for (j = 0; j < dim + 1; j++)
          if (j != opp[t]) b *= lam[j];
        for (k = 0; k < dim; k++) ft += f[k] * elm_data->f_norm[t * dim + k];
        ft *= b;
      }
      bLoc[t] += w * ft;
    }
    // P1 components (local P1 dofs = local vertices)
    l = dpeB;
    for (k = 0; k < dim; k++)
      for (t = 0; t < dim + 1; t++, l++) bLoc[l] += w * f[dim + k] * lam[t];
    // pressure (P0): -(g, q)
    for (t = 0; t < dpeP; t++) bLoc[nvel + t] += w * f[2 * dim];
  }

  if (stokes_test != 2) return;

  /* Traction on the faces F of T in Gamma_N: <g_N, v>_F, or, if stokes_br_pr,
   *   <g_N.n, (Pi_h v).n>_F + <g_N - (g_N.n) n, v>_F,
   * i.e. only the normal component, which carries the pressure, is tested
   * with Pi_h v.  Pi_h v = v for P1; for bubbles (Pi_h b_F').n = 0 on F if
   * F' != F, (Pi_h b_F).n_F = c_F/|F| = cF on F and b_F has no tangential part.
   * Face quadrature with the nquad of input.dat: nquad-point Gauss on [-1,1]
   * (2D, at most 5) or the nquad x nquad rule on the reference triangle (3D,
   * at most 7).  With rhs=pr the pressure part of (f, Pi_h v) cancels the
   * term <p n, Pi_h v>_{Gamma_N} only up to quadrature error, which enters
   * the velocity divided by mu; a face rule less accurate than the volume
   * rule leaves that error at small mu. */
  br_ref_init(dim, elm_data->quad_local->nq1d);
  INT nfq = br_ref.nfq;
  REAL* gp = br_ref.fgp;
  REAL* gw = br_ref.fgw;
  REAL lam[4], x[3], nout[3], gN[3];

  for (INT F = 0; F < dpeB; F++) {
    REAL xm[3] = {0.0, 0.0, 0.0};
    for (k = 0; k < dim; k++)
      for (j = 0; j < dim; j++) xm[j] += elm_data->xv[ef[F][k] * dim + j] / dim;
    if (!on_gamma_N(xm, dim)) continue;
    for (j = 0; j < dim; j++) nout[j] = sF[F] * elm_data->f_norm[F * dim + j];

    for (INT fq = 0; fq < nfq; fq++) {
      // barycentric coordinates on F; reference weights sum to 2 (2D), 1/2 (3D)
      REAL wf = elm_data->f_area[F] * ((dim == 2) ? 0.5 * gw[fq] : 2.0 * gw[fq]);
      for (j = 0; j < dim + 1; j++) lam[j] = 0.0;
      if (dim == 2) {
        lam[ef[F][0]] = 0.5 * (1.0 - gp[fq]);
        lam[ef[F][1]] = 0.5 * (1.0 + gp[fq]);
      } else {
        lam[ef[F][0]] = 1.0 - gp[fq] - gp[nfq + fq];
        lam[ef[F][1]] = gp[fq];
        lam[ef[F][2]] = gp[nfq + fq];
      }
      for (j = 0; j < dim; j++) {
        x[j] = 0.0;
        for (k = 0; k < dim + 1; k++) x[j] += lam[k] * elm_data->xv[k * dim + j];
      }
      traction(gN, x, nout, dim);

      // bubbles: only b_F is nonzero on F and it is normal to F
      REAL gn = 0.0;
      for (k = 0; k < dim; k++) gn += gN[k] * elm_data->f_norm[F * dim + k];
      REAL bF = cF;
      if (!stokes_br_pr) {
        bF = (dim == 2) ? 4.0 : 8.0;
        for (k = 0; k < dim; k++) bF *= lam[ef[F][k]];
      }
      bLoc[F] += wf * gn * bF;
      // P1 components (local P1 dofs = local vertices)
      l = dpeB;
      for (k = 0; k < dim; k++)
        for (t = 0; t < dim + 1; t++, l++) bLoc[l] += wf * gN[k] * lam[t];
    }
  }
}

/*!
 * \fn void local_assembly_stokes_br(REAL* ALoc, REAL* bLoc, REAL* u_local, simplex_local_data* elm_data, fe_local_data* fe_data, void (*rhs)(REAL*, REAL*, REAL, void*), void (*coeff)(REAL*, REAL*, REAL, void*), REAL time)
 *
 * \brief Local matrix for BR/P0 Stokes (and local rhs if bLoc != NULL).
 *        Velocity is written as w = sum_t c_t b_t + sum_k sum_s u_{k,s} phi_s e_k;
 *        for each quadrature point we form, for every velocity basis function,
 *        its symmetric gradient E[n][i*dim+j] = eps_ij(w), then
 *        a = 2 mu sum_ij E_ij E_ij (symmetric, upper triangle mirrored).
 *        The integrands are polynomials of degree <= 2d-2, so the matrix uses
 *        the reference rule with d points per direction (exact to degree 2d-1)
 *        instead of the nquad rule of elm_data, which is kept for the rhs.
 *        With stokes_br_condense only the diagonal of the bubble-bubble block
 *        is computed: condense_bubbles uses D = diag(A_bb) and nothing else.
 */
void local_assembly_stokes_br(REAL* ALoc, REAL* bLoc, REAL* u_local,
                              simplex_local_data* elm_data, fe_local_data* fe_data,
                              void (*rhs)(REAL*, REAL*, REAL, void*),
                              void (*coeff)(REAL*, REAL*, REAL, void*),
                              REAL time)
{
  INT dim = elm_data->dim;
  INT n = fe_data->n_dof;
  INT iP = dim + 1;
  INT dpeB = fe_data->n_dof_per_space[0];
  INT dpeP = fe_data->n_dof_per_space[iP];
  INT nvel = n - dpeP;  // velocity dofs come first in the local ordering
  INT dd = dim * dim;
  INT q, k, t, T, i, j;

  if (ALoc != NULL) {
    br_ref_init(dim, elm_data->quad_local->nq1d);
    INT nq = br_ref.nq;
    REAL wscl = elm_data->vol * ((dim == 2) ? 2.0 : 6.0);
    REAL* dlam = elm_data->gradlams;  // constant on T
    REAL cb = (dim == 2) ? 4.0 : 8.0;
    REAL tmu = 2.0 * stokes_mu;
    REAL sF[4];
    INT ef[4][3], opp[4];
    br_face_map(elm_data, ef, opp, sF);

    // E[t][i*dim+j] = eps_ij(phi_t), div[t] = div(phi_t); P1 rows are constant
    // on T, bubble rows are integrated: Eb = int eps(b_t), divb = int div(b_t)
    REAL E[16 * 9], dv[16], Eq[4 * 9], dvq[4], gb[4][3];
    INT row = dpeB;
    for (k = 0; k < dim; k++)
      for (t = 0; t < dim + 1; t++, row++) {
        // phi = lam_t e_k: G_ij = delta_ik d_j lam_t
        for (i = 0; i < dd; i++) E[row * dd + i] = 0.0;
        for (j = 0; j < dim; j++) {
          E[row * dd + k * dim + j] += 0.5 * dlam[t * dim + j];
          E[row * dd + j * dim + k] += 0.5 * dlam[t * dim + j];
        }
        dv[row] = dlam[t * dim + k];
      }
    for (t = 0; t < dpeB; t++) {
      dv[t] = 0.0;
      for (i = 0; i < dd; i++) E[t * dd + i] = 0.0;
    }

    for (q = 0; q < nq; q++) {
      REAL w = wscl * br_ref.mw[q];
      // bubbles: b_t = cb prod_{j != o} lam_j n_t, G_ij = n_i d_j(cb prod)
      for (t = 0; t < dpeB; t++) {
        REAL* nt = elm_data->f_norm + t * dim;
        for (j = 0; j < dim; j++) {
          REAL g = 0.0;
          for (k = 0; k < dim + 1; k++) g += br_ref.mgr[q][opp[t]][k] * dlam[k * dim + j];
          gb[t][j] = cb * g;
        }
        dvq[t] = 0.0;
        for (i = 0; i < dim; i++) {
          dvq[t] += nt[i] * gb[t][i];
          for (j = 0; j < dim; j++) Eq[t * dd + i * dim + j] = 0.5 * (nt[i] * gb[t][j] + nt[j] * gb[t][i]);
        }
        dv[t] += w * dvq[t];
        for (i = 0; i < dd; i++) E[t * dd + i] += w * Eq[t * dd + i];
      }
      // bubble-bubble: 2 mu (eps(b_T), eps(b_t)), upper triangle
      for (t = 0; t < dpeB; t++)
        for (T = t; T < dpeB; T++) {
          if (stokes_br_condense && T != t) continue;
          REAL aij = 0.0;
          for (i = 0; i < dd; i++) aij += Eq[t * dd + i] * Eq[T * dd + i];
          ALoc[t * n + T] += w * tmu * aij;
        }
    }

    // bubble-P1 (int eps(b_t) : eps(phi_T)) and P1-P1 (|T| eps : eps)
    for (t = 0; t < nvel; t++)
      for (T = (t < dpeB) ? dpeB : t; T < nvel; T++) {
        REAL aij = 0.0;
        for (i = 0; i < dd; i++) aij += E[t * dd + i] * E[T * dd + i];
        ALoc[t * n + T] += tmu * ((t < dpeB) ? aij : elm_data->vol * aij);
      }
    // velocity-pressure (P0): -(p, div v) and -(div u, q)
    for (t = 0; t < nvel; t++) {
      REAL bij = -((t < dpeB) ? dv[t] : elm_data->vol * dv[t]);
      for (T = 0; T < dpeP; T++) {
        ALoc[t * n + nvel + T] += bij;
        ALoc[(nvel + T) * n + t] += bij;
      }
    }
    // lower triangle of the velocity block
    for (t = 0; t < nvel; t++)
      for (T = 0; T < t; T++) ALoc[t * n + T] = ALoc[T * n + t];
  }
  if (bLoc != NULL && rhs != NULL)
    local_rhs_br(bLoc, elm_data, fe_data, rhs, time);
  return;
}

/*!
 * \fn void meanzero_pressure_br(block_dCSRmat* A, dvector* b, scomplex* sc)
 *
 * \brief Append the constraint int p = 0 (P0 pressure, last block) as an
 *        extra row and column of A, and an extra zero entry of b.
 */
void meanzero_pressure_br(block_dCSRmat* A, dvector* b, scomplex* sc)
{
  b->val = (REAL*)realloc(b->val, (b->row + 1) * sizeof(REAL));
  b->val[b->row] = 0.0;
  b->row++;
  bdcsr_extend(A, sc->fem->el_vol, sc->fem->el_vol, A->brow - 1, 1.0, 1.0);
}

/*!
 * \fn void condense_bubbles(block_dCSRmat* A, dvector* b, block_dCSRmat* Ac, dvector* bc, REAL* Dinv)
 *
 * \brief Static condensation of the bubbles (block 0 of A, b).  With
 *        D = diag(A_00) and i,j >= 1:
 *          Ac_ij = A_ij - A_i0 D^{-1} A_0j,   bc_i = b_i - A_i0 D^{-1} b_0.
 *        On output Dinv = D^{-1} (length of block 0).
 */
void condense_bubbles(block_dCSRmat* A, dvector* b, block_dCSRmat* Ac, dvector* bc, REAL* Dinv)
{
  INT ns = A->brow, nc = ns - 1;
  INT i, j, k, nb = A->blocks[0]->row;

  dvector diag;
  dcsr_getdiag(nb, A->blocks[0], &diag);
  for (k = 0; k < nb; k++) Dinv[k] = 1.0 / diag.val[k];
  dvec_free(&diag);

  bdcsr_alloc(nc, nc, Ac);
  for (i = 1; i < ns; i++) {
    dCSRmat* Ai0 = A->blocks[i * ns];
    // S = A_i0 D^{-1}
    dCSRmat S = dcsr_create(0, 0, 0);
    if (Ai0 != NULL) {
      dcsr_alloc(Ai0->row, Ai0->col, Ai0->nnz, &S);
      dcsr_cp(Ai0, &S);
      for (k = 0; k < S.nnz; k++) S.val[k] *= Dinv[S.JA[k]];
    }
    for (j = 1; j < ns; j++) {
      dCSRmat* Aij = A->blocks[i * ns + j];
      dCSRmat* A0j = A->blocks[j];
      dCSRmat* C = Ac->blocks[(i - 1) * nc + (j - 1)];
      if (Ai0 != NULL && A0j != NULL) {
        dCSRmat T = dcsr_create(0, 0, 0);
        dcsr_mxm(&S, A0j, &T);
        dcsr_add(Aij, 1.0, &T, -1.0, C);
        dcsr_free(&T);
      } else if (Aij != NULL) {
        dcsr_alloc(Aij->row, Aij->col, Aij->nnz, C);
        dcsr_cp(Aij, C);
      } else {
        free(C);
        Ac->blocks[(i - 1) * nc + (j - 1)] = NULL;
      }
    }
    if (Ai0 != NULL) dcsr_free(&S);
  }

  // Right hand side
  // block sizes from the bubble row: A_0j has as many columns as block j
  INT* off = (INT*)calloc(ns + 1, sizeof(INT));
  for (i = 0; i < ns; i++) off[i + 1] = off[i] + A->blocks[i]->col;
  dvector t = dvec_create(nb);
  for (k = 0; k < nb; k++) t.val[k] = Dinv[k] * b->val[k];
  dvector bnew = dvec_create(b->row - nb);
  for (k = 0; k < bnew.row; k++) bnew.val[k] = b->val[nb + k];
  for (i = 1; i < ns; i++)
    if (A->blocks[i * ns] != NULL)
      dcsr_aAxpy(-1.0, A->blocks[i * ns], t.val, bnew.val + off[i] - nb);
  *bc = bnew;
  dvec_free(&t);
  free(off);
}

/*!
 * \fn void recover_bubbles(block_dCSRmat* A, dvector* b, REAL* Dinv, REAL* xc, REAL* x)
 *
 * \brief Full solution x = (x_b, xc) with x_b = D^{-1} (b_0 - sum_j A_0j xc_j).
 */
void recover_bubbles(block_dCSRmat* A, dvector* b, REAL* Dinv, REAL* xc, REAL* x)
{
  INT ns = A->brow;
  INT j, k, nb = A->blocks[0]->row;
  INT n = b->row;
  for (k = 0; k < n - nb; k++) x[nb + k] = xc[k];
  for (k = 0; k < nb; k++) x[k] = b->val[k];
  INT offj = nb;
  for (j = 1; j < ns; j++) {
    dcsr_aAxpy(-1.0, A->blocks[j], x + offj, x);
    offj += A->blocks[j]->col;
  }
  for (k = 0; k < nb; k++) x[k] *= Dinv[k];
}

/*!
 * \fn void stokes_br_errors(REAL* err, REAL* sol, block_fespace* FE, scomplex* sc, qcoordinates* cq)
 *
 * \brief Errors of the BR/P0 solution (bubble + P1 summed into one velocity):
 *        err[0] = ||u - u_h||_L2,  err[1] = |u - u_h|_H1,
 *        err[2] = min_c ||p - p_h - c||_L2 (pure Dirichlet), ||p - p_h||_L2 (traction),
 *        err[3] = max_T |int_T div(u - u_h)|.
 */
void stokes_br_errors(REAL* err, REAL* sol, block_fespace* FE, scomplex* sc,
                      qcoordinates* cq)
{
  INT dim = sc->dim;
  INT iP = dim + 1;
  INT elm, q, s, k, t, j;
  REAL eu = 0.0, egu = 0.0, ep = 0.0, pint = 0.0, vol = 0.0, maxdiv = 0.0;
  REAL uh[3], guh[9], ue[3], gue[9];

  INT* off = (INT*)calloc(FE->nspaces + 1, sizeof(INT));
  for (s = 0; s < FE->nspaces; s++) off[s + 1] = off[s] + FE->var_spaces[s]->ndof;

  simplex_local_data elm_data;
  fe_local_data fe_data;
  memset(&elm_data, 0, sizeof(simplex_local_data));
  memset(&fe_data, 0, sizeof(fe_local_data));
  initialize_localdata_elm(&elm_data, &fe_data, sc, FE, cq->nq1d);

  for (elm = 0; elm < sc->fem->ns_leaf; elm++) {
    get_elmlocaldata(&elm_data, sc, elm);
    get_felocaldata_elm(&fe_data, FE, NULL, elm);
    INT dofoff = 0;
    for (s = 0; s < FE->nspaces; s++) {
      for (k = 0; k < fe_data.n_dof_per_space[s]; k++)
        fe_data.u_local[dofoff + k] = sol[fe_data.local_dof[dofoff + k] + off[s]];
      dofoff += fe_data.n_dof_per_space[s];
    }

    REAL intdiv = 0.0;
    for (q = 0; q < elm_data.quad_local->nq; q++) {
      REAL* x = elm_data.quad_local->x + q * dim;
      REAL w = elm_data.quad_local->w[q];
      for (s = 0; s < FE->nspaces; s++)
        get_FEM_basis_at_quadpt(&elm_data, &fe_data, s, q);

      for (k = 0; k < dim * dim; k++) guh[k] = 0.0;
      for (k = 0; k < dim; k++) uh[k] = 0.0;
      // bubbles
      INT dpeB = fe_data.n_dof_per_space[0];
      for (t = 0; t < dpeB; t++) {
        REAL c = fe_data.u_local[t];
        for (k = 0; k < dim; k++) {
          uh[k] += c * fe_data.phi[0][t * dim + k];
          for (j = 0; j < dim; j++)
            guh[k * dim + j] += c * fe_data.dphi[0][t * dim * dim + j * dim + k];
        }
      }
      // P1 components
      dofoff = dpeB;
      for (k = 0; k < dim; k++) {
        for (t = 0; t < fe_data.n_dof_per_space[1 + k]; t++) {
          REAL c = fe_data.u_local[dofoff + t];
          uh[k] += c * fe_data.phi[1 + k][t];
          for (j = 0; j < dim; j++) guh[k * dim + j] += c * fe_data.dphi[1 + k][t * dim + j];
        }
        dofoff += fe_data.n_dof_per_space[1 + k];
      }
      // pressure
      REAL ph = 0.0;
      for (t = 0; t < fe_data.n_dof_per_space[iP]; t++)
        ph += fe_data.u_local[dofoff + t] * fe_data.phi[iP][t];

      exact_u(ue, x, dim);
      Dexact_u(gue, x, dim);
      REAL divh = 0.0;
      for (k = 0; k < dim; k++) {
        eu += w * (uh[k] - ue[k]) * (uh[k] - ue[k]);
        divh += guh[k * dim + k];
      }
      for (k = 0; k < dim * dim; k++) egu += w * (guh[k] - gue[k]) * (guh[k] - gue[k]);
      REAL dp = ph - exact_p(x, dim);
      ep += w * dp * dp;
      pint += w * dp;
      vol += w;
      intdiv += w * (divh - exact_div(x, dim));
    }
    if (fabs(intdiv) > maxdiv) maxdiv = fabs(intdiv);
  }

  free_simplexlocaldata(&elm_data);
  free_felocaldata(&fe_data);
  free(off);

  if (stokes_test != 2) ep -= pint * pint / vol;
  err[0] = sqrt(eu);
  err[1] = sqrt(egu);
  err[2] = sqrt(ep > 0.0 ? ep : 0.0);
  err[3] = maxdiv;
}
