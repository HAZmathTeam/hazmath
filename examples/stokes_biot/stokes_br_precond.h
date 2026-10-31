/*! \file examples/stokes_biot/stokes_br_precond.h
 *
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief Block lower triangular preconditioner for the BR/P0 Stokes system
 *        with blocks [ v_0, ..., v_{m-1}, p ] (velocity blocks first,
 *        P0 pressure last):
 *
 *          P = [ A_00                     ]
 *              [ A_10  A_11               ]
 *              [ ...         ...          ]
 *              [ B_0   B_1  ...   -S      ],   S = M_p/mu - A_pp,
 *
 *        where M_p is the (diagonal) P0 mass matrix and A_pp is the pressure
 *        block of the system (zero for full BR, -B_b D^{-1} B_b^T after
 *        condensation of the bubbles).  Every diagonal block is inverted
 *        approximately by one AMG cycle with forward/backward Gauss-Seidel
 *        smoothing, so P^{-1} is a fixed linear operator:
 *          BR_AMG_RS: classical Ruge-Stuben AMG (src/amg_classical), V-cycle;
 *          BR_AMG_UA: unsmoothed aggregation AMG, AMLI-, W- or V-cycle
 *                     (the AMLI-cycle is the polynomial, i.e. linear, one).
 */

#define BR_AMG_RS 0
#define BR_AMG_UA 1

typedef struct {
  INT nb;            // number of blocks
  block_dCSRmat* A;  // the system
  INT* off;          // block offsets in the global vector
  AMG_data** mgl;    // AMG hierarchies for the diagonal blocks (pressure: S)
  AMG_param amgparam;
  INT amg;           // BR_AMG_RS or BR_AMG_UA
  INT nvcyc;         // cycles per application
  REAL* t;           // work vector
} precond_br_data;

/*!
 * \fn static void br_amg_setup(AMG_data** mgl, dCSRmat* A, AMG_param* param, INT amg)
 *
 * \brief AMG hierarchy for one block.
 */
static void br_amg_setup(AMG_data** mgl, dCSRmat* A, AMG_param* param, INT amg)
{
  *mgl = amg_data_create(param->max_levels);
  if (amg == BR_AMG_RS) {
    rs_amg_build_hierarchy(*mgl, param, A);
  } else {
    (*mgl)[0].A = dcsr_create(A->row, A->col, A->nnz);
    dcsr_cp(A, &(*mgl)[0].A);
    (*mgl)[0].b = dvec_create(A->row);
    (*mgl)[0].x = dvec_create(A->row);
    // amg_setup_ua allocates the AMLI coefficients; keep one copy
    if (param->amli_coef) { free(param->amli_coef); param->amli_coef = NULL; }
    amg_setup_ua(*mgl, param);
  }
}

/*!
 * \fn static void br_amg_cycle(AMG_data* mgl, AMG_param* param, INT amg, REAL* r, REAL* z)
 *
 * \brief z = one AMG cycle applied to A z = r starting from z = 0.
 */
static void br_amg_cycle(AMG_data* mgl, AMG_param* param, INT amg, REAL* r, REAL* z)
{
  INT m = mgl[0].A.row;
  if (amg == BR_AMG_RS) {
    rs_amg_vcycle_precond(mgl, param, r, z);
    return;
  }
  array_cp(m, r, mgl[0].b.val);
  dvec_set(m, &mgl[0].x, 0.0);
  if (param->cycle_type == AMLI_CYCLE) amli(mgl, param, 0);
  else mgcycle(mgl, param);
  array_cp(m, mgl[0].x.val, z);
}

/*!
 * \fn static void br_amg_apply(AMG_data* mgl, AMG_param* param, INT amg, INT nvcyc, REAL* r, REAL* z)
 *
 * \brief z = nvcyc AMG cycles applied to A z = r starting from z = 0.
 */
static void br_amg_apply(AMG_data* mgl, AMG_param* param, INT amg, INT nvcyc, REAL* r, REAL* z)
{
  INT c, k, m = mgl[0].A.row;
  br_amg_cycle(mgl, param, amg, r, z);
  if (nvcyc > 1) {
    REAL* res = (REAL*)calloc(m, sizeof(REAL));
    REAL* dz = (REAL*)calloc(m, sizeof(REAL));
    for (c = 1; c < nvcyc; c++) {
      array_cp(m, r, res);
      dcsr_aAxpy(-1.0, &mgl[0].A, z, res);
      br_amg_cycle(mgl, param, amg, res, dz);
      for (k = 0; k < m; k++) z[k] += dz[k];
    }
    free(res);
    free(dz);
  }
}

/*!
 * \fn void precond_br_setup(precond_br_data* pd, block_dCSRmat* A, scomplex* sc, REAL mu, INT amg, INT cycle, INT nvcyc)
 *
 * \brief Builds the AMG hierarchies for the velocity blocks and for S.
 *        amg = BR_AMG_RS (V-cycle) or BR_AMG_UA with
 *        cycle = AMLI_CYCLE, W_CYCLE or V_CYCLE.
 */
void precond_br_setup(precond_br_data* pd, block_dCSRmat* A, scomplex* sc, REAL mu,
                      INT amg, INT cycle, INT nvcyc)
{
  INT i, nb = A->brow, np = sc->fem->ns_leaf;
  pd->nb = nb;
  pd->A = A;
  pd->amg = amg;
  pd->nvcyc = nvcyc;

  param_amg_init(&pd->amgparam);
  if (amg == BR_AMG_RS) {
    pd->amgparam.cycle_type = V_CYCLE;
    pd->amgparam.strong_coupled = 0.25;
  } else {
    pd->amgparam.AMG_type = UA_AMG;
    pd->amgparam.cycle_type = cycle;
  }
  pd->amgparam.smoother = SMOOTHER_GS;
  pd->amgparam.presmooth_iter = 2;
  pd->amgparam.postsmooth_iter = 2;
  pd->amgparam.max_levels = 20;
  pd->amgparam.coarse_dof = 100;
  pd->amgparam.print_level = 0;

  // block sizes from the last block column (B_j^T and A_pp)
  pd->off = (INT*)calloc(nb + 1, sizeof(INT));
  for (i = 0; i < nb - 1; i++) pd->off[i + 1] = pd->off[i] + A->blocks[i * nb + i]->row;
  pd->off[nb] = pd->off[nb - 1] + np;
  pd->t = (REAL*)calloc(pd->off[nb], sizeof(REAL));

  pd->mgl = (AMG_data**)calloc(nb, sizeof(AMG_data*));
  for (i = 0; i < nb - 1; i++) br_amg_setup(&pd->mgl[i], A->blocks[i * nb + i], &pd->amgparam, amg);

  // S = M_p/mu - A_pp
  dvector md = dvec_create(np);
  for (i = 0; i < np; i++) md.val[i] = sc->fem->el_vol[i] / mu;
  dCSRmat M = dcsr_create_diagonal_matrix(&md);
  dCSRmat S = dcsr_create(0, 0, 0);
  dCSRmat* App = A->blocks[nb * nb - 1];
  if (App != NULL && App->nnz > 0) dcsr_add(&M, 1.0, App, -1.0, &S);
  else { dcsr_alloc(M.row, M.col, M.nnz, &S); dcsr_cp(&M, &S); }
  br_amg_setup(&pd->mgl[nb - 1], &S, &pd->amgparam, amg);
  dcsr_free(&S);
  dcsr_free(&M);
  dvec_free(&md);
}

/*!
 * \fn void precond_br_lower(REAL* r, REAL* z, void* data)
 *
 * \brief z = P^{-1} r by forward block substitution.
 */
void precond_br_lower(REAL* r, REAL* z, void* data)
{
  precond_br_data* pd = (precond_br_data*)data;
  INT i, j, k, nb = pd->nb;
  INT* off = pd->off;
  REAL* t = pd->t;

  // velocity blocks: z_i = A_ii^{-1} (r_i - sum_{j<i} A_ij z_j)
  for (i = 0; i < nb - 1; i++) {
    array_cp(off[i + 1] - off[i], r + off[i], t);
    for (j = 0; j < i; j++)
      if (pd->A->blocks[i * nb + j] != NULL)
        dcsr_aAxpy(-1.0, pd->A->blocks[i * nb + j], z + off[j], t);
    br_amg_apply(pd->mgl[i], &pd->amgparam, pd->amg, pd->nvcyc, t, z + off[i]);
  }
  // pressure: z_p = S^{-1} (sum_j B_j z_j - r_p)
  INT np = off[nb] - off[nb - 1];
  for (k = 0; k < np; k++) t[k] = -r[off[nb - 1] + k];
  for (j = 0; j < nb - 1; j++)
    if (pd->A->blocks[(nb - 1) * nb + j] != NULL)
      dcsr_aAxpy(1.0, pd->A->blocks[(nb - 1) * nb + j], z + off[j], t);
  br_amg_apply(pd->mgl[nb - 1], &pd->amgparam, pd->amg, pd->nvcyc, t, z + off[nb - 1]);
}

/*!
 * \fn void precond_br_free(precond_br_data* pd)
 */
void precond_br_free(precond_br_data* pd)
{
  INT i;
  // the AMLI coefficients are shared by all blocks: freed once below
  AMG_param ptmp = pd->amgparam;
  ptmp.cycle_type = V_CYCLE;
  for (i = 0; i < pd->nb; i++) {
    if (pd->amg == BR_AMG_RS) {
      rs_amg_free(pd->mgl[i]);
    } else {
      amg_data_free(pd->mgl[i], &ptmp);
      free(pd->mgl[i]);
    }
  }
  if (pd->amgparam.amli_coef) free(pd->amgparam.amli_coef);
  free(pd->mgl);
  free(pd->off);
  free(pd->t);
}

/*!
 * \fn INT solve_br_gmres(block_dCSRmat* A, dvector* b, scomplex* sc, REAL mu, INT amg, INT cycle, dvector* x)
 *
 * \brief Right preconditioned GMRES with the block lower triangular
 *        preconditioner (one AMG cycle per block).  For Dirichlet BC on the
 *        whole boundary the system is singular (constant pressures) and
 *        consistent; on return the pressure is shifted to int p = 0.
 *
 * \return Number of GMRES iterations (negative on failure).
 */
INT solve_br_gmres(block_dCSRmat* A, dvector* b, scomplex* sc, REAL mu, INT amg, INT cycle,
                   dvector* x)
{
  INT k, np = sc->fem->ns_leaf;
  precond_br_data pd;
  precond_br_setup(&pd, A, sc, mu, amg, cycle, 1);
  precond pc;
  pc.data = &pd;
  pc.fct = precond_br_lower;

  // Dirichlet BC on the whole boundary: the kernel of A^T is the constant
  // pressure; remove the mean of the pressure rhs (quadrature of -(g,q))
  // so that the system is consistent
  if (stokes_test != 2) {
    REAL* bp = b->val + b->row - np;
    REAL bm = 0.0;
    for (k = 0; k < np; k++) bm += bp[k];
    for (k = 0; k < np; k++) bp[k] -= bm / np;
  }

  dvec_set(x->row, x, 0.0);
  INT its = bdcsr_pvgmres(A, b, x, &pc, 1e-10, 500, 100, STOP_REL_RES, 0);
  precond_br_free(&pd);

  if (stokes_test == 2) return its;  // traction: the pressure is unique
  REAL* p = x->val + x->row - np;
  REAL pint = 0.0, vol = 0.0;
  for (k = 0; k < np; k++) {
    pint += sc->fem->el_vol[k] * p[k];
    vol += sc->fem->el_vol[k];
  }
  for (k = 0; k < np; k++) p[k] -= pint / vol;
  return its;
}

/*!
 * \fn INT solve_br_umfpack(block_dCSRmat* A, dvector* b, scomplex* sc, dvector* x)
 *
 * \brief UMFPACK with the constraint int p = 0 appended to A and b
 *        (A and b are extended in place); no constraint with traction BC.
 */
INT solve_br_umfpack(block_dCSRmat* A, dvector* b, scomplex* sc, dvector* x)
{
  INT k;
  if (stokes_test != 2) meanzero_pressure_br(A, b, sc);
  dvector y = dvec_create(b->row);
  dvec_set(y.row, &y, 0.0);
  INT flag = block_directsolve_HAZ(A, b, &y, 0);
  for (k = 0; k < x->row; k++) x->val[k] = y.val[k];
  dvec_free(&y);
  return flag;
}
