/*
 * build_hierarchy.c - Build multilevel AMG hierarchy
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Classical Ruge-Stuben coarsening using rs_strength, rs_coarsening,
 * rs_standard_interpolation from rs_classical.c.
 *
 * Uses AMG_data from HAZmath; RS-specific data stored in wdata.
 */
#include "hazmath.h"

/* Compute l1 row norms: d_i = sum_j |a_ij| */
static void compute_l1_diag(const dCSRmat* A, REAL** out) {
  INT n = A->row;
  REAL* d = (REAL*)malloc(n * sizeof(REAL));
  for (INT i = 0; i < n; i++) {
    REAL s = 0.0;
    for (INT k = A->IA[i]; k < A->IA[i + 1]; k++)
      s += fabs(A->val[k]);
    d[i] = s;
  }
  *out = d;
}

/* Allocate and zero-init RS auxiliary data for a level */
static rs_level_aux* rs_aux_create(void) {
  rs_level_aux* aux = (rs_level_aux*)calloc(1, sizeof(rs_level_aux));
  return aux;
}

/* ======================================================================
 *  Classical Ruge-Stuben hierarchy
 * ====================================================================== */
void rs_amg_build_hierarchy(AMG_data* mgl, AMG_param* param, const dCSRmat* A) {
  REAL threshold = param->strong_coupled;
  INT min_size   = param->coarse_dof;
  INT max_levels = param->max_levels;

  /* Copy A to level 0 */
  dcsr_alloc(A->row, A->col, A->nnz, &mgl[0].A);
  dcsr_cp((dCSRmat*)A, &mgl[0].A);

  rs_level_aux* aux0 = rs_aux_create();
  compute_l1_diag(&mgl[0].A, &aux0->l1_diag);
  mgl[0].wdata = aux0;

  fprintf(stderr,
    "Building RS hierarchy (threshold=%.2f, min_size=%d)\n",
    threshold, min_size);
  fprintf(stderr, "  Level 1: n = %d, nnz = %d\n", A->row, A->nnz);

  INT nlev = 1;
  const dCSRmat* Ac = &mgl[0].A;

  while (Ac->row >= min_size && nlev < max_levels) {
    INT n_fine = Ac->row;
    rs_level_aux* cur_aux = RS_AUX(mgl, nlev - 1);

    /* Step 1: Strength of connection */
    rs_strength(Ac, threshold, &cur_aux->A_filtered);

    /* Step 2: Classical RS coarsening */
    INT* cf = (INT*)calloc(n_fine, sizeof(INT));
    INT n_coarse = 0;
    rs_coarsening(&cur_aux->A_filtered, cf, &n_coarse);

    if (n_coarse == 0 || n_coarse >= n_fine) {
      fprintf(stderr, "  Coarsening stalled at level %d (n_coarse=%d)\n",
              nlev, n_coarse);
      free(cf);
      break;
    }

    /* Populate mis[] and isolated[] */
    cur_aux->mis = (INT*)calloc(n_fine, sizeof(INT));
    cur_aux->isolated = (INT*)calloc(n_fine, sizeof(INT));
    for (INT i = 0; i < n_fine; i++)
      cur_aux->mis[i] = (cf[i] == 1) ? 1 : 0;

    /* Build CF ordering: C-points first, then F-points */
    cur_aux->cf_order = (INT*)malloc(n_fine * sizeof(INT));
    INT idx = 0;
    for (INT i = 0; i < n_fine; i++)
      if (cf[i] == 1) cur_aux->cf_order[idx++] = i;
    for (INT i = 0; i < n_fine; i++)
      if (cf[i] != 1) cur_aux->cf_order[idx++] = i;

    /* Step 3: Standard interpolation */
    rs_standard_interpolation(&mgl[nlev - 1].A, &cur_aux->A_filtered, cf, &mgl[nlev - 1].P);

    free(cf);

    /* Step 4: Galerkin coarse-grid operator Ac = P^T A P */
    dCSRmat Pt;
    dcsr_alloc(mgl[nlev - 1].P.col, mgl[nlev - 1].P.row, mgl[nlev - 1].P.nnz, &Pt);
    dcsr_transz(&mgl[nlev - 1].P, NULL, &Pt);

    dCSRmat Ac_new;
    dcsr_rap(&Pt, &mgl[nlev - 1].A, &mgl[nlev - 1].P, &Ac_new);
    dcsr_free(&Pt);

    /* Symmetrize: Ac = (Ac + Ac') / 2 */
    {
      dCSRmat Ac_t;
      dcsr_alloc(Ac_new.col, Ac_new.row, Ac_new.nnz, &Ac_t);
      dcsr_transz(&Ac_new, NULL, &Ac_t);

      dCSRmat Ac_sum;
      dcsr_add(&Ac_new, 0.5, &Ac_t, 0.5, &Ac_sum);

      dcsr_free(&Ac_new);
      dcsr_free(&Ac_t);
      Ac_new = Ac_sum;
    }

    fprintf(stderr, "  Level %d: n = %d, nnz = %d (ratio %.2f)\n",
            nlev + 1, n_coarse, Ac_new.nnz, (REAL)n_fine / n_coarse);

    /* Store coarse level */
    mgl[nlev].A = Ac_new;
    rs_level_aux* next_aux = rs_aux_create();
    compute_l1_diag(&mgl[nlev].A, &next_aux->l1_diag);
    mgl[nlev].wdata = next_aux;

    nlev++;
    Ac = &mgl[nlev - 1].A;
  }

  mgl[0].num_levels = nlev;

  /* Coarsest-level solver: try ichol, fall back to direct if it fails */
  rs_level_aux* coarsest_aux = RS_AUX(mgl, nlev - 1);
  memset(&coarsest_aux->L_ichol, 0, sizeof(dCSRmat));
  ichol_compute(&mgl[nlev - 1].A, &coarsest_aux->L_ichol);
  if (coarsest_aux->L_ichol.nnz > 0) {
    fprintf(stderr, "RS hierarchy complete: %d levels (coarsest ichol nnz=%d)\n\n",
            nlev, coarsest_aux->L_ichol.nnz);
  } else {
    fprintf(stderr, "RS hierarchy complete: %d levels (coarsest: direct solve)\n\n",
            nlev);
  }
}

/* ======================================================================
 *  Rebuild hierarchy values for a spectrally equivalent matrix
 * ====================================================================== */
void rs_amg_rebuild_values(AMG_data* mgl, AMG_param* param, const dCSRmat* A_new) {
  INT nlev = mgl[0].num_levels;
  (void)param;

  /* Free old ichol on coarsest level */
  dcsr_free(&RS_AUX(mgl, nlev - 1)->L_ichol);

  /* Replace level 0 matrix */
  dcsr_free(&mgl[0].A);
  dcsr_alloc(A_new->row, A_new->col, A_new->nnz, &mgl[0].A);
  dcsr_cp((dCSRmat*)A_new, &mgl[0].A);
  rs_level_aux* aux0 = RS_AUX(mgl, 0);
  if (aux0->l1_diag) { free(aux0->l1_diag); aux0->l1_diag = NULL; }
  compute_l1_diag(&mgl[0].A, &aux0->l1_diag);

  fprintf(stderr,
    "Rebuilding AMG values (reusing hierarchy structure)\n");
  fprintf(stderr, "  Level 1: n = %d, nnz = %d\n", A_new->row, A_new->nnz);

  for (INT k = 0; k < nlev - 1; k++) {
    rs_level_aux* cur_aux = RS_AUX(mgl, k);
    INT n = mgl[k].A.row;

    /* Reconstruct cf array from mis */
    INT* cf = (INT*)malloc(n * sizeof(INT));
    for (INT i = 0; i < n; i++)
      cf[i] = cur_aux->mis[i] ? 1 : -1;

    /* Recompute P values (sparsity pattern is preserved) */
    rs_standard_interpolation_values(&mgl[k].A, &cur_aux->A_filtered, cf, &mgl[k].P);

    free(cf);

    /* Recompute coarse-grid operator: Ac = P^T A P */
    dCSRmat Pt;
    dcsr_alloc(mgl[k].P.col, mgl[k].P.row, mgl[k].P.nnz, &Pt);
    dcsr_transz(&mgl[k].P, NULL, &Pt);

    dCSRmat Ac_new;
    dcsr_rap(&Pt, &mgl[k].A, &mgl[k].P, &Ac_new);
    dcsr_free(&Pt);

    /* Symmetrize: Ac = (Ac + Ac') / 2 */
    {
      dCSRmat Ac_t;
      dcsr_alloc(Ac_new.col, Ac_new.row, Ac_new.nnz, &Ac_t);
      dcsr_transz(&Ac_new, NULL, &Ac_t);

      dCSRmat Ac_sum;
      dcsr_add(&Ac_new, 0.5, &Ac_t, 0.5, &Ac_sum);

      dcsr_free(&Ac_new);
      dcsr_free(&Ac_t);
      Ac_new = Ac_sum;
    }

    fprintf(stderr, "  Level %d: n = %d, nnz = %d\n",
            k + 2, Ac_new.row, Ac_new.nnz);

    /* Replace coarse level matrix */
    dcsr_free(&mgl[k + 1].A);
    mgl[k + 1].A = Ac_new;
    rs_level_aux* next_aux = RS_AUX(mgl, k + 1);
    if (next_aux->l1_diag) { free(next_aux->l1_diag); next_aux->l1_diag = NULL; }
    compute_l1_diag(&mgl[k + 1].A, &next_aux->l1_diag);
  }

  /* Recompute coarsest-level solver */
  rs_level_aux* coarsest_aux = RS_AUX(mgl, nlev - 1);
  if (coarsest_aux->L_ichol.nnz > 0) {
    dcsr_free(&coarsest_aux->L_ichol);
    ichol_compute(&mgl[nlev - 1].A, &coarsest_aux->L_ichol);
  }
  fprintf(stderr, "AMG values rebuilt: %d levels\n", nlev);
}

/* ====================================================================== */
void rs_amg_free(AMG_data* mgl) {
  INT nlev = mgl[0].num_levels;
  for (INT k = 0; k < nlev; k++) {
    dcsr_free(&mgl[k].A);
    if (k < nlev - 1)
      dcsr_free(&mgl[k].P);

    rs_level_aux* aux = RS_AUX(mgl, k);
    if (aux) {
      if (k < nlev - 1)
        dcsr_free(&aux->A_filtered);
      if (aux->mis) free(aux->mis);
      if (aux->isolated) free(aux->isolated);
      if (aux->cf_order) free(aux->cf_order);
      if (aux->l1_diag) free(aux->l1_diag);
      if (k == nlev - 1) dcsr_free(&aux->L_ichol);
      free(aux);
      mgl[k].wdata = NULL;
    }
  }
  mgl[0].num_levels = 0;
  free(mgl);
}
