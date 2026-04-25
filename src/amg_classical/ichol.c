/*
 * ichol.c - Incomplete factorizations and AMG+ichol preconditioner
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Incomplete Cholesky (SPD matrices):
 *   ichol_compute    - ichol(0), keeps sparsity of lower triangle of A
 *   icholt_compute   - ichol(tau), drops entries with |l_ij| < tau * ||a_i||
 *   ichol_solve      - triangular solve L L' x = b
 *
 * Incomplete LU (general matrices):
 *   ilu_compute      - ILU(0), keeps sparsity pattern of A
 *   ilut_compute     - ILU(tau), drops entries with |l_ij|,|u_ij| < tau * ||a_i||
 *   ilu_solve        - triangular solve L U x = b
 *
 * Combined AMG+ichol preconditioner:
 *   (1) u1 = S * g              (backslash AMG cycle)
 *   (2) u2 = u1 + L\(L'\(g - A*u1))  (ichol correction)
 *   (3) x  = u2 + S^T * (g - A*u2)   (fwdslash AMG cycle)
 *
 * Uses dCSRmat and AMG_data from HAZmath.
 */
#include "hazmath.h"

/* Sort column indices within each row via double transpose */
static void dcsr_sort_columns(dCSRmat* A) {
  dCSRmat At;
  dcsr_alloc(A->col, A->row, A->nnz, &At);
  dcsr_transz(A, NULL, &At);
  dcsr_free(A);
  dcsr_alloc(At.col, At.row, At.nnz, A);
  dcsr_transz(&At, NULL, A);
  dcsr_free(&At);
}

void ichol_compute(const dCSRmat* Ain, dCSRmat* L) {
  INT n = Ain->row;
  if (n <= 0) return;
  const dCSRmat* A = Ain;

  /* Count lower-triangular entries (including diagonal) */
  INT nnz_lower = 0;
  for (INT i = 0; i < n; i++)
    for (INT k = A->IA[i]; k < A->IA[i + 1]; k++)
      if (A->JA[k] <= i) nnz_lower++;

  dcsr_alloc(n, n, nnz_lower, L);

  /* Extract lower triangle */
  INT pos = 0;
  for (INT i = 0; i < n; i++) {
    L->IA[i] = pos;
    for (INT k = A->IA[i]; k < A->IA[i + 1]; k++) {
      INT j = A->JA[k];
      if (j <= i) {
        L->JA[pos] = j;
        L->val[pos] = A->val[k];
        pos++;
      }
    }
  }
  L->IA[n] = pos;

  /* ichol(0) algorithm with dense work vector */
  REAL* w = (REAL*)calloc((size_t)n, sizeof(REAL));
  INT* w_marker = (INT*)calloc((size_t)n, sizeof(INT));
  for (INT i = 0; i < n; i++) w_marker[i] = -1;

  /* Find diagonal positions */
  INT* diag_pos = (INT*)malloc((size_t)n * sizeof(INT));
  for (INT i = 0; i < n; i++) {
    for (INT k = L->IA[i]; k < L->IA[i + 1]; k++) {
      if (L->JA[k] == i) { diag_pos[i] = k; break; }
    }
  }

  for (INT i = 0; i < n; i++) {
    /* Load row i of L into w */
    for (INT k = L->IA[i]; k < L->IA[i + 1]; k++) {
      w[L->JA[k]] = L->val[k];
      w_marker[L->JA[k]] = i;
    }

    /* For each off-diagonal entry L(i,k) with k < i */
    for (INT kk = L->IA[i]; kk < L->IA[i + 1]; kk++) {
      INT k = L->JA[kk];
      if (k >= i) continue;

      REAL s = w[k];
      for (INT jj = L->IA[k]; jj < L->IA[k + 1]; jj++) {
        INT j = L->JA[jj];
        if (j >= k) break;
        if (w_marker[j] == i) {
          s -= w[j] * L->val[jj];
        }
      }
      REAL Lkk = L->val[diag_pos[k]];
      w[k] = s / Lkk;
    }

    /* Diagonal: L(i,i) = sqrt(w[i] - sum_{j<i} w[j]^2) */
    REAL diag_val = w[i];
    for (INT kk = L->IA[i]; kk < L->IA[i + 1]; kk++) {
      INT j = L->JA[kk];
      if (j >= i) break;
      diag_val -= w[j] * w[j];
    }
    if (diag_val <= 0.0) diag_val = 1e-14;
    w[i] = sqrt(diag_val);

    /* Write back to L */
    for (INT kk = L->IA[i]; kk < L->IA[i + 1]; kk++) {
      L->val[kk] = w[L->JA[kk]];
    }

    /* Clear w */
    for (INT kk = L->IA[i]; kk < L->IA[i + 1]; kk++) {
      w[L->JA[kk]] = 0.0;
      w_marker[L->JA[kk]] = -1;
    }
  }

  free(w); free(w_marker); free(diag_pos);
}

/* Solve L * L' * x = b  =>  L*y = b (fwd), L'*x = y (bwd) */
void ichol_solve(const dCSRmat* L, const REAL* b, REAL* x) {
  INT n = L->row;
  REAL* y = (REAL*)calloc((size_t)n, sizeof(REAL));

  /* Forward solve: L * y = b */
  for (INT i = 0; i < n; i++) {
    REAL s = b[i];
    REAL diag = 1.0;
    for (INT k = L->IA[i]; k < L->IA[i + 1]; k++) {
      INT j = L->JA[k];
      if (j < i)
        s -= L->val[k] * y[j];
      else if (j == i)
        diag = L->val[k];
    }
    y[i] = s / diag;
  }

  /* Backward solve: L' * x = y */
  array_cp(n, y, x);

  /* Find diagonal values */
  REAL* Ldiag = (REAL*)calloc((size_t)n, sizeof(REAL));
  for (INT i = 0; i < n; i++) {
    for (INT k = L->IA[i]; k < L->IA[i + 1]; k++) {
      if (L->JA[k] == i) { Ldiag[i] = L->val[k]; break; }
    }
  }

  for (INT i = n - 1; i >= 0; i--) {
    x[i] /= Ldiag[i];
    for (INT k = L->IA[i]; k < L->IA[i + 1]; k++) {
      INT j = L->JA[k];
      if (j < i) x[j] -= L->val[k] * x[i];
    }
  }

  free(y); free(Ldiag);
}

/* ======================================================================
 *  ichol with threshold: ichol(tau)
 *  Drops entries with |l_ij| < tau * ||a_i||_1 / n_i
 *  Allows fill-in (unlike ichol(0)).
 * ====================================================================== */
void icholt_compute(const dCSRmat* Ain, REAL tau, dCSRmat* L) {
  INT n = Ain->row;
  if (n <= 0) return;
  /* Sort column indices */
  dCSRmat As;
  dcsr_alloc(n, Ain->col, Ain->nnz, &As);
  dcsr_cp((dCSRmat*)Ain, &As);
  dcsr_sort_columns(&As);
  const dCSRmat* A = &As;

  /* Row norms for drop tolerance */
  REAL* row_norm = (REAL*)calloc((size_t)n, sizeof(REAL));
  for (INT i = 0; i < n; i++) {
    for (INT k = A->IA[i]; k < A->IA[i + 1]; k++)
      row_norm[i] += fabs(A->val[k]);
    INT row_len = A->IA[i + 1] - A->IA[i];
    if (row_len > 0) row_norm[i] /= row_len;
  }

  /* Dense work vector: w[j] holds the current row being factored */
  REAL* w = (REAL*)calloc((size_t)n, sizeof(REAL));
  INT* w_nz = (INT*)malloc((size_t)n * sizeof(INT)); /* nonzero indices in w */
  char* w_flag = (char*)calloc((size_t)n, sizeof(char)); /* 1 if w[j] is active */
  REAL* diag = (REAL*)calloc((size_t)n, sizeof(REAL));

  /* Dynamic arrays for L */
  INT nnz_alloc = A->nnz;
  INT* L_IA = (INT*)malloc(((size_t)n + 1) * sizeof(INT));
  INT* L_JA = (INT*)malloc(nnz_alloc * sizeof(INT));
  REAL* L_val = (REAL*)malloc(nnz_alloc * sizeof(REAL));
  INT nnz = 0;

  for (INT i = 0; i < n; i++) {
    L_IA[i] = nnz;

    /* Scatter lower triangle of row i of A into w */
    INT nw = 0;
    for (INT k = A->IA[i]; k < A->IA[i + 1]; k++) {
      INT j = A->JA[k];
      if (j <= i) {
        w[j] = A->val[k];
        if (!w_flag[j]) { w_flag[j] = 1; w_nz[nw++] = j; }
      }
    }

    /* Eliminate columns j = 0, 1, ..., i-1 in order.
     * Only process columns that appear in the stored L rows (no fill-in
     * in the elimination scan — fill occurs only through updates). */
    for (INT j = 0; j < i; j++) {
      if (!w_flag[j]) continue;
      REAL wj = w[j];
      if (wj == 0.0) continue;

      REAL dj = diag[j];
      if (dj <= 0.0) { w[j] = 0.0; continue; }

      REAL lij = wj / dj;
      w[j] = lij;

      /* Update: w[k] -= lij * L(j,k) for off-diagonal k in row j of L */
      for (INT p = L_IA[j]; p < L_IA[j + 1]; p++) {
        INT k = L_JA[p];
        if (k >= j) break;
        w[k] -= lij * L_val[p];
        if (!w_flag[k]) { w_flag[k] = 1; w_nz[nw++] = k; }
      }
    }

    /* Diagonal: L(i,i) = sqrt(a_ii - sum_{j<i} l_ij^2) */
    REAL dval = w[i];
    for (INT kk = 0; kk < nw; kk++) {
      INT j = w_nz[kk];
      if (j < i) dval -= w[j] * w[j];
    }
    /* Diagonal compensation: add dropped entries back to diagonal */
    REAL drop_tol = tau * row_norm[i];
    REAL dropped_sum = 0.0;
    for (INT kk = 0; kk < nw; kk++) {
      INT j = w_nz[kk];
      if (j < i && fabs(w[j]) < drop_tol)
        dropped_sum += w[j] * w[j];
    }
    dval += dropped_sum;

    if (dval <= 0.0) {
      REAL aii = 0.0;
      for (INT k = A->IA[i]; k < A->IA[i + 1]; k++)
        if (A->JA[k] == i) { aii = fabs(A->val[k]); break; }
      dval = (aii > 0.0) ? 1e-2 * aii : 1e-10;
    }
    dval = sqrt(dval);
    diag[i] = dval;

    /* Store off-diagonal entries that pass threshold (sorted by column) */
    for (INT j = 0; j < i; j++) {
      if (!w_flag[j]) continue;
      if (fabs(w[j]) < drop_tol) continue;
      if (nnz >= nnz_alloc) {
        nnz_alloc = nnz_alloc * 3 / 2 + n;
        L_JA = (INT*)realloc(L_JA, nnz_alloc * sizeof(INT));
        L_val = (REAL*)realloc(L_val, nnz_alloc * sizeof(REAL));
      }
      L_JA[nnz] = j;
      L_val[nnz] = w[j];
      nnz++;
    }
    /* Diagonal always kept */
    if (nnz >= nnz_alloc) {
      nnz_alloc = nnz_alloc * 3 / 2 + n;
      L_JA = (INT*)realloc(L_JA, nnz_alloc * sizeof(INT));
      L_val = (REAL*)realloc(L_val, nnz_alloc * sizeof(REAL));
    }
    L_JA[nnz] = i;
    L_val[nnz] = dval;
    nnz++;

    /* Clear w */
    for (INT kk = 0; kk < nw; kk++) {
      w[w_nz[kk]] = 0.0;
      w_flag[w_nz[kk]] = 0;
    }
  }
  L_IA[n] = nnz;

  /* Build output dCSRmat */
  dcsr_alloc(n, n, nnz, L);
  memcpy(L->IA, L_IA, (n + 1) * sizeof(INT));
  memcpy(L->JA, L_JA, nnz * sizeof(INT));
  memcpy(L->val, L_val, nnz * sizeof(REAL));

  free(L_IA); free(L_JA); free(L_val);
  free(w); free(w_nz); free(w_flag); free(diag); free(row_norm);
  dcsr_free(&As);
}

/* ======================================================================
 *  ILU(0): Incomplete LU with no fill
 *  Keeps the sparsity pattern of A.  A = L * U  (L unit lower, U upper)
 *  L and U stored separately as dCSRmat.
 * ====================================================================== */
void ilu_compute(const dCSRmat* Ain, dCSRmat* L, dCSRmat* U) {
  INT n = Ain->row;
  if (n <= 0) return;
  /* Sort column indices */
  dCSRmat As;
  dcsr_alloc(n, Ain->col, Ain->nnz, &As);
  dcsr_cp((dCSRmat*)Ain, &As);
  dcsr_sort_columns(&As);
  const dCSRmat* A = &As;

  /* Count lower (including diagonal for U) and upper entries */
  INT nnz_L = 0, nnz_U = 0;
  for (INT i = 0; i < n; i++)
    for (INT k = A->IA[i]; k < A->IA[i + 1]; k++) {
      if (A->JA[k] < i) nnz_L++;
      else               nnz_U++;
    }
  nnz_L += n;  /* unit diagonal in L */

  dcsr_alloc(n, n, nnz_L, L);
  dcsr_alloc(n, n, nnz_U, U);

  /* Dense row workspace */
  REAL* w = (REAL*)calloc((size_t)n, sizeof(REAL));
  char* w_flag = (char*)calloc((size_t)n, sizeof(char));

  /* Diagonal of U for back-references */
  INT* u_diag_pos = (INT*)calloc((size_t)n, sizeof(INT));

  INT posL = 0, posU = 0;

  for (INT i = 0; i < n; i++) {
    L->IA[i] = posL;
    U->IA[i] = posU;

    /* Load row i of A into w */
    for (INT k = A->IA[i]; k < A->IA[i + 1]; k++) {
      w[A->JA[k]] = A->val[k];
      w_flag[A->JA[k]] = 1;
    }

    /* Elimination: process columns j = 0, 1, ..., i-1 in order */
    for (INT j = 0; j < i; j++) {
      if (!w_flag[j] || w[j] == 0.0) continue;

      REAL u_jj = U->val[u_diag_pos[j]];
      if (fabs(u_jj) < 1e-30) continue;
      REAL mult = w[j] / u_jj;
      w[j] = mult;

      /* w[col] -= mult * U(j, col) for col > j, only if col in pattern */
      for (INT p = u_diag_pos[j] + 1; p < U->IA[j + 1]; p++) {
        INT col = U->JA[p];
        if (w_flag[col])
          w[col] -= mult * U->val[p];
      }
    }

    /* Store L: strictly lower part (sorted) + unit diagonal */
    for (INT k = A->IA[i]; k < A->IA[i + 1]; k++) {
      INT j = A->JA[k];
      if (j < i) {
        L->JA[posL] = j;
        L->val[posL] = w[j];
        posL++;
      }
    }
    L->JA[posL] = i;
    L->val[posL] = 1.0;
    posL++;

    /* Store U: diagonal + strictly upper part */
    u_diag_pos[i] = posU;
    for (INT k = A->IA[i]; k < A->IA[i + 1]; k++) {
      INT j = A->JA[k];
      if (j >= i) {
        U->JA[posU] = j;
        U->val[posU] = w[j];
        posU++;
      }
    }

    /* Clear w */
    for (INT k = A->IA[i]; k < A->IA[i + 1]; k++) {
      w[A->JA[k]] = 0.0;
      w_flag[A->JA[k]] = 0;
    }
  }
  L->IA[n] = posL;
  U->IA[n] = posU;

  free(w); free(w_flag); free(u_diag_pos);
  dcsr_free(&As);
}

/* Solve L * U * x = b  =>  L*y = b (fwd), U*x = y (bwd)
 * L is unit lower triangular, U is upper triangular. */
void ilu_solve(const dCSRmat* L, const dCSRmat* U, const REAL* b, REAL* x) {
  INT n = L->row;
  REAL* y = (REAL*)calloc((size_t)n, sizeof(REAL));

  /* Forward solve: L * y = b  (L has unit diagonal) */
  for (INT i = 0; i < n; i++) {
    REAL s = b[i];
    for (INT k = L->IA[i]; k < L->IA[i + 1]; k++) {
      INT j = L->JA[k];
      if (j < i) s -= L->val[k] * y[j];
    }
    y[i] = s;
  }

  /* Backward solve: U * x = y */
  for (INT i = n - 1; i >= 0; i--) {
    REAL s = y[i];
    REAL diag = 1.0;
    for (INT k = U->IA[i]; k < U->IA[i + 1]; k++) {
      INT j = U->JA[k];
      if (j == i)      diag = U->val[k];
      else if (j > i)  s -= U->val[k] * x[j];
    }
    x[i] = s / diag;
  }

  free(y);
}

/* ======================================================================
 *  ILUT: Incomplete LU with threshold dropping
 *  Drops fill entries with |val| < tau * ||a_i||_1 / n_i
 * ====================================================================== */
void ilut_compute(const dCSRmat* Ain, REAL tau, dCSRmat* L, dCSRmat* U) {
  INT n = Ain->row;
  if (n <= 0) return;
  /* Sort column indices */
  dCSRmat As;
  dcsr_alloc(n, Ain->col, Ain->nnz, &As);
  dcsr_cp((dCSRmat*)Ain, &As);
  dcsr_sort_columns(&As);
  const dCSRmat* A = &As;

  /* Row norms for drop tolerance */
  REAL* row_norm = (REAL*)calloc((size_t)n, sizeof(REAL));
  for (INT i = 0; i < n; i++) {
    for (INT k = A->IA[i]; k < A->IA[i + 1]; k++)
      row_norm[i] += fabs(A->val[k]);
    INT row_len = A->IA[i + 1] - A->IA[i];
    if (row_len > 0) row_norm[i] /= row_len;
  }

  /* Dense workspace */
  REAL* w = (REAL*)calloc((size_t)n, sizeof(REAL));
  char* w_flag = (char*)calloc((size_t)n, sizeof(char));
  INT* w_nz = (INT*)malloc((size_t)n * sizeof(INT));

  /* Dynamic arrays for L and U */
  INT alloc_L = A->nnz, alloc_U = A->nnz;
  INT* LIA = (INT*)malloc(((size_t)n + 1) * sizeof(INT));
  INT* LJA = (INT*)malloc(alloc_L * sizeof(INT));
  REAL* Lval = (REAL*)malloc(alloc_L * sizeof(REAL));
  INT* UIA = (INT*)malloc(((size_t)n + 1) * sizeof(INT));
  INT* UJA = (INT*)malloc(alloc_U * sizeof(INT));
  REAL* Uval = (REAL*)malloc(alloc_U * sizeof(REAL));
  INT posL = 0, posU = 0;

  /* U diagonal positions for back-reference */
  INT* u_diag_pos = (INT*)calloc((size_t)n, sizeof(INT));

  for (INT i = 0; i < n; i++) {
    LIA[i] = posL;
    UIA[i] = posU;

    /* Load row i of A into w */
    INT nw = 0;
    for (INT k = A->IA[i]; k < A->IA[i + 1]; k++) {
      INT j = A->JA[k];
      w[j] = A->val[k];
      w_flag[j] = 1;
      w_nz[nw++] = j;
    }

    /* Elimination: process j = 0, 1, ..., i-1 in order */
    for (INT j = 0; j < i; j++) {
      if (!w_flag[j] || w[j] == 0.0) continue;

      REAL u_jj = Uval[u_diag_pos[j]];
      if (fabs(u_jj) < 1e-30) continue;
      REAL mult = w[j] / u_jj;
      w[j] = mult;

      /* Update: w[col] -= mult * U(j,col) for col > j */
      for (INT p = u_diag_pos[j] + 1; p < UIA[j + 1]; p++) {
        INT col = UJA[p];
        if (w_flag[col]) {
          w[col] -= mult * Uval[p];
        } else {
          /* Fill-in */
          w[col] = -mult * Uval[p];
          w_flag[col] = 1;
          w_nz[nw++] = col;
        }
      }
    }

    REAL drop_tol = tau * row_norm[i];

    /* Store L: strictly lower (sorted) + unit diagonal */
    for (INT j = 0; j < i; j++) {
      if (!w_flag[j]) continue;
      if (fabs(w[j]) >= drop_tol) {
        if (posL >= alloc_L) {
          alloc_L = alloc_L * 3 / 2 + n;
          LJA = (INT*)realloc(LJA, alloc_L * sizeof(INT));
          Lval = (REAL*)realloc(Lval, alloc_L * sizeof(REAL));
        }
        LJA[posL] = j;
        Lval[posL] = w[j];
        posL++;
      }
    }
    /* Unit diagonal */
    if (posL >= alloc_L) {
      alloc_L = alloc_L * 3 / 2 + n;
      LJA = (INT*)realloc(LJA, alloc_L * sizeof(INT));
      Lval = (REAL*)realloc(Lval, alloc_L * sizeof(REAL));
    }
    LJA[posL] = i;
    Lval[posL] = 1.0;
    posL++;

    /* Store U: diagonal (always kept) + upper entries passing threshold */
    u_diag_pos[i] = posU;
    if (posU >= alloc_U) {
      alloc_U = alloc_U * 3 / 2 + n;
      UJA = (INT*)realloc(UJA, alloc_U * sizeof(INT));
      Uval = (REAL*)realloc(Uval, alloc_U * sizeof(REAL));
    }
    UJA[posU] = i;
    Uval[posU] = (fabs(w[i]) > 1e-30) ? w[i] : 1e-14;
    posU++;

    for (INT j = i + 1; j < n; j++) {
      if (!w_flag[j]) continue;
      if (fabs(w[j]) >= drop_tol) {
        if (posU >= alloc_U) {
          alloc_U = alloc_U * 3 / 2 + n;
          UJA = (INT*)realloc(UJA, alloc_U * sizeof(INT));
          Uval = (REAL*)realloc(Uval, alloc_U * sizeof(REAL));
        }
        UJA[posU] = j;
        Uval[posU] = w[j];
        posU++;
      }
    }

    /* Clear w */
    for (INT kk = 0; kk < nw; kk++) {
      w[w_nz[kk]] = 0.0;
      w_flag[w_nz[kk]] = 0;
    }
  }
  LIA[n] = posL;
  UIA[n] = posU;

  /* Build output */
  dcsr_alloc(n, n, posL, L);
  memcpy(L->IA, LIA, (n + 1) * sizeof(INT));
  memcpy(L->JA, LJA, posL * sizeof(INT));
  memcpy(L->val, Lval, posL * sizeof(REAL));

  dcsr_alloc(n, n, posU, U);
  memcpy(U->IA, UIA, (n + 1) * sizeof(INT));
  memcpy(U->JA, UJA, posU * sizeof(INT));
  memcpy(U->val, Uval, posU * sizeof(REAL));

  free(LIA); free(LJA); free(Lval);
  free(UIA); free(UJA); free(Uval);
  free(w); free(w_flag); free(w_nz);
  free(u_diag_pos); free(row_norm);
  dcsr_free(&As);
}

/* ======================================================================
 *  Combined AMG + ichol preconditioner
 * ====================================================================== */

/* Data bundle passed through void* to PCG */
typedef struct {
  AMG_data* mgl;
  AMG_param* param;
  const dCSRmat* A;
  const dCSRmat* L;
} amg_ichol_data;

static void amg_ichol_apply(REAL* g, REAL* x, void* data) {
  const amg_ichol_data* d = (const amg_ichol_data*)data;
  INT n = d->mgl[0].A.row;
  INT nu = d->param->presmooth_iter;

  REAL* u1 = (REAL*)calloc((size_t)n, sizeof(REAL));
  REAL* r1 = (REAL*)calloc((size_t)n, sizeof(REAL));
  REAL* ichol_corr = (REAL*)calloc((size_t)n, sizeof(REAL));
  REAL* r2 = (REAL*)calloc((size_t)n, sizeof(REAL));
  REAL* corr = (REAL*)calloc((size_t)n, sizeof(REAL));

  /* Step 1: u1 = backslash(g) */
  rs_amg_backslash(d->mgl, d->param, 0, g, u1, nu);

  /* Step 2: r1 = g - A*u1; u2 = u1 + ichol_solve(r1) */
  array_cp(n, g, r1); dcsr_aAxpy(-1.0, (dCSRmat*)d->A, u1, r1);
  ichol_solve(d->L, r1, ichol_corr);
  /* u2 stored in x: x = u1 + ichol_corr */
  array_cp(n, u1, x);
  array_axpy(n, 1.0, ichol_corr, x);

  /* Step 3: r2 = g - A*u2; x = u2 + fwdslash(r2) */
  array_cp(n, g, r2); dcsr_aAxpy(-1.0, (dCSRmat*)d->A, x, r2);
  rs_amg_fwdslash(d->mgl, d->param, 0, r2, corr, nu);
  array_axpy(n, 1.0, corr, x);

  free(u1); free(r1); free(ichol_corr); free(r2); free(corr);
}

void rs_amg_ichol_precond(AMG_data* mgl, AMG_param* param,
                          const dCSRmat* A, const dCSRmat* L,
                          const REAL* g, REAL* x) {
  amg_ichol_data d;
  d.mgl = mgl;
  d.param = param;
  d.A = A;
  d.L = L;
  amg_ichol_apply((REAL*)g, x, &d);
}
