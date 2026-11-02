/* Shared by the reproducers: a test matrix and an AMG setup. */
#include "hazmath.h"

static dCSRmat poisson_2d(INT n)
{
  const INT N = n * n;
  dCSRmat A = dcsr_create(N, N, 5 * N);
  INT k = 0;
  A.IA[0] = 0;
  for (INT i = 0; i < N; ++i)
  {
    const INT x = i % n, y = i / n;
    if (y > 0)
    {
      A.JA[k] = i - n;
      A.val[k++] = -1.0;
    }
    if (x > 0)
    {
      A.JA[k] = i - 1;
      A.val[k++] = -1.0;
    }
    A.JA[k] = i;
    A.val[k++] = 4.0;
    if (x < n - 1)
    {
      A.JA[k] = i + 1;
      A.val[k++] = -1.0;
    }
    if (y < n - 1)
    {
      A.JA[k] = i + n;
      A.val[k++] = -1.0;
    }
    A.IA[i + 1] = k;
  }
  A.nnz = k;
  return A;
}

/* AMG hierarchy for A with the given parameters (UA setup). */
static AMG_data* setup_ua(const dCSRmat* A, AMG_param* param)
{
  AMG_data* mgl
      = amg_data_create(param->max_levels > 0 ? param->max_levels : 1);
  mgl[0].A = dcsr_create(A->row, A->col, A->nnz);
  dcsr_cp((dCSRmat*)A, &mgl[0].A);
  mgl[0].b = dvec_create(A->row);
  mgl[0].x = dvec_create(A->row);
  amg_setup_ua(mgl, param);
  return mgl;
}
