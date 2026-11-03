/* Step 3 check (issue 15) for all six HX apply functions:
 *   - the caller's r is unchanged after an apply;
 *   - apply(r, r) (aliased) gives the same z as apply(r, z);
 *   - a read-only r (mprotect) works;
 *   - the factories leave the caller's coarse_solver alone. */
#include <sys/mman.h>
#include <unistd.h>
#include "poisson.h"
precond* create_precond_hxcurl(dCSRmat*, dCSRmat*, dCSRmat*, SHORT, AMG_param*);
precond* create_precond_hxdiv_3D(dCSRmat*, dCSRmat*, dCSRmat*, dCSRmat*, SHORT,
                                 AMG_param*);
precond* create_precond_hxdiv_2D(dCSRmat*, dCSRmat*, dCSRmat*, SHORT,
                                 AMG_param*);
void free_precond_hxcurl(precond*);
void free_precond_hxdiv(precond*);

static dCSRmat identity(INT n)
{
  dCSRmat I = dcsr_create(n, n, n);
  for (INT i = 0; i < n; ++i)
  {
    I.IA[i] = i;
    I.JA[i] = i;
    I.val[i] = 1.0;
  }
  I.IA[n] = n;
  return I;
}

static void check(const char* name, precond* pc, INT n)
{
  REAL* r = malloc(n * sizeof(REAL));
  REAL* r_orig = malloc(n * sizeof(REAL));
  REAL* z = calloc(n, sizeof(REAL));
  REAL* w = malloc(n * sizeof(REAL));
  for (INT i = 0; i < n; ++i) r[i] = r_orig[i] = 1.0 + (i % 7);

  pc->fct(r, z, pc->data);
  REAL dr = 0, dz = 0, nz = 0;
  for (INT i = 0; i < n; ++i) dr = fmax(dr, fabs(r[i] - r_orig[i]));

  for (INT i = 0; i < n; ++i) w[i] = r_orig[i];
  pc->fct(w, w, pc->data); /* aliased */
  for (INT i = 0; i < n; ++i)
  {
    dz = fmax(dz, fabs(w[i] - z[i]));
    nz = fmax(nz, fabs(z[i]));
  }

  /* read-only r */
  const size_t page = (size_t)sysconf(_SC_PAGESIZE);
  const size_t bytes = ((n * sizeof(REAL) + page - 1) / page) * page;
  REAL* r_ro = mmap(NULL, bytes, PROT_READ | PROT_WRITE,
                    MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  for (INT i = 0; i < n; ++i) r_ro[i] = r_orig[i];
  mprotect(r_ro, bytes, PROT_READ);
  REAL* z2 = calloc(n, sizeof(REAL));
  pc->fct(r_ro, z2, pc->data);
  REAL dro = 0;
  for (INT i = 0; i < n; ++i) dro = fmax(dro, fabs(z2[i] - z[i]));
  munmap(r_ro, bytes);

  printf("%-14s max|z| = %.6e  r changed by %g  aliased differs by %g  read-only r differs by %g\n",
         name, nz, dr, dz, dro);
  free(r); free(r_orig); free(z); free(w); free(z2);
}

int main(void)
{
  dCSRmat A = poisson_2d(16);
  const INT n = A.row;
  dCSRmat I1 = identity(n), I2 = identity(n), I3 = identity(n);
  AMG_param param;

  param_amg_init(&param);
  param.coarse_solver = SOLVER_UMFPACK;
  precond* pc = create_precond_hxcurl(&A, &I1, &I2, PREC_HX_CURL_A, &param);
  printf("coarse_solver before %d, after the setup %d\n", SOLVER_UMFPACK, param.coarse_solver);
  check("HX-curl add", pc, n);
  free_precond_hxcurl(pc);

  param_amg_init(&param);
  pc = create_precond_hxcurl(&A, &I1, &I2, PREC_HX_CURL_M, &param);
  check("HX-curl mult", pc, n);
  free_precond_hxcurl(pc);

  param_amg_init(&param);
  pc = create_precond_hxdiv_3D(&A, &I1, &I2, &I3, PREC_HX_DIV_A, &param);
  check("HX-div 3D add", pc, n);
  free_precond_hxdiv(pc);

  param_amg_init(&param);
  pc = create_precond_hxdiv_3D(&A, &I1, &I2, &I3, PREC_HX_DIV_M, &param);
  check("HX-div 3D mult", pc, n);
  free_precond_hxdiv(pc);

  param_amg_init(&param);
  pc = create_precond_hxdiv_2D(&A, &I1, &I2, PREC_HX_DIV_A, &param);
  check("HX-div 2D add", pc, n);
  free_precond_hxdiv(pc);

  param_amg_init(&param);
  pc = create_precond_hxdiv_2D(&A, &I1, &I2, PREC_HX_DIV_M, &param);
  check("HX-div 2D mult", pc, n);
  free_precond_hxdiv(pc);

  dcsr_free(&A);
  dcsr_free(&I1);
  dcsr_free(&I2);
  dcsr_free(&I3);
  return 0;
}
