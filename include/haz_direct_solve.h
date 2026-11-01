/*! \file include/haz_direct_solve.h
 *
 * \brief UMFPACK calls that match the width of INT.
 *
 * INT 32-bit (HAZ_INT64==0) -> umfpack_di_*  (int indices)
 * INT 64-bit (HAZ_INT64==1) -> umfpack_dl_*  (SuiteSparse_long indices)
 *
 * Use haz_umfpack_* everywhere instead of umfpack_di_* / umfpack_dl_*,
 * and pass the INT arrays of a dCSRmat directly (no (int *) casts).
 * Include after "hazmath.h" (for INT, HAZ_INT64) and <umfpack.h>.
 */
#ifndef HAZ_DIRECT_SOLVE_H
#define HAZ_DIRECT_SOLVE_H

#if HAZ_INT64
_Static_assert(sizeof(INT) == sizeof(SuiteSparse_long),
               "HAZ_INT64=1 requires sizeof(INT)==sizeof(SuiteSparse_long)");
typedef SuiteSparse_long haz_umf_int;
#define HAZ_UMF_(name) umfpack_dl_##name
#else
_Static_assert(sizeof(INT) == sizeof(int),
               "HAZ_INT64=0 requires sizeof(INT)==sizeof(int)");
typedef int haz_umf_int;
#define HAZ_UMF_(name) umfpack_di_##name
#endif

static inline void haz_umfpack_defaults(double *Control)
{ HAZ_UMF_(defaults)(Control); }

static inline INT haz_umfpack_symbolic(INT nrow, INT ncol, const INT *Ap,
                                       const INT *Ai, const double *Ax,
                                       void **Symbolic, const double *Control,
                                       double *Info)
{
  return (INT)HAZ_UMF_(symbolic)((haz_umf_int)nrow, (haz_umf_int)ncol,
                                 (const haz_umf_int *)Ap,
                                 (const haz_umf_int *)Ai, Ax, Symbolic,
                                 Control, Info);
}

static inline INT haz_umfpack_numeric(const INT *Ap, const INT *Ai,
                                      const double *Ax, void *Symbolic,
                                      void **Numeric, const double *Control,
                                      double *Info)
{
  return (INT)HAZ_UMF_(numeric)((const haz_umf_int *)Ap,
                                (const haz_umf_int *)Ai, Ax, Symbolic,
                                Numeric, Control, Info);
}

static inline INT haz_umfpack_solve(INT sys, const INT *Ap, const INT *Ai,
                                    const double *Ax, double *X,
                                    const double *B, void *Numeric,
                                    const double *Control, double *Info)
{
  return (INT)HAZ_UMF_(solve)((int)sys, (const haz_umf_int *)Ap,
                              (const haz_umf_int *)Ai, Ax, X, B, Numeric,
                              Control, Info);
}

static inline void haz_umfpack_free_symbolic(void **Symbolic)
{ HAZ_UMF_(free_symbolic)(Symbolic); }

static inline void haz_umfpack_free_numeric(void **Numeric)
{ HAZ_UMF_(free_numeric)(Numeric); }

static inline void haz_umfpack_report_status(const double *Control, INT status)
{ HAZ_UMF_(report_status)(Control, (int)status); }

static inline void haz_umfpack_report_info(const double *Control,
                                           const double *Info)
{ HAZ_UMF_(report_info)(Control, Info); }

#endif /* HAZ_DIRECT_SOLVE_H */
