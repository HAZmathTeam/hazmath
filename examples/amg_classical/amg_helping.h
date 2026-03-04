/*
 * AMG + Incomplete Cholesky Preconditioner - Header
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Classical Ruge-Stuben AMG hierarchy with selectable smoothers.
 *
 * Uses AMG_data and AMG_param from HAZmath.
 */

#ifndef AMG_H
#define AMG_H

#include <time.h>
#include <math.h>

/* HAZmath main header: types, constants, and function declarations */
#include "hazmath.h"

/* Smoother types: uses SMOOTHER_* from hazmath/macro.h
 *   SMOOTHER_JACOBI   (1)  - damped Jacobi
 *   SMOOTHER_GS       (2)  - lexicographic Gauss-Seidel
 *   SMOOTHER_L1DIAG   (10) - l1-Jacobi: D_l1 = diag(|A|*1)
 *   SMOOTHER_GS_CF    (14) - CF-ordered Gauss-Seidel
 */

/* ======================================================================
 * AMG data structures
 * ====================================================================== */

/* ======================================================================
 * AMG_param field mapping for RS-AMG:
 *
 *   param->strong_coupled   = strength threshold (default 0.25)
 *   param->coarse_dof       = min coarsest-level size
 *   param->max_levels       = max AMG levels
 *   param->smoother         = smoother type (RS_SMOOTHER_*)
 *   param->presmooth_iter   = number of pre-smoothing sweeps
 *   param->postsmooth_iter  = number of post-smoothing sweeps
 *   param->relaxation       = Jacobi damping factor (default 0.5)
 * ====================================================================== */

/* ======================================================================
 * Function prototypes
 * ====================================================================== */

/* --- hierarchy construction (build_hierarchy.c) --- */
void rs_amg_build_hierarchy(AMG_data *mgl, AMG_param *param, const dCSRmat* A);
void rs_amg_rebuild_values(AMG_data *mgl, AMG_param *param, const dCSRmat* A_new);
void rs_amg_free(AMG_data *mgl);

/* --- V-cycle (amg_cycle.c) --- */
void rs_amg_backslash(AMG_data *mgl, AMG_param *param, INT lev,
                      const REAL* b, REAL* x, INT nu);
void rs_amg_fwdslash(AMG_data *mgl, AMG_param *param, INT lev,
                     const REAL* b, REAL* x, INT nu);
void rs_amg_vcycle_precond(AMG_data *mgl, AMG_param *param,
                           const REAL* g, REAL* x);

/* --- ichol (ichol.c) --- */
void ichol_compute(const dCSRmat* A, dCSRmat* L);
void ichol_solve(const dCSRmat* L, const REAL* b, REAL* x);

/* --- combined preconditioner (amg_ichol_precond.c) --- */
void rs_amg_ichol_precond(AMG_data *mgl, AMG_param *param,
                          const dCSRmat* A, const dCSRmat* L,
                          const REAL* g, REAL* x);

/* --- classical RS AMG (rs_classical.c) --- */
void rs_strength(const dCSRmat* A, REAL theta, dCSRmat* S);
void rs_coarsening(const dCSRmat* S, INT* cf, INT* n_coarse);
void rs_standard_interpolation_sparsity(const dCSRmat* S,
                                        const INT* cf, dCSRmat* P);
void rs_standard_interpolation_values(const dCSRmat* A, const dCSRmat* S,
                                      const INT* cf, dCSRmat* P);
void rs_standard_interpolation(const dCSRmat* A, const dCSRmat* S,
                               const INT* cf, dCSRmat* P);

#endif /* AMG_H */
