/*! \file src/approximation/haz_brasil.c
 *
 * \brief This file contains the BRASIL algorithm for best uniform
 *        rational approximation, as well as general support code
 *        for barycentric rational functions.
 *
 * see:
 * C. Hofreither, "An algorithm for best rational approximation based on
 * barycentric rational interpolation." Numerical Algorithms, 2021
 * https://doi.org/10.1007/s11075-020-01042-0
 *
 * This code is a port of the Python version contained in the baryrat
 * package: https://github.com/c-f-h/baryrat
 *
 * Created by Clemens Hofreither on 20210618.
 */
#include "hazmath.h"

/* struct representing a rational function in barycentric form */
typedef struct {
  INT nn;               /* number of nodes - degree is nn-1 */
  REAL *z;              /* array of nodes */
  REAL *w;              /* array of weights */
  REAL *f;              /* array of values */
} barycentric_t;

/* initialize a barycentric function with given number of nodes */
static void bary_init(barycentric_t *bary, INT nn)
{
  bary->nn = nn;
  bary->z = calloc(nn, sizeof(bary->z[0]));
  bary->w = calloc(nn, sizeof(bary->w[0]));
  bary->f = calloc(nn, sizeof(bary->f[0]));
}

/* free allocated memory for barycentric function */
static void bary_free(barycentric_t *bary)
{
  bary->nn = 0;
  free(bary->z);
  free(bary->w);
  free(bary->f);
}

/* evaluate a barycentric rational function at x */
static REAL bary_eval(const barycentric_t *bary, REAL x)
{
  REAL numer = 0.0, denom = 0.0;
  for (INT k = 0; k < bary->nn; ++k) {
    const REAL z_k = bary->z[k];
    const REAL w_k = bary->w[k];
    const REAL f_k = bary->f[k];

    if (x == z_k) {     // are we exactly on a node? (avoid division by zero)
      if (w_k != 0.0)   // if weight w_k is nonzero, r interpolates f_k at x_k
        return f_k;
      else
        continue;       // weight is 0, so no contribution from this node
    } else {            // not on a node, use standard barycentric formula
      const REAL s = w_k / (x - z_k);
      numer += f_k * s;
      denom += s;
    }
  }
  // NB: if there were no contributions (e.g. nn=0), we should return 1.0
  return numer / denom;
}
/* Set the (already initialized to degree n = nn - 1) barycentric function
 * `bary` to interpolate the given nodes and values.
 *
 * nodes and values must be arrays of length 2*n + 1. It is recommended
 * to have the nodes array be in strictly increasing order.
 *
 * The resulting rational function will have degree of numerator and
 * denominator both at most n.
 *
 * ref: (Knockaert 2008), doi:10.1109/LSP.2007.913583
 */
static void bary_interpolate(barycentric_t *bary, const REAL *nodes, const REAL *values, REAL *L)
{
  const INT n = bary->nn - 1;
  if (n < 0) {
    fprintf(stderr, "ERROR: barycentric_t is not initialized\n");
    return;
  }

  /* allocate Loewner matrix (size n x (n + 1))
   *
   * Added an additional last row set to 0 for now because svdgeneral() gives
   * us the singular vector for the n-th (nonzero) singular value otherwise; what
   * we want is the vector spanning the null space (belonging to sigma_{n+1} = 0).
   */
  // compute the Loewner matrix
  for (INT i = 0; i < n; ++i) {
    /* The given nodes are divided into "primary" nodes (with even indices,
     * being n + 1 in number; columns of L) which end up as nodes of the
     * barycentric function and "secondary" nodes (with odd indices, being n in
     * number; rows of L) which are used to set up the constraints for the
     * weight vector.
     */
    const INT b = 2 * i + 1;                    // secondary node index
    for (INT j = 0; j < n + 1; ++j) {
      const INT a = 2 * j;                      // primary node index
      L[i * (n+1) + j] = (values[b] - values[a]) / (nodes[b] - nodes[a]);
    }
  }

  // weight vector = vector in the null space of L
  REAL smin = 0.0;
  INT info = svdgeneral(n + 1, n + 1, L, &smin, bary->w);
  if (info) {
    fprintf(stderr, "\n%% *** HAZMATH WARNING*** IN %s: SVD-INFO IS NOT ZERO; INFO=%lld;\n", __FUNCTION__, (long long )info);
  }
  // set up nodes and values
  for (INT j = 0; j < n + 1; ++j) {
    bary->z[j] = nodes[2 * j];
    bary->f[j] = values[2 * j];
  }
}
// gm=(3-sqrtl(5))*0.5e0;
#define GOLDEN_MEAN 0.3819660112501051

// golden section search to find the maximum of a real function
static REAL golden_search_max(REAL (*f)(REAL, void*), void *data, REAL a, REAL b, REAL *f_max, int num_iter)
{
  // compute first two internal nodes
  REAL x1 = b - (1.0 - GOLDEN_MEAN) * (b - a);
  REAL x2 = a + (1.0 - GOLDEN_MEAN) * (b - a);
  REAL fx1 = f(x1, data);
  REAL fx2 = f(x2, data);

  for (INT k = 0; k < num_iter; ++k) {
    if (fx1 < fx2) {
      a = x1;
      //if (stopping_rule(a, b, tolerance)) break;
      x1 = x2;
      fx1 = fx2;
      x2 = b - GOLDEN_MEAN * (b - a);
      fx2 = f(x2, data);
    } else {
      b = x2;
      //if (stopping_rule(a, b, tolerance)) break;
      x2 = x1;
      fx2 = fx1;
      x1 = a + GOLDEN_MEAN * (b - a);
      fx1 = f(x1, data);
    }
  }
  *f_max = fx1;
  return x1;
}

/* bisection search for the case when the maximum may be on the boundary of the interval */
static REAL boundary_search_max(REAL (*f)(REAL, void*), void *data, REAL a, REAL b, REAL *f_max, int num_iter)
{
  REAL x[2] = { a, b };
  REAL xvals[2] = { f(a, data), f(b, data) };
  const INT max_side = (xvals[0] >= xvals[1] ? 0 : 1);
  const INT other_side = 1 - max_side;

  for (INT k = 0; k < num_iter; ++k) {
    const REAL xm = 0.5 * (x[0] + x[1]);
    const REAL fm = f(xm, data);
    if (fm < xvals[max_side]) {
      // no new maximum found; shrink interval and iterate
      x[other_side] = xm;
      xvals[other_side] = fm;
    } else {
      // found a bracket for the minimum
      return golden_search_max(f, data, x[0], x[1], f_max, num_iter - k);
    }
  }
  *f_max = xvals[max_side];
  return x[max_side];
}

typedef struct {
  REAL16 (*f)(REAL16 x, void *param);   // function to approximate
  void *f_param;                        // user data for f
  barycentric_t *bary;                  // rational approximant
} errfun_data_t;

// computes absolute error | f(x) - r(x) |
static REAL errfun(REAL x, void *data)
{
  errfun_data_t *efdata = (errfun_data_t*)data;
  return fabs(efdata->f(x, efdata->f_param) - bary_eval(efdata->bary, x));
}

// comparator function for sorting REALs via qsort
/* int compare_REAL(const void *a, const void *b) */
/* { */
/*   if (*(REAL*)a < *(REAL*)b) */
/*     return -1; */
/*   else if (*(REAL*)a > *(REAL*)b ) */
/*     return 1; */
/*   else */
/*     return 0; */
/* } */

/*
 * Compute best uniform rational approximation to a scalar function
 * by the BRASIL algorithm.
 *
 * C. Hofreither, "An algorithm for best rational approximation based on
 * barycentric rational interpolation." Numerical Algorithms, 2021
 * https://doi.org/10.1007/s11075-020-01042-0
 *
 * Relatively robust default parameter settings are the following:
 *   init_steps=100
 *   maxiter=1000
 *   step_factor=0.1
 *   max_step_size=0.1
 *   tol=1e-4
 *
 * To use this, first initialize bary to the proper degree (nn = degree+1)
 * and then pass it into this function.
 *
 * Returns:
 *   the computed maximum error
 */
static REAL bary_brasil(
			REAL16 (*f)(REAL16, void*), // function to approximate
			void *param,                // user parameter for f
			REAL a,                     // start point of interval
			REAL b,                     // end point of interval
			barycentric_t *bary,        // output function - should already be
			// initialized with the proper degree
			INT init_steps,             // how many steps of the initialization alg
			INT maxiter,                // maximum number of iterations
			REAL step_factor,           // step factor
			REAL max_step_size,         // maximum allowed step size
			REAL tol,                   // maximum allowed deviation from equioscillation
			INT *iter_brasil,           // number of iterations performed;
			INT print_level             // level of verbosity
			)
{
  // number of golden section search steps to find local maxima
  const INT NUM_GOLDEN_STEPS = 30;
  
  const INT deg = bary->nn - 1;
  assert(deg >= 0);
  const INT nn = 2 * deg + 1;   // number of interpolation nodes
  const INT ni = nn + 1;        // number of intervals
  REAL return_value = -1.0;
  REAL *x  = malloc(nn * sizeof(REAL));  // interpolation nodes
  REAL *fx = malloc(nn * sizeof(REAL));  // f(x) at interpolation nodes
  REAL *local_max_x  = malloc(ni * sizeof(REAL));   // abscissae of local maxima
  REAL *local_max    = malloc(ni * sizeof(REAL));   // values of local maxima
  REAL *intv_lengths = malloc(ni * sizeof(REAL));   // length of intervals

  *iter_brasil=maxiter+1;

  // set up data for the error function |f - r|
  errfun_data_t errfun_data;
  errfun_data.f = f;
  errfun_data.f_param = param;
  errfun_data.bary = bary;

  // initialize interpolation nodes to Chebyshev nodes of 1st kind
  for (INT i = 0; i < nn; ++i)
    x[i] = (1.0 - cos((2*(i+1) - 1.0) / (2*nn) * PI)) / 2 * (b - a) + a;

  REAL *L = calloc((bary->nn + 1) * (bary->nn + 1), sizeof(REAL));
  INT iter=maxiter+init_steps+1;
  for (iter = 0; iter < init_steps + maxiter; ++iter) {
    ///////////////////////////////////////////////
    *iter_brasil=iter;
    //////////////////////////////////////////////
    // construct current rational interpolant
    for (INT i = 0; i < nn; ++i)
      fx[i] = f(x[i], param);
    bary_interpolate(bary, x, fx,L);

    // compute abscissae and values of local maxima in all intervals
    local_max_x[0] = boundary_search_max(errfun, &errfun_data,
        a, x[0], &local_max[0], NUM_GOLDEN_STEPS);
    for (INT k = 0; k < nn - 1; ++k) {
      local_max_x[k+1] = golden_search_max(errfun, &errfun_data,
          x[k], x[k+1], &local_max[k+1], NUM_GOLDEN_STEPS);
    }
    local_max_x[nn] = boundary_search_max(errfun, &errfun_data,
        x[nn-1], b, &local_max[nn], NUM_GOLDEN_STEPS);

    // find intervals with largest and smallest errors
    INT max_intv = darray_max(ni, local_max);
    INT min_intv = darray_min(ni, local_max);

    // compute deviation from equioscillation
    REAL deviation = local_max[max_intv] / local_max[min_intv] - 1.0;
    if (print_level > 2)
      printf("iter=%lld deviation=%e error=%e\n", (long long )iter, deviation, local_max[max_intv]);

    // check for convergence
    int converged = (deviation <= tol);
    if (converged || (iter == init_steps + maxiter - 1)) {
      // converged or reached maximum iterations
      if (!converged)
        fprintf(stderr, "%%\t****warning: BRASIL did not converge; deviation=%e error=%e\n",
            deviation, local_max[max_intv]);
      else if (print_level > 0)
        printf("%%%%BRASIL converged (%lld iter): deviation=%e error=%e\n", (long long )iter, deviation, local_max[max_intv]);
      return_value = local_max[max_intv];
      break;
    }

    if (iter < init_steps) {
      // PHASE 1:
      // move an interpolation node to the point of largest error
      REAL max_err_x = local_max_x[max_intv];
      // we can't move a node to the boundary, so check for that case
      // and move slightly inwards
      if (max_err_x == a)
        max_err_x = (3 * a + x[0]) / 4.0;
      else if (max_err_x == b)
        max_err_x = (x[nn-1] + 3 * b) / 4.0;

      // find the node to move (neighboring the interval with smallest error)
      INT min_j;
      if (min_intv == 0)
        min_j = 0;
      else if (min_intv == nn)
        min_j = nn - 1;
      else {
        // of the two nodes on this interval, choose the farther one
        if (fabs(max_err_x - x[min_intv-1]) < fabs(max_err_x - x[min_intv]))
          min_j = min_intv;
        else
          min_j = min_intv - 1;
      }
    
      // shift node min_j to max_err_x and re-sort the nodes
      x[min_j] = max_err_x;
      dsi_sort(nn,x);
      //      qsort(x, nn, sizeof(x[0]), compare_REAL);

    } else {
      // PHASE 2:
      // global interval size adjustment

      // compute mean error over all intervals
      REAL mean_err = 0.0;
      for (INT k = 0; k < ni; ++k)
        mean_err += local_max[k];
      mean_err /= ni;

      // compute maximum deviation from mean error
      REAL max_dev = 0.0;
      for (INT k = 0; k < ni; ++k)
        max_dev = MAX(max_dev, fabs(local_max[k] - mean_err));

      REAL stepsize = MIN(max_step_size, step_factor * max_dev / mean_err);

      REAL total_length = 0.0;
      for (INT k = 0; k < ni; ++k) {
        // normalize local_max errors to [-1, 1]
        REAL normalized_dev = (local_max[k] - mean_err) / max_dev;
        REAL scaling = pow(1.0 - stepsize, normalized_dev);

        // determine current length of k-th interval
        REAL old_length;
        if (k == 0)         old_length = x[0] - a;
        else if (k == nn)   old_length = b - x[nn-1];
        else                old_length = x[k] - x[k-1];

        // rescale k-th interval
        intv_lengths[k] = scaling * old_length;
        total_length += intv_lengths[k];
      }

      // compute scaling factor such that intervall adds up to b-a
      REAL scaling_factor = (b - a) / total_length;

      // compute new, rescaled nodes
      x[0] = a + intv_lengths[0] * scaling_factor;
      for (INT k = 1; k < ni - 1; ++k)
        x[k] = x[k-1] + intv_lengths[k] * scaling_factor;
    }
  } // end iter loop:
  // cleanup
  free(L);
  free(x);
  free(fx);
  free(local_max_x);
  free(local_max);
  free(intv_lengths);
  return return_value;
}

/**********************************************************************/
/*!
 * \fn REAL get_rpzwf_brasil(REAL16 (*f)(REAL16, void*), void *param,
 *            REAL **rpzwf, REAL a, REAL b, INT deg, INT init_steps,
 *            INT maxiter, REAL step_factor, REAL max_step_size,
 *            REAL tol, INT print_level)
 *
 * \brief Uses the BRASIL algorithm to compute residues, poles, nodes, weights
 *        and values (rpzwf) for the best uniform rational approximation to f().
 *
 * \param f             the function to approximate
 * \param param         user parameters for f
 * \param rpzwf         output arrays; should have space for five REAL*
 * \param a             start point of interval
 * \param b             end point of interval
 * \param deg           desired degree of rational approximation (both numerator and denominator)
 * \param init_steps    how many steps of the initialization algorithm to perform (default: 100)
 * \param maxiter       maximum number of iterations (default: 1000)
 * \param step_factor   step factor (default: 0.1)
 * \param max_step_size maximum allowed step size (default: 0.1)
 * \param tol           maximum allowed deviation from equioscillation (default: 1e-4)
 * \param print_level   level of verbosity
 *
 * \return The maximum error of the computed rational approximation. Its parameters
 *         are returned in the five arrays of rpzwf: rpzwf[0:1] = residues,
 *         rpzwf[2:3] = poles, rpzwf[4] = nodes, rpzwf[5] = weights, rpzwf[6] = values.
 *         These arrays each have length degree + 1 and are allocated as one big
 *         block, i.e., the caller should only call free(rpzwf[0]).
 *
 * \note Please cite https://doi.org/10.1007/s11075-020-01042-0 when using this
 *       code in published research.
 */
REAL get_rpzwf_brasil(
    REAL16 (*f)(REAL16, void*), // function to approximate
    void *param,                // user parameter for f
    REAL **rpzwf,               // output - see haz_aaa.c
    REAL a,                     // start point of interval
    REAL b,                     // end point of interval
    INT deg,                    // desired degree of rational approximation
    INT init_steps,             // how many steps of the initialization alg
    INT maxiter,                // maximum number of iterations
    REAL step_factor,           // step factor
    REAL max_step_size,         // maximum allowed step size
    REAL tol,                   // maximum allowed deviation from equioscillation
    INT *iter_brasil,
    INT print_level             // level of verbosity
)
{
  barycentric_t bary;
  bary_init(&bary, deg + 1);
  *iter_brasil=maxiter+1;
  const REAL error = bary_brasil(f, param, a, b,			\
				 &bary, init_steps, maxiter, step_factor, max_step_size, \
				 tol,iter_brasil,print_level);

  // allocate and fill output data structure
  /* rpzwf[0] = calloc(7 * bary.nn, sizeof(**rpzwf)); */
  /* rpzwf[1] = rpzwf[0] + bary.nn; */
  /* rpzwf[2] = rpzwf[1] + bary.nn; */
  /* rpzwf[3] = rpzwf[2] + bary.nn; */
  /* rpzwf[4] = rpzwf[3] + bary.nn; */

  rpzwf[0]=calloc(7*(bary.nn), sizeof(REAL)); // real (resid) and all;
  rpzwf[1]=rpzwf[0] + bary.nn;// imag (resid)
  rpzwf[2]=rpzwf[1] + bary.nn;// real (poles)
  rpzwf[3]=rpzwf[2] + bary.nn;// real(resid)
  rpzwf[4]=rpzwf[3] + bary.nn;
  rpzwf[5]=rpzwf[4] + bary.nn;
  rpzwf[6]=rpzwf[5] + bary.nn;

  // copy nodes, weights and values
  memcpy(rpzwf[4], bary.z, bary.nn * sizeof(**rpzwf));
  memcpy(rpzwf[5], bary.w, bary.nn * sizeof(**rpzwf));
  memcpy(rpzwf[6], bary.f, bary.nn * sizeof(**rpzwf));

  // compute residues and poles from nodes/weights/values
  //  residues_poles(bary.nn, rpzwf[2], rpzwf[3], rpzwf[4], rpzwf[0], rpzwf[1]);
  residues_poles(bary.nn,rpzwf[4],rpzwf[5],rpzwf[6],rpzwf[0],rpzwf[1],rpzwf[2],rpzwf[3]);
  /* should this be here?
  // this is copied from haz_aaa.c for consistency:
  // rotate the last residue to the front
  */
  /* REAL rswp; */
  /* INT i,m1 = bary.nn - 1; */
  /* rswp = rpzwf[0][m1]; */
  /* for(i = m1; i>0; --i) */
  /*   rpzwf[0][i] = rpzwf[0][i-1]; */
  /* rpzwf[0][0] = rswp;   */
  REAL rswpr,rswpi;
  INT i,m1 = bary.nn - 1;
  rswpr=rpzwf[0][m1];
  rswpi=rpzwf[1][m1];
  for(i=m1;i>0;i--){
    rpzwf[0][i]=rpzwf[0][i-1];
    rpzwf[1][i]=rpzwf[1][i-1];
  }
  rpzwf[0][0]=rswpr;
  rpzwf[1][0]=rswpi;
  bary_free(&bary);
  return error;
}
