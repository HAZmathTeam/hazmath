/* Item 2: linear_solver_dcsr_krylov_amg with the AMLI cycle. Its local
 * precond_data got a copy of the AMLI coefficients from param_amg_to_prec
 * that was never freed. Under LSan no leak should be reported. */
#include "poisson.h"

int main(void)
{
  dCSRmat A = poisson_2d(32);
  const INT n = A.row;
  dvector b = dvec_create(n), x = dvec_create(n);
  dvec_set(n, &b, 1.0);
  for (int k = 0; k < 3; ++k)
  {
    AMG_param amgparam;
    param_amg_init(&amgparam);
    amgparam.cycle_type = AMLI_CYCLE;
    amgparam.amli_degree = 2;
    amgparam.print_level = PRINT_NONE;
    linear_itsolver_param itparam;
    param_linear_solver_init(&itparam);
    itparam.linear_print_level = PRINT_NONE;
    dvec_set(n, &x, 0.0);
    INT its = linear_solver_dcsr_krylov_amg(&A, &b, &x, &itparam, &amgparam);
    printf("solve %d: %lld iterations\n", k, (long long)its);
    free(amgparam.amli_coef); /* the caller's copy, if the solver left it */
  }
  dvec_free(&b);
  dvec_free(&x);
  dcsr_free(&A);
  return 0;
}
