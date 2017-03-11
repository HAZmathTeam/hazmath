
#ifdef WITH_SUITESPARSE
#undef WITH_SUITESPARSE
#endif

#ifndef M_PI
#define M_PI 3.141592653589793e+00;
#endif

// PDE Coefficients
void diffusion_coeff(REAL *val,REAL* x, REAL t) {
  // a(x)
  *val = 1.e-4;
  return;
}

// Exact Solution (if you have one)
// Change as needed for different dimensions
void exactsol(REAL *val,REAL* x, REAL t) {
  // 1D
  //*val = sin(M_PI*x[0])*exp(-M_PI*M_PI*time);
  // 2D
  *val = x[0]; //sin(M_PI*x[0])*sin(M_PI*x[1]);
  // 3D
  //  *val=x[2]*x[2];
  //  *val=10.*x[0]*(1.-x[0]) *x[1]*(1.-x[1])*x[2];
  // *val= x[2]+x[0];
  return;
}

void advection(REAL *val, REAL *x,REAL t) {
  val[0] = -1.; // or beta_1(x)
  val[1] = 0.; // or beta_2(x)
  val[2] = 0.; // or beta_3(x)
  return;
}

// Right-hand Side
void f_rhs(REAL *val,REAL* x, REAL t) {
  *val = 1.0;
}

// ANY Boundary Conditions
void bc_any(REAL *val, REAL* x, REAL t) {
  *val= x[0];
  return;
}
