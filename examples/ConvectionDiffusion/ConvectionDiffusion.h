#ifdef WITH_SUITESPARSE
#undef WITH_SUITESPARSE
#endif

REAL eps0=1e-6;
REAL pi=3.141592653589793e0;

// PDE Coefficients
void diffusion_coeff(REAL *val,REAL* x) {
  // a(x)
  *val = 1.0;
}

// Exact Solution (if you have one)
// Change as needed for different dimensions
void exactsol(REAL *val,REAL* x0) {
  // 1D
  //*val = sin(M_PI*x[0])*exp(-M_PI*M_PI*time);
  // 2D
  *val = sin(M_PI*x[0])*sin(M_PI*x[1])*exp(-2.0*M_PI*M_PI*time);
  // 3D
  ///  *val=x0[2]*x0[2];
  //  *val=10.*x0[0]*(1.-x0[0]) *x0[1]*(1.-x0[1])*x0[2];
  // *val= x0[2]+x0[0];
  //  *val=exp(-x0[2])*sin(pi*x0[0])*sin(pi*x0[1]);
  //*val = sin(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2])*exp(-3*M_PI*M_PI*time);
}

void advection(const REAL* x0, REAL *advcoeff) {
  advcoeff[0] = 1e2; // or beta_1(x0)
  advcoeff[1] = 0.; // or beta_2(x0)
  advcoeff[2] = -1.; // or beta_3(x0)
}

REAL bernoulli(const REAL z)
{
  // returns B(z) = z/(exp(z)-1)
  double tolb=1e-12,zlarge=256e+0;  
  if (fabs(z) < tolb)
    return (1.-z*0.5); // around 0 this is the asymptotic;
  else if(z<zlarge)
    return (z/(exp(z)-1.));
  else //z>tlarge this is zero pretty much
    return 0.;
}
// Right-hand Side
void myrhs(REAL *val,REAL* x) {
  *val = 0.0;
}

// Boundary Conditions
void bc(REAL *val, REAL* x) {
  *val= 0.0;
}

REAL f_rhs(REAL* x0){
  //  return -2.*eps0+2*x0[2]; // -\Delta u-eps0*u_{tt}
  //  return 2.*x0[2]; // pretend eps0 is 0 and u_{tt} is small. 
  //  return 20.*x0[1]*(1.-x0[1]) + 20.*x0[0]*(1.-x0[0])*x0[2];
  //   return (2.*pi*pi-1.)*exp(-x0[2])*sin(pi*x0[0])*sin(pi*x0[1]);
    return 1e+00;
}
////////////////////////////////////////////////////////////
void DiffusionTensor(const REAL* x0, REAL *alpha) {
  /* 
     alpha is a full matrix, in C stored by rows:
     alpha[i][j] <-- alpha[i*n+j]
   */
  INT i,n=x0.row;
  for (i=0;i<n*n;i++) alpha[i]=0.;
  alpha[0*n+0]=1.; //(1,1)
  alpha[1*n+1]=1.; //(2,2)
  alpha[2*n+2]=1.; //(3,3)
}


