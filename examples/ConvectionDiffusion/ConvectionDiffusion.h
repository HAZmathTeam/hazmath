#ifdef WITH_SUITESPARSE
#undef WITH_SUITESPARSE
#endif

#ifndef M_PI
#define M_PI 3.141592653589793e+00;
#endif

#ifndef EPS0
#define EPS0 1e00
#endif


// PDE Coefficients
void diffusion_coeff(REAL *val,REAL* x, REAL t) {
  // a(x)
  //  fprintf(stdout,"%f %s\n",x[0],__FUNCTION__);
  *val = EPS0;
  return;
}

// Exact Solution (used also to set boundary conditions)
void exactsol(REAL *val,REAL* x, REAL t) {
  //*val = sin(2*M_PI*x[0])*sin(2*M_PI*x[1]);
  *val = (x[0]-x[0]*x[0])*(x[1]-x[1]*x[1]);
  //  *val = x[0]*x[1];
  //  *val = x[0];
  return;
}

void advection(REAL *val, REAL *x,REAL t) {
  val[0] = 1.; // or beta_1(x)
  val[1] = 1.; // or beta_2(x)
  val[2] = 1e10; // or beta_3(x)
  return;
}

void coeff_low_order(REAL *val,
                    REAL *x,
                    REAL time)
{
    *val = 1.0;
    return;
}
/************************************************************************/
// Right-hand Side
void f_rhs(REAL *val,REAL* x, REAL t) {
  REAL divj=-1e20,gg=-1e20,uu=-1e20, ad0[3],gruu[3];
  advection(ad0,x,0.0);
  exactsol(&uu,x,0.0);
  coeff_low_order(&gg,x,0.0);
  gruu[0]=(1.-2.*x[0])*(x[1]-x[1]*x[1]);
  gruu[1]=(1.-2.*x[1])*(x[0]-x[0]*x[0]);
  gruu[0]=x[1];
  gruu[1]=x[0];
  REAL divad=0.;
  divj=EPS0*(2.*x[0])*(x[1]-x[1]*x[1]) + EPS0*(2.*x[1])*(x[0]-x[0]*x[0]) \
    - ad0[0]*gruu[0] - ad0[1]*gruu[1] - divad*uu;
   divj= 0. - ad0[0]*gruu[0] - ad0[1]*gruu[1] - divad*uu;
  *val = divj+gg*uu;
  //  *val = gg*uu;
  //  fprintf(stdout,"%f %f %s\n",gg,uu,__FUNCTION__);
  return;
}
// ANY Boundary Conditions
void bc_any(REAL *val, REAL* x, REAL t) {
  REAL uu=-1e20, ad0[3],sj[3];
  //  fprintf(stdout,"%f %s\n",x[0],__FUNCTION__);
  advection(ad0,x,0.0);
  exactsol(&uu,x,0.0);
  sj[0]=EPS0*(1.-2.*x[0])*(x[1]-x[1]*x[1])+ad0[0]*uu;
  sj[1]=EPS0*(1.-2.*x[1])*(x[0]-x[0]*x[0])+ad0[1]*uu;
  sj[0]=EPS0*(x[1])+ad0[0]*uu;
  sj[1]=EPS0*(x[0])+ad0[1]*uu;
  if(fabs(x[0])< 1e-10){
    // n=(-1,0)
    *val=-sj[0];
    exactsol(val,x,0.0);
  } else if(fabs(x[0]-1.) < 1e-10){
    // n=(1,0)
    *val=sj[0];
  } else if(fabs(x[1]) < 1e-10){
    //n=(0,-1)
    *val=-sj[1];
  } else if(fabs(x[1]-1.) < 1e-10){
    //n=(0,1)
    *val=sj[1];
  } else {
    fprintf(stderr,"\n*** Warning: evaluating boundary value at an interior point %f %f\n",x[0],x[1]);
  }
  //  fprintf(stdout,"%e %e %f %s\n",x[0],x[1],*val,__FUNCTION__);
  //  *val=0.;
  /*Dirichlet:*/
  //  exactsol(val,x,0.0);
  return;
}
//void mgraph_wrap(dCSRmat A, dvector rhs, dvector *sol);
void mgraph_wrap(INT idoilu, INT nrow, INT *ia, INT *ja, REAL *a, REAL *rhs, REAL *sol, INT *jareb, REAL *areb, INT *ka);
