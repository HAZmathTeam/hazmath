#ifdef WITH_SUITESPARSE
#undef WITH_SUITESPARSE
#endif

#ifndef M_PI
#define M_PI 3.141592653589793e+00;
#endif

#ifndef EPS0
#define EPS0 1e-10
#endif


// PDE Coefficients
void diffusion_coeff(REAL *val,REAL *x, REAL t,void *param) {
  // K = 3x3 matrix
  // ordered row by row
  INT elem_code = *(INT *) param; //element code
  // one integer and 3 reals: the surface elevation and  
  //scale[k]=(scalex,scaley,scalez,scalet,K_z/K_x), k=1:5;
  memset(val,0,9*sizeof(REAL));
  switch(elem_code) {
    // smaller code, closer to the surface; larger code going towards bottom.
  case 0:
    val[0] = 1e-2;  
    break;
  case 1:
    val[0] = 1e-4;  
    break;
  case 2:
    val[0] = 1e-5;  
    break;
  case 3:
    val[0] = 1e-6;  
    break;
  case 4:
    val[0] = 1e-6;  
    break;
  case 5:
    val[0] = 1e-6;  
    break;
  default:    
    val[0] = 1e0;
    //    fprintf(stderr,"\n***WARNING(in %s): unexpected conductivity value:%i\n",__FUNCTION__,elem_code);
    break;
  }
  val[4] = val[0]; //K_y=K_x=K_{xy}};
  val[8] = val[0]; //K_z=K_{xy}
  //  print_full_mat(3,3,val,__FUNCTION__);
  return;
}

// Exact Solution (used also to set boundary conditions)
void exactsol(REAL *val,REAL* x, REAL t,void *param) {
  //*val = sin(2*M_PI*x[0])*sin(2*M_PI*x[1]);
  *val = (x[0]-x[0]*x[0])*(x[1]-x[1]*x[1]);
  //  *val = x[0]*x[1];
  //  *val = x[0];
  return;
}

void advection(REAL *val, REAL *x,REAL t,void *param) {
  val[0] = 1.; // or beta_1(x)
  val[1] = 1.; // or beta_2(x)
  val[2] = 1.; // or beta_3(x)
  return;
}

void coeff_low_order(REAL *val,
                    REAL *x,
                    REAL time,
                     void *param)
{
    *val = 1.0;
    return;
}
/************************************************************************/
// Right-hand Side
void f_rhs(REAL *val,REAL* x, REAL t,void *param) {
  REAL divj=-1e20,gg=-1e20,uu=-1e20, ad0[3],gruu[3];
  advection(ad0,x,0.0,param);
  exactsol(&uu,x,0.0,param);
  coeff_low_order(&gg,x,0.0,param);
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
void bc_any(REAL *val, REAL* x, REAL t,void *param) {
  REAL uu=-1e20, ad0[3],sj[3];
  //  fprintf(stdout,"%f %s\n",x[0],__FUNCTION__);
  advection(ad0,x,0.0,param);
  exactsol(&uu,x,0.0,param);
  sj[0]=EPS0*(1.-2.*x[0])*(x[1]-x[1]*x[1])+ad0[0]*uu;
  sj[1]=EPS0*(1.-2.*x[1])*(x[0]-x[0]*x[0])+ad0[1]*uu;
  sj[0]=EPS0*(x[1])+ad0[0]*uu;
  sj[1]=EPS0*(x[0])+ad0[1]*uu;
  if(fabs(x[0])< 1e-10){
    // n=(-1,0)
    *val=-sj[0];
    exactsol(val,x,0.0,param);
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

void eafe1(dCSRmat *ain, dvector *rhs,scomplex *sc);
