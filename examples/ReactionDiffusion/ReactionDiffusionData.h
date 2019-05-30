
/******** Data Input Reaction-Diffusion****************************************/
// PDE Coefficients

  //Scalar array of dim^2 x 1, A(x)
void diffusion_coeff(REAL *val,REAL* x,REAL time,void *dim) {
  val[0] = 1.0;
  val[1] = 0.0;
  val[2] = 0.0;
  val[3] = 1.0;
  
}

  //scalar fuction, c(x)
void reaction_coeff(REAL *val,REAL* x,REAL time,void *dim) {
	*val = 0.0;
}

//combines coeff into one struct. Don't need to change anything here
void pde_coeff(REAL *val,REAL* x,REAL time,void *dim){
int d = *((int *) dim);
diffusion_coeff(val, x, time, &d);
reaction_coeff(val+(d)*(d),x,time, &d);
}

// Exact Solution (if you have one)
void exactsol_1D(REAL *val,REAL* x,REAL time,void *param) {
  *val = sin(M_PI*x[0]);
}
void exactsol_2D(REAL *val,REAL* x,REAL time,void *param) {
  *val = sin(M_PI*x[0])*sin(M_PI*x[1]);
}
void exactsol_3D(REAL *val,REAL* x,REAL time,void *param) {
  *val = sin(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
}

// Derivative of Exact Solution (if you have one)
void D_exactsol_1D(REAL *val,REAL* x,REAL time,void *param) {
  *val = M_PI*cos(M_PI*x[0]);
}
void D_exactsol_2D(REAL *val,REAL* x,REAL time,void *param) {
  val[0] = M_PI*cos(M_PI*x[0])*sin(M_PI*x[1]);
  val[1] = M_PI*sin(M_PI*x[0])*cos(M_PI*x[1]);
}
void D_exactsol_3D(REAL *val,REAL* x,REAL time,void *param) {
  val[0] = M_PI*cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
  val[1] = M_PI*sin(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]);
  val[2] = M_PI*sin(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*x[2]);
}

// Right-hand Side
void rhs_1D(REAL *val,REAL* x,REAL time,void *param) {
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu=-666.6;
  reaction_coeff(&myc,x,time,param);
  diffusion_coeff(&mya,x,time,param);
  exactsol_1D(&myu,x,time,param);
  
  *val = (mya*M_PI*M_PI + myc)*myu;		// A=I

}
void rhs_2D(REAL *val,REAL* x,REAL time,void *param) {
  REAL myc=-666.6;
  REAL myu=-666.6;
  reaction_coeff(&myc,x,time,param);
  exactsol_2D(&myu,x,time,param);
  
  *val = (2*M_PI*M_PI + myc)*myu; // A=I
 
}
void rhs_3D(REAL *val,REAL* x,REAL time,void *param) {
  REAL myc=-666.6;
  REAL myu=-666.6;
  reaction_coeff(&myc,x,time,param);
  exactsol_3D(&myu,x,time,param);
  
  *val = (3*M_PI*M_PI + myc)*myu; // A=I

}
// Boundary Conditions
void bc_1D(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu;
  exactsol_1D(&myu,x,time,param);
  *val= myu;
}
void bc_2D(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu;
  exactsol_2D(&myu,x,time,param);
  *val= myu;
}
void bc_3D(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu;
  exactsol_3D(&myu,x,time,param);
  *val= myu;
}