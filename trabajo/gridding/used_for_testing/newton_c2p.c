#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <getopt.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>
#include "hazmath.h"
#ifndef INT
#define INT int
#endif
#ifndef REAL
#define REAL double
#endif
#ifndef DIM
#define DIM 3
#endif
#ifndef PI
#define PI 3.14159265358979323851281e+00
#endif
void p2c(INT dim, REAL *px, REAL *cx)
{
  // polar is r, theta1,...theta[n-1]; cart is x[0]...x[n-1] px are
  // n-coordinates in polar coord system converts polar coordnates to
  // cartesian in d dimensions.  origin is set at 0,0,0 so if it is
  // different, translation needs to be done after return from here.
  INT i,j;
  if(dim==1) return;
  //  if(polar->type!=1) return;
  REAL rho = px[0];
  REAL cend=rho; 
  memset(cx,0,dim*sizeof(REAL));
  for(i=0;i<(dim-1);i++){
    cx[i]=rho*cos(px[i+1]);
    for(j=0;j<i;j++){
      cx[i]*=sin(px[j+1]);
      //      fprintf(stdout,"\niiiii=%d,jjjjj=%d",i,j+1);
    }
    //    print_full_mat(1,dim,cx,"c1");
    cend*=sin(px[i+1]);    
  }
  cx[dim-1]=cend;
  /* print_full_mat(1,dim,cx,"c11"); */
  /* cx[0]=rho*cos(px[1]); */
  /* cx[1]=rho*sin(px[1])*cos(px[2]); */
  /* cx[2]=rho*sin(px[1])*sin(px[2]); */
  /* print_full_mat(1,dim,cx,"c22"); */
  return;
}
static void jacobip(INT dim, REAL *p, REAL *ajac)
{
  /* 
     computes the jacobian matrix: polar coordinates to cartesian. p
     is the current value of the polar coordinates. ajac is the matrix
     with the jacobian. 
  */
  INT i,j,k;
  REAL pio2=0.5*(PI),r=-1e19,angle=-1e20;
  for(i=0;i<dim;i++) ajac[i*dim+i]=1e0;
  if (p[0]==0) return;
  REAL *c=calloc(dim,sizeof(REAL));
  memset(ajac,0,dim*dim*sizeof(REAL));
  memset(c,0,dim*sizeof(REAL));
  r=p[0];p[0]=1.;
  p2c(dim,p,c);
  //  print_full_mat(dim,1,c,"c2");
  p[0]=r;
  for(j=0;j<dim;j++)ajac[j*dim]=c[j];
  REAL above_d=-p[0];
  for(k=1;k<dim;k++){
    angle=p[k];
    p[k]=pio2-p[k];
    p2c(dim,p,c);
    for(j=k-1;j<dim;j++){
      ajac[j*dim+k]=c[j];
      //      fprintf(stdout,"\nj=%d,k=%d,v=%e",j+1,k+1,c[j]);
    }
    p[k]=angle;  
    above_d*=sin(p[k]);
    ajac[(k-1)*dim+k]=above_d;
  }
  /* REAL *ajac1=calloc(dim*dim,sizeof(REAL)); */
  /* memset(ajac1,0,dim*dim*sizeof(REAL)); */
  /* if(dim==2){ */
  /*   ajac1[0]=cos(p[1]); ajac1[1]=-r*sin(p[1]); */
  /*   ajac1[2]=sin(p[1]); ajac1[3]= r*cos(p[1]); */
  /* } else if(dim==3){ */
  /*   ajac1[0]=cos(p[1]);           ajac1[1]=-r*sin(p[1]);          ajac1[2]=0.; */
  /*   ajac1[3]=sin(p[1])*cos(p[2]); ajac1[4]=r*cos(p[1])*cos(p[2]); ajac1[5]=-r*sin(p[1])*sin(p[2]); */
  /*   ajac1[6]=sin(p[1])*sin(p[2]); ajac1[7]=r*cos(p[1])*sin(p[2]); ajac1[8]=r*sin(p[1])*cos(p[2]); */
  /* } */
  /* REAL cmax=-1e0,cdiff; */
  /* for(i=0; i<dim*dim;i++){ */
  /*   cdiff=fabs(ajac[i]-ajac1[i]); */
  /*   if(cdiff>cmax) cmax=cdiff; */
  /* } */
  /* fprintf(stdout,"\ncmax=%e",cmax); */
  /* if(ajac1) free(ajac1); */
  if(c) free(c);
  return;
}
static void prange(INT dim, REAL *p)
{
  // puts the angles in polar coordinates in [-pi/2,pi/2] and [-pi,pi].
  INT i,dimm1=dim-1;
  REAL pi2=2e0*((REAL )PI);
  /* for(i=1;i<dimm1;i++) */
  /*   p[i]=p[i]-(PI)*(floor(p[i]/(PI))) - (0.5*PI); */
  for(i=1;i<dimm1;i++)
    p[i]=acos(cos(p[i]));
  //  fprintf(stdout,"\nfib=%f",p[dimm1]*180./(PI));
  if(p[dimm1]<0e0)
    p[dimm1]+=pi2*(floor(fabs(p[dimm1])/pi2));
  else
    p[dimm1]-=pi2*(floor(p[dimm1]/pi2));
  //  fprintf(stdout,"\nfia=%f",p[dimm1]*180./(PI));
    
return;
}
INT c2p(INT dim, REAL *c,REAL *p)
{
  // converts cartesian to polar coords using Newton's method; origin
  // is assumed at (0,0,0);
  INT isolv,iter,i,j,k,maxitr=100;
  REAL pi2=2e0*((REAL )PI);
  REAL e=1e20;
  REAL *cw=calloc((3+dim)*dim,sizeof(REAL));
  REAL *f=cw+dim;
  REAL *ajac=f+dim;
  REAL *piv=ajac+dim*dim;
  INT *ip=calloc(dim,sizeof(INT));
  p[0]=0e0;
  for(i=0;i<dim;i++){
    p[0]+=(c[i]*c[i]);
    if(i>0) p[i]=0.25*(PI);
  }
  p[0]=sqrt(p[0]);
  iter=0;
  while(iter<maxitr){
    p2c(dim,p,cw);
    e=0e0;
    for(i=0;i<dim;i++){
      f[i]=c[i]-cw[i];
      e+=f[i]*f[i];
    }
    e=sqrt(e);
    //    fprintf(stdout,"\niter=%4d,e=%16.8e",iter,e);
    if(e<1e-14) break;
    jacobip(dim,p,ajac);
    //    print_full_mat(dim,1,p,"p0");
    //    print_full_mat(dim,dim,ajac,"ajac0");
    //INT solve_pivot(INT dopivot, INT n, REAL *A, REAL *b, INT *p,REAL *piv)
    isolv=solve_pivot(1,dim,ajac,f,ip,piv);
    for(i=0;i<dim;i++) {
      if(i) p[i]+=f[i];
    }
    prange(dim,p);
    iter++;
  }
  if(ip)free(ip);
  if(cw)free(cw);
  if(iter<maxitr){
    fprintf(stdout,"\ne(c2p)=%e",e);
    return 0;
  } else
    return 1;
}
