/*! \file src/amr/amr_coords.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note containing all routines for polar to cartesian; cartesian to
 *  polar, etc. SHould work in dimension n. 
 *
 */
#include "hazmath.h"
/**********************************************************************/
static void coord_perm(SHORT type, INT n,void *x, size_t elsize)
{
  /*
    permutes coordinate arrays (or any array whose elements are of size): 
    if type=1: x[0],x[1]...x[n-1]--> x[1]...x[n-1],x[0] 
    if type=0(inverse): y[0],y[1]...y[n-1]--> y[n-1],y[0]...y[n-2] 
  */
  INT i;  
  void *xx0n=(void *)calloc(1,elsize*sizeof(void));
  if(type){
    memcpy(xx0n,x,elsize);
    for(i=0;i<n;i++)
      memcpy(x+i*elsize,x+(i+1)*elsize,elsize);
    memcpy(x+(n-1)*elsize,xx0n,elsize);
  } else {
    memcpy(xx0n,x+(n-1)*elsize,elsize);
    for(i=(n-1);i>=1;i--)
      memcpy(x+i*elsize,x+(i-1)*elsize,elsize);
    memcpy(x,xx0n,elsize);
  }
  if(xx0n)free(xx0n);
  return;
}
/************************************************************************/
void polar2cart(INT dim, REAL *px, REAL *cx)
{
  // polar is r, theta1,...theta[n-1]; cart is x[0]...x[n-1] px are
  // n-coordinates in polar coord system converts polar coordnates to
  // cartesian in d dimensions.  origin is set at 0,0,0 so if it is
  // different, translation needs to be done after return from here.
  INT i,j;
  if(dim==1) return;
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
  /* cx[0]=rho*cos(px[1]); */
  /* cx[1]=rho*sin(px[1])*cos(px[2]); */
  /* cx[2]=rho*sin(px[1])*sin(px[2]); */
  // the conversion above puts cx[n-1] first, so put it back at the end.  
  coord_perm(1,dim,cx,sizeof(REAL));
  return;
}
/************************************************************************/
INT cart2polar(INT dim, REAL *c,REAL *p)
{
  INT i,dimm1=dim-1;
  REAL rl,r;
  // first put c[n-1] first to agree with the polar ordering;
  coord_perm(0,dim,c,sizeof(REAL));
  r=0.;
  for(i=0;i<dim;i++){
    r+=(c[i]*c[i]);
    p[i]=0e0;
  }  
  if(fabs(r)<1e-14){
    for(i=1;i<dim;i++) p[i]=-1e20;
    return 1;
  }
  r=sqrt(r);
  rl=r;
  INT flag=1;
  for(i=1;i<dim;i++){
    p[i]=acos(c[i-1]/rl);
    if(fabs(sin(p[i]))<1e-14){
      flag=0;
      break;
    }    
    rl/=sin(p[i]);
  }
  if(flag) p[dimm1]=atan2(c[dimm1],c[dimm1-1]);
  p[0]=r;
  // permute c back;
  coord_perm(1,dim,c,sizeof(REAL));
  return 0;
}
/************************************************************************/
