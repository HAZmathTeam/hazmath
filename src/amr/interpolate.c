/*! \file src/amr/interpolate.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note iterpolating routines to map from the reference
 *  PARALLELEPIPED (2^n vertices) to a physical paralelepiped. Also
 *  routine to interpolate any data on a uniform grid.
 *
 */
#include "hazmath.h"
/**********************************************************************/
/*!
 * \fn REAL interp4(cube2simp *c2s, REAL *u, REAL *xhat)
 *
 * \brief Interpolate d-linearly in d dimensions on the UNIT cube
 *        given a point with coordinates xhat[j],j=0:dim-1 in the unit
 *        cube. and function values u[k] at the vertices of the unit
 *        cube (k=0:2^{dim}-1) this returns the value of the
 *        d-linear interpolant at xhat[].
 *
 * \param 
 *
 * \return
 *
 * \note
 *
 */
REAL interp4(cube2simp *c2s, REAL *u, REAL *xhat)
{
  INT dim=c2s->n,k,i,kdim;
  REAL phik;
  REAL s=0.;
  for(k = 0;k<c2s->nvcube;k++){
    kdim=k*dim;
    phik=1e0;
    for(i = 0;i<dim;++i){      
      if(!c2s->bits[kdim+i])
	phik*=(1e0-xhat[i]);
      else	
	phik*=xhat[i];
    }
    s+=u[k]*phik;
  }
  return s;
}
/**********************************************************************/
/*!
 * \fn REAL interp8(cube2simp *c2s, REAL *u, REAL *ue, REAL *xhat) 
 *
 * \brief Interpolate quadratically in every dimension on the UNIT
 *        cube. It needs improvement for the polar coordinates
 *        interpolation.
 *
 * \param 
 *
 * \return
 *
 * \note
 *
 */
REAL interp8(cube2simp *c2s, REAL *u, REAL *ue, REAL *xhat)
{
  INT dim=c2s->n,k,i,k1,k2,kdim;
  REAL c1mid_guess=0.5/((REAL )dim);
  REAL c2mid_guess=1.-c1mid_guess;
  REAL phik,zmid,psimid;
  REAL s=0.;
  for(k = 0;k<c2s->nvcube;k++){
    kdim=k*dim;
    phik=1e0;
    for(i = 0;i<dim;++i){      
      if(!c2s->bits[kdim+i]){
	phik*=(1e0-xhat[i]);
      } else {
	phik*=xhat[i];
      }
    }
    psimid=0.;zmid=0.;
    for(i=0;i<dim;++i){
      if(!c2s->bits[kdim+i]){
	psimid-=xhat[i];
	zmid+=c1mid_guess;//was 0.25e00;
      } else {
	psimid+=xhat[i];
	zmid-=c2mid_guess;//was 0.75
      }
    }
    phik*=2e00*(zmid+psimid);
    s+=u[k]*phik;
  }  
  REAL se,phie;
  se=0.;
  for(k=0;k<c2s->ne;k++){
    k2=c2s->edges[2*k];    
    k1=c2s->edges[2*k+1];
    //    fprintf(stdout,"\nmid=(%d,%d);",k1,k2);
    phie=1e0;
    for(i = 0;i<dim;++i){      
      if(!c2s->bits[k1*dim+i]){
	phie*=(1e0-xhat[i]);
      } else {
	phie*=xhat[i];
      }
    }
    for(i = 0;i<dim;++i){      
      if(c2s->bits[k2*dim+i]==c2s->bits[k1*dim+i]) continue;
      if(!c2s->bits[k2*dim+i]){
	phie*=(1e0-xhat[i]);
      } else {
	phie*=xhat[i];
      }
    }
    phie*=4e00;
    se+=ue[k]*phie;
  }
  return (s+se);
}
/**********************************************************************/
/*!
 * \fn void data_transform(const INT nv, const int m, REAL *data, 
 *                         REAL *xodst, REAL *xndst)
 *
 * \brief Transform each column of a REAL data given in (nv x m)
  *        matrix to the intervals xodst[j],xndst[j], j=1:m
*
  * \param 
  *
  * \return
  *
  * \note
  *
  */
void data_transform(const INT nv, const INT m,			\
		    REAL *data, REAL *xodst, REAL *xndst)
{
  INT i,j;
  REAL *xo=(REAL *)calloc(2*m,sizeof(REAL));
  REAL *xn=xo+m;
  REAL dxsrc,a,b;
  for(j=0;j<m;j++){
    xn[j]=xo[j]=data[j];
    for(i=1;i<nv;i++){
      if(data[m*i+j]<xo[j]) xo[j]=data[m*i+j];
      if(data[m*i+j]>xn[j]) xn[j]=data[m*i+j];
    }
  }  
  for(j=0;j<m;j++){
    dxsrc=xn[j]-xo[j];
    a=(xndst[j]-xodst[j])/dxsrc;
    b=(xodst[j]*xn[j]-xndst[j]*xo[j])/dxsrc;    
    for(i=0;i<nv;i++){
      data[m*i+j] = a*data[m*i+j] +b;
    }
  }
  if(xo) free(xo);
  return; 
}
/**********************************************************************/
/*!
 * \fn void interp1(const INT dimbig, REAL *fi, unigrid *ug,REAL *x,
 * const INT nvert, INT *mask)
 *
 * \brief interpolates using data from lattice grid. Anything outside
 *  the lattice is moved to the closest lattice boundary in L_1 and
 *  interpolated too.
 *
 * \param 
 *
 * \return
 *
 * \note
 *
 */
void interp1(const INT dimbig, REAL *fi, unigrid *ug,	\
	     REAL *x, const INT nvert, INT *mask)
{  
  /*  
      find the values of ug->data at x by interpolation.  ug->data is
      given on a uniform grid ug.
  */
  INT dim = ug->n,i,j,k,kdim,kf;
  REAL xoj;
  //  unsigned int bi,bi1;
  REAL *xo = ug->xo,*dx = ug->dx;// *xn = ug->xn;
  unsigned INT *bits=ug->bits;
  INT nvcube=ug->nvcube;
  //    fprintf(stderr,"\nnvcube=%d\n",nvcube);
  cube2simp *c2s=cube2simplex(dim);
  INT *nd=ug->ndiv;
  INT *mo=(INT *)calloc(dim,sizeof(INT));
  REAL *u=(REAL *)calloc(nvcube,sizeof(REAL));
  REAL *xhat=(REAL *)calloc(dim,sizeof(REAL));
  /*NOTE: nd[] are the number of DIVISIONS (not number of vertices)
    for every direction */
  unsigned INT found;
  REAL eps0=1e-8;
  for (i = 0; i<nvert;i++){
    if(mask != NULL) {
      if(mask[i]) continue; /* mask false means do not mask */
    }
    found=FALSE;
    for(j=0;j<dim;j++){ 
      xhat[j] = (x[dimbig*i+j]-xo[j])/dx[j];
      mo[j] = (INT )floor(xhat[j]-eps0);
      if(mo[j]<0 && fabs(xhat[j])<eps0*2.) {mo[j]=0;}
      else if(mo[j]<0) {mo[j]=0; continue;}
      if(mo[j]>nd[j]) {mo[j]=nd[j]; continue;}
      found=TRUE;
    }   
    if(!found){
      fprintf(stdout,"\nvertex=%i NOT found. x (xhat)=",i); fflush(stdout);
      for(j=0;j<dim;j++){ 
	fprintf(stdout," %12.6e (%10.5f)",x[dimbig*i+j],xhat[j]);
      }
      fprintf(stdout,"\n");
      //    }else{
      //      fprintf(stdout,"\nvertex=%i found.",i); fflush(stdout);
    }
    //    if(!found)continue;
    for(j=0;j<dim;j++){
      xoj=xo[j]+((REAL )mo[j])*dx[j];
      xhat[j] = (x[dimbig*i+j]-xoj)/dx[j];
    }
    /* for(j=0;j<dim;j++){ */
    /*   fprintf(stdout,"xo=%10.5e; dx=%12.3e ; xhat=%10.5f (%10.5f)\n ",xo[j],dx[j],xhat[j],x[dim*i+j]); fflush(stdout); */
    /* } */
    for(k = 0;k<nvcube;k++){
      kdim=k*dim;	
      for (j=0;j<dim; j++){mo[j]+=(INT )(!bits[kdim+j]);}
      /* compute the global index */
      kf=num_lattice(mo,dim,nd);
      /* return to the previous state */
      for (j=0;j<dim; j++){mo[j]-=(INT )(!bits[kdim+j]); }
      /* set the local values at the vertices */
      u[k]=ug->data[kf];
    }
    /* for(j=0;j<nvcube;j++){ */
    /*   fprintf(stdout,"u=%10.5e ",u[j]); fflush(stdout); */
    /* } */
    fi[i] = interp4(c2s,u,xhat);
    /*    fprintf(stdout,"Err(%i)=%12.3e;\n",i+1,fi[i]-ff(dim,(x+dim*i)));*/
  }
  if(u) free(u);
  if(mo) free(mo);
  if(xhat) free(xhat);
  return;
}
/*EOF*/
