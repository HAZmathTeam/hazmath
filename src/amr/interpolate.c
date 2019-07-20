/*! \file src/amr/amr_coords.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note iterpolating routines to map from the reference
 *  PARALLELEPIPED (2^n vertices) to a physical paralelepiped.
 *
 */
#include "hazmath.h"
/************************************************************************/
REAL interp4(cube2simp *c2s, REAL *u, REAL *xhat)
{
  /*INTerpolate d-linearly in d dimensions on the UNIT cube */
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
/************************************************************************/
REAL interp8(cube2simp *c2s, REAL *u, REAL *ue, REAL *xhat)
{
  /*INTerpolate quadratically in d dimensions on the UNIT cube */
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
	//	fprintf(stdout,"\nold:%d:vert:(1-x[%d])",k,i);
      } else {
	phik*=xhat[i];
	//	fprintf(stdout,"\nold:%d:vert:x[%d]",k,i);
      }
    }
    psimid=0.;zmid=0.;
    //    fprintf(stdout,"\nvertex=%d; normal=",k);
    for(i=0;i<dim;++i){
      if(!c2s->bits[kdim+i]){
	psimid-=xhat[i];
	zmid+=c1mid_guess;//was 0.25e00;
	//	fprintf(stdout,"\n%d:vert:(-x[%d])",k,i);
      } else {
	psimid+=xhat[i];
	zmid-=c2mid_guess;//was 0.75
	//	fprintf(stdout,"\n%d:vert:x[%d]",k,i);
      }
    }
    //    fprintf(stdout,"\n%d:zmid=%5.1f",k,zmid);
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
	//	fprintf(stdout,"\n1:(1-x[%d])",i);
      } else {
	phie*=xhat[i];
	//	fprintf(stdout,"\n1:x[%d]",i);
      }
    }
    for(i = 0;i<dim;++i){      
      if(c2s->bits[k2*dim+i]==c2s->bits[k1*dim+i]) continue;
      if(!c2s->bits[k2*dim+i]){
	phie*=(1e0-xhat[i]);
	//	fprintf(stdout,"\n2:(1-x[%d])=%e",i,1.-xhat[i]);
      } else {
	phie*=xhat[i];
	//	fprintf(stdout,"\n2:x[%d]=%e",i,xhat[i]);
      }
    }
    phie*=4e00;
    //    fprintf(stdout,"\nedge=(%d,%d),coord=%e,phie=%e",k1,k2,ue[k],phie);
    se+=ue[k]*phie;
  }
  //  fprintf(stdout,"\n");
  //  fprintf(stdout,"\n***************** xhat=(%e,%e):%e",xhat[0],xhat[1],se);
  //  fprintf(stdout,"\ns=%e",s);
  //  fprintf(stdout,"\n");
  //  exit(99);
  return (s+se);
}
/**********************************************************************/
