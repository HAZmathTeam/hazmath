/*! \file src/amr/marking.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *   \note routines to mark simplices for refinement
*/

#include "hazmath.h"

INT xins(INT n, INT *nodes, REAL *xs, REAL *xstar)
{
  /* checks if the point given with coordinates xstar[] is in the
  simplex defined by "nodes[] and xs[]" in dimension "n" construct the
  map A from reference simplex and solves linear systems with partial pivoting to determine this */  
  INT n1=n+1,i,j,l0n,ln,j1;
  INT *p=NULL;
  REAL *A=NULL,*xhat=NULL, *piv=NULL;
  A=(REAL *)calloc(n*n,sizeof(REAL));
  xhat=(REAL *)calloc(n,sizeof(REAL));
  piv=(REAL *)calloc(n,sizeof(REAL));
  p=(INT *)calloc(n,sizeof(INT));
  /* fprintf(stdout,"\nj=%d; vertex=%d\n",0+1,nodes[0]+1); fflush(stdout); */
  /* fprintf(stdout,"\nvertex=%d\n",nodes[0]+1); */
  /* for(i=0;i<n;i++){ */
  /*   fprintf(stdout,"xyz=%f ",xs[l0n+i]); */
  /* } */
  /* fprintf(stdout,"\n"); */
  l0n=nodes[0]*n;
  for (j = 1; j<n1;j++){
    /* grab the vertex */
    ln=nodes[j]*n;
    j1=j-1;
    for(i=0;i<n;i++){
      /*A_{ij} = x_{i,nodes(j)}-x0_{i,nodes(j)},j=1:n (skip 0)*/
      A[i*n+j1] = xs[ln+i]-xs[l0n+i];
    }
  }
  //  fflush(stdout);
  for(i=0;i<n;i++){
    xhat[i] = xstar[i]-xs[l0n+i];
  }
  //print_full_mat(n,1,xstar,"xstar");
  //  print_full_mat(n,n,A,"A");
  //  print_full_mat(n,1,xhat,"xhat0");
  solve_pivot(1, n, A, xhat, p, piv);
  //  print_full_mat(n,1,xhat,"xhat");
  /* check the solution if within bounds */
  //  fprintf(stdout,"\nSolution:\n");
  //  print_full_mat(n,1,xhat,"xhat");
  REAL xhatn=1e0,eps0=1e-10,xmax=1e0+eps0;
  INT flag = 0;
  for(j=0;j<n;j++){
    if((xhat[j]<-eps0) || (xhat[j]>xmax)){
      flag=(j+1);
      //      fprintf(stdout,"\nNOT FOUND: xhat(%d)=%e\n\n",flag,xhat[j]);
      break;
    }
    xhatn -= xhat[j];
    if((xhatn<-eps0) || (xhatn>xmax)) {
      flag=n+1;
      //          fprintf(stdout,"\nNOT FOUND: xhat(%d)=%e\n\n",flag,xhatn);
      break;
    }
  }
  if(A) free(A);
  if(xhat) free(xhat);
  if(p) free(p);
  if(piv) free(piv);
  return flag;
}
void marks(INT level,scomplex *sc)
{
  /* mark simplices depending on te value of an estimator */
  /* the estimator here is the aspect ratio of the simplex */
  INT n=sc->n,n1=n+1,ns=sc->ns,nv=sc->nv;
  INT kbadel,ke,i,j,j1,k,p,q,ni,node1,node2,mj,mk;
  INT ne=(INT )((n*n1)/2);
  REAL slmin,slmax,asp,aspmax=-10.;;
  REAL *sl=(REAL *)calloc(ne,sizeof(REAL));
  INT kbad=0;
  for(i = 0;i<ns;i++){
    if(sc->gen[i] <level) continue;
    ni=n1*i;
    ke=0;
    for (j = 0; j<n;j++){      
      mj=sc->nodes[ni+j]*n;
      j1=j+1;
      for (k = j1; k<n1;k++){
	mk=sc->nodes[ni+k]*n;
	sl[ke]=0e0;
	for(p=0;p<n;p++){
	  sl[ke]+=(sc->x[mj+p]-sc->x[mk+p])*(sc->x[mj+p]-sc->x[mk+p]);
	}
	sl[ke]=sqrt(sl[ke]);
	ke++;
      }
    }
    slmin=sl[0];
    slmax=sl[0];
    for(j=1;j<ne;j++){
      if(sl[j]<slmin) slmin=sl[j];
      if(sl[j]>slmax) slmax=sl[j];
    }
    asp=slmax/slmin;
    if(asp>1e1){
      sc->marked[i]=1;
      kbad++;
      //      fprintf(stdout,"\nlev=%d, gen=%d, asp= %e(%e/%e)\n",level,sc->gen[i],asp,slmax,slmin);
      if(asp>aspmax){kbadel=i;aspmax=asp;}
    }
  }
  //  fprintf(stdout,"\nbad:%d, aspectmax=%e (at el=%d)\n",kbad,aspmax,kbadel);
  //  exit(33);
  if(sl)free(sl);
  return;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
void markstar(INT level,scomplex *sc, INT nstar, REAL *xstar)
{
  /* 
     from marked simplices, remove any simplex that does not contain a
     point from xstar[...]. simplices which are initially unmarked are
     left unmarked
  */
  //  fprintf(stdout,"\nNSTAR=%d\n",nstar);fflush(stdout);
  INT n=sc->n,n1=n+1,ns=sc->ns;
  INT istar,jstar,j=-1;
  if(!nstar){
    /* mark everything and bail out. */
    for(j=0;j<ns;j++){
      sc->marked[j]=TRUE;
    }
    return;
  }
  INT *scnjn=NULL,mrkno=0,mrktst=0,flagmrk=0;  
  REAL *xstar0=NULL;
  for(j=0;j<ns;j++){
    flagmrk=0;
    if(sc->child0[j]>=0) {continue;}
    if(!sc->marked[j]) {continue;}
    mrktst++;
    scnjn = sc->nodes+j*n1;
    for(jstar=0;jstar<nstar;jstar++){
      //      fprintf(stdout,"\nel=%d\n",jstar+1);
      xstar0=xstar+jstar*n;
      if(!xins(n,scnjn,sc->x,xstar0)){
	//	fprintf(stdout,"\nel=%d, found: %d; NOmrk=%d",j,jstar+1,mrkno);
	mrkno--;
	sc->marked[j]=jstar+1;
	flagmrk=1;
	break;
      }
    }
    if(flagmrk) continue;
    mrkno++;
    sc->marked[j]=FALSE;
  }
  /* fprintf(stdout,							\ */
  /* 	  "\nat level %d: UN-MARKING %d elements out of %d tested.",	\ */
  /* 	  level, mrkno,mrktst); */
  return;
}
