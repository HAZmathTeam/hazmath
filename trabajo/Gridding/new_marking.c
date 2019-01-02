/*! \file  marking.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *   \note routines to mark simplices for refinement
*/

#include "hazmath.h"
#include "grid_defs.h"
void marks(scomplex *sc)
{
  /* mark simplices depending on the value of an estimator */
  /* the estimator here is the aspect ratio of the simplex */
  INT n=sc->n,n1=n+1,ns=sc->ns,nv=sc->nv,level=sc->level;
  INT kbadel,ke,i,j,j1,k,p,q,ni,node1,node2,mj,mk;
  INT ne=(INT )((n*n1)/2);
  REAL slmin,slmax,asp,aspmax=-10.;;
  REAL *sl=(REAL *)calloc(ne,sizeof(REAL));
  INT kbad=0;
  for(i = 0;i<ns;i++){
    if(sc->gen[i] < level) continue;
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
void markstar(scomplex *sc,dvector *w)
{
  /* 
     from marked simplices, remove any simplex that does not contain a
     point from w[...]. simplices which are initially unmarked are
     left unmarked.
  */
  if(!(w->row) || !(w->val)) return;
  INT n=sc->n,n1=n+1,ns=sc->ns,nstar=w->row;
  //  fprintf(stdout,"\nNSTAR=%d\n",nstar);fflush(stdout);
  INT istar,jstar,j=-1;
  REAL *xstar = w->val;
  INT *scnjn=NULL,mrkno=0,mrktst=0,flagmrk=0;  
  REAL *xstar0=NULL;
  for(j=0;j<ns;j++){
    flagmrk=0;
    if(sc->child0[j]>=0) {continue;}
    if(!sc->marked[j]) {continue;}
    mrktst++;
    scnjn = sc->nodes+j*n1; //beginning of local vertex numbering for
			    //simplex j.
    for(jstar=0;jstar<nstar;jstar++){
      //      fprintf(stdout,"\nel=%d\n",jstar+1);
      xstar0=xstar+jstar*n; 
      if(!xins(n,scnjn,sc->x,xstar0)){
	//fprintf(stdout,"\nel=%d, found: %d; NOmrk=%d",j,jstar+1,mrkno);
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
  //  for(j=0;j<sc->ns;j++){
  //    if(sc->marked[j])
      //      fprintf(stdout,"\nj,mark=%d %d; gen=%d",j,sc->marked[j],sc->gen[j]);
  //  }
    /* fprintf(stdout,							\  */
    /* 	  "\nat level %d: UN-MARKING %d elements out of %d tested.",	\  */
    /* 	  level, mrkno,mrktst);  */
  return;
}
unsigned int markeql(scomplex *sc, dvector *w)
{
  /* input: simplicial complex sc and a vector with error estimate on
     each simplex. In general the values of w should be
     nonnegative. We take the absolute value of the elements of w
     below to make sure we have nonnegative values*/
  /* ON RETURN: a nonzero value indicates wrong input and all simplices are marked */
  /* mark simplices equilibrating errors given by a dvector w*/
  /* w must have length nslast: the number of simplices on the last grid otherwise */
  unsigned int flag;
  if((w->row<=0) || !(w->val)) {
    /* for(i=0;i<ns;i++){ */
    /*   if(sc->gen[i] < level) continue; */
    /*   sc->marked[i]=TRUE; */
    /* } */
    return 2;
  }
  INT n=sc->n,n1=n+1,ns=sc->ns,nv=sc->nv,level=sc->level;
  INT i,nslast;
  INT ntotal=w->row;
  REAL *wval=w->val;
  nslast=0;
  REAL avew=0e0;
  // find the average of the error 
  for(i = 0;i<ns;i++){
    // we look at the elements on the last level:
    if(sc->childn[i] < 0) {
      //    if(sc->gen[i] < level) continue;
      avew+=fabs(wval[nslast]);
      nslast++;
      if(nslast > ntotal) { fprintf(stderr,"\n***ERROR: nslast > ntotal (%d > %d) in function %s;\n",nslast, ntotal, __FUNCTION__); return 1; }
    }
  }
  avew=avew/nslast;
  nslast=0;
  //mark for refinement every unrefined and unmarked simplex for which error > average error
  for(i = 0;i<ns;i++){
    //    if(sc->gen[i] < level) continue;    
    if(!(sc->marked[i]) && sc->childn[i] <0) {
      if(fabs(wval[nslast])>avew) {
	sc->marked[i]=TRUE;
      } else {
	sc->marked[i]=FALSE;
      }
      nslast++;
    }
  }
  return 0;
}
void markall(scomplex *sc){
  /* mark all simplices for refinement */
  INT ns=sc->ns,level=sc->level;
  INT i,nslast;
  /* mark everything and bail out. */
  for(i=0;i<ns;i++)
    sc->marked[i]=TRUE;
  return;
}
INT xins(INT n, INT *nodes, REAL *xs, REAL *xstar)
{
  /* 
     In dimension "n" constructs the map from reference simplex to
     simplex with coordinates xs[0..n].  Then solves a linear system
     with partial pivoting to determine if a point given with
     coordinates xstar[0..n-1] is in the (closed) simplex defined by
     "nodes[0..n] and xs[0..n]"
  */  
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
  //  fprintf(stdout,"\nSolution:\n");
  //  print_full_mat(n,1,xhat,"xhat");
  REAL xhatn=1e0,eps0=1e-10,xmax=1e0+eps0;
  /* check the solution if within bounds */
  INT flag = 0;
  for(j=0;j<n;j++){
    if((xhat[j] < -eps0) || (xhat[j] > xmax)){
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
