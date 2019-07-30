/*! \file  marking.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *   \note routines to mark simplices for refinement
*/

#include "hazmath.h"
#include "grid_defs.h"
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
  INT jstar,j=-1;
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
  /* w is a vector with errors on the last level...*/
  if((w->row<=0) || !(w->val)) {
    return 2;
  }
  INT ns=sc->ns,i,nslast=0,ntotal=w->row;
  REAL *wval=w->val;
  REAL avew=0e0;
  // find the average of the error 
  for(i = 0;i<ns;i++){
    // we look at the elements on the last level:
    if(sc->childn[i] < 0 || sc->childn[i] < 0) {
      avew+=fabs(wval[i]);
      nslast++;
      if(nslast > ntotal) { fprintf(stderr,"\n***ERROR: nslast > ntotal (%d > %d) in function %s;\n",nslast, ntotal, __FUNCTION__); return 1; }
    }
  }
  avew=avew/nslast;
  fprintf(stdout,"\n%s: average=%e",__FUNCTION__,avew);
  nslast=0;
  for(i = 0;i<ns;i++){
    if((sc->gen[i]!=sc->level) && (sc->child0[i]>=0 && sc->childn[i]>=0)) continue;
    sc->marked[i]=0;
    /*    fprintf(stdout,"\n%s: BEFORE: simplex %d (level=%d,mark=%d, gen=%d, c0=%d)",__FUNCTION__,i,sc->level,sc->marked[i],sc->gen[i],sc->child0[i]);
*/
  }
  INT nmark=0;
  //mark for refinement every unrefined and unmarked simplex for which error > average error  
  for(i = 0;i<ns;i++){ 
    if(!(sc->marked[i]) && (sc->childn[i] <0||sc->childn[i]<0)) {
      if(fabs(wval[i])>avew) {
        fprintf(stdout,"\n%s: marking simplex %d (mark=%d, c0=%d, e=%e)", \
		__FUNCTION__,i,sc->marked[i],sc->child0[i],wval[nslast]);
	sc->marked[i]=1;
	nmark++;
      } else {
	sc->marked[i]=0;
      }
      nslast++;
    }
  }
  return 0;
}
void markall(scomplex *sc, const INT amark){
  /* 
     mark everything with "amark"; amark could be 0 or 1 or true or
     false. if 1 (true) then every simplex is marked for refinement. If 0 no simplex is marked for refinement. 
 */
  INT ns=sc->ns;
  INT i;
  for(i=0;i<ns;i++)
    sc->marked[i]=amark;
  return;
}
