/****************************************************************/
dvector *exmpl_solve(scomplex *sc, void *all)
{
  dvector *pts;/**/
  /* number of points, here used to define the solfem as this is only
     an example not involving solution of anything! */
  INT npts;
  npts=1;
  pts=(dvector *)malloc(sizeof(dvector));
  pts[0]=dvec_create(sc->n*npts);
  memset(pts->val,0,pts->row*sizeof(REAL));// this is set to be the origin;
  /*SET UP A SOLUTION:*/
  dvector *solfem=malloc(sizeof(dvector));
  solfem[0]=dvec_create(sc->ns);
  /* 
   *  Sets the solution equal to 1 at select simplices.
   */
  dvec_set_amr(1.,sc,pts,solfem->val);  
  /*free*/
  dvec_free(pts);
  /*SOLVE: as output we have the fem solution solfem as a dvector
    (value on every simplex).*/
  return solfem;
}
dvector *exmpl_estimate(scomplex *sc, dvector *solfem,void *all)
{
  /* if all is null there is nothing to estimate here */  
  /* extract the numerical and the "exact" solution and allocate space for the estimator */
  dvector *exact,*estimator;
  exact=malloc(1*sizeof(dvector));
  exact[0]=dvec_create(solfem->row);  
  estimator=malloc(1*sizeof(dvector));
  estimator[0]=dvec_create(sc->ns);  
  /* 
   * Compute the calues of the error estimator. In this example: the
   * estimator on every simplex equals to (the absolute value of the
   * difference between the exact solution and the numerical
   * solution).
   *
   * WARNING: THIS IS JUST AN EXAMPLE, in real applications the error
   *          estimator used could be much more involved.
  */
  INT j;
  for(j=0;j<sc->ns;j++)
    estimator->val[j]=fabs(exact->val[j]-solfem->val[j]);
  dvec_free(exact);// not needed any more;
  return estimator;
}
/**********************************************************************/
ivector *exmpl_mark(scomplex *sc,dvector *estimator,void *all)
{
  INT i,k=0;
  REAL errmax,gamma=0.75;
  ivector *marked=malloc(1*sizeof(ivector));
  marked[0]=ivec_create(sc->ns);
  errmax=estimator->val[0];
  for(i=1;i<sc->ns;i++)
    if(fabs(estimator->val[i])>errmax)
      errmax=estimator->val[i];
  for(i=0;i<sc->ns;i++)
    if(fabs(estimator->val[i])>gamma*errmax){
      marked->val[i]=TRUE;
      k++;
    }    else
      marked->val[i]=FALSE;
  //  fprintf(stdout,"\nnumber of simplices to be refined=%d\n",k);
  return marked;
}
/*************************************************************************/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
ivector *mark_near_points(scomplex *sc, INT nstar, REAL *xstar)
{
  /* 
     from marked simplices, remove any simplex that does not contain a
     point from xstar[...]. unmarked simplices are left unmarked 
  */
//  fprintf(stdout,"\nNSTAR=%d\n",nstar);fflush(stdout);
  INT n=sc->n,n1=n+1,ns=sc->ns;
  INT istar,jstar,j=-1;
  /* mark everything */
  ivector *marked=malloc(1*sizeof(ivector));
  marked[0]=ivec_create(sc->ns);
  for(j=0;j<ns;j++) marked->val[j]=TRUE;
  /* bail out if no points */
  if(!nstar || (xstar==NULL)) return marked;
  INT *scnjn=NULL,mrkno=0,mrktst=0,flagmrk=0;  
  REAL *xstar0=NULL;
  for(j=0;j<ns;j++){
    flagmrk=0;
    // check if we have been marked:
    if(!marked->val[j]) {continue;}
    mrktst++;
    scnjn = sc->nodes+j*n1;
    for(jstar=0;jstar<nstar;jstar++){
      //      fprintf(stdout,"\nel=%d\n",jstar+1);
      xstar0=xstar+jstar*n;
      // check if the point is inside the simplex. 
      if(!xins(n,scnjn,sc->x,xstar0)){
	//	fprintf(stdout,"\nel=%d, found: %d; NOmrk=%d",j,jstar+1,mrkno);
	mrkno--;
	marked->val[j]=TRUE;
	flagmrk=1;
	break;
      }
    }
    if(flagmrk) continue;
    mrkno++;
    marked->val[j]=FALSE;
  }
  return marked;
}
