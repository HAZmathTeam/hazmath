/*! \file examples/amr_grids/amr_grids_supporting.h
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 2019/01/09.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program generates simplicial grids in 2,3,4... dimension.
 *
 * \note mock code for solve-estimate-mark and refinement. 
 */
/**********************************************************************/
/*!
 * \fn INT isxnears(INT n, INT *splx_nodes, REAL *splx_coord, REAL *x, REAL threshold)
 *
 * \brief check if the point x is close to a given simplex (within the given distance)
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
INT isxnears(INT n, INT *splx_nodes, REAL *splx_coord, REAL *x, REAL threshold)
{

  //printf("dimension: %d\n", n);

  // local variables
  INT i,j, node_start;
  REAL *temp = (REAL *)calloc(n,sizeof(REAL));
  REAL distance = BIGREAL;

  // return variable
  INT flag = 0;

  //printf("loop to compute the barycenter\n");

  // compute the barycenter of the simplex
  for (i=0; i<n+1; i++) // loop over nodes
  {
      // get the start index for node i
      //printf("i=%d; node index = %d\n", i, splx_nodes[i]*n);
      node_start = splx_nodes[i]*n;

      // compute node[i] - x
      array_axpyz(n, -1.0, x, splx_coord+node_start, temp);

      // compute distance and keep the minimum one
      distance = MIN(distance, array_norm2(n,temp));

  }

  //printf("distance = %f\n", distance);

  // check
  if (distance > threshold)
  {
    flag = n+1;
  }

  return flag;
}


/****************************************************************/
dvector exmpl_solve(scomplex *sc, void *all)
{
  REAL *pts;/**/
  /* number of points, here used to define the solfem as this is only
     an example not involving solution of anything! */
  /**********************************************************************/
  INT npts;
  npts=1;
  pts=(REAL *)calloc(sc->n*npts,sizeof(REAL)); // npts x n vector of
					       // points used to mock
					       // marking
  memset(pts,0,npts*sizeof(REAL));// this is set to be the origin;
  /**********************************************************************/
  /*SET UP A SOLUTION:*/
  dvector solfem=dvec_create(sc->ns);
  /*
   *  Sets the solution equal to 1 at select simplices.
   */
  dvec_set_amr(1.,sc,npts,pts,solfem.val);
  /*free*/
  free(pts);
  /*SOLVE: as output we have the fem solution solfem as a dvector
    (value on every simplex).*/
  return solfem;
}
dvector exmpl_estimate(scomplex *sc, dvector *solfem,void *all)
{
  /* if all is null there is nothing to estimate here */
  /* extract the numerical and the "exact" solution and allocate space for the estimator */
  dvector exact,estimator;
  exact=dvec_create(solfem->row);
  estimator=dvec_create(sc->ns);
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
    estimator.val[j]=fabs(exact.val[j]-solfem->val[j]);
  dvec_free(&exact);// not needed any more;
  return estimator;
}
/**********************************************************************/
ivector exmpl_mark(scomplex *sc,dvector *estimator,void *all)
{
  INT i,k=0;
  REAL errmax,gamma=0.75;
  ivector marked=ivec_create(sc->ns);
  errmax=estimator->val[0];
  for(i=1;i<sc->ns;i++)
    if(fabs(estimator->val[i])>errmax)
      errmax=estimator->val[i];
  for(i=0;i<sc->ns;i++)
    if(fabs(estimator->val[i])>gamma*errmax){
      marked.val[i]=TRUE;
      k++;
    }    else
      marked.val[i]=FALSE;
  //  fprintf(stdout,"\nnumber of simplices to be refined=%d\n",k);
  return marked;
}
/*************************************************************************/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
ivector mark_near_points(scomplex *sc, INT nstar, REAL *xstar, REAL threshold)
{
  /*
     from marked simplices, remove any simplex that does not contain a
     point from xstar[...]. unmarked simplices are left unmarked
  */
//  fprintf(stdout,"\nNSTAR=%d\n",nstar);fflush(stdout);
  INT n=sc->n,n1=n+1,ns=sc->ns;
  INT istar,jstar,j=-1;
  /* mark everything */
  ivector marked=ivec_create(sc->ns);
  for(j=0;j<ns;j++) marked.val[j]=TRUE;
  /* bail out if no points */
  if(!nstar || (xstar==NULL)) return marked;
  INT *scnjn=NULL,mrkno=0,mrktst=0,flagmrk=0;
  REAL *xstar0=NULL;
  INT flag = FALSE;

  for(j=0;j<ns;j++){
    flagmrk=0;
    // check if we have been marked:
    if(!marked.val[j]) {continue;}
    mrktst++;
    scnjn = sc->nodes+j*n1;
    for(jstar=0;jstar<nstar;jstar++){
      //      fprintf(stdout,"\nel=%d\n",jstar+1);
      xstar0=xstar+jstar*n;
      // check if the point is inside the simplex.
      //if(!xins(n,scnjn,sc->x,xstar0)){
      flag = isxnears(n, scnjn, sc->x, xstar0, threshold); 
      //flag = xins(n,scnjn,sc->x,xstar0);
      //printf("flad = %d\n", flag);
      //sleep(5);
      if(!flag){
	    //fprintf(stdout,"\nel=%d, found: %d; NOmrk=%d",j,jstar+1,mrkno);
	    mrkno--;
	    marked.val[j]=TRUE;
	    flagmrk=1;
	    break;
      }
    }
    if(flagmrk) continue;
    mrkno++;
    marked.val[j]=FALSE;
  }
  return marked;
}
