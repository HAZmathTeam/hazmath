/*! \file examples/amr_grids/amr_grids_supporting.h
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 2019/01/09.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program generates simplicial grids in 2,3,4... dimension.
 *
 * \note mock code for solve-estimate-mark and refinement. 
 */
/* #ifndef MAX_NODES_PER_SIMPLEX */
/* #define MAX_NODES_PER_SIMPLEX 10 */
/* #endif */
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
      //flag = isxnears(n, scnjn, sc->x, xstar0, threshold); 
      flag = xins(n,scnjn,sc->x,xstar0);
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
/**************************************************************************/
static ivector mark_around_pts(scomplex *sc, scomplex *scglobal, INT nstar, REAL *xstar, iCSRmat *node_ins, const INT max_nodes)
{
  /* scglobal is the global sc that contains all refinements */
  iCSRmat ins_node;
  dCSRmat xins_node,xnode_ins;
  INT *ii,*jj;
  INT dim=scglobal->n,dim1=dim+1,n=sc->n,n1=n+1,ns=sc->ns;
  INT nzw=-10,istar,jstar,jk,p,pj,cj0,cjn,j=-1,iwa,iwb,keep;
  /* un-mark everything */
  ivector marked=ivec_create(sc->ns);// this is only created on the top level....
  for(j=0;j<ns;j++) marked.val[j]=FALSE;
  //  for(j=0;j<ns;j++) marked.val[j]=TRUE;
  //  bail out if no points
  if(!nstar || (xstar==NULL)) return marked;
  INT *scnjn=NULL;
  REAL *xstar0=NULL;
  //  fprintf(stdout,"\nlevel=%d",scglobal->level);
  if(!(scglobal->level)){
    // on the coarsest level we need to set the matrix as for a first time:
    node_ins->row=nstar;
    node_ins->col=scglobal->ns;
    node_ins->IA=calloc((nstar+1),sizeof(INT));
    nzw=0;
    node_ins->IA[0]=nzw;
    for(jstar=0;jstar<nstar;++jstar){
      xstar0=xstar+jstar*dim;
      for(j=0;j<scglobal->ns;j++){
	scnjn = scglobal->nodes+j*dim1;
	// check if the point is inside the simplex.
	if(!xins(dim,scnjn,scglobal->x,xstar0)){
	  nzw++;
	}
      }  
      node_ins->IA[jstar+1]=nzw;
    }
    node_ins->nnz=node_ins->IA[nstar];
    // now run again, this time filling in the points:
    node_ins->JA=calloc(node_ins->nnz,sizeof(INT));
    nzw=0;
    for(jstar=0;jstar<nstar;++jstar){
      xstar0=xstar+jstar*dim;
      for(j=0;j<scglobal->ns;j++){
	scnjn = scglobal->nodes+j*dim1;
	// check if the point is inside the simplex.
	if(!xins(dim,scnjn,scglobal->x,xstar0)){	
	  node_ins->JA[nzw]=j;
	  marked.val[j]=TRUE;
	  nzw++;	
	}
      }  
    }
    //
  } else {
    //    haz_scomplex_print(scglobal,0,__FUNCTION__);
    ins_node=icsr_create(node_ins->col,node_ins->row,node_ins->nnz);
    xins_node=dcsr_create(0,0,0);
    xins_node.row=ins_node.row; xins_node.col=ins_node.col;  xins_node.nnz=ins_node.nnz;
    xins_node.IA=ins_node.IA;  xins_node.JA=ins_node.JA; xins_node.val=NULL;
    // input:
    xnode_ins=dcsr_create(0,0,0);
    xnode_ins.row=node_ins->row;  xnode_ins.col=node_ins->col;  xnode_ins.nnz=node_ins->nnz;
    xnode_ins.IA=node_ins->IA; xnode_ins.JA=node_ins->JA;  xnode_ins.val=NULL;
    // transpose
    dcsr_transz(&xnode_ins,NULL,&xins_node);
    //    fprintf(stdout,"\nNS=%d",scglobal->ns);fflush(stdout);
    //    icsr_write_icoo("ins_node.txt",&ins_node);
    /* INT nzw=0; */
    INT pjiter;
    nzw=0;
    // first run compute nonzeroes:
    for(j=0;j<scglobal->ns;j++){
      // we do the finest level only!
      if(scglobal->child0[j]<0||scglobal->childn[j]<0){
	pj=scglobal->parent[j];
	p=pj;
	if(pj<0) pj=j;
	//	fprintf(stdout,"\nj=%d,p=%d;pj=%d,ins_node.row=%d",j,p,pj,ins_node.row);fflush(stdout);
	pjiter=0;
	while(pj>=ins_node.row) {
	  pj=scglobal->parent[pj];
	  pjiter++;
	  if(pjiter>scglobal->level){
	  // stop w err if this is not working. 
	    fprintf(stderr,"%%%%****ERROR (%s): ancestor not found in refinement tree for element %lld on level=%lld\n%%%%****Exiting....",__FUNCTION__,(long long int)j,(long long int)scglobal->level);
	    exit(16);
	  }
	}
	/* if(pjiter>0){ */
	/*   fprintf(stdout,"\nNew:pj[%d]==%d,pjiter=%d",j,pj,pjiter);fflush(stdout); */
	/* } */
	iwa=ins_node.IA[pj];
	iwb=ins_node.IA[pj+1];
	if((iwb-iwa)<=0) {
	  //empty simplex do nothing
	  continue;
	}
	scnjn = scglobal->nodes+j*n1;
	for(jk=iwa;jk<iwb;jk++){
	  jstar=ins_node.JA[jk];
	  xstar0=xstar+jstar*n;
	  // check if the point is inside the simplex.
	  if(!xins(n,scnjn,scglobal->x,xstar0)){
	    // this is in, so we count it. 
	    nzw++;
	  } 
	}
      }
    }
    // reallocate as some nodes may belong to more than 1 simplex... otherwise this can be simplified.
    if(nzw>ins_node.nnz){
      ins_node.val=realloc(ins_node.val,nzw*sizeof(INT));
      ins_node.JA=realloc(ins_node.JA,nzw*sizeof(INT));
    }
    ins_node.nnz=nzw;
    if(nzw>node_ins->nnz){
      node_ins->JA=realloc(node_ins->JA,nzw*sizeof(INT));
    }
    node_ins->nnz=nzw;
    for(j=0;j<ins_node.nnz;++j){
      ins_node.val[j]=-1;
      node_ins->JA[j]=-1;
    }
    // now run again to get the JA in val.
    jj=node_ins->JA;
    ii=ins_node.val;
    nzw=0;
    for(j=0;j<scglobal->ns;j++){
      /* // check if we have marked and do nothing if we are marked: */
      /* if(marked.val[j]) {continue;} */
      // we do the finest level only!
      if(scglobal->child0[j]<0||scglobal->childn[j]<0){
	pj=scglobal->parent[j];
	if(pj<0) pj=j;
	pjiter=0;
	while(pj>=ins_node.row) {
	  pj=scglobal->parent[pj];
	  pjiter++;
	  if(pjiter>scglobal->level){
	  // stop w err if this is not working. 
	    fprintf(stderr,"%%%%****ERROR (%s): ancestor not found in refinement tree for element %lld on level=%lld\n%%%%****Exiting....",__FUNCTION__,(long long int)j,(long long int)scglobal->level);
	    exit(16);
	  }
	}
	iwa=ins_node.IA[pj];
	iwb=ins_node.IA[pj+1];
	if((iwb-iwa)<=0) {
	  //empty simplex do nothing
	  continue;
	}
	scnjn = scglobal->nodes+j*n1;
	for(jk=iwa;jk<iwb;jk++){
	  jstar=ins_node.JA[jk];
	  xstar0=xstar+jstar*n;
	  // check if the point is inside the simplex.
	  if(!xins(n,scnjn,scglobal->x,xstar0)){
	    // this is in, so we add it.
	    ii[nzw]=j; // note: ii is alias for ins_node.val
	    jj[nzw]=jstar; // note: jj is alias for node_ins->JA (node_ins->val does not exist).
	    nzw++;
	  } 
	}//end-for points in j-th simplex
      }//if on finest level
    }// endfor over all simplices
    //now we need to convert and clean up (use dcoo_2_dcsr for our particular case with tmp->val used as jj.
    int i, iind, jind;
    ins_node.row=scglobal->ns;
    ins_node.col=nstar;
    ins_node.nnz=nzw;
    INT *ind = calloc((ins_node.row+1),sizeof(INT));
    ins_node.IA = realloc(ins_node.IA,(ins_node.row+1)*sizeof(INT));
    // initins_node.IAlize
    memset(ind, 0, sizeof(INT)*(ins_node.row+1));
    // count number of nonzeros in each row
    for (i=0; i<nzw; ++i) ind[ii[i]+1]++;
    // set row pointer
    ins_node.IA[0] = 0;
    for (i=1; i<=ins_node.row; ++i) {
        ins_node.IA[i] = ins_node.IA[i-1]+ind[i];
        ind[i] = ins_node.IA[i];
    }

    // set column index and values
    for (i=0; i<nzw; ++i) {
        iind = ii[i];
        jind = ind[iind];
        ins_node.JA[jind] = jj[i];
        ind[iind] = ++jind;
    }
    free(ind);
    // now all should be ready we transpose and we are done:
    ins_node.nnz=nzw;
    xins_node.row=ins_node.row;xins_node.col=xins_node.col;xins_node.nnz=ins_node.nnz;
    node_ins->col=ins_node.row;//transpose should also have this, 
    node_ins->row=ins_node.col;
    node_ins->nnz=ins_node.nnz;
    // here we need to assign the right pointers:
    xins_node.IA=ins_node.IA;xins_node.JA=ins_node.JA;
    xnode_ins.IA=node_ins->IA;xnode_ins.JA=node_ins->JA;
    dcsr_transz(&xins_node,NULL,&xnode_ins);
    //    icsr_write_icoo("ins_node.txt",&ins_node);
    //    icsr_write_icoo("node_ins.txt",node_ins);
    for(j=0;j<sc->ns;j++){
      marked.val[j]=FALSE;
    }
    nzw=0; //working
    for(j=0;j<scglobal->ns;j++){
      if(scglobal->child0[j]<0||scglobal->childn[j]<0){
	nzw++;
	if((ins_node.IA[j+1]-ins_node.IA[j])>=max_nodes ){
	  /* if(scglobal->level>50){ */
	  /*   fprintf(stdout,"\nj=%d; nodes=%d (",j,ins_node.IA[j+1]-ins_node.IA[j]); */
	  /*   for(i=ins_node.IA[j];i<ins_node.IA[j+1];++i){ */
	  /*     fprintf(stdout,"%d ",ins_node.JA[i]); */
	  /*   } */
	  /*   fprintf(stdout,")\n"); */
	  /* } */
	  p=abs((scglobal->child0[j]+1));
	  //	  fprintf(stdout,"\np=%d",p);fflush(stdout);
	  marked.val[p]=TRUE;
	}
      }
    }
    icsr_free(&ins_node);
  }
  /* if(scglobal->level > 256){ */
  /*   fprintf(stdout,"\nnode_ins(row)=%d,node_ins(col)=%d,node_ins(nnz)=%d\n",node_ins->row,node_ins->col,node_ins->nnz);fflush(stdout); */
  /*   icsr_write_icoo("z123.txt",node_ins); */
  /* } */
  return marked;
}
/************************************************************/
static char *fname_set(const char *dir,const char *fname_in)
{
  // combines names: fname_in[0]=dir/fname_in[0]
  size_t ldir0=strlen(dir)+1;
  size_t lfname=strlen(fname_in);
  char *fname=strndup(dir,ldir0);
  fname=realloc(fname,(lfname+ldir0+1)*sizeof(char));
  strncat(fname,fname_in,lfname);
  return fname;
}
/************************************************************/
