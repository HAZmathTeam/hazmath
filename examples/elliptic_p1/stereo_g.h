/*! \file src/amr/stereographic.h
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note Just plotting, using stereographic projection for 2d,3d,4d
 *
 *  \note: modified by ltz on 20210602
 *
 */
//
// we first create a boundary simplicial complex which we will plot
// usiing stereographic projection
scomplex *sc_bndry(scomplex *sc)
{
  // finds the boundary simplicial complex of a simplicial complex. 
  INT ns = sc->ns,nv=sc->nv, dim=sc->n;
  scomplex *dsc=malloc(sizeof(scomplex));//boundary scomplex
  INT dim1=dim+1,iii,i,j,k,l,m,isn1,is,nnzbf,ns_b1,in1,jn1;
  INT ns_b=0;
  for(i=0;i<ns;i++){
    for(j=0;j<dim1;j++){
      if(sc->nbr[i*dim1+j]<0)
	ns_b++;
    }
  }
  dsc->ns=ns_b;
  // end computation of boundary faces. 
  INT *fnodes=calloc(ns_b*dim,sizeof(INT));
  ns_b1=0;
  for(i=0;i<ns;i++){
    for(j=0;j<dim1;j++){
      if(sc->nbr[i*dim1+j]<0) {
	k=0;
	for(m=0;m<dim1;m++){
	  if(m==j) continue;
	  fnodes[ns_b1*dim+k]=sc->nodes[i*dim1+m];
	  //	  fprintf(stdout,"\nnodes=%d",sc->nodes[i*dim1+m]);fflush(stdout);
	  k++;
	}
	ns_b1++;
      }
    }
  }
  if(ns_b!=ns_b1){
    fprintf(stderr,"\n%%***ERROR(1): num. bndry faces mismatch (ns_b=%d .ne. ns_b=%d) in %s",ns_b1,ns_b,__FUNCTION__);
    exit(65);
  }
  /* fprintf(stdout,"\nelnodes111=["); */
  /* for(i=0;i<ns_b;++i){ */
  /*   for(j=0;j<dim;++j){ */
  /*     fprintf(stdout,"%4i ",fnodes[dim*i+j]); */
  /*   } */
  /*   fprintf(stdout,";\n"); */
  /* } */
  /* fprintf(stdout,"]\n"); */
  //
  // FIX numbering, there is global numbering of nodes and local numbering of nodes:
  INT *indx=calloc(sc->nv,sizeof(INT));
  INT *indxinv=calloc(sc->nv,sizeof(INT));
  memset(indx,0,sc->nv*sizeof(INT));
  memset(indxinv,0,sc->nv*sizeof(INT));
  for(i=0;i<sc->nv;++i) indx[i]=-1;
  for(i=0;i<ns_b*dim;++i) indx[fnodes[i]]++;// for interior pts this does not get incremented
  INT nv_b=0;
  for(i=0;i<sc->nv;++i){
    if(indx[i]<0) continue;// interior pt
    indx[i]=nv_b;
    indxinv[nv_b]=i;
    nv_b++;
  }
  fprintf(stdout,"\n%%number of boundary vertices=%d (total nv=%d)\n",nv_b,sc->nv);
  /////////// init the scomplex;
  dsc=haz_scomplex_init((dim-1),ns_b,nv_b,dim);
  dsc->nbig=sc->n;
  dsc->n=(dsc->nbig-1);
  ////////////////
  if(dsc->nbig>dsc->n){
    fprintf(stdout,"\nSimplicial complex of dimension %d embedded in sc of dimension %d\n\n",dsc->n,dsc->nbig);
  }
  // there are things we dont need:
  free(dsc->nodes);
  dsc->nodes=fnodes;
  if(dsc->nv<sc->nv)
    indxinv=realloc(indxinv,dsc->nv*sizeof(INT));
  // now we can init the complex and then remove the unnecessary stuff:
  for(i=0;i<dsc->ns*(dsc->n+1);++i)
    fnodes[i]=indx[fnodes[i]];
  // set x, sc->bndry and so on:
  for(i=0;i<dsc->nv;i++){
    in1=i*dsc->nbig;
    j=indxinv[i];
    jn1=j*sc->n; // original coords:
    memcpy((dsc->x+in1),sc->x+jn1,sc->n*sizeof(REAL));
  }
  free(indx);
  free(indxinv);
  // we now have the numbering of simplices in dsc;
  // find neighbors:
  /* fprintf(stdout,"\nelnodes=["); */
  /* for(i=0;i<dsc->ns;++i){ */
  /*   for(j=0;j<(dsc->n+1);++j){ */
  /*     fprintf(stdout,"%4i ",dsc->nodes[(dsc->n+1)*i+j]+1); */
  /*   } */
  /*   fprintf(stdout,";\n"); */
  /* } */
  /* fprintf(stdout,"]\n"); */
  //
  find_nbr(dsc->ns,dsc->nv,dsc->n,dsc->nodes,dsc->nbr);
  //
  /* fprintf(stdout,"\nelnbr=["); */
  /* for(i=0;i<dsc->ns;++i){ */
  /*   //    fprintf(stdout,"\nelnbr[%d]=(",i); */
  /*   for(j=0;j<((dsc->n+1));++j){ */
  /*     fprintf(stdout,"%4i ",dsc->nbr[(dsc->n+1)*i+j]+1); */
  /*   } */
  /*   fprintf(stdout,";\n"); */
  /* } */
  /* fprintf(stdout,"]\n"); */
  return dsc;
}
//////////////////////////////////////////////////////////////////////////
INT stereo_g(scomplex *dsc)
{
// dsc should be a boundary complex.
  if(dsc->nbig != (dsc->n+1)){
    fprintf(stdout,"\nERR: Wrong dimensions of the boundary simplicial complex");
    return 1;
  }
  REAL tol=1e-10;
  INT i,j,node,n=dsc->n,n1=dsc->n+1,ns=dsc->ns,nv=dsc->nv;  
  INT nbig=dsc->nbig,nbig1=dsc->nbig+1;
  // mark for removal all vertices with last coordinate WITHIN TOL OF (1) and all simplices attached to them;
  INT *indv=calloc(nv,sizeof(INT));
  for(i=0;i<nv;i++) indv[i]=-1;
  REAL xn;
  INT nvnew=0;
  for(i=0;i<nv;i++){
    xn=dsc->x[nbig*i+nbig-1];
    if(fabs(xn-1e0)<tol) continue;
    indv[i]=nvnew;
    nvnew++;
  }
  INT nsnew=0,remove;
  for(i=0;i<ns;i++){
    remove=0;
    for(j=0;j<n1;j++){
      node=dsc->nodes[i*n1+j];
      if(indv[node]<0){
	remove=1;
	break;
      }
    }
    if(remove) continue;
    //    fprintf(stdout,"\nnews=%d: ",nsnew);fflush(stdout);
    for(j=0;j<n1;j++){
      node=dsc->nodes[i*n1+j];
      dsc->nodes[nsnew*n1+j]=indv[node];
      //      fprintf(stdout,"(%d->%d)",node,indv[node]);fflush(stdout);
    }
    dsc->flags[nsnew]=0;
    nsnew++;
  }
  for(i=0;i<nv;i++){
    node=indv[i];
    if(node>=0) {
      // copy coordinates after stereographic projection:
      xn=1e0/(1e0-dsc->x[nbig*i+nbig-1]);
      //      fprintf(stdout,"\ni=%d,node=%d,xnold=%e,xn=%e",i,node,dsc->x[nbig*i+nbig-1],xn);
      for(j=0;j<dsc->n;j++){
	dsc->x[node*n+j]=dsc->x[i*nbig+j]*xn;
      }
      dsc->bndry[node]=0;
    }
  }
  //  fprintf(stdout,"\n");
  free(indv);
  // adjust sizes;
  dsc->nbig=n;
  dsc->ns=nsnew;
  dsc->nv=nvnew;
  dsc->nodes=realloc(dsc->nodes,dsc->ns*(dsc->n+1)*sizeof(INT));
  dsc->x=realloc(dsc->x,(dsc->nv*dsc->n)*sizeof(REAL));
  fprintf(stdout,"\nStereographic projection: simplices=%d,vertices=%d\n\n",	\
	  dsc->ns,dsc->nv);fflush(stdout);
  //haz_scomplex_print(dsc,0,__FUNCTION__);
  return 0;
}
  
