// should be included after "hazmath.h"
void find_cc_bndry_cc(scomplex *sc)
{
  // find connected components of the simplex->simplex graph and the
  // number of connected components on the boundary. sets all boundary codes
  // on every connected component to be different. the arrays nodes[]
  // and nbr[] should be set.
  //
  INT ns = sc->ns,nv=sc->nv, dim=sc->n;
  INT dim1=dim+1,iii,i,j,k,l,m,isn1,is,nbf,nnzbf;
  iCSRmat *neib=malloc(1*sizeof(iCSRmat));
  neib[0]=icsr_create(ns,ns,dim1*ns+ns);
  INT *ibnd=neib->val;// use as working space
  nbf=0;
  nnzbf=0;
  iii=0;
  neib->IA[0]=iii;
  for(i=0;i<ns;i++) ibnd[i]=-1;
  for(i=0;i<ns;i++){
    neib->JA[iii]=i;
    iii++;
    isn1=i*dim1;
    for(j=0;j<dim1;j++){
      is=sc->nbr[isn1+j];
      if(is>=0){
	neib->JA[iii]=is;
	iii++;
      } else {
	nbf++;
	nnzbf+=dim; 
      }
    }
    neib->IA[i+1]=iii;
  }
  INT *jblk=calloc(2*ns+2,sizeof(INT));
  INT *iblk=jblk+ns+1;
  sc->cc=-10;
  iCSRmat *blk_dfs=run_dfs(ns,neib->IA, neib->JA);
  sc->cc=blk_dfs->row;
  iblk=blk_dfs->IA;
  jblk=blk_dfs->JA;
  //  dfs00_(&ns,neib->IA, neib->JA,&sc->cc,iblk,jblk);
  // set up simplex flags= connected component number:
  for(i=0;i<sc->cc;++i){
    for(k=iblk[i];k<iblk[i+1];++k){
      j=jblk[k];
      sc->flags[j]=i+1;
    }
  }
  fprintf(stdout,"\n%%number of connected components in the bulk=%d\n",sc->cc);
  fprintf(stdout,"\n%%number of boundary faces=%d (nnzbf=%d)\n",nbf,nnzbf);
  ///////////////////////////////////////////////////////////////////////
  // now working on the boundary:
  jblk=realloc(jblk,(2*nbf+2)*sizeof(INT));
  iblk=jblk + nbf + 1;
  INT *fnodes=calloc(nbf*dim,sizeof(INT));
  INT *fnbr=calloc(nbf*dim,sizeof(INT));
  INT nbfnew=0;
  for(i=0;i<ns;i++){
    for(j=0;j<dim1;j++){
      if(sc->nbr[i*dim1+j]<0) {
	k=0;
	for(m=0;m<dim1;m++){
	  if(m==j) continue;
	  fnodes[nbfnew*dim+k]=sc->nodes[i*dim1+m];
	  k++;
	}
	nbfnew++;
      }
    }
  }
  if(nbf!=nbfnew){
    fprintf(stderr,"\n%%***ERROR(1): num. bndry faces mismatch (nbf=%d .ne. nbfnew=%d) in %s",nbf,nbfnew,__FUNCTION__);
    exit(65);
  }
  // FIX numbering:
  INT *indx=calloc(2*sc->nv,sizeof(INT));
  INT *indxinv=indx+sc->nv;
  memset(indx,0,sc->nv*sizeof(INT));
  for(i=0;i<sc->nv;++i) indx[i]=-1;
  for(i=0;i<nbf*dim;++i) indx[fnodes[i]]++;
  INT nvbnd=0;
  for(i=0;i<sc->nv;++i){
    if(indx[i]<0) continue;
    indx[i]=nvbnd;
    indxinv[nvbnd]=i;
    nvbnd++;
  }
  fprintf(stdout,"\n%%number of boundary vertices=%d (total nv=%d)\n",nvbnd,sc->nv);
  if(nvbnd<sc->nv)
    indx=realloc(indx,(sc->nv+nvbnd)*sizeof(INT));
  for(i=0;i<nbf*dim;++i)
    fnodes[i]=indx[fnodes[i]];
  // end fix numbering
  /* fprintf(stdout,"\nelnodes=["); */
  /* for(i=0;i<nbf;++i){ */
  /*   //    fprintf(stdout,"\nelnodes[%d]=(",i); */
  /*   for(j=0;j<dim;++j){ */
  /*     fprintf(stdout,"%4i ",fnodes[dim*i+j]+1); */
  /*   } */
  /*   fprintf(stdout,";\n"); */
  /* } */
  /* fprintf(stdout,"]\n"); */
  find_nbr(nbf,nvbnd,(dim-1),fnodes,fnbr);
  /* fprintf(stdout,"\nelnbr=["); */
  /* for(i=0;i<nbf;++i){ */
  /*   //    fprintf(stdout,"\nelnbr[%d]=(",i); */
  /*   for(j=0;j<dim;++j){ */
  /*     fprintf(stdout,"%4i ",fnbr[dim*i+j]+1); */
  /*   } */
  /*   fprintf(stdout,";\n"); */
  /* } */
  /* fprintf(stdout,"]\n"); */
  /* fprintf(stdout,"\n\nnbf=%d,nnzbf=%d;nnzbf+nbf=%d\n\n",nbf,nnzbf,nnzbf+nbf);fflush(stdout); */
  neib->IA=realloc(neib->IA,(nbf+1)*sizeof(INT));
  neib->JA=realloc(neib->JA,(nnzbf+nbf)*sizeof(INT));
  neib->row=nbf;
  neib->col=nbf;
  iii=0;
  neib->IA[0]=iii;
  for(i=0;i<nbf;i++){
    neib->JA[iii]=i;
    iii++;
    for(j=0;j<dim;++j){
      is=fnbr[i*dim+j];
      if(is>=0){
	//	fprintf(stdout,"\ni=%d,j=%d,nnz=%d",i,is,iii);fflush(stdout);
  	neib->JA[iii]=is;
  	iii++;
      }
    }
    neib->IA[i+1]=iii;
  }
  //  fprintf(stdout,"\n\nrows=%d,nnzOOO=%d;nnzXXX=%d\n\n",neib->row,iii,neib->IA[neib->row]);fflush(stdout);
  neib->nnz=neib->IA[neib->row];
  //  icsr_print_matlab(stdout,neib);
  /* fprintf(stdout,"\nnbr00=["); */
  /* for(i=0;i<nbf;i++){ */
  /*   for(k=neib->IA[i];k<neib->IA[i+1];++k){ */
  /*     j=neib->JA[k]; */
  /*     fprintf(stdout,"%d %d %d\n",i+1,j+1,1); */
  /*   } */
  /* } */
  /* fprintf(stdout,"];\n"); */
  sc->bndry=(INT *)calloc(sc->nv,sizeof(INT));
  sc->bndry_cc=-10;
  //  dfs00_(&nbf,neib->IA, neib->JA,&sc->bndry_cc,iblk,jblk);  
  icsr_free(blk_dfs); blk_dfs=NULL;
  blk_dfs=run_dfs(nbf,neib->IA, neib->JA);
  sc->bndry_cc=blk_dfs->row;
  iblk=blk_dfs->IA;
  jblk=blk_dfs->JA;
  icsr_free(neib);
  // set up bndry flags= connected component number;
  for(i=0;i<sc->bndry_cc;++i){
    for(k=iblk[i];k<iblk[i+1];++k){
      j=jblk[k];
      for(m=0;m<dim;m++){
	sc->bndry[indxinv[fnodes[dim*j+m]]]=i+1;
      }
    }
  }
  icsr_free(blk_dfs);
  /* for(j=0;j<sc->nv;j++){ */
  /*   fprintf(stdout,"\ncode[%d]=%d",j,sc->bndry[j]); */
  /* } */
  fprintf(stdout,"%%number of connected components in the bulk=%d\n",sc->cc);
  //  fprintf(stdout,"%%number of boundary faces=%d (nnzbf=%d)\n",nbf,nnzbf);
  fprintf(stdout,"%%number of connected components on the boundary=%d\n",sc->bndry_cc);
  free(fnodes);
}
//////////////////////////////////////////////////////////////////////////
