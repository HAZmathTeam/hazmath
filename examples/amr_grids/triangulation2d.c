#include "hazmath.h"
void find_cc_bndry_cc(scomplex *sc)
{
  // find connected components of the simplex->simplex graph and the
  // number of connected components on the boundary. sets all boundary codes
  // on every connected component to be different. the arrays nodes
  // and nbr should be set.
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
  dfs00_(&ns,neib->IA, neib->JA,&sc->cc,iblk,jblk);
  // set up simplex flags= connected component number:
  for(i=0;i<sc->cc;++i){
    for(k=iblk[i];k<iblk[i+1];++k){
      j=jblk[k];
      sc->flags[j]=i+1;
    }
  }
  /* fprintf(stdout,"\n%%number of connected components in the bulk=%d\n",sc->cc); */
  /* fprintf(stdout,"\n%%number of boundary faces=%d (nnzbf=%d)\n",nbf,nnzbf); */
  ///////////////////////////////////////////////////////////////////////
  // now working on the boundary:
  jblk=realloc(jblk,(2*nbf+2)*sizeof(INT));
  iblk=jblk+nbf+1;
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
    fprintf(stderr,"\n%%***ERROR: num. bndry faces mismatch (nbf=%d .ne. nbfnew=%d) in %s",nbf,nbfnew,__FUNCTION__);
    exit(65);
  }
  /* fprintf(stdout,"\nelnodes=["); */
  /* for(i=0;i<nbf;++i){ */
  /*   //    fprintf(stdout,"\nelnodes[%d]=(",i); */
  /*   for(j=0;j<dim;++j){ */
  /*     fprintf(stdout,"%4i ",fnodes[dim*i+j]+1); */
  /*   } */
  /*   fprintf(stdout,";\n"); */
  /* } */
  /* fprintf(stdout,"]\n"); */
  find_nbr(nbf,nv,(dim-1),fnodes,fnbr);
  /* fprintf(stdout,"\nelnbr=["); */
  /* for(i=0;i<nbf;++i){ */
  /*   //    fprintf(stdout,"\nelnbr[%d]=(",i); */
  /*   for(j=0;j<dim;++j){ */
  /*     fprintf(stdout,"%4i ",fnbr[dim*i+j]+1); */
  /*   } */
  /*   fprintf(stdout,";\n"); */
  /* } */
  /* fprintf(stdout,"]\n"); */
  neib->IA=realloc(neib->IA,(nbf+1)*sizeof(INT));
  neib->JA=realloc(neib->JA,(nnzbf+nbf)*sizeof(INT));
  iii=0;
  neib->IA[0]=iii;
  for(i=0;i<nbf;i++){
    neib->JA[iii]=i;
    iii++;
    for(j=0;j<dim;++j){
      is=fnbr[i*dim+j];
      if(is>=0){
	//	fprintf(stdout,"\ni=%d,j=%d",i,is);
  	neib->JA[iii]=is;
  	iii++;
      }
    }
    neib->IA[i+1]=iii;
  }
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
  dfs00_(&nbf,neib->IA, neib->JA,&sc->bndry_cc,iblk,jblk);  
  icsr_free(neib);
  // set up bndry flags= connected component number;
  for(i=0;i<sc->bndry_cc;++i){
    for(k=iblk[i];k<iblk[i+1];++k){
      j=jblk[k];
      for(m=0;m<dim;m++){
	sc->bndry[fnodes[dim*j+m]]=i+1;
      }
    }
  }
  /* for(j=0;j<sc->nv;j++){ */
  /*   fprintf(stdout,"\ncode[%d]=%d",j,sc->bndry[j]); */
  /* } */
  fprintf(stdout,"%%number of connected components in the bulk=%d\n",sc->cc);
  //  fprintf(stdout,"%%number of boundary faces=%d (nnzbf=%d)\n",nbf,nnzbf);
  fprintf(stdout,"%%number of connected components on the boundary=%d\n",sc->bndry_cc);
}
//////////////////////////////////////////////////////////////////////////
INT main()
{
  // construction from a list of points :
  INT ncir=-20,nsq=-20;
  INT np=-20,neold=-20,ne=-20,dim,ii,jj,i,j,k;
  REAL pi=M_PI, hmax=1e0, shape_bound=0.125e0;
  REAL r, a, h, sx,sy, twopi=2e0*pi;
  //////////////////////////////////////////////////////////
  dim=2; /// dimension
  REAL *c=calloc(dim,sizeof(REAL)); // seed point which is inside
					 // a hole of the domain; for
					 // a circular "hole" this is
  					 // the center;
  FILE *fp=stdin; //fp=fopen("hole_grid.input","r");
  //read input:
  for(j=0;j<dim;j++){
    fscanf(fp,"%lg",(c+j));
  }
  fscanf(fp,"%lg",&r);// radius
  fscanf(fp,"%i",&ncir);//number of pts on circle
  fscanf(fp,"%lg",&sx); //length of square edge parallel to x
  fscanf(fp,"%lg",&sy); //length of square edge parallel to y
  fclose(fp);
  //end read input
  nsq=4;// points on the square
  np=ncir+nsq; //total number of points to describe the domain
  REAL *x_in=(REAL *)calloc(np*dim,sizeof(REAL));// points on the
						 // boundary of the
						 // domain.
  h=twopi/((REAL ) ncir);
  a=0.;
  // points(circle)
  for(i=0;i<ncir;++i){
    ii=dim*i;
    x_in[ii]  =c[0] + r*cos(a);
    x_in[ii+1]=c[1] + r*sin(a);
    a += h;
  }
  // we need half lengths of the sides of the rectangle
  sx=sx*0.5;
  sy=sy*0.5;
  // (1,1);(0,1);(0,0);(1,0)
  REAL *xx=x_in+ncir*dim;
  // Corners of the square centered at the center of the circle:
  // NE
  xx[0] = c[0]+sx; xx[1]= c[1]+sy;
  // NW
  xx[2] = c[0]-sx; xx[3]= c[1]+sy;
  // SW
  xx[4] = c[0]-sx; xx[5]= c[1]-sy;
  // SE
  xx[6] = c[0]+sx; xx[7]= c[1]-sy;
  //
  /* for(i=0;i<np;++i){ */
  /*   fprintf(stdout,"\nx[%d]=(",i); */
  /*   fprintf(stdout,"%.4f,%.4f)",x_in[dim*i],x_in[dim*i+1]); */
  /* } */
  /* fprintf(stdout,"\n"); */
  //edges(circle)
  neold=np; // as many edges as points (2 loops around circle and
  // around the square,
  // array for edge constraints:
  INT *edges=(INT *)calloc(2*neold,sizeof(INT));
  ne=0;
  for(i=0;i<(ncir-1);++i){
    edges[2*ne]=i;
    edges[2*ne+1]=i+1;
    ne++;
  }
  // ncir-1 edges are done, we add the last one:
  edges[2*ne]=0;
  edges[2*ne+1]=ncir-1;
  ne++;
  // edges(square):
  for(i=0;i<(nsq-1);++i){
    ii=ncir+i;
    edges[2*ne]=ii;
    edges[2*ne+1]=ii+1;
    ne++;
  }
  edges[2*ne]=ncir;
  edges[2*ne+1]=ncir+nsq-1;
  ne++;
  //  fprintf(stdout,"\nCheck:neold-ne=%d\nEdges:\n",neold-ne);
  shape_bound=0.125;
  REAL d,dxy;
  // first edge lenngth:
  ii=dim*edges[0];jj=dim*edges[1];
  d=(x_in[ii]-x_in[jj])*(x_in[ii]-x_in[jj])		\
    + (x_in[ii+1]-x_in[jj+1])*(x_in[ii+1]-x_in[jj+1]);
  hmax=sqrt(d);
  // compute the min edge length to put it as a constraint for the
  // edge length in the triangulation
  for(i=0;i<ne;++i){
    ii=dim*edges[2*i];
    jj=dim*edges[2*i+1];
    /* fprintf(stdout,"\ne(%d)=(%d,%d)=[(",i,ii/dim,jj/dim); */
    /* fprintf(stdout,"%.3f,%.3f)--(%.3f,%.3f)]",	\ */
    /* 	    x_in[ii],x_in[ii+1],		\ */
    /* 	    x_in[jj],x_in[jj+1]);	       */
    d=0e0;
    for(j=0;j<dim;++j){
      dxy=x_in[ii+j]-x_in[jj+j];
      d+=dxy*dxy;
    }
    if(hmax>sqrt(d)) hmax=sqrt(d);
  }
  //  fprintf(stdout,"\nhmax=%f\n\n",hmax);
  //scomplex *sc=(scomplex *)malloc(sizeof(scomplex));
  scomplex *sc=(scomplex *)haz_scomplex_init(dim,0,0);    
  //    exit(55);
  sc->n=dim;
  sc->nbig=dim;
  INT dim1=dim+1;
  sc->nv=0;
  sc->ns=0;
  hmax=0.5;
  c_do_tri_2d(x_in,edges,					\
	      c,						\
	      np, ne, dim,					\
	      shape_bound,					\
	      hmax,
	      &sc->nodes,&sc->x,	\
	      &sc->nv,&sc->ns);
  free(c);
  free(x_in);
  free(edges);
  //
  sc->nbr=(INT *)calloc(dim1*sc->ns,sizeof(INT));
  find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
  haz_scomplex_init_part(sc);
  find_cc_bndry_cc(sc);
  // add space at the end but 12 is one byte less so space is
  // overwritten with \0 as it should.
  char *fname=strndup("hole_XXX.yyy ",12*sizeof(char));
  size_t nbyt=(strlen(fname)+1)*sizeof(char);
  snprintf(fname,nbyt,"hole%03i.haz",ncir);
  //  fprintf(stdout,"\n\nFNAME=|%s|\n",fname);
  /// output-haz
  hazw(fname,sc,0);
  /// output-vtu
  snprintf(fname,nbyt,"hole%03i.vtu",ncir);
  vtkw(fname,sc,0,1.);
  /*FREE*/
  haz_scomplex_free(sc);
  free(fname);
  return 0;
}
