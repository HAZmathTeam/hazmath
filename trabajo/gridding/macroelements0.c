/*! \file src/amr/macroelements.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note all files related to constructing mesh of macroelements and
 *  splitting it into simplices
 *
 */
#include "hazmath.h"
/************************************************************************/
macrocomplex *set_mmesh(input_grid *g0,					\
		cube2simp *c2s,						\
		INT *wrk)
{
  /*prepare macro element mesh (mmesh) for passing to the mesh generator*/
  /*wrk is working integer array should have at least size 2*nvcube+2 */
  INT i,j0,j1,kel,jel,ke;
  INT nvcube=c2s->nvcube,nvface=c2s->nvface;
  INT nel0,je,kj,k2,iel2v,jel2v,k1,kface,kbnd,found;
  INT *p=wrk; 
  INT *mnodes=p+nvcube+1;
  /*macro complex creation:)*/
  macrocomplex *mc=malloc(1*sizeof(macrocomplex));
  mc->nel=nel0=g0->nel; //important to set;
  mc->cc=mc->bndry_cc=1;
  mc->nf=mc->nfi=mc->nfb=-1;
  /*macro complex allocation*/
  /********************************************************************/
  mc->flags=calloc(mc->nel,sizeof(INT));
  /**/
  mc->nd=NULL;//later.
  mc->elneib=calloc(mc->nel,sizeof(INT *));
  mc->el2fnum=calloc(mc->nel,sizeof(INT *)); 
  mc->iindex=calloc(mc->nel,sizeof(INT *)); 
  for(i=0;i<g0->nel;i++){
    mc->elneib[i]=calloc(c2s->nf,sizeof(INT)); /* to hold the neighbors */
    mc->el2fnum[i]=calloc(c2s->nf,sizeof(INT)); /* to hold the face
						   numbers for every
						   element */
    mc->iindex[i]=NULL;//later, when nd is known;
  }
  /**/
  INT **elneib=mc->elneib;
  INT **el2fnum=mc->el2fnum;
  mc->etree=calloc((mc->nel+1),sizeof(INT));
  INT *etree=mc->etree;
  mc->bcodesf=NULL; /* later when num faces is known*/
  mc->isbface=NULL; /* later when num faces is known*/
  mc->bfs=NULL; /* later when bfs is called*/
  mc->fullel2el=NULL; /* later when el2el is formed*/
  /********************************************************************/
  /* let us also allocate all the rest here */
  /*  print_full_mat_int(g0->nel,(c2s->nvcube+1),g0->mnodes,"ex");*/
  ilexsort(g0->nf, (c2s->nvface+1),g0->mfaces,p);
  ilexsort(g0->nel,(c2s->nvcube+1),g0->mnodes,p);
  /*-------------------------------------------------------------------*/
  iCSRmat *el2v=malloc(1*sizeof(iCSRmat));
  el2v[0]=icsr_create(g0->nel,g0->nv,nvcube*g0->nel);
  el2v->IA[0]=0;
  for(kel=0;kel<nel0;kel++){
    memcpy((el2v->JA+el2v->IA[kel]),(g0->mnodes+kel*(nvcube+1)),nvcube*sizeof(INT));
    el2v->IA[kel+1]=el2v->IA[kel]+nvcube;
  }
  for(i=0;i<el2v->IA[nel0];i++)
    el2v->val[i]=1;
  iCSRmat *v2el=malloc(1*sizeof(iCSRmat));
  icsr_trans(el2v,v2el);	
  iCSRmat *el2el=malloc(1*sizeof(iCSRmat));
  icsr_mxm(el2v,v2el,el2el);
  icsr_free(v2el);
  /* create fullel2el for later and work here with the copy */
  iCSRmat *fullel2el=malloc(1*sizeof(iCSRmat));
  fullel2el[0]=icsr_create(el2el->row,el2el->col,el2el->nnz);  
  memcpy(fullel2el->IA,el2el->IA,(fullel2el->row+1)*sizeof(INT));
  memcpy(fullel2el->JA,el2el->JA,(fullel2el->nnz)*sizeof(INT));
  memcpy(fullel2el->val,el2el->val,(fullel2el->nnz)*sizeof(INT));
  /***********************************************************/  
  icsr_tri(fullel2el,'l');
  icsr_nodiag(fullel2el);
  /***********************************************************/  
  mc->fullel2el=fullel2el;
  /***********************************************************/  
  /*   shrink the el2el matrix by removing any entry with value not
       equal to nvface; */
  el2el->nnz=el2el->IA[0];
  for(kel=0;kel<nel0;kel++){    
    j0=el2el->IA[kel];
    j1=el2el->IA[kel+1];
    el2el->IA[kel]=el2el->nnz;
    for(ke=j0;ke<j1;ke++){
      jel=el2el->JA[ke];
      if(jel==kel) {
	el2el->JA[el2el->nnz]=kel;
	el2el->val[el2el->nnz]=1;
	el2el->nnz++;
	continue;
      }
      if(el2el->val[ke]!=nvface) continue;
      el2el->JA[el2el->nnz]=jel;
      el2el->val[el2el->nnz]=el2el->val[ke];
      el2el->nnz++;
    }
  }  
  el2el->IA[nel0]=el2el->nnz;
  el2el->JA=realloc(el2el->JA,el2el->nnz*sizeof(INT));
  el2el->val=realloc(el2el->val,el2el->nnz*sizeof(INT));
  /* find the connected components:*/
  INT *iblk=calloc((nel0+1),sizeof(INT));
  INT *jblk=calloc((nel0),sizeof(INT));
  // number of connected domains and boundaries;
  dfs00_(&nel0,el2el->IA, el2el->JA,&mc->cc,iblk,jblk);
  iblk=realloc(iblk,(mc->cc+1)*sizeof(INT));
  /*  prepare for bfs: remove diagonal in el2el */
  icsr_nodiag(el2el);
  mc->bfs=bfscc(mc->cc,iblk,jblk,el2el,etree);
  /* construct element to face to element maps for macroelements */
  /*FACES******************************************************/
  mc->nf=g0->nel*c2s->nf-(el2el->nnz/2);
  mc->nfi=(INT )(el2el->nnz/2);
  mc->nfb=mc->nf-mc->nfi;
  // face to vertex:
  iCSRmat *f2v=malloc(sizeof(iCSRmat));
  f2v[0]=icsr_create(mc->nf,g0->nv,nvface*mc->nf);
  /**/
  mc->bcodesf=calloc(mc->nf,sizeof(INT));
  mc->isbface=calloc(mc->nf,sizeof(INT));
  INT *bcodesf=mc->bcodesf;
  INT *isbface=mc->isbface;
  /*init*/
  /*-------------------------------------------------------------------*/
  for(kel=0;kel<g0->nel;kel++){
    for(i=0;i<c2s->nf;i++){
      elneib[kel][i]=-1;
      el2fnum[kel][i]=-1;
    }
  }
  /*-------------------------------------------------------------------*/
  f2v->IA[0]=0;
  for(i=0;i<mc->nf;i++){
    f2v->IA[i+1]=f2v->IA[i]+nvface;
    isbface[i]=0;
  }
  INT *facei=calloc(2*nvface,sizeof(INT));
  INT *facej=facei+nvface;
  kface=0;
  for(kel=0;kel<nel0;kel++){
    iel2v=el2v->IA[kel];
    for (ke=0;ke<c2s->nf;ke++){
      found=0;
      k1=ke*nvface;
      for(i=0;i<nvface;i++){
	facei[i]=el2v->JA[iel2v+c2s->faces[k1+i]];
      }
      j0=el2el->IA[kel];
      j1=el2el->IA[kel+1];
      for(kj=j0;kj<j1;kj++){
	jel=el2el->JA[kj];
	jel2v=el2v->IA[jel];
	for (je=0;je<c2s->nf;je++){
	  k2=je*nvface;
	  for(i=0;i<nvface;i++){
	    facej[i]=el2v->JA[jel2v+c2s->faces[k2+i]];
	  }
	  if(aresame(facei,facej,nvface)){
	    if(kel<jel){
	      /* fprintf(stdout,"\nptr=%d;FOUND: fi(ei),fj(ej)=%d(%d),%d(%d);",kface,ke,kel,je,jel); */
	      /* print_full_mat_int(1,nvface,facei,"facei"); */
	      /* print_full_mat_int(1,nvface,facej,"facej"); */
	      /* fprintf(stdout,"\nadd nonzero at:%d ",f2v->IA[kface]); */
	      memcpy((f2v->JA+f2v->IA[kface]),facei,nvface*sizeof(INT));
	      elneib[kel][ke]=jel;
	      elneib[jel][je]=kel;
	      el2fnum[kel][ke]=kface;
	      el2fnum[jel][je]=kface;
	      kface++;// pointer for the face2 vertex matrix
	    }
	    found=1;
	  }
	}
      }
      if(!found){
	/* fprintf(stdout,"\nptr=%d;NOT FOUND: fi(ei)=%d(%d);",kface,ke,kel); */
	/* print_full_mat_int(1,nvface,facei,"facei"); */
	/* fprintf(stdout,"\nadd nonzero at:%d",f2v->IA[kface]); */
	isbface[kface]=1;
	memcpy((f2v->JA+f2v->IA[kface]),facei,nvface*sizeof(INT));
/* 
in this way bcodesf[1:elneib[kel][ke]] gives us the code of the corresponding face;; this can also be linked to f2v and so on. 
*/
	el2fnum[kel][ke]=kface;
	kface++;
      }     
    }
  }
  /*****************************************************/
  for(i=0;i<f2v->IA[mc->nf];i++)
    f2v->val[i]=1;
  /*******************************************************************/    
  //  fprintf(stdout,"\nkface=%d  ? = ? %d\n",kface,f2v->nnz);
  /*ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ*/
  INT cfbig=((INT )MARKER_BOUNDARY_NO)+1;
  for(i=0;i<mc->nf;i++){
    bcodesf[i]=-cfbig;
    j1=f2v->IA[i+1]-f2v->IA[i];
    //    fprintf(stdout,"\nnnz(%d)=%d",i,j1); 
    if(j1>0){
      memcpy(facei,(f2v->JA+f2v->IA[i]),j1*sizeof(INT));
      for(je=0;je<g0->nf;je++){
	memcpy(facej,(g0->mfaces+je*(nvface+1)),nvface*sizeof(INT));
	kbnd=g0->mfaces[je*(nvface+1)+nvface];
	if(aresame(facei,facej,nvface)){
	  // if the face was in the list, take its code.
	  bcodesf[i]=kbnd;
	  //	  fprintf(stdout,"\ncode:%d ",bcodesf[i]);
	  //	  print_full_mat_int(1,nvface,facei,"facei");
	  //	  print_full_mat_int(1,nvface,facej,"facej");
	  break;
	}
      }
    }
  }
  // set all interior faces faces with no code to 0 code and all
  // boundary faces with no code to dirichlet code 1.
  for(i=0;i<mc->nf;i++){
    if(bcodesf[i]<(1-cfbig)){
      if(isbface[i])bcodesf[i]=1;
      else bcodesf[i]=0;
    }
  }
  for(i=0;i<mc->nf;i++){
    if(isbface[i] && (bcodesf[i]==0))bcodesf[i]=1;
    //    fprintf(stdout,"\n[%d]=%d",i,bcodesf[i]);
  }
  /***********************************************************************/
  /*Connected components on the boundaries. First find face to face map*/
  iCSRmat *v2f=malloc(1*sizeof(iCSRmat));
  icsr_trans(f2v,v2f);
  iCSRmat *f2f=malloc(1*sizeof(iCSRmat));
  /*******************************************************************/    
  icsr_mxm(f2v,v2f,f2f);
  /* 
     now remove all rows in f2f that correspond to interior faces and
     all entries that are not nvface/2, i.e. the number of vertices in
     a (n-2) cube;
  */
  f2f->nnz=f2f->IA[0];
  for(i=0;i<f2f->row;i++){    
    j0=f2f->IA[i];
    j1=f2f->IA[i+1];
    f2f->IA[i]=f2f->nnz;
    //    if(isbface[i]==0) continue;
    for(ke=j0;ke<j1;ke++){
      je=f2f->JA[ke];
      if(je==i){
	f2f->JA[f2f->nnz]=i;
	f2f->val[f2f->nnz]=0;
	f2f->nnz++;
	continue;
      }
      if((isbface[je]==0)||(f2f->val[ke]!=(nvface/2))) continue;
      //      if((je==i)||(isbface[je]==0)||(f2f->val[ke]!=(nvface/2))) continue;
      f2f->JA[f2f->nnz]=je;
      f2f->val[f2f->nnz]=f2f->val[ke];
      f2f->nnz++;
    }
  }    
  f2f->IA[f2f->row]=f2f->nnz;
  f2f->JA=realloc(f2f->JA,f2f->nnz*sizeof(INT));
  f2f->val=realloc(f2f->val,f2f->nnz*sizeof(INT));
  /*******************************************************************/    
  icsr_free(v2f);
  /*connected comps on the boundary*/
  iblk=realloc(iblk,(mc->nf+1)*sizeof(INT));
  jblk=realloc(jblk,(mc->nf)*sizeof(INT));
  dfs00_(&mc->nf,f2f->IA, f2f->JA,&mc->bndry_cc,iblk,jblk);
  mc->bndry_cc-=mc->nfi;
  //icsr_nodiag(f2f);
  icsr_free(f2f);
  free(iblk); 
  free(jblk);   
  /* /\*****************************************************\/     */
  /* fprintf(stdout,"\n");  */
  /* /\*****************************************************\/ */
  icsr_free(el2v);
  icsr_free(f2v);
  icsr_free(el2el);
  mc->elneib=elneib;
  mc->el2fnum=el2fnum;
  return mc;
}
/*******************************************************************/
/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
void scomplex_merge1(const INT nvall,		\
		     const INT nsall,		\
		     macrocomplex *mc,		\
		     scomplex **sc0,		\
		     cube2simp *c2s)
{
  /* combines an array of simplicial complexes constructed from a
     macro complex together. the mc->iindex array should contain all
     vertex numbers without repetitions
*/
  if(mc->nel==1) return;
  scomplex *sc=sc0[0];
  INT n1=(sc->n+1),ns0,nv=nvall,ns=nsall;
  INT kel,i,ii,j,in1,iin1,newv;
  sc->marked=realloc(sc->marked,ns*sizeof(INT));
  sc->gen=realloc(sc->gen,ns*sizeof(INT));
  sc->nbr=realloc(sc->nbr,ns*n1*sizeof(INT));
  sc->parent=realloc(sc->parent,ns*sizeof(INT));
  sc->child0=realloc(sc->child0,ns*sizeof(INT));
  sc->childn=realloc(sc->childn,ns*sizeof(INT));
  sc->nodes=realloc(sc->nodes,ns*n1*sizeof(INT));
  sc->bndry=realloc(sc->bndry,nv*sizeof(INT));
  sc->csys=realloc(sc->csys,nv*sizeof(INT));/* coord sys: 1 is polar, 2
					    is cyl and so on */
  /*connected components*/
  sc->cc=mc->cc;sc->bndry_cc=mc->bndry_cc;
  sc->flags=(INT *)realloc(sc->flags,ns*sizeof(INT));
  sc->x=(REAL *)realloc(sc->x,nv*(sc->n)*sizeof(REAL));
  sc->vols=(REAL *)realloc(sc->vols,ns*sizeof(REAL));
  sc->fval=(REAL *)realloc(sc->fval,nv*sizeof(REAL)); // function values at every vertex; not used in general;
  fprintf(stdout,"\nnsall=%d,nvall=%d",nsall,nvall);fflush(stdout);
  for(kel=1;kel<mc->nel;kel++){
    //    fprintf(stdout,"\n*********YYYYYYYYYY nv[%d]=%d\n",kel,sc0[kel]->nv);
    //    haz_scomplex_print(sc,0,"Z"); fflush(stdout);
    ns0=sc->ns;
    for (ii = 0;ii<sc0[kel]->ns;ii++) {
      i=ii+ns0;
      sc->marked[i] = sc0[kel]->marked[ii];
      sc->gen[i] = sc0[kel]->gen[ii];
      sc->parent[i]=sc0[kel]->parent[ii];
      sc->child0[i]=sc0[kel]->child0[ii];
      sc->childn[i]=sc0[kel]->childn[ii];
      sc->flags[i]=sc0[kel]->flags[ii];
      sc->vols[i]=sc0[kel]->vols[ii];
      in1=i*n1;
      iin1=ii*n1;
      for(j=0;j<n1;j++){
	newv=mc->iindex[kel][sc0[kel]->nodes[iin1+j]];
	sc->nodes[in1+j]=newv;
	sc->nbr[in1+j]=sc0[kel]->nbr[iin1+j]+ns0;
      }
    }
    for (ii = 0;ii<sc0[kel]->nv;ii++) {
      i=mc->iindex[kel][ii];
      //      fprintf(stdout,"\n*********YYYYYYYYYY nv[%d]=%d:::%d-->%d\n",kel,sc0[kel]->nv,ii,i);
      sc->bndry[i]=sc0[kel]->bndry[ii];
      sc->csys[i]=sc0[kel]->csys[ii];
      sc->fval[i]=sc0[kel]->fval[ii];
      in1=i*sc->n;
      iin1=ii*sc->n;
      for(j=0;j<sc->n;j++)
	sc->x[in1+j]=sc0[kel]->x[iin1+j];
    }
    sc->ns+=sc0[kel]->ns;
    haz_scomplex_free(sc0[kel]);
  }
  sc->nv=nvall;
  //  fprintf(stdout,"\nsc->nv=%d,sc->ns=%d; nvall=%d,nsall=%d\n",sc->nv,sc->ns,nvall,nsall);fflush(stdout);
  return;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
INT locate1(INT *b,
	    INT *a,INT n,			\
	    INT *a2,INT n2,INT m2)
{  
  /* 
     locates the elements of an array a in a two dimensional array
     a2. a has n elements and a2 is n2 by m2.  returns the the number
     of the rows(a2) that contain ALL elements of a; and b contains
     the indices of these rows; in case there is an element of a not
     contained in any row returns 0; in case no row contains ALL
     elements, returns 0; b should have size n2*sizeof(INT);
  */
  INT i,j,im,nb,aj,bi;
  for(i=0;i<n2;i++) b[i]=i;
  for(j=0;j<n;j++){    
    aj=a[j];
    nb=0;
    //    fprintf(stdout,"\na[%d]=%d; in:  ",j,aj); 
    for(i=0;i<n2;i++){
      if(b[i]<0) continue;
      im=i*m2;
      bi=locate0(aj,(a2+im),m2);
      if(bi<0){b[i]=bi;}
      else {b[i]=i;nb++;}
      //      fprintf(stdout,"\nnb=%d,n=%d bi=%d: ",nb,n,bi); 
      //      print_full_mat_int(1,m2,(a2+im),"a2");
      //      print_full_mat_int(1,n2,b,"b2");
    }
    if(!nb) return 0;// one of the elements was not found
  }
  //  fprintf(stdout,"\nnb=%d,n=%d: ",nb,n); 
  nb=0;
  for(i=0;i<n2;i++){
    if(b[i]<0) continue;
    b[nb]=b[i];
    //    fprintf(stdout,"XZXZ=b[%d]=%d;",nb,b[nb]); 
    nb++;
    b[i]=-1;
  }
  //  fprintf(stdout,"\n");
  return nb;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
void macrocomplex_free(macrocomplex *mc)
{
  INT i;
  for(i=0;i<mc->nel;i++){
    free(mc->nd[i]);
    free(mc->elneib[i]);
    free(mc->el2fnum[i]);
    free(mc->iindex[i]);
  }
  free(mc->nd);
  free(mc->elneib);
  free(mc->el2fnum);
  free(mc->iindex);
  free(mc->isbface);
  free(mc->bcodesf);
  free(mc->flags);
  free(mc->etree);
  icsr_free(mc->bfs);
  icsr_free(mc->fullel2el);
  return;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
void fix_grid(macrocomplex *mc,		\
	       scomplex **scin,			\
	       cube2simp *c2s,			\
	       input_grid *g0)
{
  /****************************************************************/
  /* 
     this piece of code removes all repeated vertices in macro
     elements making the numbering of the  vertices a valid numbering;
  */
  if(mc->nel<=1) return;
  scomplex *scp;
  INT dim=c2s->n,dim1=c2s->n+1,nvface=c2s->nvface,nvcube=c2s->nvcube;
// at most dim (n-1)dimensional faces may intersect to form a vertex
  INT *mp = (INT *)calloc(dim1,sizeof(INT));
  INT *m = (INT *)calloc(dim1,sizeof(INT));
  INT *mi = (INT *)calloc(2*dim1,sizeof(INT));
  INT *mip = (INT *)calloc(2*dim1,sizeof(INT));
  INT *vertk = (INT *)calloc(dim1,sizeof(INT));
  INT *vertj = (INT *)calloc(dim1,sizeof(INT));
  INT *facesk = (INT *)calloc(2*dim*nvface,sizeof(INT));
  INT *facesj = (INT *)calloc(2*dim1*nvface,sizeof(INT));
  INT *ti=vertj,*tip=vertk;
  //scalars
  INT nv=-10,kel, jel,nk,nj,flag;
  INT numv,kf,kfp,kz,kdim,kj;
  INT i,j,k,iaa, iab;
  INT *nodesj=NULL,*nodesk=NULL;
  // place holders (short hand)
  iCSRmat *fel2el=mc->fullel2el;
  INT **nd=mc->nd;
  //  INT **elneib=mc->elneib;
  //  INT **el2fnum=mc->el2fnum;
  INT *iindex,*iindexp;
  //  for(kz=(dim-1);kz>=0;kz--){
  //serching for intersection of kz-dimensional faces.
  //    kdim=(1<<kz);
  //  for(knnz=0;knnz<mc->bfs->nnz;knnz++){    
  //    kel=mc->bfs->JA[knnz];
  INT neg,nsall,nvall,nvold;

  nvall=0;nsall=0;

  print_full_mat_int(g0->nel,c2s->nvcube+1,g0->mnodes,"mel0");
  
  for(kel=0;kel<mc->nel;kel++){
    // we have not been here, so let us set the initial indexing to be the original indexing; 
    nvold=nvall;
    for(i=0;i<scin[kel]->nv;i++){
      mc->iindex[kel][i]=i+nvall;// this is the global number if there are no removals of vertices. 
    }
    neg=0;
    iaa=fel2el->IA[kel];
    iab=fel2el->IA[kel+1];
    if((iab-iaa)<=0){
      if(g0->print_level>4)
	fprintf(stdout,"\nmacroelement=%d; vertices=%d; overlaps=%d;",kel,scin[kel]->nv,neg);
      nvall+=scin[kel]->nv;      
      nsall+=scin[kel]->ns;
      continue;
    } 
    scp=scin[kel];// sc complex visited earlier.
    iindexp=mc->iindex[kel];
    nodesk=(g0->mnodes+kel*(nvcube+1));
    for(i=0;i<c2s->nf;i++){
      for(j=0;j<nvface;j++){
	k=c2s->faces[i*nvface+j];
	facesk[i*nvface+j]=nodesk[k];
      }
    }
    for (kj=iaa;kj<iab;kj++){
      jel=fel2el->JA[kj];
      nodesj=(g0->mnodes+jel*(nvcube+1));
      for(i=0;i<c2s->nf;i++){
	for(j=0;j<nvface;j++){
	  k=c2s->faces[i*nvface+j];
	  facesj[i*nvface+j]=nodesj[k];
	}
      }
      numv=0;
      //	if(fel2el->val[kj]!=kdim) continue;
      kdim=fel2el->val[kj];
      kz=((INT )floor(log2(((REAL) kdim)-1e-3)))+1;
      //      scomplex *sc=scin[jel];
      iindex=mc->iindex[jel];
      //      fprintf(stdout,"\nkel=%d,jel=%d,kz=%d,kdim=%d\n",kel,jel,kz,kdim);
      // find now how many times these two elements intersect:
      for(k=0;k<c2s->nvcube;k++){
	i=nodesk[k];
	i=locate0(i,nodesj,c2s->nvcube);
	if(i<0) continue;
	//	vertk[numv]=k;
	//	vertj[numv]=i; // this uses that the positions are the same, we need to change this. 
	vertk[numv]=nodesk[k];
	vertj[numv]=nodesj[i];
	numv++;
      }
      //      nk=locate1(mip,vertk,numv,c2s->faces,c2s->nf,c2s->nvface);
      //      nj=locate1(mi,vertj,numv,c2s->faces,c2s->nf,c2s->nvface);
      nk=locate1(mip,vertk,numv,facesk,c2s->nf,c2s->nvface);
      nj=locate1(mi,vertj,numv,facesj,c2s->nf,c2s->nvface);
      if((nk == nj) && nk==(dim-kz) && nj==(dim-kz)) {
	if(nk==1){
	  fprintf(stdout,"\nkel=%d,jel=%d,mik=%d,mij=%d;",kel,jel,mip[0],mi[0]);
	  print_full_mat_int(1,nvcube,nodesk,"nodesk");
	  print_full_mat_int(1,numv,vertk,"vertk");
	  print_full_mat_int(c2s->nf,nvface,facesk,"facesk");
	  print_full_mat_int(1,dim,mip,"mik");
	  print_full_mat_int(1,nvcube,nodesj,"nodesj");
	  print_full_mat_int(1,numv,vertj,"vertj");
	  print_full_mat_int(c2s->nf,nvface,facesj,"facesj");
	  print_full_mat_int(1,dim,mi,"mij");
 	}
	for(j=0;j<nj;j++){
	  /* 
	     for two macroelements kel and jel, from faces mi[] and
	     mip[] whose intersection forms the intersection of kel and jel
	     macroelements get the ti[] and tip[] arrays which
	     describe on which boundary using the number of
	     divisions we have the corresponding vertices lying
	  */
	  if(mi[j]<dim){
	    mi[j]=dim-(mi[j]+1);// which place in the array m[] this face corresponds to.
	    ti[j]=0;
	  } else {
	    mi[j]=dim-((mi[j]%dim)+1);
	    ti[j]=nd[jel][mi[j]];
	  }
	  if(mip[j]<dim){
	    mip[j]=dim-(mip[j]+1);
	    tip[j]=0;
	  } else {
	    mip[j]=dim-((mip[j]%dim)+1);
	    tip[j]=nd[kel][mip[j]];
	  }	    
	}
	nv=0;
	for(kfp=0;kfp<scp->nv;kfp++){
	  coord_lattice(mp,dim,kfp,scp->nv,nd[kel]);
	  flag=1;
	  for(j=0;j<nj;j++) if(mp[mip[j]]!=tip[j]){flag=0;break;}
	  if(flag){
	    memcpy(m,mp,dim*sizeof(INT));
	    for(k=0;k<nk;k++)m[mi[k]]=ti[k];
	    kf=num_lattice(m,dim,nd[jel]);
	    if(kz==2){
	      fprintf(stdout,"\nkel=%d; jel=%d,kfp=%d,iindexp=%d,iindex=%d",kel,jel,kfp,iindexp[kfp],iindex[kf]);fflush(stdout);
	    }
	    if(iindexp[kfp]>=nvold){
	      iindexp[kfp]=iindex[kf];
	      neg++;
	    }
	  } else {
	    if(kz==2){
	      fprintf(stdout,"\nkel=%d; jel=%d,kfp=%d",kel,jel,kfp);fflush(stdout);
	    }
	    if(iindexp[kfp]>=nvold) {
	      iindexp[kfp]=nv+nvold;
	      nv++;
	    }
	  }
	}
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
      } else {
	fprintf(stderr,"\n*** An error in counting overlaps in %s ***\n",__FUNCTION__);
      }
    }
    if(g0->print_level>0)
      fprintf(stdout,"\nmacroelement=%d; total=%d; overlaps:%d",kel,scp->nv,neg);
    nvall+=nv;
    nsall+=scin[kel]->ns;
    //fprintf(stdout,"\nGLOBALLY:v_total=%d; s_total=%d",nvall,nsall);fflush(stdout);
  }
  if(g0->print_level>5){
    for(kel=0;kel<mc->nel;kel++){
      fprintf(stdout,"\nelement:%d;",kel);
      for(i=0;i<scin[kel]->nv;i++){
	coord_lattice(m,dim,i,scin[kel]->nv,nd[kel]);
	fprintf(stdout,"\n\t[%d=( ",i);
	for(j=0;j<dim;j++){
	  fprintf(stdout,"%d ",m[j]);
	}
	fprintf(stdout,")]<-->%d",mc->iindex[kel][i]);
	fprintf(stdout,"; x=(");
	for(j=0;j<dim;j++)
	  fprintf(stdout,"%.3f ", scin[kel]->x[i*dim+j]);
	fprintf(stdout,")");
      }
    }
    fprintf(stdout,"\n");
  }
  free(m);
  free(mp);
  free(mi);
  free(mip);
  free(vertk);
  free(vertj);
  //  free(ti); //same as vertk
  //  free(tip); //same as vertj
  /* nvall=0;nsall=0; */
  /* INT nvloc,nsloc; */
  /* for(kel=0;kel<g0->nel;kel++){ */
  /*   nvloc=1;nsloc=1; */
  /*   for(i=0;i<g0->dim;i++){ */
  /*     nvloc*=(nd[kel][i]+1); */
  /*     nsloc*=nd[kel][i]; */
  /*   } */
  /*   nvall+=nvloc; */
  /*   nsall+=nsloc*(scin[0]->factorial); */
  /* } */
  /* scomplex_merge(scin,			\ */
  /* 		 nsall, nvall,			\ */
  /* 		 scin[0]->cc, scin[0]->bndry_cc,	\ */
  /* 		 g0,c2s); */
  scomplex_merge1(nvall,nsall,mc,scin,c2s);
  return;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ*/
scomplex *generate_grid(input_grid *g0)
{
  /* 
     From an input grid loops over the macroelements and generates
     meshes depending on the division given in every dimension. 
  */
  input_grid *g;
  INT i,j,kel,pmem;
  /***************************************************************************/
  cube2simp *c2s=cube2simplex(g0->dim);
  /***************************************************************************/
  INT nvcube=c2s->nvcube;
  scomplex **sc=malloc(g0->nel*sizeof(scomplex *));
  pmem=2*g0->nel*c2s->nf*c2s->nvface;
  if(pmem<2*g0->ne) pmem=2*g0->ne;
  if(pmem<2*g0->nel) pmem=2*g0->nel;
  if(pmem<2*g0->nf) pmem=2*g0->nf;
  INT *p=calloc(pmem,sizeof(INT));
  //  print_full_mat_int(g0->nel,(c2s->nvcube+1),g0->mnodes,"ex");
  ilexsort(g0->nf, (c2s->nvface+1),g0->mfaces,p);
  ilexsort(g0->nel,(c2s->nvcube+1),g0->mnodes,p);
  /* for(kel=0;kel<g0->nel;kel++){ */
  /*   print_full_mat_int(1,c2s->nvcube+1,g0->mnodes+kel*(c2s->nvcube+1),"mnodes33"); */
  /* } */
  /*-------------------------------------------------------------------*/
  //  INT *efound=calloc(c2s->ne*(g0->nel),sizeof(INT));  
  INT **nd=calloc(g0->nel,sizeof(INT*));
  for(i=0;i<g0->nel;i++){
    nd[i]=calloc(c2s->nf,sizeof(INT)); /* to hold the number of
 					 divisions in every coordinate
					 direction */
  }
  for(kel=0;kel<g0->nel;kel++)
    for(i=0;i<c2s->n;i++){
      nd[kel][i]=-1;
    }
  g=malloc(1*sizeof(input_grid));
  /**/
  g->title=g0->title;
  g->dgrid=g0->dgrid;
  g->fgrid=g0->fgrid;
  g->dvtu=g0->dvtu;
  g->fvtu=g0->fvtu;
  g->print_level=g0->print_level;
  g->ref_type=g0->ref_type;
  g->nref=g0->nref;
  g->err_stop=g0->err_stop;
  /**/
  g->dim=c2s->n;
  g->ncsys=g0->ncsys;
  g->nv=c2s->nvcube;
  g->nf=c2s->nf;
  g->ne=c2s->ne;
  g->nel=1;
  input_grid_arrays(g);
  /* reassign this as these are the same as g0 */
  free(g->systypes);   g->systypes=g0->systypes;
  free(g->syslabels);   g->syslabels=g0->syslabels;
  free(g->ox); g->ox=g0->ox;
  INT chng=1,iter=0,maxiter=1024;  
  //  INT je,kj,k2,iel2v,jel2v,k1,kface,kbnd,found;
  set_edges(g0,c2s);  
  while(chng && (iter<maxiter)){
    iter++;
    // make the divisions in g0->seg consistent;
    chng=set_ndiv_edges(g,g0,c2s,nd,iter);
  }
  // get the macroelement mesh in a structure
  macrocomplex *mc=set_mmesh(g0,c2s,p);
  /*PLACE HOLDERS*/
  /* iCSRmat *fullel2el=mc->fullel2el; */
  /* INT **elneib=mc->elneib; */
  //  INT *etree=mc->etree;
  INT **el2fnum=mc->el2fnum;
  INT *isbface=mc->isbface;
  INT *bcodesf=mc->bcodesf;
  iCSRmat *bfs0=mc->bfs;
  fprintf(stdout,"\n%%DFS(domains): %d connected components",mc->cc);
  fprintf(stdout,"\n%%DFS(boundaries): %d connected components",mc->bndry_cc);
  if(g0->print_level>3){
    fprintf(stdout,"\n%%Input faces with codes: %d; Total number of faces:%d\n",g0->nf,mc->nf);
  }else if(g0->print_level>10){
    fprintf(stdout,"\nbfs0=[");
    icsr_print_matlab_val(stdout,bfs0);
    fprintf(stdout,"];");
    fprintf(stdout,"\nbfs=sparse(bfs0(:,1),bfs0(:,2),bfs0(:,3));\n");
    print_full_mat_int(1,mc->nf,mc->isbface,"isbface");
    print_full_mat_int(1,mc->nf,mc->bcodesf,"bcodesf");
    fprintf(stdout,"\n*****   nf=%d (%d)******* \n",g0->nf,mc->nf);
    fprintf(stdout,"\nfel2el0=[");
    icsr_print_matlab_val(stdout,mc->fullel2el);
    fprintf(stdout,"];");
    fprintf(stdout,"\nfel2el=sparse(fel2el0(:,1),fel2el0(:,2),fel2el0(:,3));\n");
    print_full_mat_int(1,mc->nf,mc->isbface,"isbface");
    print_full_mat_int(1,mc->nf,mc->bcodesf,"bcodesf");
  }
  /***********************************************************************/    
  /* set the divisions on every edge now; since they are consistent we
   have: */
  if(set_ndiv_edges(g,g0,c2s,nd,0)) {
    fprintf(stderr,"\n\n***ERR in %s: the divisions of the edges cannod be inconsistent during second call of set_ndiv_edges()\n\n",__FUNCTION__);
    exit(4);
  }
  mc->nd=nd;
  for(kel=0;kel<g0->nel;kel++){
    for(i=0;i<c2s->n;i++){
      if(nd[kel][i]<=0) nd[kel][i]=1;
    }
  }
  // use bfs to split elements:
  if(g0->print_level>5) input_grid_print(g0);
  INT nsall,nvall;
  nsall=0;nvall=0;
  /* now nd is known, let us allocate iindex */
  for(kel=0;kel<g0->nel;kel++){
    nvall=1;
    for(i=0;i<g0->dim;i++)
      nvall*=(nd[kel][i]+1);
    mc->iindex[kel]=calloc(nvall,sizeof(INT));
    // TRUE    mc->flags[kel]=g0->mnodes[nvcube+kel*(nvcube+1)];
    mc->flags[kel]=2*kel+1;
  }
  INT intype=g0->ref_type-1,kj,jel;  
  INT *codef=calloc(c2s->nf,sizeof(INT));
  INT *isbndf=calloc(c2s->nf,sizeof(INT));
  INT kjj;
  for(kj=0;kj<bfs0->row;kj++){
    if(g0->ref_type>=0) intype++;
    else intype=-1;
    for(kjj=bfs0->IA[kj];kjj<bfs0->IA[kj+1];kjj++){
      jel=bfs0->JA[kjj];    
      if((intype>=0) && (mc->etree[jel]<0)) intype=g0->ref_type; //reset reftype;
      //    print_full_mat_int(1,c2s->nf,elneib[kel],"neib");
      memcpy(g->mnodes,(g0->mnodes+jel*(nvcube+1)),(nvcube+1)*sizeof(INT));
      for(i=0;i<c2s->nvcube;i++){
	j=g->mnodes[i];// vertex number (global)
	g->csysv[i]=g0->csysv[j]; 
	g->labelsv[i]=g0->csysv[j];
	memcpy((g->xv+i*g->dim),(g0->xv+j*g0->dim),g->dim*sizeof(REAL));
      }
      /***************************************************/
      /*element code is in mc->flags[]; we now do the face codes:*/      
      for(i=0;i<c2s->nf;i++){
	codef[i]=bcodesf[el2fnum[jel][i]];
	isbndf[i]=isbface[el2fnum[jel][i]];       
      }      
      sc[jel]=umesh(g->dim,nd[jel],c2s,					\
		    isbndf,codef,mc->flags[jel],			\
		    intype);
      nsall+=sc[jel]->ns;
      nvall+=sc[jel]->nv;
      //    if((mc->fullel2el->IA[jel+1]-mc->fullel2el->IA[jel])<=0){
      for(i=0;i<sc[jel]->nv;i++){
	mc->iindex[jel][i]=i;
      }
      //}
      //      fprintf(stdout,"\n%%Mapping back to the macroelement...\n");
      map2mac(sc[jel],c2s,g);
      fprintf(stdout,"\n%%mesh(macroelement=%d): nv=%d; nsimp=%d",jel,sc[jel]->nv,sc[jel]->ns);      
    }
  }
  //  fprintf(stdout,"\n%%Removing overlaps...");
  fix_grid(mc,sc,c2s,g0);
  fprintf(stdout,"..\tdone");
  if(g0->print_level>0){
    fprintf(stdout,"\n%%merged(macroelements=[%d..%d]): vertices=%d; simplices=%d", \
	    0,mc->nel-1,sc[0]->nv,sc[0]->ns);  
    fprintf(stdout," ..done.\n\n");
  }
  free(p);
  free(isbndf);
  free(codef);
  macrocomplex_free(mc);
  cube2simp_free(c2s);
  /* prepare for adaptive refinement */

  //  haz_scomplex_print(sc[0],0,"TTT");
  return sc[0];
  //    find_nbr(sc[0]->ns,sc[0]->nv,sc[0]->n,sc[0]->nodes,sc[0]->nbr);
  //  INT *wrk1=calloc(5*(sc[0]->n+2),sizeof(INT));
  /* construct bfs tree for the dual graph */
  //    abfstree(0,sc[0],wrk1,g0->print_level);
  //  free(wrk1);
  //  return sc[0];  
}
/****************************************************************************************/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
void scomplex_merge(scomplex **sc0,			\
		    const INT nsall, const INT nvall,	\
		    const INT cc, const INT bndry_cc,	\
		    input_grid *g0,cube2simp *c2s)
{
  /* combains an array of simplicial complexes together */
  if(g0->nel==1) return;
  scomplex *sc=sc0[0];
  INT n1=(sc->n+1),nv0,ns0,nv=nvall,ns=nsall;
  INT kel,i,ii,j,in1,iin1;
  sc->marked=realloc(sc->marked,ns*sizeof(INT));
  sc->gen=realloc(sc->gen,ns*sizeof(INT));
  sc->nbr=realloc(sc->nbr,ns*n1*sizeof(INT));
  sc->parent=realloc(sc->parent,ns*sizeof(INT));
  sc->child0=realloc(sc->child0,ns*sizeof(INT));
  sc->childn=realloc(sc->childn,ns*sizeof(INT));
  sc->nodes=realloc(sc->nodes,ns*n1*sizeof(INT));
  sc->bndry=realloc(sc->bndry,nv*sizeof(INT));
  sc->csys=realloc(sc->csys,nv*sizeof(INT));/* coord sys: 1 is polar, 2
					    is cyl and so on */
  /*connected components*/
  sc->cc=cc;sc->bndry_cc=bndry_cc;
  sc->flags=(INT *)realloc(sc->flags,ns*sizeof(INT));
  sc->x=(REAL *)realloc(sc->x,nv*(sc->n)*sizeof(REAL));
  sc->vols=(REAL *)realloc(sc->vols,ns*sizeof(REAL));
  sc->fval=(REAL *)realloc(sc->fval,nv*sizeof(REAL)); // function values at every vertex; not used in general;
  fprintf(stdout,"\nnsall=%d,nvall=%d",nsall,nvall);fflush(stdout);
  for(kel=1;kel<g0->nel;kel++){
    ns0=sc->ns;nv0=sc->nv;
    for (ii = 0;ii<sc0[kel]->ns;ii++) {
      i=ii+ns0;
      sc->marked[i] = sc0[kel]->marked[ii];
      sc->gen[i] = sc0[kel]->gen[ii];
      sc->parent[i]=sc0[kel]->parent[ii];
      sc->child0[i]=sc0[kel]->child0[ii];
      sc->childn[i]=sc0[kel]->childn[ii];
      sc->flags[i]=sc0[kel]->flags[ii];
      sc->vols[i]=sc0[kel]->vols[ii];
      in1=i*n1;
      iin1=ii*n1;
      for(j=0;j<n1;j++){
	sc->nodes[in1+j]=sc0[kel]->nodes[iin1+j]+nv0;
	sc->nbr[in1+j]=sc0[kel]->nbr[iin1+j]+ns0;
      }
    }
    for (ii = 0;ii<sc0[kel]->nv;ii++) {
      i=ii+nv0;
      sc->bndry[i]=sc0[kel]->bndry[ii];
      sc->csys[i]=sc0[kel]->csys[ii];
      sc->fval[i]=sc0[kel]->fval[ii];
      in1=i*sc->n;
      iin1=ii*sc->n;
      for(j=0;j<sc->n;j++)
	sc->x[in1+j]=sc0[kel]->x[iin1+j];
    }
    sc->nv+=sc0[kel]->nv;
    sc->ns+=sc0[kel]->ns;
    haz_scomplex_free(sc0[kel]);
  }
   fprintf(stdout,"\nsc->nv=%d,sc->ns=%d; nvall=%d,nsall=%d\n",sc->nv,sc->ns,nvall,nsall);fflush(stdout);
  return;
}
/*ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ*/
/****************************************************************************************/
