/*! \file src/amr/macroelements.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note all files related to constructing mesh of macroelements and
 *  splitting it into simplices
 *  \date 20170715 
 *  \modified 20190715 (ltz);
 *  \note: modified by ltz on 20210813
 */
#include "hazmath.h"
/**************************************************************************/
/*!
   * \fn void icsr_add (iCSRmat *A, iCSRmat *B, iCSRmat *C)
   *
   * \brief symbolic sparse matrix addition C = alpha*A + beta*B
   *
   * \param A         Pointer to the iCSRmat matrix A
   * \param B         Pointer to the iCSRmat matrix B
   * \param C         Pointer to iCSRmat matrix with structure A + B; 
   */
/**********************************************************************/
static INT icsr_add(iCSRmat *A,iCSRmat *B,iCSRmat *C)
{
  INT i,j,k,l;
  INT count=0, added, countrow;
  INT status = SUCCESS;
  // both matrices A and B are NULL
  if (A == NULL && B == NULL) {
    printf("%%%% ERROR: both matrices are null in %s\n", __FUNCTION__);
    status = ERROR_MAT_SIZE;
    goto FINISHED;
  }
  // dimensions  mismatch!
  if ((A->row != B->row) || (A->col != B->col)) {
    printf("### ERROR HAZMATH DANGER: Dimensions of matrices do not match!!! %s\n", __FUNCTION__);
    status = ERROR_MAT_SIZE;
    goto FINISHED;
  }
  // Both matrices A and B are neither NULL or empty
  C->row=A->row; C->col=A->col;
  C->IA=(INT*)calloc((C->row+1),sizeof(INT));
  // allocate work space for C->JA and C->val
  C->JA=(INT *)calloc((A->nnz+B->nnz),sizeof(INT));
  C->val=(INT *)calloc((A->nnz+B->nnz),sizeof(INT));
  // initialize C->IA
  memset(C->IA, 0, sizeof(INT)*(C->row+1));
  for (i=0; i<(A->nnz+B->nnz); ++i) {
    C->JA[i]=-1;
  }
  for (i=0; i<A->row; ++i) {
    countrow = 0;
    for (j=A->IA[i]; j<A->IA[i+1]; ++j) {
      //      C->val[count] = alpha * A->val[j];
      C->val[count] = A->val[j];// only in A. 
      C->JA[count] = A->JA[j];
      C->IA[i+1]++;
      count++;
      countrow++;
    } // end for js

    for (k=B->IA[i]; k<B->IA[i+1]; ++k) {
      added = 0;
      for (l=C->IA[i]; l<C->IA[i]+countrow+1; l++) {
        if (B->JA[k] == C->JA[l]) {
	  // getting to here means this is in the intersection of the patterns. We just skip this...
          added = 1;
	  //	  C->val[l] = C->val[l] + B->val[k];
	  // here, in this particular application A->val must be equal to b->val
	  //	  fprintf(stdout,"\nA(%d,%d)=%d .eq. %d=B(%d,%d)",i,C->JA[l],C->val[l],B->val[k],i,B->JA[k]);
	  break;
        }
      } // end for l

      if (added == 0) {
	//        C->val[count] = beta * B->val[k];
        C->val[count] = B->val[k]; // only in B, just add it to the pattern.
        C->JA[count] = B->JA[k];
        C->IA[i+1]++;
        count++;
      }
    } // end for k
    C->IA[i+1] += C->IA[i];
  }

  C->nnz = count;
  C->JA  = (INT *)realloc(C->JA, (count)*sizeof(INT));
  C->val = (INT *)realloc(C->val, (count)*sizeof(INT));

FINISHED:
  return status;
}
/**************************************************************************/
/*!
 * \fn static INT ilog2(const INT k)
 *
 * \brief integer log for an integer. For example 1=ilog2(2);
 *
 * \param 
 *
 * \return
 *
 * \note
 *
 */
static INT ilog2(const INT k)
{
  if(k>0)
    return (((INT )floor(log2(((REAL) k)-1e-3)))+1);
  else
    return -(1<<30);
}
/**********************************************************************/
/*!
 * \fn void align_lattice(INT *pd,const INT nkj, INT *nodes0,
 *  		   INT *nodes1, cube2simp *c2s)
 *
 * \brief constructs permutation of the vertices of an element based
 *        on the numbering in a neighboring element. 
 *
 * \param 
 *
 * \return
 *
 * \note 
 *
 */
void align_lattice(INT *pd, const INT nkj,	\
		   INT *nodes0,			\
		   INT *nodes1,	\
		   cube2simp *c2s)
{
  INT dim=c2s->n,i;
  //fprintf(stdout,"\n****%d ^^^\n",nkj);
  for(i=0;i<dim;i++)pd[i]=i;  
  if(nkj<2){
    return; // meaning only one common vertex, so we should do nothing;
  }
  INT j,d0,d1,k0[2],k1[2];  
  for(i=0;i<dim;i++) pd[i]=i+1;
  for(i=0;i<c2s->ne;i++){
    k0[0]=nodes0[c2s->edges[2*i]];
    k0[1]=nodes0[c2s->edges[2*i+1]];
    d0=ilog2(abs(c2s->edges[2*i+1]-c2s->edges[2*i]));
    for(j=0;j<c2s->ne;j++){
      k1[0]=nodes1[c2s->edges[2*j]];
      k1[1]=nodes1[c2s->edges[2*j+1]];
      d1=ilog2(abs(c2s->edges[2*j+1]-c2s->edges[2*j]));      
      if(k0[0]==k1[0]&&k0[1]==k1[1]){
	pd[d0]=(d1+1);
	break;
      }else if(k0[0]==k1[1]&&k0[1]==k1[0]){
	// reversed order of the common edge ends
	pd[d0]=-(d1+1);
	break;
      }
    }
  }
  return;
}
/**********************************************************************/
/*!
 * \fn macrocomplex set_mmesh(input_grid *g0,cube2simp *c2s,INT *wrk) 
 *
 * \brief prepare macro element mesh (mmesh) for passing to the mesh
 *        generator. wrk is working integer array should have at least
 *        size 2*nvcube+2
 *
 * \param wrk (working array of size 2*nvcube +2)
 *
 * \return
 *
 * \note
 *
 */
macrocomplex *set_mmesh(input_grid *g0,cube2simp *c2s,INT *wrk)
{
  INT i,j,j0,j1,kel,jel,ke;
  INT nvcube=c2s->nvcube,nvface=c2s->nvface;
  INT nel0,je,kj,k2,iel2v,jel2v,k1,kface,kbnd,found;
  INT *p=wrk; 
  INT *mnodes=p+nvcube+1;
  /*macro complex creation:)*/
  macrocomplex *mc=malloc(sizeof(macrocomplex));
  mc->nel=nel0=g0->nel; //important to set;
  mc->cc=mc->bndry_cc=1;
  mc->nf=mc->nfi=mc->nfb=-1;
  /*macro complex allocation*/
  /********************************************************************/
  mc->flags=calloc(mc->nel,sizeof(INT));
  /**/
  mc->nd=calloc(mc->nel,sizeof(INT*));
  for(i=0;i<mc->nel;i++){
    mc->nd[i]=calloc(c2s->nf,sizeof(INT)); /* to hold the number of
					      divisions in every coordinate
					      direction.*/
  }
  for(kel=0;kel<mc->nel;kel++)
    for(i=0;i<c2s->n;i++){
      mc->nd[kel][i]=-1;
    }
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
  mc->etree=NULL;// array
  mc->bcodesf=NULL; /* array: later when num faces is known*/
  mc->isbface=NULL; /* array: later when num faces is known*/
  mc->bfs=NULL; /* icsr_mat later when bfs is called*/
  mc->fullel2el=NULL; /* later when el2el is formed*/
  /********************************************************************/
  /* let us also allocate all the rest here */
  /*  print_full_mat_int(g0->nel,(c2s->nvcube+1),g0->mnodes,"ex");*/
  ilexsort(g0->nf, (c2s->nvface+1),g0->mfaces,p);
  //  ilexsort(g0->nel,(c2s->nvcube+1),g0->mnodes,p);
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
  icsr_free(v2el);free(v2el);
  /* create mc->fullel2el for later and work here with the copy */
  mc->fullel2el=malloc(sizeof(iCSRmat));
  mc->fullel2el[0]=icsr_create(el2el->row,el2el->col,el2el->nnz);  
  memcpy(mc->fullel2el->IA,el2el->IA,(mc->fullel2el->row+1)*sizeof(INT));
  memcpy(mc->fullel2el->JA,el2el->JA,(mc->fullel2el->nnz)*sizeof(INT));
  memcpy(mc->fullel2el->val,el2el->val,(mc->fullel2el->nnz)*sizeof(INT));
  /***********************************************************/  
  icsr_tri(mc->fullel2el,'l');
  icsr_nodiag(mc->fullel2el);
  /***********************************************************/  
  /*   shrink the el2el matrix by removing any entry with value not
       equal to nvface=(number of vertices per face); */
  el2el->nnz=el2el->IA[0];
  for(kel=0;kel<nel0;kel++){    
    j0=el2el->IA[kel];
    j1=el2el->IA[kel+1];
    el2el->IA[kel]=el2el->nnz;
    for(ke=j0;ke<j1;ke++){
      jel=el2el->JA[ke];
      //      if(jel==kel) {
      // do not need the diagonal.
      //	el2el->JA[el2el->nnz]=kel;
      //	el2el->val[el2el->nnz]=1;
      //	el2el->nnz++;       
      //      continue;
      //      }
      if(jel==kel) continue;
      if(el2el->val[ke]!=nvface) continue;
      el2el->JA[el2el->nnz]=jel;
      el2el->val[el2el->nnz]=el2el->val[ke];
      el2el->nnz++;
    }
  }  
  el2el->IA[nel0]=el2el->nnz;
  el2el->JA=realloc(el2el->JA,el2el->nnz*sizeof(INT));
  el2el->val=realloc(el2el->val,el2el->nnz*sizeof(INT));
  iCSRmat *blk_dfs=run_dfs(nel0,el2el->IA, el2el->JA);
  ////////////////////////////////////////////////////////////////////////////
  mc->cc=blk_dfs->row;
  ivector roots,anc;
  roots.row=mc->cc;
  roots.val=calloc(roots.row,sizeof(INT));
  for(kel=0;kel<mc->cc;++kel){
    roots.val[kel]=blk_dfs->JA[blk_dfs->IA[kel]];
  }
  //
  mc->bfs=run_bfs(mc->nel,el2el->IA,el2el->JA,&roots,&anc,(mc->nel+1));
  //
  ivec_free(&roots);
  mc->etree=anc.val;// put this pointer in place;
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
  //////////////////////////////
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
	      memcpy((f2v->JA+f2v->IA[kface]),facei,nvface*sizeof(INT));
	      elneib[kel][ke]=jel;
	      elneib[jel][je]=kel;
	      el2fnum[kel][ke]=kface;
	      el2fnum[jel][je]=kface;
	      //	      fprintf(stdout,"\nxkel=%d,xke=%d;xjel=%d,xje=%d:::kface=%d",kel,ke,jel,je,kface);
	      kface++;// pointer for the face2vertex matrix
	    }
	    found=1;
	    //	    print_full_mat_int(1,nvface,facei,"facei");
	    //	    print_full_mat_int(1,nvface,facej,"facej");
	  }
	}
      }
      if(!found){
	isbface[kface]=1;
	memcpy((f2v->JA+f2v->IA[kface]),facei,nvface*sizeof(INT));
	/* 
	   in this way bcodesf[1:elneib[kel][ke]] gives us the code of the corresponding face;; this can also be linked to f2v and so on. 
	*/
	//	fprintf(stdout,"\nykel=%d,yke=%d:::kface=%d",kel,ke,kface);
	el2fnum[kel][ke]=kface;
	kface++;
      }     
    }
  }
  /********************************************************************/
  for(i=0;i<f2v->IA[mc->nf];i++) f2v->val[i]=1;
  /*******************************************************************/   
  INT cfmax,cfmin;
  cfmax=g0->mfaces[nvface]; // first face: get the code;
  cfmin=cfmax; 
  for(je=1;je<g0->nf;je++){
    kbnd=g0->mfaces[je*(nvface+1) + nvface]; // get the code;
    //    fprintf(stdout,"g0_code(%d)=%%d;cfmax=%d\n",je,kbnd,cfmax);fflush(stdout);
    if(kbnd>cfmax) cfmax=kbnd;
    if(kbnd<cfmin) cfmin=kbnd;
  }
  cfmax++; cfmin--;
  if(cfmax<2) cfmax=2;
  if(cfmin>=0) cfmin=-1;
  //  fprintf(stdout,"CFMAX=%d;CFMIN=%d\n",cfmax,cfmin);fflush(stdout);
  for(i=0;i<mc->nf;i++){
    bcodesf[i]=cfmax;// something that is not a prescribed code:
    j1=f2v->IA[i+1]-f2v->IA[i];
    if(j1>0){
      memcpy(facei,(f2v->JA+f2v->IA[i]),j1*sizeof(INT));
      for(je=0;je<g0->nf;je++){
	memcpy(facej,(g0->mfaces+je*(nvface+1)),nvface*sizeof(INT));
	kbnd=g0->mfaces[je*(nvface+1)+nvface];
	if(aresame(facei,facej,nvface)){
	  /* if the face was in the list, take its code. In this way
	     we only take codes of faces that are actually
	     macroelement faces. */
	  bcodesf[i]=kbnd;
	  break;
	}
      }
    }
  }
  /* 
   * set all interior faces faces with no code to 0 code and all
   * boundary faces with no code to Dirichlet code 1. If a boundary
   * face has a code 0 it is set to cmin; cmin<0 always;
   */
  for(i=0;i<mc->nf;i++){
    if(bcodesf[i]!=cfmax) continue;
    if(isbface[i])
      bcodesf[i]=1;
    else
      bcodesf[i]=0;
    //      fprintf(stdout,"\n[%d]=%d",i,bcodesf[i]);
  }
  //  fprintf(stdout,"\n");fflush(stdout);
  for(i=0;i<mc->nf;i++){
    if(isbface[i] && (bcodesf[i]==0)){
      bcodesf[i]=cfmin; // faces on the boundary with code 0. 
    }
  }
  /*******************************************************************/   
  /* for(i=0;i<mc->nf;i++){ */
  /*   fprintf(stdout,"\nFace(%d; code=%d)=[",i,bcodesf[i]); */
  /*   for(je=f2v->IA[i];je<f2v->IA[i+1];++je){ */
  /*     fprintf(stdout,"%d ",f2v->JA[je]); */
  /*   } */
  /*   fprintf(stdout,"]"); */
  /* } */
  /* fprintf(stdout,"\n");fflush(stdout); */
  /***********************************************************************/
  /*Connected components on the boundaries. First find face to face map*/
  iCSRmat *v2f=malloc(sizeof(iCSRmat));
  icsr_trans(f2v,v2f);
  iCSRmat *f2f=malloc(sizeof(iCSRmat));
  /*******************************************************************/    
  icsr_mxm(f2v,v2f,f2f);
  /*******************************************************************/    
  icsr_free(v2f);free(v2f);
  /*******************************************************************/    
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
    if(isbface[i]==0) continue;
    for(ke=j0;ke<j1;ke++){
      je=f2f->JA[ke];
      if(je==i){
	//	f2f->JA[f2f->nnz]=i;
	//	f2f->val[f2f->nnz]=0;
	//	f2f->nnz++;
      	continue;
      }
      if((isbface[je]==0)||(f2f->val[ke]!=(nvface/2))) continue;
      //      if((je==i)||(isbface[je]==0)||(f2f->val[ke]!=(nvface/2))) continue;
      f2f->JA[f2f->nnz]=je;
      f2f->val[f2f->nnz]=f2f->val[ke];
      f2f->nnz++;
    }
  }
  /*******************************************************************/    
  /* connected comps on the boundary */
  icsr_free(blk_dfs);free(blk_dfs);
  blk_dfs=run_dfs(mc->nf,f2f->IA, f2f->JA);
  mc->bndry_cc=0;
  for(i=0;i<blk_dfs->row;++i){
    found=blk_dfs->IA[i+1]-blk_dfs->IA[i];
    if(found>1) mc->bndry_cc++;
  }
  icsr_free(f2f);free(f2f);
  icsr_free(blk_dfs);free(blk_dfs);
  /*now use the bfs*/
  INT lvl,keok,swp,keswp;
  /* for(kj=0;kj<mc->nel;++kj){ */
  /*   print_full_mat_int(1,c2s->nf,elneib[kj],"elneib"); */
  /* } */
  /* for(kj=0;kj<mc->nel;++kj){ */
  /*    print_full_mat_int(1,c2s->nf,el2fnum[kj],"el2fnum"); */
  /* } */
  for(lvl=0;lvl<mc->bfs->row;lvl++){
    j0=mc->bfs->IA[lvl];
    j1=mc->bfs->IA[lvl+1];
    for(kj=j0;kj<j1;kj++){
      jel=mc->bfs->JA[kj];
      kel=mc->etree[jel];// ancestor, this stays unchanged
      if(kel>=0){
  	je=locate0(jel,elneib[kel], c2s->nf);
  	ke=locate0(kel,elneib[jel], c2s->nf);
  	keok=(je+c2s->n)%c2s->nf;
	//////////////////////////////////////////////
	/* fprintf(stdout,"\nKEL=%d,JEL=%d;  KE=%d;KEOK=%d;JE=%d",kel,jel,ke,keok,je);fflush(stdout); */
	/* print_full_mat_int(1,(nvface),(c2s->faces + je*nvface),"FACES00(KEL)"); */
	/* print_full_mat_int(1,(nvcube+1),(g0->mnodes+kel*(nvcube+1)),"NODES00(KEL)");	 */
	/* print_full_mat_int(1,(nvface),(c2s->faces+ke*nvface),"FACES00(JEL)"); */
	/* print_full_mat_int(1,(nvcube+1),(g0->mnodes+jel*(nvcube+1)),"NODES00(JEL)"); */
	///////////////////////////////////////////////////////////////
  	if(keok!=ke){
	  /*************************************************/	  
	  /*************************************************/
	  /* print_full_mat_int(1,c2s->nf,elneib[kel],"elneib1(kel)"); */
	  /* print_full_mat_int(1,c2s->nf,el2fnum[kel],"el2fnum1(kel)"); */
	  /* print_full_mat_int(1,c2s->nf,elneib[jel],"elneib1(jel)"); */
	  /* print_full_mat_int(1,c2s->nf,el2fnum[jel],"el2fnum1(jel)"); */
	  /*************************************************/	  
	  /*************************************************/
  	  //	  swap in jel:
  	  swp=elneib[jel][ke];
  	  elneib[jel][ke]=elneib[jel][keok];
  	  elneib[jel][keok]=swp;
  	  swp=el2fnum[jel][ke];
  	  el2fnum[jel][ke]=el2fnum[jel][keok];
  	  el2fnum[jel][keok]=swp;
	  /*************************************************/
	  /* print_full_mat_int(1,c2s->nf,elneib[kel],"elneib2(kel)"); */
	  /* print_full_mat_int(1,c2s->nf,el2fnum[kel],"el2fnum2(kel)"); */
	  /* print_full_mat_int(1,c2s->nf,elneib[jel],"elneib2(jel)"); */
	  /* print_full_mat_int(1,c2s->nf,el2fnum[jel],"el2fnum2(jel)"); */
	  /*************************************************/	  
  	  // we now need to swap vertices in g0->mnodes
  	  for(i=0;i<nvface;i++){
  	    facei[i]=c2s->faces[ke*nvface+i];
  	    keswp=(ke+c2s->n)%c2s->nf;
  	    facei[i+nvface]=c2s->faces[keswp*nvface+i];
  	    // use mnodes as work space here;
  	    mnodes[i]=c2s->faces[keok*nvface+i];
  	    keswp=(keok+c2s->n)%c2s->nf;
  	    mnodes[i+nvface]=c2s->faces[keswp*nvface+i];
  	  }
  	  for(i=0;i<nvcube;i++)
  	    el2v->JA[i]=g0->mnodes[jel*(nvcube+1)+facei[i]];
  	  for(i=0;i<nvcube;i++)
  	    g0->mnodes[jel*(nvcube+1)+mnodes[i]]=el2v->JA[i];
  	}       
	//////////////////////////////////////////////
	/* print_full_mat_int(1,(nvcube+1),(g0->mnodes+kel*(nvcube+1)),"NODES11(KEL)");	 */
	/* print_full_mat_int(1,(nvcube+1),(g0->mnodes+jel*(nvcube+1)),"NODES11(JEL)"); */
	///////////////////////////////////////////////////////////////
      }
    }
  }
  //  input_grid_print(g0);
  /* for(kj=0;kj<mc->nel;++kj){ */
  /*   print_full_mat_int(1,c2s->nf,elneib[kj],"elneib1"); */
  /* } */
  /* for(kj=0;kj<mc->nel;++kj){ */
  /*    print_full_mat_int(1,c2s->nf,el2fnum[kj],"el2fnum1"); */
  /* } */
  ////////////////////////// 
  //  input_grid_print(g0);
  //////////////////////////
  if(g0->print_level>10){
    fprintf(stdout,"\nbfs tree(in %s):\n",__FUNCTION__);
    for(lvl=0;lvl<mc->bfs->row;lvl++){
      j0=mc->bfs->IA[lvl];
      j1=mc->bfs->IA[lvl+1];
      for(kj=j0;kj<j1;kj++){
  	jel=mc->bfs->JA[kj];
  	kel=mc->etree[jel];// ancestor
  	fprintf(stdout,"** (%d--%d)",kel,jel);
      }
    }
    fprintf(stdout,"\n");
  }
  // BEGIN: REBUILD THE ELNEIB LIST BECAUSE SOME permutation of vertices may have made the elneib out of date
  /*-------------------------------------------------------------------*/
  /*initialize again*/
  for(kel=0;kel<mc->nel;kel++){
    for(i=0;i<c2s->nf;i++){
      elneib[kel][i]=-1;
    }
  }
  INT *f_save=calloc(c2s->nf,sizeof(INT));
  for(lvl=0;lvl<mc->bfs->row;lvl++){
    j0=mc->bfs->IA[lvl];
    j1=mc->bfs->IA[lvl+1];
    for(kj=j0;kj<j1;kj++){
      jel=mc->bfs->JA[kj];
      kel=mc->etree[jel];// ancestor and already established neighbor
      //      fprintf(stdout,"** (%d-->%d(tochange))",kel,jel);
      if(kel<0) continue;
      // save the faces in the element jel.
      memcpy(f_save,el2fnum[jel],c2s->nf*sizeof(INT));
      for (je=0;je<c2s->nf;je++){
	found=-1;
	k1=je*nvface;
  	for(j=0;j<nvface;j++){
  	  facej[j]=c2s->faces[je*nvface+j];
  	  facej[j]=g0->mnodes[jel*(nvcube+1)  + facej[j]];
  	}
	for (ke=0;ke<c2s->nf;ke++){
	  k2=ke*nvface;
	  for(i=0;i<nvface;i++){
	    facei[i]=c2s->faces[ke*nvface+i];
	    facei[i]=g0->mnodes[kel*(nvcube+1)  + facei[i]];
	  }
	  if(aresame(facei,facej,nvface)){
	    elneib[kel][ke]=jel;
	    elneib[jel][je]=kel;
	    //	    fprintf(stdout,"%%%%\nCommon face (ke=%d,je=%d):(elements %d-%d)",ke,je,kel,jel);fflush(stdout);
	    found=ke;
	    break;
	  }
	}
	if(found>=0)
	  break;
      }     
      if(found<0) {
	fprintf(stderr,"%%%%***ERROR: elements %d and %d are in the neighboring list but they do not share a face****\n\n",kel,jel);
	exit(12);
      }
      /* //SEE IF WE NEED TO run again around the faces in jel and rebuild el2fnum[jel][] */
      /* for (je=0;je<c2s->nf;je++){ */
      /* 	found=-1; */
      /* 	k1=je*nvface; */
      /* 	for(j=0;j<nvface;j++){ */
      /* 	  facej[j]=c2s->faces[je*nvface+j]; */
      /* 	  facej[j]=g0->mnodes[jel*(nvcube+1)  + facej[j]]; */
      /* 	} */
      /* 	for(ke=0;ke<c2s->nf;ke++){ */
      /* 	  k2=f_save[ke]; */
      /* 	  for(j=f2v->IA[k2];j<f2v->IA[k2+1];j++){ */
      /* 	    facei[j-f2v->IA[k2]]=f2v->JA[j]; */
      /* 	  } */
      /* 	  if(aresame(facei,facej,nvface)){	     */
      /* 	    found=ke; */
      /* 	    break; */
      /* 	  } */
      /* 	} */
      /* 	if(found<0){ */
      /* 	  fprintf(stderr,"%%%%***ERROR: FACE %d not found in macroelement %d****\n\n",k2,jel); */
      /* 	  exit(13); */
      /* 	} else { */
      /* 	  el2fnum[jel][je]=k2; */
      /* 	} */
      /* } */
      /* ////////////////////////////////////////////// */
      /* /\* print_full_mat_int(1,(nvcube+1),(g0->mnodes+kel*(nvcube+1)),"NODES33(KEL)"); *\/ */
      /* /\* print_full_mat_int(1,(nvcube+1),(g0->mnodes+jel*(nvcube+1)),"NODES33(JEL)"); *\/ */
      /* /////////////////////////////////////////////////////////////// */
    }    
  }
  // END REBUILD REBUILD
  /**** FINAL REORDER: make the vertices in shared faces the same!!! ***/
  for(lvl=0;lvl<mc->bfs->row;lvl++){
    j0=mc->bfs->IA[lvl];
    j1=mc->bfs->IA[lvl+1];
    for(kj=j0;kj<j1;kj++){
      jel=mc->bfs->JA[kj];
      kel=mc->etree[jel];// 
      //      if(kel<0){
      //	fprintf(stdout,"\n%%splitting element=%d",jel);
      //      } else {
      if(kel>=0){
  	je=locate0(jel,elneib[kel], c2s->nf);
  	ke=locate0(kel,elneib[jel], c2s->nf);
  	/*
	  in kel, we have the face je; in jel we have the face ke. we
	  want to make ke in jel same as je in kel; also the opposite
	  face needs to be reordered.
  	*/
	/* fprintf(stdout,"\nKEL=%d,JEL=%d,KE=%d; JE=%d;",kel,jel,ke,je);fflush(stdout); */
	/* print_full_mat_int(1,c2s->nf,elneib[kel],"Zelneib(kel)"); */
	/* print_full_mat_int(1,c2s->nf,elneib[jel],"Zelneib(jel)"); */
	/* print_full_mat_int(1,(nvface),(c2s->faces + je*nvface),"FACES(KEL)"); */
	/* print_full_mat_int(1,(nvcube+1),(g0->mnodes+kel*(nvcube+1)),"NODES(KEL)"); */
	/* print_full_mat_int(1,(nvface),(c2s->faces+ke*nvface),"FACES(JEL)"); */
	/* print_full_mat_int(1,(nvcube+1),(g0->mnodes+jel*(nvcube+1)),"NODES(JEL)"); */
  	for(i=0;i<nvface;i++){
  	  facei[i]=c2s->faces[ke*nvface+i];
  	  facei[i]=g0->mnodes[jel*(nvcube+1)  + facei[i]];
  	  mnodes[i]=c2s->faces[je*nvface+i];
  	  mnodes[i]=g0->mnodes[kel*(nvcube+1) + mnodes[i]];
  	}
  	k1=aresamep(mnodes,facei,nvface,p);
  	if(!k1){
  	  for(i=0;i<nvface;i++)p[i]=i;
  	}	      
	keswp=(ke+c2s->n)%c2s->nf;
	/*KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK*/
	/* fprintf(stdout,"\nKESWP=%d;k1=%d",keswp,k1);fflush(stdout); */
	/* print_full_mat_int(1,nvface,p,"permute"); */
	/* print_full_mat_int(1,nvface,facei,"facei(ke)"); */
	/* print_full_mat_int(1,nvface,mnodes,"facej(je)"); */
	/* print_full_mat_int(1,c2s->nf,elneib[kel],"elneib(kel)"); */
	/* print_full_mat_int(1,c2s->nf,el2fnum[kel],"el2fnum(kel)"); */
	/* print_full_mat_int(1,c2s->nf,elneib[jel],"elneib(jel)"); */
	/* print_full_mat_int(1,c2s->nf,el2fnum[jel],"el2fnum(jel)"); */
	/*KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK*/
  	for(i=0;i<nvface;i++){
	  //	  fprintf(stdout,"\n11111111:::%d,%d",facei[c2s->faces[ke*nvface+p[i]]],c2s->faces[ke*nvface+i]);
  	  facei[c2s->faces[ke*nvface+p[i]]]=c2s->faces[ke*nvface+i];
	  //  	  keswp=(ke+c2s->n)%c2s->nf; it is independent of the loop. 
	  //	  fprintf(stdout,"\n22222222:::%d,%d",facei[c2s->faces[keswp*nvface+p[i]]],c2s->faces[keswp*nvface+i]);
  	  facei[c2s->faces[keswp*nvface+p[i]]]=c2s->faces[keswp*nvface+i];
  	}
  	for(i=0;i<nvcube;i++)
  	  el2v->JA[i]=g0->mnodes[jel*(nvcube+1)+facei[i]];
  	for(i=0;i<nvcube;i++)
  	  g0->mnodes[jel*(nvcube+1)+i]=el2v->JA[i];
	/* print_full_mat_int(1,(nvface),(c2s->faces+ke*nvface),"FACES22(JEL)"); */
	/* print_full_mat_int(1,(nvcube+1),(g0->mnodes+jel*(nvcube+1)),"NODES22(JEL)"); */      
	///////////////////////////////AGAIN REBUILD AGAIN....      
	///////////////////////////////////////////////////////////////
	//	print_full_mat_int(1,(nvcube+1),(g0->mnodes+jel*(nvcube+1)),"NODES88(JEL)");
	///////////////////////////////////////////////////////////////
	memcpy(f_save,el2fnum[jel],c2s->nf*sizeof(INT));
	for (je=0;je<c2s->nf;je++){
	  found=-1;
	  k1=je*nvface;
	  for(j=0;j<nvface;j++){
	    facej[j]=c2s->faces[je*nvface+j];
	    facej[j]=g0->mnodes[jel*(nvcube+1)  + facej[j]];
	  }
	  //	print_full_mat_int(1,nvface,facej,"faceJJJJ33()");
	  //	print_full_mat_int(1,nvface,facei,"faceIIII33()");
	  for(ke=0;ke<c2s->nf;ke++){
	    k2=f_save[ke];
	    for(j=f2v->IA[k2];j<f2v->IA[k2+1];j++){
	      facei[j-f2v->IA[k2]]=f2v->JA[j];
	    }
	    if(aresame(facei,facej,nvface)){	    
	      //	      print_full_mat_int(1,nvface,facej,"facej33()");
	      //	      print_full_mat_int(1,nvface,facei,"facei33()");
	      found=ke;
	      break;
	    }
	  }
	  if(found<0){
	    fprintf(stderr,"%%%%***ERROR: FACE %d not found in macroelement %d****\n\n",k2,jel);
	    exit(13);
	  } else {
	    el2fnum[jel][je]=k2;
	    //	  if(el2fnum[jel][je]!=k2){
	    //	    fprintf(stdout,"\nel2fnum(%d,%d)=%d;",jel,je,el2fnum[jel][je]);fflush(stdout);	    
	    //	  }	  
	  }
	}	
      }      
    }          
  }
  /// AGAIN REBUILD::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  
  free(f_save);
  /// END REBUILD::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //  input_grid_print(g0);
 /*****************************************************/
  //YZYZYZYZ
  /* for(i=0;i<mc->nf;i++){ */
  /*   j1=f2v->IA[i+1]-f2v->IA[i]; */
  /*   //    memcpy(facei,(f2v->JA+f2v->IA[i]),j1*sizeof(INT)); */
  /*   fprintf(stdout,"code(%d)=%d; face_data(%d)=[",i,bcodesf[i],i); */
  /*   for (j1=f2v->IA[i];j1<f2v->IA[i+1]-1;++j1){ */
  /*     fprintf(stdout,"%d,",f2v->JA[j1]); */
  /*   } */
  /*   j1=f2v->IA[i+1]-1; */
  /*   fprintf(stdout,"%d]\n",f2v->JA[j1]); */
  /* } */
  // free
  icsr_free(el2v);free(el2v);
  icsr_free(f2v);free(f2v);
  icsr_free(el2el);free(el2el);
  free(facei);
  mc->elneib=elneib;
  mc->el2fnum=el2fnum;
  return mc;
}
/**********************************************************************/
/*!
 * \fn static void scomplex_merge1(const INT nvall, const INT nsall,
 * 		     macrocomplex *mc, scomplex **sc0, cube2simp *c2s)
 *
 * \brief combines an array of simplicial complexes together. REMOVES
 *        the common vertices.  Combines an array of simplicial
 *        complexes constructed from a macro complex together. the
 *        mc->iindex array should contain all global vertex numbers
 *        without repetitions.
 *
 * \param 
 *
 * \return
 *
 * \note
 *
 */
static void scomplex_merge1(const INT nvall,		\
		     const INT nsall,		\
		     macrocomplex *mc,		\
		     scomplex **sc0,		\
		     cube2simp *c2s)
{
  scomplex *sc=sc0[0];
  INT kel,i,ii,j,in1,iin1,newv,nnz;
  iCSRmat bndry_v1;//,bndry_v2;
  INT n1=(sc->n+1),ns0,nv=nvall,ns=nsall;
  if(mc->nel!=1) {
    for(kel=0;kel<mc->nel;++kel){
      sc0[kel]->bndry_v->col=nvall;
    }
    sc->marked=realloc(sc->marked,ns*sizeof(INT));
    sc->gen=realloc(sc->gen,ns*sizeof(INT));
    sc->nbr=realloc(sc->nbr,ns*n1*sizeof(INT));
    sc->parent=realloc(sc->parent,ns*sizeof(INT));
    sc->child0=realloc(sc->child0,ns*sizeof(INT));
    sc->childn=realloc(sc->childn,ns*sizeof(INT));
    sc->nodes=realloc(sc->nodes,ns*n1*sizeof(INT));
    sc->bndry=realloc(sc->bndry,nv*sizeof(INT));
    /*  
	every vertex can be on at most "dim" boundaries in one macroelement. if every
	boundary on every macroelement has a different code, then these
	are at most mc->nel*2*dim different codes as every macroelement
	has at most 2*dim faces. 
    */
    //  icsr_realloc(nv,mc->nel*2*sc->n,nv*sc->n,sc->bndry_v); // no need of this here. 
    /* 
     * coord sys: 1 is polar, 2 is cyl and so on: not fully implemented
     * and tested yet
     */
    /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
    /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
    sc->csys=realloc(sc->csys,nv*sizeof(INT));
    /*connected components*/
    sc->cc=mc->cc;sc->bndry_cc=mc->bndry_cc;
    sc->flags=(INT *)realloc(sc->flags,ns*sizeof(INT));
    sc->x=(REAL *)realloc(sc->x,nv*(sc->n)*sizeof(REAL));
    sc->vols=(REAL *)realloc(sc->vols,ns*sizeof(REAL));
    //  sc->fval=(REAL *)realloc(sc->fval,nv*sizeof(REAL)); // function values at every vertex; not used in general;
    //
    /* for(i=0;i<sc->bndry_v->row;++i){ */
    /*   if((sc->bndry_v->IA[i+1]-sc->bndry_v->IA[i])){ */
    /*     fprintf(stdout,"\n0size(row=%d)=%d; 0entries=[ ",i,sc->bndry_v->IA[i+1]-sc->bndry_v->IA[i]); */
    /*     for(j=sc->bndry_v->IA[i];j<sc->bndry_v->IA[i+1];++j){ */
    /* 	//	fprintf(stdout,"%d(Xc=%d,Xb=%d) ",sc->bndry_v->JA[j],sc->bndry_v->val[j],sc->bndry_v->val[nnz+j]); */
    /* 	fprintf(stdout,"%d(0c=%d) ",sc->bndry_v->JA[j],sc->bndry_v->val[j]); */
    /*     }       */
    /*     fprintf(stdout,"]"); fflush(stdout); */
    /*   } */
    /* }   */
    for(kel=1;kel<mc->nel;kel++){
      ns0=sc->ns;
      /* fprintf(stdout,"\nElement=%d; ns0=%d",kel,ns0);       */
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
	sc->bndry[i]=sc0[kel]->bndry[ii];
	sc->csys[i]=sc0[kel]->csys[ii];
	//      sc->fval[i]=sc0[kel]->fval[ii];
	in1=i*sc->n;
	iin1=ii*sc->n;
	for(j=0;j<sc->n;j++)
	  sc->x[in1+j]=sc0[kel]->x[iin1+j];
      }
      nnz=sc0[kel]->bndry_v->nnz;
      for(i=0;i<sc0[kel]->bndry_v->row;++i){
	for(j=sc0[kel]->bndry_v->IA[i];j<sc0[kel]->bndry_v->IA[i+1];++j){
	  ii=sc0[kel]->bndry_v->JA[j];// this is vertex number;
	  //	fprintf(stdout,"\nkel=%d:ii=%d,newindex=%d",kel,ii,mc->iindex[kel][ii]);fflush(stdout);
	  sc0[kel]->bndry_v->JA[j]=mc->iindex[kel][ii];
	}
      }
      bndry_v1=icsr_create(sc->bndry_v->row,sc->bndry_v->col,sc->bndry_v->nnz);
      memcpy(bndry_v1.IA,sc->bndry_v->IA,(bndry_v1.row+1)*sizeof(INT));
      memcpy(bndry_v1.JA,sc->bndry_v->JA,bndry_v1.nnz*sizeof(INT));
      memcpy(bndry_v1.val,sc->bndry_v->val,bndry_v1.nnz*sizeof(INT));
      nnz=sc->bndry_v->nnz;
      /*** MOVE POINTERS AND ADD BNDRY CODES ***/
      // free, and use as adding
      icsr_free(sc->bndry_v);
      //add once
      icsr_add(&bndry_v1,sc0[kel]->bndry_v,sc->bndry_v); //
      //    sc->bndry_v->val=realloc(sc->bndry_v->val,2*sc->bndry_v->nnz*sizeof(INT));
      /*END MOVE POINTERS AND ADD BNDRY CODES*/
      sc->ns+=sc0[kel]->ns;
      haz_scomplex_free(sc0[kel]);
      // very wasteful
      //    free(bndry_v2.val);
      icsr_free(&bndry_v1);
    }
    ///////////////////////////////////////////
    sc->nv=nvall;
  }
  ///////////////////////////////////////////  
  /* for(i=0;i<sc->bndry_v->row;++i){ */
  /*   if((sc->bndry_v->IA[i+1]-sc->bndry_v->IA[i])){ */
  /*     fprintf(stdout,"\nXsize(row=%d)=%d; Xentries=[ ",i,sc->bndry_v->IA[i+1]-sc->bndry_v->IA[i]); */
  /*     for(j=sc->bndry_v->IA[i];j<sc->bndry_v->IA[i+1];++j){ */
  /* 	//	fprintf(stdout,"%d(Xc=%d,Xb=%d) ",sc->bndry_v->JA[j],sc->bndry_v->val[j],sc->bndry_v->val[nnz+j]); */
  /* 	fprintf(stdout,"%d(Xc=%d) ",sc->bndry_v->JA[j],sc->bndry_v->val[j]); */
  /*     } */
  /*     fprintf(stdout,"]"); fflush(stdout); */
  /*   } */
  /* } */
  ////////////////////////////////////////////////////////////////////////////
  /*Now we transpose to obtain the vertex->face correspondence*/
  bndry_v1=icsr_create(sc->bndry_v->row,sc->bndry_v->col,sc->bndry_v->nnz);
  memcpy(bndry_v1.IA,sc->bndry_v->IA,(bndry_v1.row+1)*sizeof(INT));
  memcpy(bndry_v1.JA,sc->bndry_v->JA,bndry_v1.nnz*sizeof(INT));
  memcpy(bndry_v1.val,sc->bndry_v->val,bndry_v1.nnz*sizeof(INT));
  icsr_free(sc->bndry_v);
  icsr_trans(&bndry_v1,sc->bndry_v);
  icsr_free(&bndry_v1);  
  nnz=sc->bndry_v->nnz;
  sc->bndry_v->val=realloc(sc->bndry_v->val,2*nnz*sizeof(INT));
  for(i=0;i<sc->bndry_v->row;++i){
    if((sc->bndry_v->IA[i+1]-sc->bndry_v->IA[i])){
      for(ii=sc->bndry_v->IA[i];ii<sc->bndry_v->IA[i+1];++ii){
	j=sc->bndry_v->JA[ii];	
	//sc->bndry_v->val[ii] must be equal to the mc->bcodesf[j].
	if(j<mc->nf && j>=0){	    
	  sc->bndry_v->val[ii+nnz]=mc->isbface[j];	  
	}	
      }
    }
  }
  return;
}
/**********************************************************************/
/*!
 * \fn INT locate1(INT *b,INT *a,INT n,INT *a2,INT n2,INT m2) 
 *
 * \brief locates the elements of an array a in a two dimensional
 *     array a2. a has n elements and a2 is n2 by m2.  returns the the
 *     number of the rows(a2) that contain ALL elements of a; and b
 *     contains the indices of these rows; in case there is an element
 *     of a not contained in any row returns 0; in case no row
 *     contains ALL elements, returns 0; b should have size
 *     n2*sizeof(INT);
 *
 * \param 
 *
 * \return
 *
 * \note
 *
 */
INT locate1(INT *b,INT *a,INT n,INT *a2,INT n2,INT m2)
{  
  INT i,j,im,nb,aj,bi;
  for(i=0;i<n2;i++) b[i]=i;
  for(j=0;j<n;j++){    
    aj=a[j];
    nb=0;
    for(i=0;i<n2;i++){
      if(b[i]<0) continue;
      im=i*m2;
      bi=locate0(aj,(a2+im),m2);
      if(bi<0){b[i]=bi;}
      else {b[i]=i;nb++;}
    }
    if(!nb) return 0;// one of the elements was not found
  }
  nb=0;
  for(i=0;i<n2;i++){
    if(b[i]<0) continue;
    b[nb]=b[i];
    nb++;
  }
  return nb;
}
/**********************************************************************/
/*!
 * \fn void macrocomplex_free(macrocomplex *mc) 
 *
 * \brief free the structure macroelements. 
 *
 * \param 
 *
 * \return
 *
 * \note
 *
 */
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
  if(mc->bfs) {
    icsr_free(mc->bfs);free(mc->bfs);
  }
  if(mc->fullel2el) {
    icsr_free(mc->fullel2el);free(mc->fullel2el);
  }
  if(mc) free(mc);
  return;
}
/**********************************************************************/
/*!
 * \fn void fix_grid(macrocomplex *mc, scomplex **scin,cube2simp *c2s,
 *                   input_grid *g0)
 *
 * \brief A loop over all macroelements and generates array
 *        mc->iindex[kel][] for kel=1:nel (number of
 *        macroelements=nel). Each of the entries in iindex gives the
 *        global number of a local vertex.
 *
 * \param scin is an array of simplicial complexes each generated
 *        independently in a macroelement. this program combines them
 *
 * \return
 *
 * \note
 *
 */
void fix_grid(macrocomplex *mc,		\
	       scomplex **scin,			\
	       cube2simp *c2s,			\
	       input_grid *g0)
{
  INT nsall,nvall;
  if(mc->nel<=1){
    nvall=scin[0]->nv;
    nsall=scin[0]->ns;
    scomplex_merge1(nvall,nsall,mc,scin,c2s);
    return;
  }
  scomplex *scp;
  INT dim=c2s->n,dim1=c2s->n+1,nvface=c2s->nvface,nvcube=c2s->nvcube;
// at most dim (n-1)dimensional faces may intersect to form a vertex
  INT *mp = (INT *)calloc(dim1,sizeof(INT));
  INT *m = (INT *)calloc(dim1,sizeof(INT));
  INT *mwrk = (INT *)calloc(dim1,sizeof(INT));
  INT *mi = (INT *)calloc(2*dim1,sizeof(INT));
  INT *mip = (INT *)calloc(2*dim1,sizeof(INT));
  INT *vertk = (INT *)calloc(dim1,sizeof(INT));
  INT *vertj = (INT *)calloc(dim1,sizeof(INT));
  INT *facesk = (INT *)calloc(2*dim*nvface,sizeof(INT));
  INT *facesj = (INT *)calloc(2*dim1*nvface,sizeof(INT));
  INT *pd=(INT *)calloc(dim,sizeof(INT));
  INT *ti=vertj,*tip=vertk;
  //scalars
  INT nv=-10,kel, jel,nk,nj,flag;
  INT numv,kf,kfp,kz,kdim,kj;
  INT i,j,k,iaa, iab;
  INT *nodesj=NULL,*nodesk=NULL;
  // place holders (short hand)
  //  iCSRmat *fel2el=mc->fullel2el;
  //  INT **nd=mc->nd;
  //  INT **elneib=mc->elneib;
  //  INT **el2fnum=mc->el2fnum;
  INT *iindex,*iindexp;
  //  for(kz=(dim-1);kz>=0;kz--){
  //serching for intersection of kz-dimensional faces.
  //    kdim=(1<<kz);
  //  for(knnz=0;knnz<mc->bfs->nnz;knnz++){    
  //    kel=mc->bfs->JA[knnz];
  INT neg,nvold;

  nvall=0;nsall=0;

  //  print_full_mat_int(g0->nel,c2s->nvcube+1,g0->mnodes,"mel0");
  
  for(kel=0;kel<mc->nel;kel++){
    // we have not been here, so let us set the initial indexing to be the original indexing; 
    nvold=nvall;
    for(i=0;i<scin[kel]->nv;i++){
      mc->iindex[kel][i]=i + nvall;// this is the global number if there are no removals of vertices. 
    }
    neg=0;
    iaa=mc->fullel2el->IA[kel];
    iab=mc->fullel2el->IA[kel+1];
    if((iab-iaa)<=0){
      if(g0->print_level>5){
	fprintf(stdout,"\nin %s: macroelement=%d; vertices=%d; overlaps=%d;",__FUNCTION__,kel,scin[kel]->nv,neg);fflush(stdout);
      }
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
      jel=mc->fullel2el->JA[kj];
      nodesj=(g0->mnodes+jel*(nvcube+1));
      for(i=0;i<c2s->nf;i++){
	for(j=0;j<nvface;j++){
	  k=c2s->faces[i*nvface+j];
	  facesj[i*nvface+j]=nodesj[k];
	}
      }
      numv=0;
      // find now how many times these two elements intersect:
      kdim=mc->fullel2el->val[kj];
      // calculate the dimension of the intersection:
      kz=((INT )floor(log2(((REAL) kdim)-1e-3)))+1;
      iindex=mc->iindex[jel];
      for(k=0;k<c2s->nvcube;k++){
	i=nodesk[k];
	i=locate0(i,nodesj,c2s->nvcube);
	if(i<0) continue;
	vertk[numv]=nodesk[k];
	vertj[numv]=nodesj[i];
	numv++;
      }
      nk=locate1(mip,vertk,numv,facesk,c2s->nf,c2s->nvface);
      nj=locate1(mi,vertj,numv,facesj,c2s->nf,c2s->nvface);
      if((nk == nj) && nk==(dim-kz) && nj==(dim-kz)) {
	align_lattice(pd,numv,nodesj,nodesk,c2s);
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
	    ti[j]=mc->nd[jel][mi[j]];
	  }
	  if(mip[j]<dim){
	    mip[j]=dim-(mip[j]+1);
	    tip[j]=0;
	  } else {
	    mip[j]=dim-((mip[j]%dim)+1);
	    tip[j]=mc->nd[kel][mip[j]];
	  }	    
	}
	nv=0;
	for(kfp=0;kfp<scp->nv;kfp++){
	  coord_lattice(mp,dim,kfp,scp->nv,mc->nd[kel]);
	  flag=1;
	  for(j=0;j<nj;j++) if(mp[mip[j]]!=tip[j]){flag=0;break;}
	  if(flag){
	    // find out which (i,j) is (ip,jp) in kel
	    memcpy(mwrk,mp,dim*sizeof(INT)); //	    
	    // put everything in sync:
	    for(k=0;k<dim;k++){
	      if(pd[k]<0){
		mwrk[abs(pd[k])-1]=mc->nd[jel][k]-mwrk[abs(pd[k])-1];
		m[k]=mwrk[abs(pd[k])-1]; // added later??? why???
	      }
	    }
	    for(k=0;k<dim;k++){
	      if(pd[k]>0){
		m[k]=mwrk[abs(pd[k])-1];
	      }
	    }
	    for(k=0;k<dim;k++){
	      if(pd[k]==0){
		m[k]=mwrk[abs(pd[k])];
	      }
	    }
	    for(k=0;k<nk;k++) m[mi[k]]=ti[k];
	    kf=num_lattice(m,dim,mc->nd[jel]);
	    if(iindexp[kfp]>=nvold){
	      iindexp[kfp]=iindex[kf];
	      neg++;
	    }
	  } else {
	    if(iindexp[kfp]>=nvold) {
	      iindexp[kfp]=nv+nvold;
	      nv++;
	    }
	  }
	}
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/	
      } else {
	fprintf(stderr,"\n*** An error in counting overlaps in %s ***\n",__FUNCTION__);
	exit(128);
      }
    }
    if(g0->print_level>5)
      fprintf(stdout,"\nin %s: macroelement=%d; total=%d; overlaps=%d",__FUNCTION__,kel,scp->nv,neg);
    nvall+=nv;
    nsall+=scin[kel]->ns;
    //fprintf(stdout,"\nGLOBALLY:v_total=%d; s_total=%d",nvall,nsall);fflush(stdout);
  }
  /* INT nnz,iloc,iglob; */
  /* for(jel=0;jel<mc->nel;++jel){ */
  /*   fprintf(stdout,"\n(Again)Element=%d;",jel); */
  /*   for(kf=0;kf<scin[jel]->nv;++kf){ */
  /*     fprintf(stdout,"\noldv=%d-->newv=%d",kf,mc->iindex[jel][kf]); */
  /*   } */
  /*   nnz=scin[jel]->bndry_v->nnz; */
  /*   for(kf=0;kf<scin[jel]->bndry_v->row;++kf){ */
  /*     if((scin[jel]->bndry_v->IA[kf+1]-scin[jel]->bndry_v->IA[kf])){ */
  /* 	fprintf(stdout,"\nTsize(Trow=%d)=%d; Tentries=[ ",kf,scin[jel]->bndry_v->IA[kf+1]-scin[jel]->bndry_v->IA[kf]); */
  /* 	for(j=scin[jel]->bndry_v->IA[kf];j<scin[jel]->bndry_v->IA[kf+1];++j){ */
  /* 	  iloc=scin[jel]->bndry_v->JA[j]; */
  /* 	  iglob=mc->iindex[jel][iloc]; */
  /* 	  fprintf(stdout,"%d(Tc=%d,Tb=%d) ",iglob,scin[jel]->bndry_v->val[j],scin[jel]->bndry_v->val[nnz+j]); */
  /* 	} */
  /* 	fprintf(stdout,"]"); fflush(stdout); */
  /* 	fprintf(stdout,"\nZsize(Zrow=%d)=%d; Zentries=[ ",kf,scin[jel]->bndry_v->IA[kf+1]-scin[jel]->bndry_v->IA[kf]); */
  /* 	for(j=scin[jel]->bndry_v->IA[kf];j<scin[jel]->bndry_v->IA[kf+1];++j){ */
  /* 	  fprintf(stdout,"%d(Zc=%d,Zb=%d) ",scin[jel]->bndry_v->JA[j],scin[jel]->bndry_v->val[j],scin[jel]->bndry_v->val[nnz+j]); */
  /* 	} */
  /* 	fprintf(stdout,"]"); fflush(stdout); */
  /*     } */
  /*   } */
  /*   haz_scomplex_print(scin[jel],0,"JEL"); */
  /* }   */
  if(g0->print_level>5){
    for(kel=0;kel<mc->nel;kel++){
      fprintf(stdout,"\nelement{%d}=[",kel);
      for (i = 0;i<(nvcube+1);i++){
	fprintf(stdout,"%d ",g0->mnodes[kel*(nvcube+1)+i]);
      }
      fprintf(stdout,"];");
      for(i=0;i<scin[kel]->nv;i++){
	coord_lattice(m,dim,i,scin[kel]->nv,mc->nd[kel]);
	fprintf(stdout,"\n     [%d=( ",i);
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
  free(mwrk);
  free(mp);
  free(mi);
  free(mip);
  free(facesk);
  free(facesj);
  free(vertk);
  free(vertj);
  free(pd);
  scomplex_merge1(nvall,nsall,mc,scin,c2s);
  return;
}
/**********************************************************************/
/*!
 * \fn scomplex *generate_initial_grid(input_grid *g0) 
 *
 * \brief From the data from input_grid g0 generates the global
 *        simplicial complex with the grid based on the macroelements
 *        and the divisions given by *g0.
 *
 * \param 
 *
 * \return
 *
 * \note
 *
 */
scomplex **generate_initial_grid(input_grid *g0)
{
  /* 
     From an input grid loops over the macroelements and generates
     meshes depending on the number of divisions  in every dimension. 
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
  //  ilexsort(g0->nel,(c2s->nvcube+1),g0->mnodes,p);
  /* for(kel=0;kel<g0->nel;kel++){ */
  /*   print_full_mat_int(1,c2s->nvcube+1,g0->mnodes+kel*(c2s->nvcube+1),"mnodes33"); */
  /* } */
  //  input_grid_print(g0);
  /*-------------------------------------------------------------------*/
  //  INT *efound=calloc(c2s->ne*(g0->nel),sizeof(INT));  
  g=malloc(sizeof(input_grid));  //temp grid for one macroelement// we need to free this at the end
  /**/
  g->title=strndup(g0->title,(strlen(g0->title)+1)*sizeof(char));
  g->fgrid=strndup(g0->fgrid,(strlen(g0->fgrid)+1)*sizeof(char));
  g->fvtu=strndup(g0->fvtu,(strlen(g0->fvtu)+1)*sizeof(char));
  /**/
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
  g->num_refine_points=g0->num_refine_points;
  input_grid_arrays(g);
  /* copy these as these are the same as g0 */
  memcpy(g->systypes,g0->systypes,abs(g0->ncsys)*sizeof(INT));
  memcpy(g->syslabels,g0->syslabels,abs(g0->ncsys)*sizeof(INT));
  memcpy(g->ox,g0->ox,g0->dim*abs(g0->ncsys)*sizeof(REAL));
  INT chng=1,iter=0,maxiter=1024;  
  // get the macroelement mesh in a structure
  macrocomplex *mc=set_mmesh(g0,c2s,p);
  /****************************************************************/
  set_edges(g0,c2s);  
  while(chng && (iter<maxiter)){
    iter++;
    // make the divisions in g0->seg consistent;
    chng=set_ndiv_edges(g,g0,c2s,mc->nd,iter);
  }
  /*PLACE HOLDERS*/
  INT **el2fnum=mc->el2fnum;
  INT *isbface=mc->isbface;
  INT *bcodesf=mc->bcodesf;
  iCSRmat *bfs0=mc->bfs;
  if(g0->print_level>0){
    fprintf(stdout,"\n%%DFS(domains): %d connected components",mc->cc);
    fprintf(stdout,"\n%%DFS(boundaries): %d connected components",mc->bndry_cc);
  }else if (g0->print_level>3){
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
    /* fprintf(stdout,"\nbfs0=["); */
    /* icsr_print_matlab_val(stdout,bfs0); */
    /* fprintf(stdout,"];"); */
    /* fprintf(stdout,"\nbfs=sparse(bfs0(:,1),bfs0(:,2),bfs0(:,3));\n"); */
    /* fprintf(stdout,"\n%%%%*****   g0_nf=%d; mc_nf=%d\n",g0->nf,mc->nf); */
    /* fprintf(stdout,"\nfel2el0=["); */
    /* icsr_print_matlab_val(stdout,mc->fullel2el); */
    /* fprintf(stdout,"];"); */
    /* fprintf(stdout,"\nfel2el=sparse(fel2el0(:,1),fel2el0(:,2),fel2el0(:,3));\n"); */
    /* print_full_mat_int(1,mc->nf,mc->isbface,"isbface"); */
    /* print_full_mat_int(1,mc->nf,mc->bcodesf,"bcodesf"); */
  /* set the divisions on every edge now; since they are consistent we have: */
  if(set_ndiv_edges(g,g0,c2s,mc->nd,0)) {
    fprintf(stderr,"\n\n***ERR in %s: the divisions of the edges cannod be inconsistent during second call of set_ndiv_edges()\n\n",__FUNCTION__);
    exit(4);
  }
  for(kel=0;kel<mc->nel;kel++){
    for(i=0;i<c2s->n;i++){
      if(mc->nd[kel][i]<=0) mc->nd[kel][i]=1;
    }
  }
  if(g0->print_level>15) input_grid_print(g0);
  INT nsall,nvall,nnz;
  nsall=0;nvall=0;
  /* now mc->nd is known, let us allocate iindex */
  for(kel=0;kel<g0->nel;kel++){
    nvall=1;
    for(i=0;i<g0->dim;i++)
      nvall*=(mc->nd[kel][i]+1);
    mc->iindex[kel]=calloc(nvall,sizeof(INT));
    // TRUE    mc->flags[kel]=g0->mnodes[nvcube+kel*(nvcube+1)];
    mc->flags[kel]=2*kel+1;
  }
  INT intype=g0->ref_type-1,kj,jel;  
  INT *codef=calloc(c2s->nf,sizeof(INT));
  INT *labelf=calloc(c2s->nf,sizeof(INT));
  INT *isbndf=calloc(c2s->nf,sizeof(INT));
  iCSRmat bndry_v1,bndry_v2;// local vertex/bface relation. Later combined in sc->bndry_v;
  INT *tmp_ptr; // to store isbface
  INT kjj;
  for(kj=0;kj<bfs0->row;kj++){
    if(g0->ref_type>=0) intype++;
    else intype=-1;
    for(kjj=bfs0->IA[kj];kjj<bfs0->IA[kj+1];kjj++){
      jel=bfs0->JA[kjj];    
      /* fprintf(stdout,"\nAAAElement=%d;",jel);       */
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
	labelf[i]=el2fnum[jel][i];  // this is the face global number;
	codef[i]=bcodesf[labelf[i]];// code associated with this macroelement face
	isbndf[i]=isbface[labelf[i]];// is this macroelement face on the boundary
      }
      sc[jel]=umesh(g->dim,mc->nd[jel],c2s,		\
		    labelf,isbndf,codef,mc->flags[jel],	\
		    intype);
      // now we make the boundary matrix global.... transpose it so it is "face"->"vertex"
      sc[jel]->bndry_v->col=mc->nf;//
      nnz=sc[jel]->bndry_v->nnz;
      icsr_trans(sc[jel]->bndry_v,&bndry_v1);
      tmp_ptr=sc[jel]->bndry_v->val;
      sc[jel]->bndry_v->val += nnz;
      icsr_trans(sc[jel]->bndry_v,&bndry_v2);
      sc[jel]->bndry_v->val = tmp_ptr;
      nnz=sc[jel]->bndry_v->nnz;
      /* Now we copy the face->vertex correspondence over sc->bndry_v */
      sc[jel]->bndry_v->col=bndry_v1.col;//col
      sc[jel]->bndry_v->row=bndry_v1.row;//row
      sc[jel]->bndry_v->nnz=bndry_v1.nnz;//nnz
      sc[jel]->bndry_v->IA=realloc(sc[jel]->bndry_v->IA,(sc[jel]->bndry_v->row+1)*sizeof(INT));      
      memcpy(sc[jel]->bndry_v->IA,bndry_v1.IA,(bndry_v1.row+1)*sizeof(INT));
      memcpy(sc[jel]->bndry_v->JA,bndry_v1.JA,bndry_v1.nnz*sizeof(INT));
      memcpy(sc[jel]->bndry_v->val,bndry_v1.val,bndry_v1.nnz*sizeof(INT));
      memcpy((sc[jel]->bndry_v->val+nnz),bndry_v2.val,bndry_v2.nnz*sizeof(INT));
      icsr_free(&bndry_v1);
      icsr_free(&bndry_v2);
      //////////////////////////////////
      nsall+=sc[jel]->ns;
      nvall+=sc[jel]->nv;
      //    if((mc->fullel2el->IA[jel+1]-mc->fullel2el->IA[jel])<=0){
      for(i=0;i<sc[jel]->nv;i++){
	mc->iindex[jel][i]=i;
      }
      //}
      //      fprintf(stdout,"\n%%Mapping back to the macroelement...\n");
      map2mac(sc[jel],c2s,g);
      //      fprintf(stdout,"\n%%in %s: mesh(macroelement=%d): nv=%d; nsimp=%d",__FUNCTION__,jel,sc[jel]->nv,sc[jel]->ns); fflush(stdout);
    }
  }
  fix_grid(mc,sc,c2s,g0);
  /*initialize the parent_v matrix*/
  icsr_free(sc[0]->parent_v);
  sc[0]->parent_v[0]=icsr_create(sc[0]->nv,sc[0]->nv,sc[0]->nv);
  sc[0]->parent_v->IA[0]=0;
  for(i=0;i<sc[0]->parent_v->row;++i){
    sc[0]->parent_v->JA[sc[0]->parent_v->IA[i]]=i;
    sc[0]->parent_v->IA[i+1]=i+1;
  }  
  if(g0->print_level>4){
    fprintf(stdout,"\n%%merged(macroelements=[%d..%d]): vertices=%d; simplices=%d", \
	    0,mc->nel-1,sc[0]->nv,sc[0]->ns);  
    fprintf(stdout," ..done.\n\n");fflush(stdout);
  }
  free(p);
  free(isbndf);
  free(codef);
  free(labelf);
  input_grid_free(g);
  macrocomplex_free(mc);
  cube2simp_free(c2s);
  /* order simplex2vertex array as required by the adaptive refinement */
  find_nbr(sc[0]->ns,sc[0]->nv,sc[0]->n,sc[0]->nodes,sc[0]->nbr);
  //
  INT *wrk1=calloc(5*(sc[0]->n+2),sizeof(INT));
  /*   construct bfs tree for the dual graph*/
  sc_vols(sc[0]);
  abfstree(0,sc[0],wrk1,g0->print_level);
  free(wrk1);
  //  sc=realloc(sc,sizeof(scomplex *));
  return sc;  
}
/*EOF*/
