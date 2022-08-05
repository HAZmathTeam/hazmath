/*! \file src/graphs/dfs.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note depth first search and related routines. 
 *
 */
#include "hazmath.h"
/***********************************************************************************************/
/*!
 * \fn iCSRmat **lex_bfs(INT n,INT *ia, INT *ja, ivector *anc, const REAL *guide)
 *
 * \brief lexicographical breath first search
 *
 * \param n                number of vertices
 * \param (ia,ja):         adjacency structure of the graph; diagonels are allowed
 * \param anc[n]:           a vector containing the ancestors of the vertices:
 *                         v<->anc[v] is an edge in the bfs tree. 
 *
 * \param guide[n]         an array to identify the "initial points"
 *                         for BFS. In every connected component this 
 *                         will pick as roots the minimum and the maximum
 *                         of "guide". If guide is NULL than the root
 *                         in a connected component is picked
 *                         randomly.
 *
 * \return pointer to iCSRmat of size n,n with n nonzeroes, which
 *         contains the bfs ordering permutation of the vertices
 *         (component after component).
 *
 * \note Following the algol-like pseudocode from:
 *  Rose, Donald J. ; Tarjan, R. Endre ; Lueker, George S.
 *  Algorithmic aspects of vertex elimination on graphs.  SIAM
 *  J. Comput. 5 (1976), no. 2, 266–283 (MR0408312).
 *  https://doi.org/10.1137/0205021,
 *
 * Created by ltz1 on 20190621.
 */
iCSRmat **lex_bfs(INT n,INT *ia, INT *ja,ivector *inv_ord,ivector *anc, const REAL *guide)
{
  INT c,h,v,w,p;  
  INT ibeg,iend,i,iv,nfx,ijk,i0;
  iCSRmat **blk;
  /* 
   *  two lists: queue and sets of vertices. the queue is described by
   *  head[] and backp[] and the sets are described by next[] and
   *  back[]. Each cell "c" is either a header cell or describes a
   *  set of vertices that have the same label.
   *
  */
  /*----------------------------------------------------------------*/
  // find all connected components:
  blk=(iCSRmat **)malloc(2*sizeof(iCSRmat *));
  blk[0]=run_dfs(n,ia,ja);
  /* Move to first place in every connected component, the max abs(guide)*/
  if(guide!=NULL){
    INT iaa,iab;
    REAL emax, etest;
    for(i=0;i<blk[0]->row;++i){
      iaa=blk[0]->IA[i];
      iab=blk[0]->IA[i+1];
      if((iab-iaa)<2) continue;
      emax=fabs(guide[blk[0]->JA[iaa]]);
      iv=iaa;
      iaa++;
      for(ijk=iaa;ijk<iab;++ijk){
	etest=fabs(guide[blk[0]->JA[ijk]]);
	if(emax < etest){
	  emax = etest;
	  iv=ijk;
	}
      }
      // now swap the first vertex in this connected component with the emax one.
      iaa--;
      fprintf(stdout,"\nconnected_component=%d;max_node=%d,emax=%e,iaa=%d\n",i,blk[0]->JA[iv],emax,iaa);
      if(iv==iaa) continue;
      ijk=blk[0]->JA[iaa];
      blk[0]->JA[iaa] = blk[0]->JA[iv];
      blk[0]->JA[iv]=ijk;
    }
  }
  blk[1]=malloc(sizeof(iCSRmat));
  *blk[1]=icsr_create(blk[0]->row,blk[0]->col,blk[0]->nnz);
  /*----------------------------------------------------------------*/  
  memcpy(blk[1]->JA,blk[0]->JA,blk[0]->nnz*sizeof(INT));
  INT *iord=blk[1]->JA;// this makes ordering following connected components.
  INT *level=blk[1]->val;// this is the level (distance) from every 
  INT *tree=anc->val;
  INT *iord1=inv_ord->val;
  //////////////////////////////////////////////////////
  INT nmem=n;
  // allocate memory;
  INT *cell=calloc(n,sizeof(INT));
  for (i=0;i<n;i++) cell[i]=-1;
  INT *flag=calloc(nmem,sizeof(INT));
  INT *head=calloc(nmem,sizeof(INT));
  INT *next=calloc(nmem,sizeof(INT));
  INT *back=calloc(nmem,sizeof(INT));
  INT *fixlst=calloc(nmem,sizeof(INT));
  for (i=0;i<nmem;i++){
    fixlst[i]=-1;
    head[i]=-1;
    back[i]=-1;
    next[i]=-1;
    flag[i]=-1;
  }
  /*     assign empty label to all vertices */
  head[0] = 1; // the first cell is the head of the queue; it does not
	       // describe any set of vertices, but points to the
	       // cell=1 which will describe the first set of
	       // vertices.
  back[1] = 0; // the second cell points to the first cell as its predecessor
  /**/
  head[1] = -1;   
  back[0] = -1;
  next[0] = -1;
  flag[0] = -1;
  flag[1] = -1;
  // cell 0 is the beginning; cell 1 is header cell for the first set;
  c=2; // thus the first empty cell is cell 2. 
  for(iv=0;iv<n;iv++){
    v = iord[iv]; // vertex label
    head[c] = v;  //head(cell)=element;
    cell[v] = c;  //cell(vertex)=cell.
    next[c-1] = c;
    flag[c] = 1;
    back[c] = c-1;
    c++;
    if(c>=nmem){
      nmem=c+1;
      flag=realloc(flag,nmem*sizeof(INT));
      head=realloc(head,nmem*sizeof(INT));
      next=realloc(next,nmem*sizeof(INT));
      back=realloc(back,nmem*sizeof(INT));
      head[c]=-1; back[c]=-1; next[c]=-1; flag[c]=-1;
    }
    iord1[v] = -1;
    level[v] = -1;
    tree[v] = -1;
  }
  next[c-1] = -1;
  for (i = n-1;i>=0;i--){
    /*********************************************************************/
    nfx=0;
    //    print_cells("%%%%first:",nfx,i,cell,flag,head,next,back,fixlst);
    //    if(i<n-10) break;
    //C  Skip empty sets
    while(next[head[0]] < 0){
      head[0] = head[head[0]];
      back[head[0]]=0;
    }
    //C  pick the cell for the next vertex to number
    p=next[head[0]];
    //C     C ADDITION BY ltz
    v = head[p];
    //C  delete cell of vertex from set
    next[head[0]] = next[p];
    next[back[cell[v]]]=next[cell[v]];
    if(next[cell[v]]>=0){
      back[next[cell[v]]]=back[cell[v]];
    }
    //C assign number i to v
    iord[i] = v;
    iord1[v] = i;
    if(level[v]<0) level[v]=0;
    nfx=0;
    ibeg = ia[v]-1;
    iend = ia[v+1]-1;
    for(ijk = iend;ijk>ibeg;ijk--){
      w = ja[ijk];
      if((tree[w]<0)&&(level[w]<0)) tree[w]=v;
      if(level[w]<0) level[w]=level[v]+1;
      // if vertex is not numbered:
      if(iord1[w] < 0) {
	// delete cell of w from the set (this will go into a new set). 
	next[back[cell[w]]]=next[cell[w]]; // skip cell[w] in the element list
	if(next[cell[w]] >=0){
	  back[next[cell[w]]]=back[cell[w]]; // if there was a cell
					     // pointed to as next by
					     // cell(w), then put its
					     // prev. cell to be the
					     // prev(cell(w)) and not
					     // cell(w);
	}
	h = back[flag[cell[w]]];  // this is the header set cell which
				  // is a predecessor of the set header cell
				  // containing w
	// if h is an old set, then we create a new set (c=c+1)
	if(flag[h]<0) {//
	  head[c] = head[h];// insert c between h and head[h]
	  head[h] = c;
	  back[head[c]] = c;
	  back[c] = h;
	  flag[c] = 0;
	  next[c] = -1;
	  // add the new set to fixlst:
	  fixlst[nfx] = c;
	  nfx++;
	  if(nfx>nmem) fixlst=realloc(fixlst,nfx*sizeof(INT));
	  h=c;
	  c++;
	  if(c>=nmem){
	    nmem=c+1;
	    flag=realloc(flag,nmem*sizeof(INT));
	    head=realloc(head,nmem*sizeof(INT));
	    next=realloc(next,nmem*sizeof(INT));
	    back=realloc(back,nmem*sizeof(INT));
	    head[c]=-1; back[c]=-1; next[c]=-1; flag[c]=-1;
	  }
	}
	// add cell of w to the new set
	next[cell[w]] = next[h];
	if (next[h] >= 0) {
	  back[next[h]] = cell[w];
	}
	flag[cell[w]] = h;
	back[cell[w]] = h;
	next[h] = cell[w];
      }//	end if
    }   //end for
    for(i0 = 0;i0< nfx;i0++){
      flag[fixlst[i0]] = -1;
    }// 
  }  //end do
  free(fixlst);
  free(flag);
  free(head);
  free(next);
  free(back);
  free(cell);
  //
  blk[0]->IA[blk[0]->row]=blk[0]->nnz;
  blk[1]->IA[blk[1]->row]=blk[1]->nnz;
  /* INT swp; */
  /* for(i=0;i<((INT ) n/2);i++){ */
  /*   swp=iord[i]; */
  /*   iord[i]=iord[n-i-1]; */
  /*   iord[n-i-1]=swp; */
  /* } */
  /**************************************set up level array*************************************/
  // backup the level array;
  /* for(i=0;i<n;++i) */
  /*   iord1[i]=level[i]; */
  /* // put the level array back setting the level numbers correctly */
  /* for(i=0;i<n;++i) */
  /*   level[i]=iord1[blk[1]->JA[i]]; */
  /* /\**************************************set up ancestor array*************************************\/ */
  /* // backup the ancestor array */
  /* for(i=0;i<n;++i) */
  /*   iord1[i]=anc->val[i]; */
  /* // put the level array back setting the level numbers correctly */
  /* for(i=0;i<n;++i) */
  /*   anc->val[i]=iord1[blk[1]->JA[i]]; */
  /* // */
  //finally, define the inverse permutation
  //////////////////////////////////////////////////
  /* for(i=0;i<n;i++) */
  /*   iord1[iord[i]]=i; */
  return blk;
}
/*****************************************************************************************/
/*!
 * \fn iCSRmat *run_bfs(INT n,INT *ia, INT *ja,ivector *roots,ivector *anc,const INT lvlmax)
 *
 * \brief recurrsive bfs search of a graph
 *
 * \param n                number of vertices
 * \param (ia,ja):         adjacency structure of the graph;
 * \param roots[]:         a vector containing starting points for BFS. Aiming
 *                         at having one starting point in every
 *                         connected component.
 * \param anc[]:           a vector containing the ancestors of the vertices:
 *                         v<->anc[v] is an edge in the bfs tree. 
 *
 * \return pointer to iCSRmat of size n,n with n nonzeroes, which
 *         contains the bfs ordering permutation of the vertices
 *         (component after component). The values in the matrix are
 *         the distance from the root vertex in every connected
 *         component: for the vertex v bfs->val[v]=distance to the
 *         root[j], j=1:ncomponents
 *
 * \author Ludmil Zikatanov
 * \date   20210516
 */
/*****************************************************************************/
iCSRmat *run_bfs(INT n,INT *ia, INT *ja,	\
		 ivector *roots,		\
		 ivector *anc,			\
		 const INT lvlmax)
{
  /* n, ia, ja are the CSR matrix elements. */
  /* anc[v] is the ancestor of v; */
  /* roots[] is input; */
  INT i,k,q,v,vi,qbeg,qend,lvl;  
  iCSRmat *blk=malloc(sizeof(iCSRmat));
  blk[0]=icsr_create(n,n,n);
  /* Not true that we need that: blk->IA=realloc(blk->IA,(n+2)*sizeof(INT));*/
  anc->row=n;
  anc->val=(INT *)calloc(anc->row,sizeof(INT));
  for(i=0;i<blk->nnz;++i) blk->val[i]=0;
  for(i=0;i<anc->row;++i) anc->val[i]=-1;
  lvl=0;
  blk->IA[lvl]=0;
  k=blk->IA[lvl];
  //  print_full_mat_int(1,roots->row,roots->val,"roots");
  if(roots->row<=0){
    /* take the first vertex as root if none are given as input */
    roots->row=1;
    roots->val=(INT *)realloc(roots->val,roots->row*sizeof(INT));
    roots->val[0]=0;
  }
  /* Now roots are set as they are either input or root[0]=0 */
  for(i=0;i<roots->row;++i){
    blk->val[roots->val[i]]=lvl+1;
    blk->JA[k]=roots->val[i];
    k++;
  }
  blk->IA[lvl+1]=k;
  if(n<=1)
    return blk;
  /* n>1 ...  */
  while(1){
    qbeg=blk->IA[lvl];
    qend=blk->IA[lvl+1];
    for(q=qbeg;q<qend;++q){
      v = blk->JA[q];
      for(vi=ia[v];vi<ia[v+1];++vi){
	i=ja[vi];
	if (!blk->val[i]){
	  blk->val[i]=lvl+1;// mark as visited;
	  blk->JA[k]=i;
	  k++;
	  anc->val[i]=v; // ancestor;
	}
	/* fprintf(stdout,"\nlvl=%d,v=%d; nbr=%d,blk->val=%d",lvl,v,i,blk->val[i]);fflush(stdout); */
      }
    }
    lvl++;
    if(k<=qend) break;
    /* fprintf(stdout,"\nn=%d,(lvl+1)=%d,k=%d",n,lvl+1,k); */
    blk->IA[lvl+1]=k;    
  }
  blk->row=lvl;
  blk->IA=(INT *)realloc(blk->IA,(blk->row+1)*sizeof(INT));
  /* icsr_print_rows(stdout,blk,"BLK"); */
  //end
  return blk;
}
/*****************************************************************************/
/*!
 * \fn iCSRmat *bfs_di(void *a, const char c,ivector *roots,
 *                     ivector *anc,const INT lvlmax)
 *
 * \brief bfs on graphs given by INT or REAL CSR matrix.
 *
 * \param a:                  The CSR matrix 
 * \param c:                  a char ('r' or 'i' or 'R' or 'I')
 *                            indicating whether this is a REAL or INT matrix;
 *
 * \return returns the output of bfs_di(a->row,a->IA,a->JA)
 *
 * \author Ludmil Zikatanov
 * \date   20210516
 */
iCSRmat *bfs_di(void *a, const char c,		\
		ivector *roots,			\
		ivector *anc,			\
		const INT lvlmax)
{
  /* 
     do a breadth first search for INT or REAL matrix. It does not use
     a->val, so upon entry these can be null 
  */
  dCSRmat *ad=NULL;
  iCSRmat *ai=NULL;
  if(c=='R' || c=='r'){
    ad=(dCSRmat *)a;
    return run_bfs(ad->row, ad->IA, ad->JA, roots, anc, lvlmax);  
  } else if(c=='I' || c=='i'){
    ai=(iCSRmat *)a;
    return run_bfs(ai->row, ai->IA, ai->JA, roots, anc, lvlmax);  
  } else {
    fprintf(stderr,"### ERROR: Wrong value of c in %s (c=%c should be \'r\' or \'R\')\n",__FUNCTION__,c);
    exit(ERROR_INPUT_PAR);
  }
}
/***********************************************************************************************/
/*!
 * \fn static void dfs00_(INT *nin, INT *ia, INT *ja, INT *nblko,INT *iblk, INT *jblk)
 *
 * \brief non-recurrsive dfs search of a graph.  This function
 *        implements Tarjan's DFS algorithm for finding the strongly
 *        connected components of a di-graph.  It is the FRED
 *        GUSTAVSON'S non-recurrsive version of the algorithm.
 *
 * \param n                number of vertices
 * \param (ia,ja):         adjacency structure of the graph; all diagonal
 *                         elements must be nonzero in ia,ja.
 * \param nblko:           number of connected components
 * \param (iblk,jblk):     CSR matrix containing the permutation: vertices
 *                         are ordered ccomponent by compoent.
 *
 * \return The new numbering is in CSR format: the pointers are in
 *         IBLK[], and the actual numbering is in JBLK[].
 *
 * \author Ludmil Zikatanov
 * \date   20210516
 */
static void dfs00_(INT *nin, INT *ia, INT *ja, INT *nblko,INT *iblk, INT *jblk)
{
  INT n,n1,nb,nblk,count,k1,sp,vp,v=0,wp=0,w=0,v1=0,i,ii,myflag;
  INT *numb=NULL,*iedge=NULL,*lowlink=NULL;
  /*---------------------------------------------------------------------
    ---------------------------------------------------------------------*/
  /*C...  Initialization.*/
  n=*nin;
  n1 = n+1;
  iedge = (INT *) calloc(n1,sizeof(INT));
  numb = (INT *) calloc(n1,sizeof(INT));
  lowlink= (INT *) calloc(n1,sizeof(INT));
  /* Check */
  if(!(numb && iedge && lowlink)) {
    fprintf(stderr,"\nCould not allocate local mem in dfs\n");
    exit(19);
  }

  nblk = -1;
  nb = -1;
  count  = nb+1;
  k1 = count;
  vp = n;//
  sp = n;// may be n???
  for (i=0; i<n;i++) {
    iedge[i] = ia[i];
    numb[i]=-1;
    lowlink[i]=-1;
  }
  numb[n] = -1;
  lowlink[n] = -1;
  myflag=10;
  while (1) { //10 continue
    //    fprintf(stdout," myflag=%i %i, %i %i\n",myflag,nblk,count,n);fflush(stdout);
    if(myflag==10){
  /*
    C...  Get out when the  renumbering is done;
    C...  Otherwise begin at vertex K1
  */
      if(count==n) {
	nblk++; //not sure why....
	iblk[nblk] = n;
	*nblko=nblk;
	if(iedge) free(iedge);
	if(numb) free(numb);
	if(lowlink) free(lowlink);
	return;
      }
      for (ii=k1;ii<n;ii++) {
	i=ii;
	if(numb[ii]<0) {
	  myflag=30;
	  break; //go to 30
	}
      }
      if(myflag !=30){
	fprintf (stderr,"There is an error in DEPTH FIRST SEARCH:  %i == %i\n",i,n);
	exit(254);
      } else {
	v = i;
	k1 = v + 1;
	myflag=50;
      }
    }
    /*C...  :::*/
    if(myflag==40) {
      vp--;
      iblk[vp] = v;
      v=w;
      myflag=50;
    }
    if(myflag==50) {
      nb++;
      numb[v] = nb;
      lowlink[v] = numb[v];
      sp--;
      jblk[sp] = v;
      v1 = v+1;
      myflag=60;
    }
    /*...  */
    while(60) { // 60   continue;
      if(myflag == 60) {
	wp = iedge[v];
	//	fprintf(stdout,"\nv=%d myflag=%i wp=%i, count=%i n=%i\n",v,myflag,wp,count,n);fflush(stdout);
	w = ja[wp];
	iedge[v] = wp+1;
	if(numb[w] >= numb[v])  {
	  myflag=70;
	} else if (numb[w]<0) {
	  myflag=40;
	  break;
	} else {
	  if(lowlink[v]>numb[w])
	    lowlink[v]=numb[w];
	  myflag=70;
	}
      }
      if(iedge[v] < ia[v1]) {
	myflag=60;
	continue;
      }
      /*...*/
      if(lowlink[v] >= numb[v]) {//don't {go to 90}
	nblk++;
	iblk[nblk] = count;
	/**/
	while(80) {     //    80 continue;
	  w = jblk[sp];
	  numb[w] = n;
	  sp++;
	  jblk[count] = w;
	  count++;
	  if(v == w) break; //{don't go to 80}
	}
	/*C...  */
	if(sp==n) {
	  myflag=10;
	  break;
	}
      }
      myflag=70;
      w = v;
      v = iblk[vp];
      vp++;
      v1 = v + 1;
      if(lowlink[v]>lowlink[w])
	lowlink[v]=lowlink[w];
    }
  }
}
/***********************************************************************************************/
/*!
 * \fn iCSRmat *run_dfs(INT n, INT *ia, INT *ja)
 *
 * \brief recurrsive dfs search of a graph
 *
 * \param n                number of vertices
 * \param (ia,ja):         adjacency structure of the graph;
 *
 * \return pointer to iCSRmat of size n,n with n nonzeroes, which
 *         contains the connected components and the permutation of
 *         the vertices (component after component). The values in the
 *         matrix are connected component number for the vertex
 *         dfs->val[i]=connected component number(dfs->JA[i])
 *
 * \author Ludmil Zikatanov
 * \date   20210516
 */
iCSRmat *run_dfs(INT n, INT *ia, INT *ja)
{
  INT i,pos,flag;
  //  INT j,k;
  // short hand;
  iCSRmat *dfs=malloc(sizeof(iCSRmat));
  dfs[0]=icsr_create(n,n,n);
  INT *iawrk=calloc(n+1,sizeof(INT));
  INT *jawrk=calloc(ia[n]+n,sizeof(INT)); // use it to add a diagonal.
  INT nnzwrk=0;
  iawrk[0]=nnzwrk;
  for(i=0;i<n;++i){
    flag=0;
    for(pos=ia[i];pos<ia[i+1];++pos){
      jawrk[nnzwrk]=ja[pos];
      if(ja[pos]==i) flag=1;// diagonal was found. 
      nnzwrk++;
    }
    if(!flag){
      jawrk[nnzwrk]=i;
      nnzwrk++;
    }
    iawrk[i+1]=nnzwrk;
  }
  //  fprintf(stdout,"\n*******:::: %d and %d\n",ia[n],iawrk[n]);
  jawrk=realloc(jawrk,iawrk[n]*sizeof(INT));
  dfs->row=-10;
  dfs00_(&n,iawrk,jawrk,&dfs->row,dfs->IA,dfs->JA);
  dfs->IA=realloc(dfs->IA,(dfs->row+1)*sizeof(INT));
  free(iawrk);
  free(jawrk);
  INT ipstrt=-10,ipend=-10;// INT lp=-10,swp;
  for(i=0;i<dfs->row;++i){
    ipstrt=dfs->IA[i];
    ipend=dfs->IA[i+1];
    //    lp=(INT )((ipend-ipstrt)/2);
    for(pos=ipstrt;pos<ipend;++pos){
      dfs->val[pos]=i+1;
    }
    /* for(pos=0;pos<lp;++pos){ */
    /*   swp=dfs->JA[ipstrt+pos]; */
    /*   dfs->JA[ipstrt+pos]=dfs->JA[ipend-pos-1]; */
    /*   dfs->JA[ipend-pos-1]=swp; */
    /* } */
  }
  //  icsr_print_matlab(stdout,dfs);fflush(stdout);
  /* for(i=0;i<dfs->row;++i){ */
  /*   for(pos=dfs->IA[i];pos<dfs->IA[i+1];++pos){ */
  /*     dfs->val[pos]=i+1; */
  /*   } */
  /* } */
  //  icsr_print_matlab(stdout,dfs);fflush(stdout);
  return dfs;
}
/*****************************************************************************/
/*!
 * \fn iCSRmat *dfs_di(void *a, const char c)
 *
 * \brief dfs on graphs given by INT or REAL CSR matrix.
 *
 * \param a:                  The CSR matrix 
 * \param c:                  a char ('r' or 'i' or 'R' or 'I')
 *                            indicating whether this is a REAL or INT matrix;
 *
 * \return returns the output of run_dfs(a->row,a->IA,a->JA)
 *
 * \author Ludmil Zikatanov
 * \date   20210516
 */
iCSRmat *dfs_di(void *a, const char c)
{
  /* 
     do a depth first search for INT or REAL matrix. It does not use
     a->val, so upon entry this can be null 
  */
  dCSRmat *ad=NULL;
  iCSRmat *ai=NULL;
  if(c=='R' || c=='r'){
    ad=(dCSRmat *)a;
    return run_dfs(ad->row, ad->IA, ad->JA);  
  } else if(c=='I' || c=='i'){
    ai=(iCSRmat *)a;
    return run_dfs(ai->row, ai->IA, ai->JA);  
  } else {
    fprintf(stderr,"### ERROR: Wrong value of c in %s: c=%c\n",__FUNCTION__,c);
    exit(ERROR_INPUT_PAR);
  }
}
/*********************************************************************************************************/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx*/
INT check0(weights *elem1, weights *elem2)
{
  if ( elem1->val < elem2->val)
    return -1;
   else if (elem1->val > elem2->val)
      return 1;
   else
      return 0;
}

INT check1(iweights *elem1, iweights *elem2)
{
  //descending order. 
  if  (elem1->mask > elem2->mask)
    return -1;
  else if(( elem1->mask == elem2->mask) && (elem1->val > elem2->val))
    return -1;
  else if(elem1->mask < elem2->mask)
    return 1;
  else if(( elem1->mask == elem2->mask) && (elem1->val < elem2->val))
    return 1;
  else
    return 0;
}

void getp(INT *ie, INT *je, REAL *w, INT ne, INT *p)
{
  INT k;
  weights *tosort=NULL;  
  tosort = (weights *)calloc(ne,sizeof(weights));
  for (k=0; k<ne;k++)
    {
      tosort[k].val=w[k];
      tosort[k].id=k;
    }
  //Sort now.
  qsort((void *) tosort,ne,sizeof(weights),(testit )check0 );                  
  for (k=0; k<ne;k++)
    p[k]=tosort[k].id;

  if(tosort) free(tosort);  
  return;
}
void getpz(REAL *z, INT nv, INT *p)
{
  INT k;
  weights *tosort=NULL;  
  tosort = (weights *)calloc(nv,sizeof(weights));
  for (k=0; k<nv;k++)
    {
      tosort[k].val=-z[k];
      tosort[k].id=k;
    }
  qsort((void *) tosort,nv,sizeof(weights),(testit )check0 );                  
  // after things are sorted,  we get the permutation
  for (k=0; k<nv;k++)
    p[k]=tosort[k].id;
  if(tosort) free(tosort);  
  return;
}
void getpi(INT *iz, INT *maskv, INT nv, INT *p)
{
  INT k;
  iweights *tosort=NULL;  
  tosort = (iweights *)calloc(nv,sizeof(iweights));
  for (k=0; k<nv;k++)
    {
      tosort[k].mask=maskv[k];
      tosort[k].val=iz[k];
      tosort[k].id=k;
    }
  qsort((void *) tosort,nv,sizeof(iweights),(testit )check1 );                  
  // after things are sorted,  we have the permutation
  for (k=0; k<nv;k++)
    p[k]=tosort[k].id;
  if(tosort) free(tosort);  
  return;
}
