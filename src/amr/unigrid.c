/*! \file src/amr/unigrid.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *   \note routines to initialize construct and destruct uniform grids
 *  in d-dimensions for general d. 
*/
#include "hazmath.h"
/**********************************************************************/
/*!
 * \fn void coord_lattice(INT *m,const INT dim,const INT kf, 
 *                        const INT nall, const INT *nd)
 *
 * \brief given a global number kf on a lattice grid with
 *        lexicographical ordering, this returns the n-tuple of latice
 *        coordinates, m[i],i=0:n-1). nall is the total number of
 *        vertices in the lattice, nd[i] is the number of divisions in
 *        every direction;
 *
 * \param 
 *
 * \return
 *
 * \note
 *
 */
void coord_lattice(INT *m,const INT dim,const INT kf, \
		   const INT nall, const INT *nd)
{
  INT i,j,k;
  /*
  */
  j=nall;  k=kf;
  for(i=dim;i>0;i--){
    j=j/(nd[i-1]+1);  m[i-1] = k/j;   k=k-m[i-1]*j;
  }
  return;
}
/**********************************************************************/
/*!
 * \fn INT num_lattice(INT *m,const INT dim,INT *nd)
 *
 * \brief given the lattice coordinates m[i],i=0:n-1 of a vertex, in
 *        dimension n, finds the global number kf of a vertex on a
 *        lexicographically ordered lattice grid. Divisions in each
 *        directions are stored in nd[].
 *
 * \param 
 *
 * \return
 *
 * \note
 *
 */
INT num_lattice(INT *m,const INT dim,INT *nd)
{
  INT i,kf;
  /*
  */
  kf=m[dim-1];
  //      kf = m[dim-1];
  for (i=dim-1; i>0; i--){
    kf=kf*(nd[i-1]+1)+m[i-1];
  }
  return kf;
}

void binary1(const INT dim, unsigned int *bits, INT *nvloc)
{
  // coordinates of the vertices of the unit dim-cube as arrays of 0/1 
  INT i,k,kdim=-10,nbits=dim-1;
  *nvloc = (1 << dim);
  for(k = 0;k<(*nvloc);k++){
    kdim=k*dim;
    for (i=nbits ; i >=0; --i){
      bits[kdim+i] = (unsigned int )(!(k >> i & 1));
    }
  }
  return;
}

/**********************************************************************/
/*!
 * \fn scomplex *umesh(const INT dim,INT *nd, cube2simp *c2s,\ 
 *                    INT *labelf,INT *isbndf, INT *codef,INT elflag, \
 *                    const INT intype)
 *
 * \brief Uniform simplicial mesh of the unit cube in dimension dim.nd
 *        is the number of grid points in each dimension.  ordering is
 *        lexicographically by name=(x[0],...,x[n]).  more than 3D is
 *        not fully tested xmacro[] are the coordinates of a domain
 *        isomorphic to the cube via a bilinear or "Q2" change of
 *        coordinates. output is a simplicial complex sc.
 *
 * \param intype: type of uniform grid. 
 *
 *    if(intype == -2) use unirefine() function (from unigrid.c in
 *   src/amr) if(intype == -1)construct grid using diagonals pointing
 *   0-7(0...0)-->(1...1).  if (intype>0) starting with intype the
 *   mesh is constructed like criss-cross grid. This works in 2D and
 *   3D, and is unclear whether it works in d>3.
 *
 * \return scomplex containing the simplicial grid. 
 *
 * \note
 *
 */
scomplex *umesh(const INT dim,					\
		INT *nd, cube2simp *c2s,			\
		INT *labelf,INT *isbndf, INT *codef,INT elflag, \
		const INT intype)
{
  /*   */
  INT iz1;
  INT jperm,i,j,flag,kf,type;
  INT dim1 = dim+1;
  // m is dim+1 so that we can handle even dimensions
  INT *m = (INT *)calloc(dim1,sizeof(INT));
  INT *mm = (INT *)calloc(dim1,sizeof(INT));
  INT *cnodes = (INT *)calloc(c2s->nvcube,sizeof(INT));  
  //  INT *icycle = (INT *)calloc(dim+1,sizeof(INT));
  INT nv=1,ns=1;
  for(i=0;i<dim;i++){
    nv*=(nd[i]+1);
    ns*=nd[i];
  }
  ns*=c2s->ns; /*multiply by the number of simplices in the unit cube
		 (2 in 2D and 6 in 3d and 24 in 4d*/
  scomplex *sc = (scomplex *)haz_scomplex_init(dim,ns,nv,dim); 
  //  fprintf(stdout,"\nFaces=(%d,%d)=(face,face_parent)\n",face,face_parent);fflush(stdout);
  for(kf=0;kf<sc->nv;kf++){
    coord_lattice(m,dim,kf,sc->nv,nd);
    for(i=0;i<dim;i++){
      /* THIS HERE HAS THE X COORD FIRST WHICH IS NOT WHAT ONE HAS IF
	 USING THE BIJECTION BETWEEN BINARY NUMBERS AND THE
	 COORDINATES OF VERTICES IN THE UNIT CUBE> SO WE REVERSE THE
	 ORDERING OF DIVISIONS SO THAT WE PARTITION FIRST X and so
	 on. so basically we come here with the last coordinate
	 first. That is why we also have (dim-i-1) instead of i*/
      sc->x[kf*dim+(dim-i-1)]=((REAL )m[i])/((REAL )nd[i]);
      /* OLD: sc->x[kf*dim+i]=((REAL )m[i])/((REAL )nd[i]); */
    }    
    //      print_full_mat_int(1,dim,m,"m1=");
    //      fprintf(stdout,"; iglobal=%d;",kf);
  }
  ns=0;
  for(kf=0;kf<sc->nv;kf++){
    coord_lattice(m,dim,kf,sc->nv,nd);
    flag=0;
    for(i=0;i<dim;i++){
      if(m[i]==nd[i]) {flag=1; break;}
    }
    if(flag) continue;
    if(intype==-1) {
      type=0;
    }else{
      /*criss-cross in any D*/
      /* determine type; this is not fully rigorously justified, but
       works in d=2,3*/
      type=(m[0]+intype)%2;
      for(i=1;i<dim-1;i++){
	type+=2*(m[i]%2);
      }
      // this is a hack here to work in 2D. unclear how to do in 2D yet or 4D. 
      if(dim==2){type=(abs(m[1]-m[0])+intype)%dim;}
      if((m[dim-1]%2)) type=dim-type;
      if((!(dim%2)) && (type>=(dim))) {type%=(dim);}
      if(dim%2 && (type>(dim+1))) {type%=(dim+1);}
    }
    //    type=0;
    /*depending on the type, split a cube in simplices*/
    //    fprintf(stdout,"\ntype=%d\n",type+1);
    for(j=0;j<c2s->nvcube;j++){
      //            fprintf(stdout,"j:%d; ",j+1);
      for(i=0;i<dim;i++){
    	mm[i]=m[i]+(c2s->bits[dim*j+i]);
	//		fprintf(stdout,"mm[%d]=%d; ",i+1,mm[i]+1);
      }
      cnodes[j]=num_lattice(mm,dim,nd);
      //            fprintf(stdout,"\n"); fflush(stdout);
    }   
    //    fprintf(stdout,"\n"); fflush(stdout);
    //  
    for(i=0;i<c2s->ns;i++){
      for(j=0;j<dim1;j++){
	iz1=c2s->nodes[i*dim1+j];
	jperm=c2s->perms[type*c2s->nvcube+iz1];
	//	fprintf(stdout,"\ntype=%d,ns=%d,jperm=%d,iz1=%d",type,ns,jperm,iz1);
	sc->nodes[ns*dim1+j]=cnodes[jperm];
      }
      sc->flags[ns]=elflag;
      ns++;      
    }    
  }
  INT cfbig=((INT )MARKER_BOUNDARY_NO)+100;
  INT facei,bf,cf,mi;
  //  INT kfp,ijk,mi,mip,toskip,toadd;
  /******************************************************************/
  /*  
   *  when we come here, all boundary faces have codes and they are
   *  non-zero. All interior faces should have a code zero.
   */
  /******************************************************************/
  for(kf=0;kf<sc->nv;kf++) sc->bndry[kf]=cfbig;
  /*  icsr_print_matlab(stdout,sc->parent_v);*/
  /* for(facei=0;facei<c2s->nf;facei++){ */
  /*   if(facei<dim){ */
  /*     mi=dim-(facei+1); */
  /*     bf=0; */
  /*   } else{ */
  /*     mi=dim-((facei%dim)+1); */
  /*     bf=nd[mi]; */
  /*   } */
  /*   cf=codef[facei]; */
  /*   if(!isbndf[facei]){ */
  /*     // first pass: set the interior faces; */
  /*     for(kf=0;kf<sc->nv;kf++){ */
  /* 	coord_lattice(m,dim,kf,sc->nv,nd); */
  /* 	if(m[mi]==bf){ */
  /* 	  if(sc->bndry[kf]>cf && (cf !=0)) sc->bndry[kf]=cf; */
  /* 	} */
  /*     } */
  /*   } */
  /* } */
  /******************************************************************/
  // second pass: set boundaries, so that the bondaries are the ones
  // that we care about:
  /******************************************************************/
  icsr_realloc(sc->nv,sc->n,sc->n*sc->nv,sc->bndry_v); // a vertex belongs to at most dim and also for the codes.
  sc->bndry_v->val=realloc(sc->bndry_v->val,2*sc->bndry_v->nnz*sizeof(INT));// 2 values per vertex per dimension
  // init the column indices to negative
  for(j=0;j<sc->bndry_v->nnz;++j)
    sc->bndry_v->JA[j]=-1;
  memset(sc->bndry_v->val,0,2*sc->bndry_v->nnz*sizeof(INT));
  sc->bndry_v->IA[0]=0;
  for(kf=0;kf<sc->nv;++kf){
    /* from 1 to n this holds the address of the beginning of the previous row */
    sc->bndry_v->IA[kf+1]=sc->bndry_v->IA[kf]+sc->n; 
  }
  /**/
  for(facei=0;facei<c2s->nf;facei++){
    if(facei<dim){
      mi=dim-(facei+1);
      bf=0;
    } else{
      mi=dim-((facei%dim)+1);
      bf=nd[mi];
    }
    cf=codef[facei];
    if(isbndf[facei]){
      for(kf=0;kf<sc->nv;kf++){
	coord_lattice(m,dim,kf,sc->nv,nd);
	if(m[mi]==bf){
	  if(sc->bndry[kf]>cf) sc->bndry[kf]=cf;
	}
      } 
    }
    for(kf=0;kf<sc->nv;kf++){
      coord_lattice(m,dim,kf,sc->nv,nd);
      if(m[mi]==bf){
	j=sc->bndry_v->IA[kf];
	sc->bndry_v->JA[j]=facei;
	sc->bndry_v->IA[kf]++;
	//	fprintf(stdout,"\nvertex=%d; face=%d; j=%d; bf=%d; cf=%d; bndry=%d,label=%d",kf,facei,j,bf,cf,isbndf[facei],labelf[facei]);
      } 
    }
  }
  fflush(stdout);
  //return IA to previous state.
  sc->bndry_v->IA[0]=0;  
  for(kf=0;kf<sc->nv;kf++)
    sc->bndry_v->IA[kf+1]=sc->bndry_v->IA[kf]+sc->n;
  INT row_begin=sc->bndry_v->IA[0],nnz=0;  
  for(kf=0;kf<sc->nv;kf++){
    for(j=row_begin;j<sc->bndry_v->IA[kf+1];++j){
      if(sc->bndry_v->JA[j]<0) continue;
      sc->bndry_v->JA[nnz]=sc->bndry_v->JA[j];
      nnz++;
    }
    row_begin=sc->bndry_v->IA[kf+1];
    sc->bndry_v->IA[kf+1]=nnz;
  }
  // set up the final sc->bndry_v for this uniform mesh
  sc->bndry_v->JA=realloc(sc->bndry_v->JA,nnz*sizeof(INT));
  sc->bndry_v->val=realloc(sc->bndry_v->val,2*nnz*sizeof(INT));
  sc->bndry_v->nnz=nnz;
  for(kf=0;kf<sc->nv;kf++){
    for(j=sc->bndry_v->IA[kf];j<sc->bndry_v->IA[kf+1];++j){
      facei=sc->bndry_v->JA[j];
      sc->bndry_v->JA[j]=labelf[facei];
      //      sc->bndry_v->JA[j]=facei;
      sc->bndry_v->val[j]=codef[facei];
      sc->bndry_v->val[nnz+j]=isbndf[facei];
    }
  }
  /******************************************************************/
  // Only interior points should now be left; set them to 0:
  for(kf=0;kf<sc->nv;kf++)
    if(sc->bndry[kf]>=cfbig) sc->bndry[kf]=0;      
  /*clean up*/
  if(m) free(m);
  if(mm) free(mm);
  if(cnodes) free(cnodes);
  /**/
  return sc;
}
/**********************************************************************/
/*!
 * \fn void unirefine(INT *nd,scomplex *sc)  
 *
 * \brief refine uniformly l levels, where 2^l>max_m nd[m] using the
 *        generic algorithm for refinement.  Works in the following
 *        way: first construct a grid with refinements up to 2^{l}
 *        such that 2^{l}>max_m nd[m]. then remove all x such that
 *        x[k]>nd[k]*2^{-l} and then remap to the unit square.
 *
 * \param 
 *
 * \return
 *
 * \note  not used; (20180718)--ltz
 *
 */
void unirefine(INT *nd,scomplex *sc)  
{
/* 
*/
  INT ndmax=-1,i=-1,j=-1;
  for(i=0;i<sc->n;i++)
    if(ndmax<nd[i]) ndmax=nd[i];
  //  fprintf(stdout,"\nmax split=%d",ndmax);
  REAL sref=log2((REAL )ndmax);
  if(sref-floor(sref)<1e-3)
    sref=floor(sref);
  else
    sref=floor(sref)+1.;
  INT ref_levels= sc->n*((INT )sref);
  //  fprintf(stdout,"\nlog2 of the max=%e, l=%d",log2((REAL )ndmax)+1,ref_levels);
  find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
  //  haz_scomplex_print(sc,0,__FUNCTION__);  fflush(stdout);
  INT *wrk=calloc(5*(sc->n+2),sizeof(INT));
  /* construct bfs tree for the dual graph */
  abfstree(0,sc,wrk,0);
  free(wrk);
  ref_levels=0;
  if(ref_levels<=0) return;
  INT nsold,print_level=0;//ns,nvold,level;
  if(!sc->level){
    /* form neighboring list; */
    find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
    //    haz_scomplex_print(sc,0,__FUNCTION__);  fflush(stdout);
    /* construct bfs tree for the dual graph */
    abfstree(0,sc,wrk,print_level=0);
    //    haz_scomplex_print(sc,0,__FUNCTION__);fflush(stdout);
    //    exit(100);
  }
  //INT n=sc->n, n1=n+1,level=0;
  fprintf(stdout,"refine: ");
  while(sc->level < ref_levels && TRUE){
    nsold=sc->ns;
    //    nvold=sc->nv;
    for(j = 0;j < nsold;j++)sc->marked[j]=TRUE;
    for(j = 0;j < nsold;j++)
      if(sc->marked[j] && (sc->child0[j]<0||sc->childn[j]<0))
	haz_refine_simplex(sc, j, -1);
    /* new mesh */
    //    ns=sc->ns; 
    sc->level++;
    fprintf(stdout,"u%du",sc->level);//,nsold,ns,nv);
  }
  fprintf(stdout,"\n");
  scfinalize(sc,(INT )0);
  return;
}
/**********************************************************************/
/*!
 * \fn unigrid *ugrid_init(INT n, INT *nd, REAL *xo, REAL *xn)
 *
 * \brief init structure unigrid
 *
 * \param 
 *
 * \return
 *
 * \note
 *
 */
unigrid *ugrid_init(INT n, INT *nd, REAL *xo, REAL *xn)
{
  INT j;
  unigrid *ug=(unigrid *)malloc(sizeof(unigrid));
  ug->n=n; 
  ug->nvcube=(1<<n); 
  ug->ndiv=nd; 
  ug->bits=(unsigned int *)calloc(n*(ug->nvcube),sizeof(unsigned int));
  binary1(n,ug->bits,&(ug->nvcube));
 if(!xo || !xn){
   ug->xo=(REAL *)calloc(n,sizeof(REAL)); /* coordinates of the origin
					     xo[dim]*/
   ug->xn=(REAL *)calloc(n,sizeof(REAL)); /* coordinates of the max
					     corner(NE in 2D)
					     xn[dim]*/   
   for(j=0;j<n;j++){
     ug->xo[j]=0.;
     ug->xn[j]=1.;
   }
 } else {
   ug->xo=xo;
   ug->xn=xn;
 }
 ug->dx = (REAL *)calloc(n,sizeof(REAL));
 ug->nall=1;
 for(j=0;j<n;j++){
   ug->dx[j]=(ug->xn[j]-ug->xo[j])/((REAL )ug->ndiv[j]);
   ug->nall *= (ug->ndiv[j]+1);
 }
 /* allocate space for the data */
 ug->data=(REAL *)calloc(ug->nall,sizeof(REAL));
 if(!ug->data) {
   fprintf(stdout,"\n***ERROR: could not allocate memory for a uniform grid with data in %s\n",__FUNCTION__);fflush(stdout);
   exit(7);
 }
 ug->ugtype = 0; 
 return ug; 
}
void ugrid_free(unigrid *ug)
{
  if(ug->bits) free(ug->bits);
  if(ug->dx) free(ug->dx);
  if(ug->data) free(ug->data);
  ug->data=NULL;
  if(ug) free(ug);
  ug=NULL;
  return; 
}
/**********************************************************************/
/*!
 * \fn void ugrid_transform(const INT n,unigrid *ug,REAL *xodst, REAL *xndst)
 *
 * \brief transform the unform grid ug to be over a parallelepiped
 *        with corners xodst and xndst
 *
 * \param 
 *
 * \return
 *
 * \note
 *
 */
void ugrid_transform(const INT n,unigrid *ug,REAL *xodst, REAL *xndst)
{
  /* transform the unform grid ug to be over a parallelepiped with
     corners xodst and xndst */
  INT j;
  //  REAL a,b,dxsrc,dxdst;
  ug->xo=xodst;
  ug->xn=xndst;
  for(j=0;j<n;j++){
    /*    
	  dxsrc=ug->xn[j]-ug->xo[j];
	  dxdst=xndst[j]-xodst[j];
	  a=dxdst/dxsrc; b=(xodst[j]*ug->xn[j]-xndst[j]*ug->xo[j])/dxsrc;
    */
    ug->dx[j]=(ug->xn[j]-ug->xo[j])/((REAL )ug->ndiv[j]);
  }
  return;
}
/*EOF*/
