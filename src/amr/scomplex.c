/*! \file src/amr/scomplex.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note containing all essentials routines for mesh refinement
 *
 */
#include "hazmath.h"
/**********************************************************************/
/*!
 * \fn REAL chk_sign(const int it, const int nbrit)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
REAL chk_sign(const int it, const int nbrit)
{
/*
  nbrit is a neighbor of it and this picks the sign of the normal
  vector for the face shared by both. Used to construct dcsrmat atf
  and local matrices to build aff. If the face is on the boundary,
  i.e. nbrit <0, we always return 1 (i.e. the outward normal).
*/
  if(it>nbrit) return 1e0;
  return -1e0;
}
/**********************************************************************/
/*!
 * \fn scomplex *haz_scomplex_init(const INT n,INT ns, INT nv, const INT nbig)
 *
 * \brief Initialize simplicial complex in dimension n with ns
 *        simplices and nv vertices.
 *
 * \param n is the dimension of the simplicial complex; 
 * \param ns is the number of simplices
 * \param nv is the number of vertices
 * \param nbig is the dimension of the space where the simplicial
 *        complex is embedded;
 *
 * \return initialized structure of type scomplex
 *
 * \note
 *
 */
scomplex *haz_scomplex_init(const INT n,INT ns, INT nv,const INT nbig)
{
  /*
     n = dimension of the simplicial complex;
     nbig = dimension of the space where the complex is embedded;
     Future work: we should think of making this for different
     dimensions, e.g. 2-homogenous complex in 3d
  */
  scomplex *sc=(scomplex *) malloc(sizeof(scomplex));
  sc->nbig=nbig; sc->n=n;
  if(sc->nbig<=0)sc->nbig=n;
  INT n1=sc->n+1,i,j,in1;
  sc->factorial=1.;
  sc->level=0;
  for (j=2;j<n1;j++) sc->factorial *= ((REAL )j);
  // fprintf(stdout,"\nIMPORTANT: NS=%d (%d!)=%f",ns,n,sc->factorial);
  sc->marked=(INT *) calloc(ns,sizeof(INT));
  sc->gen=(INT *) calloc(ns,sizeof(INT));
  sc->nbr=(INT *) calloc(ns*n1,sizeof(INT));
  sc->parent=(INT *)calloc(ns,sizeof(INT));
  sc->child0=(INT *)calloc(ns,sizeof(INT));
  sc->childn=(INT *)calloc(ns,sizeof(INT));
  sc->nodes=(INT *)calloc(ns*n1,sizeof(INT));
  sc->bndry=(INT *)calloc(nv,sizeof(INT));
  sc->csys=(INT *)calloc(nv,sizeof(INT));/* coord sys: 1 is polar, 2
					    is cyl and so on */
  sc->flags=(INT *)calloc(ns,sizeof(INT));
  sc->x=(REAL *)calloc(nv*nbig,sizeof(REAL));
  sc->vols=(REAL *)calloc(ns,sizeof(REAL));
  sc->fval=(REAL *)calloc(nv,sizeof(REAL)); // function values at every vertex; not used in general;
  for (i = 0;i<ns;i++) {
    sc->marked[i] = FALSE; // because first is used for something else.
    sc->gen[i] = 0;
    sc->parent[i]=-1;
    sc->child0[i]=-1;
    sc->childn[i]=-1;
    sc->flags[i]=-1;
    in1=i*n1;
    for(j=0;j<n1;j++){
      sc->nodes[in1+j]=-1;
      sc->nbr[in1+j]=-1;
    }
  }
  for (i = 0;i<nv;i++) {
    sc->bndry[i]=0;
    sc->csys[i]=0;
    sc->fval[i]=0.;
  }
  sc->nv=nv;
  sc->ns=ns;
  sc->bndry_cc=1; // one connected component on the boundary for now.
  sc->cc=1; // one connected component for now.
  return sc;
}
/**********************************************************************/
/*!
 * \fn scomplex *haz_scomplex_init_part(scomplex *sc)
 *
 * \brief Initialize simplicial complex which has already been
 *        allocated partially. Assumption is: nv,ns,nodes,nbr, cc,
 *        bndry_cc are known. in dimension n with ns simplices and nv
 *        vertices.
 *
 * \param sc pointer to partially allocated simplician complex; 
 *
 * \return full structure of type scomplex with all allocations done
 *
 * \note
 *
 */
void haz_scomplex_init_part(scomplex *sc)
{
  INT nv=sc->nv,ns=sc->ns,n1=sc->n+1,i,j;
  sc->factorial=1.;
  sc->level=0;
  for (j=2;j<n1;j++) sc->factorial *= ((REAL )j);
  // fprintf(stdout,"\nIMPORTANT: NS=%d (%d!)=%f",ns,n,sc->factorial);
  sc->marked=(INT *) calloc(ns,sizeof(INT));
  sc->gen=(INT *) calloc(ns,sizeof(INT));
  sc->parent=(INT *)calloc(ns,sizeof(INT));
  sc->child0=(INT *)calloc(ns,sizeof(INT));
  sc->childn=(INT *)calloc(ns,sizeof(INT));
  sc->bndry=(INT *)calloc(nv,sizeof(INT));
  sc->csys=(INT *)calloc(nv,sizeof(INT));/* coord sys: 1 is polar, 2
					    is cyl and so on */
  sc->flags=(INT *)calloc(ns,sizeof(INT)); // element flags
  sc->vols=(REAL *)calloc(ns,sizeof(REAL));// simplex volumes
  sc->fval=(REAL *)calloc(nv,sizeof(REAL));// simplex volumes
  for (i = 0;i<ns;i++) {
    sc->marked[i] = FALSE; // because first is used for something else.
    sc->gen[i] = 0;
    sc->parent[i]=-1;
    sc->child0[i]=-1;
    sc->childn[i]=-1;
    sc->flags[i]=-1;
  }
  for (i = 0;i<nv;i++) {
    sc->bndry[i]=0;
    sc->csys[i]=0;
    sc->fval[i]=0.;
  }
  return;
}
/**********************************************************************/
/*!
 * \fn void vol_simplex(INT dim, REAL fact, REAL *xf, REAL *volt, void *wrk)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void vol_simplex(INT dim, REAL fact, REAL *xf, REAL *volt, void *wrk)
{
  /*
     computes the volume of a simplex; wrk should be at least
     dim*(dim+1) REALS and dim integers.
  */
  INT dim1 = dim+1,i,j,ln,ln1;
  REAL *bt=(REAL *)wrk;
  REAL *piv=bt+dim*dim;
  INT *p = (INT *)(wrk+(dim*dim + dim)*sizeof(REAL));
  // construct bt using xf;
  for (j = 1;j<dim1;j++){
    ln=j*dim; ln1=ln-dim;
    for(i=0;i<dim;i++){
      bt[ln1+i] = xf[ln+i]-xf[i];
    }
  }
  //  print_full_mat(dim,dim,bt,"bt");
  if(lufull(1, dim, volt, bt,p,piv))
    *volt=0e0;
  else
    *volt=fabs(*volt)/fact;
  return;
}
/**********************************************************************/
/*!
 * \fn scomplex *haz_scomplex_read(FILE *fp)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
scomplex *haz_scomplex_read(FILE *fp,INT print_level)
{
  INT i,ns,nv,n,dummy;
  i=fscanf(fp,"%d %d %d %d",&ns,&nv,&n,&dummy);
  INT n1=n+1,j,k,n1kj=-10,nbig=n;// we can only read same dimension complexes now. 
  scomplex *sc=(scomplex *)haz_scomplex_init(n,ns,nv,n);
  for (j=0;j<n1;j++) {
    for (k=0;k<ns;k++){
      n1kj=n1*k+j;
      dummy=fscanf(fp," %d ", sc->nodes+n1kj);
      /* shift if needed ; this should not be here: later CHANGE */
      sc->nodes[n1kj]=sc->nodes[n1kj]-1;
    }
  }
  for (k=0;k<ns;k++){
    dummy=fscanf(fp," %d ", sc->flags+k);
  }
  for(j=0;j<nbig;j++){
    for(i=0;i<nv;i++){
      dummy=fscanf(fp,"%lg",sc->x+i*nbig+j);
    }
  }
  for(i=0;i<nv;i++){
    dummy=fscanf(fp,"%i",sc->bndry+i);
    /* fprintf(stdout,"%i: %i\n",i,sc->bndry[i]); */
  }
  for(i=0;i<nv;i++){
    dummy=fscanf(fp,"%lg",sc->fval+i);
    if(dummy<0 && (print_level>5)){
      fprintf(stderr,"***WARNING(in %s): failed reading the function value at node %d\n                 Continuing with sc->fval=0 for all points\n",__FUNCTION__,i);
      break;
    }
  }
  /*************************************************************/
  return sc;
}
/**********************************************************************/
/*!
 * \fn void haz_scomplex_print(scomplex *sc, const INT ns0,const char *infor)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void haz_scomplex_print(scomplex *sc, const INT ns0,const char *infor)
{
  // print simplicial complex, starting with ns0.
  INT i,j,in,in1;
  INT n=sc->n,n1=n+1,ns=sc->ns,nv=sc->nv,nbig=sc->nbig;
  if (ns0 < 0 || ns0>ns) return;
  fprintf(stdout,"\nN=%d,NBIG=%d, NS=%d, NV=%d\n",sc->n,sc->nbig,sc->ns,sc->nv);fflush(stdout);
  fprintf(stdout,"\n%s printout: %s\n",__FUNCTION__,infor);
  fprintf(stdout,"\nNODES list:\n");
  if(sc->parent){
    for(i=ns0;i<ns;i++){
      in1=i*n1;
      fprintf(stdout,"Element: %d ; vol=%e, Parent=%d; NODES=",i-ns0,sc->vols[i],sc->parent[i-ns0]);
      for(j=0;j<n1;j++)
	fprintf(stdout,"%d  ",sc->nodes[in1+j]);
      fprintf(stdout,"\n");
    }
  } else {
    for(i=ns0;i<ns;i++){
      in1=i*n1;
      fprintf(stdout,"Element: %d ; vol=%e, NODES=",i-ns0,sc->vols[i]);
      for(j=0;j<n1;j++)
	fprintf(stdout,"%d  ",sc->nodes[in1+j]);
      fprintf(stdout,"\n");
    }
  }
  fprintf(stdout,"\nNBR list:\n");
  if(sc->gen){
    for(i=ns0;i<ns;i++){
      in1=i*n1;
      fprintf(stdout,"Element: %d (%d) ; NBR=",i-ns0,sc->gen[i-ns0]);
      for(j=0;j<n1;j++)
	fprintf(stdout,"%d  ",sc->nbr[in1+j]-ns0);
      fprintf(stdout,"\n");
    }
  } else { 
    for(i=ns0;i<ns;i++){
      in1=i*n1;
      fprintf(stdout,"Element: %d ; NBR=",i-ns0);
      for(j=0;j<n1;j++)
	fprintf(stdout,"%d  ",sc->nbr[in1+j]-ns0);
      fprintf(stdout,"\n");
    }
  }
  //
  if(sc->bndry){
    for(i=0;i<nv;i++){
      in=i*nbig;
      fprintf(stdout,"Node: %d ; Code: %d ; COORDS=",i,sc->bndry[i]);
      for(j=0;j<nbig;j++){
	fprintf(stdout,"%e  ",sc->x[in+j]);
      }
      fprintf(stdout,"\n");
    }
  } else {
    for(i=0;i<nv;i++){
      in=i*nbig;
      fprintf(stdout,"Node: %d ; COORDS=",i);
      for(j=0;j<nbig;j++){
	fprintf(stdout,"%e  ",sc->x[in+j]);
      }
      fprintf(stdout,"\n");
    }
  }
  fflush(stdout);
  return;
}
/**********************************************************************/
/*!
 * \fn void haz_scomplex_free(scomplex *sc)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void haz_scomplex_free(scomplex *sc)
{
  if(sc->gen) free(sc->gen);
  if(sc->parent) free(sc->parent);
  if(sc->child0) free(sc->child0);
  if(sc->childn) free(sc->childn);
  if(sc->bndry) free(sc->bndry);
  if(sc->flags) free(sc->flags);
  if(sc->nodes) free(sc->nodes);
  if(sc->x) free(sc->x);
  if(sc->vols) free(sc->vols);
  //  fprintf(stdout,"9\n");fflush(stdout);
  if(sc->fval) free(sc->fval);
  if(sc->nbr) free(sc->nbr);
  if(sc) free(sc);
  return;
}
/**********************************************************************/
/*!
 * \fn void faces_cnt(subscomplex *subsc)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void faces_cnt(subscomplex *subsc)
{
  /*
     dim1 here is dim + 1: number of vertices in a simplex This
     function finds face count and creates element to face map. It
     preserves the ordering in sc->nbr, which means that in element
     el, face[el*dim1+k] is oposize to the vertex
     sc->nodes[ele*dim1+k] and is shared with the element number
     sitting at sc->nbr[ele*dim1+k]
     Such ordering (for sc->nbr) is needed for the refinement of an
n-dimensional simplicial grid so it is also used here to construct
"face-face", etc.
*/
  INT dim=subsc->nbig, dim1=dim+1;
  INT is,it,ir,di,dj,j,k;
  scomplex *sc=subsc->parent;
  INT ns=sc->ns;
  INT *nbr=sc->nbr, *elf=subsc->elf;
  INT nf=0;
  /*
      fprintf(stdout,"\n------Elements: vertices, n=%d, %d, %d\n",subsc->parent->ns,subsc->parent->nv,sc->n);fflush(stdout);
  */
  for(it=0;it<ns;it++){
    di=dim1*it;
    for(j=0;j<dim1;j++){
      is=nbr[di+j];
      if(it>is){
	elf[di+j]=nf;
	/*
	   in the neighboring element, find the place of i in the
	   neighboring list of j and put the face number in that spot
	   in elf[].
	*/
	if(is>=0){
	  dj=dim1*is;
	  ir=-1;
	  for(k = 0;k<dim1;k++){
	    ir=nbr[dj+k];
	    if(ir==it) { elf[dj+k]=nf; break; }
	  }
	  if(ir < 0) {
	    fprintf(stderr,						\
		    "\n\n*** ERROR in %s: incompatible neighbors: %d is neighbor of %d but %d is not a neighbor of %d\n", \
		    __FUNCTION__,is,it,it,is);
	    exit(127);
	  }
	}
	nf++;
      }
    }
  }
  // set the number of simplices in the subsc (faces)
  subsc->ns=nf;
  return;
}
/**********************************************************************/
/*!
 * \fn void area_face(INT dim, REAL fact, REAL *xf, REAL *sn,  REAL
 *	       *areas,REAL *volt,  void *wrk)
 *
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void area_face(INT dim, REAL fact, REAL *xf, REAL *sn,	\
	       REAL *areas,REAL *volt,			\
	       void *wrk)
{
  /*
     computes areas of all faces (n-1) dimensional simplices and their
     normal vectors for a given simplex. it also computes the volume
     of the simplex.

     work space: wrk should be at least dim*(dim+1) REALS and dim
     integers.
  */
  INT dim1 = dim+1,i,j,j1,ln,ln1;
  REAL *bt=(REAL *)wrk;
  REAL *piv=bt+dim*dim;
  INT *p = (INT *)(wrk+(dim*dim + dim)*sizeof(REAL));
  // construct bt using xf;
  for (j = 1;j<dim1;j++){
    ln=j*dim; ln1=ln-dim;
    for(i=0;i<dim;i++){
      bt[ln1+i] = xf[ln+i]-xf[i];
    }
  }
  //  print_full_mat(dim,dim,bt,"bt");
  if(lufull(1, dim, volt, bt,p,piv))
    *volt=0e0;
  else
    *volt=fabs(*volt)/fact;
  memset(sn,0,dim1*dim*sizeof(REAL));
  for (j = 1; j<dim1;j++){
    j1=j-1;
    ln=j*dim;
    sn[ln+j1]=-1.;
    solve_pivot(0, dim, bt, (sn+ln), p,piv);
    //areas[] contains the norm of the normal vectors first (this is 1/altitude);
    areas[j]=0.;
    for(i=0;i<dim;i++){
      sn[i]-=sn[ln+i];
      areas[j]+=sn[ln+i]*sn[ln+i];
    }
    areas[j]=sqrt(areas[j]);
  }
  areas[0]=0.; for(i=0;i<dim;i++) areas[0]+=sn[i]*sn[i];
  areas[0]=sqrt(areas[0]);
  //  print_full_mat(dim1,dim,sn,"snsn22");
  //rescale all normal vectors:
  for(j=0;j<dim1;j++){
    ln=j*dim;
    for(i=0;i<dim;i++) { sn[ln+i]/=areas[j]; }
    areas[j]=fabs((*volt)*areas[j])*((REAL )dim);
  }
  //  print_full_mat(dim1,1,areas,"areas");
  //  print_full_mat(dim1,dim,sn,"snsn33");
  return;
}
/**********************************************************************/
/*!
 * \fn void faces_attr(subscomplex *subsc)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void faces_attr(subscomplex *subsc)
{
  INT dim=subsc->nbig, dim1=dim+1;
  size_t nbits=dim*sizeof(REAL);
  scomplex *sc=subsc->parent;
  INT *nbr=sc->nbr, *elf=subsc->elf;
  INT is,it,di,i,j,k;
  INT *fflags=subsc->flags;
  INT l,node,bflag,j1,jf;
  //  find the minimal vertex bounday flag:
  INT minvflag=sc->bndry[0];
  for(i=1;i<sc->nv;i++)
    if(minvflag>sc->bndry[i])  minvflag=sc->bndry[i];
  for(jf=0;jf<subsc->ns;jf++){fflags[jf]=minvflag-1;}
  //  fprintf(stdout,"\nmin vertex flag=%d\n",minvflag);
  /* face_vertex map*/
  REAL *xf=(REAL *)calloc(dim*dim1,sizeof(REAL));
  REAL *snsn=(REAL *)calloc(dim*dim1,sizeof(REAL));
  REAL *areas=(REAL *)calloc(dim1,sizeof(REAL));
  // length of the work space needed : 3 matrices(dim*dim), mass
  // center(dim); pivots(dim),integer permutation(dim) ;
  INT nlength=3*dim*dim*sizeof(REAL)+2*dim*sizeof(REAL)+dim*sizeof(INT);
  void *wrk=(void *)calloc(nlength,sizeof(char));
  REAL fact=sc->factorial;
  REAL s;
  for(it=0;it<sc->ns;it++){
    di=dim1*it;
    //    fprintf(stdout,"\nel=%d, faces:",it);fflush(stdout);
    for(j=0;j<dim1;j++){
      is=nbr[di+j];
      jf=elf[di+j];
      //      fprintf(stdout,"%d(%d) ",jf,fflags[jf]);fflush(stdout);
      j1=jf*((subsc->n+1));
      // vertices in it which are opposite to jf
      if(fflags[jf]>=minvflag)
	{
	  fflags[jf]=0; //this is an interior face because we have
			//already been here.
	  continue;
	}
      l=0;
      bflag=minvflag;
      for(k=0;k<dim1;k++){
	node=sc->nodes[di+k];
	if(k==j) continue;
	subsc->nodes[j1+l]=node;
	if(abs(bflag)<abs(sc->bndry[node])) bflag=sc->bndry[node];
	l++;
      }
      fflags[jf]=bflag;
    }
    // use nbig below in the future!
    for(j=0;j<dim1;j++){
      node=sc->nodes[di+j];
      memcpy((xf+j*dim),(sc->x+node*dim),nbits);
    }
    /* xf[9]=0.; xf[10]=0.; xf[11]=0.; */
    /* xf[6]=1.; xf[7]=0.; xf[8]=0.; */
    /* xf[0]=0.; xf[1]=1.; xf[2]=0.; */
    /* xf[3]=0.; xf[4]=0.; xf[5]=1.; */
    area_face(dim, fact, xf, snsn,		\
	      areas,(sc->vols+it),		\
	      wrk);
    // get the signs of the normals to agree is it>is;
    for(j=0;j<dim1;j++){
      is=nbr[di+j];
      //      jf=elf[di+j];
      s=chk_sign(it,is);
      for(i=0;i<dim;i++){snsn[j*dim+i]*=s;}
      //      if(fflags[jf]==22) fprintf(stdout,"\n(%d,%d); s=%f",it,is,s);
    }
    //    print_full_mat(dim1,dim,snsn,"snsn");
    for(j=0;j<dim1;j++){
      jf=elf[di+j];
      memcpy(subsc->normals+jf*dim,(snsn+j*dim),nbits);
      subsc->areas[jf]=areas[j];
    }
    //    exit(11);
  }
  //  haz_scomplex_print(sc,0,"after attr");
  //  haz_subscomplex_print(subsc,0,__FUNCTION__);
  return;
}
/**********************************************************************/
/*!
 * \fn subscomplex *haz_subscomplex_init(scomplex *sc)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
subscomplex *haz_subscomplex_init(scomplex *sc)
{
  /*
     initialize s subcomplex; since we do not know how many simplices
     are in the complex, we need to find this out. But we know how to
     allocate memory for subsc->elf because this comes from the parent
     sc;
  */
  subscomplex *subsc=(subscomplex *) malloc(sizeof(subscomplex));
  subsc->nbig=sc->n;
  subsc->n=sc->n-1;
  subsc->parent = sc;
  /* allocate */
  subsc->elf=(INT *)calloc(sc->ns*(sc->n+1),sizeof(INT));
  /*
      here faces are counted and elf is constructed and these are
      stored in subsc->ns and subsc->elf;
  */
  faces_cnt(subsc);
  //  fprintf(stdout,"\nNUMBER of faces=%d dim_faces=%d\n",subsc->ns,subsc->n);fflush(stdout);
  /* allocate other arrays for the subcomplex */
  subsc->nodes=calloc(subsc->ns*(subsc->n+1),sizeof(INT));
  subsc->flags=calloc(subsc->ns,sizeof(INT));
  subsc->areas=calloc(subsc->ns,sizeof(REAL));
  subsc->normals=calloc((subsc->parent->n)*subsc->ns,sizeof(REAL));
  // fill in the rest of data; this below also computes the volume of
  // the elements and fills also sc->vols;
  faces_attr(subsc);
  return subsc;
}
/*
   read a simplicial complex (a mesh) in haz format and initialize the
   rest of the structure scomplex.
*/
/**********************************************************************/
/*!
 * \fn void haz_subscomplex_print(subscomplex *subsc, const INT ns0, const char *infor)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void haz_subscomplex_print(subscomplex *subsc, const INT ns0, const char *infor)
{
  // print simplicial subcomplex, starting with ns0.
  INT i,j,in1;
  scomplex *sc=subsc->parent;
  if (ns0 < 0 || ns0>subsc->ns) return;
  fprintf(stdout,"\nSubsc info: %s\n",infor);fflush(stdout);
  fprintf(stdout,"\nFace list:\n");
  for(i=ns0;i<sc->ns;i++){
    in1=(sc->n+1)*i;
    fprintf(stdout,"\nt= %d ; vol(el)=%e, faces=",	\
  	    i-ns0,sc->vols[i]); fflush(stdout);
    for(j=0;j<(sc->n+1);j++){
      fprintf(stdout,"%d  ",subsc->elf[in1+j]);fflush(stdout);
    }
  }
  fprintf(stdout,"\n");
  fprintf(stdout,"\n%d %d vertex list:\n",subsc->ns,sc->nv);fflush(stdout);
  for(i=ns0;i<subsc->ns;i++){
    in1=i*(subsc->n+1);
    fprintf(stdout,"f= %d (flag=%d) area(f)=%e; vert=",	\
  	    i-ns0,subsc->flags[i-ns0],subsc->areas[i-ns0]);  fflush(stdout);
    for(j=0;j<(subsc->n+1);j++){
      fprintf(stdout,"%d  ",subsc->nodes[in1+j]);  fflush(stdout);
    }
    fprintf(stdout,"\n");  fflush(stdout);
  }
  return;
}
/**********************************************************************/
/*!
 * \fn
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void haz_subscomplex_free(subscomplex *subsc)
{
  if(subsc->elf) free(subsc->elf);
  if(subsc->flags) free(subsc->flags);
  if(subsc->nodes) free(subsc->nodes);
  if(subsc->normals) free(subsc->normals);
  if(subsc->areas) free(subsc->areas);
  if(subsc) free(subsc);
  return;
}
/**********************************************************************/
/*!
 * \fn static unsigned int cmp_simplex(INT n, INT sim1, INT sim2, INT
 *			 *sv1, INT *sv2, INT *stos1, INT *stos2)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
static unsigned int cmp_simplex(INT n, INT sim1, INT sim2,	\
			 INT *sv1, INT *sv2, INT *stos1, INT *stos2)
{
  //sim1 and sim2 are pointers to the neighboring lists of two
  //simplices stos1 and stos2 are pointers to rows in
  //simplex-to-simplex incidence; sv1 and sv2 are pointers to rows in
  //simplex-vertex incidence for two neighboring simplices.
  //
  //rearrange to meet the structural condition.  compare two sets of
  //length n. they should overlap at exactly n-1 members.  If they do
  //not, 0 is returned, otherwise 1. no check if the members of the
  //sets are distinct (they should be for the purpose of comparing two
  //simplices.
  unsigned int fnd;
  unsigned int nf=0;
  INT i1,k1=-10,i2,k2=-10;
  INT i,j;
  for (i=0;i<n;i++){
    fnd=0;
    i1=sv1[i];
    for (j = 0;j < n;j++){
      if(i1 == sv2[j]){
	fnd=1;
	break;
      }
    }
    if(fnd) {
      continue;
    } else {
      // not found
      nf++;
      if(nf > 1) return 0;
      k1=i;
    }
  }
  // same with sv1 and sv2 swapped.
  nf=0;
  for (i=0;i<n;i++){
    fnd=0;
    i2=sv2[i];
    for (j=0;j<n;j++){
      if(i2 == sv1[j]){
	fnd=1;
	break;
      }
    }
    if(fnd) {
      continue;
    } else {
      // not found
      nf++;
      if(nf > 1) return 0;
      k2=i;
    }
  }
  /* NOW put the neightbors at the right places ******* */
  if(k1<0||k2<0){
    fprintf(stderr,"\n***ERROR in %s; k1=%d,k2=%d\n\n",__FUNCTION__,k1,k2);
    exit(65);
  } else {
    stos1[k1]=sim2;
    stos2[k2]=sim1;
    return 1;
  }
}
/**********************************************************************/
/*!
 * \fn void find_nbr(INT ns,INT nv,INT n,INT *sv,INT *stos)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void find_nbr(INT ns,INT nv,INT n,INT *sv,INT *stos)
{
  // find neighboring list
  INT* ivs=NULL, *jvs=NULL, *stosi=NULL, *stosk=NULL, *svi=NULL, *svk=NULL;
  INT i,j,k,jp,jia,jib,iabeg,iaend,ibbeg,ibend,kn1;
  INT n1 = n+1,nv1=nv+1,nsv=ns*n1;
  /* allocate */
  ivs=(INT *) calloc(nv1,sizeof(INT));
  jvs=(INT *) calloc(nsv,sizeof(INT));
  /* transpose to get the vertex--simplex incidence */
  for (i = 0; i < nv1; ++i)
    ivs[i] = 0;
  //
  for (k = 0; k < ns; ++k) {
    kn1=k*n1;
    for (i = 0; i < n1; i++) {
      j = sv[kn1+i] + 2;
      if (j < nv1) ivs[j]++;
    }
  }
  ivs[0] = 0;
  ivs[1] = 0;
  if (nv != 1) {
    for (i= 2; i < nv1; ++i) {
      ivs[i] += ivs[i-1];
    }
  }
  for (i = 0; i < ns; ++i) {
    iabeg = i*n1;
    iaend = i*n1+n1;
    for (jp = iabeg; jp < iaend; ++jp) {
      j = sv[jp]+1;
      k = ivs[j];
      jvs[k] = i;
      ivs[j] = k+1;
    }
  }
  /* for (i = 0; i < nv; ++i) { */
  /*   iabeg = ivs[i]; */
  /*   iaend = ivs[i+1]; */
  /*   fprintf(stdout,"row %d: ",i+1); */
  /*   for (jia = iabeg; jia < iaend; ++jia) { */
  /*     fprintf(stdout,"%d ",jvs[jia]+1); */
  /*   } */
  /*   fprintf(stdout,"\n"); */
  /* } */
  /*
   */
  INT *icp=(INT *) calloc(ns,sizeof(INT));
  for (i = 0; i < ns; ++i) icp[i] = -1;
  for (i = 0; i < nsv; ++i) stos[i] = -1;
  for (i = 0; i < ns; ++i) {
      iabeg = i*n1;
      iaend = iabeg + n1;
      stosi=stos+iabeg; svi=sv+iabeg;
      for (jia = iabeg; jia < iaend; ++jia) {
	j = sv[jia]; // vertex number
	ibbeg = ivs[j];
	ibend = ivs[j+1];
	// loop over all simplices with this j as a vertex.
	for (jib = ibbeg; jib < ibend; ++jib) {
	  k = jvs[jib];
	  if(k<=i) continue;
	  if (icp[k] != i) {
	    icp[k] = i;
	    kn1=k*n1;
	    stosk=stos+kn1; svk=sv+kn1;
	    if(!cmp_simplex(n1,i,k,svi,svk,stosi,stosk)) continue;
	  }//if
	} //for
      } //for
    } //for (i...
    if (icp) free(icp);
    if (ivs) free(ivs);
    if (jvs) free(jvs);
    //  }
    return;
}
/**********************************************************************/
/*!
 * \fn INT haz_add_simplex(INT is, scomplex *sc,REAL *xnew,INT
 *		    ibnew,INT csysnew,INT nsnew, INT nvnew)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
INT haz_add_simplex(INT is, scomplex *sc,REAL *xnew,INT ibnew,INT csysnew, \
		    INT nsnew, INT nvnew)
{
  /* adds nodes and coords as well */
  INT n=sc->n, nbig=sc->nbig, n1 = n+1,nv=sc->nv;//ns=sc->ns;
  INT ks0=sc->child0[is], ksn=sc->childn[is];
  INT isc0=ks0*n1, iscn=ksn*n1;
  //  INT *dsti,*srci;
  INT j,j0,jn;
  REAL *dstr;
  /* nodes  AND neighbors */
  sc->nbr=realloc(sc->nbr,(nsnew*n1)*sizeof(INT));
  sc->nodes=realloc(sc->nodes,(nsnew*n1)*sizeof(INT));
  for(j=0;j<n1;j++){
    j0=isc0+j;
    jn=iscn+j;
    //    fprintf(stdout,"\nNSN1 %d, NSN1J %d\n",j0,jn);   fflush(stdout);
    sc->nodes[j0]=-1;
    sc->nbr[j0]=-1;
    sc->nodes[jn]=-1;
    sc->nbr[jn]=-1;
  }
  //new vertex (if any!!!)
  if(nvnew != nv) {
    sc->x=realloc(sc->x,(nvnew*nbig)*sizeof(REAL));
    dstr=(sc->x+nv*nbig); memcpy(dstr,xnew,n*sizeof(REAL));
    if(xnew) free(xnew);
    sc->bndry=realloc(sc->bndry,(nvnew)*sizeof(INT));
    sc->bndry[nv]=ibnew;
    sc->csys=realloc(sc->csys,(nvnew)*sizeof(INT));
    sc->csys[nv]=csysnew;
    sc->fval=realloc(sc->fval,(nvnew)*sizeof(REAL));
    sc->fval[nv]=0.;
  }
  //generation
  sc->gen=realloc(sc->gen,(nsnew)*sizeof(INT));
  sc->gen[ks0]=sc->gen[is]+1;
  sc->gen[ksn]=sc->gen[is]+1;
  //marked
  sc->marked=realloc(sc->marked,(nsnew)*sizeof(INT));
  sc->marked[ks0]=sc->marked[is];
  sc->marked[ksn]=sc->marked[is];
  //flags
  sc->flags=realloc(sc->flags,(nsnew)*sizeof(INT));
  sc->flags[ks0]=sc->flags[is];
  sc->flags[ksn]=sc->flags[is];
  //parents
  sc->parent=realloc(sc->parent,(nsnew)*sizeof(INT));
  sc->parent[ks0]=is;
  sc->parent[ksn]=is;
  //child0
  sc->child0=realloc(sc->child0,(nsnew)*sizeof(INT));
  sc->child0[ks0]=-1;
  sc->child0[ksn]=-1;
  //childn
  sc->childn=realloc(sc->childn,(nsnew)*sizeof(INT));
  sc->childn[ks0]=-1; sc->childn[ksn]=-1;
  //scalars
  sc->ns=nsnew;
  sc->nv=nvnew;
  return 0;
}
/**********************************************************************/
/*!
 * \fn INT haz_refine_simplex(scomplex *sc, const INT is, const INT it)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
INT haz_refine_simplex(scomplex *sc, const INT is, const INT it)
{
  INT n=sc->n, nbig=sc->nbig,ns=sc->ns,nv=sc->nv;
  INT nsnew=ns,nvnew=nv;
  INT itype,nodnew=-10;
  INT n1=n+1,j,i,p,p0,pn,isn1,snbri,snbrp,snbrn,snbr0;
  INT jt,jv0,jvn,ks0,ksn,s0nbri,snnbri;//,isn;
  REAL *xnew;
  INT csysnew,ibnew;
  if(is<0) return 0;
  if(sc->child0[is] >= 0) return 0;
  //  isn=is*n;
  isn1=is*n1;
  haz_scomplex_print(sc,-10,__FUNCTION__);
  for (i=1;i<n;i++){
    snbri=sc->nbr[isn1+i] ; // the on-axis neighbor.
    if(snbri<0) continue; //this is a boundary
    if (sc->gen[snbri]<sc->gen[is]){//this was wrong in the traxler's paper
      haz_refine_simplex(sc,snbri,-1);
      nsnew=sc->ns;
      nvnew=sc->nv;
      haz_scomplex_print(sc,-10,"AXIS NEIGHBOR");
    }
  }
  if (it<0){
    xnew = (REAL *)calloc(nbig,sizeof(REAL));
    jv0=sc->nodes[isn1];
    jvn=sc->nodes[isn1+n];
    // here we can store the edge the new vertex comes from
    for(j=0;j<nbig;j++){
      xnew[j] = 0.5*(sc->x[jv0*nbig+j]+sc->x[jvn*nbig+j]);
    }
    if(sc->bndry[jv0] > sc->bndry[jvn])
      ibnew=sc->bndry[jv0];
    else
      ibnew=sc->bndry[jvn];
    if(sc->csys[jv0] < sc->csys[jvn])
      csysnew=sc->csys[jv0];
    else
      csysnew=sc->csys[jvn];
    /* we have added a vertex, let us indicate this */
    nodnew = nvnew;
    nvnew++;
    /* fprintf(stdout,"\nnew vertex %d: is=%d  it=%d, edge (%d,%d)\n",nodnew,is,it,jv0,jvn); fflush(stdout); */
  } else {
    //    jx=    newvertex=t->child0->vertex[1];
    jt=sc->child0[it]; // child simplex number
    jv0 = sc->nodes[jt*n1+1]; // global vertex number of vertex 1.
    xnew = (sc->x + jv0*nv); /* xnew is the same pointer, the vertex
				has already been added. */
    nodnew=jv0;
    ibnew=sc->bndry[nodnew];
    csysnew=sc->csys[nodnew];
    /* fprintf(stdout,"\n no new vertex: is=%d  it=%d; vertex=%d\n",is,it,jv0); fflush(stdout); */
  }
  ks0=nsnew;
  sc->child0[is]=ks0; // child0 simplex number
  nsnew++;
  ksn=nsnew;
  sc->childn[is]=ksn; // childn simplex number
  nsnew++;
  /*
    Add two new simplices and initialize their parents, etc
  */
  haz_add_simplex(is,sc,xnew,ibnew,csysnew,nsnew,nvnew);
  /*
    Initialize all vertex pointers of the children according to the
    scheme. Also initialize all pointers to bring simplices as long
    as they do not require a recursive call for subdivision.  Always
    remember: The mesh is supposed to meet the structural condition when
    this routine is called, and it will meet the structural condition
    again after this routine has terminated.
  */
  haz_scomplex_print(sc,-10,"AFTER ADD");
  INT isc0=-100,iscn=-100;
  //ks0 = sc->child0[is];
  //ksn = sc->childn[is];
  isc0=ks0*n1;
  iscn=ksn*n1;
  sc->nodes[isc0+0]=sc->nodes[isn1];   // nodes[is][0]
  sc->nodes[iscn+0]=sc->nodes[isn1+n]; // nodes[is][n1] nodes is (ns x n1)
  /*backup:   sn->nodes[1]=s0->nodes[1]=nodnew;*/
  sc->nodes[iscn+1]=sc->nodes[isc0+1]=nodnew;
  sc->nbr[isc0]=sc->childn[is];
  sc->nbr[iscn]=sc->child0[is];
  /* fprintf(stdout,"\nOOOOOOOOOO: Refining: %d ; child0=%d child1=%d\n",is+1,ks0+1,ksn+1);  */
  /* fprintf(stdout,"\nOOOOOOOOOO: sc->ns=%d ; nsnew=%d\n",sc->ns,nsnew);  */
  /* fprintf(stdout,"\nOOOOOOOOOO: i0=%d;in=%d\n",sc->nodes[isn1]+1,sc->nodes[isn1+n]+1); fflush(stdout); */
  /*
    We know the structure of the nbrs children, if existent already,
    thanks to the structural condition.
  */
  snbr0=sc->nbr[isn1];
  snbrn=sc->nbr[isn1+n];
  if(snbrn>=0){
    if (sc->child0[snbrn]>=0){
      //      s0->nbr[1]=sc->child0[snbrn];
      sc->nbr[isc0+1]=sc->child0[snbrn];
      /* fprintf(stdout,"\nchild: %d\n", sc->child0[snbrn]+1); */
    } else {
      //      s0->nbr[1]=snbrn;
      sc->nbr[isc0+1]=snbrn;
      /* fprintf(stdout,"\nNO CHILD: %d\n", snbrn+1); */
    }
  }
  if(snbr0>=0){
    if(sc->childn[snbr0]>=0)
      //      sn->nbr[1]=sc->childn[snbr0];
      sc->nbr[iscn+1]=sc->childn[snbr0];
    else
      //      sn->nbr[1]=snbr0;
      sc->nbr[iscn+1]=snbr0;
  }
  haz_scomplex_print(sc,-10,"CHECK 1");
  /* Compute the simplex type. */
  itype=(sc->gen[is]) % n;
  // for each vertex p=1..n-1 of the parent simplex S
  for (p=1;p<n;p++)
    {
      /*
	 p0 is the index of S->vertex[p] in child0
	 pn is the index of S->vertex[p] in childn
      */
      pn = p+1;
      p0 = p+1;
      if (p > itype) pn=n-p+itype+1;
      /*       YYYYYYYYYYYYYYYYYYYYYY */
      /* initialize children % vertex according to structural cond */
      //      sn->nodes[pn]=s0->nodes[p0]=sc->nodes[isn1+p];
      sc->nodes[isc0+p0]=sc->nodes[isn1+p];
      sc->nodes[iscn+pn]=sc->nodes[isn1+p];
      snbrp=sc->nbr[isn1+p]; /* s->nbr[p] */
      if(snbrp<0) continue;
      if (sc->child0[snbrp]>=0) {
	/*
	   S-> nbr[p] is subdivided. The corresponding nbr pointers of the
	   children of S should then point to S nbr[p] -> children (structural
	   condition). It might however be that the vertices 0 and n of S have
	   exchanged local indices in the nbr. That has to be tested.
	*/
	if (sc->nodes[snbrp*n1+0]==sc->nodes[isn1+0]) {
	  //	  s0->nbr[p0]=sc->child0[snbrp];
	  //	  sn->nbr[pn]=sc->childn[snbrp];
	  sc->nbr[isc0+p0]=sc->child0[snbrp];
	  sc->nbr[iscn+pn]=sc->childn[snbrp];
	} else {
	  //	  s0->nbr[p0]=sc->childn[snbrp];
	  //	  sn->nbr[pn]=sc->child0[snbrp];
	  sc->nbr[isc0+p0]=sc->childn[snbrp];
	  sc->nbr[iscn+pn]=sc->child0[snbrp];
	}
	/* fprintf(stdout,"is=%d (node0,node1,node2)=(%d,%d,%d)\n",is+1,sc->nodes[isn1]+1,sc->nodes[isn1+1]+1,sc->nodes[isn1+2]+1); */
	/* fprintf(stdout,"is_nbr=%d (node0,node1,node2)=(%d,%d,%d)\n",snbrp+1,sc->nodes[snbrp*n1]+1,sc->nodes[snbrp*n1+1]+1,sc->nodes[snbrp*n1+2]+1); */
	/* fprintf(stdout,"xxxxxxxxxx, (ks0,ksn)=(%d,%d), (nod_nbrp(0),nod_is(0))=(%d,%d)\n",ks0+1,ksn+1,sc->nodes[snbrp]+1,sc->nodes[isn1]+1); */
	/* fprintf(stdout,"we are herezzzzzzzzzz, (ks0,ksn)=(%d,%d), (p0,pn)=(%d,%d), (nbrc0p0,nbrcnpn)=(%d,%d)\n",ks0+1,ksn+1,p0,pn,sc->nbr[isc0+p0]+1,sc->nbr[iscn+pn]+1); */
	/* fprintf(stdout,"we are here1111111; is=%d, is_nbr[p]=%d:(c0,cn)=(%d,%d)\n",is+1,snbrp+1,sc->child0[snbrp]+1,sc->childn[snbrp]+1); */
      } else {
	/*
	  s->nbr[p] is not subdivided. The corresponding neighbor pointers
	  of the children of s are now set to s->nbr[p], which will be
	  corrected later, when this simplex is divided as well.
	*/
	//	sn->nbr[pn]=snbrp;
	//	s0->nbr[p0]=snbrp;
	sc->nbr[isc0+p0]=snbrp;
	sc->nbr[iscn+pn]=snbrp;
      }
    }
  //  fprintf(stdout,"\nis,it=(%d,%d); snbrn=%d,snbrn0=%d\n",is+1,it+1,snbrn+1,snbr0+1);
  haz_scomplex_print(sc,-10,"CHECK 2"); fflush(stdout);
  for(i=0;i<n1;i++) {
     /*
	The children OF NEIGHBORS, if existent, still point to s as
	one of their nbrs. This contradicts the structural condition
	and is corrected here.
     */
    //    s0nbri=s0->nbr[i];
    //    snnbri=sn->nbr[i];
    //    fprintf(stdout,"\nisc0i=%d, iscni=%d",isc0+i,iscn+i);   fflush(stdout);
    s0nbri=sc->nbr[isc0+i];    /*s->child0->neighbor[i]*/
    snnbri=sc->nbr[iscn+i]; /*s->childn->neighbor[i]*/
    //    fprintf(stdout,"\n\nXXXnsnew,s0nbri,snnbri,ns: %i %i %i %i\n\n",nsnew,snnbri,s0nbri,sc->ns); fflush(stdout);
    if(s0nbri>=0){
      if(s0nbri >=sc->ns) {
	fprintf(stdout,"\n\nSTOPPING: nsnew,s0nbri,snnbri,ns: %i %i %i %i\n\n",nsnew,snnbri,s0nbri,sc->ns); fflush(stdout);
	exit(222);
      }
      //      if(sc->gen[s0nbri]==s0->gen)
      if(sc->gen[s0nbri]==sc->gen[ks0])
	/* sc->nbr[s0nbri*n1+i] = s->child0->neighbor[i]->neighbor[i]*/
	sc->nbr[s0nbri*n1+i]=ks0;
    }
    if(snnbri>=0){
      if(snnbri >=sc->ns) {
	fprintf(stdout,"\n\nSTOPPING2: s0nbri,snnbri,ns: %i %i %i %i\n",nsnew,snnbri,s0nbri,sc->ns); fflush(stdout);
	exit(223);
      }
      //      if(sc->gen[snnbri]==sn->gen)
      if(sc->gen[snnbri]==sc->gen[ksn])
	/* sc->nbr[snnbri*n1+i] = s->childn->neighbor[i]->neighbor[i]*/
	sc->nbr[snnbri*n1+i]=ksn;
    }
  }
  haz_scomplex_print(sc,-10,"CHECK 3");
  /*
     NOW call the on-axis nbrs for refinement, passing to them a
     pointer to this simplex S so they find our new vertex.  Skip the
     neighbors opposite to x0 and xn, they do not have to be refined
     and refine the rest of the "on-axis" neighbors */
  for(i=1 ; i < n; i++){
    //    fprintf(stdout,"\n%%trying to refine also (%i) coming from  (%i)\n ",sc->nbr[isn1+i]+1, is+1); fflush(stdout);
    haz_refine_simplex(sc, sc->nbr[isn1+i], is);
        haz_scomplex_print(sc,-10,"\nEND OF haz_refine\n");
  }
  return 0;
}
/******************************************************************/
/*!
 * \fn void refine(const INT ref_levels, scomplex *sc,ivector *marked)
 *
 * \brief Refines a simplicial complex.
 *
 * \param sc: scomplex containing the whole hierarchy of refinements
 *
 * \param marked: input ivector containing all simplices from the
 *               finest level which are marked for refinement; If
 *               marked is null then uniformly renfines the grid
 *               ref_levels.
 *
 *
 * \param ref_levels number of refinements. If marked is not null,
 * then this is ignored and only one refinement is done;
 *
 * \return void
 *
 */
void refine(const INT ref_levels, scomplex *sc,ivector *marked)
{
  if(ref_levels<=0) return;
  /*somethind to be done*/
  INT j=-1,i,nsold,print_level=0,nsfine=-1;
  if(!sc->level){
    /* sc->level this is set to 0 in haz_scomplex_init */
    /* form neighboring list on the coarsest level */
    find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
    //    haz_scomplex_print(sc,0,__FUNCTION__);  fflush(stdout);
    INT *wrk=calloc(5*(sc->n+2),sizeof(INT));
    /* construct bfs tree for the dual graph */
    abfstree(0,sc,wrk,print_level);
    //    haz_scomplex_print(sc,0,__FUNCTION__);fflush(stdout);
    //    exit(100);
    if(wrk) free(wrk);
  }
  if((!marked)){
    // just refine everything that was not refined:
    for (i=0;i<ref_levels;i++){
      nsold=sc->ns;
      for(j = 0;j < nsold;j++)
	if((sc->child0[j]<0||sc->childn[j]<0))
	  haz_refine_simplex(sc, j, -1);
      sc->level++;
    }
    for(j=0;j<sc->ns;j++) sc->marked[j]=TRUE; // not sure we need this.
    // we are done here;
    return;
  } else if(!sc->level){
    // we have not refined anything yet and marked is set, so
    if((marked->row)&&(marked->val))
      for(j=0;j<sc->ns;j++) sc->marked[j]=marked->val[j];
    else {
      //      issue a warning and mark everything for refinement;
      for(j=0;j<sc->ns;j++) sc->marked[j]=TRUE;
    }
  } else {
    /* we come here if we have refined few times and in such case we
       need to re-mark our simplices on the finest level using the
       values of marked at abs(sc->child0[j]+1) which were set from
       the previous visit here */
    for(j=0;j<sc->ns;j++) {
      if(sc->child0[j]<0||sc->childn[j]<0)
	nsfine=abs(sc->child0[j]+1);
      // if(nsfine>sc->ns) issue an error and exit;
	sc->marked[j]=marked->val[nsfine];
    }
  }
  /*
   * refine everything that is marked on the finest level and is
   * not yet refined: (marked>0 and child<0)
   */
  nsold=sc->ns;
  for(j = 0;j < nsold;j++)
    if(sc->marked[j] && (sc->child0[j]<0||sc->childn[j]<0))
      haz_refine_simplex(sc, j, -1);
  sc->level++;
  //  fprintf(stdout,"\n.%d.\n",sc->level);
  //  fprintf(stdout,"\n");
  //  haz_scomplex_print(sc,0,__FUNCTION__);  fflush(stdout);
}
/******************************************************************/
/*!
 * \fn void sc2mesh(scomplex *sc,mesh_struct *mesh)
 *
 * \brief copies scomplex structure to mesh struct (not all but the
 *        bare minimum needed to define mesh_struct.
 *
 * \param scomplex sc;
 *
 * \param mesh_struct mesh;
 *
 *
 */
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
mesh_struct *sc2mesh(scomplex *sc)
{
  /*********************************************************************/
  /* INPUT: pointer to a simplicial complex sc; returns mesh_struct
     mesh.
 
     Store the finest mesh in mesh_struct structure.  sc has all the
    hierarchy, mesh_struct will have only the last mesh.
  */
  /*********************************************************************/
  /* copy the final grid at position 1*/
  INT ns=0,nv=sc->nv,n1=sc->n+1,dim=sc->n,dimbig=sc->nbig;
  if(dimbig!=dim)
    return NULL;
  INT jk=-10,k=-10,j=-10,i=-10;
  mesh_struct *mesh=malloc(1*sizeof(mesh_struct));
  initialize_mesh(mesh);
  /*
    count the number of elements on the last level of refinement.
    On the last grid are all simplices that were not refined, so
    these are the ones for which child0 and childn are not set.
  */
  ns=0;
  for (j=0;j<sc->ns;j++)
    if(sc->child0[j]<0 || sc->childn[j]<0) ns++;
  /*Update mesh with known quantities*/
  mesh->dim = sc->n;
  mesh->nelm = ns; //do not ever put sc->ns here
  mesh->nv=nv; 
  mesh->nconn_reg = sc->cc; // dummy argument yet.
  mesh->nconn_bdry = sc->bndry_cc;// dummy argument
  mesh->v_per_elm = n1;
  mesh->nelm=ns;
  /*Allocate all pointers needed to init the mesh struct*/
  mesh->el_flag = (INT *)calloc(ns,sizeof(INT));
  mesh->el_vol = (REAL *)calloc(ns,sizeof(REAL));
  mesh->v_flag = (INT *)calloc(nv,sizeof(INT));
  mesh->cv=allocatecoords(nv,dim);
  mesh->el_v=(iCSRmat *)malloc(1*sizeof(iCSRmat));
  mesh->el_v->row=mesh->nelm;
  mesh->el_v->col=mesh->nv;
  mesh->el_v->nnz=mesh->nelm*mesh->v_per_elm;
  mesh->el_v->IA = (INT *)calloc((ns+1),sizeof(INT));
  mesh->el_v->JA = (INT *)calloc(ns*n1,sizeof(INT));
  mesh->el_v->val=NULL;
  INT chk=(INT )(!(mesh &&					  \
		 mesh->el_flag && mesh->el_vol && mesh->v_flag && \
		 mesh->cv && mesh->cv->x && \
		 mesh->el_v && mesh->el_v->IA && mesh->el_v->JA
		   ));
  if(chk){
    fprintf(stderr,"\nCould not allocate memory for mesh in %s\n *** RETURNING a NULL pointer to mesh", \
	    __FUNCTION__); 
    return NULL;
  }
  /********************************************************************/
  /*element flag and element volumes; volumes are recomputed later in
    build_mesh_all()*/
  INT nsfake=0;// this one must be ns at the end. 
  for (j=0;j<sc->ns;j++){
    if(sc->child0[j]<0 || sc->childn[j]<0){
      mesh->el_flag[nsfake]=sc->flags[j];
      mesh->el_vol[nsfake]=sc->vols[j];
      nsfake++;
    }
  }
  /*boundary flags*/
  mesh->nbv=0;
  for (j=0;j<mesh->nv;j++){
    if(sc->bndry[j]!=0) mesh->nbv++;
    mesh->v_flag[j]=sc->bndry[j];
  }
  /*Coordinates*/
  for(j=0;j<mesh->dim;j++){
    for(i=0;i<nv;i++){
      mesh->cv->x[j*nv+i]=sc->x[i*sc->n+j];
    }
  }
  // el_v
  mesh->el_v->IA[0]=0;
  for (j=0;j<mesh->nelm;j++)
    mesh->el_v->IA[j+1]=mesh->el_v->IA[j]+n1;
  jk=0;
  for (j=0;j<sc->ns;j++){
    /*  copy el_v map using only the top grid;    */
    if(sc->child0[j]<0 || sc->childn[j]<0){
      for (k=0;k<n1;k++)
        memcpy(mesh->el_v->JA+jk*n1,sc->nodes+j*n1,n1*sizeof(INT));
      jk++;
    }
  }
  fprintf(stdout,"\n%%After %d levels of refinement:\tvertices=%d ; simplices=%d in dim=%d;\n components(bdry):%d\n",   sc->level,mesh->nv,mesh->nelm,mesh->dim,mesh->nconn_bdry); fflush(stdout);
  return mesh;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*EOF*/
