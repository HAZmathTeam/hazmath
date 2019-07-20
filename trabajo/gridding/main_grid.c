#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <getopt.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>
#include "hazmath.h"
#include "grid_defs.h"
#ifndef INT
#define INT int
#endif
#ifndef REAL
#define REAL double
#endif
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
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
  INT ii,i,j,im,nb,aj,bi;
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
  }
  //  fprintf(stdout,"\n");
  return nb;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
scomplex *umesh(const INT dim,				\
		INT *nd, cube2simp *c2s,		\
		INT *isbndf, INT *codef,		\
		INT elflag,				\
		const INT intype);
void unirefine(INT *nd,scomplex *sc);
macrocomplex *set_mmesh(input_grid *g0,		\
			cube2simp *c2s,		\
			INT *wrk);
void set_edges(input_grid *g0,cube2simp *c2s);
INT set_ndiv_edges(input_grid *g,		\
		   input_grid *g0,		\
		   cube2simp *c2s,		\
		   INT **nd,			\
		   const INT iter);
void map2mac(scomplex *sc,cube2simp *c2s, input_grid *g);
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
void *fix_grid(macrocomplex *mc,		\
	       scomplex **scin,			\
	       cube2simp *c2s,			\
	       input_grid *g0)
{
  /****************************************************************/
  /* 
     this piece of code removes all repeated vertices in macro
     elements making the numbering of the  vertices a valid numbering;
  */
  if(mc->nel==1) return;
  scomplex *sc,*scp;
  INT dim=c2s->n,nel=mc->nel,dim1=c2s->n+1;
// at most dim (n-1)dimensional faces may intersect to form a vertex
  INT *mp = (INT *)calloc(dim1,sizeof(INT));
  INT *m = (INT *)calloc(dim1,sizeof(INT));
  INT *mi = (INT *)calloc(dim1,sizeof(INT));
  INT *mip = (INT *)calloc(dim1,sizeof(INT));
  INT *vertk = (INT *)calloc(dim1,sizeof(INT));
  INT *vertj = (INT *)calloc(dim1,sizeof(INT));
  INT *mm = (INT *)calloc(dim1,sizeof(INT));
  INT *ti=vertj,*tip=vertk;
  //scalars
  INT nv,kel, jel,nk,nj,flag;
  INT numv,kf,kfp,kz,kdim,kj,ijk;
  INT i,j,k,iaa, iab;
  INT isbf,bf,cf,knnz;
  INT *nodesj=NULL,*nodesk=NULL;
  // place holders (short hand)
  iCSRmat *fel2el=mc->fullel2el;
  INT **nd=mc->nd;
  INT **elneib=mc->elneib;
  INT **el2fnum=mc->el2fnum;
  INT *iindex,*iindexp;
  INT nvc=c2s->nvcube,nvc1=c2s->nvcube+1;
  for(kz=(dim-1);kz>=0;kz--){
    //serching for intersection of kz-dimensional faces.
    kdim=(1<<kz);
    //    fprintf(stdout,"\nkdim=%d",kdim);
    for(knnz=0;knnz<mc->bfs->nnz;knnz++){    
      kel=mc->bfs->JA[knnz];
      scp=scin[kel];
      iindexp=mc->iindex[kel];
      nodesk=(g0->mnodes+kel*nvc1);
      iaa=fel2el->IA[kel];
      iab=fel2el->IA[kel+1];
      for (kj=iaa;kj<iab;kj++){
	jel=fel2el->JA[kj];
	nodesj=(g0->mnodes+jel*nvc1);
	numv=0;
	if(fel2el->val[kj]!=kdim) continue;
	sc=scin[jel];
	iindex==mc->iindex[jel];
	fprintf(stdout,"\nkel=%d,jel=%d,kz=%d\n",kel,jel,kz);
	// find now how many times these two elements intersect:
	for(k=0;k<c2s->nvcube;k++){
	  i=nodesk[k];
	  i=locate0(i,nodesj,c2s->nvcube);
	  if(i<0) continue;
	  vertk[numv]=k;
	  vertj[numv]=i;
	  numv++;
	}
	nk=locate1(mip,vertk,numv,c2s->faces,c2s->nf,c2s->nvface);
	nj=locate1(mi,vertj,numv,c2s->faces,c2s->nf,c2s->nvface);
	if((nk == nj) && nk==(dim-kz) && nj==(dim-kz)) {
	  /* print_full_mat_int(1,nvc,nodesk,"nodesk");  */
	  /* print_full_mat_int(1,numv,vertk,"vertk");  */
	  /* print_full_mat_int(1,nvc,nodesj,"nodesj");  */
	  /* print_full_mat_int(1,numv,vertj,"vertj");  */
	  for(j=0;j<nj;j++){
	    /* 
	       for two macroelements kel and jel, from faces mi[] and
	       mip[] whose intersection forms the intersection kel and jel
	       macroelements get the ti[] and tip[] arrays which
	       describe on which boundary using the number of
	       divisions we have the corresponding vertices lying
	    */
	    if(mi[j]<dim){
	      mi[j]=dim-(mi[j]+1);
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
	  for(kf=0;kf<sc->nv;kf++){
	    coord_lattice(m,dim,kf,sc->nv,nd[jel]);
	    flag=1;
	    for(j=0;j<nj;j++) if(m[mi[j]]!=ti[j]){flag=0;break;}
	    if(flag){
	      memcpy(mp,m,dim*sizeof(INT));
	      for(k=0;k<nk;k++)mp[mip[k]]=tip[k];
	      kfp=num_lattice(mp,dim,nd[kel]);
	      sc->bndry[kf]=100;
	      fprintf(stdout,"Skipping %d<-->%d",kf,kfp);
	      fprintf(stdout,": oldx=(");
	      for(ijk=0;ijk<dim;ijk++)
		fprintf(stdout,"%.5e ", scp->x[kfp*dim+ijk]);
	      sc->bndry[kfp]=-100;
	      fprintf(stdout,"); newx=[");
	      for(ijk=0;ijk<dim;ijk++)
		fprintf(stdout,"%.5e ", sc->x[kf*dim+ijk]);
	      //	iindex[kf]=-abs(iindex_parent[kfp]);
	      fprintf(stdout,"]\n");
	    } else {
	      //	iindex[kf]=nv+(*nvall0);
	      nv++;
	    }
	  }
	  // now we have numv vertices which we want to match to faces. 
	  /* print_full_mat_int(1,nk,mip,"mip");  */
	  /* print_full_mat_int(1,nk,tip,"tip");  */
	  /* print_full_mat_int(1,nj,mi,"mi");  */
	  /* print_full_mat_int(1,nj,ti,"ti");  */
	} else {
	  fprintf(stderr,"\n*** An error in counting overlaps in %s ***\n",__FUNCTION__);
	}
      }
    }
  }
  //  fprintf(stdout,"\nface(%d) in a divided el(%d) is also face(%d) in the current el(%d)",je,kel,ke,jel);
  free(m);
  free(mp);
  free(mi);
  free(mip);
  free(vertk);
  free(vertj);
  //  free(ti); //same as vertk
  //  free(tip); //same as vertj
  return;
}
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
  //  fprintf(stdout,"\nnsall=%d,nvall=%d",nsall,nvall);fflush(stdout);
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
  //  fprintf(stdout,"\nsc->nv=%d,sc->ns=%d; nvall=%d,nsall=%d\n",sc->nv,sc->ns,nvall,nsall);fflush(stdout);
  return;
}
/*ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ*/
scomplex *macro_split(input_grid *g0,cube2simp *c2s)
{
  /* 
     From an input grid loops over the macroelements and generates
     meshes depending on the division given in every dimension. 
  */
  input_grid *g;
  INT i,j,j0,j1,kel,ke,pmem;
  INT nvcube=c2s->nvcube,nvface=c2s->nvface;
  scomplex **sc=malloc(g0->nel*sizeof(scomplex *));
  pmem=2*g0->nv;
  if(pmem<2*g0->ne) pmem=2*g0->ne;
  if(pmem<2*g0->nel) pmem=2*g0->nel;
  if(pmem<2*g0->nf) pmem=2*g0->nf;
  INT *p=calloc(pmem,sizeof(INT));
  //  print_full_mat_int(g0->nel,(c2s->nvcube+1),g0->mnodes,"ex");
  ilexsort(g0->nf, (c2s->nvface+1),g0->mfaces,p);
  ilexsort(g0->nel,(c2s->nvcube+1),g0->mnodes,p);
  /*-------------------------------------------------------------------*/
  //  INT *efound=calloc(c2s->ne*(g0->nel),sizeof(INT));
  INT **nd=calloc(g0->nel,sizeof(INT *));
  for(i=0;i<g0->nel;i++){
    nd[i]=calloc(c2s->n,sizeof(INT)); /* to hold the number of
 					 divisions in every coordinate
					 direction */
  }
  for(kel=0;kel<g0->nel;kel++)
    for(i=0;i<c2s->nf;i++){
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
  //print_full_mat_int(g0->ne,3,g0->seg,"newseg");
  // get the macroelement mesh in a structure
  macrocomplex *mc=set_mmesh(g0,c2s,p);
  /*PLACE HOLDERS*/
  INT **elneib=mc->elneib;
  INT **el2fnum=mc->el2fnum;
  INT *isbface=mc->isbface;
  INT *bcodesf=mc->bcodesf;
  INT *etree=mc->etree;
  iCSRmat *bfs0=mc->bfs;
  iCSRmat *fullel2el=mc->fullel2el;
  fprintf(stdout,"\n%%DFS(domains): %d connected components",mc->cc);
  fprintf(stdout,"\n%%DFS(boundaries): %d connected components",mc->bndry_cc);
  /* fprintf(stdout,"\nbfs0=["); */
  /* icsr_print_matlab_val(stdout,bfs0); */
  /* fprintf(stdout,"];"); */
  /* fprintf(stdout,"\nbfs=sparse(bfs0(:,1),bfs0(:,2),bfs0(:,3));\n"); */
  /* print_full_mat_int(1,mc->nf,mc->isbface,"isbface"); */
  /* print_full_mat_int(1,mc->nf,mc->bcodesf,"bcodesf"); */
  /* fprintf(stdout,"\n*****   nf=%d (%d)******* \n",g0->nf,mc->nf); */
  /*full el2el*/
  fprintf(stdout,"\nfel2el0=[");
  icsr_print_matlab_val(stdout,mc->fullel2el);
  fprintf(stdout,"];");
  fprintf(stdout,"\nfel2el=sparse(fel2el0(:,1),fel2el0(:,2),fel2el0(:,3));\n");
  /* print_full_mat_int(1,mc->nf,mc->isbface,"isbface"); */
  /* print_full_mat_int(1,mc->nf,mc->bcodesf,"bcodesf"); */
  fprintf(stdout,"\n*****   nf=%d (%d)******* \n",g0->nf,mc->nf);
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
  INT lvl,intype=-1,kj,je,jel;  
  INT *codef=calloc(c2s->nf,sizeof(INT));
  INT *isbndf=calloc(c2s->nf,sizeof(INT));
  // here this can be also done without bfs
  for(kj=0;kj<bfs0->nnz;kj++){    
    jel=bfs0->JA[kj];
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
		  isbndf,codef,mc->flags[jel],				\
		  intype);
    nsall+=sc[jel]->ns;
    nvall+=sc[jel]->nv;
    //      fprintf(stdout,"\n%%Mapping back to the macroelement...\n");
    map2mac(sc[jel],c2s,g);
    fprintf(stdout,"\n%%mesh(macroelement=%d): nv=%d; nsimp=%d",jel,sc[jel]->nv,sc[jel]->ns);      
    haz_scomplex_print(sc[jel],0,"HAHA");
  }
  fix_grid(mc,			\
	   sc,			\
	   c2s,			\
	   g0);
  scomplex_merge(sc,			\
		 nsall, nvall,			\
		 mc->cc, mc->bndry_cc,		\
		 g0,c2s);
  fprintf(stdout,"\n%%");
  fprintf(stdout,"\n%%merged(macroelements=%d:%d): nv=%d; nsimp=%d",	\
	  0,g0->nel-1,sc[0]->nv,sc[0]->ns);  
  fprintf(stdout,"\n%%\n");
  //  haz_scomplex_print(sc[0],0,"MERGED");
  //  for(i=0;i<g0->nel;i++){
  //    
  //  }
  free(p);
  free(isbndf);
  free(codef);
  macrocomplex_free(mc);
  return sc[0];  
}
/****************************************************************/
INT main(INT argc, char **argv)
{
  //  INT i=-1;
  input_grid *g=parse_input_grid("grid.input");
  INT dim=g->dim;
  cube2simp *c2s=cube2simplex(dim);
  scomplex *sc=macro_split(g,c2s);
  if(dim==2||dim==3) {
    fprintf(stdout,"Writing vtk file...\n");
    vtkw("newmesh.vtu",sc,0,0,1.);
  }
  /*FREE*/
  haz_scomplex_free(sc);
  input_grid_free(g); 
  cube2simp_free(c2s);
  fprintf(stdout,"\nDone.\n");
  return 0;
}
/*********************EOF**********************************/
