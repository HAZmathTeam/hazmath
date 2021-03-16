/*! \file examples/basic_elliptic/basic_elliptic.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 2019/01/09.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program generates simplicial grids in 2,3,4... dimension.
 *
 * \note This example highlights some of the features of the simple
 * mesh generator included with HAZmath. It is only to illustrate how
 * to use the mesh refinement.
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
/*********************************************************************/
#include "solve_estimate_mark.h"
/*********************************************************************/
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
  fprintf(stdout,"\nelnodes=[");
  for(i=0;i<nbf;++i){
    //    fprintf(stdout,"\nelnodes[%d]=(",i);
    for(j=0;j<dim;++j){
      fprintf(stdout,"%4i ",fnodes[dim*i+j]+1);
    }
    fprintf(stdout,";\n");
  }
  fprintf(stdout,"]\n");
  find_nbr(nbf,nv,(dim-1),fnodes,fnbr);
  fprintf(stdout,"\nelnbr=[");
  for(i=0;i<nbf;++i){
    //    fprintf(stdout,"\nelnbr[%d]=(",i);
    for(j=0;j<dim;++j){
      fprintf(stdout,"%4i ",fnbr[dim*i+j]+1);
    }
    fprintf(stdout,";\n");
  }
  fprintf(stdout,"]\n");
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
// an interior vertex is defined as the vertex such that all faces attached to it are interior.
//
INT main(INT   argc,   char *argv[])
{
  INT i;
  FILE *fp=stdin;
  //  fp=HAZ_fopen("2d_square.input","r");
  /*
    PARSE THE INPUT.
  */
  input_grid *g=parse_input_grid(fp);
  fclose(fp);
  /*
    GENERATE INITIAL GRID AND DECLARE VARIABLES.
  */
  scomplex *sc=generate_initial_grid(g);
  scomplex *sctop=NULL;
  INT ref_levels=g->nref, amr_marking_type=g->mark_type,j;
  dvector *solfem=NULL,*estimator=NULL;
  ivector *marked=NULL;
  void *all=NULL;
  REAL *xstar=NULL;
  INT nstar,dim=sc->n;
  if(amr_marking_type==0){
    // refine ref_levels;
    refine(ref_levels,sc,NULL);
  } else if(amr_marking_type==33){
    REAL h = 1.0/128;  // step distance of points
    REAL threshold = h; // threshold for close to the points or not
    INT nstep = 0;
    nstar= 2 + nstep*4; // refining near several points: (even number )
    xstar=(REAL *)calloc(nstar*dim,sizeof(REAL));
    xstar[0*dim+0]=1.666667e-1;
    xstar[0*dim+1]=6.666667e-1;
    xstar[1*dim+0]=8.333333e-1;
    xstar[1*dim+1]=6.666667e-1;
    for (i=1; i<nstep; i++){
      // points near (1/6, 2/3)
      xstar[(2+(i-1)*4)*dim+0]=1.666667e-1 + i*h;
      xstar[(2+(i-1)*4)*dim+1]=6.666667e-1;
      xstar[(3+(i-1)*4)*dim+0]=1.666667e-1 - i*h;
      xstar[(3+(i-1)*4)*dim+1]=6.666667e-1;

      // points near (5/6, 2/3)
      xstar[(4+(i-1)*4)*dim+0]=8.333333e-1 + i*h;
      xstar[(4+(i-1)*4)*dim+1]=6.666667e-1;
      xstar[(5+(i-1)*4)*dim+0]=8.333333e-1 - i*h;
      xstar[(5+(i-1)*4)*dim+1]=6.666667e-1;
    }
    for(j=0;j<ref_levels;j++){
      /*
       * SELECT the finest grid:
       */
      sctop=scfinest(sc);
      /* MARK: marked is an ivector with num.rows=the number of
       *       simplices; its componenets are nonzero if the simplex
       *       is marked for refinement and 0 if it is not marked for
       *       refinement.
       */
      marked=mark_near_points(sctop,nstar,xstar, threshold);
      refine(1,sc,marked);
      /* free */
      haz_scomplex_free(sctop);
    }
    ivec_free(marked);
    free(xstar);
  } else {
    /*
       Use "all" here can pass data around. Below we make 4 dvectors
       and one ivector, just as example. A good example for using the
       array will be to pass the whole hierarchy, not only the fine
       grid via all, e.g. all=(coid *)sc
    */
    all=(void *)malloc(5*sizeof(dvector)+sizeof(ivector));
    /**/
    for(j=0;j<ref_levels;j++){
      /*
       * SELECT the finest grid:
       */
      sctop=scfinest(sc);
      /*
       * SOLVE on the finest (for now) grid
       */
      solfem=(dvector *)exmpl_solve(sctop,all);
      /*
       * ESTIMATE using the numerical solution and the data stored in *all
       */
      estimator=(dvector *)exmpl_estimate(sctop,solfem,all);
      /*
       * MARK: marked is an ivector with num.rows=the number of
       *       simplices; its componenets are nonzero if the simplex
       *       is marked for refinement and 0 if it is not marked for
       *       refinement.
       */
      marked=exmpl_mark(sctop,estimator,all);
      /*
       *  Refine the grid. this always refines 1 time, but since we
       *  are in a loop, it will refine ref_levels times;
       */
      refine(1,sc,marked);
      /* free */
      haz_scomplex_free(sctop);
      dvec_free(solfem);
      dvec_free(estimator);
      ivec_free(marked);
    }
    free(all);
  }
  /*  MAKE sc to be the finest grid only */
  scfinalize(sc);
  find_cc_bndry_cc(sc);
  /* write the output mesh file:    */
  hazw(g->fgrid,sc,0);
  /* WRITE THE OUTPUT vtu file for paraview:    */
  vtkw(g->fvtu,sc,0,1.);
  /*FREE*/
  input_grid_free(g);
  haz_scomplex_free(sc);
  return 0;
}
