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
void get_faces(scomplex *sc)
{
  INT i,j,k,nfi,nfb,nf,nbr,dim=sc->n,nv=sc->nv,ns=sc->ns;
  INT dim1=dim+1;
  // for a constructed simplicial complex, finds face to vertex map,
  // dfs ordering of the elements, connected components, connected
  // components on the boundary.
  sc->f2s=malloc(1*sizeof(iCSRmat));
  // nfi interior faces, nfb boundary faces. 
  nfi=0;
  nfb=0;
  nf=0;
  for(k=0;k<dim1;k++){
    for(j=0;j<sc->ns;j++){
      nbr=sc->nbr[j*dim1+k];
      if(nbr>j) {
	nf++;
	nfi++;
	fprintf(stdout,"\n: j=%d, nbr=%d",j,nbr);
      } else if(nbr<0) {
	nfb++;
	nf++;
	fprintf(stdout,"\n: j=%d, nbr=%d",j,nbr);
      } 
    }
  }
  fprintf(stdout,"\n: interior faces=%d, boundary faces=%d,all faces=%d",nfi,nfb,nf);
  sc->f2s->IA=icsr_create_p(nf,ns,2*nf)
  // we now get the interior faces and then the boundary faces
  return; 
}
SHORT area_face0(INT dim, REAL factorial, REAL *xf, REAL *sn,	\
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
  //  print_full_mat(dim,1,piv,"piv");
  if(lufull(1, dim, volt, bt,p,piv)) {
    //print_full_mat(dim,dim,bt,"bt");
    volt[0]=0.;
    return 1; // this is a degenerate simplex
  } else
    volt[0]=fabs(volt[0])/factorial;
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
  return 0;
}
/*********************************************************************/
// an interior vertex is defined as the vertex such that all faces attached to it are interior.
//
INT main(INT   argc,   char *argv[])
{
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
  if(amr_marking_type){
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
       * SOLVE on the finest (for now grid)
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
  } else {
    // refine ref_levels;
    refine(ref_levels,sc,NULL);
  }
  /*  MAKE sc to be the finest grid only */
  scfinalize(sc);
  get_faces(sc);
  /* write the output mesh file:    */
  hazw(g->fgrid,sc,0);
  /* WRITE THE OUTPUT vtu file for paraview:    */
  vtkw(g->fvtu,sc,0,1.);
  /*FREE*/
  input_grid_free(g);
  free(all);
  haz_scomplex_free(sc);
  return 0;
}
