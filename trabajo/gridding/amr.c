#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "hazmath.h"
void markall(scomplex *sc,const int amark);
void markstar(scomplex *sc, dvector *w);
unsigned int markeql(scomplex *sc, dvector *w);

/* refinement */
void n_refine(INT ref_type, INT ref_levels, scomplex *sc,	\
	      dvector *errors,					\
	      void (*solving)(INT, scomplex *, void *),		\
	      void (*estimating)(INT , scomplex *, void *),	\
	      void (*marking)(INT , scomplex *, void *));
/******************************************************************/
/*!
 * \fn INT dvec_set_amr(const REAL value, scomplex *sc, dvector *pts,
 * dvector *toset)
 *
 * \brief given a dvector of size sc->ns sets the values of a vector
 *        to equal value at every simplex that is on the last level of
 *        refinement (not refined simplex) and contains a point from
 *        dvector pts (note that the size of pts->val should be
 *        sc->n*pts->row)
 *
 * \param dvector toset; toset->val is (re)-allocated if needed
 *
 * \return number of simplices where the value was assigned
 *
 */
INT dvec_set_amr(const REAL value, scomplex *sc, dvector *pts, dvector *toset)
{
  INT j,k,jpts,n=sc->n,n1=sc->n+1,ns=sc->ns;
  REAL *pval0=NULL; /*place holder*/
  INT *scnjn=NULL; /*place holder*/
  k=0;
  toset->row=sc->ns;
  toset->val=realloc(toset->val,toset->row*sizeof(REAL));
  for(j=0;j<ns;j++) {
    if((sc->gen[j]!=sc->level) && (sc->child0[j]>=0 && sc->childn[j]>=0)){
      /* this was refined*/
      continue;
    }
    scnjn = sc->nodes+j*n1; /* beginning of local vertex numbering for
			       simplex j.*/
    for(jpts=0;jpts<pts->row;jpts++){
      pval0=pts->val+jpts*n; 
      if(!xins(n,scnjn,sc->x,pval0)){
	fprintf(stdout,"\nel=%d, found: %d",j,jpts);
	toset->val[j]=value;
	k++;
	break;
      }
    }
  }
  return k;
}
/*****************************************************************/
void exmpl_solve(INT amr_ref_type, scomplex *sc, void *all){
  INT dim=sc->n,ns=sc->ns,i;
  dvector *exact,*solfem,*pts;// place holders for data in all.
  INT npts=1;
  /*solve phase; as output we have a fem solution femsol (value on
    every simplex) */  
  switch(amr_ref_type){
  case 1: case 2: case 3:
    /* 
       case 1 and we can do 2 and 3 later; 
       refine around a given set of points or a point. use the void
       pointer to carry data around.
    */
    if(!sc->level){
      // if no refinement was done yet, initialize
      all=(void *)malloc(3*sizeof(dvector));
      pts=(dvector *)all;
    /* ALT for this below is: 
       solfem=(dvector *)(all +sizeof(dvector)); 
    */
      solfem=pts+1;
      pts[0]=dvec_create(dim*npts);
      memset(pts->val,0,pts->row*sizeof(REAL));// always near the origin
    /*    
	  in this example scenario: exact solution is zero and the fem
	  solution equals 1 in every element that contains the points in
	  pts;
    */
      solfem->row=0; solfem->val=NULL; //allocated when "solving"
    }
    dvec_set_amr(1.,sc,pts,solfem);
    //exact solution is 0 in this example. 
    break;
  default:
    /* 
       Do nothing, amr_ref_type corresponds to a uniform
       refinement
    */
    break;
  }
  return;
}
void exmpl_estimate(INT amr_ref_type, scomplex *sc, void *all){
  INT j,lvl_prev=(sc->level-1),dim=sc->n;
  REAL gamma=5e-1;// reduction factor. 
  dvector *errs=all;
  switch(amr_ref_type){
  case 1:
    /*
      we estimate here by setting the errors equal to 1 at every
      simplex that contains points from pts; 
    */
    if(!sc->level){
      /*on the coarsest grid generate random errors */
      INT ntotal=errs->row;
      dvec_set(ntotal,errs,0e0);
      dvec_rand_true(sc->ns,errs);
      memset(sc->marked,0x0,sizeof(INT)*sc->ns);
      errs->row=ntotal;
      /* in this mock scenario we just have random error on 0th level */
      for(j=0;j<sc->ns;j++){
	errs->val[j]=fabs(errs->val[j]);
	fprintf(stdout,"\n: err[%d]=%e;",	\
		j,errs->val[j]);    
      }
    } else {
      /*    
	    assuming that every refined simplex reduces the error by a
	    factor of gamma
      */
      for(j=0;j<sc->ns;j++){
	if(sc->child0[j]<0 || sc->childn[j]<0 || sc->gen[j] != lvl_prev)continue;
	errs->val[sc->childn[j]]=errs->val[sc->child0[j]]=gamma*errs->val[j];
	fprintf(stdout,"\ngeneration: %d;j=%d,e=%e (child0[%d]=%e; childn[%d]=%e", \
		sc->gen[j],						\
		j,errs->val[j],						\
		sc->child0[j],errs->val[sc->child0[j]],			\
		sc->childn[j],errs->val[sc->childn[j]]);
      }
      for(j=0;j<sc->ns;j++)
	fprintf(stdout,"\nMMMerrs[%d]: %e",j,errs->val[j]);
    }
    break;
  default:
    // do nothing this is a regular refinement;
    break;
  }
  fprintf(stdout,"\n%s: %d",__FUNCTION__,sc->level);
  return;
}
void exmpl_mark(INT amr_ref_type, scomplex *sc,void *all){
  dvector *errs=(dvector *)all;
  //  markeql(sc, errs);
  switch(amr_ref_type){
  case 1: case 2:
    //mareql(sc,all);
    break;
  default:
    markall(sc,1);
    break;
  }
  //  markstar(sc,errs);
  //  print_full_mat_int(sc->ns,1,sc->marked,"mark");
  //  fprintf(stdout,"\nXXXXXXXXXX%s: %d\n",__FUNCTION__,sc->level);
  return;
}
/*****************************************************************************/
INT main(INT   argc,   char *argv[])
{
  FILE *fp=stdin;     
  //  fp=HAZ_fopen("3d_cube.input","r"); 
  input_grid *g=parse_input_grid(fp);
  fclose(fp);
  //  input_grid_print(g);
  scomplex *sc=generate_grid(g);
  INT ref_levels=g->nref,dim=g->dim,dim1=g->dim+1;
  input_grid_free(g);
  if(ref_levels){
    /* lots=number of elements if every element was refined ref_level times*/
    long int lots = (1<<ref_levels)*sc->ns;
    //fprintf(stdout,"\n*** NNN=%li %d\n",lots,(1<<ref_levels));
    /* create a vector to contain the errors of all refinements*/
    dvector *errest=(dvector *)malloc(sizeof(dvector));
    *errest=dvec_create((INT )lots);
    /*******************************************/
    INT amr_ref_type=0;
    n_refine(amr_ref_type,ref_levels,sc,errest,	\
	     &exmpl_solve,			\
	     &exmpl_estimate,			\
	     &exmpl_mark);  
  }
  // write the output mesh file:    
  //  hazw(nameout,sc,0,0);
  fprintf(stdout,"\n\n%%Writing a vtk file...\n");
  vtkw("newamr.vtu",sc,0,0,1.);
  /*FREE*/
  haz_scomplex_free(sc);
  return 0;
}
