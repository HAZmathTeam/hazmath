#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "hazmath.h"
#include "grid_defs.h"
#include "grid_params.h"
/*****************************************************************/
/*****************************************************************/
void mock_solve(INT ref_type, scomplex *sc, void *all){
  if(!ref_type) return;
  fprintf(stdout,"\n%s: %d",__FUNCTION__,sc->level);
  //  INT j,k;
  //  INT dim=sc->n,nv=sc->nv,ns=sc->ns;
  //  dvector *u=all;
  /* 
     if(u)
     for(j=0;j<ns;j++)
     fprintf(stdout,"\nu[%d]: %f",j,u[j]);    
  */
  return;
}
void mock_estimate(INT ref_type, scomplex *sc, void *all){
  INT j,lvl_prev=(sc->level-1);
  REAL gamma=5e-1;
  dvector *errs=all;
  if(ref_type<2)return;
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
      fprintf(stdout,"\n: err[%d]=%e;",		\
	      j,errs->val[j]);    
    }
  } else {
    /*    assuming that every refined simplex reduces the error by a
	  factor of gamma */
    for(j=0;j<sc->ns;j++){
      if(sc->child0[j]<0 || sc->childn[j]<0 || sc->gen[j] != lvl_prev)continue;
      errs->val[sc->childn[j]]=errs->val[sc->child0[j]]=gamma*errs->val[j];
      fprintf(stdout,"\ngeneration: %d;j=%d,e=%e (child0[%d]=%e; childn[%d]=%e", \
      	      sc->gen[j],					\
      	      j,errs->val[j],					\
      	      sc->child0[j],errs->val[sc->child0[j]],		\
      	      sc->childn[j],errs->val[sc->childn[j]]);
    }
    for(j=0;j<sc->ns;j++)
      fprintf(stdout,"\nMMMerrs[%d]: %e",j,errs->val[j]);
  }
  fprintf(stdout,"\n%s: %d",__FUNCTION__,sc->level);
  return;
}
void mock_mark(INT ref_type, scomplex *sc,void *all){
  dvector *errs=all;
  //  markeql(sc, errs);
  markall(sc,1);
  //  markstar(sc,errs);
    //  print_full_mat_int(sc->ns,1,sc->marked,"mark");
    //  fprintf(stdout,"\nXXXXXXXXXX%s: %d\n",__FUNCTION__,sc->level);
  return;
}
/*****************************************************************************/
INT main(INT   argc,   char *argv[])
{
  INT use_features=0;
  FILE *fp=stdin;     
  //  fp=HAZ_fopen("3d_cube.input","r"); 
  input_grid *g=parse_input_grid(fp);
  fclose(fp);
  //  input_grid_print(g);
  scomplex *sc=generate_grid(g);
  input_grid_free(g);
  INT ref_levels=g->nref,dim=g->dim,dim1=g->dim+1;
  if(ref_levels){
    /* lots=number of elements if every element was refined ref_level times*/
    long int lots = (1<<ref_levels)*sc->ns;
    //fprintf(stdout,"\n*** NNN=%li %d\n",lots,(1<<ref_levels));
    /* create a vector to contain the errors of all refinements*/
    dvector *errest=(dvector *)malloc(sizeof(dvector));
    *errest=dvec_create((INT )lots);
    /*******************************************/
    INT amr_ref_type=2; /*default is refine equilibrium*/
    if(use_features){
    /* if we have to refine only elements containing points from an
       array called "features array" */
      features *feat=malloc(sizeof(features));
      feat->nbig=dim;
      feat->n=dim;
      REAL vvv=0.;
      features_r(dim,use_features,feat,vvv);
      errest->row=feat->nf;
      errest->val=feat->x;
      amr_ref_type=1;
    } 
    n_refine(amr_ref_type,ref_levels,sc,errest,	\
	     &mock_solve,			\
	     &mock_estimate,			\
	     &mock_mark);  
    scfinalize(sc);
  }
  // write the output mesh file:    
  //  hazw(nameout,sc,0,0);
  fprintf(stdout,"\n\n%%Writing a vtk file...\n");
  vtkw("newamr.vtu",sc,0,0,1.);
  /*FREE*/
  haz_scomplex_free(sc);
  return 0;
  free(sc);
  return 0;
}
