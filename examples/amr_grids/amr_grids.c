#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "hazmath.h"
#include "solve_estimate_mark.h"
INT main(INT   argc,   char *argv[])
{
  FILE *fp=stdin;     
  //  fp=HAZ_fopen("2d_square.input","r"); 
  input_grid *g=parse_input_grid(fp);
  fclose(fp);
  //    input_grid_print(g);
  //    exit(44);
  scomplex *sc=generate_initial_grid(g);
  scomplex *sctop=NULL;
  INT ref_levels=g->nref, amr_marking_type=g->mark_type,j; 
  dvector *solfem=NULL,*estimator=NULL;
  ivector *marked=NULL;
  void *all=NULL;
  if(amr_marking_type){
    // all here can pass data around. so we make it 4 dvectors worth and one ivector (as an example)
    all=(void *)malloc(5*sizeof(dvector)+sizeof(ivector)); 
    /**/    
    for(j=0;j<ref_levels;j++){
      /*
       * select the finest level: 
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
       *       refinement
       */
      marked=exmpl_mark(sctop,estimator,all);
      /* 
       *  Refine the grid. this always refines 1 time;
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
  //  shrink the sc to hold the finest grid only
  scfinalize(sc);
  // write the output mesh file:    
  //  hazw(g->fgrid,sc,0,0);
  fprintf(stdout,"\n\n%%Writing a vtk on file...%s\n",g->fvtu);
  vtkw(g->fvtu,sc,0,0,1.);
  /*FREE*/
  input_grid_free(g);
  free(all);
  haz_scomplex_free(sc);
  return 0;
}
