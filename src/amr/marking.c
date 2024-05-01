/*! \file src/amr/marking.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20240430.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note contains functions for marking strategies used in AMR.
 *
 */
#include "hazmath.h"

/*!
* \fn ivector amar_maximal_mark(scomplex *sc,dvector *estimator, REAL gamma)
*
* \brief mark elements for adaptive refinement using maximal marking strategy
*           mark all elements with error estimator greater than gamma*max(estimator)
*
* \param sc 	       Pointer to the simplicial complex
* \param estimator   Pointer to error estimator
* \param gamma       Fraction of max error estimator to mark (eg. 0 -> mark all, 1 -> mark none)
*
* \return marked     ivector of elements to be marked, TRUE for marked
*
*/
ivector amr_maximal_mark(scomplex *sc, REAL *estimator, REAL gamma) 
{
  
  INT i,k=0;
  INT ns = sc->ns;
  REAL errmax;
  ivector marked=ivec_create(sc->ns);
  errmax=estimator[0];
  for(i=1;i<ns;i++) {
    if(fabs(estimator[i])>errmax) {
      errmax=estimator[i];
    }
  }
  for(i=0;i<ns;i++) {
    if(fabs(estimator[i])>gamma*errmax){
      marked.val[i]=TRUE;
      k++;
    } else {
      marked.val[i]=FALSE;
    }
  }
  
  printf("\nMarking elements for AMR using Maximal Mark strategy.\n\tMarking elements with %3.2f%% of the max error: --> %d of %d elements were marked to be refined.\n",gamma*100,k,ns);
  return marked;
}
