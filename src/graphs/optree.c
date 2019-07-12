/*! \file src/graphs/optree.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190115.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note Kruskal algorithm;
 *
 */
/*  Kruskal algorithm for min weighted tree input is the signed
      edge<->vertex incidence matrix and a permutation vector p which
      shows in what order the edges have to be examined in the Kruskal
      algorithm */

#include "hazmath.h"

INT optree(INT *gradi, INT *gradj, INT *p, INT ne, INT n, INT *et, INT *net)
{
  INT ket=-1, kea=-1, nm1=0, ie=0,je=0,iee=0,jee=0;
  INT i,nem1,nm2;
  INT *pa=NULL; /* working */
  nm1=n-1;
  nm2=nm1-1;
  nem1=ne-1;
  pa = (INT *)calloc(n,sizeof(INT));
  /* init pa=parent[i] -- each poINT is a root of a tree at the beginning. */
  for (i=0;i<n;i++)  pa[i]=i;
  while(kea< nem1 && ket < nm2) {
    kea++;
    ie=gradi[p[kea]];
    je=gradj[p[kea]];
    iee = ie;
    jee = je;
    while(pa[iee] != iee)  iee=pa[iee];
    while(pa[jee] !=jee) jee=pa[jee]; // fprintf(stderr,"%i\n",iee); }
    /* the possibilities now are that iee ~=jee or iee=jee */
    if(iee < jee) {
      pa[jee]=iee;
    }
    else if(iee > jee) {
      pa[iee]=jee;
    }
    else
      continue;
    ket++;
    et[ket]=p[kea];    
  }
  *net=ket+1;
  /*  fprintf(stdout,"\nNodes and Trees\n");
  for (i=0;i<n;++i)
    {
      fprintf(stdout,"Node=%i; Tree=%i\n",i,pa[i]);
    }
  */
  if(pa) free(pa);
  return 0; 
}

