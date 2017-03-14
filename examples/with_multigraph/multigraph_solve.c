/*! \file examples/ConvectionDiffusion/multigraph_solve.c
 *
 *  Created by Xiaozhe Hu, James Adler, and Ludmil Zikatanov 20170309
 *
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This is a wrapper to call the
 *      multigraph package (see http://ccom.ucsd.edu/~reb/software.html)
 *      by R.E.Bank for the solution of general linear system Ax = b
 *
 *      
 *
 *        
 */
#if defined (__cplusplus)
extern "C" {
#endif

#include "hazmath.h"
#include "multigraph_solve.h"

  //  void mgraph_wrap(dCSRmat A, dvector f, dvector *u)
  void mgraph_wrap(INT doilu, INT nrow, INT *ia, INT *ja, REAL *a, REAL *rhs, REAL *sol, INT *jareb, REAL *areb, INT *ka)
  {
    /* all indices start from 0 here so we do too many shifts: TODO
       less */
    INT n=nrow,n1=n+1;
    INT i,nnzlu=-16;
    //    INT  *jareb=NULL;
    //    REAL *areb=NULL;
    //    FILE *fp;
    csrreb(&n,&n, &nnzlu,ia,ja,a,&jareb,&areb);
    INT maxja = MAXJA_MULT*(n1 + jareb[n]-jareb[0]);
    if(maxja < 32*256*256*256) maxja=32*256*256*256;
    INT maxa = MAXA_MULT*maxja;
    if(doilu){
      fprintf(stdout,"\n\nmaxja=%d maxa=%d\n",maxja,maxa);
      fflush(stdout);
      jareb = (INT *)realloc(jareb,maxja*sizeof(INT));
      areb = (REAL *)realloc(areb,maxa*sizeof(REAL));
      INT  nblock = 1;
      INT iblk[2]={0,0};
      iblk[nblock]=n1;
      iblk[nblock-1]=1;
      ka = (INT *) calloc(10*(maxlvl+1),sizeof(INT));
      //shift and do mginit
      for (i=0;i<n1+nnzlu;++i) jareb[i]+=1;
      mginit_(&n, &ispd, &nblock, iblk,			\
	    &maxja, jareb, &maxa, areb,			\
	      &ncfact, &maxlvl, &maxfil, ka,		\
	      &lvl, &dtol, &method,			\
	      &iflag );
    }
    REAL hist[LENGTH_HIST];
    mg_(&ispd, &lvl, &mxcg, &eps1, jareb, areb,	\
	sol, rhs, ka, &relerr, &iflag, hist );
    /* unshift */
    for (i=0;i<n1+nnzlu;++i) jareb[i]-=1;
    /* */
    fprintf(stdout,"\n Rel. Err.=%12.5e\n\n",relerr);   
    //ka is the array which contains all info. 
    INT iend=ITNUM_FROM_HIST, itnum= (INT )hist[iend];
    if(itnum<iend) iend=itnum;
    fprintf(stdout,"\nMultigraph history (iter_num=%5i)\n",itnum);    
    for (i=0;i<iend;++i){    
      fprintf(stdout,"iter =  %5i; res %12.4e\n",itnum-iend+i+1,hist[i]);    
    }
    return;
  }

#if defined (__cplusplus)
}
#endif

