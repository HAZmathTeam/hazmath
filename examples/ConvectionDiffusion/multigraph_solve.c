/*! \file examples/ConvectionDiffusion/multigraph_solve.c
 *
 *  Created by Xiaozhe Hu, James Adler, and Ludmil Zikatanov 20170309
 *
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This is a wrapper to call the
 *      multigraph package (see http://ccom.ucsd.edu/~reb/software.html)
 *      by R.E.Bank for the solution of 
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
  
  void mgraph_wrap(dCSRmat A, REAL *rhs, REAL *sol)
  {
    INT n=A.row,*ia=A.IA,*ja=A.JA,n1=n+1;
    REAL *a=A.val;
    INT  *jareb=NULL;
    REAL *areb=NULL;
    INT i,nnzlu=-16;
    FILE *fp;
    csrreb(&n,&n, &nnzlu,ia,ja,a,&jareb,&areb);
    // guess for the storage at all levels.
    // triple the memory, because this is to say
    // that the storage we need for all levels triples. 
    INT maxja = MAXJA_MULT*(n1 + jareb[n]-jareb[0]);
    // this is because it does not work for a diagonal matrix
    if(maxja < 7*n) maxja=7*n;
    INT maxa = MAXA_MULT*maxja;
    jareb = (INT *)realloc(jareb,maxja*sizeof(INT));
    areb = (REAL *)realloc(areb,maxa*sizeof(REAL));
    INT  nblock = 1;
    INT iblk[2]={0,0};
    iblk[nblock]=n1;
    iblk[nblock-1]=1;
    INT iflag=-16;  
    /// outputs
    INT lvl=-16;
    INT *ka = (INT *) calloc(10*(maxlvl+1),sizeof(INT));
    /*  
	INT lenz = 2*(2*maxja+5*n1);
	REAL *z = (REAL *)calloc(lenz,sizeof(double));
    */
    /*  fprintf(stdout,"\n************ nnzlu+n1 =  %i; other = %i\n\n", 
	nnzlu+n1,n1+jareb[n]-jareb[0]);
    */
    ///
    //shift and do mginit
    for (i=0;i<n1+nnzlu;++i) jareb[i]+=1;
    mginit_(&n, &ispd, &nblock, iblk,			\
	    &maxja, jareb, &maxa, areb,			\
	    &ncfact, &maxlvl, &maxfil, ka,
	    &lvl, &dtol, &method, 
	    &iflag );
    INT j, ij;
    REAL eps1=1e-8;
    INT mxcg=1000;
    INT ns=ka[10*(lvl+1-1)+(2-1)]-1;
    REAL hist[LENGTH_HIST];
    REAL relerr=1e0;  
    mg_(&ispd, &lvl, &mxcg, &eps1, jareb, areb,	\
	sol, rhs, ka, &relerr, &iflag, hist );
    /* unshift */
    for (i=0;i<n1+nnzlu;++i) jareb[i]-=1;
    if(exsol){
      REAL err0=fabs(sol[0]-exsol[0]);
      for(j=1; j<n;++j){
	if(err0<fabs(sol[j]-exsol[j])){
	  err0=fabs(sol[j]-exsol[j]);
	}
      }
      fprintf(stdout,"\n Rel. Err.=%12.5e, err_infty=%12.5e\n\n",relerr, err0);
    }else{
      fprintf(stdout,"\n Rel. Err.=%12.5e\n\n",relerr);
    }
    //ka is the array which contains all info. 
    INT iend=20, itnum= (INT )hist[iend];
    if(itnum<iend) iend=itnum;
    fprintf(stdout,"\nMultigraph history (iter_num=%5i)\n",itnum);    
    for (i=0;i<iend;++i){    
      fprintf(stdout,"iter =  %5i; res %12.4e\n",itnum-iend+i+1,hist[i]);    
    }
    //    fprintf(stdout,"iflag=%i ; lvl= %i\n",iflag,lvl);
    //  for (i=0;i<n1+nnzlu;++i)jareb[i]-=1;
    ///  fprintf(stdout,"\n");
    if(ka) free(ka);
    //  if(z) free(z);
    if(jareb) free(jareb); /* these are all integers */
    if(areb) free(areb); /*this should be all  reals */
    return;
  }

#if defined (__cplusplus)
}
#endif

