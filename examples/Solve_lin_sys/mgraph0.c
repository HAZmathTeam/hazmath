#if defined (__cplusplus)
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#define INT int
#define REAL double

void  mginit_(INT *n, INT *ispd, INT *nblock, INT *iblock,		\
	      INT *maxja, INT *jareb, INT *maxa, REAL *areb,		\
	      INT *ncfact, INT *maxlvl,	INT *maxfil, INT *ka,
	      INT *lvl, REAL *dtol, INT *method, 
	      INT *iflag );
  /*        subroutine mginit(n,ispd,nblock,ib,maxja,ja,maxa,a,ncfact,
	    +      maxlvl,maxfil,ka,lvl,dtol,method,iflag)
  */
  ///  call mg(ispd,lvl,mxcg,eps1,jareb,areb,z,rhs,
  ///	  >     ka,relerr,iflag,z1,hist)
void mg_(INT *ispd, INT *lvl, INT *mxcg, REAL *eps1,			\
	 INT *jareb, REAL *areb,		\
	 REAL *sol, REAL *rhs, INT *ka, REAL *relerr,		\
	 INT *iflag, REAL *hist );
  /*   
       subroutine mg(ispd,lvl,mxcg,eps1,ja,a,dr,br,ka,relerr,
       +      iflag,hist)
  */
  // io routines for debugging. 


void mgraph(INT *ia, INT *ja, REAL *a, INT *nrow,	\
	    REAL *rhs, REAL *sol, REAL *exsol)
{
  INT    *jareb=NULL;
  REAL *areb=NULL;
  INT i,nnzlu=-16,n=*nrow,n1=*nrow+1;
  FILE *fp;
  csrreb(&n,&n, &nnzlu,ia,ja,a,&jareb,&areb);
  INT maxja = 5*(n1 + jareb[n]-jareb[0]);
  if(maxja < 7*n) maxja=7*n;
  INT maxa = 5*maxja;
  jareb=(INT *)realloc(jareb,maxja*sizeof(INT));
  areb=(REAL *)realloc(areb,maxa*sizeof(INT));
  INT  nblock = 1;
  INT iblk[2]={0,0};
  iblk[nblock]=n1;
  iblk[nblock-1]=1;
  INT ispd=0;
  INT ncfact=4;
  INT maxlvl=20;
  INT maxfil=32;
  INT method=0;
  REAL dtol=1e-2;
  INT iflag=-16;  
  /// outputs
  INT lvl=-16;
  INT *ka = (INT *) calloc(10*(maxlvl+1),sizeof(INT));
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
  REAL hist[22];
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


