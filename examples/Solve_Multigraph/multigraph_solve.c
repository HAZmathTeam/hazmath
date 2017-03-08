#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "hazmat.h"
#include "multigraph_solve.h"

int main (int argc, char* argv[])
{
  dCSRmat A; //matrix
  dCOOmat B; //matrix
  dvector rhs,sol,exsol; //RHS,solution,exact solution
  INT idummy[3],i=3,n=-16,ncol=-16;
  // read the matrix
  FILE* fp=HAZ_fopen("symm.coo","r");
  rveci_(fp,idummy, &i);
  n=A.row=B.row=idummy[0];
  ncol=A.col=B.col=idummy[1];
  INT m=ncol,n1=n+1;
  A.nnz=B.nnz=idummy[2];
  B.rowind=(INT *) calloc(B.nnz,sizeof(INT));
  B.colind=(INT *) calloc(B.nnz,sizeof(INT));
  B.val=(REAL *) calloc(B.nnz,sizeof(REAL));
  
  rveci_(fp,B.rowind, &B.nnz);
  rveci_(fp,B.colind, &B.nnz);
  rvecd_(fp,B.val, &B.nnz);
  //rhs and exact sol
  sol.row=rhs.row=exsol.row=n;
  rhs.val=(REAL *) calloc(rhs.row,sizeof(REAL));
  sol.val=(REAL *) calloc(sol.row,sizeof(REAL));
  exsol.val=(REAL *) calloc(exsol.row,sizeof(REAL));
  rvecd_(fp,rhs.val, &n);
  rvecd_(fp,exsol.val, &n);
  for(i=0;i<sol.row;i++) sol.val[i]=0.;
  fclose(fp);
  //shift
  for(i=0;i<B.nnz;i++) B.rowind[i]--;
  for(i=0;i<B.nnz;i++) B.colind[i]--;
  //A in csr format
  A.IA=(INT *) calloc((A.row+1),sizeof(INT));
  A.JA=(INT *) calloc(A.nnz,sizeof(INT));
  A.val=(REAL *) calloc(A.nnz,sizeof(REAL));
  dcoo_2_dcsr (&B,&A);
  if(B.rowind) free(B.rowind); if(B.colind) free(B.colind); if(B.val) free(B.val);
  fprintf(stdout,"\nrows=%i,cols=%i,nnz=%i\n",A.row,A.col,A.nnz);
  // use now multigraph
  clock_t clk_mesh_start = clock(); // Time
  INT    *jareb=NULL;
  REAL *areb=NULL;
  INT nnzlu=-16;  
  fprintf(stdout," *** Starting multigraph Solver\n");
  csrreb(&n,&n, &nnzlu,A.IA,A.JA,A.val,&jareb,&areb);
  INT maxja = 4*(n1 + jareb[n]-jareb[0]);
  if(maxja < 7*n) maxja=7*n;
  INT maxa = 5*maxja;
  jareb=(INT *)realloc(jareb,maxja*sizeof(INT));
  areb=(REAL *)realloc(areb,maxa*sizeof(REAL));
  INT  nblock = 1;
  INT iblk[2]={0,0};
  iblk[nblock]=n1;
  iblk[nblock-1]=1;
  INT *ka = (INT *) calloc(10*(maxlvl+1),sizeof(INT));
  //shift and do mginit
  for (i=0;i<n1+nnzlu;++i)
    jareb[i]+=1;
  mginit_(&n, &ispd, &nblock, iblk,			\
  	  &maxja, jareb, &maxa, areb,			\
  	  &ncfact, &maxlvl, &maxfil, ka,
  	  &lvl, &dtol, &method, 
  	  &iflag );
  fprintf(stdout,"\n\n*** mginit: flag=%i\n\n",iflag);
  fflush(stdout);
  INT j, ij;
  INT ns=ka[10*(lvl+1-1)+(2-1)]-1;
  REAL *hist= (REAL *) calloc(22,sizeof(REAL));
  mg_(&ispd, &lvl, &mxcg, &eps1, jareb, areb,	\
       sol.val, rhs.val, ka, &relerr, &iflag, hist );
  fprintf(stdout,"\n\n*** mg_: flag=%i\n\n",iflag);
  fflush(stdout);
  /* unshift */
  for (i=0;i<n1+nnzlu;++i) jareb[i]-=1;
  INT iend=hist[22], itnum= (INT )hist[iend];
  if(itnum<iend) iend=itnum;
  fprintf(stdout,"\nMultigraph history (iter_num=%5i)\n",itnum);    
  for (i=0;i<iend;++i){    
    fprintf(stdout,"iter =  %5i; res %12.4e\n",itnum-iend+i+1,hist[i]);    
  }
  //  for (i=0;i<n;++i){    
  //    fprintf(stdout,"%22.16e\n",sol.val[i]);    
  //  }
  if(ka) free(ka);
  if(hist) free(hist);
  if(jareb) free(jareb);
  if(areb) free(areb); 
  return 0;
}
#if defined (__cplusplus)
}
#endif


