#if defined (__cplusplus)
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "hazmath.h"

void csrreb(INT *nrow, INT *ncol, INT *nnzluo, \
	    INT *ia, INT *ja, REAL *a,       \
	    INT **jareb, REAL **areb)
{
  INT ipntr,i1,i,j,jk,itstrt,itend,istrt,iend,ni,nii,nnzu=0,nnzlu=0;
  INT n=*nrow,m=*ncol,n1=*nrow+1;
  INT *ikeep=NULL,*itmp=NULL,*jtmp=NULL;
  REAL *x=NULL,*utmp=NULL,*adiag=NULL;
  ///
  INT *ib=NULL,*jb=NULL;
  REAL *b=NULL;
  /* converts a CSR matrix A=(ia,ja,a) with possibly nonsymmetric
     pattern to rebsm format (R.E.Bank and R.K.Smith) and the
     multigraph package (see http://ccom.ucsd.edu/~reb/software.html)
     by R.E.Bank. 
    
     The rebsm format is for a sparse matrices with symmetric
     pattern. Let A be one such matrix.  We write A=D(A) + U(A) +
     U(A^T)^T, where D(A) is the diagonal of A, U(A) is the strict
     lower triangle of A and U(A^T) is the strict lower triangle of
     A^T.  This format is similar to the CSR format. It is described
     by two arrays, i.e. A=(ia,a), but in terms of pointers, we can
     put ja=ia+n+1 and have
     
     ia = array of pointers to beginning of the rows of U(A^T) as they
     are stored in ja. For a CSR storage we have ia(here) =
     ia(csr)+n+1;
     
     ja(1:nnzlu) = column indices in U(A^T) (ordered in increasing order in
              this implementation); 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//no ordering needed in multigraph!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       a[0:n-1] =  D;

      a[n] is arbitrary;

      a[n+1:n+nnzlu] = values of U(A^T) by columns, or, equivalently, the upper
                 triangle of A by rows)

      a[n+nnzlu:n+nnzlu*2]= values of U(A) by columns.

      In CSR terms: U(A^T)=(ia,ja,a+n1) and U(A)=(ia,ja,a+n1+nnzlu).

      The input here is an n by n CSR A(ia,ja,a) and the output is
      B=(ib,jb,b), with jb=ib+n+1. If A originally does not have a symmetric
      sparsity pattern, then it is padded with zeroes so that B has a
      symmetric pattern pattern. the storage for B is allocated here. 
*/
  itmp=(INT *) calloc(n1,sizeof(INT)); // pointers to rows (working)
  /*reset */
  for (i=0; i < n1; ++i) itmp[i] = 0;
  /*  count nonzeroes*/
  /* 1. find nnz(U(A^T)) */
  nnzu = 0;
  for(i=0;i<n;i++){
    istrt  = ia[i];
    iend  = ia[i+1];
    if(iend>istrt)     
      for (jk = istrt;jk<iend; jk++){
	j = ja[jk];
	if( j < i )  {itmp[j+1]++; nnzu++;}
      }
  }
  //  fprintf(stdout,"\nnonzeroes in L^T=U(A^T) in %s: nnzu=%i\n",__FUNCTION__,nnzu);
  //  fflush(stdout);
  jtmp = (INT *)calloc(nnzu,sizeof(INT));
  /* this is initial allocation for the B; if the sparsity is
     originally symmetric, there is no need to reallocate
     reallocated */
  adiag = (REAL *)calloc(n1 + 2*nnzu,sizeof(double));
  utmp = adiag + n1; 
  ni = itmp[1];
  itmp[0]=0;
  itmp[1]=0;
  for (i = 1; i < n; ++i){
    i1=i+1;
    nii = itmp[i1];
    itmp[i1]=itmp[i]+ni;
    ni=nii;
  }
  /* store the L^T=U(A^T) as (itmp,jtmp,utmp)*/
  ipntr=0;
  for(i=0;i<n;i++){
    i1=i+1;
    istrt  = ia[i];
    iend  = ia[i+1];
    if(iend>istrt){
      for (jk = istrt;jk<iend; jk++){
	j = ja[jk];	
	if(j < i) {
	  ipntr = itmp[j+1]++;
	  jtmp[ipntr]=i;
	  utmp[ipntr]=a[jk];
	} else if (i==j) {
	  adiag[i]=a[jk];
	}
      }
    }
  }
  /****************************************** */
  /* 2. Now we count the nonzeroes in U(A) + U(A^T) */
  /* Allocate working array for keeping record of who is in U(A^T) and
     also in U(A) */
  ikeep=(INT *) calloc(m,sizeof(int)); /* (working array: who is in
					  both U(A) and U(A^T)*/
  /* reset */
  for (i=0; i < m; ++i) ikeep[i] = -1;
  /* count first nonzeroes in U(A) + U(A^T) */
  nnzlu = 0;
  for (i=0; i< n; ++i) {
    i1=i+1;
    istrt = ia[i];
    iend = ia[i1];
    if (iend > istrt) {
      for (jk = istrt; jk < iend; ++jk) {
	j = ja[jk];
	if( j <= i) continue;
	++nnzlu;
	ikeep[j] = i;
      }
    }
    itstrt = itmp[i];
    itend  = itmp[i1];
    if (itend > itstrt) {
      for (jk = itstrt; jk < itend; ++jk) {
	j = jtmp[jk];
	if (ikeep[j] != i) {
	  ++nnzlu;
	}
      }
    }
  }
  /* up to here we only computed the number of nonzeroes in U(A)+U(A^T) */
  fprintf(stdout,"nnz(U(A+A^T)) vs. nnz(U(A)+U(A^T)) =  %i vs. %i\n",nnzlu,nnzu); fflush(stdout);
   if(nnzlu>nnzu){
     jtmp=(INT *)realloc(jtmp,nnzlu*sizeof(int));
     adiag=(REAL *)realloc(adiag,(2*nnzlu+n1)*sizeof(double));
   }
  /* still have L^T=U(A^T) in (itmp,jtmp,utmp) but now with enough space
     for storing the U(A) + U(A^T) */
  ib=(INT *)calloc(n1+nnzlu,sizeof(int));
  jb=ib + n1;
  b = adiag + n1; // because of realloc this needs to be here;
		  // otherwise b=utmp will work
  utmp = b + nnzlu;
  /* rearrange to put utmp at the end */
  //  for(i=0;i<nnzu;++i) utmp[i] = b[i];
  /* reset */
  for(i = 0;i<m;++i) ikeep[i]=-1;
  /* Store the nnz structure of U(A) + U(A^T) in (ib,jb) */
  ipntr = 0;
  for (i=0; i < n; ++i) {
    ib[i] = ipntr;
    i1=i+1;
    istrt = ia[i];
    iend = ia[i1];
    if (iend > istrt) {
      for (jk = istrt; jk < iend; ++jk) {
	j = ja[jk];
	if(j > i) {//Use only U(A)
	  jb[ipntr] = j;
	  ++ipntr;
	  ikeep[j] = i;
	}
      }
    }
    itstrt = itmp[i];
    itend = itmp[i1];
    if (itend > itstrt) {
      for (jk = itstrt; jk < itend; ++jk) {
	j = jtmp[jk];
	if (ikeep[j] != i) {
	  jb[ipntr] = j;
	  ++ipntr;
	}
      }
    }
  }  
  ib[n] = ipntr;
  if(ikeep) free(ikeep);
  /* now we have the indexing for U(A)+U(A^T) in ib and jb and
     the old indexing for B in itmp,jtmp; 
     The B holds the entries
     of B (not padded with zeros)*/
  x=(REAL *)calloc(m,sizeof(double));
  /* get U(A^T) padded with zeroes */
  for (i=0;i<n;++i) {
    i1=i+1;
    istrt = ib[i];
    iend = ib[i1];
    if (iend > istrt) {
      for (jk = istrt;jk<iend;++jk) {
	j=jb[jk];
	x[j] = 0e+00;
      }
      /* only do B */
      itstrt = itmp[i];
      itend = itmp[i1];
      for (jk = itstrt;jk<itend;++jk) {
	j = jtmp[jk];
	x[j] = b[jk];
      }
      for (jk = istrt;jk<iend;++jk) {
	j=jb[jk];
	utmp[jk] = x[j];
	//	fprintf(stdout,"(%i,%i)=%f\n",i+1,j+1,utmp[jk]);
      }
    } // if (icstrt > icend)...
  } // loop i=0; i< n
  //
  if(itmp) free(itmp);
  if(jtmp) free(jtmp);
  /*  clean up */
  /* final step is to store U(A) padded with zeroes where U(A^T) has nonzeroes */
  /* U(A) is pointed by utmp */
  
  for (i=0;i<n;++i) {
    i1=i+1;
    istrt = ib[i];
    iend = ib[i1];
    if (iend > istrt) {
      for (jk = istrt;jk<iend;++jk) {
	j=jb[jk];
	x[j] = 0e+00;
      }
      /* only do B */
      itstrt = ia[i];
      itend = ia[i1];
      for (jk = itstrt;jk<itend;++jk) {
	j = ja[jk];
	if( j > i ) x[j] = a[jk];
      }
      for (jk = istrt;jk<iend;++jk) {
	j=jb[jk];
	b[jk] = x[j];
      }
    }
  }
  jb=ib+n1;
  /* finalize */
  for (i=0;i<n1;++i) ib[i] = ib[i]+n1;  
  *areb=adiag;
  *jareb=ib;
  *nnzluo=nnzlu;
  return;
}

#if defined (__cplusplus)
}
#endif

