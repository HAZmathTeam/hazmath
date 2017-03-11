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
  void mg_(INT *ispd, INT *lvl, INT *mxcg, REAL *eps1,	\
	   INT *jareb, REAL *areb,				\
	   REAL *sol, REAL *rhs, INT *ka, REAL *relerr,		\
	   INT *iflag, REAL *hist );
  /*   
       subroutine mg(ispd,lvl,mxcg,eps1,ja,a,dr,br,ka,relerr,
       +      iflag,hist)
  */
  // io routines for debugging. 
  void zwrijv(INT *ia, INT *ja, REAL *a,		\
	      INT *nrow, INT *ncol,const char *fname)
  {
    /* prints a csr matrix in matlab output*/
    INT i,j1,j2,j,nr=*nrow,nc=*ncol,nnz=ia[nr]-1; 
    INT shift_flag = 1; /* INDEXING STARTS AT 0 so for output we need to shift */
    /* write everything in three long rows */
    FILE *fp=fopen(fname,"w+");
    if(!fp)
      fp=stdout;
    for(i=0;i<nr;i++) {
      j1 = ia[i];
      j2 = ia[i+1];
      for(j=j1;j<j2;j++) {
	fprintf(fp," %i",i+shift_flag);
      }
    }
    fprintf(fp,"\n");
    for(i=0;i<nr;i++) {
      j1 = ia[i];
      j2 = ia[i+1];
      for(j=j1;j<j2;j++) {
	fprintf(fp," %i",ja[j]+shift_flag);
      }
    }
    fprintf(fp,"\n");
    for(i=0;i<nr;i++) {
      j1 = ia[i];
      j2 = ia[i+1];
      for(j=j1;j<j2;j++) {
	fprintf(fp," %22.16e",a[j]);
      }
    }
    fprintf(fp,"\n");
    fclose(fp);
    return;
  }
  INT swrveci (FILE *out, INT *vec, INT *n,const INT sht);
  INT swrvecd (FILE *out, REAL *vec, INT *n, const INT sht);
  INT rveci (FILE *inp, INT *vec, INT *n, const INT sht);
  INT rvecd (FILE *inp, REAL *vec, INT *n, const INT sht);

  ///////////////////////////////////////////////////////////////////////////

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
    adiag = (REAL *)calloc(n1 + 2*nnzu,sizeof(REAL));
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
    //  fprintf(stdout,"\nXXXXXXXXXXnonzeroes in U(A)+U(A^T) in %s: nnzLU=%i\n\n",__FUNCTION__,ib[n]);
    //  fflush(stdout);
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
    //  fprintf(stdout," %i %i\n",n,m);
    //  for (i=0;i<nnzlu;++i) {utmp[i]=0e0;}
  
    for (i=0;i<n;++i) {
      i1=i+1;
      istrt = ib[i];
      iend = ib[i1];
      if (iend > istrt) {
	for (jk = istrt;jk<iend;++jk) {
	  j=jb[jk];
	  //	fprintf(stdout," %i %i\n",i,j);
	  x[j] = 0e+00;
	}
	//      fflush(stdout);
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
    //  zwrijv(ia, ja,a,&n,&n,"AX.dat");
    //  zwrijv(ib, jb,b,&n,&n,"AU2.dat");
    //  zwrijv(ib, jb,b+1,&n,&n,"AU1.dat");
    //  fprintf(stdout,"\nibibibibib=%i\n",ib[n]);
    //  fprintf(stdout,"\nbububu=%e\n",b[n1-1]);
    //  fprintf(stdout,"\nblblblb=%e\n",b[n1-2]);
    /* finalize */
    for (i=0;i<n1;++i) ib[i] = ib[i]+n1;  
    *areb=adiag;
    *jareb=ib;
    *nnzluo=nnzlu;
    return;
  }
  void mgraph(INT *ia, INT *ja, REAL *a, INT *nrow,	\
	      REAL *rhs, REAL *sol, REAL *exsol)
  {
    INT    *jareb=NULL;
    REAL *areb=NULL;
    INT i,nnzlu=-16,n=*nrow,n1=*nrow+1;
    FILE *fp;
    csrreb(&n,&n, &nnzlu,ia,ja,a,&jareb,&areb);
    ///
    /*   INT mdata[3]={n,n,n1+jareb[n]-jareb[0]}; */
    /*   INT nout=3; */
    /*   fp = fopen("af77.dat","w"); */
    /*   swrveci(fp,mdata,&nout,0); */
    /*   nout=n1+jareb[n]-jareb[0]; */
    /*   swrveci(fp,jareb,&nout,1); */
    /*   nout = n1+2*(jareb[n]-jareb[0]); */
    /*   swrvecd(fp,areb,&nout,1); */
    /*   nout=n; */
    /*   swrvecd(fp,rhs,&nout,1); */
    /*   swrvecd(fp,exsol,&nout,1); */
    /*   fclose(fp); */
  
    //  fprintf(stdout,"\nn,nnzlu maxja= %i %i %i\n",n,nnzlu,jareb[n]-jareb[0]);
    // guess for the storage at all levels.
    // triple the memory, because this is to say
    // that the storage we need for all levels triples. 
    INT maxja = 5*(n1 + jareb[n]-jareb[0]);
    // this is because it does not work for a diagonal matrix
    if(maxja < 7*n) maxja=7*n;
    INT maxa = 5*maxja;
    jareb=(INT *)realloc(jareb,maxja*sizeof(INT));
    areb=(REAL *)realloc(areb,maxa*sizeof(REAL));
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
    //  fprintf(stdout,"\n\n*** ns=%i", ns);
    REAL hist[22];
    REAL relerr=1e0;  
    /*  lenz = 11*n+2*ns;
	z = (REAL *)realloc(z,lenz*sizeof(double));
    */
    /* fprintf(stdout,"lvl=%i; ns = %i\n",lvl,ns); */
    /* for (i=0;i<10;++i){ */
    /*   for (j=0;j<=lvl;++j){ */
    /*     ij=i+j*10; */
    /*     fprintf(stdout,"  ka(%i,%i)=%i ; ",i+1,j+1,ka[ij]); */
    /*   } */
    /*   fprintf(stdout,"\n"); */
    /* } */
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

