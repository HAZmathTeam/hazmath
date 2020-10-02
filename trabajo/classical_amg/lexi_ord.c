#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#ifndef INT
#define INT int
#endif

#ifndef REAL
#define REAL double
#endif

/**
 * \struct dCSRmat
 * \brief Sparse matrix of REAL type in CSR format
 *
 * CSR Format (IA,JA,A) in REAL
 *
 * \note The starting index of A is 0.
 */
typedef struct dCSRmat{

    //! row number of matrix A, m
    INT row;

    //! column of matrix A, n
    INT col;

    //! number of nonzero entries
    INT nnz;

    //! integer array of row pointers, the size is m+1
    INT *IA;

    //! integer array of column indexes, the size is nnz
    INT *JA;

    //! nonzero entries of A
    REAL *val;

} dCSRmat; /**< Sparse matrix of REAL type in CSR format */

/**
 * \struct dCOOmat
 * \brief Sparse matrix of REAL type in COO (or IJ) format
 *
 * Coordinate Format (I,J,A)
 *
 * \note The starting index of A is 0.
 */
typedef struct dCOOmat{

    //! row number of matrix A, m
    INT row;

    //! column of matrix A, n
    INT col;

    //! number of nonzero entries
    INT nnz;

    //! integer array of row indices, the size is nnz
    INT *rowind;

    //! integer array of column indices, the size is nnz
    INT *colind;

    //! nonzero entries of A
    REAL *val;

} dCOOmat; /**< Sparse matrix of REAL type in COO format */


/**
 * \struct dLILmat
 * \brief Sparse matrix of REAL type in Linked Lists D. Knuth format
 *
 * Linked List Format (dCOOmat A,NEXTR,NEXTC,RSTRT,NSTRT) in REAL
 *
 * \note
 */
typedef struct dLILmat{
  
  //! dLIL mat has a dCOOmat as a sub-struct
  dCOOmat *acoo;
  //
  INT *nextr; // array next[] for rows
  
  INT *nextc; // array next[] for columns

  INT *rstrt; // beginning of each row

  INT *cstrt; // beginning of the linked list for each column

} dLILmat; /**< Sparse matrix of REAL type in CSR format */

/*csr2lil*/
dLILmat *csr2lil(dCSRmat *a)
// convert csr matrix to d.knuth's storage scheme which uses linked lists and is a bit expensive but efficient in adding/removing elements, etc. 
{  
  INT nnz=a->nnz;
  INT n=a->row;
  INT m=a->col;
  dLILmat *alil=malloc(sizeof(alil));
  alil->acoo=malloc(sizeof(dCOOmat));
  dCOOmat *acoo=alil->acoo;
  void *coo=calloc(nnz*sizeof(REAL)/sizeof(char) + (4*nnz + 2*n*sizeof(INT))/sizeof(char),sizeof(char));
  acoo->val=(REAL *)coo;
  acoo->rowind=(INT *)(coo + nnz*sizeof(REAL));
  acoo->colind=acoo->rowind + nnz;
  alil->nextr=acoo->colind   + nnz;
  //
  alil->nextc=alil->nextr + nnz;
  alil->rstrt=alil->nextc + nnz;
  alil->cstrt=alil->rstrt + n;
  void *cooend=(void *)(alil->rstrt + m);
  REAL *aa=(REAL *)(coo + 2*nnz*sizeof(INT));
  return alil;
}
/**********************************************************************/
void fill(INT n, INT *ia, INT *ja, INT *iord, INT *iord1)
{
  /*
    calculates the fill-in in the LU decomposition for a given ordering iord. Algorithm is
    linear time algorithm in number of edges in the elimination graph.
  */
  return;
}
/**********************************************************************/
void lexi(INT n, INT *ia, INT *ja, INT *iord, INT *iord1)
{  /*
     Output: iord is the lexicographical ordering which for
     triangulated graph is an ordering withouut fill in
   */
  /* 
   *    Lexicographical (recurrsive) breath first search.       
   *	Following the algol-like pseudocode from 
   *	Rose, Donald J. ; Tarjan, R. Endre ; Lueker, George S. 
   *	Algorithmic aspects of vertex elimination on graphs. 
   *	SIAM J. Comput. 5 (1976), no. 2, 266–283 (MR0408312).
  */
  INT c,h,v,w,p;  
  INT ibeg,iend,i,iv,nfx,ijk,i0;
  /* 
   *  two lists: queue and sets of vertices. the queue is described by
   *  head[] and backp[] and the sets are described by next[] and
   *  back[]. Each cell "c" is either a header cell or describes a
   *  set of vertices that have the same label.
   *
  */
  INT n2=2*n;
  // allocate memory;
  INT *mbegin=calloc(6*n2,sizeof(INT));
  memset(mbegin,0,6*n2*sizeof(INT));
  // set pointers;
  INT *fixlst = mbegin;
  INT *head   = fixlst + n2;
  INT *back   = head   + n2;
  INT *next   = back   + n2;
  INT *cell   = next   + n2;
  INT *flag   = cell   + n2;
  INT *mend   = flag   + n2;
  /*     assign empty label to all vertices */
  head[0] = 1; // the first cell is the head of the queue; it does not
	       // describe any set of vertices, but points to the
	       // cell=1 which will describe the first set of
	       // vertices.
  back[1] = 0; // the second cell points to the first cell as its predecessor
  /**/
  head[1] = -1;   
  back[0] = -1;
  next[0] = -1;
  flag[0] = -1;
  flag[1] = -1;
  // cell 0 is the beginning; cell 1 is header cell for the first set;
  c=2; // thus the first empty cell is cell 2. 
  for(iv = 0;iv<n;iv++){
    v = iord[iv]; // vertex label
    head[c] = v;  //head(cell)=element;
    cell[v] = c;  //cell(vertex)=cell.
    next[c-1] = c;
    flag[c] = 1;
    back[c] = c-1;
    c++;
    iord1[v] = -1;
  }
  next[c-1] = -1;
  for (i = n-1;i>=0;i--){
    //C  Skip empty sets
    while(next[head[0]] < 0){
      head[0] = head[head[0]];
      back[head[0]]=0;
    }
    //C  pick next vertex to number
    p=next[head[0]];
    //C     C ADDITION BY ltz
    v = head[p];
    //    fprintf(stdout,"\np=%d,v=%d,cell[v]=%d",p,v,cell[v]);fflush(stdout);
    //C  delete cell of vertex from set
    next[head[0]] = next[p];
    next[back[cell[v]]]=next[cell[v]];
    if(next[cell[v]]>=0){
      back[next[cell[v]]]=back[cell[v]];
    }
    //C assign number i to v
    iord[i] = v;
    iord1[v] = i;
    nfx=0;
    ibeg = ia[v]-1;
    iend = ia[v+1]-1;
    for(ijk = iend;ijk>ibeg;ijk--){
      w = ja[ijk];
      if(iord1[w] < 0) {
	// delete cell of w from the set (this will go into new set. 
	next[back[cell[w]]]=next[cell[w]];
	if(next[cell[w]] >=0){
	  back[next[cell[w]]]=back[cell[w]];
	}
	h = back[flag[cell[w]]];
	// if h is an old set, then we create a new set (c=c+1)
	if(flag[h]<0) {
	  head[c] = head[h];
	  head[h] = c;
	  back[head[c]] = c;
	  back[c] = h;
	  flag[c] = 0;
	  next[c] = -1;
	  // add the new set to fixlst:
	  nfx++;
	  fixlst[nfx] = c;
	  h=c;
	  fprintf(stdout,"\nc=%d;fixlst[%d]=%d",c,nfx,fixlst[nfx]);fflush(stdout);
	  c++;
	}
	// add cell of w to the new set
	next[cell[w]] = next[h];
	if (next[h] >= 0) {
	  back[next[h]] = cell[w];
	}
	flag[cell[w]] = h;
	back[cell[w]] = h;
	next[h] = cell[w];
      }//	end if
    }   //end for
    for(i0 = 0;i0< nfx;i0++){
      //      h = fixlst[i0];
      //      flag[h] = -1;
      flag[fixlst[i0]] = -1;
    }// 
  }  //end do
  // free
  free(mbegin);
  return;
}
/*********************************************************************/
INT main(INT argc, char **argv)
{
  INT n,nnz,k,iread=-1,i0=-1,i1=-1;
  FILE *fp   = NULL;
  INT  *ia   = NULL, *ja    = NULL;
  INT  *iord = NULL, *iord1 = NULL;
  //
  fprintf(stdout,"\nData type? ");
  fscanf(stdin,"%d",&iread);
  if(iread) {
    fp=fopen("fort.15","r");
    fscanf(fp,"%d",&n);
    ia=calloc(n+1,sizeof(INT));
    for(k=0;k<n;k++){
      fscanf(fp,"%d",ia+k);
      ia[k]--;
    }
    i0 = ia[0]+1;
    ia[0]  = 0;
    for(k=0;k<n;k++){
      i1 = ia[k+1] + 1;
      ia[k+1] = ia[k] + i0;
      i0 = i1;
      fprintf(stdout,"\nia[%d]=%d",k,ia[k]);fflush(stdout);
    }
    nnz = ia[n];
    fprintf(stdout,"\nn=%d,nnz=%d\n",n,nnz);fflush(stdout);
    ja=calloc(nnz,sizeof(INT));
    for(k=0;k<nnz;k++){
      fscanf(fp,"%d",ja+k);
      fprintf(stdout,"\nja[%d]=%d",k,ja[k]);fflush(stdout);
      ja[k]--;
    }
  } else {
    fp=fopen("abc","r");
    fscanf(fp,"%d",&n);
    ia=calloc(n+1,sizeof(INT));
    for(k=0;k<=n;k++){
      fscanf(fp,"%d",ia+k);
      fprintf(stdout,"\nia[%d]=%d",k,ia[k]);fflush(stdout);
      ia[k]--;
    }
    nnz = ia[n];
    ja=calloc(nnz,sizeof(INT));
    fprintf(stdout,"\nn=%d,nnz=%d\n",n,nnz);fflush(stdout);
    for(k=0;k<nnz;k++){
      fscanf(fp,"%d",ja+k);
      fprintf(stdout,"\nja[%d]=%d",k,ja[k]);fflush(stdout);
      ja[k]--;
    }
  }//end if
  fclose(fp);
  iord=calloc(n,sizeof(INT));
  for(k=0;k<n;k++){
    iord[k] = k;
  }
  iord1=calloc(n,sizeof(INT));
  memset(iord1,0,n*sizeof(INT));
  // memory:
  //
  lexi(n,ia,ja,iord, iord1);
  //
  fprintf(stdout,"\n");
  for(k=0;k<n;k++){
    fprintf(stdout,"%d ",iord[k]);
  }
  fprintf(stdout,"\n");
  for(k=0;k<n;k++){
    fprintf(stdout,"%d ",iord1[k]);
  }
  fprintf(stdout,"\n");
  free(iord);
  free(iord1);
  free(ia);
  free(ja);
  return 0;
}
//
