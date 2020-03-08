#include "hazmath.h"
/* coord comparison function */
static INT realcmp(const void *x0, const void *x1)
{
  /* for sorting: coord comparison: returns negative if lexicographically 
     x0< x1 and positive if x0 > x1,   and on rare occasions returns 0*/
  const REAL *a = (const REAL *)x0;
  const REAL *b = (const REAL *)x1;
/* a and b must be vectors with 2 components.*/
  if(a[0]>b[0]) return 1;
  else if (a[0]<b[0]) return -1;
  else if (a[1]>b[1]) return 1;
  else if (a[1]<b[1]) return -1;
  else return 0;
}
/* 
   file to support reading coordinates of points of interest around
   which we may refine. it requires reading/sorting/etc. 
   reads file with features (coordinates of poINTs to refine around
     them). allocates memory 
*/
typedef struct /* segments for constructing initial grid */
{
  INT dim; /* spatial dimension */
  INT np; /* number of points */
  INT ns; /* number of segments */
  iCSRmat *adj; /* adjacency matrix constructed from segments */ 
  REAL *x; /*(np x dim) array to hold the coordinates of the
	     vertices */  
  INT *bndry; /* boundary codes for vertices */
  INT *sbndry; /* boundary codes for segments */
  INT *polarsys; /* polar coord systems for some of the points */
} segments;

typedef struct /* n-homogenous grid: elements isomorphic to
		  parallelepipeds */
{
  INT nbig; /* spatial dimension in which the cubic grid is embedded*/
  INT n; /* dimension of the grid */
  INT np; /* number of vertices */
  INT ns; /* number of n-dimensional parallelepipeds */
  REAL *x; /*(nv x n) array to hold the coordinates of vertices */
  INT *nodes; /*element-vertex incidence (ns by 2^n) */
  INT *flags; /*flag of the hex */
} cubegrid;

INT seg2cube(segments *s,cubegrid *c)
{
  INT ns=s->ns,np=s->np;
  INT *iwrk=2*nseg;
  for(i=0;i<nseg;
  if(use_features){
    char *fname=(char *)calloc(FILENAMESIZE,sizeof(char));  
    size_t lfname=0,lfnamei=0;
    /* read a file with coordinates on the surface where we refine */
    lfname=strlen((char *)FEATURES_DIR);
    lfnamei=strlen((char *)FEATURES_FILE_IN);
    strcpy(fname,(char *)FEATURES_DIR);
    strcpy((fname+lfname),(char *)FEATURES_FILE_IN);
    fname[lfnamei+lfname]='\0';
    feat->fpf=fopen(fname,"r");
    feat->nbig=dim_orig; feat->n=2;/* dim_orig or dim_orig- 1 */
    feat->fill=vfill;
    if(fname) free(fname);
  } else {
    feat->nf=0;
    feat->x=NULL;
    return 5;
  }
  INT dim=feat->n, dimbig=feat->nbig;
/***************************************************/
  INT ichk,k=0,count=0,i,j,m;
  REAL xtmp;
  char ch;
  if(!feat->fpf) {
    feat->nf=0;
    feat->x=NULL;
    fprintf(stderr,"****ERROR:Could not open file for reading in %s; no features read",__FUNCTION__);
    return 3;
  }
  /* 
     Read csv file. When exported via excel such file has 3-4 control
     chars at the beginning.
  */ 
  /*read special chars if any*/
  ch=fgetc(feat->fpf); while((INT )ch < 0){ch=fgetc(feat->fpf);count++;}
  if(count){
    fprintf(stdout,"%%Read: %d control chars...", count);
    fseek(feat->fpf,count*sizeof(char),SEEK_SET);
  } else rewind(feat->fpf);
  k=0;
  while(1){
    if(feof(feat->fpf) || k>40000) break;
    for (j=0;j<dim;j++){
      ichk=fscanf(feat->fpf,"%lg", &xtmp);
      if(ichk<0) break;
      k++;
    }
  }
  k=k/dim;
  /* read now k numbers from the file */ 
  if(count){
    fseek(feat->fpf,count*sizeof(char),SEEK_SET);
  } else {
    rewind(feat->fpf);
  }
  feat->nf=k;
  /* we allocate always the max of dim or dimbig */
  if(dimbig>dim) {
    feat->x=(REAL *)calloc(dimbig*(feat->nf),sizeof(REAL));
    //    fprintf(stdout,"\ndimbig=%d, dim=%d",dimbig,dim);fflush(stdout);
  }else{
    feat->x=(REAL *)calloc(dim*(feat->nf),sizeof(REAL));
  }
  for (i=0;i<k;i++){
    for(j=0;j<dim;j++){
      count = fscanf(feat->fpf,"%lg", (feat->x+dim*i+j));
    }
  }
  fprintf(stdout,"Read %i coord pairs\n",k);
  fclose(feat->fpf);
  //  fclose(feat->fpf);
  /* sort so that duplicates */
  qsort(feat->x,(feat->nf), dim*sizeof(REAL), realcmp);
  /* Clean up the data by removing all the duplicates */
  k=feat->nf-1;
  i=0;j=0;
  REAL *xc=(REAL *)calloc(dim, sizeof(REAL));
  REAL dli,dli1; //l1 distance
  while (i<k){
    if(fabs(feat->x[dim*i])<1e-6) {i++;continue;}
    if(j==0) {
      for(m=0;m<dim;m++){
	feat->x[m]=feat->x[2*i+m];
      }
    }
    dli=0.;
    for(m=0;m<dim;m++){
      xc[m]=feat->x[2*j+m];
      dli=dli+fabs(xc[m]);
    }
    while(1 && i<k) {
      dli1=0.;
      for(m=0;m<dim;m++){
	dli1+=fabs(xc[m]-feat->x[dim*i+dim+m]);
      }
      dli1=dli1/dli;
      if(dli1>1e-6){
	j++;i++;
	for(m=0;m<dim;m++){feat->x[dim*j+m]=feat->x[dim*i+m];}
	break;
      }
      i++;
      //      fprintf(stdout,"i=%i\n",i);
    }
    //    fprintf(stdout,"i=%i, j=%i\n",i,j);
  }
  i++;j++; feat->nf=j;
  return;

void a123()
{
  
  /* Clean up the data by removing all the duplicates */
  k=feat->nf-1;
  i=0;j=0;
  REAL *xc=(REAL *)calloc(dim, sizeof(REAL));
  REAL dli,dli1; //l1 distance
  isi_sortp(nseg, iseg,p,invp);
  for(i=0;i<nseg,i++)jseg[i]=jseg[invp[i]];
  while (i<k){
    if(fabs(feat->x[dim*i])<1e-6) {i++;continue;}
    if(j==0) {
      for(m=0;m<dim;m++){
	feat->x[m]=feat->x[2*i+m];
      }
    }
    dli=0.;
    for(m=0;m<dim;m++){
      xc[m]=feat->x[2*j+m];
      dli=dli+fabs(xc[m]);
    }
    while(1 && i<k) {
      dli1=0.;
      for(m=0;m<dim;m++){
	dli1+=fabs(xc[m]-feat->x[dim*i+dim+m]);
      }
      dli1=dli1/dli;
      if(dli1>1e-6){
	j++;i++;
	for(m=0;m<dim;m++){feat->x[dim*j+m]=feat->x[dim*i+m];}
	break;
      }
      i++;
      //      fprintf(stdout,"i=%i\n",i);
    }
    //    fprintf(stdout,"i=%i, j=%i\n",i,j);
  }
  i++;j++; feat->nf=j;
  //  fprintf(stdout,"i=%i, j=%i\n",i,j);
  for(m=0;m<dim;m++){feat->x[dim*j+m]=feat->x[dim*i+m];}
return
}
/***********************************************************************************************/
/*!
 * \fn SHORT 
 *
 * \brief Transform a dCOOmat matrix to a dCSRmat format.
 *
 * \param A   Pointer to dCOOmat matrix
 * \param B   Pointer to dCSRmat matrix
 *
 * \return    SUCCESS if successed; otherwise, error information.
 *
 */
void (npts,nseg,iseg,jseg,nd,hexc,
{   
  // sort the segments    
    // local variable
    INT *ia = B->IA;
    INT *ja = B->JA;
    REAL *Bval = B->val;
    INT *row_idx = A->rowind;
    INT *col_idx = A->colind;
    REAL *Aval = A->val;
    INT i, iind, jind;
    
    INT *ind = (INT *) calloc(m+1,sizeof(INT));
    
    // initialize
    memset(ind, 0, sizeof(INT)*(m+1));
    
    // count number of nonzeros in each row
    for (i=0; i<nnz; ++i) ind[row_idx[i]+1]++;
    
    // set row pointer
    ia[0] = 0;
    for (i=1; i<=m; ++i) {
        ia[i] = ia[i-1]+ind[i];
        ind[i] = ia[i];
    }
    
    // set column index and values
    for (i=0; i<nnz; ++i) {
        iind = row_idx[i];
        jind = ind[iind];
        ja[jind] = col_idx[i];
        Bval[jind] = Aval[i];
        ind[iind] = ++jind;
    }
    
    if (ind) free(ind);
    
    return SUCCESS;
}

/********************************  END  ********************************************************/
