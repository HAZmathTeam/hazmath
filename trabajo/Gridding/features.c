#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "hazmath.h"
#include "grid_defs.h"
#include "grid_params.h"
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
INT features_r(const INT dim_orig,const INT use_features,features *feat, REAL vfill)
{
  /* dimbig is the spatial dimension; dim can be either dimbig or
     something smaller as which will indicate the features are on the
     boundary. In all cases, feat->x[] here is dimbig by f->n array */
  /*inputs*/
  if(use_features){
    feat->fpf=NULL;
    char *fname=(char *)calloc(FILENAMESIZE,sizeof(char));  
    size_t lfname=0,lfnamei=0;
    /* read a file with coords where we can refine */
    lfname=strlen((const char *)FEATURES_DIR);
    lfnamei=strlen((const char *)FEATURES_FILE_IN);
    strncpy(fname,(const char *)FEATURES_DIR,lfname);
    strncpy((fname+lfname),FEATURES_FILE_IN,lfnamei);
    feat->fpf=fopen(fname,"r");
    //    if(feat->fpf)
    //      {      fprintf(stdout,"\nZZZZZZZZZZZZvfill=%f; lfname,lfnamei=%li ;  %li ; name=%s\n",vfill,lfname,lfnamei,fname);fflush(stdout);}
    //    else
    //      {fprintf(stderr,"\n*** ERROR***\n");fflush(stderr);}
    /*feat nbig or feat*/
    //    feat->nbig=dim_orig; feat->n=dim_orig;/* dim_orig or dim_orig- 1 */
    feat->fill=vfill;
    if(fname) free(fname);
  } else {
    //    fprintf(stdout,"\nnbig==%d\n",feat->nbig);fflush(stdout);
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
     Read (csv) file. When exported via excel such file has 3-4 control
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
  //  fprintf(stdout,"\nfeatures=%d ; %d ;  %d\n",feat->nf,feat->n,feat->nbig);fflush(stdout);
  /* we allocate always the max of dim or dimbig */
  if(dimbig>dim) {
    feat->x=(REAL *)calloc(dimbig*(feat->nf),sizeof(REAL));
    //    fprintf(stdout,"\ndimbig=%d, dim=%d",dimbig,dim);fflush(stdout);
  }else{
    feat->x=(REAL *)calloc(dim*(feat->nf),sizeof(REAL));
  }
  for (i=0;i<k;i++){
    //    fprintf(stdout,"\nnode=%d ; ",i);
    for(j=0;j<dim;j++){
      count = fscanf(feat->fpf,"%lg", (feat->x+dim*i+j));
      //      fprintf(stdout,"%g ; ",feat->x[dim*i+j]);
    }
  }
  fprintf(stdout,"\nRefinement around points: Read %i coordinate %d-tuples\n",k,dim);
  fclose(feat->fpf);

  /* sort so that no duplicates are present */
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
  //  fprintf(stdout,"i=%i, j=%i\n",i,j);
  for(m=0;m<dim;m++){feat->x[dim*j+m]=feat->x[dim*i+m];}
  /* fprintf(stdout,"\nSorted Coords:\n");  */
  /* for (i=0;i<feat->nf;i++){  */
  /*   fprintf(stdout,"%17.12g %17.12g\n",feat->x[2*i],feat->x[2*i+1]);  fflush(stdout); */
  /* } */
  /* if dimbig is larger than dim, i.e. we have read a 2d array but we are in 3D we need to rearrange feat->x so that it is dimbig by feat->nf */
  /* fprintf(stdout,"\nCoords:\n");  */
  /* for (i=0;i<feat->nf;i++){  */
  /*   fprintf(stdout,"%17.12g %17.12g\n",feat->x[2*i],feat->x[2*i+1]);  fflush(stdout); */
  /* } */
  if(dimbig > dim) {
    k=feat->nf-1;
    for(i=k;i>0;i--){
      for(m=dim-1;m>=0;m--){
	feat->x[dimbig*i+m]=feat->x[dim*i+m];
	//	fprintf(stdout,"(%d,%d): %d  %d\n",i,m,dimbig*i+m,dim*i+m);  fflush(stdout);
      }
      for(m=dim;m<dimbig;m++){
	feat->x[dimbig*i+m]=feat->fill;
      }
    }
  }
  /* fprintf(stdout,"\nCoords:\n"); */
  /* for (i=0;i<feat->nf;i++){  */
  /*   fprintf(stdout,"%17.12g %17.12g %10.4g\n",feat->x[3*i],feat->x[3*i+1],feat->x[3*i+2]);  fflush(stdout); */
  /* }  */
  return 0;
}
INT features_w(features *feat,REAL *extra)
{
  INT dimbig=feat->nbig, dim=feat->n;
  INT i,j;
  size_t lfname=strlen((char *)FEATURES_DIR);
  size_t lfnamei=strlen((char *)FEATURES_FILE_OUT);
  char *fname=(char *)calloc(lfname+lfnamei+1,sizeof(char));
  strcpy(fname,(char *)FEATURES_DIR);
  strcpy((fname+lfname),(char *)FEATURES_FILE_OUT);
  fname[lfnamei+lfname]='\0';
  FILE *fpo=fopen(fname,"w");
  //  FILE *fpo=fopen("w.w","w");  
  if(!fpo) {
    fprintf(stderr,"****ERROR:Could not open %s for writing; no features read",	\
	    fname);
    return 3;
  } else {
    fprintf(stdout,"\nWriting %d coords to %s;\n",feat->nf,fname);
  }
  /* if extra coordinate comes from interpolation, write it here*/
  if(extra){
    //    fprintf(fpo,"%d\n",feat->nf);
    fprintf(fpo,"%10i %10i\n",feat->nf,feat->n+1);
    for (i=0;i<feat->nf;i++){
      for(j=0;j<dimbig;j++){
	fprintf(fpo,"%23.16g ",feat->x[dimbig*i+j]);
      }
      fprintf(fpo,"%23.16g\n",extra[i]);
    }
  } else {
    for (i=0;i<feat->nf;i++){
      for(j=0;j<dimbig;j++){
	fprintf(fpo,"%23.16g ",feat->x[dimbig*i+j]);
      }
      fprintf(fpo,"\n");
    }
  }
  fclose(fpo);
  if(fname) free(fname);
  return 0;
}


