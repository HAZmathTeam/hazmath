#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef REAL
#define REAL double
#endif

#ifndef INT
#define INT int
#endif

void d_sort(INT n, REAL *a)
{
//implements STRAIGHT  INSERT sort for integers
  INT i,j;
  REAL aj;
  unsigned int lt=0;
  for (j = 1; j < n; ++j) {
    aj = a[j];
    i = j - 1;
    if(aj>a[i]) lt=0;
    else lt=1;
    //      while ((i >=0) && (aj<a[i])){
    while ((i >=0) && (lt)){
      a[i + 1] = a[i];
      --i;
      if(i<0) break;
      if(aj>a[i]) lt=0;
      else lt=1;
    }
    a[i + 1] = aj; 
  }
  return;
}
void lex_sort(const INT nr, const INT nc, REAL *a, REAL *aj)
{
/*
  implements STRAIGHT INSERT sorting to order lexicographically nr
  names with nc components each.  the array a is overwritten. aj[] on
  input should be allocated. aj is a working double array with nc
  elements.
*/
  INT i,j,k,k1;
  unsigned int lt=0;
  for (j = 1; j < nr; ++j) {
    //    for(k=0;k<nc;k++)aj[k] = a[j*nc+k];
    memcpy(aj,(a+j*nc),nc*sizeof(REAL));
    i = j - 1;
    lt=0;
    for(k=0;k<nc;k++) {
      if(aj[k]>a[i*nc+k]){lt=0;break;}
      else if(aj[k]<a[i*nc+k]){lt=1;break;}
      else {continue;}	
    }
    //    while ((i >=0) && (aj<a[i])){
    while ((i>=0) && (lt)){
      memcpy((a+(i+1)*nc),(a+i*nc),nc*sizeof(REAL));
      //      for(k=0;k<nc;k++) a[(i + 1)*nc+k] = a[i*nc+k];
      --i;
      if(i<0) break;
      lt=0;
      for(k=0;k<nc;k++) {
	if(aj[k]>a[i*nc+k]){lt=0;break;}
	else if(aj[k]<a[i*nc+k]){lt=1;break;}
	else{continue;}
      }
    }
    memcpy((a+(i+1)*nc),aj,nc*sizeof(REAL));
    //    for(k=0;k<nc;k++) a[(i + 1)*nc+k] = aj[k]; 
  }
  return;
}
INT main(int argc, char **argv){
  INT i,j,k,n=16,dim=4;
  REAL *a=calloc(2*(n*dim),sizeof(REAL));
  REAL *wrk=a+n*dim;
  FILE *fp=fopen("test.input","r");
  for (i=0;i<n;i++){
    fprintf(stdout,"a(%i,:)=",i);
    for (j=0;j<dim;j++){
      k=fscanf(fp,"%lg", (a+i*dim+j));    
      fprintf(stdout,"%g ",a[i*dim+j]);
    }
    fprintf(stdout,";\n");    
  }
  fprintf(stdout,"\n");    
      lex_sort(n,dim,a,wrk);
      //  d_sort(n,a);
  for (i=0;i<n;i++){
    fprintf(stdout,"\nasort(%i,:)=",i);
    for (j=0;j<dim;j++){
      fprintf(stdout,"%g ",a[i*dim+j]);
    }
  }
  fprintf(stdout,"\n");    
  return 0;
}
