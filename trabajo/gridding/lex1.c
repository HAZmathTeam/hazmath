#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef REAL
#define REAL double
#endif

#ifndef INT
#define INT int
#endif

void lexsort(const INT nr, const INT nc,REAL *a,INT *p)
{
  /*
    implements STRAIGHT INSERT sorting to order lexicographically nr
    names with nc components each.  the array a is overwritten. aj[] on
    input should be allocated. aj is a working double array with nc
    elements. used to sort coordinates of the vertices of the
    macroelements lexicographically.  on output, a is ordered, p is the
    permutation used to order a. The ogiginal a is recovered with inf
    permutation aorig[]=a[invp[i]];
  */
  INT i,j,k,k1,pj;
  unsigned int lt=0;
  REAL *aj=(REAL *)calloc(nc,sizeof(REAL));
  for (i = 0; i < nr; i++){p[i]=i;} 
  for (j = 1; j < nr; ++j) {
    //    for(k=0;k<nc;k++)aj[k] = a[j*nc+k];
    memcpy(aj,(a+j*nc),nc*sizeof(REAL)); pj = *(p + j);
    i = j - 1;
    lt=0;
    for(k=0;k<nc;k++) {
      if(aj[k]>a[i*nc+k]){lt=0;break;}
      else if(aj[k]<a[i*nc+k]){lt=1;break;}
      else {continue;}	
    }
    //    while ((i >=0) && (aj<a[i])){
    while ((i>=0) && (lt)){
      memcpy((a+(i+1)*nc),(a+i*nc),nc*sizeof(REAL)); *(p+i+1)=*(p+i);
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
    memcpy((a+(i+1)*nc),aj,nc*sizeof(REAL)); *(p+i+1)=pj;
  }
  if(aj) free(aj);
  return;
}
INT main(int argc, char **argv){
  INT i,j,k,n=17,dim=4;
  REAL *a=calloc((n*dim),sizeof(REAL));
  INT *p=(INT *)calloc(n,sizeof(INT));
  FILE *fp=fopen("test4.input","r");
  for (i=0;i<n;i++){
    fprintf(stdout,"a(%i,:)=",i);
    for (j=0;j<dim;j++){
      k=fscanf(fp,"%lg", (a+i*dim+j));    
      fprintf(stdout,"%g ",a[i*dim+j]);
    }
    fprintf(stdout,";\n");    
  }
  fprintf(stdout,"\n");    
  lexsort(n,dim,a,p);
      //  d_sort(n,a);
  for (i=0;i<n;i++){
    fprintf(stdout,"\nasort(%3i(%3i),:)=",i,p[i]);
    for (j=0;j<dim;j++){
      fprintf(stdout,"%g ",a[i*dim+j]);
      //      fprintf(stdout,"%g(%g) ",aord[i*dim+j],aord[i*dim+j]-a[p[i]*dim+j]);
    }
  }
  fprintf(stdout,"\n");
  if(a) free(a);
  if(p) free(p);
  return 0;
}
