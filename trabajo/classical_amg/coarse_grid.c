#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "hazmath.h"
/*********************************************************************/
INT main(INT argc, char **argv)
{
  INT dummy,ish,n,nnz,k,iread;
  FILE *fp   = NULL;
  INT  *ia   = NULL, *ja    = NULL;
  dCSRmat a;
  ivector anc,inv_ord;
  iCSRmat *lbfs;
  //
  //  fprintf(stdout,"\n%%%%argc=%3d\n",argc);
  //  for(k=0;k<argc;k++)
  //    fprintf(stdout,"\n%%%%arg(%3d)=\'%s\'",k,argv[k]);
  //  fprintf(stdout,"\n%%%%==============\n");
  if(argc>1){
    //    fp=fopen(argv[1],"r");    
    dcoo_read_dcsr(argv[1],&a);
    n=a.row;
    nnz=a.nnz;
    ia=a.IA;
    ja=a.JA;
    //    fclose(fp);
  } else {
    fp=fopen("abc","r");
    iread=fscanf(fp,"%d %d %d",&n,&dummy,&nnz);
    ia=calloc(n+1,sizeof(INT));
    for(k=0;k<=n;k++) iread=fscanf(fp,"%d",ia+k);
    // shift if ia[0] is not 0.
    ish=0;
    if(ia[0]>0) ish=-ia[0];  
    nnz = ia[n]+ish;
    ja=calloc(nnz,sizeof(INT));
    for(k=0;k<nnz;k++) iread=fscanf(fp,"%d",ja+k);
    fclose(fp);
  }
  fprintf(stdout,"\n%%%%n=%3d,nnz=%3d,shift=%3d\n",n,nnz,ish);fflush(stdout);
  if(ish!=0){
    for(k=0;k<=n;k++) ia[k]+=ish;
    for(k=0;k<nnz;k++) ja[k]+=ish;
  }
  //
  anc=ivec_create(n);
  inv_ord=ivec_create(n);
  lbfs=lex_bfs(n,ia,ja,&inv_ord,&anc);
  //
  fprintf(stdout,"\np=[");
  for(k=0;k<n;k++){
    fprintf(stdout,"%3d ",lbfs->JA[k]+1);
  }
  fprintf(stdout,"];\n");
  fprintf(stdout,"\npinv=[");
  for(k=0;k<n;k++){
    fprintf(stdout,"%3d ",inv_ord.val[k]+1);
  }
  fprintf(stdout,"];\n");
  fprintf(stdout,"\nlevel=[");
  for(k=0;k<n;k++){
    fprintf(stdout,"%3d ",lbfs->val[lbfs->JA[k]]+1);
  }
  fprintf(stdout,"];\n");
  fprintf(stdout,"\ntree=[");
  for(k=0;k<n;k++){
    fprintf(stdout,"%3d ",anc.val[lbfs->JA[k]]+1);
  }
  fprintf(stdout,"];\n");
  icsr_free(lbfs);
  ivec_free(&anc);
  ivec_free(&inv_ord);
  if(argc>1){
    dcsr_free(&a);
    ivec_free(&inv_ord);
    ivec_free(&anc);
  }else{
    free(ia);
    free(ja);
  }
    return 0;
}
//
