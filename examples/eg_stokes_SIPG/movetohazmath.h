static void csr2matlab_code(FILE *fp,dCSRmat *a, const char *varname)
{
  char *s;
  if((varname!="") && (varname))
    s=strndup(varname,8);
  else
    s=strdup("a");
  //
  fprintf(fp,"\n%%\n%s=[",s);
  csr_print_matlab(fp,a);
  fprintf(fp,"\n];\n%s=sparse(%s(:,1),%s(:,2),%s(:,3),%d,%d);",s,s,s,s,a->row,a->col);
  fprintf(fp,"\n%%\n");
  free(s);
  return;
}
/*********************************************************/
static void bdcsr_extend(block_dCSRmat *ba,				\
			 REAL *vcol, REAL *vrow,			\
			 const INT iblock,				\
			 const REAL c1, const REAL c2)
{
  /*
     if vrow and vcol are NULL, then adds one row and one column of
     zeros at block (iblock,iblock) also reallocs the corresponding
     blocks->IA, blocks->JA, and blocks->val. They are not reallocated
     (except blocks->IA) if vrow and vcol are null because we add
     column and row of zeros. --ltz 20200405
  */

  //
  INT k,j,n,nnz;
  REAL *vals=NULL;
  dCSRmat *atmp=NULL,*ctmp=NULL,*add=NULL;
  INT ptype=-1;
  FILE *fp;///=stdout;
  if(vcol==NULL && vrow==NULL){
    for(k=0;k<ba->brow;k++){
      // just add rows and columns of zeros to the matrix
      atmp=ba->blocks[iblock*ba->bcol+k];
      atmp->row++;
      atmp->IA=realloc(atmp->IA,sizeof(INT)*(atmp->row+1));
      atmp->IA[atmp->row]=atmp->nnz;
      atmp=ba->blocks[k*ba->bcol+iblock];
      atmp->col++;
    }
    return;
  }
  if(vcol==NULL){
    vals=vrow;
    ptype=2;
  } else if(vrow==NULL){
    vals=vcol;
    ptype=1;
  } else
    ptype=0;
  fprintf(stdout,"\n**** ptype=%d\n",ptype);
  if(ptype){
    // block rows
    atmp=ba->blocks[iblock*ba->bcol+iblock];
    n=atmp->row;
    nnz=n;
    // new matrix to add to the existing one
    add=dcsr_create_p(n,n,nnz);
    add->IA[0]=0;
    for(k=0;k<n;k++){
      add->JA[add->IA[k]]=k;
      add->val[add->IA[k]]=vals[k];
      add->IA[k+1]=add->IA[k]+1;
    }
  } else {
    for(k=0;k<ba->brow;k++){
      // block pressure rows
      atmp=ba->blocks[iblock*ba->bcol+k];
      atmp->row++;
      atmp->IA=realloc(atmp->IA,sizeof(INT)*(atmp->row+1));
      atmp->IA[atmp->row]=atmp->nnz;
      atmp=ba->blocks[k*ba->bcol+iblock];
      atmp->col++;
    }
    atmp=ba->blocks[iblock*ba->bcol+iblock];
    n=atmp->row;
    nnz=2*(n-1);
    /* fp=fopen("debug/x123.m","w"); */
    /* csr2matlab_code(fp,atmp,"xtmp"); */
    /* fclose(fp); */
    /* fprintf(stdout,"\nnnz=%d;ncol=%d,nrow=%d\n\n",atmp->nnz,atmp->col,atmp->row);       */
    /* fprintf(stdout,"\nn=%d;n=%d,nnz=%d\n\n",n,n,nnz); */
    /* fflush(stdout); */
    add=dcsr_create_p(n,n,nnz);
    add->IA[0]=0;
    add->val[0]=vcol[0];
    for(k=0;k<(n-1);k++){
      add->JA[add->IA[k]]=(n-1);
      add->val[add->IA[k]]=vcol[k];
      add->IA[k+1]=add->IA[k]+1;
    }
    add->IA[n]=nnz;
    //    fprintf(stdout,"\nn=%d,add->IA[n,n-1]=(%d,%d); (c1,c2)=(%e,%e)\n\n",n,add->IA[n],add->IA[n-1],c1,c2);
    j=0;
    for(k=add->IA[n-1];k<(add->IA[n]);k++){
      add->JA[k]=j;
      add->val[k]=vrow[j];
      j++;
    }
  }
  // now we add to atmp and then copy back to atmp;
  ctmp=malloc(1*sizeof(dCSRmat));
  // ctmp=c1*atmp+c2*add;
  dcsr_add(atmp,c1,add,c2,ctmp);
  free(add);
  atmp->IA=(INT *)realloc(atmp->IA,sizeof(INT)*(ctmp->row+1));
  atmp->JA=(INT *)realloc(atmp->JA,sizeof(INT)*ctmp->nnz);
  atmp->val=(REAL *)realloc(atmp->val,sizeof(REAL)*ctmp->nnz);
  // copy atmp=ctmp
  dcsr_cp(ctmp,atmp);
  /////////////////// print to debug/*
  /* fp=fopen("debug/app.m","w"); */
  /* csr2matlab_code(fp,atmp,"app"); */
  /* fclose(fp); */
  /* ///////////////////// */
  /* fprintf(stdout,"\nnnz-c=%d;nnz-a=%d\n",ctmp->nnz,atmp->nnz); fflush(stdout); */
  /* fp=fopen("debug/a.data","w"); */
  /* bdcsr_print_matlab(fp,ba); */
  /* fclose(fp); */
  //////////////////end print for debug
  // free the working memory.
  free(ctmp);
  return;
}
