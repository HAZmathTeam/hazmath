/*! \file src/amr/input_grid.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note containing all essential routines for input for the mesh
 *  generation mesh refinement
 *
 */
#include "header_input.h"
/*********************************************************************/
void r2c1(const INT n, const INT m, const size_t sizeel, void *x)
{
  INT i,j,ji,nms=n*m*sizeel;
  void *y=(void *)malloc(nms);
  for (i=0;i<n;i++){
    for (j=0;j<m;j++){
      ji=sizeel*(n*j+i);
      memcpy((y+ji),x,sizeel);
      x+=sizeel;
    }
  }
  memcpy((x-nms),y,nms);
  if(y) free(y);
  return;
}
/*************************************************/
char **input_strings(INT *nall_out)
{
  /* take as input the strings in INPUT_GRID_DATA and on output has an
     array with all such strings which are used in parsing the input
     file for grid generation. */
  INT k,nall;
  const char *indata_const[]= { INPUT_GRID_DATA_ ,	\
				"\0"		\
  };
  nall=0;
  while(strlen(indata_const[nall]))nall++;
  //  fprintf(stdout,"\nthe nall is: %d\n",nall);
  char **indata=malloc(nall*sizeof(char *));
  for(k=0;k<nall;k++){
    indata[k]=strndup(indata_const[k],strlen(indata_const[k]));   
  }
  *nall_out=nall;
  return indata;
}
/*********************************************************************/
void input_grid_arrays(input_grid *g)
{
  //allocate the arrays used in input grid
  INT nvcube=(1 << g->dim),nvface=(1 << (g->dim-1));
  g->ox=(REAL *)calloc(g->dim*g->ncsys,sizeof(REAL));
  g->systypes=(INT *)calloc(g->ncsys,sizeof(INT)); 
  g->syslabels=(INT *)calloc(g->ncsys,sizeof(INT)); 
  g->csysv=(INT *)calloc(g->nv,sizeof(INT)); 
  g->labelsv=(INT *)calloc(g->nv,sizeof(INT)); 
  g->bcodesv=(INT *)calloc(g->nv,sizeof(INT)); 
  g->xv=(REAL *)calloc(g->dim*g->nv,sizeof(REAL)); 
  g->xe=(REAL *)calloc(g->dim*g->ne,sizeof(REAL)); 
  g->seg=(INT *)calloc(3*g->ne,sizeof(INT));
  g->mnodes=(INT *)calloc(g->ne*(nvcube+1),sizeof(INT));
  g->mfaces=(INT *)calloc(g->nf*(nvface+1),sizeof(INT));
  return;
}
/******************************************************************/
void input_grid_free(input_grid *g)
{
  //free input grid structure. 
  if(g->title) free(g->title);
  if(g->dgrid) free(g->dgrid);
  if(g->fgrid) free(g->fgrid);
  if(g->dvtu) free(g->dvtu);
  if(g->fvtu) free(g->fvtu);
  if(g->ox) free(g->ox);
  if(g->systypes) free(g->systypes);
  if(g->syslabels) free(g->syslabels);
  if(g->csysv) free(g->csysv);
  if(g->bcodesv) free(g->bcodesv);
  if(g->xv) free(g->xv);
  if(g->xe) free(g->xe);
  if(g->seg) free(g->seg);
  if(g->mnodes) free(g->mnodes);
  if(g->mfaces) free(g->mfaces);
  if(g) free(g);
  return;
}
/**********************************************************************/
void input_grid_print(input_grid *g)
{
  // prints what was read from the inpput grid file. 
  INT i,j,dim=g->dim;
  fprintf(stdout,"\n\nTITLE: %s",g->title);
  fprintf(stdout,"\ndimension=%d",g->dim);
  fprintf(stdout,"\ndir_grid=%s",g->dgrid);
  fprintf(stdout,"\ndir_vtu=%s",g->dvtu);
  fprintf(stdout,"\nfile_grid=%s",g->fgrid);
  fprintf(stdout,"\nfile_vtu=%s",g->fvtu);
  fprintf(stdout,"\nprint_level=%d",g->print_level);
  fprintf(stdout,"\nnum_refinements=%d",g->nref);
  fprintf(stdout,"\nrefinement_type=%d",g->ref_type);
  fprintf(stdout,"\nerr_stop_amr=%.3g",g->err_stop);
  /*ARRAYS*/
  fprintf(stdout,"\n\nnum_coordsystems=%d",g->ncsys);
  for(i=0;i<g->ncsys;i++){
    fprintf(stdout,"\nlabel=%d,type=%d, origin(",g->syslabels[i],g->systypes[i]);
    for(j=0;j<g->dim;j++) fprintf(stdout," %6.2f ",g->ox[i*dim+j]);
    fprintf(stdout,")");
  }
  fprintf(stdout,"\n\nnum_vertices=%d",g->nv);
  for(i=0;i<g->nv;i++){
    fprintf(stdout,"\nvertex=%d, coord_system=%d, bcode=%d, coords(",i,g->csysv[i],g->bcodesv[i]);
    if(g->systypes[g->csysv[i]]==1){
      fprintf(stdout," %6.2f ",g->xv[i*dim]);
      for(j=1;j<g->dim;j++)	
	fprintf(stdout," %6.2f ",(g->xv[i*dim+j])/((REAL )PI)*180.);
    }else{
      for(j=0;j<g->dim;j++) fprintf(stdout," %6.2f ",g->xv[i*dim+j]);
    }
    fprintf(stdout,")");
  }
  fprintf(stdout,"\n\nnum_edges=%d\n",g->ne);
  for(i=0;i<g->ne;i++){
    fprintf(stdout,"\nedge=(%d,%d) div=%d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
  }
  INT nvcube=(1<<g->dim),nvcube1=nvcube+1;
  fprintf(stdout,"\n\nnum_macroelements=%d\n",g->nel);
  for(i=0;i<g->nel;i++){
    fprintf(stdout,"\nmacroel=%d; code=%d; vertices=(",i,g->mnodes[i*nvcube1+nvcube]);
    for(j=0;j<nvcube;j++){
      fprintf(stdout,"%d ",g->mnodes[nvcube1*i+j]);
    }
    fprintf(stdout,")");
  }
  INT nvface=(1<<(g->dim-1)),nvface1=nvface+1;
  fprintf(stdout,"\n\nnum_faces=%d\n",g->nf);
  for(i=0;i<g->nf;i++){
    fprintf(stdout,"\nmacroface=%d; code=%d; vertices=(",i,g->mfaces[i*nvface1+nvface]);
    for(j=0;j<nvface;j++){
      fprintf(stdout,"%d ",g->mfaces[nvface1*i+j]);
    }
    fprintf(stdout,")");
  }
  fprintf(stdout,"\n\n");fflush(stdout);  
  return;
}
/**********************************************************************/
char **splits(char *s, const char *d, INT *num)
{
  // splits a string in an array. splitting delimiter is d
  INT k;
  char **w=calloc(strlen(s)+1,sizeof(char *));
  char *work=strtok(s,d);
  k=0;
  while(work!= NULL){	
    w[k]=calloc(strlen(work)+1,sizeof(char));
    memcpy(w[k],work,(strlen(work)+1)*sizeof(char)); // +1 to copy \0
						     // at the end;
    work = strtok(NULL,d);
    k++;
  }
  w[k]=NULL;
  *num=k;
  return w;
}
/*==========================================================*/
void *read_mixed_data(INT nrec, INT ni, INT nr, char *the_string)
{
  /* read from a string n records. Each record has */
  /* strings which describe n_i INTs and n_r REALs. Store the output
     in idata[] and rdata[] */
  char **w;
  INT i,iread,count,cni,cnr,k,ki,kr,j,num;
  void *out=(void *)malloc(nrec*(ni*sizeof(INT)+nr*sizeof(REAL)));  
  INT *idata;
  REAL *rdata;
  if(ni>0)
    idata=(INT *)out;
  else
    idata=NULL;
  if(nr>0)
    rdata=(REAL *)(out + nrec*ni*sizeof(INT));
  else
    rdata=NULL;  
  if(idata==NULL && rdata==NULL) return NULL;
  w=splits(the_string," ",&num);
  //  for(k=0;k<num;k++)
  //    {    fprintf(stdout,"\n[%d]:%s",k,w[k]);fflush(stdout);}
  k=0;
  for(count=0;count<nrec;count++){
    if(w[k]==NULL) break;
    if(idata!=NULL) {
      cni=count*ni;
      for(j=0;j<ni;j++){
	iread=sscanf(w[k],"%d",(idata+cni+j));
	//	fprintf(stdout,"\nc=%d;z=%d ",count,idata[cni+j]);fflush(stdout);
	if(iread<0) iread=0;
	free(w[k]);
	k++;
      }
    }
    if(rdata!=NULL) {
      //      fprintf(stdout,"\n[%d]:yyy=%s\n",k,w[k]);fflush(stdout);
      cnr=count*nr;
      for(j=0;j<nr;j++){
	iread=sscanf(w[k],"%lg",(rdata+cnr+j));
	if(iread<0) iread=0;
	free(w[k]);
	k++;
      }
    }
  }
  free(w);
  return out;
}
/*==========================================================*/
void  read_data1(char **clndata,input_grid *g)
{
  // read first the data related to coord systems.
  // correspondence is below: very important
  /* title=clndata[0]; */
  /* dir_grid=clndata[1]; */
  /* dir_vtu=clndata[2]; */
  /* file_grid=clndata[3]; */
  /* file_vtu=clndata[4]; */  
  /* data_coordsystems=clndata[5]; */
  /* data_vertices=clndata[6]; */
  /* data_edges=clndata[7]; */
  /* data_macroelements=clndata[8]; */
  /* data_macrofaces=clndata[9]; */
  /* dimension=clndata[10];		 */
  /* num_coordsystems=clndata[11]; */
  /* num_vertices=clndata[12]; */
  /* num_edges=clndata[13]; */
  /* num_macroelements=clndata[14]; */
  /* num_macrofaces=clndata[15]; */
  /* num_refinements=clndata[16]; */
  /* refinement_type=clndata[17]; */
  /* err_stop_refinement=clndata[18]; */
  /* print_level=clndata[19]; */
  //  
  INT i,iread,count,k,j,num;
  void *mdata;
  INT *idata=NULL;
  /********************* coord_systems*****************/
  mdata=read_mixed_data(g->ncsys,2,g->dim,clndata[5]);
  if(mdata!=NULL){
    idata=(INT *)mdata;
    r2c1(g->ncsys,2,sizeof(INT),idata);// by rows. 
    memcpy(g->ox,(mdata+g->ncsys*2*sizeof(INT)),g->ncsys*g->dim*sizeof(REAL));
    /* for(count=0;count<g->ncsys;count++){ */
    /*   fprintf(stdout,"\nrec=%d (",count); */
    /*   for(j=0;j<g->dim;j++){ */
    /* 	fprintf(stdout,"%f ",g->ox[count*g->dim + j]); */
    /*   } */
    /*   fprintf(stdout,")"); */
    /* } */
    /* fprintf(stdout,"\n"); */
    memcpy(g->syslabels,idata,g->ncsys*sizeof(INT));
    memcpy(g->systypes,(idata+g->ncsys),g->ncsys*sizeof(INT));
    free(mdata); mdata=NULL;
  }
  /***** vertices: same as coordinate systems *****/
  mdata=read_mixed_data(g->nv,2,g->dim,clndata[6]);
  if(mdata!=NULL){
    idata=(INT *)mdata;
    /* for(count=0;count<g->nv;count++){ */
    /*   fprintf(stdout,"\nrec=%d (",count); */
    /*   for(j=0;j<2;j++){ */
    /* 	fprintf(stdout,"%d ",idata[count*2 + j]); */
    /*   } */
    /*   fprintf(stdout,")"); */
    /* } */
    /* fprintf(stdout,"\n"); */
    r2c1(g->nv,2,sizeof(INT),idata);// vertex labels and coord systems by rows.
    memcpy(g->xv,(mdata+g->nv*2*sizeof(INT)),g->nv*g->dim*sizeof(REAL));
    /* for(count=0;count<g->nv;count++){ */
    /*   fprintf(stdout,"\nrec=%d (",count); */
    /*   for(j=0;j<g->dim;j++){ */
    /* 	fprintf(stdout,"%f ",g->xv[count*g->dim + j]); */
    /*   } */
    /*   fprintf(stdout,")"); */
    /* } */
    /* fprintf(stdout,"\n"); */
    memcpy(g->labelsv,idata,g->nv*sizeof(INT));
    memcpy(g->csysv,(idata+g->nv),g->nv*sizeof(INT));
    memset(g->bcodesv,0,g->nv*sizeof(INT));
    free(mdata); mdata=NULL;
    /* fprintf(stdout,"\n"); */
    /* for(j=0;j<g->nv;j++){ */
    /*   fprintf(stdout,"\nl=%d;t=%d;",g->labelsv[j],g->csysv[j]); */
    /* } */
    /* fprintf(stdout,"\n");fflush(stdout); */
  }
  /* convert the degree coordinates to radian coordinates */
  for(count=0;count<g->nv;count++){
    if(g->systypes[g->csysv[count]]==1){
      for(j=1;j<g->dim;j++){
  	//	fprintf(stdout,"\n(%d%d) x=%f",count,j,g->xv[count*g->dim + j]);
  	g->xv[count*g->dim + j]*=(((REAL )PI)/180.);
      }
    }
  }
  /* /\***** edges *****\/ */
  mdata=read_mixed_data(g->ne,3,0,clndata[7]);
  if(mdata!=NULL){
    memcpy(g->seg,mdata,3*g->ne*sizeof(INT));
    free(mdata); mdata=NULL;
  }
  INT ne = 0,iri,ici,ndd;
  // no selfedges
  for(i=0;i<g->ne;i++){
    iri=g->seg[3*i];ici=g->seg[3*i+1];ndd=g->seg[3*i+2];
    if(iri==ici) continue;
    if(iri<ici){
      g->seg[3*ne]=iri;
      g->seg[3*ne+1]=ici;
    } else {
      g->seg[3*ne]=ici;
      g->seg[3*ne+1]=iri;
    }
    g->seg[3*ne+2]=ndd;
    ne++;
  }
  if(ne<g->ne)
    g->seg=realloc(g->seg,3*ne*sizeof(INT));
  g->ne=ne;
  /* end edges */
  INT nvcube=(1<<g->dim);
  mdata=read_mixed_data(g->nel,(nvcube+1),0,clndata[8]);
  if(mdata!=NULL){
    memcpy(g->mnodes,mdata,g->nel*(nvcube+1)*sizeof(INT));
    free(mdata); mdata=NULL;
  }
  INT nvface=(1<<(g->dim-1));
  mdata=read_mixed_data(g->nf,(nvface+1),0,clndata[9]);
  if(mdata!=NULL){
    memcpy(g->mfaces,mdata,g->nf*(nvface+1)*sizeof(INT));
    free(mdata); mdata=NULL;
  }/*  else { */
  /*   free(g->mfaces); g->mfaces=NULL; */
  /* } */
  input_grid_print(g);
  return;
}
/********************************************************************/
void get_out(const char *pattern, size_t le)
{
  /* prints a string cutting it at the closest blank space <= le*/
  int i;
  fprintf(stderr, "\n\n\n           *** ERROR::::   \n     UNBALANCED \"{}\" near or before \"");
  for (i=0;i<(le-1);++i)
    fprintf(stderr, "%c",*(pattern+i));
  fprintf(stderr, "\"\n");  
  exit(12);
}
/***********************************************************************/
char *make_string_from_file(FILE *the_file, size_t *length_string)
{
  char *everything;
  char ch, ch_next;
  int count=0,i,j,flag;
  while(feof(the_file)==0) 
    {
      ch = fgetc(the_file);
      if(ch) count++;
    }
  count--;
  i = count*sizeof(char);
  everything = (char *)malloc(i);
  rewind(the_file);
  //  fprintf(stderr,"\nNumber of characters in the file %i\n", count);
  i = 0;
  j = 0;
  flag = 1;
  while(j < count) {
    if( flag ) {
      ch = fgetc(the_file);
      ++j;
    }
    if( ch == '%' ){
      /* delete to the end of line or to the end of file, whatever it is */
      do{
	ch = fgetc(the_file);
	++j;
      }  while(ch != '\n' && ch != EOF);
      if( ch == '\n' ) 
	flag=1;
      else 
	break; 
    } else {
      if(ch == '\n' || ch == '\t') ch = ' ';
      if(ch == ' ' || ch == '{' || ch ==  '}'){
	do{
	  ch_next = fgetc(the_file);
	  if(ch_next == '\n' || ch_next == '\t') ch_next = ' ';
	  ++j;
	}  while(ch_next == ' ');
	if(ch == ' ' && (ch_next ==  '{' || ch_next ==  '}')){
	  ch = ch_next;
	  flag=0;
	} else {
	  /*		  printf(" %i ",ch);*/
	  *(everything+i) = ch;
	  ++i;
	  ch = ch_next;
	  flag=0;
	}
      } else if( ch != EOF ) {
	/*	      printf(" %i ",ch);*/
	*(everything+i) = ch;
	++i;
	flag=1;
      }
    }
    /*           printf("\n i,j, %i %i  %i\n",i,j,j-count);*/
  }
  if(i) 
    if (*(everything+i-1)== ' ') i--;
  /* end of string */
  *(everything+i) = '\0';
  *length_string = i;
  everything = (char *)realloc(everything,i*sizeof(char));
  /* for(j=0;j<i;j++){ */
  /*   if(!everything[j]) continue; */
  /*   everything[j]=toupper(tolower(everything[j])); */
  /* } */
  /* fprintf(stderr,"\nNumber of characters in the supressed string %li\n",strlen(everything)); */
  /* fprintf(stdout,"\nString=%s\n",everything); */
  return everything;
}
/********************************************************************/
char *get_substring(const char *pattern,		\
		    size_t *length_substring,	\
		    char *the_string)
{
  /* 
     get substring from a string matching a pattern; it can probably
     be done with strtok() but let us use this for now. 
  */  
  size_t le;
  INT i;
  char *found, *wrk;
  //  char ch;
  le = strlen(pattern);
  found = (char *) strstr(the_string,pattern);
  if(found != NULL){
    //    found = &(*(found + le));
    found += le;
    wrk = (char *) strstr(found,"}");
    if(wrk == NULL ){ get_out(pattern,le);}
    *length_substring=strlen(found)-strlen(wrk);
    *wrk = '\t';
    wrk = (char *) strstr(found,"\t");
    i = strlen(found)-strlen(wrk);
    if(i != *length_substring ){ get_out(pattern,le); }
  }else{
    fprintf(stderr, "\n\n\n *** WARNING: \"" );
    for (i=0;i<le-1;++i)
      fprintf(stderr, "%c",*(pattern+i));
    fprintf(stderr, "\" has not been found in the input file\n\n");
    found = (char *)malloc(sizeof(char));
    *length_substring=0;
  }
  return found;
}
/********************************************************************/
input_grid *parse_input_grid(const char *input_file_grid)
{
  INT iread,numel_data,k;
  FILE *the_file;
  char *everything;
  size_t length_string=0;
  char **indata;
  char **clndata;
  size_t *lengths;
  the_file = fopen(input_file_grid,"r");
  everything = make_string_from_file(the_file, &length_string);
  fclose(the_file);
  //  fprintf(stdout,"\n%s\n",everything);
  indata=input_strings(&numel_data);
  clndata=malloc(numel_data*sizeof(char *));
  lengths=calloc(numel_data,sizeof(size_t));
  //
  /* get all substrings */
  for(k=0;k<numel_data;k++){
    clndata[k] = get_substring(indata[k],(lengths+k), everything);
  }
  for(k=0;k<numel_data;k++)
    clndata[k][lengths[k]] = '\0';
  /* initialize */
  input_grid *g=malloc(1*sizeof(input_grid));    
  /* ... PARSE ... */
  /* strings */
  g->title=strdup(clndata[0]);
  g->dgrid=strdup(clndata[1]);
  g->fgrid=strdup(clndata[2]);
  g->dvtu=strdup(clndata[3]);
  g->fvtu=strdup(clndata[4]);
  /*INTEGERS*/
  iread=sscanf( clndata[10],"%d",&g->dim); // the dimension of the problem.
  iread=sscanf( clndata[11],"%d",&g->ncsys);//
  iread=sscanf( clndata[12],"%d",&g->nv);//
  iread=sscanf( clndata[13],"%d",&g->ne);//
  iread=sscanf( clndata[14],"%d",&g->nel);//
  iread=sscanf(clndata[15],"%d",&g->nf);//  
  iread=sscanf(clndata[16],"%d",&g->nref);//
  iread=sscanf(clndata[17],"%d",&g->ref_type);//
  iread=sscanf(clndata[18],"%lg",&g->err_stop);//
  iread=sscanf(clndata[19],"%hd",&g->print_level);//
  if(iread<0) iread=0;
  input_grid_arrays(g);
  read_data1(clndata,g);  
  /*FREE*/
  if(everything)free(everything);
  for(k=0;k<numel_data;k++){
    free(indata[k]);
  }
  free(indata);
  free(clndata);
  return g;
}
/********************************************************************/