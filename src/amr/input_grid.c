/*! \file src/amr/input_grid.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note containing all essential routines for input for the mesh
 *  generation mesh refinement
 *
 */
#include "hazmath.h"
/**********************************************************************/
void input_grid_free(input_grid *g)
{
  if(g->title) free(g->title);
  if(g->dgrid) free(g->dgrid);
  if(g->fgrid) free(g->fgrid);
  if(g->dvtu) free(g->dvtu);
  if(g->fvtu) free(g->fvtu);
  if(g->ox) free(g->ox);
  if(g->systypes) free(g->systypes);
  if(g->syslabels) free(g->syslabels);
  if(g->labels) free(g->labels);
  if(g->bcodes) free(g->bcodes);
  if(g->x) free(g->x);
  if(g) free(g);
  icoo_free(g->seg);
  return;
}
/**********************************************************************/
void input_grid_print(input_grid *g)
{
  // prints input grid. 
  INT i,j,dim=g->dim;
  fprintf(stdout,"\n\nTITLE: %s",g->title);
  fprintf(stdout,"\ndimension=%d",g->dim);
  fprintf(stdout,"\nprint_level=%d",g->print_level);
  fprintf(stdout,"\ndir_grid=%s",g->dgrid);
  fprintf(stdout,"\ndir_vtu=%s",g->dvtu);
  fprintf(stdout,"\nfile_grid=%s",g->fgrid);
  fprintf(stdout,"\nfile_vtu=%s",g->fvtu);
  fprintf(stdout,"\n\nnum_coordsystems=%d",g->ncsys);
  for(i=0;i<g->ncsys;i++){
    fprintf(stdout,"\nlabel=%d,type=%d, origin(",g->syslabels[i],g->systypes[i]);
    for(j=0;j<g->dim;j++) fprintf(stdout," %6.2f ",g->ox[i*dim+j]);
    fprintf(stdout,")");
  }
  fprintf(stdout,"\n\nnum_vertices=%d",g->nv);
  for(i=0;i<g->nv;i++){
    fprintf(stdout,"\nvertex=%d, coord_system=%d, bcode=%d, coords(",i,g->labels[i],g->bcodes[i]);
    if(g->systypes[g->labels[i]]==1){
      fprintf(stdout," %6.2f ",g->x[i*dim]);
      for(j=1;j<g->dim;j++)	
	fprintf(stdout," %6.2f ",(g->x[i*dim+j])/((REAL )PI)*180.);
    }else{
      for(j=0;j<g->dim;j++) fprintf(stdout," %6.2f ",g->x[i*dim+j]);
    }
    fprintf(stdout,")");
  }
  fprintf(stdout,"\n\nnum_edges=%d\n",g->ne);
  /* INT iaa,iab; */
  /* for(i=0;i<(g->seg->row);i++){ */
  /*   iaa=g->seg->IA[i]; */
  /*   iab=g->seg->IA[i+1]; */
  /*   for(j=iaa;j<iab;j++){ */
  /*     fprintf(stdout,"\nedge=(%d,%d) div=%d",i,g->seg->JA[j],g->seg->val[j]); */
  /*   } */
  /* } */
  /* fprintf(stdout,"\n\n"); */
  for(i=0;i<(g->seg)->nnz;i++){
    fprintf(stdout,"\nedge=(%d,%d) div=%d",g->seg->rowind[i],g->seg->colind[i],g->seg->val[i]);
  }
  fprintf(stdout,"\n\n");
  
  return;
}
/*---------------------------------------------------------------------*/
void coo2csr(INT nrow,INT ncol,INT nnz,					\
	     INT *row_idx,INT *col_idx, void *aval,			\
	     INT *ia,INT *ja, void *bval,				\
	     size_t elsize)
{
  // converts (i,j,value) to (ia,ja,bval) for matrices matrices whose
  // elements have elsize bytes
  //
  // get dimensions of the matrix
  //
  const INT m=nrow, n=ncol;
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
    memcpy((bval+jind*elsize),(aval+i*elsize),elsize);
    ind[iind] = ++jind;
  }    
  if (ind) free(ind);    
  return;
}
/*---------------------------------------------------------------------*/
char **splits(char *s, const char *d, INT *num)
{
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
void read_data(char *data_coordsystems,		\
	       char *data_vertices,		\
	       char *data_edges,		\
	       input_grid *g)
{
  // read first the data related to coord systems. 
  char **w;
  INT i,iread,count,k,j,ksys,num;
  w=splits(data_coordsystems," ",&num);
  k=0;
  for(count=0;count<g->ncsys;count++){
    if(w[k]==NULL) break;
    iread=sscanf(w[k],"%d",&g->syslabels[count]);
    free(w[k]);
    k++;
    for(j=0;j<g->dim;j++){
      iread=sscanf(w[k],"%lg",(g->ox +count*g->dim+j));
      free(w[k]);
      k++;
    }
    iread=sscanf(w[k],"%d",&g->systypes[count]);
    free(w[k]);
    k++;
  }
  if(w) free(w);
  /***** vertices *****/
  w=splits(data_vertices," ",&num);
  k=0;
  for(count=0;count<g->nv;count++){
    if(w[k]==NULL) break;
    for(j=0;j<g->dim;j++){
      iread=sscanf(w[k],"%lg",(g->x + count*g->dim+j));
      if(w[k]) free(w[k]);
      k++;
    }
    iread=sscanf(w[k],"%d",&g->labels[count]);
    if(w[k]) free(w[k]);
    k++;
    iread=sscanf(w[k],"%d",&g->bcodes[count]);
    if(w[k]) free(w[k]);
    k++;
  }
  if(w) free(w);
  for(count=0;count<g->nv;count++){
    if(g->systypes[g->labels[count]]==1){
      for(j=1;j<g->dim;j++){
	//	fprintf(stdout,"\n(%d%d) x=%f",count,j,g->x[count*g->dim + j]);
	g->x[count*g->dim + j]*=(((REAL )PI)/180.);
      }
    }
  }
  /***** edges *****/
  w=splits(data_edges," ",&num);
  INT *ir=(INT *)calloc(g->ne,sizeof(INT));
  INT *ic=(INT *)calloc(g->ne,sizeof(INT));
  INT *ndiv=(INT *)calloc(g->ne,sizeof(INT));
  k=0;
  for(count=0;count<g->ne;count++){
    if(w[k]==NULL) break;
    iread=sscanf(w[k],"%d",ir+count);
    if(w[k]) free(w[k]);
    k++;
    iread=sscanf(w[k],"%d",ic+count);
    if(w[k]) free(w[k]);
    k++;
    iread=sscanf(w[k],"%d",ndiv+count);
    if(w[k]) free(w[k]);
    k++;
  }
  INT ne = 0,iri,ici;
  // no selfedges
  for(i=0;i<g->ne;i++){
    iri=ir[i];ici=ic[i];
    if(iri==ici) continue;
    if(iri<ici){
      ir[ne]=iri;
      ic[ne]=ici;
    } else {
      ir[ne]=ici;
      ic[ne]=iri;
    }
    ndiv[ne]=ndiv[i];
    ne++;
  }
  if(ne<g->ne){
    ir=realloc(ir,ne*sizeof(INT));
    ic=realloc(ic,ne*sizeof(INT));
    ndiv=realloc(ndiv,ne*sizeof(INT));
    g->ne=ne;
  }
  /*uncomment the block below if CSR is used later. */
  // transpose:
  /* if(ne<g->ne){ */
  /*   ir=realloc(ir,2*ne*sizeof(INT)); */
  /*   ic=realloc(ic,2*ne*sizeof(INT)); */
  /*   ndiv=realloc(ndiv,2*ne*sizeof(INT)); */
  /*   g->ne=ne;  */
  /* } */
  /* for(i=0;i<g->ne;i++){ */
  /*   ir[g->ne+i]=ic[i]; */
  /*   ic[g->ne+i]=ir[i]; */
  /*   ndiv[g->ne+i]=ndiv[i]; */
  /*   ne++; */
  /* } */
  g->ne=ne;
  g->seg=malloc(sizeof(iCOOmat));
  (g->seg)->row=(g->seg)->col=g->nv;
  (g->seg)->nnz=g->ne;
  (g->seg)->rowind=ir;
  (g->seg)->colind=ic;
  (g->seg)->val=ndiv;
  /* g->seg=malloc(sizeof(iCSRmat)); */
  /* (g->seg)->row=(g->seg)->col=g->nv; */
  /* (g->seg)->nnz=g->ne; */
  /* (g->seg)->IA=calloc((g->seg)->row+1,sizeof(INT)); */
  /* (g->seg)->JA=calloc((g->seg)->nnz,sizeof(INT)); */
  /* (g->seg)->val=calloc((g->seg)->nnz,sizeof(INT)); */
  /* coo2csr((g->seg)->row,			\ */
  /* 	  (g->seg)->col,			\ */
  /* 	  (g->seg)->nnz,			\ */
  /* 	  ir,ic,ndiv,				\ */
  /* 	  (g->seg)->IA,				\ */
  /* 	  (g->seg)->JA,				\ */
  /* 	  (g->seg)->val,			\ */
  /* 	  sizeof(INT)); */
  /* fix the ndiv; if in a edge */
  /* for(i=0;i<g->seg->nnz;i++){ */
  /*   fprintf(stdout,"\ndiff=%d; (%d,%d)",ic[i]-ir[i],ir[i],ic[i]); */
  /* } */
  /* fprintf(stdout,"\n"); */
  if(w) free(w);
  return;
}
/********************************************************************/
void get_out(char *pattern, size_t le)
{
  /* prints a string cutting it at the closest blank space <= le*/
  int i;
  fprintf(stderr, "\n\n\n           *** ERROR::::   \n     UNBALANCED \"{}\" near:::   ");
  for (i=0;i<le;++i)
    fprintf(stderr, "%c",*(pattern+i));
  fprintf(stderr, "...}\n");
  exit(128);
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
  /*  fprintf(stderr,"\nNumber of characters in the supressed string
      %i\n", i); */
  /* fprintf(stderr,"\nNumber of characters in the supressed string
     %li\n",strlen(everything)); */
  return everything;
}
/********************************************************************/
char *get_substring(char *pattern,		\
		    size_t *length_substring,	\
		    char *the_string)
{
  /* 
     get substring from a string matching a pattern; it can probably
     be done with strtok()
  */
  
  size_t le;
  INT i;
  char *found, *wrk;
  char ch;
  le = strlen(pattern);
  found = (char *) strstr(the_string,pattern);
  if(found != NULL){
    found = &(*(found + le));
    wrk = (char *) strstr(found,"}");
    if(wrk == NULL ){ get_out(pattern,le);}
    *length_substring=strlen(found)-strlen(wrk);
    *wrk = '\t';
    wrk = (char *) strstr(found,"\t");
    i = strlen(found)-strlen(wrk);
    if(i != *length_substring ){ get_out(pattern,le); }
  }else{
    fprintf(stderr, "\n\n\n *** WARNING:::: " );
    for (i=0;i<le;++i)
      fprintf(stderr, "%c",*(pattern+i));
    fprintf(stderr, "...} has not been found in the input file\n\n");
    found = (char *)malloc(sizeof(char));
    *length_substring=0;
  }
  return found;
}
input_grid *parse_input_grid(const char *input_file_grid)
{
  INT iread,i,j;
  FILE *the_file;
  char *everything;
  size_t length_string=0;
  char *title=NULL, 
    *dimension=NULL, 
    *print_level=NULL, 
    *dir_grid=NULL,
    *dir_vtu=NULL,
    *file_grid=NULL,
    *file_vtu=NULL,
    *num_coordsystems=NULL, 
    *data_coordsystems=NULL,
    *num_vertices=NULL,
    *data_vertices=NULL,
    *num_edges=NULL,
    *data_edges=NULL;
  size_t length_info_file;
  size_t length_title=0, 
    length_dimension=0, 
    length_print_level=0, 
    length_dir_grid=0,
    length_dir_vtu=0,
    length_file_grid=0,
    length_file_vtu=0,
    length_num_coordsystems=0, 
    length_data_coordsystems=0,
    length_num_vertices=0,
    length_data_vertices=0,
    length_num_edges=0,
    length_data_edges=0;
  the_file = fopen(input_file_grid,"r");
  everything = make_string_from_file(the_file, &length_string);
  fclose(the_file);
  //  fprintf(stdout,"\n%s\n",everything);
  /* get all substrings */
  title  = get_substring("title{",&length_title, everything);
  dimension  = get_substring("dimension{",&length_dimension, everything);
  print_level        = get_substring("print_level{",&length_print_level, everything);
  dir_grid      = get_substring("dir_grid{",   &length_dir_grid, everything);
  dir_vtu       = get_substring("dir_vtu{",&length_dir_vtu, everything);
  file_grid      = get_substring("file_grid{",&length_file_grid, everything);
  file_vtu         = get_substring("file_vtu{",&length_file_vtu, everything);
  num_coordsystems       = get_substring("num_coordsystems{",&length_num_coordsystems, everything);
  data_coordsystems       = get_substring("data_coordsystems{",&length_data_coordsystems, everything);
  num_vertices        = get_substring("num_vertices{",&length_num_vertices, everything);
  data_vertices     = get_substring("data_vertices{",&length_data_vertices, everything);
  num_edges  = get_substring("num_edges{",&length_num_edges, everything);
  data_edges = get_substring("data_edges{",&length_data_edges, everything);
  /* ... */
  //  fprintf(stdout,"\n***^^^TITLE:%s\n\n",title);fflush(stdout);
  *(title + length_title) = '\0';
  *(dimension + length_dimension) = '\0';
  *(print_level + length_print_level) = '\0';
  *(dir_grid + length_dir_grid) = '\0';
  *(dir_vtu + length_dir_vtu) = '\0';
  *(file_grid + length_file_grid) = '\0';
  *(file_vtu + length_file_vtu) = '\0';
  *(num_coordsystems + length_num_coordsystems) = '\0';
  *(data_coordsystems + length_data_coordsystems) = '\0';
  *(num_vertices + length_num_vertices) = '\0';
  *(data_vertices + length_data_vertices) = '\0';
  *(num_edges + length_num_edges) = '\0';
  *(data_edges + length_data_edges) = '\0';
  /* ... PARSE ... */
  input_grid *g=malloc(1*sizeof(input_grid));
  iread=sscanf(dimension,"%d",&g->dim); // the dimension of the problem.
  iread=sscanf(print_level,"%hd",&g->print_level);//
  g->title=(char *)calloc(strlen(title),sizeof(char));
  g->dgrid=(char *)calloc(strlen(dir_grid),sizeof(char));
  g->fgrid=(char *)calloc(strlen(file_grid),sizeof(char));
  g->dvtu=(char *)calloc(strlen(dir_vtu),sizeof(char));
  g->fvtu=(char *)calloc(strlen(file_vtu),sizeof(char));
  strncpy(g->title,title,strlen(title));  /**< grid dir name */
  strncpy(g->dgrid,dir_grid,strlen(dir_grid));  /**< grid dir name */
  strncpy(g->fgrid,file_grid,strlen(file_grid));  /**< grid file name */
  strncpy(g->dvtu,dir_vtu,strlen(dir_vtu));  /**< vtu grid */
  strncpy(g->fvtu,file_vtu,strlen(file_vtu));  /**< vtu file */
  // integers;
  iread=sscanf(num_coordsystems,"%d",&g->ncsys);//
  iread=sscanf(num_vertices,"%d",&g->nv);//
  iread=sscanf(num_edges,"%d",&g->ne);//
  //mixed
  g->ox=(REAL *)calloc(g->dim*g->ncsys,sizeof(REAL));
  g->systypes=(INT *)calloc(g->ncsys,sizeof(INT)); 
  g->syslabels=(INT *)calloc(g->ncsys,sizeof(INT)); 
  g->labels=(INT *)calloc(g->nv,sizeof(INT)); 
  g->bcodes=(INT *)calloc(g->nv,sizeof(INT)); 
  g->x=(REAL *)calloc(g->dim*g->nv,sizeof(REAL)); 
  /*  iCSRmat *graph;// icsrmat thing for the macroelement graph. */
  /* read arrays and convert data */
  read_data(data_coordsystems,data_vertices, data_edges,g);
  if(everything)free(everything);
  return g;
}
/********************************************************************/
