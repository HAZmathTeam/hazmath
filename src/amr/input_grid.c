/*! \file src/amr/input_grid.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note containing all essential routines for input for the mesh
 *  generation and mesh refinement
 *
 *  \note: modified by ltz on 20190327
 *  \note: modified by ltz on 20190728
 *  \note: modified by ltz on 20210813
 *
 */
/*********************************************************************/
#include "hazmath.h"
/*********************************************************************/
/*!
 * \fn char **input_strings(INT *nall_out)
 *
 * \brief take as input the strings defining the initial grid and on
 *        output has an array with all such strings which are used in
 *        parsing the input file for grid generation.
 *
 * \param nall_out=how many strings are in the INPUT_GRID_DATA;
 *
 * \return
 *
 * \note
 */
static char **input_strings(INT *nall_out)
{
  INT k,nall;
  const char *indata_const[]={			\
			      "title{",		\
			      "file_grid{",	\
			      "file_vtu{",	\
			      "data_coordsystems{",	\
			      "data_vertices{",		\
			      "data_edges{",		\
			      "data_macroelements{",	\
			      "data_macrofaces{",	\
			      "dimension{",		\
			      "num_coordsystems{",	\
			      "num_vertices{",		\
			      "num_edges{",		\
			      "num_macroelements{",	\
			      "num_macrofaces{",	\
			      "num_refinements{",	\
			      "refinement_type{",	\
			      "amr_marking_type{",	\
			      "err_stop_refinement{",	\
			      "print_level{",		\
			      "num_refine_points{",	\
			      "data_refine_points{",	\
			      "\0" };
  
  nall=0;
  while(strlen(indata_const[nall]))nall++;
  //  fprintf(stdout,"\nthe nall is: %d\n",nall);
  char **indata=malloc(nall*sizeof(char *));
  for(k=0;k<nall;k++)
    indata[k]=strndup(indata_const[k],strlen(indata_const[k]));
  *nall_out=nall;
  return indata;
}
/*********************************************************************/
/*!
 * \fn void input_grid_arrays(input_grid *g)
 *
 * \brief Allocate the arrays used in input grid. Some of the values
 *        have abs() because they could be negative if their values
 *        were missing or were not read correctly from the input file.
 *
 * \param initial macroelement grid g.
 *
 * \return
 *
 * \note
 *
 */
void input_grid_arrays(input_grid *g)
{
  INT nvcube=(1 << g->dim),nvface=(1 << (g->dim-1));
  g->ox=(REAL *)calloc(g->dim*abs(g->ncsys),sizeof(REAL));
  g->systypes=(INT *)calloc(g->ncsys,sizeof(INT));
  g->syslabels=(INT *)calloc(g->ncsys,sizeof(INT));
  g->csysv=(INT *)calloc(g->nv,sizeof(INT));
  g->labelsv=(INT *)calloc(g->nv,sizeof(INT));
  g->bcodesv=(INT *)calloc(g->nv,sizeof(INT));
  g->xv=(REAL *)calloc(g->dim*g->nv,sizeof(REAL));
  g->xe=(REAL *)calloc(g->dim*abs(g->ne),sizeof(REAL));
  g->seg=(INT *)calloc(3*abs(g->ne),sizeof(INT));
  g->mnodes=(INT *)calloc(g->nel*(nvcube+1),sizeof(INT));
  g->mfaces=(INT *)calloc(abs(g->nf)*(nvface+1),sizeof(INT));
  if(g->num_refine_points){
    g->data_refine_points=(REAL *)calloc(g->num_refine_points*g->dim,sizeof(REAL));
    memset(g->data_refine_points,0,g->num_refine_points*g->dim*sizeof(REAL));
  } else {
    g->data_refine_points=NULL;
  }
  //init
  memset(g->ox,0,g->dim*abs(g->ncsys)*sizeof(REAL));
  memset(g->systypes,0,g->ncsys*sizeof(INT));
  memset(g->syslabels,0,g->ncsys*sizeof(INT));
  memset(g->csysv,0,g->nv*sizeof(INT));
  memset(g->labelsv,0,g->nv*sizeof(INT));
  memset(g->bcodesv,0,g->nv*sizeof(INT));
  memset(g->xv,0,g->dim*g->nv*sizeof(REAL));
  memset(g->xe,0,g->dim*abs(g->ne)*sizeof(REAL));
  memset(g->seg,0,3*abs(g->ne)*sizeof(INT));
  memset(g->mnodes,0,g->nel*(nvcube+1)*sizeof(INT));
  memset(g->mfaces,0,abs(g->nf)*(nvface+1)*sizeof(INT));
  return;
}
/******************************************************************/
/*!
 * \fn void input_grid_free(input_grid *g)
 *
 * \brief frees all arrays in a structure input_grid *
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void input_grid_free(input_grid *g)
{
  if(g->title) free(g->title);
  if(g->fgrid) free(g->fgrid);
  if(g->fvtu) free(g->fvtu);
  if(g->ox) free(g->ox);
  if(g->systypes) free(g->systypes);
  if(g->syslabels) free(g->syslabels);
  if(g->labelsv) free(g->labelsv);
  if(g->csysv) free(g->csysv);
  if(g->bcodesv) free(g->bcodesv);
  if(g->xv) free(g->xv);
  if(g->xe) free(g->xe);
  if(g->seg) free(g->seg);
  if(g->mnodes) free(g->mnodes);
  if(g->mfaces) free(g->mfaces);
  if(g->data_refine_points) free(g->data_refine_points);
  if(g) free(g);
  return;
}
/**********************************************************************/
/*!
 * \fn void input_grid_print(input_grid *g)
 *
 * \brief prints what was read from the input grid file.
 *
 * \param input grid of macroelements g (read usually from a file)
 *
 * \return
 *
 * \note
 *
 */
void input_grid_print(input_grid *g)
{
  INT i,j,dim=g->dim;
  fprintf(stdout,"\n\nTITLE: %s",g->title);fflush(stdout);
  fprintf(stdout,"\ndimension=%d",g->dim);fflush(stdout);
  fprintf(stdout,"\nfile_grid=%s",g->fgrid);fflush(stdout);
  fprintf(stdout,"\nfile_vtu=%s",g->fvtu);fflush(stdout);
  fprintf(stdout,"\nprint_level=%d",g->print_level);fflush(stdout);
  fprintf(stdout,"\nnum_refinements=%d",g->nref);fflush(stdout);
  fprintf(stdout,"\nrefinement_type=%d",g->ref_type);fflush(stdout);
  fprintf(stdout,"\namr_marking_type=%d",g->mark_type);fflush(stdout);
  fprintf(stdout,"\nerr_stop_amr=%.3g",g->err_stop);fflush(stdout);
  /*ARRAYS*/
  fprintf(stdout,"\n\nnum_coordsystems=%d",g->ncsys);
  for(i=0;i<g->ncsys;i++){
    fprintf(stdout,"\nlabel=%d,type=%d, origin(",g->syslabels[i],g->systypes[i]);fflush(stdout);
    for(j=0;j<g->dim;j++) fprintf(stdout," %6.2f ",g->ox[i*dim+j]);
    fflush(stdout);
    fprintf(stdout,")");
  }
  fprintf(stdout,"\n\nnum_vertices=%d\n",g->nv);fflush(stdout);
  for(i=0;i<g->nv;i++){
    fprintf(stdout,"\nvertex=%d, coord_system=%d, bcode=%d, coords(",i,g->csysv[i],g->bcodesv[i]);fflush(stdout);
    if(g->systypes[g->csysv[i]]==1){
      fprintf(stdout," %6.2f ",g->xv[i*dim]);
      for(j=1;j<g->dim;j++) fprintf(stdout," %6.2f ",(g->xv[i*dim+j])/((REAL )PI)*180.);
      fflush(stdout);
    }else{
      for(j=0;j<g->dim;j++) fprintf(stdout," %6.2f ",g->xv[i*dim+j]);
      fflush(stdout);
    }
    fprintf(stdout,")");
  }
  fprintf(stdout,"\n\nnum_edges=%d\n",g->ne);
  for(i=0;i<g->ne;i++)
    fprintf(stdout,"\nedge=(%d,%d) div=%d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
  fflush(stdout);
  INT nvcube=(1<<g->dim),nvcube1=nvcube+1;
  fprintf(stdout,"\n\nnum_macroelements=%d\n",g->nel);fflush(stdout);
  for(i=0;i<g->nel;i++){
    fprintf(stdout,"\nmacroel=%d; code=%d; vertices=(",i,g->mnodes[i*nvcube1+nvcube]);
    for(j=0;j<nvcube;j++)
      fprintf(stdout,"%d ",g->mnodes[nvcube1*i+j]);
    fprintf(stdout,")");fflush(stdout);
  }
  INT nvface=(1<<(g->dim-1)),nvface1=nvface+1;
  fprintf(stdout,"\n\nnum_faces=%d\n",g->nf);
  for(i=0;i<g->nf;i++){
    fprintf(stdout,"\nmacroface=%d; code=%d; vertices=(",i,g->mfaces[i*nvface1+nvface]);
    fflush(stdout);
    for(j=0;j<nvface;j++)
      fprintf(stdout,"%d ",g->mfaces[nvface1*i+j]);
    fprintf(stdout,")");
    fflush(stdout);
  }

  // print refine points
  fprintf(stdout,"\n\nnum_refine_points=%d\n",g->num_refine_points);fflush(stdout);
  for(i=0;i<g->num_refine_points;i++){
    fprintf(stdout,"\n refine point=%d, coord_system=%d, coords(",i,g->csysv[i]);fflush(stdout);
    if(!(g->num_refine_points)) continue;
    if(g->systypes[g->csysv[i]]==1){
      fprintf(stdout," %6.2f ",g->data_refine_points[i*dim]);
      for(j=1;j<g->dim;j++) fprintf(stdout," %6.2f ",(g->data_refine_points[i*dim+j])/((REAL )PI)*180.);
      fflush(stdout);
    }else{
      for(j=0;j<g->dim;j++) fprintf(stdout," %6.2f ",g->data_refine_points[i*dim+j]);
      fflush(stdout);
    }
    fprintf(stdout,")");
  }

  fprintf(stdout,"\n\n");fflush(stdout);
  return;
}
/**********************************************************************/
/*!
 * \fn void input_grid_example_file(FILE *fp,input_grid *g)
 *
 * \brief For a given input grid g, it outputs to stderr the file
 *        which could have generated such macroelement grid. It is
 *        used to print an example of an input file for the partition
 *        of unit cube in 3D when errors are found in the user input
 *        (basically fail-safe strategy)
 *
 * \param input grid g.
 *
 * \return
 *
 * \note
 *
 */
void input_grid_example_file(FILE *fp,input_grid *g)
{
  INT i,j,dim=g->dim;
  fprintf(fp,"\n%s\n%s","%%%%%% An *example of an input file* follows.",
	  "For a lattice grid on the unit cube in 3D:\nCopy and paste in a file the text between %=== and %---):");
  fprintf(fp,"\n%%%s%s\n%%%s",					\
	  "===================================================",	\
	  "===================================================",	\
	  "%%%%BEGIN(input grid file)");
  fprintf(fp,"\ntitle{%s}",g->title);
  fprintf(fp,"\ndimension{%d}",g->dim);
  fprintf(fp,"\nfile_grid{%s}",g->fgrid);
  fprintf(fp,"\nfile_vtu{%s}",g->fvtu);
  fprintf(fp,"\nprint_level{%d}",g->print_level);
  fprintf(fp,"\nnum_refinements{%d}",g->nref);
  fprintf(fp,"\nrefinement_type{%d}",g->ref_type);
  fprintf(fp,"\namr_marking_type{%d}",g->mark_type);
  fprintf(fp,"\nerr_stop_amr{%.3g}",g->err_stop);
  //
  /*COORDSYSTEMS*/
  fprintf(fp,"\nnum_coordsystems{%d}",g->ncsys);
  fprintf(fp,"\ndata_coordsystems{");
  fprintf(fp,"%d %d ",g->syslabels[0],g->systypes[0]);
  for(j=0;j<g->dim;j++) fprintf(fp," %.4e ",g->ox[j]);
  for(i=1;i<g->ncsys;i++){
    fprintf(fp,"\n%d %d ",g->syslabels[i],g->systypes[i]);
    for(j=0;j<g->dim;j++) fprintf(fp," %.4e ",g->ox[i*dim+j]);
  }
  fprintf(fp,"}\n");
  /*VERTICES */
  fprintf(fp,"\nnum_vertices{%d}",g->nv);
  fprintf(fp,"\ndata_vertices{");
  fprintf(fp,"%d %d ",0,g->csysv[0]);
  if(g->systypes[g->csysv[0]]==1){
    fprintf(fp," %.4e ",g->xv[0]);
    for(j=1;j<g->dim;j++) fprintf(fp," %.4e ",g->xv[j]/(PI)*180e00);
  }else{
    for(j=0;j<g->dim;j++) fprintf(fp," %.4e ",g->xv[j]);
  }
  for(i=1;i<g->nv;i++){
    fprintf(fp,"\n%d %d ",i,g->csysv[i]);
    if(g->systypes[g->csysv[i]]==1){
      fprintf(fp," %.4e ",g->xv[i*dim]);
      for(j=1;j<g->dim;j++) fprintf(fp," %.4e ",g->xv[i*dim+j]/(PI)*180);
    }else{
      for(j=0;j<g->dim;j++) fprintf(fp," %.4e ",g->xv[i*dim+j]);
    }
  }
  fprintf(fp,"}\n");
  /*EDGES*/
  fprintf(fp,"\nnum_edges{%d}",g->ne);
  fprintf(fp,"\ndata_edges{%d %d %d",g->seg[0],g->seg[1],g->seg[2]);
  for(i=1;i<g->ne;i++)
    fprintf(fp,"\n%d %d   %d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
  fprintf(fp,"}\n");
  /*MACROELEMENTS*/
  INT nvcube=(1<<g->dim),nvcube1=nvcube+1;
  fprintf(fp,"\nnum_macroelements{%d}\n",g->nel);
  fprintf(fp,"\ndata_macroelements{");
  for(j=0;j<nvcube1;j++)
    fprintf(fp,"%d ",g->mnodes[j]);
  for(i=1;i<g->nel;i++){
    fprintf(fp,"\n");
    for(j=0;j<nvcube1;j++)
      fprintf(fp,"%d ",g->mnodes[nvcube1*i+j]);
  }
    fprintf(fp,"}\n");
  /*MACROFACES */
  INT nvface=(1<<(g->dim-1)),nvface1=nvface+1;
  fprintf(fp,"\nnum_macrofaces{%d}\n",g->nf);
  fprintf(fp,"\ndata_macrofaces{");
  for(j=0;j<nvface1;j++)
    fprintf(fp,"%d ",g->mfaces[j]);
  for(i=1;i<g->nf;i++){
    fprintf(fp,"\n");
    for(j=0;j<nvface1;j++)
      fprintf(fp,"%d ",g->mfaces[i*nvface1+j]);
  }
  fprintf(fp,"}");fflush(stdout);
  fprintf(fp,"\n%%%%%s\n%%%s%s",					\
	  "END(input grid file)",					\
	  "---------------------------------------------------",	\
	  "---------------------------------------------------\n");
  return;
}
/**********************************************************************/
/*!
 * \fn static char **splits(char *s, const char *d, INT *num)
 *
 * \brief splits a string s into pieces using as delimiter d.
 *
 * \param num (number of pieces)
 *
 * \return the array of strings after splitting s.
 *
 * \note
 *
 */
static char **splits(char *s, const char *d, INT *num)
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
/**********************************************************************/
/*!
 * \fn void *read_mixed_data(INT nrec, INT ni, INT nr, char *the_string)
 *
 * \brief Read nrec records from a string the_string. Each record has
 *        strings which describe n_i INTs and n_r REALs. Store the
 *        output in idata[] and rdata[] all together in a void array
 *        which is the return value of the function.
 *
 * \param
 *
 * \return void array with the data read (mixed, int and double)
 *
 * \note
 *
 */
void *read_mixed_data(INT nrec, INT ni, INT nr, char *the_string)
{
  if(nrec<0||ni<0||nr<0 ||(the_string==NULL)) return NULL;
  char **w;
  INT iread,count,cni,cnr,k,j,num;
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
  k=0;
  for(count=0;count<nrec;count++){
    if(idata!=NULL) {
      cni=count*ni;
      for(j=0;j<ni;j++){
	iread=sscanf(w[k],"%d",(idata+cni+j));
	if(iread<0) break;
	k++;
      }
      if(iread<0) break;
    }
    if(rdata!=NULL) {
      cnr=count*nr;
      for(j=0;j<nr;j++){
	iread=sscanf(w[k],"%lg",(rdata+cnr+j));
	if(iread<0) break;
	k++;
      }
      if(iread<0) break;
    }
  }
  //  fprintf(stdout,"\n%%%%%%k=%d,nrec=%d,num=%d",k,nrec,num);fflush(stdout);
  for(j=0;j<num;++j) free(w[j]);    
  free(w);
  return out;
}
/*****************************************************************/
/*!
 * \fn static INT  read_data(char **clndata,input_grid *g)
 *
 * \brief reads data from an array of strings which was created from a
 *        input file for the grid generator.
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
static INT  read_data(char **clndata,input_grid *g)
{
  /* title=clndata[0]; */
  /* file_grid=clndata[1]; */
  /* file_vtu=clndata[2]; */
  /* data_coordsystems=clndata[3]; */
  /* data_vertices=clndata[4]; */
  /* data_edges=clndata[5]; */
  /* data_macroelements=clndata[6]; */
  /* data_macrofaces=clndata[7]; */
  /* dimension=clndata[8];		 */
  /* num_coordsystems=clndata[9]; */
  /* num_vertices=clndata[10]; */
  /* num_edges=clndata[11]; */
  /* num_macroelements=clndata[12]; */
  /* num_macrofaces=clndata[13]; */
  /* num_refinements=clndata[14]; */
  /* refinement_type=clndata[15]; */
  /* amr_marking_type=clndata[16]; */
  /* err_stop_refinement=clndata[17]; */
  /* print_level=clndata[18]; */
  /* num_refinements = clndata[19]; */
  /* data_refine_points = clndata[20]; */
  //
  INT i,count,j,status=0;
  void *mdata;
  INT *idata=NULL;
  /********************* coord_systems*****************/
  //  fprintf(stdout,"%s:RMD(3)",__FUNCTION__);fflush(stdout);
  mdata=read_mixed_data(g->ncsys,2,g->dim,clndata[3]);
  if(mdata!=NULL){
    idata=(INT *)mdata;
    r2c(g->ncsys,2,sizeof(INT),idata);// by rows.
    memcpy(g->ox,(mdata+g->ncsys*2*sizeof(INT)),g->ncsys*g->dim*sizeof(REAL));
    memcpy(g->syslabels,idata,g->ncsys*sizeof(INT));
    memcpy(g->systypes,(idata+g->ncsys),g->ncsys*sizeof(INT));
    free(mdata); mdata=NULL;
  } else {
    status--;
    fprintf(stderr,"\n%%%%%% WARNING (in %s):     ",__FUNCTION__);
    fprintf(stderr,"%s\n%s\n%s",					\
	    "*** Coord systems are read incorrectly or they are not specified in the input file;", \
	    "*** Switching to default (one coord system):",		\
	    "***           origin(0,0,...,0);systype=0(cartesian); syslabel=0");
    g->ncsys=1;
    g->ox=realloc(g->ox,g->dim*sizeof(REAL));
    g->systypes=realloc(g->systypes,g->ncsys*sizeof(INT));
    g->syslabels=realloc(g->syslabels,g->ncsys*sizeof(INT));
    for(i=0;i<g->dim;i++)g->ox[i]=0;
    g->syslabels[0]=0;g->systypes[0]=0;
  }

  /***** vertices: same as coordinate systems *****/
  //  fprintf(stdout,"%s:RMD(4)",__FUNCTION__);fflush(stdout);
  mdata=read_mixed_data(g->nv,2,g->dim,clndata[4]);
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
    r2c(g->nv,2,sizeof(INT),idata);// vertex labels and coord systems by rows.
    memcpy(g->xv,(mdata+g->nv*2*sizeof(INT)),g->nv*g->dim*sizeof(REAL));
    //    print_full_mat(g->nv,g->dim,g->xv,"x1");
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
    /* convert the degree coordinates to radian coordinates (polar); also convert all angles to be in [0,2*PI] */
    for(count=0;count<g->nv;count++){
      if(g->systypes[g->csysv[count]]==1){
	for(j=1;j<g->dim;j++){
	  //	fprintf(stdout,"\n(%d%d) x=%f",count,j,g->xv[count*g->dim + j]);
	  g->xv[count*g->dim + j]=zero_twopi_deg(g->xv[count*g->dim + j]);
	}
      }
    }
  } else {
    // nothing to save here bc vertices are needed.
    return 12;
  }
  /* /\***** edges *****\/ */
  //  fprintf(stdout,"%s:RMD(5)",__FUNCTION__);fflush(stdout);
  mdata=read_mixed_data(g->ne,3,0,clndata[5]);
  if(mdata!=NULL){
    memcpy(g->seg,mdata,3*g->ne*sizeof(INT));
    free(mdata); mdata=NULL;
  } else {
    status--;
    fprintf(stderr,"\n%%%%%% WARNING (in %s):     ",__FUNCTION__);
    fprintf(stderr,"\n%s\n%s\n%s",					\
	    "*** Edges are read incorrectly or they are not specified in the input file;", \
	    "*** Switching to default (one edge):",		\
	    "***          edge=(0,1);number of divisions = 1");
    g->ne=1;
    g->seg=realloc(g->seg,3*g->ne*sizeof(INT));
    g->seg[0]=0;g->seg[1]=1;g->seg[2]=1;
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
  //  fprintf(stdout,"%s:RMD(6)",__FUNCTION__);fflush(stdout);
  mdata=read_mixed_data(g->nel,(nvcube+1),0,clndata[6]);
  if(mdata!=NULL){
    memcpy(g->mnodes,mdata,g->nel*(nvcube+1)*sizeof(INT));
    free(mdata); mdata=NULL;
  } else {
    return 14;
  }

  // faces
  INT nvface=(1<<(g->dim-1));
  //  fprintf(stdout,"%s:RMD(7)",__FUNCTION__);fflush(stdout);
  mdata=read_mixed_data(g->nf,(nvface+1),0,clndata[7]);
  if(mdata!=NULL){
    memcpy(g->mfaces,mdata,g->nf*(nvface+1)*sizeof(INT));
    free(mdata); mdata=NULL;
  } else {
    status--;
    fprintf(stderr,"\n%%%%%% WARNING (in %s):     ",__FUNCTION__);
    fprintf(stderr,"%s\n%s\n%s",					\
	    "*** Faces and their codes are read incorrectly or they are not specified in the input file;", \
	    "*** Switching to default (one face):",		\
	    "***          face=(0,1,...,dim); bndry code=1 (Dirichlet)");
    g->nf=1;
    g->mfaces=realloc(g->mfaces,(nvface+1)*sizeof(INT));
    for(i=0;i<(nvface);i++) g->mfaces[i]=i;
    g->mfaces[nvface]=1;
  }

  // refine points
  //  fprintf(stdout,"%s:RMD(20)",__FUNCTION__);fflush(stdout);
  mdata=read_mixed_data(g->num_refine_points,2,g->dim,clndata[20]);
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
    if(g->num_refine_points){
      r2c(g->num_refine_points,2,sizeof(INT),idata);// vertex labels and coord systems by rows.
      memcpy(g->data_refine_points,(mdata+g->num_refine_points*2*sizeof(INT)),g->num_refine_points*g->dim*sizeof(REAL));
    }
    //    print_full_mat(g->nv,g->dim,g->xv,"x1");
    /* for(count=0;count<g->nv;count++){ */
    /*   fprintf(stdout,"\nrec=%d (",count); */
    /*   for(j=0;j<g->dim;j++){ */
    /* 	fprintf(stdout,"%f ",g->xv[count*g->dim + j]); */
    /*   } */
    /*   fprintf(stdout,")"); */
    /* } */
    /* fprintf(stdout,"\n"); */
    free(mdata); mdata=NULL;
    /* fprintf(stdout,"\n"); */
    /* for(j=0;j<g->nv;j++){ */
    /*   fprintf(stdout,"\nl=%d;t=%d;",g->labelsv[j],g->csysv[j]); */
    /* } */
    /* fprintf(stdout,"\n");fflush(stdout); */
    /* convert the degree coordinates to radian coordinates (polar); also convert all angles to be in [0,2*PI] */
    for(count=0;count<g->nv;count++){
      //      fprintf(stdout,"\ncount=%d,g->dim=%d,(systype)=%d(%d)",count,g->dim,g->systypes[g->csysv[count]],g->csysv[count]);fflush(stdout);
      if(g->systypes[g->csysv[count]]==1){
	for(j=1;j<g->dim;j++){
	  if(!(g->num_refine_points)) continue;
	  fprintf(stdout,"\ncount=%d;j=%d;indx=%d, g_xv=%f",count,j,count*g->dim+j,g->xv[count*g->dim + j]); fflush(stdout);
	  g->data_refine_points[count*g->dim + j]=zero_twopi_deg(g->data_refine_points[count*g->dim + j]);
	}
      }
    }
  } else {
    // nothing to save here bc vertices are needed.
    return 20;
  }


  if(status<0)
    fprintf(stdout,"\n%d warnings were issued;",status);
  //  input_grid_print(g);
  return status;
}
/**********************************************************************/
/*!
 * \fn void x_out(const char *pattern, size_t le)
 *
 * \brief prints a string cutting it at the closest blank space
 *        with location <=le.
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void x_out(const char *pattern, size_t le)
{
  int i;
  fprintf(stderr, "\n\n\n *** ERROR(%s)::::   \n     UNBALANCED \"{}\" near or before \"",__FUNCTION__);
  for (i=0;i<(le-1);++i)
    fprintf(stderr, "%c",*(pattern+i));
  fprintf(stderr, "\"\n");
  exit(12);
}
/**********************************************************************/
/*!
 * \fn char *make_string_from_file(FILE *the_file,size_t *length_string)
 *
 * \brief puts the content of a file in a string, ignoring all comments which start with %.
 *
 * \param the_file is the input channel we read from.
 *
 * \return a string formed by at most of
 *         MAX_CHARS_INPUT_GRID_FILE(macro.h) are read with all empty
 *         lines consequtive spaces and comments ignored.
 *
 * \note
 *
 */
static char *make_string_from_file(FILE *the_file, size_t *length_string)
{
  char *file2str;
  char ch, ch_next;
  INT count=0,i,j,flag,maxcount=((INT )MAX_CHARS_INPUT_GRID_FILE); //maxcount=(1<<15-1);
  while(feof(the_file)==0)
    {
      ch = fgetc(the_file);
      if(ch) count++;
      if(count>maxcount) {
	fprintf(stderr,"\n%%%%%% WARNING (in %s): The grid input file is too large. Truncating the input to %d chars\n",__FUNCTION__,count);
	break;
      }
    }
  count--;
  //  fprintf(stdout,"\n%%%%%%count=%d",count);
  //  if(count<0) count=0;
  //  else if(count>maxcount) {count=maxcount;}
  //  fprintf(stdout,"\n%%%%%%*** count=%d",count);
  i = count*sizeof(char);
  file2str = (char *)malloc(i);
  memset(file2str,0,i);
  rewind(the_file);
  //fprintf(stderr,"\nNumber of characters in the file %i\n", count);
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
      //      if(ch == ' ' || ch == '{' || ch ==  '}'){
      if(ch == ' ' || ch == '{' || ch ==  '}' || ch == '='){
	do{
	  ch_next = fgetc(the_file);
	  if(ch_next == '\n' || ch_next == '\t') ch_next = ' ';
	  ++j;
	}  while(ch_next == ' ');
	//	if((ch == ' ')&& (ch_next ==  '{' || ch_next ==  '}')){
	if((ch == ' ' || ch=='=')&& (ch_next ==  '{' || ch_next ==  '}')){
	  ch = ch_next;
	  flag=0;
	} else {
	  /*		  printf(" %i ",ch);*/
	  *(file2str+i) = ch;
	  ++i;
	  ch = ch_next;
	  flag=0;
	}
      } else if( ch != EOF ) {
	/*	      printf(" %i ",ch);*/
	*(file2str+i) = ch;
	++i;
	flag=1;
      }
    }
    /*           printf("\n i,j, %i %i  %i\n",i,j,j-count);*/
  }
  if(i && (file2str[i-1]==' '))   i--;
  /* end of string */
  INT any=0;
  /**/
  for(j=i;j>=0;j--) if(file2str[j]=='}'){any=1;break;}
  /**/
  if(any) i=j+1;  else i=0;
  file2str = (char *)realloc(file2str,(i+1)*sizeof(char));
  file2str[i] = '\0';
  *length_string = strlen(file2str)+1;
  //  fprintf(stdout,"\ncount=%d:%s",i,file2str);fflush(stdout);
  /* for(j=0;j<i;j++){ */
  /*   if(!file2str[j]) continue; */
  /*   file2str[j]=toupper(tolower(file2str[j])); */
  /* } */
  /*
      fprintf(stderr,"\nNumber of characters in the supressed string %li\n",strlen(file2str));
      fprintf(stdout,"\nString=%s\n",file2str);
  */
  return file2str;
}
/**********************************************************************/
/*!
 * \fn char *get_substring(const char *pattern, size_t
 * *length_substring, char *the_string)
 *
 * \brief Finds the  pattern in the_string.
 *
 * \param
 *
 * \return returns a pointer to the part of the string after
 * pattern. if the pattern is not found returns '\0'
 *
 * \note
 *
 */
static char *get_substring(const char *pattern,	\
		    size_t *length_substring,	\
		    char *the_string)
{
  /*
     get substring from a string matching a pattern; it can probably
     be done with strtok() more efficiently...
  */
  size_t le;
  INT i;
  char *found=NULL, *wrk=NULL;
  //  char ch;
  le = strlen(pattern);
  found = (char *) strstr(the_string,pattern);
  if(found != NULL){
    found += le;
    wrk = (char *) strstr(found,"}");
    if(wrk == NULL ){ x_out(pattern,le);}
    *length_substring=strlen(found)-strlen(wrk);
    *wrk = '\t';
    wrk = (char *) strstr(found,"\t");
    i = strlen(found)-strlen(wrk);
    if(i != *length_substring ){ x_out(pattern,le); }
  }else{
    found=strndup("",1);
    *length_substring=strlen(found);
  }
  return found;
}
/**********************************************************************/
/*!
 * \fn char *safe_parse(const char *sinp, const char *warn0,
 *		 const char *default_s, const INT max_length)
 *
 * \brief If sinp is not empty returns a copy of it up to max_length
 *        chars; if sinp is empty then returns default_s and issues a
 *        standard warning
 *
 * \param default_s default string
 *
 * \param warn0 string to be included in the warning message.
 *
 * \param sinp is the input string
 *
 * \return either the  default string (if sinp is null) or the input string.
 *
 * \note
 *
 */
static char *safe_parse(const char *sinp,		\
		 const char *warn0,		\
		 const char *default_s,		\
		 const INT max_length)
{
  // if the input string sinp is empty, return the default string and issue
  // a warning; else return the input string;
  char *s;
  if(strlen(sinp)<=0){
    s=strndup(default_s,max_length);
    fprintf(stderr,"\n%%%%%%WARNING (input file): %s is not set correctly. Setting it to the default value: %s",warn0,s);
  } else {
    s=strndup(sinp,max_length);
    //  fprintf(stdout,"\n%%%%%%%s={%s}",warn0,s);
  }
  return s;
}
/**********************************************************************/
/*!
 * \fn static INT check_input(char * file2str, input_grid *g,
		char **indata, char **notes, INT numel_data)
 *
 * \brief reads all mixed data (strings, int and double) from the
 *        string file2str. If the data is read OK, then the values are
 *        stored in the input_grid structure g. If success, returns 0;
 *        if there is a problem reading the data a warning is issued
 *        and the function returns nonzero value
 *
 * \param file2str a string with input data;
 *
 * \param indata are the patterns we look for to extract the
 *         corresponding values
 *
 * \param notes are strings used in the warning messages.
 *
 * \return 0 on successful reading of all data; non- there is a problem.
 *
 * \note
 *
 */
static INT check_input(char * file2str, input_grid *g,	\
		char **indata, char **notes,	\
		INT numel_data)
{
  INT k,status=0;
  size_t *lengths=calloc(numel_data,sizeof(size_t));
  char **clndata=malloc(numel_data*sizeof(char *));
  INT *iread=calloc(numel_data,sizeof(INT));
  /*clndata are only addresses in file2str containing the information to be read later*/
  for(k=0;k<numel_data;k++){
    clndata[k] = get_substring(indata[k],(lengths+k), file2str);
  }
  for(k=0;k<numel_data;k++)
    clndata[k][lengths[k]] = '\0';
  /* initialize */
  /* ... PARSE ... strings */
  g->title=safe_parse(clndata[0],"Title","Untitled",256);
  g->fgrid=safe_parse(clndata[1],"Filename for the grid file","mesh.haz",256);
  g->fvtu=safe_parse(clndata[2],"Filename for the VTU file","mesh.vtu",256);
  /*INTEGERS*/
  iread[8]=sscanf( clndata[8],"%d",&g->dim); // the dimension of the problem.
  if(iread[8]<0) return 8;
  //
  iread[9]=sscanf( clndata[9],"%d",&g->ncsys);//
  if(iread[9]<0) g->ncsys=-1; // this is fixable fixable;
  //
  iread[10]=sscanf( clndata[10],"%d",&g->nv);//
  if(iread[10]<0) return 10;
  //
  iread[11]=sscanf( clndata[11],"%d",&g->ne);//
  if(iread[11]<0) g->ne=-1; // this is fixable;
  //
  iread[12]=sscanf( clndata[12],"%d",&g->nel);//
  if(iread[12]<0) return 12;// no
  //
  iread[13]=sscanf(clndata[13],"%d",&g->nf);//
  if(iread[13]<0) g->nf=-1;//ok with some check
  //
  iread[14]=sscanf(clndata[14],"%d",&g->nref);//
  if(iread[14]<0)g->nref=0;//ok
  //
  iread[15]=sscanf(clndata[15],"%d",&g->ref_type);//
  if(iread[15]<0)g->ref_type=-1;//ok
  //
  iread[16]=sscanf(clndata[16],"%d",&g->mark_type);//
  if(iread[16]<0)g->mark_type=0;//ok
  //
  iread[17]=sscanf(clndata[17],"%lg",&g->err_stop);//
  if(iread[17]<0)g->err_stop=-1e-10;//ok
  //
  iread[18]=sscanf(clndata[18],"%hd",&g->print_level);//
  if(iread[18]<0)g->print_level=0;//ok
  //
  iread[19]=sscanf(clndata[19],"%d",&g->num_refine_points);//
  if(iread[19]<0)g->num_refine_points=0;//ok
  /*FREE*/
  free(iread);
  /*FIX*/
  if(g->ncsys<=0){g->ncsys=-1;}
  if(g->ne<=0){g->ne=-1;}
  if(g->nf<=0){g->nf=-1;}
  if(g->num_refine_points<=0){g->num_refine_points=0;}
    //
  /* now read the data*/
  input_grid_arrays(g); // allocation only
  /**/
  status=read_data(clndata,g);
  /*FREE*/
  for(k=0;k<numel_data;k++){
    if(lengths[k]) continue;
    //    fprintf(stdout,"\nlengths[%d]=%ld,clndata=%s\n",k,lengths[k],clndata[k]);fflush(stdout);
    if(clndata[k]) free(clndata[k]);
  }
  free(lengths);
  free(clndata);
  return status;
}
/**********************************************************************/
/*!
 * \fn input_grid *parse_input_grid(FILE *the_file)
 *
 * \brief parses the input for grid generation from a file on the
 *        stream "the_file".
 *
 * \param  the_file is a stream open with fopen.
 *
 * \return returns a pointer to a structure describing the
 *         macroelement grid.
 *
 * \note
 *
 */
input_grid *parse_input_grid(FILE *the_file)
{
  input_grid *g=malloc(sizeof(input_grid));
  INT numel_data,k,knext;
  char *file2str=NULL;
  size_t ll=-1,length_string=0;
  char **indata,**notes;
  indata=input_strings(&numel_data);
  notes=malloc(numel_data*sizeof(char *));
  for(k=0;k<numel_data;k++){
    /* fprintf(stdout,"\n%d %ld: %s",k,ll,indata[k]); fflush(stdout); */
    ll=strlen(indata[k]);
    notes[k]=strndup(indata[k],ll);
    /*remove '{' at the end of indata*/
    if(ll>0) notes[k][ll-1]='\0'; else notes[k][0]='\0';
  }
  /* get all substrings */
  //  size_t *lengths;
  file2str=make_string_from_file(the_file,&length_string);
  //  fprintf(stdout,"\neve0[%ld]=%s\n",length_string,file2str);
  length_string=strlen(file2str);
  file2str=realloc(file2str,(length_string+1)*sizeof(char));
  file2str[length_string]='\0';
  //  fprintf(stdout,"\neve1[%ld]=%s\n",length_string,file2str);
  k=check_input(file2str,g,indata,notes,numel_data);
  switch(k){
  case 8: case 10: case 12:
    // issue a warning message and switch to the unit cube in 3D.
    file2str = strndup( DEFAULT_GRID_DATA_ , 1024);
    length_string=strlen(file2str);
    file2str=realloc(file2str,(length_string+1)*sizeof(char));
    file2str[length_string]='\0';
    knext=check_input(file2str,g,indata,notes,numel_data);
    break;
  case 0:
    if(g->print_level>0)
      fprintf(stdout,"\nEXAMPLE: %s\n",g->title);
    //    input_grid_example_file(g);
    for(k=0;k<numel_data;k++) free(indata[k]);
    free(indata);
    for(k=0;k<numel_data;k++) free(notes[k]);
    free(notes);
    free(file2str);
    return g;
  default:
    break;
  }
  for(k=0;k<numel_data;k++) free(indata[k]);
  free(indata);
  for(k=0;k<numel_data;k++) free(notes[k]);
  free(notes);
  free(file2str);
  if(knext){
    //issue an error message and exit
    fprintf(stderr,"\n%s\n%s\n","%%%% **** ERROR: the input is not as expected.",
                   "   Consult examples of input files and re-run. ***");
    exit(12);
  }  else {
    fprintf(stderr,"\n\n%s\n%s\n%s\n","%%%%%% ****ERROR(input file): One of the following", \
	    "%%%%%%     (dimension), (data_macroelements), (data_vertices)", \
	    "%%%%%%     could is incorrect read. EXITING.");
    //    if(g->print_level>3)
    //      input_grid_print(g);
    input_grid_example_file(stderr,g);
    exit(16);
    //    return g;
  }
}
/**********************************************************************/
/*!
 * \fn void set_input_grid(INT *nd,input_grid *g,cube2simp *c2s)
 *
 * \brief Every edge is put into a subset, i.e. two edges
 *        (v(i1),v(i2)) and (v(j1),v(j2)) are considered equivalent
 *        iff (i2-i1)=(j2-j1).  The number of divisions in an
 *        equivalent set of edges is taken to be the largest from the
 *        equivalence class.  OUTPUT array is a "dim" array and for
 *        each direction gives the number of partitions.
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void set_input_grid(INT *nd,input_grid *g,cube2simp *c2s)
{
  INT i,j,k,iri,ici,pmem;
  pmem=2*g->nv;
  if(pmem<2*g->ne) pmem=2*g->ne;
  if(pmem<2*g->nel) pmem=2*g->nel;
  if(pmem<2*g->nf) pmem=2*g->nf;
  //  fprintf(stdout,"pmem=%d",pmem);fflush(stdout);  
  INT *p=calloc(pmem,sizeof(INT));// permutation and inverse permutation;
  //
  memset(p,0,pmem);
  for (i=0;i<g->ne;i++){
    iri=g->seg[3*i];
    ici=g->seg[3*i+1];
    if(iri<ici){
      g->seg[3*i]=iri;
      g->seg[3*i+1]=ici;
    } else {
      g->seg[3*i]=ici;
      g->seg[3*i+1]=iri;
    }
    /* set up divisions */
    j=g->seg[3*i+1]-g->seg[3*i]; // should be always positive;
    //    fprintf(stdout,"\n%%z123=%d:(%d,%d);%d",i,3*i,3*i+1,g0->seg[3*efound[i]+2]);
    if(g->seg[3*i+2]>p[j])
      p[j]=g->seg[3*i+2];
  }
  for (i=0;i<g->ne;i++){
    j=g->seg[3*i+1]-g->seg[3*i];
    g->seg[3*i+2]=p[j];
    //    fprintf(stdout,"\n[%d,%d]:div=%d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
  }
  for (i=0;i<g->ne;i++){
    j=g->seg[3*i+1]-g->seg[3*i];
    g->seg[3*i+2]=p[j];
    //    fprintf(stdout,"\n[%d,%d]:div=%d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
  }
  /*ORDER*/
  ilexsort(g->ne, 3,g->seg,p);
  k=0;
  for (i=0;i<g->ne;i++){
    if(g->seg[3*i]) continue;
    //    j=g->seg[3*i+1]-g->seg[3*i]-1;
    p[k]=g->seg[3*i+2];
    k++;
    //    fprintf(stdout,"\n[%d,%d]:div=%d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
  }
  ////////////////////////////////////////
  //  p=realloc(p,g->dim*sizeof(INT)); // realloc to dimension g->dim
  ///////////////////////////////////
  //  for (i=0;i<g->dim;i++){
  //    fprintf(stdout,"\ndirection:%d; div=%d",i,p[i]);
  //  }
  //  input_grid_print(g);
  //  print_full_mat_int(g->ne,3,g->seg,"med");
  //  print_full_mat_int(g->nf,(c2s->nvface+1),g->mfaces,"mf");
  //  print_full_mat_int(g->nel,(c2s->nvcube+1),g->mnodes,"mel");
  ///newnew
  memcpy(nd,p,g->dim*sizeof(int));
  free(p);
  return;
}
/**********************************************************************/
/*!
 * \fn void set_edges(input_grid *g0,cube2simp *c2s)
 *
 * \brief adds missing edges to input grid g0 so that all edges are
 *        present.
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void set_edges(input_grid *g0,cube2simp *c2s)
{
  /*adds missing edges to input grid*/
  INT newne,i,j,k,kel,ke,swp,cols;
  INT j01[2],k01[2];
  INT nvcube=c2s->nvcube;
  cols=3;
  INT *newseg=calloc(cols*c2s->ne*g0->nel,sizeof(INT));
  INT *mnodes=calloc((c2s->nvcube+1),sizeof(INT));
  newne=0;
  INT found=0;
  INT nseg_chk=0;
  ////check if a set macroelement segment is an interior edge for the cube/
  for(ke=0;ke<g0->ne;ke++){
    k01[0]=g0->seg[cols*ke];
    k01[1]=g0->seg[cols*ke+1];
    for(kel=0;kel<g0->nel;kel++){
      memcpy(mnodes,(g0->mnodes+kel*(nvcube+1)),(nvcube+1)*sizeof(INT));
      found=0;
      for(i=0;i<c2s->ne;i++){
	j01[0]=mnodes[c2s->edges[2*i]];
	j01[1]=mnodes[c2s->edges[2*i+1]];
	//      fprintf(stdout,"\n*YYY*YYY j01=(%d,%d)",j01[0],j01[1]);
	if(j01[0]>j01[1]){swp=j01[0];j01[0]=j01[1];j01[1]=swp;}
	if((k01[0]==j01[0])&&(k01[1]==j01[1])){
	  found=1;break;
	}	
      }
      if(found) break;
    }
    if(found) {
      memcpy(g0->seg+cols*nseg_chk,g0->seg+cols*ke,cols*sizeof(INT));
      nseg_chk++;
    } else if(g0->print_level>1){
      fprintf(stdout,							\
	      "\n%%%%%% WARNING (in %s): Removing input segment=(%d,%d) as it is NOT on the boundary of any macroelement\n", \
	      __FUNCTION__,k01[0],k01[1]);
    }
  }
  g0->ne=nseg_chk;
  g0->seg=realloc(g0->seg,g0->ne*cols*sizeof(INT));
  /////////////////////////////////////////////////////
  for(kel=0;kel<g0->nel;kel++){
    memcpy(mnodes,(g0->mnodes+kel*(nvcube+1)),(nvcube+1)*sizeof(INT));
    for(i=0;i<c2s->ne;i++){
      j01[0]=mnodes[c2s->edges[2*i]];
      j01[1]=mnodes[c2s->edges[2*i+1]];
      if(j01[0]>j01[1]){swp=j01[0];j01[0]=j01[1];j01[1]=swp;}
      found=0;
      for(ke=0;ke<g0->ne;ke++){
	k01[0]=g0->seg[cols*ke];
	k01[1]=g0->seg[cols*ke+1];
	//	fprintf(stdout,"\n%s(%d):*YYY*YYY k01=(%d,%d)?=?(%d,%d)",__FUNCTION__,newne,k01[0],k01[1],j01[0],j01[1]);fflush(stdout);
	if((k01[0]==j01[0])&&(k01[1]==j01[1])){
	  found=1;break;
	}
      }
      if(!found){
	newseg[cols*newne]=j01[0];
	newseg[cols*newne+1]=j01[1];
	newseg[cols*newne+2]=1;
	newne++;
      }
    }
  }
  //  fprintf(stdout,"\n%%%%%%* newne=%d",newne); fflush(stdout);
  if(!newne){
    free(newseg);
    free(mnodes);
    return;
  }
  //  fprintf(stdout,"\n%s:newne=%d",__FUNCTION__,newne);fflush(stdout);
  newseg=realloc(newseg,cols*newne*sizeof(INT));
  INT *p=calloc(newne,sizeof(INT));
  ilexsort(newne, cols,newseg,p);
  free(p);
  // remove dupps
  INT m,li1,ic[cols];
  k=newne-1;
  i=0;j=0;
  while (i<k){
    if(j==0) {for(m=0;m<cols;m++) {newseg[m]=newseg[cols*i+m];}}
    for(m=0;m<cols;m++)  {ic[m]=newseg[cols*j+m];}
    while(i<k) {
      li1=0;
      for(m=0;m<cols;m++){li1+=abs(ic[m]-newseg[cols*i+cols+m]);}
      if(li1>0){
  	j++;i++;
  	for(m=0;m<cols;m++){newseg[cols*j+m]=newseg[cols*i+m];}
  	break;
      }
      i++;
      //      fprintf(stdout,"i=%i\n",i);
    }
    //    fprintf(stdout,"i=%i, j=%i\n",i,j);
  }
  i++;j++; newne=j;
  g0->seg=realloc(g0->seg,(cols*(g0->ne+newne))*sizeof(INT));
  memcpy((g0->seg+cols*g0->ne),newseg,cols*newne*sizeof(INT));
  g0->ne+=newne;
  free(newseg);
  free(mnodes);
  return;
}
/**********************************************************************/
/*!
 * \fn INT set_ndiv_edges(input_grid *g, input_grid *g0,
 *         cube2simp *c2s, INT **nd, const INT iter)
 *
 * \brief For a given global input grid g0 creates local input grids
 *        for every macroelement and computes the divisions in each
 *        direction for it. It is used iteratively in macro_split to
 *        set the correct divisions for every macroelement.  The
 *        input_grid *g0 should be all set, and the input_grid *g
 *        should have all its scalar values set.  nd is the array with
 *        the divisions, it must be g0->nel by c2s->n.
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
INT set_ndiv_edges(input_grid *g, input_grid *g0,		\
		   cube2simp *c2s, INT **nd, const INT iter)
{
  INT kel0,i,j0,j1,swp,kel,ke,k0,k1,ndiv;
  INT nel0=g0->nel,nvcube=c2s->nvcube;
  /*for easier reference*/
  //  INT *efound=calloc(c2s->ne*(g0->nel+1),sizeof(INT));
  INT *efound=malloc(c2s->ne*(g0->nel+1)*sizeof(INT));
  INT *e0found=efound + c2s->ne;
  for(i=0;i<g0->ne;i++){
    /* fprintf(stdout,"\n%s:i=%d,totals=%d,c2s_ne=%d,g0_ne=%d,g0_nel=%d",	\ */
    /* 	    __FUNCTION__,i,c2s->ne*g0->nel,c2s->ne,g0->ne,g0->nel);fflush(stdout); */
    e0found[i]=-1;
  }
  // make all divisions > 0
  for(ke=0;ke<g0->ne;ke++){
    ndiv=abs(g0->seg[3*ke+2]);
    if(ndiv<=0)ndiv=1;
    g0->seg[3*ke+2]=ndiv;
  }
  for(ke=0;ke<g0->ne;ke++)
    e0found[ke]=g0->seg[3*ke+2];
  //  print_full_mat_int(g0->ne,3,g0->seg,"seg0");
  for(kel0=0;kel0<nel0;kel0++){
    if((iter%2)) kel=nel0-kel0-1; else kel=kel0;
    // macroelement by macroelement try to find the edge divisions
    for(i=0;i<c2s->ne;i++){
      g->seg[3*i]=c2s->edges[2*i];
      g->seg[3*i+1]=c2s->edges[2*i+1];
      g->seg[3*i+2]=-1;
      efound[i]=-1;
    }
    memcpy(g->mnodes,(g0->mnodes+kel*(nvcube+1)),(nvcube+1)*sizeof(INT));
    for(i=0;i<c2s->ne;i++){
      j0=g->mnodes[c2s->edges[2*i]];
      j1=g->mnodes[c2s->edges[2*i+1]];
      if(j0>j1){swp=j0;j0=j1;j1=swp;}
      for(ke=0;ke<g0->ne;ke++){
	k0=g0->seg[3*ke];
	k1=g0->seg[3*ke+1];
	if((k0==j0)&&(k1==j1)){
	  g->seg[3*i+2]=g0->seg[3*ke+2];
	  efound[i]=ke;
	}
      }
      //      if(iter==0){
      //	fprintf(stdout,"\n%%iter=%d;Element:%d, edge=(%d,%d);",iter,kel,j0,j1);fflush(stdout);
      //}
    }
    //    if(iter==0){
    //      input_grid_print(g);
    //      fflush(stdout);
    //}
    ////////////////////////////////////////////////////////////
    set_input_grid(nd[kel],g,c2s);
    for(i=0;i<g->ne;i++){
      if(efound[i]<0)
	continue;
      ke=efound[i];
      g0->seg[3*ke+2]=g->seg[3*i+2];
    }
  }
  INT chng=0;
  for(ke=0;ke<g0->ne;ke++){
    k1=abs(e0found[ke]-g0->seg[3*ke+2]);
    if(k1>chng)chng=k1;
  }
  //  print_full_mat_int(g0->ne,3,g0->seg,"seg1");
  //  print_full_mat_int(1,c2s->n,nd[kel],"ndnd");
  //  fprintf(stderr,"\nchng=%d",chng);
  free(efound);
  return chng;
}
/*EOF*/
