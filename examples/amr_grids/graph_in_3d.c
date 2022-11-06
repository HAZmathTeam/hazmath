/*! \file examples/amr_grids/amr_grids.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 2019/01/09.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This example exports a 1d simplicial complex embedded in 3d
 * to a vtk file
 *
 */
/*********** HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
/**********************************************************************/
static void iread_eof(const char *fname,ivector *ivec)
/* reads integers until eof. skips any control chars at the beginning
   of the file. result is stored in an ivector ivec; ivec->val is
   allocated here*/
{
  INT i,k,count=0,ichk=-1;
  long long lint0;
  char ch;
  FILE *fp=fopen(fname,"r");
  ch=fgetc(fp); while((INT )ch < 0){ch=fgetc(fp);count++;}
  if(count){
    fprintf(stdout,"%%Read: %lld control chars...", (long long int)count);
    fseek(fp,count*sizeof(char),SEEK_SET);
  } else rewind(fp);
  k=0;
  while(1){
    if(feof(fp)) break;
    ichk=fscanf(fp,"%lld", &lint0);
    if(ichk<0) break;
    k++;
  }
  /* read now k numbers from the file */
  if(count){
    fseek(fp,count*sizeof(char),SEEK_SET);
  } else {
    rewind(fp);
  }
  /* we allocate k array*/
  ivec->row=k;
  ivec->val=(INT *)calloc(ivec->row,sizeof(INT));
  rveci_(fp,ivec->val,&ivec->row);
  fprintf(stdout,"Read %lld integers\n",(long long int)ivec->row);
  fclose(fp);
  return;
}
/**********************************************************************************************/
INT main(INT   argc,   char *argv[])
{
  INT i,j,k;
  FILE *fp;
  char *data_coords,*data_seg;
  features feat;
  feat.n=0;feat.nbig=0;feat.x=NULL;feat.fill=-1e20;feat.fpf=NULL;
  // nodes coords of the 1d graph embedded in 3d. 
  data_coords=strdup("./try_features.txt");
  data_seg=strdup("./try_seg.txt");
  //
  //
  feat.fpf=fopen(data_coords,"r");
  feat.nbig=3;// this is the number of coordinates for the feature points:
  feat.n=3; //this is the same as bnbig here as every point has 3 coordinates to be read. 
  feat.fill=-1e20;
  //
  scomplex *sc=NULL;
  j=features_r(&feat,sc,(INT )0); free(data_coords);
  //
  ivector ivec;
  ivec.row=0;ivec.val=NULL;
  iread_eof(data_seg,&ivec);  free(data_seg);
  INT dim,dimbig,ns,nv;
  //
  dim=1; dimbig=3;
  ns=(INT )(ivec.row/(dim+1));  nv=feat.nf;
  //
  sc=haz_scomplex_init(dim,ns,nv,dimbig);
  free(sc->nodes);
  free(sc->x);
  sc->nodes=ivec.val;// this is freed as a member of sc;
  sc->x=feat.x;// this is also freed as a member of sc.
  /* WRITE THE OUTPUT vtu file for paraview:    */
  vtkw("output/1d_graph.vtu",sc,0,1.);
  /*FREE: the input grid is freed here, because it has the filenames in it*/
  haz_scomplex_free(sc);
  return 0;
}
