#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "hazmath.h"
/*********************************************************************/
static void print_cells(const char *s,INT nfx, INT num,	\
			INT *cell,			\
			INT *flag,			\
			INT *head,			\
			INT *next,			\
			INT *back,			\
			INT *fixlst)
{
  INT hc,hcnext,f,v,b,c,strt;
  hc=0;
  //  fprintf(stdout,"\n%s:number=%d:",s,num);
  fprintf(stdout,"\n%%%%%%\tHEAD(ALL):|%3d|%3d|%3d|%3d|;",flag[0],head[0],next[0],back[0]);fflush(stdout);
  while(1){
    hcnext=head[hc];
    if(hcnext<0) break;
    strt=next[hcnext];
    //    fprintf(stdout,"\n\theader_cell=%3d; list_strt=%3d",hcnext,strt);
    //    fprintf(stdout,"\n");
    if(strt>=0) {
      fprintf(stdout,"\n%%%%%%\theader_cell=%3d:|%3d|%3d|%3d|%3d|* ",hcnext,flag[hcnext],head[hcnext],next[hcnext],back[hcnext]);fflush(stdout);
    }
    while(strt>=0){
      v=head[strt];
      c=cell[v];
      f=flag[c];
      b=back[c];
      fprintf(stdout,"%6d:|%3d|v=%3d|%3d|%3d|;",strt,f,v,next[strt],b);fflush(stdout);
      strt=next[strt];
    }
    //    jjj++;
    //    if(jjj>200) break;
    //    fprintf(stdout,"\njjj=%d\n",jjj);
    hc=hcnext;
  }
  fprintf(stdout,"\n");
  /* INT i0; */
  /* for(i0 = 0;i0< nfx;i0++){ */
  /*   fprintf(stdout,"\n%%%%%%nfx=%3d,xFIXLIST[%d]=%d,flag=%d",nfx,i0,fixlst[i0],flag[fixlst[i0]]);fflush(stdout); */
  /* } */
  return;
}
/**********************************************************************/
void lexi(INT n, INT *ia, INT *ja, INT *iord, INT *iord1, \
	  INT *level, INT *tree)
{  /*
     Output: iord is the lexicographical ordering which for
     triangulated graph is an ordering withouut fill in
   */
  /* 
   *    Lexicographical (recurrsive) breath first search.       
   *	Following the algol-like pseudocode from 
   *	Rose, Donald J. ; Tarjan, R. Endre ; Lueker, George S. 
   *	Algorithmic aspects of vertex elimination on graphs. 
   *	SIAM J. Comput. 5 (1976), no. 2, 266–283 (MR0408312).
  */
  INT c,h,v,w,p;  
  INT ibeg,iend,i,iv,nfx,ijk,i0;
  /* 
   *  two lists: queue and sets of vertices. the queue is described by
   *  head[] and backp[] and the sets are described by next[] and
   *  back[]. Each cell "c" is either a header cell or describes a
   *  set of vertices that have the same label.
   *
  */
  INT nmem=n;
  // allocate memory;
  INT *cell=calloc(n,sizeof(INT));
  for (i=0;i<n;i++) cell[i]=-1;
  INT *flag=calloc(nmem,sizeof(INT));
  INT *head=calloc(nmem,sizeof(INT));
  INT *next=calloc(nmem,sizeof(INT));
  INT *back=calloc(nmem,sizeof(INT));
  INT *fixlst=calloc(nmem,sizeof(INT));
  for (i=0;i<nmem;i++){
    fixlst[i]=-1;
    head[i]=-1;
    back[i]=-1;
    next[i]=-1;
    flag[i]=-1;
  }
  /*     assign empty label to all vertices */
  head[0] = 1; // the first cell is the head of the queue; it does not
	       // describe any set of vertices, but points to the
	       // cell=1 which will describe the first set of
	       // vertices.
  back[1] = 0; // the second cell points to the first cell as its predecessor
  /**/
  head[1] = -1;   
  back[0] = -1;
  next[0] = -1;
  flag[0] = -1;
  flag[1] = -1;
  // cell 0 is the beginning; cell 1 is header cell for the first set;
  c=2; // thus the first empty cell is cell 2. 
  //  fprintf(stdout,"\n");fflush(stdout);
  for(iv=0;iv<n;iv++){
    v = iord[iv]; // vertex label
    //    fprintf(stdout,"%d ",v);fflush(stdout);
    head[c] = v;  //head(cell)=element;
    cell[v] = c;  //cell(vertex)=cell.
    next[c-1] = c;
    flag[c] = 1;
    back[c] = c-1;
    c++;
    if(c>=nmem){
      nmem=c+1;
      flag=realloc(flag,nmem*sizeof(INT));
      head=realloc(head,nmem*sizeof(INT));
      next=realloc(next,nmem*sizeof(INT));
      back=realloc(back,nmem*sizeof(INT));
      head[c]=-1; back[c]=-1; next[c]=-1; flag[c]=-1;
    }
    iord1[v] = -1;
    level[v] = -1;
    tree[v] = -1;
  }
  next[c-1] = -1;
  for (i = n-1;i>=0;i--){
    /*********************************************************************/
    nfx=0;
    //    print_cells("%%%%first:",nfx,i,cell,flag,head,next,back,fixlst);
    //    if(i<n-10) break;
    //C  Skip empty sets
    while(next[head[0]] < 0){
      head[0] = head[head[0]];
      back[head[0]]=0;
    }
    //C  pick the cell for the next vertex to number
    p=next[head[0]];
    //C     C ADDITION BY ltz
    v = head[p];
    /* fprintf(stdout,"\n%%%%%%AAA:p=%3d,v=%3d,cell[v]=%3d",p,v,cell[v]);fflush(stdout); */
    //C  delete cell of vertex from set
    next[head[0]] = next[p];
    next[back[cell[v]]]=next[cell[v]];
    if(next[cell[v]]>=0){
      back[next[cell[v]]]=back[cell[v]];
    }
    //C assign number i to v
    iord[i] = v;
    iord1[v] = i;
    if(level[v]<0) level[v]=0;
    // fprintf(stdout,"\n*** NUMBERING at (%d) v=%d: inv(%d)=%d",i,v,v,iord1[v]);
    nfx=0;
    ibeg = ia[v]-1;
    iend = ia[v+1]-1;
    for(ijk = iend;ijk>ibeg;ijk--){
      w = ja[ijk];
      if((tree[w]<0)&&(level[w]<0)) tree[w]=v;
      if(level[w]<0) level[w]=level[v]+1;
      //fprintf(stdout,"\n%%%% %3d:v=%3d,w=%3d,h[0]=%3d,n(h(0))=%3d",ijk-ibeg,v,w,head[0],next[head[0]]);
      //      fprintf(stdout,"\n%%%%ibeg=%3d,iend=%3d,ijk=%3d,v=%3d,w=%3d,iord(v)=%3d,iord(w)=%3d",ibeg,iend,ijk,v,w,iord[v],iord[w]);fflush(stdout);
      if(iord1[w] < 0) {
	// delete cell of w from the set (this will go into a new set). 
	next[back[cell[w]]]=next[cell[w]]; // skip cell[w] in the element list
	if(next[cell[w]] >=0){
	  back[next[cell[w]]]=back[cell[w]]; // if there was a cell
					     // pointed to as next by
					     // cell(w), then put its
					     // prev. cell to be the
					     // prev(cell(w)) and not
					     // cell(w);
	}
	h = back[flag[cell[w]]];  // this is the header set cell which
				  // is a predecessor of the set header cell
				  // containing w
	//	fprintf(stdout,"\n%%%%i=%3d,v=%3d,w=%3d,iord1(%3d)=%3d,iord1(%3d)=%3d",i,v,w,v,iord1[v],w,iord1[w]);fflush(stdout);
	//	fprintf(stdout,"\n%%%%102:::h=%d,b=%d,f=%d,i=%3d,h=%3d, flag(h)=%3d, cell(v)=%3d,cell[w]=%3d\n",head[102],back[102],flag[102],i,h,flag[h],cell[v],cell[w]);fflush(stdout);
	// if h is an old set, then we create a new set (c=c+1)
	if(flag[h]<0) {//
	  head[c] = head[h];// insert c between h and head[h]
	  head[h] = c;
	  back[head[c]] = c;
	  back[c] = h;
	  flag[c] = 0;
	  next[c] = -1;
	  // add the new set to fixlst:
	  fixlst[nfx] = c;
	  //	  fprintf(stdout,"\n%%%%%%h=%3d,c=%3d;fixlst[%3d]=%3d",h,c,nfx,fixlst[nfx]);fflush(stdout);
	  nfx++;
	  if(nfx>nmem) fixlst=realloc(fixlst,nfx*sizeof(INT));
	  h=c;
	  c++;
	  if(c>=nmem){
	    nmem=c+1;
	    //	    fprintf(stdout,"\n1nmem=%d,c=%d,n=%3d,ratio=%f\n",nmem,c,n,((REAL )c)/((REAL )n)); fflush(stdout);
	    flag=realloc(flag,nmem*sizeof(INT));
	    head=realloc(head,nmem*sizeof(INT));
	    next=realloc(next,nmem*sizeof(INT));
	    back=realloc(back,nmem*sizeof(INT));
	    head[c]=-1; back[c]=-1; next[c]=-1; flag[c]=-1;
	  }
	}
	// add cell of w to the new set
	next[cell[w]] = next[h];
	if (next[h] >= 0) {
	  back[next[h]] = cell[w];
	}
	flag[cell[w]] = h;
	back[cell[w]] = h;
	next[h] = cell[w];
      }//	end if
    }   //end for
    /* for(ijk = iend;ijk>ibeg;ijk--){ */
    /*   w = ja[ijk]; */
    /*   fprintf(stdout,"\n%%%%%%cell[%3d]=%3d,",w,cell[w]); */
    /*   fprintf(stdout,"flag(%3d)=%3d,",cell[w],flag[cell[w]]); */
    /*   fprintf(stdout,"head(%3d)=%3d,",cell[w],head[cell[w]]); */
    /*   fprintf(stdout,"next(%3d)=%3d,",cell[w],next[cell[w]]); */
    /*   fprintf(stdout,"back(%3d)=%3d",cell[w],back[cell[w]]); */
    /* } */
    /* fprintf(stdout,"\n%%%%%%vertex v=%3d",v); */
    //    print_cells("%%secon:",nfx,i,cell,flag,head,next,back,fixlst);
    //    nfx--;
    for(i0 = 0;i0< nfx;i0++){
      //      h = fixlst[i0];
      //      flag[h] = -1;
      //      fprintf(stdout,"\nFIXLIST[%d]=%d,flag=%d",i0,fixlst[i0],flag[fixlst[i0]]);
      flag[fixlst[i0]] = -1;
      //      fprintf(stdout,"\nn=%4d;FIXLIST[%d]=%d,flag=%d",n,i0,fixlst[i0],flag[fixlst[i0]]);
      //      fprintf(stdout,"\nc=%7d,n=%7d,ratio=%f\n",c,n,((REAL )c)/((REAL )n)); fflush(stdout);
    }// 
  }  //end do
  free(fixlst);
  free(flag);
  free(head);
  free(next);
  free(back);
  free(cell);
  return;
}
/*********************************************************************/
INT main(INT argc, char **argv)
{
  INT dummy,ish,n,nnz,k,iread;
  FILE *fp   = NULL;
  INT  *ia   = NULL, *ja    = NULL;
  INT  *iord = NULL, *iord1 = NULL, *level = NULL,*tree = NULL;
  dCSRmat *a;
  //
  //  fprintf(stdout,"\n%%%%argc=%3d\n",argc);
  //  for(k=0;k<argc;k++)
  //    fprintf(stdout,"\n%%%%arg(%3d)=\'%s\'",k,argv[k]);
  //  fprintf(stdout,"\n%%%%==============\n");
  if(argc>1){
    //    fp=fopen(argv[1],"r");    
    a=malloc(1*sizeof(dCSRmat));
    // read the matrix in coo format; convert it to csr and this is a:
    dcoo_read_dcsr("graphij0",a);
    n=a->row;
    nnz=a->nnz;
    ia=a->IA;
    ja=a->JA;
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
  /* for(k=0;k<n;k++){ */
  /*   fprintf(stdout,"\n%%%%ia[%3d]=%3d",k,ia[k]);fflush(stdout); */
  /* } */
  /* for(k=0;k<nnz;k++){ */
  /*   fprintf(stdout,"\n%%%%ja[%3d]=%3d",k,ja[k]);fflush(stdout); */
  /* } */
  iord1=calloc(n,sizeof(INT));
  level=calloc(n,sizeof(INT));
  tree=calloc(n,sizeof(INT));
  iCSRmat *blk_dfs=run_dfs(n,ia,ja);
  iord=blk_dfs->JA;
  fprintf(stdout,"\nend of list: blk_dfs->IA[%d]=%d (.eq. %d \?)\n",blk_dfs->row,blk_dfs->IA[blk_dfs->row],blk_dfs->row);
  /* for(i=0;i<blk_dfs->row;i++){ */
  /*   fprintf(stdout,"\nblock[%d]; size=%d:  ",i,blk_dfs->IA[i+1]-blk_dfs->IA[i]); */
  /*   for(k=blk_dfs->IA[i];k<blk_dfs->IA[i+1];k++){  */
  /*     fprintf(stdout,"%d ",blk_dfs->JA[k]); */
  /*   }  */
  /* }   */
  fprintf(stdout,"\nblocks=%d; vertices=%d; vertices(orig)=%d\n",blk_dfs->row,blk_dfs->col,n);

  /*******************************************************************/
  /* ivector *roots=malloc(1*sizeof(ivector)); */
  /* // grab a point in every connected component: */
  /* roots->row=blk_dfs->row; */
  /* roots->val=(INT *)calloc(roots->row,sizeof(INT)); */
  /* for(i=0;i<blk_dfs->row;i++){ */
  /*   k=blk_dfs->IA[i]; */
  /*   roots->val[i]=blk_dfs->JA[k]; */
  /* } */
  lexi(n,ia,ja,iord, iord1,level,tree);
  //
  fprintf(stdout,"\np=[");
  for(k=0;k<n;k++){
    fprintf(stdout,"%3d ",iord[k]+1);
  }
  fprintf(stdout,"];\n");
  fprintf(stdout,"\npinv=[");
  for(k=0;k<n;k++){
    fprintf(stdout,"%3d ",iord1[k]+1);
  }
  fprintf(stdout,"];\n");
  fprintf(stdout,"\nlevel=[");
  for(k=0;k<n;k++){
    fprintf(stdout,"%3d ",level[iord[k]]+1);
  }
  fprintf(stdout,"];\n");
  fprintf(stdout,"\ntree=[");
  for(k=0;k<n;k++){
    fprintf(stdout,"%3d ",tree[iord[k]]+1);
  }
  fprintf(stdout,"];\n");
  icsr_free(blk_dfs);
  free(level);
  free(iord1);
  if(argc>1){
    dcsr_free(a);
  }else{
    free(ia);
    free(ja);
  }
    return 0;
}
//
