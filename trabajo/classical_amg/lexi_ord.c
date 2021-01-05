#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "hazmath.h"
/*********************************************************************/
static void print_cells(INT num,
			INT *cell,		\
			INT *flag,		\
			INT *head,		\
			INT *next,		\
			INT *back)
{
  INT hc,hcnext,f,v,b,c,strt;
  hc=0;
  fprintf(stdout,"\nnumber=%d:",num);
  while(1){
    hcnext=head[hc];
    if(hcnext<0) break;
    strt=next[hcnext];
    fprintf(stdout,"\n\theader_cell=%3d; list_strt=%3d",hcnext,strt);
    fprintf(stdout,"\n");
    while(strt>0){
      v=head[strt];
      c=cell[v];
      f=flag[c];
      b=back[c];
      fprintf(stdout,"%6d:|%3d|%3d|%3d|%3d|;",strt,f,v,next[strt],b);
      strt=next[strt];
    }
    fprintf(stdout,"\n");
    hc=hcnext;
  }
  return;
}
/**********************************************************************/
void lexi(INT n, INT *ia, INT *ja, INT *iord, INT *iord1)
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
  INT n2=2*n;
  // allocate memory;
  INT *mbegin=calloc(5*n2,sizeof(INT));
  for (i=0;i<5*n2;i++) mbegin[i]=-1;
  INT *fixlst=calloc(n2,sizeof(INT));
  for (i=0;i<n2;i++) fixlst[i]=-1;
  // set pointers;
  INT *head   = mbegin;
  INT *back   = head   + n2;
  INT *next   = back   + n2;
  INT *cell   = next   + n2;
  INT *flag   = cell   + n2;
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
  for(iv = 0;iv<n;iv++){
    v = iord[iv]; // vertex label
    head[c] = v;  //head(cell)=element;
    cell[v] = c;  //cell(vertex)=cell.
    next[c-1] = c;
    flag[c] = 1;
    back[c] = c-1;
    c++;
    iord1[v] = -1;
  }
  next[c-1] = -1;
  for (i = n-1;i>=0;i--){
    /*********************************************************************/
    //    if(i<n-3) break;
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
    nfx=0;
    ibeg = ia[v]-1;
    iend = ia[v+1]-1;
    for(ijk = iend;ijk>ibeg;ijk--){
      w = ja[ijk];
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
	/* fprintf(stdout,"\n%%%%i=%3d,v=%3d,w=%3d,iord1(%3d)=%3d,iord1(%3d)=%3d",i,v,w,v,iord1[v],w,iord1[w]);fflush(stdout); */
	fprintf(stdout,"\n%%%%102:::h=%d,b=%d,f=%d,i=%3d,h=%3d, flag(h)=%3d, cell(v)=%3d,cell[w]=%3d\n",\
		head[102],back[102],flag[102],i,h,flag[h],cell[v],cell[w]);fflush(stdout);
	// if h is an old set, then we create a new set (c=c+1)
	if(flag[h]<0) {//
	  head[c] = head[h];// insert c between h and head[h]
	  head[h] = c;
	  back[head[c]] = c;
	  back[c] = h;
	  flag[c] = 0;
	  next[c] = -1;
	  // add the new set to fixlst:
	  nfx++;
	  fixlst[nfx] = c;
	  /* fprintf(stdout,"\n%%%%%%h=%3d,c=%3d;fixlst[%3d]=%3d",h,c,nfx,fixlst[nfx]);fflush(stdout); */
	  h=c;
	  c++;
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
      print_cells(i,cell,flag,head,next,back);
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
    for(i0 = 0;i0< nfx;i0++){
      //      h = fixlst[i0];
      //      flag[h] = -1;
      flag[fixlst[i0]] = -1;
    }// 
  }  //end do
  // free
  free(fixlst);
  free(mbegin);
  return;
}
/*********************************************************************/
INT main(INT argc, char **argv)
{
  INT dummy,ish,n,nnz,k,iread;
  FILE *fp   = NULL;
  INT  *ia   = NULL, *ja    = NULL;
  INT  *iord = NULL, *iord1 = NULL;
  //
  //  fprintf(stdout,"\n%%%%argc=%3d\n",argc);
  //  for(k=0;k<argc;k++)
  //    fprintf(stdout,"\n%%%%arg(%3d)=\'%s\'",k,argv[k]);
  //  fprintf(stdout,"\n%%%%==============\n");
  if(argc>1){
    fp=fopen(argv[1],"r");    
  } else {
    fp=fopen("abc","r");
  }
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
  iord=calloc(n,sizeof(INT));
  iord1=calloc(n,sizeof(INT));
  for(k=0;k<n;k++){
    iord[k] = k;
    iord1[k]=-1;
  }
  lexi(n,ia,ja,iord, iord1);
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
  free(iord);
  free(iord1);
  free(ia);
  free(ja);
  return 0;
}
//
