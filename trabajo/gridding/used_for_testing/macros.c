#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <getopt.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>
#include "hazmath.h"
#ifndef INT
#define INT int
#endif
#ifndef REAL
#define REAL double
#endif
#ifndef DIM
#define DIM 3
#endif
/********************************FINCTIONS:*/
// form macroelements from segments;
kvertex=9;
dim=2;

xy=[...
    1, 1.5  ,3;
    2, 2.5  ,3;
    3, 2    ,1.5;
    4, 0.75 ,0.5;
    5, 2    ,3.75;
    6, 1    ,3.75;
    7, 2    ,5;
    8, 0.5, 2; 
     9, 1.5,2.1;
     10,0.5,3;
     11, 3.1,1.9...
    ];

seg=[...
        4,3;
	8,9;
	 10,1;
	 3,11;
	 9,2;
	 9,1;
	 1,5;
	 6,7;
	 6,1;
	 7,5;
	 10,8;
	 2,11;
	 9,3;
	 1,9;
	 8,4];
nv=size(xy,1);
ne = size(seg,1);
o=ones(ne,1);
A=sparse(seg(:,1),seg(:,2),o,nv,nv);
[i,j,o]=find(A+transpose(A));
ne1 = length(i);
o=ones(ne1,1);
A=sparse(i,j,o,nv,nv);

hold on;
ek=zeros(nv,1);
for k=1:nv
  ek(k)=1;
  a=A*(A*ek+ek);
  jj=find(a==dim);
  jj=jj(find(jj>k));
  ek(k)=0;
  disp([int2str(k),':  ',int2str(jj')])
end
gplot(A,xy(:,2:3)); 
plot(xy(:,2),xy(:,3),'ro'); 
for i=1:nv
set(text(xy(i,2)+0.01,xy(i,3)+0.01,int2str(i)),'FontSize',24);
end
axis off
hold off 
/*********************************/
INT main(INT argc, char **argv)
{
  INT i=-1,j=-1,k=-1,kperm,dim=DIM;
  cube2simp *c2s=cube2simplex(dim);
  INT ns=c2s->ns,dim1=c2s->n+1;
  INT *nd=(INT *)calloc(dim,sizeof(INT));  
  for(i=0;i<dim;i++) nd[i]=33+i;
  INT memx=dim;
  for(i=0;i<dim;i++) memx*=(nd[i]+1);
  scomplex *sc=(scomplex *)umesh(dim,nd,c2s);
  fprintf(stdout,"\nuniform mesh in dim=%d; vertices: %d, simplexes %d\n",dim,sc->nv,sc->ns);  
  REAL xmac[24]={-2.1,-2.1,-2.1,		\
		 1.,-1.,-1.,			\
		 -1., 1.,-1.,			\
		 1., 1.,-1.,			\
		 -1.,-1., 1.,			\
		 1.,-1., 1.,			\
		 -1., 1., 1.,			\
		 1., 1., 1.};
  
  map2mac(sc,c2s,xmac);
  //  haz_scomplex_print(sc,0,"HAHA");
  if(dim==2||dim==3) {
    vtkw("newmesh.vtu",sc,0,0,1.);
  }
  cube2simp_free(c2s);
  if(nd) free(nd);
  haz_scomplex_free(sc);
  return 0;
}
  
