/*******************LIBLIB*********************/
typedef struct /* structure to support splitting unit cube into simplices */
{ INT n; /* spatial dimension of the grid */
  INT nvcube; /* number of vertices on the unit cube in R^n=2^{n}.*/
  INT ns; /* number of n dimensional simplices in the unit
	     cube(n_factorial of them) */
  INT ne; // number of edges in the cube.
  INT nf; // number of n-1 dimensional faces in the cube.
  unsigned INT *bits; /* the binary digits of all the integers from 0
		to 2^{n-1} as an array. These are also the coordinates
		of the vertices of the unit cube in R^n*/
  INT *nodes; /* the array describing each of the n factorial simplices
		in the unit cube */
  INT *perms; /* the n by nvcube array describing the permutations
		which give different orderings of the cube's vertices
		so that we get a criss-cross grid in n
		dimensions. There are n such permutations reflecting
		the cube across the hyperplanes meeting at the vertex
		with coordinates [1,1,1,1,...,1] perms[nodes] will
		give consistent splitting of neighboring cubes.
	     */
} cube2simp;

/********************************LIBLIB */
void coord_lattice(INT *m,const INT dim,			\
		   const INT kf, const INT nall, const INT *nd)
{
  INT i,j,k;
  /*
    given a global number kf on a lattice grid with lexicographical
    ordering, this returns the n-tuple of latice coordinates,
    m[i],i=0:n-1). nall is the total number of vertices in the
    lattice, nd[i] is the number of divisions in every direction;
  */
  j=nall;  k=kf;
  for(i=dim;i>0;i--){
    j=j/(nd[i-1]+1);  m[i-1] = k/j;   k=k-m[i-1]*j;
  }
  return;
}
INT num_lattice(INT *m,const INT dim,INT *nd)
{
  INT i,j,k,kf;
  /*
    given the lattice coordinates m[i],i=0:n-1 of a vertex, in
    dimension n, finds the global number kf of a vertex on a
    lexicographically ordered lattice grid. divisions in each
    directions are stored in nd[].
  */
  kf=m[dim-1];
  //      kf = m[dim-1];
  for (i=dim-1; i>0; i--){
    kf=kf*(nd[i-1]+1)+m[i-1];
  }
  return kf;
}
INT reverse(INT *arr,INT length,size_t elsize)
{
  int i,j,k,swap,nnn=(INT)(length/2);
  //  reverses ordering in an INT array;
  k=length-1;
  for(i=0;i<nnn;i++){
    swap=arr[i];
    arr[i]=arr[k];
    arr[k]=swap;
    k--;
  }  
}
void coord_lattice1(INT *m,const INT dim,			\
		   const INT kf, const INT nall, const INT *nd)
{
  INT i,j,k,swap;
  /*
    given a global number kf on a lattice grid with lexicographical
    ordering, this returns the n-tuple of latice coordinates,
    m[i],i=0:n-1). nall is the total number of vertices in the
    lattice, nd[i] is the number of divisions in every direction;
  */
  j=nall;  k=kf;
  for(i=dim;i>0;i--){
    j=j/(nd[i-1]+1);  m[i-1] = k/j;   k=k-m[i-1]*j;
  }
  return;
}
INT num_lattice1(INT *m,const INT dim,INT *nd)
{
  INT i,j,k,kf;
  /*
    given the lattice coordinates m[i],i=0:n-1 of a vertex, in
    dimension n, finds the global number kf of a vertex on a
    lexicographically ordered lattice grid. divisions in each
    directions are stored in nd[].
  */
  kf=m[dim-1];
  //      kf = m[dim-1];
  for (i=dim-1; i>0; i--){
    kf=kf*(nd[i-1]+1)+m[i-1];
  }
  INT kf1 = m[2]*(nd[0]+1)*(nd[1]+1) + m[1]*(nd[0]+1) + m[0];
  //  fprintf(stdout,"\nkf1=%d,kf=%d",kf,kf1);
  return kf1;
}
    //    fprintf(stdout,"\n"); fflush(stdout);
    /*PRINTING ONLY HERE*/
    /* fprintf(stdout,"\nkf=%d; [",kf);fflush(stdout); */
    /* for(i=0;i<dim-1;i++){ */
    /*   fprintf(stdout,"%d,",m[i]+1); */
    /* } */
    /* fprintf(stdout,"%d]; type =%d; (",m[dim-1]+1,type+1); */
    /* for(j=0;j<c2s->nvcube;j++){ */
    /*   for(i=0;i<dim;i++){ */
    /* 	mm[i]=m[i]+(c2s->bits[dim*j+i]); */
    /*   } */
    /*   k=num_lattice(mm,dim,nd); */
    /*   fprintf(stdout," %d ",k+1); */
    /* } */
    /* fprintf(stdout,")");fflush(stdout); */
    /* fprintf(stdout,"\nkf=%d;type=%d s=",kf,type+1); */
    /* for(i=0;i<c2s->ns;i++){ */
    /*   fprintf(stdout,"["); */
    /*   for(j=0;j<dim;j++){ */
    /* 	fprintf(stdout,"%d,",cnodes[c2s->nodes[i*dim1+j]]+1); */
    /*   } */
    /*   fprintf(stdout,"%d]; ",cnodes[c2s->nodes[i*dim1+dim]]+1); */
    /* } */
    /**/
    /*   for(i=0;i<dim-1;i++){ */
    /* 	fprintf(stdout,"%.3f,",x[k*dim+i]); */
    /*   } */
    /*   fprintf(stdout,"%.3f)",x[(k+1)*dim-1]); */
    /* } */
    /* fprintf(stdout,"\n"); */
