/**************************************************************************/
#include "hazmath.h"
/*************************************************************************/
//
static REAL *local_m(const INT dim)
{
  INT dim1=dim+1,j,k;
  // all integrals divided by the volume(reference_simplex);
  REAL dd,smij,smii;
  dd=(REAL )(dim1);
  smij=1e0/(dd*dd);
  smii=2e0*smij;
  REAL *mlocal=calloc(dim1*dim1,sizeof(REAL));// local matrix stored by rows.
  for(j=0;j<dim1;j++)  {
    for(k=0;k<dim1;k++){
      if(k==j) continue;
      mlocal[j*dim1+k]=smij;
    }
    mlocal[j*dim1+j]=smii;
  }
  return mlocal;
}
/*******************************************************************************/
static void local_coords(const INT dim,REAL *xs,INT *nodes, REAL *x)
{
  // from global coords x and nodes=element-vertex for a particular element
  // correspondence extract the local coordinates from x in xs.
  // xs must be allocated earlier and should have space for dim1*dim REAL;
  INT dim1=dim+1,j,k;
  for(j=0;j<dim1;j++){
    k=nodes[j]; // pick a node in the simplex
    memcpy(&xs[j*dim],&x[k*dim],dim*sizeof(REAL)); // copy coordinates;
  }
  return;
}
/*******************************************************************************/
void main()
{
  // read the mesh on the cube or square:
  INT dim=3;
  INT i,j,k,ijdim1,idim1,jdim1,iaij;  // loop and working
  INT ns,nv,nnz,nnzp; // num simplices, vertices, num nonzeroes
  REAL volume,smij,smii,fact; //mass matrix entries and factorial.  
  char *meshfile=NULL;
  switch(dim){
  case 3:
    meshfile=strdup("../../examples/grids/3D/unitCUBE_n65.haz");
    break;
  default:
    meshfile=strdup("../../examples/grids/2D/unitSQ_n1025.haz");
  }
  FILE *fp=fopen(meshfile,"r");
  fprintf(stdout,"\nReading mesh file(dim=%d)...",dim);fflush(stdout);
  scomplex *sc=haz_scomplex_read(fp,0);// 0 is the print level; reads the mesh
  fprintf(stdout,"done;\n");fflush(stdout);
  // for simplices: number of vertices per simplex. 
  nv=sc->nv;// short hand for num vertices. 
  ns=sc->ns; // short hand for num simplices
  INT dim1=sc->n+1; //number of vertices in a simplex. 
  /*=====================================================*/
  // compute dim!=dim_factorial (needed for the volume):  
  fact=1e0;
  for (j=2;j<dim1;j++) fact *= ((REAL )j);
  //
  // local mass matrix: it is constant times the volume of an element.
  REAL *mlocal=local_m(dim);
  fprintf(stdout,"\nnum simplices=%d ; num vertices=%d",ns,nv);fflush(stdout);
  // to compute the volume and to grab the local coordinates of the vertices in the simplex we need some work space
  REAL *xs=calloc(dim*dim1,sizeof(REAL));// for the local coordinates of vertices of a simplex;
  void *wrk=calloc(dim1*dim1,sizeof(REAL));// this is used in every simplex but is allocated only once.
  ////////////////////// ASSEMBLY BEGINS HERE:
  clock_t clk_assembly_start = clock(); // End of timing for mesh and FE setup
  // create a "big" block diagonal matrix with the local mass matrices on the diagonal
  dCSRmat *Mbig=malloc(sizeof(dCSRmat));
  Mbig[0]=dcsr_create(ns*dim1,ns*dim1,ns*dim1*dim1);
  nnz=0;
  Mbig->IA[0]=nnz;  
  for(i=0;i<ns;++i){
    idim1=i*dim1;
    // we now grab the vertex coordinates and then compute the volume of the simplex.
    local_coords(dim,xs,&sc->nodes[idim1],sc->x);// sc->nodes[idim1:idim1+dim]
						 // point to global
						 // vertex numbers in
						 // element i
    // find volume of the simplex:
    vol_simplex(dim,fact, xs, &volume, wrk);
    //    fprintf(stdout,"\nelem=%d ; volume=%e",i,volume);fflush(stdout);
    for(j=0;j<dim1;j++){
      jdim1=j*dim1;
      ijdim1=idim1+jdim1;
      for(k=0;k<dim1;k++){
	//	fprintf(stdout,"\nrow=%d ; col=%d",idim1+j,idim1+k);fflush(stdout);
	Mbig->JA[nnz]=idim1+k;
	Mbig->val[nnz]=mlocal[jdim1+k]*volume;
	nnz++;
      }
      //      fprintf(stdout,"\nnext_row=%d; nnz=%d",idim1+j+1,nnz);
      Mbig->IA[idim1+j+1]=nnz;
    }
  }
  //  fprintf(stdout,"\n check nonzeroes(nnz1(Mbig)=nnz2(MBig))?:%d %d\n",Mbig->nnz,Mbig->IA[Mbig->row]);
  // create a sparse matrix for the natural inclusion of P1-continuous in P1 discontinuous
  dCSRmat *P=malloc(sizeof(dCSRmat));
  P[0]=dcsr_create(ns*dim1,nv,ns*dim1);//P needs no val, as all nonzeroes are 1. so we free the val;
  free(P->val);  P->val=NULL;// we do not need the val.
  nnzp=0;
  P->IA[0]=nnzp;
  for(i=0;i<ns;++i){
    idim1=i*dim1;
    k=idim1;// starting point for the element->vertex map;
    for(j=0;j<dim1;j++){
      P->JA[P->IA[k]]=sc->nodes[k];
      k++; nnzp++;
      P->IA[k]=nnzp;
    }
  }
  P->nnz=nnzp;  
  //  fprintf(stdout,"\n P-check nonzeroes(a=b):%d %d\n",P->nnz,P->IA[Mbig->row]);fflush(stdout);
  //  So now we have P and Mbig, the assembly is P^T*Mbig*P;
  dCSRmat *PT=malloc(sizeof(dCSRmat));
  dcsr_trans(P,PT);// transpose P;
  //  Here PT->val should be NULL;
  //  csr_print_matlab(stdout,Mbig);
  dCSRmat *M=malloc(sizeof(dCSRmat));
  dcsr_rap_agg(PT,Mbig,P,M);
  clock_t clk_assembly_end = clock(); // End of timing for mesh and FE setup
  printf("\n\nelapsed CPU time for assembling the global mass matrix = %f seconds.\n\n",
         (REAL) (clk_assembly_end - clk_assembly_start)/CLOCKS_PER_SEC);

  // no need of Mbig or anything like that;
  dcsr_free(Mbig);
  // the rest of dcsr may be needed. 

  // free all;
  free(xs);
  free(wrk);
  haz_scomplex_free(sc);  
  dcsr_free(P);
  dcsr_free(PT);
  dcsr_free(M);  
  return;
}
