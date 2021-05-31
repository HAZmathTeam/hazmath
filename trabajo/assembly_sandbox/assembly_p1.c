/**************************************************************************/
#include "hazmath.h"
/*************************************************************************/
SHORT grad_compute(INT dim, REAL factorial, REAL *xs, REAL *grad,	\
		   REAL *vols,						\
		   void *wrk)
{
  /*
    computes the (d+1 x d) matrix [B^{-T} ; -B^{-T}*ones(dim,1)]=grad(barycentric).
    Here B is the mapping such that \hat{x} = B^{-1}(x-x_0) 
    it also computes the volume of the simplex, which is det(B)/d_factorial.
    work space: wrk should be at least dim*(dim+1) REALS and dim
    integers.
  */
  INT dim1 = dim+1,i,j,j1,ln,ln1;
  REAL *bt=(REAL *)wrk;
  REAL *piv=bt+dim*dim;
  INT *p = (INT *)(wrk+(dim*dim + dim)*sizeof(REAL));
  // construct bt using xs;
  for (j = 1;j<dim1;j++){
    ln=j*dim; ln1=ln-dim;
    for(i=0;i<dim;i++){
      bt[ln1+i] = xs[ln+i]-xs[i];  // k-th row of bt is [x(k)-x(0)]^T. x(k) are coords of vertex k. 
    }
  }
  //  print_full_mat(dim,dim,bt,"bt");
  //  print_full_mat(dim,1,piv,"piv");
  // lu decomposition of bt
  if(lufull(1, dim, vols, bt,p,piv)) {
    //print_full_mat(dim,dim,bt,"bt");
    vols[0]=0.;
    return 1; // this is a degenerate simplex
  } else
    vols[0]=fabs(vols[0])/factorial;
  // in this setting grad is (dim+1) x dim matrix, one barycentric gradient per row???
  memset(grad,0,dim1*dim*sizeof(REAL));
  for (j = 1; j<dim1;j++){
    j1=j-1;
    ln=j*dim;
    grad[ln+j1]=1.;
    // compute B^{-T}*"a column(identity)"
    solve_pivot(0, dim, bt, (grad+ln), p,piv);
    for(i=0;i<dim;i++){
      //This below is:      grad[0*dim+i]-=grad[ln+i];
      grad[i]-=grad[ln+i];
    }
  }
  //  print_full_mat(dim1,dim,grad,"gradgrad22");
  return 0;
}
/*************************************************************************/
//
static void local_sm(REAL *slocal,		\
		     REAL *grad,		\
		     const INT dim,		\
		     const REAL vols)
{
  // grad is a (dim+1) x (dim) gradients of the barycentric coords.
  // slocal(dim,dim) must be alllocated before entering here
  INT dim1=dim+1,i,j,k,idim,idim1,jdim;
  REAL s;
  // peforms grd*transpose(grad) which is d+1 by d+1 matrix
  for(i=0;i<dim1;++i){
    idim=i*dim;
    idim1=i*dim1;
    for(j=0;j<dim1;++j){
      s=0e0;
      jdim=j*dim;
      for(k=0;k<dim;++k){
	//sum_k grad_{k,i}*grad_{k,j}
	  s+=grad[idim+k]*grad[jdim+k];
      }
      //      fprintf(stdout,"\ns123123=%f",s);fflush(stdout);
      slocal[idim1+j]=s*vols;
    }
  }
  return;
}
/*************************************************************************/
//
static REAL *local_mm(const INT dim)
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
  INT dim=2;
  INT i,j,k,ijdim1,idim1,jdim1,iaij;  // loop and working
  INT ns,nv,nnz,nnzp; // num simplices, vertices, num nonzeroes
  REAL volume,smij,smii,fact; //mass matrix entries and factorial.  
  char *meshfile=NULL;
  switch(dim){
  case 3:
    meshfile=strdup("../../examples/grids/3D/unitCUBE_n3.haz");
    break;
  default:
    meshfile=strdup("../../examples/grids/2D/unitSQ_n5.haz");
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
  REAL *mlocal=local_mm(dim);
  fprintf(stdout,"\nnum simplices=%d ; num vertices=%d",ns,nv);fflush(stdout);
  // to compute the volume and to grab the local coordinates of the vertices in the simplex we need some work space
  REAL *slocal=calloc(dim1*dim1,sizeof(REAL));// local stiffness matrix.
  REAL *grad=calloc(dim1*dim,sizeof(REAL));// local stiffness matrix.
  REAL *xs=calloc(dim*dim1,sizeof(REAL));// for the local coordinates of vertices of a simplex;
  void *wrk=calloc(dim1*dim1,sizeof(REAL));// this is used in every simplex but is allocated only once.
  ////////////////////// ASSEMBLY BEGINS HERE:
  clock_t clk_assembly_start = clock(); // begin assembly timing;
  // create a "big" block diagonal mass and stiffness matrices with the local matrices on the diagonal
  dCSRmat *Mbig=malloc(sizeof(dCSRmat));
  dCSRmat *Abig=malloc(sizeof(dCSRmat));// we are saving here, we will
					// use the IA, JA from mass
					// and stiffness, so stiffness
					// will only have a different
					// val
  Mbig[0]=dcsr_create(ns*dim1,ns*dim1,ns*dim1*dim1);
  //
  Abig->row=Mbig->row;
  Abig->col=Mbig->col;
  Abig->nnz=Mbig->nnz;
  Abig->IA=Mbig->IA;
  Abig->JA=Mbig->JA;
  Abig->val=calloc(Abig->nnz,sizeof(REAL));// needs to be freed at the end
  //
  // Now action:
  nnz=0;
  Mbig->IA[0]=nnz;  
  for(i=0;i<ns;++i){
    idim1=i*dim1;
    // we now grab the vertex coordinates and then compute the volume of the simplex.
    local_coords(dim,xs,&sc->nodes[idim1],sc->x);// sc->nodes[idim1:idim1+dim]
						 // point to global
						 // vertex numbers in
						 // element i
    //    vol_simplex(dim,fact, xs, &volume, wrk);
    grad_compute(dim, fact, xs, grad,		\
		 &volume,			\
		 wrk);
    local_sm(slocal,	\
	     grad,	\
	     dim,	\
	     volume);
    print_full_mat(dim1,dim,slocal,"slocal");
    fprintf(stdout,"\nelem=%d ; volume=%e",i,volume);fflush(stdout);
    for(j=0;j<dim1;j++){
      jdim1=j*dim1;
      ijdim1=idim1+jdim1;
      for(k=0;k<dim1;k++){
	//	fprintf(stdout,"\nrow=%d ; col=%d",idim1+j,idim1+k);fflush(stdout);
	Mbig->JA[nnz]=idim1+k;
	Mbig->val[nnz]=mlocal[jdim1+k]*volume;
	Abig->val[nnz]=slocal[jdim1+k];
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
  dCSRmat *M=malloc(sizeof(dCSRmat));
  dcsr_rap_agg(PT,Mbig,P,M);
  dCSRmat *A=malloc(sizeof(dCSRmat));
  dcsr_rap_agg(PT,Abig,P,A);
  clock_t clk_assembly_end = clock(); // End of timing for mesh and FE setup
  printf("\n\nelapsed CPU time for assembling the global mass matrix = %f seconds.\n\n",
         (REAL) (clk_assembly_end - clk_assembly_start)/CLOCKS_PER_SEC);
  csr_print_matlab(stdout,A);
  // no need of Mbig or anything like that;
  dcsr_free(Mbig);
  free(Abig->val);
  // free all;
  free(xs);
  free(slocal);
  free(grad);
  free(wrk);
  haz_scomplex_free(sc);  
  dcsr_free(P);
  dcsr_free(PT);
  dcsr_free(M);  
  return;
}
