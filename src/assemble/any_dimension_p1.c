/*! \file src/amr/any_dimension_p1.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note containing all essential routines that assemble the
 *  stiffness matrix for the Laplace equation and the mass matrix in
 *  any spatial dimensions (1,2,3,...) on any simplicial grid
 *
 *  \note: modified by ltz on 20210531
 *
 */
/*********************************************************************/
#include "hazmath.h"
/*********************************************************************/
/*!
 * \fn grad_compute(INT dim, REAL factorial, REAL *xs, REAL *grad, 
 *                  REAL *vols,void *wrk);
 *
 * \brief Computes the (d+1 x d) matrix [B^{-T} ; -B^{-T}*ones(dim,1)]=grad(barycentric).
 *        Here B is the mapping such that \hat{x} = B^{-1}(x-x_0) 
 *        it also computes the volume of the simplex, which is det(B)/d_factorial.
 *
 * \param dim        I: The dimension of the problem.
 * \param factorial  I: dim! (dim factorial)
 * \param xs         I: coordinates of the simplex
 * \param grad       O: (dim+1) x dim matrix of the gradients of the dim+1
 *                      barycentric coordinates.
 * \param vols       O: The volume of the simplex
 * \param wrk        W: working array of dimension 
 *                      (dim+1)*(dim*sizeof(REAL) + sizeof(INT))
 *
 * \return
 *
 * \note
 */
SHORT grad_compute(INT dim, REAL factorial, REAL *xs, REAL *grad,	\
		   void *wrk)
{
  INT dim1 = dim+1,i,j,j1,ln,ln1;
  REAL vols;
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
  if(ddense_lu(1, dim, &vols, bt,p,piv)) {
    //print_full_mat(dim,dim,bt,"bt");
    //    vols=0.;
    return 1; // this is a degenerate simplex
  }// else
   // vols=fabs(vols)/factorial;
  // in this setting grad is (dim+1) x dim matrix, one barycentric gradient per row???
  memset(grad,0,dim1*dim*sizeof(REAL));
  for (j = 1; j<dim1;j++){
    j1=j-1;
    ln=j*dim;
    grad[ln+j1]=1.;
    // compute B^{-T}*"a column(identity)"
    ddense_solve_pivot(0, dim, bt, (grad+ln), p,piv);
    for(i=0;i<dim;i++){
      //This below is:      grad[0*dim+i]-=grad[ln+i];
      grad[i]-=grad[ln+i];
    }
  }
  //  print_full_mat(dim1,dim,grad,"gradgrad22");
  return 0;
}
/*************************************************************************/
/*!
 * \fn static void local_sm(REAL *slocal,REAL *grad,
 *                          const INT dim,const REAL vols)
 *
 * \brief Computes (grad lambda_i,grad lambda_j) on a simplex
 *
 * \param slocal  O: (dim+1) x (dim+1) local matrix corresponding to
 *                  (grad lamnbda_i,grad lambda_j)
 * \param grad    I: (dim+1) x dim matrix of the gradients of the dim+1
 *                   barycentric coordinates.
 * \param dim     I: The dimension of the problem.
 * \param vols    I: The volume of the simplex
 *                   (dim+1)*(dim*sizeof(REAL) + sizeof(INT))
 *
 * \return
 *
 * \note
 */
void local_sm(REAL *slocal,		\
	      REAL *grad,			\
	      const INT dim,			\
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
/*!
 * \fn static REAL *local_mm(const INT dim)
 *
 * \brief Computes (lambda_i,lambda_j)/(volume) on the reference
 *        element. To compute the mass matrix on any other element one
 *        can use this and multiply by the volume of the physical
 *        element
 *
 * \param dim     I: The dimension of the problem.
 *
 * \return       The mass matrix on the reference element. 
 *
 * \note
 */
//
REAL *local_mm(const INT dim)
{
  INT dim1=dim+1,j,k;
  // all integrals divided by the volume(reference_simplex);
  REAL dd,smij,smii;
  dd=(REAL )(dim1);
  smij=1e0/(dd+1e0)/dd;
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
/*!
 * \fn void local_coords(const INT dim,REAL *xs,INT *nodes,REAL *x)
 *
 * \brief grabs the coordinates of the vertices of a simplex and
 *        stores them in an array xs[(dim+1)*dim].
 *
 * \param dim     I: The dimension of the problem.
 * \param xs      O: (dim+1) x dim matrix with coordinates of the
 *                   vertices of the simplex stored by rows.
 *
 * \param nodes I: (dim+1) vector containing the global numbers (as
 *                 vertices of the mesh) of the vertices of the
 *                 simplex. This is just a row of the
 *                 element_to_vertex matrix.
 *
 * \param x     I: coordinates of the vertices of the mesh
 *
 * \return
 *
 * \note
 */
void local_coords(const INT dim,REAL *xs,INT *nodes, REAL *x)
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
/*************************************************************************/
/*!
 * \fn void assemble_p1(scomplex *sc, dCSRmat *A, dCSRmat *M)
 *
 * \brief Assembles the stiffness and the mass matrix for linear
 *        element in any spatial dimension. 
 *
 * \param sc      I: a simplical complex (a mesh)
 * \param A       O: the assembled stiffness matrix.
 * \param M       O: The assembled mass matrix
 *
 * \return
 *
 * \note
 */
INT assemble_p1(scomplex *sc, dCSRmat *A, dCSRmat *M)
{
  // read the mesh on the cube or square:
  INT i,j,k,idim1,jdim1;  // loop and working
  INT ns,nv,nnz,nnzp; // num simplices, vertices, num nonzeroes
  REAL volume,fact; //mass matrix entries and dim factorial.  
  // for simplices: number of vertices per simplex. 
  INT dim=0,dim1=1;
  dim=sc->n;
  dim1=sc->n+1;  
  nv=sc->nv;// shorthand for num vertices. 
  ns=sc->ns; // shorthand for num simplices
  /*=====================================================*/
  fact=sc->factorial;
  //
  // local mass matrix: it is a constant matrix times the volume of an element.
  REAL *mlocal=local_mm(dim);
  //  fprintf(stdout,"\nnum_simplices=%d ; num_vertices=%d",ns,nv);fflush(stdout);
  // to compute the volume and to grab the local coordinates of the vertices in the simplex we need some work space
  REAL *slocal=calloc(dim1*dim1,sizeof(REAL));// local stiffness matrix.
  REAL *grad=calloc(dim1*dim,sizeof(REAL));// local matrix with
					   // gradients from which
					   // many things can be
					   // computed.
  REAL *xs=calloc(dim*dim1,sizeof(REAL));// for the local coordinates of vertices of a simplex;
  void *wrk=calloc(dim1*dim1,sizeof(REAL));// this is used in every simplex but is allocated only once.
  ////////////////////// ASSEMBLY BEGINS HERE:
  // create a block diagonal mass and stiffness matrices with the local matrices on the diagonal
  dCSRmat *m_dg=malloc(sizeof(dCSRmat));
  m_dg[0]=dcsr_create(ns*dim1,ns*dim1,ns*dim1*dim1);
  /*
    Here we use the fact that we can use the asme IA and JA for the stiffness and mass matrices.
  */
  dCSRmat *a_dg=malloc(sizeof(dCSRmat));
  a_dg->row=m_dg->row;
  a_dg->col=m_dg->col;
  a_dg->nnz=m_dg->nnz;
  a_dg->IA=m_dg->IA;// same as the mass matrix;
  a_dg->JA=m_dg->JA;// same as the mass matrix;
  // needs to be freed at the end
  a_dg->val=calloc(a_dg->nnz,sizeof(REAL));
  //
  nnz=0;
  m_dg->IA[0]=nnz;  
  for(i=0;i<ns;++i){
    idim1=i*dim1;
    // grab the vertex coordinates and compute the volume of the simplex.
    local_coords(dim,xs,&sc->nodes[idim1],sc->x);// sc->nodes[idim1:idim1+dim]
						 // point to global
						 // vertex numbers in
						 // element i
    // compute gradients
    grad_compute(dim, fact, xs, grad,wrk);
    //    print_full_mat(dim1,dim,grad,"grad");fflush(stdout);
    volume=sc->vols[i];
    //    sc->vols[i]=volume;
    //    fprintf(stdout,"\nvolume[%d]=%.5e",i,sc->vols[i]);fflush(stdout);
    // copute local stiffness matrix as grad*transpose(grad);
    local_sm(slocal,grad,dim,volume);
    for(j=0;j<dim1;j++){
      jdim1=j*dim1;
      for(k=0;k<dim1;k++){
	m_dg->JA[nnz]=idim1+k;
	m_dg->val[nnz]=mlocal[jdim1+k]*volume;
	a_dg->val[nnz]=slocal[jdim1+k];
	nnz++;
      }
      m_dg->IA[idim1+j+1]=nnz;
    }
  }
  //
  // create a sparse matrix representing the natural inclusion of
  // P1-continuous in P1 discontinuous
  dCSRmat *P=malloc(sizeof(dCSRmat));
  P[0]=dcsr_create(ns*dim1,nv,ns*dim1);
  //P needs no val, as all nonzeroes are 1. so we free the val;
  free(P->val);  P->val=NULL;//
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
  //  fprintf(stdout,"\n P-check nonzeroes(a=b):%d %d\n",P->nnz,P->IA[m_dg->row]);fflush(stdout);
  //  So now we have P and m_dg, the assembly is P^T*m_dg*P;
  dCSRmat *PT=malloc(sizeof(dCSRmat));
  dcsr_trans(P,PT);// transpose P;
  //  Here PT->val should be NULL;
  dcsr_rap_agg(PT,m_dg,P,M);
  dcsr_rap_agg(PT,a_dg,P,A);
  // free the element matrices
  if(m_dg) {
    dcsr_free(m_dg);
    free(m_dg);
  }  
  if(a_dg->val) free(a_dg->val);  
  if(a_dg) free(a_dg);
  // free all the rest;
  free(xs);
  free(slocal);
  free(grad);
  free(wrk);
  free(mlocal);
  if(P) {
    dcsr_free(P);free(P);
  }
  if(PT){
    dcsr_free(PT);free(PT);
  }
  return 0;
}
/****************************************************************************/
/*END OF FILE*/
