/*! \file src/fem/eafe.c
 *
 *  Created by Xiaozhe Hu, James Adler, and Ludmil Zikatanov 20170309
 *  (originally 19940524)  
 *
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief computes the EAFE FE  stiffness matrix and right hand side for the FE approximation of the equation
 *        -div(a(x)*grad(u) - u*b) = f
 *        b is an advection vector; a(x) is a diffusion matrix. 
 *        u = 0 on the boundary.
 *        assumes that all indices in CSR start from 1. 
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
/*********************************************************************/
/*!
 * \brief Calculates boundary integral using lumped mass.
 *
 * \note For an equation in divergence form, the natural boundary condition is: \grad u . n + (b . n)u=g;
 *
 * \note A boundary condition such as  \grad u . n = g;
  is a Robin type condition when the equation is in divergence form.
 *
 * \note All matrices are assumed to be indexed at 1 in the CSR formatting.
 *
 * \param mesh           Mesh Data
 *
 * \param faketime Physical Time if time dependent (not used now)
 * \param vector_val_ad is a vector valued function giving the advection coefficient
 * \param scalar_val_bndnr is a scalar valued function for evaluating the right hand side of Neumann or Robin boundary condition
 *
 * \note there are some Robin type conditions so we need to compute the addition to the stiffness matrix. We use lumped mass quadrature for this and store the diagonal element of the boundary mass matrix in dmass[*];
 * \return dmass lumped mass matrix for the boundary.
 * \return rhs  modified rhs vector if non-zero right hand side for Neumann/Robin boundary condition
 *
*/
static void LumpMassBndry(const trimesh mesh,
		   void (*vector_val_ad)(REAL *, REAL *, REAL),	\
		   void (*scalar_val_bndnr)(REAL *, REAL *, REAL),	\
		   dvector *rhs, dvector *dmass)
{
#ifndef MARKER_ROBIN
  return;
#endif
#ifndef MARKER_NEUMANN
  return;
#endif
  INT i=-1, j=-1,k=-1,ifa=-1,ifb=-1,attri=-16;
  REAL vol=1e10,bdotn=1e10,rhsi=1e10;
  /* mesh entities: number of faces, number boundary faces and so on */
  INT  nv=mesh.nv, nf=mesh.nface,nfb=mesh.nbface,dim=mesh.dim; 
  REAL ccc=1./((REAL ) dim), ad[dim],xm[dim];
  iCSRmat *f2v=mesh.f_v; /* face to vertex map as iCSR */
  INT *fonb=mesh.f_flag; /* boundary markers for every face */
  /* get areas, normal vectors and barycenters  of faces */
  REAL *fa=mesh.f_area, *fn=mesh.f_norm, *fm=mesh.f_mid;
  dmass->row=0;  dmass->val=NULL;
  fprintf(stdout,"Num Of Bdr faces =%i\n",nfb);
  for (i = 0; i < nf; i++)   {
    attri=fonb[i];
    if(attri < MARKER_NEUMANN || attri >= MARKER_BOUNDARY_NO)
      continue; 
    else if(attri >= MARKER_ROBIN && attri < MARKER_BOUNDARY_NO) {
      if(dmass->row) continue;
      dmass->row=nv;
      continue;
    } else { /* this is now Neumann condition */
      vol=fa[i];//*ccc;
      for(k=0;k<dim;k++)
	xm[k]=fm[i*dim+k];
      scalar_val_bndnr(&rhsi,xm,0.0);
      //      fprintf(stdout,"COORDS:%22.16e, %22.16e ; vol=%22.16e\n",xm[0],xm[1],rhsi);
      rhsi = rhsi*vol*ccc;
      ifa=f2v->IA[i]-1;
      ifb=f2v->IA[i+1]-1;
      for(k=ifa;k<ifb;k++){
	j=f2v->JA[k]-1;
	rhs->val[j] += rhsi;
	//	fprintf(stdout,"i,j:%i,%i, rhs_add:%22.16e,vol=%22.16e; tot=%22.16e\n",i+1,j+1,rhsi/vol,vol,rhsi);
      }
    }
  }
  return;
  if(dmass->row) return;
  dmass->val=(REAL *)calloc(dmass->row,sizeof(REAL));
  for (i = 0; i < nf; i++)   {
    attri=fonb[i];
    if(attri < MARKER_ROBIN || attri >= MARKER_BOUNDARY_NO)
      continue;
    vol=fa[i];//*ccc;
    for(k=0;k<dim;k++)  xm[k]=fm[i*dim+k];
    vector_val_ad(xm,ad,0.0);
    bdotn=0.;
    for(k=0;k<dim;k++) bdotn += ad[k]*fn[i*dim+k];
    bdotn*=vol;
    ifa=f2v->IA[i]-1;
    ifb=f2v->IA[i+1]-1;
    for(k=ifa;k<ifb;k++){
      j=f2v->JA[k]-1;
      dmass->val[j]-= bdotn;
    }
  }
  return;
}

/*!
 *
 * \fn bernoulli(const REAL z)
 * \brief returns B(z) = z/(exp(z)-1)
 *
 */
static REAL bernoulli(const REAL z)
{
  double tolb=1e-10,zlarge=37.*2.302;  
  if (fabs(z) < tolb)
    return (1.-z*0.5); // around 0 this is the asymptotic;
  else if(z<zlarge)
    return (z/(exp(z)-1.));
  else //z>tlarge this is zero pretty much
    return 0.;
}

/*!
 *
 * \fn poisson_coeff(REAL *val,REAL* x, REAL t)
 * \brief returns 1;
 *
 */
static void poisson_coeff(REAL *val,REAL* x, REAL t,void *param) {
  *val = 1.0;
  return;
}

/*!
 * \fn  eafe(dCSRmat *A, dvector *rhs, void (*local_assembly)(REAL *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL), trimesh mesh, fespace FE, qcoordinates *cq, void (*scalar_val_d)(REAL *, REAL *, REAL), void (*scalar_val_rhs)(REAL *, REAL *, REAL), void (*vector_val_ad)(REAL *, REAL *, REAL), void (*scalar_val_bndnr)(REAL *, REAL *, REAL), REAL faketime)
 * \brief Uses Schur product and from the assembled matrix for Poisson
 * equation with natural boundary conditions makes the EAFE FE
 * discretization for the equation
 *
 *        -div(a(x)*grad(u) - u*b) = f
 *
 *        b is an advection vector; a(x) is a diffusion matrix.
 *        u = 0 on the boundary.
 *
 *  \note Details about the EAFE discretization are found in: Jinchao
 *        Xu and Ludmil Zikatanov: A monotone finite element scheme
 *        for convection-diffusion equations. Math. Comp. 68 (1999),
 *        no. 228, 1429–1446.
 * \note All matrices are assumed to be indexed at 1 in the CSR formatting.
 *
 * \param local_assembly Routine to get local matrices
 * \param FE             FE Space
 * \param mesh           Mesh Data
 * \param cq             Quadrature Nodes
 *
 *\param faketime Physical Time if time dependent (not used now)
 * \param scalar_val_d is a scalar valued function -- diffusion coefficient   (should be changed to "matrix_val_d")
 * \param scalar_val_rhs is a scalar valued function: right hand side
 * \param vector_val_ad is a vector valued function giving the advection coefficient
 * \param scalar_val_bndnr is a scalar valued function for evaluating the right hand side of Neumann or Robin boundary condition
 *
 * \return A              EAFE stiffness CSR matrix
 * \return b              RHS vector
 *
 */
void eafe(dCSRmat *A, dvector *rhs,		\
      void (*local_assembly)(REAL *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL,void *),REAL), \
	  trimesh mesh, fespace FE, qcoordinates *cq,	\
      void (*scalar_val_d)(REAL *, REAL *, REAL),			\
      void (*scalar_val_rhs)(REAL *, REAL *, REAL, void *),			\
      void (*vector_val_ad)(REAL *, REAL *, REAL),			\
      void (*scalar_val_bndnr)(REAL *, REAL *, REAL), REAL faketime)
{
  assemble_global(A,rhs,local_assembly,			\
		  &FE,&mesh,cq,				\
		  scalar_val_rhs,poisson_coeff,0.0);
  INT i,j,jk,jdim,iaa,iab,nv=mesh.nv,dim=mesh.dim;
  INT *ia=A->IA, *ja=A->JA ;
  REAL *a=A->val;
  coordinates *xyz=mesh.cv;
  dvector diag0 = dvec_create(nv);
  for (i=0;i<nv;i++) {diag0.val[i]=0.;}
  REAL ad[dim], xm[dim],te[dim];
  REAL xi,yi,zi,xj,yj,zj,bte,alpe;
  for (i = 0; i < nv; i++) {
    xi=xyz->x[i];
    yi=xyz->y[i];
    if(dim>2)
      zi=xyz->z[i];
    iaa=ia[i]-1;
    iab=ia[i+1]-1;
    for (jk=iaa; jk<iab; jk++){
      j=ja[jk]-1;
      /*      fprintf(stdout,"aaaaa %i,%i\n",i,j); */
      if(i != j){
	xj=xyz->x[j];
	yj=xyz->y[j];
	if(dim>2) 
	  zj=xyz->z[j];
	/* 
	   compute the advection field at the middle of the edge and
	   then the bernoulli function
	*/
	te[0]  = xi - xj;
	te[1]  = yi - yj;
	xm[0] = (xi + xj)*0.5e+0;
	xm[1] = (yi + yj)*0.5e+0;
	if(dim>2){
	  te[2]  = zi - zj;
	  xm[2] = (zi + zj)*0.5e+0;
	}
	vector_val_ad(ad,xm,0.0);
	bte=0.;
	for(jdim=0;jdim<dim;jdim++){
	  bte += ad[jdim]*te[jdim];
	}
	scalar_val_d(&alpe,xm,0.0);
	// alpe=a(xmid)\approx harmonic_average=|e|/(int_e 1/a);
	// should be computed by quadrature in general for 1/a(x).
	// a_{ij}=B(beta.t_e/harmonic_average)*harmonic_average*omega_e;
	//	  for (i,j):    B(beta.t_e/harmonic_average);   
	//	  for (j,i):    B(-beta.t_e/harmonic_average), the minus comes from
	//     is because t_e=xi-xj;*/
	a[jk] *= (alpe*(bernoulli(bte/alpe))); 
	//	fprintf(stdout,"\ndd=%22.14e",(alpe*(bernoulli(bte/alpe))));
	/* the diagonal is equal to the negative column sum;*/
	diag0.val[j]-=a[jk];
      }
    } 
  }
  /* Another loop to set up the diagonal equal to whatever it needs to equal. */
  for (i = 0; i < nv; i++){
    iaa=ia[i]-1;
    iab=ia[i+1]-1;
    for (jk=iaa; jk<iab; jk++){
      j = ja[jk]-1; 
      if(i == j) a[jk]=diag0.val[i];
    }
  }
  dvector dmass;  
  LumpMassBndry(mesh,vector_val_ad,scalar_val_bndnr,rhs,&dmass);
  /* If Robin condition: */
  if(dmass.row==nv){
    fprintf(stdout,"\n Robin condition\n");
    for (i = 0; i < nv; i++) {
      for (jk=iaa; jk<iab; jk++){
	j = ja[jk]-1; 
	if(i == j) a[jk]+=dmass.val[i];
      }
    }
  }
  dvec_free(&diag0);
  dvec_free(&dmass);
  return;
}





