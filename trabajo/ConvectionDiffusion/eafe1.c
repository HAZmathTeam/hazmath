#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hazmath.h"
//#include "ConvectionDiffusion.h"
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
/*!
 *
 * \fn bernoulli(const REAL z)
 * \brief returns B(z) = z/(exp(z)-1)
 *
 */

static REAL bernoulli1(const REAL z){
  REAL tolb=1e-10,zlarge=37.*2.302;  
  if (fabs(z) < tolb)
    return (1.-z*0.5); // around 0 this is the asymptotic;
  else if(z<zlarge){
    return z/(exp(z)-1.);
    //    return (REAL )(((long double )z)/(expl((long double )z)-1.));
  }
  else //z>tlarge this is zero pretty much
    return 0.;
}

/*!
 *
 * \fn diff_coeff1(REAL dim, REAL *val,REAL* x, REAL t, void *param)
 * \brief returns 1;
 *
 */
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
 * \param scalar_val_d is a scalar valued function -- diffusion coefficient   (should be changed to "matrix_val_d")
 * \param scalar_val_rhs is a scalar valued function: right hand side
 * \param vector_val_ad is a vector valued function giving the advection coefficient
 * \param scalar_val_bndnr is a scalar valued function for evaluating the right hand side of Neumann or Robin boundary condition
 *
 * \return A              EAFE stiffness CSR matrix
 * \return b              RHS vector
 *
 */

void eafe1(dCSRmat *ain, dvector *rhs,scomplex *sc){
  // assemble the laplacian matrix
  INT dim=sc->n, ns=sc->ns, nv=sc->nv,dim1=dim+1,dim2=dim*dim;
  INT i,j,k,l,jk,jdim,iaa,iab;
  INT is1,node,jd;
  REAL btmp;
  INT edge_flag=-10,issym=1;
  INT *ia=ain->IA, *ja=ain->JA ;
  REAL *a=ain->val;
  //  dvector diag0 = dvec_create(nv);
  //  for (i=0;i<nv;i++) {diag0.val[i]=0.;}
  REAL xi,yi,zi,xj,yj,zj,bte,alpe;
  // over all elements
  INT nrbits=dim*sizeof(REAL),nibits=dim1*sizeof(INT);
  INT is=(INT )(sizeof(REAL)/sizeof(char));
  void *param=calloc(is*dim2,sizeof(char));//random size but enough
  void *wrk=calloc(dim1*dim1*16*is,sizeof(char));//random size but enough
  INT *nop= (INT *)wrk;
  INT *p= nop+dim1;
  is=(INT )(sizeof(INT)/sizeof(char));
  REAL *xs = (REAL *)(wrk+2*dim1*is);
  REAL *xmass=xs+dim1*dim;
  REAL *bs=xmass+dim;
  REAL *bsinv = bs+dim1*dim1;
  REAL *diff=bsinv+dim1*dim1;
  REAL *ad=diff+dim1*dim1;
  REAL *piv=ad+dim1;
  void *wrk1=(void *)piv;
  fprintf(stdout,"\n%%%%=====1Elements: vertices, n=%d, %d, %d\n", ns,nv,dim);
  for(is=0;is<ns;is++){
    //    we take the
    is1=is*dim1;
    memcpy(nop,(sc->nodes+is1),nibits); // local numbering   
    memset(xmass,0,nrbits);
    //    fprintf(stdout,"\n s=%d; nop=",is);
    for (j=0;j<dim1;j++){
      node=nop[j]; //
      //      fprintf(stdout," %d ",node);
      // xs are coords: 
      memcpy((xs+j*dim),(sc->x+node*dim),nrbits);
    }
    for(i=0;i<dim;i++){
      for (j=0;j<dim1;j++){
	xmass[i]+=xs[j*dim+i];
      }
      xmass[i]/=((REAL )dim1);
    }
    //    now we calculate the advection vector, the diffusion matrix at the center:
    for(j=1;j<dim1;j++){
      jd=j*dim;
      for(i=0;i<dim;i++){
	//	bs[jd-dim+i]=xs[jd+i]-xs[i];
	bs[j-1+i*dim]=xs[jd+i]-xs[i]; // this is B,not BT
      }
    }
    invfull(bsinv,dim,bs,wrk1);
    //
    print_full_mat(dim1,dim,xs,"xs");
    print_full_mat(1,dim,xmass,"xmass");
    print_full_mat(dim,dim,bs,"bs");
    // if b is m by n, by rows we have bij stored in i*n+j, place i=1:m-1;j=1:n-1
    for(j=0;j<dim;j++){
      btmp=0.;
      for(i=0;i<dim;i++){
	//	bsinv[dim2+j]-=bsinv[i*dim+j];
	btmp+=bsinv[i*dim+j];
      }
      /*         dim2=dim*dim;*/
      bsinv[dim2+j]=-btmp;
    }
    // now we need to multiply 
    *((INT *)param)=sc->flags[is];
    diffusion_coeff(diff,xmass,0.,param);
    advection(ad,xmass,0.,param);
    // dummy  to make diff somewhat different
    for(i = 0;i<dim;i++){
      for(j = 0;j<dim;j++){
    	btmp=0.;
    	for (k=0;k<dim;k++){
    	  btmp+=bs[i*dim+k]*bs[j*dim+k];
    	}
    	diff[i*dim+j]=btmp;
      }
      diff[i*dim+i]+=1.;
    }
    print_full_mat(dim,dim,diff,"diff");
    // inv(B)*a*inv(B')
    if(issym){
      /* SYMMETRIC*/
      fprintf(stdout,"\n%% symmetric\n");
      for(i = 0;i<dim1;i++){
	for(j = i;j<dim1;j++){
	  btmp=0.;
	  for (k=0;k<dim;k++){
	    for (l=0;l<dim;l++){
	      btmp+=bsinv[i*dim+k]*diff[k*dim+l]*bsinv[j*dim+l];
	    }
	  }
	  piv[i*dim1+j]=btmp;
	}
      }
      for(i = 1;i<dim1;i++)
	for(j = 0;j<i;j++) piv[i*dim1+j]=piv[i+j*dim1];
    } else {
      /*NON-SYMMETRIC*/
      for(i = 0;i<dim1;i++){
	for(j = 0;j<dim1;j++){
	  btmp=0.;
	  for (k=0;k<dim;k++){
	    for (l=0;l<dim;l++){
	      btmp+=bsinv[i*dim+k]*diff[k*dim+l]*bsinv[j*dim+l];
	    }
	  }
	  piv[i*dim1+j]=btmp;
	}
      }
    }
    /* print_full_mat(1,dim,ad,"ad"); */
    /* print_full_mat(dim1,dim,bsinv,"bsinv"); */
    /* print_full_mat(dim1,dim1,piv,"diffnew"); */
    // diff is in LU below, so destroyed. and piv too. 
    solve_pivot(1, dim, diff, ad, p, piv);
     /*    print_full_mat(1,dim,ad,"ad1"); */
  }
  /* matrices c = a*b+c; a is m by n, b is n by p; c is m by p */
  //
  fprintf(stdout,"\n%%%% end of element loop...\n");  
  //  dvec_free(&diag0);
  //  dvec_free(&dmass);
  return;
}

