#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hazmath.h"

void num_assembly(INT ndof, INT *nop,INT *fflags,	\
		  dCSRmat *ain, REAL *f,		\
		  REAL *ae, REAL *fe, INT *ip)
{
  //ndof is the LOCAL number of degrees of freedom in this particular
  //element. then ae is ndof by ndof; ip should be initialized before
  //the first call of this function which is usually donein element by
  //element loop.
// nop here is pointed to by (sc->nop+element_num*ndof)
  INT *ia=ain->IA,*ja=ain->JA;
  INT i,j,k,l,ll,iii,iaa,iab,jk;
  REAL *a=ain->val;
/*-------------------------------------------------------------------- */
/*   Numerical assembly */
/*   Input:   */
/*     AE, FE - element matrix and vector, respectively, to be */
/*-------------------------------------------------------------------- */
  for(l=0;l<ndof;l++){
    i = nop[l];
    if (chk_bc(fflags[i],1)) {
    //      fprintf(stdout,"\ndofs with special bcflag: (i=%d):flag=%d",i,fflags[i]);
      continue;
    }
    k = l - ndof;
    f[i]+=fe[l];
    for (ll=0;ll<ndof;ll++){
      k+=ndof;
      j=nop[ll];
      //      fprintf(stdout,"\n%d %d; flags=%d %d",i,j,fflags[i],fflags[j]);
      if(chk_bc(fflags[j],1)) {
	//	fprintf(stdout,"\ndofs(2) with special bcflag: (%d,%d):(flagi,flagj)=%d,%d",i,j,fflags[i],fflags[j]);
	f[i]-=ae[k]*f[j];
      } else {
	ip[j] = k;
	continue;
      }
    }
    iaa = ia[i];
    iab = ia[i+1];
    iii = 0;
    for(jk=iaa;jk<iab;jk++){
      j=ja[jk];
      k = ip[j];
      //      fprintf(stdout,"\nk=%d:",k);
      if(k < 0) continue;
      //      fprintf(stdout,"\n(i,j)=%d %d; flags=%d %d",i,j,fflags[i],fflags[j]);
      a[jk]+=ae[k];
      ip[j] = -1;
      iii++;
      if(iii==ndof) break;
    }
  }
  return;
}
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
 * \fn  eafe1(dCSRmat *A, dvector *rhs, void (*local_assembly)(REAL *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL), trimesh mesh, fespace FE, qcoordinates *cq, void (*scalar_val_d)(REAL *, REAL *, REAL), void (*scalar_val_rhs)(REAL *, REAL *, REAL), void (*vector_val_ad)(REAL *, REAL *, REAL), void (*scalar_val_bndnr)(REAL *, REAL *, REAL), REAL faketime)
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
  INT is1,node,jd,id,ij;
  REAL btmp;
  INT edge_flag=-10,issym=1;
  INT *ia=ain->IA, *ja=ain->JA ;
  REAL *a=ain->val;
  //  dvector diag0 = dvec_create(nv);
  //  for (i=0;i<nv;i++) {diag0.val[i]=0.;}
  REAL bte;
  REAL *aloc=calloc(dim1*dim1,sizeof(REAL));//random size but enough
  // over all elements
  INT nrbits=dim*sizeof(REAL),nibits=dim1*sizeof(INT);
  INT is=(INT )(sizeof(REAL)/sizeof(char));
  void *param=calloc(is*dim2,sizeof(char));//random size but enough
  void *wrk=calloc(dim1*dim1*15*is,sizeof(char));//random size but enough
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
    /* dummy  to make diff somewhat different */
    /* for(i = 0;i<dim;i++){ */
    /*   for(j = 0;j<dim;j++){ */
    /* 	btmp=0.; */
    /* 	for (k=0;k<dim;k++){ */
    /* 	  btmp+=bs[i*dim+k]*bs[j*dim+k]; */
    /* 	} */
    /* 	diff[i*dim+j]=btmp; */
    /*   } */
    /*   diff[i*dim+i]+=1.; */
    /* } */
    print_full_mat(dim,dim,diff,"diff");
    /*******************************************/
    print_full_mat(1,dim,ad,"ad");
    /* ad <- D^{-1}\beta */
    solve_pivot(1, dim, diff, ad, p, piv);
    print_full_mat(1,dim,ad,"ad1");
    /* // inv(B)*a*inv(B') */
    if(issym){
      /* SYMMETRIC */
      for(i = 0;i<dim1;i++){
	for(j = i;j<dim1;j++){
	  btmp=0.;
	  for (k=0;k<dim;k++){
	    for (l=0;l<dim;l++){
	      btmp+=bsinv[i*dim+k]*diff[k*dim+l]*bsinv[j*dim+l];
	    }
	  }
	  aloc[i*dim1+j]=btmp;
	}
      }
      for(i = 1;i<dim1;i++)
	for(j = 0;j<i;j++) aloc[i*dim1+j]=aloc[i+j*dim1];
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
	  aloc[i*dim1+j]=btmp;
	}
      }
    }
    /**********************************************************************/
    for(i=0;i<dim1;i++){
      id=i*dim; // this is for xs, so only dim columns
      for(j=0;j<dim1;j++){
	if(i==j) continue;
	/* \beta\cdot \tau_e */
	ij=i*dim1+j; // this is for stiff local.
	jd=j*dim; // this is for xs, so only dim columns
	bte=0.; 
	//	if(nop[j]>nop[i]){
	  for(k=0;k<dim;k++) bte+=ad[k]*(xs[id+k]-xs[jd+k]);
	  //	  fprintf(stdout,"\nFRW: (%d--%d)bte(%d,%d)=%g",nop[i],nop[j],i,j,bte);fflush(stdout);
	  //	} else if(nop[j]!=nop[i]) {// reverse tau
	  //	  for(k=0;k<dim;k++) bte+=ad[k]*(xs[id+k]-xs[jd+k]);
	  //	  fprintf(stdout,"\nREV: (%d--%d)bte(%d,%d)=%g",nop[i],nop[j],i,j,bte);fflush(stdout);
	  //      }
	//old for the above... for(k=0;k<dim;k++) te[k]=xs[jd+k]-xs[id+k];
	//old for the above... for(k=0;k<dim;k++) te[k]=xs[id+k]-xs[jd+k];
	/*midpoints*/
	// mid points on the edges if needed!    for(k=0;k<dim;k++) xm[k]=0.5*(xs[jd+k]+xs[id+k]);
	aloc[ij] *= bernoulli1(bte); 
	//	fprintf(stdout,"\ndd=%22.14e",(alpe*(bernoulli(bte/alpe))));
	/* the diagonal is equal to the negative column sum;*/    
	piv[j]-=aloc[ij];
      }
    }
    /*DIAGONAL= negative column sum */
    for(j=0;j<dim1;j++){
      btmp=0.;
      for(i=0;i<dim1;i++){
	if(i==j) continue;
	ij=i*dim1+j; // 
	btmp+=aloc[ij];
      }
      aloc[j*dim1+j] = -btmp;
    }
    //    for (i = 0; i < dim1; i++) aloc[i+dim1+i]=piv[i];
/***********************************************************************/
    print_full_mat(dim1,dim,bsinv,"bsinv");
    print_full_mat(dim1,dim1,aloc,"diffnew");
  }/* end of element loop*/
  fprintf(stdout,"\n%%%% end of element loop...\n");  
  //  dvec_free(&diag0);
  if(wrk)free(wrk);
  if(aloc)free(aloc);
  return;
}
/*C====================================================================*/
void symba(scomplex *sc,dCSRmat *ain,INT n,INT *nnz, INT *iwk){
 /*C====================================================================*/
  INT *ia=ain->IA,*ja=ain->JA;
  INT ndof=sc->n+1,i,j,k,l,ll,iii,iaa,iab,jk;
  REAL *a=ain->val;
  iCSRmat *ieje=NULL,*iejet=NULL;
  ieje->IA=calloc(sc->ns*ndof,sizeof(INT));
    ieje->JA=nodes;
  /*
C--------------------------------------------------------------------
C...  Symbolic assembly of a general nodal assembly matrix.
C...
C...  Input:
C...    IE, JE   - mesh connectivity matrix in RRCU.
C...    IET, JET - transpose of IE and JE in RRCO.
C...    N        - number of nodes in the mesh.
C...    IWK      - array of dimension N+1 which contains the value  
C...               N+1 in the positions corresponding to Dirichlet 
C...               nodes and 0 elsewhere.
C...
C...  Output:
C...    IA, JA   - structure of the nodal assembly matrix of 
C...               dimension NxN, in RRCU form.
C...
C...  Note:
C...    N+1 is the dimension of IA.
C...    NNZ is the dimension of JA.
C--------------------------------------------------------------------
*/
      np = n + 1
      nnz = 0
      jp = 1
      do 50 i = 1, n
         jpi = jp
         if(iwk(i) .eq. np) then
            ja(jp) = i
            jp = jp + 1
            go to 40
         end if
         ieta = iet(i)
         ietb = iet(i+1) - 1
         do 30 ip = ieta,ietb
            j = jet(ip)
            iea = ie(j)
            ieb = ie(j+1) - 1
            do 20 kp = iea,ieb
               k = je(kp)
               if (iwk(k) .ge. i) go to 20
               ja(jp) = k
               jp = jp + 1
               iwk(k) = i
 20         continue
 30      continue
 40      ia(i) = jpi
 50   continue
      ia(np) = jp
      nnz = ia(np)-1
      return
      end

