/*
 *  error.c
 *  
 *  Created by James Adler and Xiaozhe Hu on 5/13/15.
 *  Copyright 2015_JXLcode__. All rights reserved.
 *
 */

// Standard Includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// Our Includes
#include "macro.h"
#include "grid.h"
#include "sparse.h"
#include "vec.h"
#include "functs.h"
#include "fem.h"


/****************************************************************************************************************************/
// Error and Norm ROUTINES
/****************************************************************************************************************************/

/***************************************************************************/
void L2norm(REAL *norm,REAL *u,fespace *FE,trimesh *mesh,qcoordinates *cq)
{

  /* Computes the L2 Norm of a FE approximation using the mass matrix assembly any type of element.
   *
   * Input:		u	Numerical Solution at DOF
   *                    FE      fespace struct
   *                    mesh    mesh struct
   *                    cq      Quadrature Nodes
   *
   * Output:		norm	L2 Norm
   *
   */
		
  INT i,j,k,rowa,rowb,jcntr,elm;
  REAL sum = 0.0;

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* MLoc = calloc(local_size,sizeof(REAL));
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  /* Loop over all Elements */
  for (i=0; i<FE->nelm; i++) {

    // Zero out local matrices
    for (j=0; j<local_size; j++) MLoc[j] = 0.0;
    
    // Find Nodes for given Element
    rowa = FE->el_dof->IA[i]-1;
    rowb = FE->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm[jcntr] = FE->el_dof->JA[j];
      jcntr++;
    }

    //Find Nodes for given Element if not H1 elements
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr=0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }
    		
    // Compute Local Stiffness Matrix for given Element
    elm = i+1;
    assemble_mass_local(MLoc,FE,mesh,cq,dof_on_elm,v_on_elm,elm,constcoeff,1.0);

    for(j=0;j<dof_per_elm;j++) {
      for(k=0;k<dof_per_elm;k++) {
	sum+=u[dof_on_elm[j]-1]*MLoc[j*dof_per_elm+k]*u[dof_on_elm[k]-1];
      }
    }
  }

  *norm = sqrt(sum);

  if(MLoc) free(MLoc);
  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);

  return;
}
/*******************************************************************************************************************************************************/

/***************************************************************************/
void L2_InnerProduct(REAL *product,REAL *u,REAL *v,fespace *FE,trimesh *mesh,qcoordinates *cq)
{

  /* Computes the L2 inner product of two FE approximations using the mass matrix assembly any type of element.
   *
   * Input:		u,v	Numerical Solutions at DOF
   *                    FE      fespace struct
   *                    mesh    mesh struct
   *                    cq      Quadrature Nodes
   *
   * Output:		product	L2 Inner Product of u and v, <u,v>
   *
   */
		
  INT i,j,k,rowa,rowb,jcntr,elm;
  REAL sum = 0.0;

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* MLoc = calloc(local_size,sizeof(REAL));
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  /* Loop over all Elements */
  for (i=0; i<FE->nelm; i++) {

    // Zero out local matrices
    for (j=0; j<local_size; j++) MLoc[j] = 0.0;
    
    // Find Nodes for given Element
    rowa = FE->el_dof->IA[i]-1;
    rowb = FE->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm[jcntr] = FE->el_dof->JA[j];
      jcntr++;
    }

    //Find Nodes for given Element if not H1 elements
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr=0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }
    		
    // Compute Local Stiffness Matrix for given Element
    elm = i+1;
    assemble_mass_local(MLoc,FE,mesh,cq,dof_on_elm,v_on_elm,elm,constcoeff,1.0);

    for(j=0;j<dof_per_elm;j++) {
      for(k=0;k<dof_per_elm;k++) {
	sum+=v[dof_on_elm[j]-1]*MLoc[j*dof_per_elm+k]*u[dof_on_elm[k]-1];
      }
    }
  }

  *product= sqrt(sum);

  if(MLoc) free(MLoc);
  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);

  return;
}
/*******************************************************************************************************************************************************/

/* /\***************************************************************************\/ */
/* void L2error(REAL *err,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT test, \ */
/* 	     REAL param,REAL time)  */
/* { */
	
/*   /\* Computes the L2 Norm of the Error using a high order quadrature or the Mass matrix if given */
/*    * */
/*    * Input:		u				Numerical Solution at DOF */
/*    *	       		xn,yn,zn		Coordinates of Nodes */
/*    *	       		iel_n,jel_n		Element to Node Map */
/*    *	       		nelm			Number of Elements */
/*    *	       		nve			Number of Vertices per element */
/*    *   			mydim			Dimension */
/*    *   			nq1d			Number of Quadrature Points per dimension */
/*    *	       		element_order	        Number of Nodes per Element */
/*    *		       	test			Which solution to compare to */
/*    * */
/*    * Output:		err				L2 Norm of the Error */
/*    * */
/*    *\/ */
	
/*   INT i,j,rowa,iq;				/\* Loop Indices *\/ */
	
/*   REAL el2 = 0;			/\* Error *\/ */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL uh,utrue;			/\* Numerical Solution and True Solution at Quadrature Nodes *\/ */
	
/*   for (i=0; i<nelm; i++) { */
		
/*     // Get Vertices (must be vertices) */
/*     rowa = iel_n[i]-1; */
/*     for (j=0; j<nve; j++) {  // First ones are always vertices  */
/*       ipv[j] = jel_n[rowa+j]; */
/*       ipn[j] = ipv[j]; */
/*     } */
/*     for (j=nve; j<element_order; j++) { */
/*       ipn[j] = jel_n[rowa+j]; */
/*     } */
/*     if (mydim==2) { */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
/*     } else if (mydim==3) { */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
/*     } else { */
/*       printf("Bad Dimension"); */
/*       return; */
/*     } */
		
/*     for (iq=0;iq<nq;iq++) { */
/*       x = xq[iq]; */
/*       y = yq[iq]; */
/*       if (mydim>2) { z = zq[iq]; } */
/*       w = wq[iq]; */
			
/*       if (element_order==1) { */
/* 	uh=u[i]; */
/*       } else { */
/* 	H1_interpolation(&uh,u,x,y,z,ipn,xn,yn,zn,0,element_order,mydim,1); */
/*       } */

/*       // get true solution at quadrature node */
/*       getknownfunction(&utrue,x,y,z,time,mydim,1,param,test); */
			
/*       el2 = el2 + ( uh - utrue )*( uh - utrue)*w; */
/*     } */
		
/*   }	 */
	
/*   el2 = sqrt ( el2 ); */
 
	
/*   *err = el2; */
	
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(ipn) free(ipn); */
/*   if(ipv) free(ipv); */
	
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\***************************************************************************\/ */
/* void H1seminormLapassemble(REAL *norm,REAL *u,CSRinc el_dof,CSRinc el_n,coordinates cn,INT nelm,INT nq1d,INT mydim,INT dof_order,INT element_order, \ */
/* 			   CSRinc dof_n,REAL* f_area,REAL* f_norm,INT elementtype)  */
/* { */

/*   /\* Computes the H1 semi-Norm of a vector using the laplacian-like matrix if given for any type of element... */
/*    *          Nodal - <grad u, grad u> - |u|_1 */
/*    *          RT - <div u, div u>      - |u|_(H(div)) */
/*    *          Nedelec - <curl u, curl u>  - |u|_(H(curl)) */
/*    * */
/*    * Input:		u		       	Numerical Solution at DOF */
/*    *                    ndof                    Number of DOF */
/*    * */
/*    * Output:		norm		       	L2 Norm */
/*    * */
/*    *\/ */
		
/*   INT i,j,k,rowa,rowb,ndof,jcntr,elm; */
/*   REAL sum = 0.0; */
/*   REAL* MLoc = calloc(dof_order*dof_order,sizeof(REAL)); */
/*   INT* ipdof = calloc(dof_order,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
/*   REAL nq; */
/*   qcoordinates cq; */

/*   /\* Get Quadrature Nodes *\/ */
/*   nq = pow(nq1d,mydim); */
/*   allocateqcoords(&cq,nq1d,nelm,mydim); */
/*   if (mydim==2) { */
/*     get_quadrature(cq.x,cq.y,NULL,cq.w,cn.x,cn.y,NULL,el_n.IA,el_n.JA,nelm,element_order,nq1d,mydim); */
/*   } else { */
/*     get_quadrature(cq.x,cq.y,cq.z,cq.w,cn.x,cn.y,cn.z,el_n.IA,el_n.JA,nelm,element_order,nq1d,mydim); */
/*   } */

/*   /\* Loop over all Elements *\/ */
/*   for (i=0; i<nelm; i++) { */
		
/*     ndof = el_dof.IA[i+1]-el_dof.IA[i]; */
		
/*     // Zero out local matrices */
/*     for (j=0; j<ndof*ndof; j++) { */
/*       MLoc[j]=0.0; */
/*     } */
		
/*     // Find DOF for given Element if not H1 elements */
/*     rowa = el_dof.IA[i]-1; */
/*     rowb = el_dof.IA[i+1]-1; */
/*     jcntr = 0; */
/*     for (j=rowa; j<rowb; j++) { */
/*       ipdof[jcntr] = el_dof.JA[j]; */
/*       jcntr++; */
/*     } */

/*     //Find Nodes for given Element if not H1 elements */
/*     rowa = el_n.IA[i]-1; */
/*     rowb = el_n.IA[i+1]-1; */
/*     jcntr=0; */
/*     for (j=rowa; j<rowb; j++) { */
/*       ipn[jcntr] = el_n.JA[j]; */
/*       jcntr++; */
/*     } */
    
		
/*     // Compute Local Stiffness Matrix for given Element */
/*     elm = i+1; */

/*     if(elementtype==0) { // H1 Elements Linears or Quadratics */
/*       H1_LapL(MLoc,cn.x,cn.y,cn.z,cq.x,cq.y,cq.z,cq.w,ipn,ndof,nq1d,mydim,elm); */
/*     } else if(elementtype==1) { // Nedelec Elements */
/*       if (mydim==2) { */
/* 	ned_curlcurlL2D(MLoc,cn.x,cn.y,cq.x,cq.y,cq.w,ipn,ipdof,dof_n.IA,dof_n.JA,ndof,element_order,nq1d,mydim,elm); */
/*       } else if (mydim==3) { */
/* 	ned_curlcurlL3D(MLoc,cn.x,cn.y,cn.z,cq.x,cq.y,cq.z,cq.w,ipn,ipdof,dof_n.IA,dof_n.JA,ndof,element_order,nq1d,mydim,elm); */
/*       } else { */
/* 	printf("Your dimension isn't 2 or 3.  Welcome to the Twilight Zone..."); */
/*       } */
/*     } else if(elementtype==2) { // Raviart-Thomas Elements */
/*       rt_divdivL(MLoc,cn.x,cn.y,cn.z,cq.x,cq.y,cq.z,cq.w,ipn,ipdof,dof_n.IA,dof_n.JA,ndof,f_area,f_norm,element_order,nq1d,mydim,elm); */
/*     } */

/*     for(j=0;j<ndof;j++) { */
/*       for(k=0;k<ndof;k++) { */
/* 	sum+=u[ipdof[j]-1]*MLoc[j*ndof+k]*u[ipdof[k]-1]; */
/*       } */
/*     } */
/*   } */

/*   if(sum<0.0) { */
/*     printf("Your H1 Semi Norm Squared is negative!  Outputting the square itself\n"); */
/*     *norm = sum; */
/*   } else { */
/*     *norm = sqrt(sum); */
/*   } */

/*   if(MLoc) free(MLoc); */
/*   if(ipn) free(ipn); */
/*   if(ipdof) free(ipdof); */
/*   freeqcoords(cq); */

/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\***************************************************************************\/ */
/* void H1semierror(REAL *err,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT test,REAL param,REAL time)  */
/* { */

/*   /\* Computes the H1 Semi-Norm of the Error using a high order quadrature  */
/*    * */
/*    * Input:		u				Numerical Solution at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				nelm			Number of Elements */
/*    *				ndof			Number of DOF */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature Points per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    *				test			Which solution to compare to (This is the gradient of the solution so it's a vector) */
/*    * */
/*    * Output:		err				H1 Norm of the Error */
/*    * */
/*    *\/ */
	
/*   INT i,j,rowa,iq;				/\* Loop Indices *\/ */

/*   REAL el = 0.0;			 */
/*   REAL elx = 0.0;			 */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0; */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* ugh = calloc(mydim,sizeof(REAL));			/\* Grad of approximation interpolated at Quadrature Nodes *\/ */
/*   REAL* ugtrue = calloc(mydim,sizeof(REAL));		/\* Grad of true solution interpolated at Quadrature Nodes *\/ */
	
/*   if (mydim==2) { */
/*     for (i=0; i<nelm; i++) { */
		
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d);		 */
		
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
			
/* 	H1Grad_interpolation(ugh,u,x,y,z,ipn,xn,yn,zn,0,element_order,mydim,1); */
			
/* 	// get true solution at quadrature node */
/* 	getknownfunction(ugtrue,x,y,z,time,mydim,2,param,test); */
			
/* 	elx = elx + ( ugh[0] - ugtrue[0] )*( ugh[0] - ugtrue[0] )*w; */
/* 	ely = ely + ( ugh[1] - ugtrue[1] )*( ugh[1] - ugtrue[1] )*w; */
/*       } */
		
/*     } */
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */

/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
				
/* 	H1Grad_interpolation(ugh,u,x,y,z,ipn,xn,yn,zn,0,element_order,mydim,1); */
				
/* 	// get true solution at quadrature node */
/* 	getknownfunction(ugtrue,x,y,z,time,mydim,2,param,test); */
				
/* 	elx = elx + ( ugh[0] - ugtrue[0] )*( ugh[0] - ugtrue[0] )*w; */
/* 	ely = ely + ( ugh[1] - ugtrue[1] )*( ugh[1] - ugtrue[1] )*w; */
/* 	elz = elz + ( ugh[2] - ugtrue[2] )*( ugh[2] - ugtrue[2] )*w; */
/*       } */
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3.  YOu are in the twilight zone"); */
/*     exit(0); */
/*   } */

	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *err = el; */
	
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(ipn) free(ipn); */
/*   if(ipv) free(ipv); */
	
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\*******************************************************************************************************************************************************\/ */
/* void L2errorRT(REAL *err,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT *iel_f,INT *jel_f,INT *if_n,INT *jf_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT face_order,INT test,REAL param,REAL time,REAL* f_area,REAL* f_norm) */
/* { */
	
/*   /\* Computes the L2 error of an approximation defined on Nedelec elements using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				iel_ed,jel_ed	Element to Edge Map */
/*    *				ied_n,jed_n		Edge to Node Map */
/*    *				nelm			Number of Elements */
/*    *				nve				Number of vertices per element */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature PoINTs per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    *				edge_order		Number of Edges per Element */
/*    * */
/*    * Output:		norm			L2 Norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
	
/*   REAL el = 0.0; */
/*   REAL elx = 0.0;			/\* Error *\/ */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0;	 */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
/*   INT* ipf = calloc(face_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* uh = calloc(mydim,sizeof(REAL));			/\* Numerical Solution interpolated at Quadrature Nodes *\/ */
/*   REAL* utrue = calloc(mydim,sizeof(REAL));		/\* True Solution at Quadrature Nodes *\/ */
	
/*   if (mydim==2) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_f[i]-1; */
/*       for (j=0; j<face_order; j++) { */
/* 	ipf[j] = jel_f[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
				
/* 	RT_interpolation(uh,u,x,y,z,ipn,ipf,if_n,jf_n,xn,yn,zn,element_order,face_order,mydim,f_area,f_norm); */
				
/* 	// get true solution at quadrature node */
/* 	getknownfunction(utrue,x,y,z,time,mydim,2,param,test); */
/* 	elx = elx + (uh[0]-utrue[0])*(uh[0]-utrue[0])*w; */
/* 	ely = ely + (uh[1]-utrue[1])*(uh[1]-utrue[1])*w; */
/*       } */
			
/*     }	 */
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_f[i]-1; */
/*       for (j=0; j<face_order; j++) { */
/* 	ipf[j] = jel_f[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
				
/* 	RT_interpolation(uh,u,x,y,z,ipn,ipf,if_n,jf_n,xn,yn,zn,element_order,face_order,mydim,f_area,f_norm); */
				
/* 	// get true solution at quadrature node */
/* 	getknownfunction(utrue,x,y,z,time,mydim,2,param,test); */
				
/* 	elx = elx + (uh[0]-utrue[0])*(uh[0]-utrue[0])*w; */
/* 	ely = ely + (uh[1]-utrue[1])*(uh[1]-utrue[1])*w; */
/* 	elz = elz + (uh[2]-utrue[2])*(uh[2]-utrue[2])*w; */
/*       } */
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3\n"); */
/*     exit(0); */
/*   } */
	
	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *err = el; */
	
/*   if(ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(ipf) free(ipf); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\*******************************************************************************************************************************************************\/ */
/* void L2errorNed(REAL *err,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT *iel_ed,INT *jel_ed,INT *ied_n,INT *jed_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT edge_order,INT test,REAL param,REAL time)   */
/* { */
	
/*   /\* Computes the L2 error of an approximation defined on Nedelec elements using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				iel_ed,jel_ed	Element to Edge Map */
/*    *				ied_n,jed_n		Edge to Node Map */
/*    *				nelm			Number of Elements */
/*    *				nve				Number of vertices per element */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature Points per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    *				edge_order		Number of Edges per Element */
/*    * */
/*    * Output:		norm			L2 Norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
/*   REAL el = 0.0; */
/*   REAL elx = 0.0;			/\* Error *\/ */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0;	 */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
/*   INT* ipe = calloc(edge_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* uh = calloc(mydim,sizeof(REAL));			/\* Numerical Solution interpolated at Quadrature Nodes *\/ */
/*   REAL* utrue = calloc(mydim,sizeof(REAL));		/\* True Solution at Quadrature Nodes *\/ */
	
/*   if (mydim==2) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
				
/* 	Ned_interpolation(uh,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
				
/* 	// get true solution at quadrature node */
/* 	getknownfunction(utrue,x,y,z,time,mydim,2,param,test); */

/* 	elx = elx + (uh[0]-utrue[0])*(uh[0]-utrue[0])*w; */
/* 	ely = ely + (uh[1]-utrue[1])*(uh[1]-utrue[1])*w; */
/*       } */
			
/*     }	 */
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
/* 	Ned_interpolation(uh,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
				
/* 	// get true solution at quadrature node */
/* 	getknownfunction(utrue,x,y,z,time,mydim,2,param,test); */
				
/* 	elx = elx + (uh[0]-utrue[0])*(uh[0]-utrue[0])*w; */
/* 	ely = ely + (uh[1]-utrue[1])*(uh[1]-utrue[1])*w; */
/* 	elz = elz + (uh[2]-utrue[2])*(uh[2]-utrue[2])*w; */
/*       } */
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3\n"); */
/*     exit(0); */
/*   } */
	
	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *err = el; */
	
/*   if (ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(ipe) free(ipe); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */


/* /\*******************************************************************************************************************************************************\/ */
/* void NedsemiCurlerror(REAL *err,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT *iel_ed,INT *jel_ed,INT *ied_n,INT *jed_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT edge_order,INT test,REAL param,REAL time)  */
/* { */
	
/*   /\* Computes the H curl Semi-Norm Error of an approximation using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				iel_ed,jel_ed	Element to Edge Map */
/*    *				ied_n,jed_n		Edge to Node Map */
/*    *				nelm			Number of Elements */
/*    *				nve				Number of vertices per element */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature PoINTs per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    *				edge_order		Number of Edges per Element */
/*    * */
/*    * Output:		norm			L2 Norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
	
/*   REAL el = 0.0; */
/*   REAL elx = 0.0;			/\* Error *\/ */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0; */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
/*   INT* ipe = calloc(edge_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* uch = calloc(mydim,sizeof(REAL));			/\* Curl of approximation interpolated at Quadrature Nodes *\/ */
/*   REAL* uctrue = calloc(mydim,sizeof(REAL));		/\* Curl of true solution interpolated at Quadrature Nodes *\/ */
	
/*   if (mydim==2 ) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
				
/* 	NedCurl_interpolation(uch,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
/* 	getknownfunction(uctrue,x,y,z,time,mydim,1,param,test); */
				
/* 	elx = elx + (uch[0]-uctrue[0])*(uch[0]-uctrue[0])*w; */
/*       } */
/*     } */
		
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
				
/* 	NedCurl_interpolation(uch,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
				
/* 	// True Solution */
/* 	getknownfunction(uctrue,x,y,z,time,mydim,2,param,test); */

				
/* 	elx = elx + (uch[0]-uctrue[0])*(uch[0]-uctrue[0])*w; */
/* 	ely = ely + (uch[1]-uctrue[1])*(uch[1]-uctrue[1])*w; */
/* 	elz = elz + (uch[2]-uctrue[2])*(uch[2]-uctrue[2])*w; */
/*       } */
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3.  What is this the Twighlight Zone??"); */
/*     exit(0); */
/*   } */
	
	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *err = el; */
	
/*   if (ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(ipe) free(ipe); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\*******************************************************************************************************************************************************\/ */
/* void L2norminterp(REAL *norm,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order)  */
/* { */
	
/*   /\* Computes the L2 Norm of an approximation using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				nelm			Number of Elements */
/*    *				ndof			Number of DOF */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature Points per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    * */
/*    * Output:		norm			L2 Norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
	
/*   REAL el2 = 0.0; */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL uh;			/\* Numerical Solution interpolated at Quadrature Nodes *\/ */
	
/*   for (i=0; i<nelm; i++) { */
		
/*     // Get Vertices (must be vertices) */
/*     rowa = iel_n[i]-1; */
/*     for (j=0; j<nve; j++) {  // First ones are always vertices  */
/*       ipv[j] = jel_n[rowa+j]; */
/*       ipn[j] = ipv[j]; */
/*     } */
/*     for (j=nve; j<element_order; j++) { */
/*       ipn[j] = jel_n[rowa+j]; */
/*     } */
		
/*     // Get quadrature nodes for given order of accuracy */
/*     if (mydim==2) { */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
/*     } else if (mydim==3) { */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
/*     } else { */
/*       printf("Bad Dimension"); */
/*       return; */
/*     } */
		
/*     for (iq=0;iq<nq;iq++) { */
/*       x = xq[iq]; */
/*       y = yq[iq]; */
/*       if (mydim>2) { z = zq[iq]; } */
/*       w = wq[iq]; */
			
/*       if (element_order==1) { */
/* 	uh=u[i]; */
/*       } else { */
/* 	H1_interpolation(&uh,u,x,y,z,ipn,xn,yn,zn,0,element_order,mydim,1); */
/*       } */
			
/*       el2 = el2 + uh*uh*w; */
/*     } */
		
/*   }	 */
	
/*   el2 = sqrt ( el2 ); */
	
/*   *norm = el2; */
	
/*   if (ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\*******************************************************************************************************************************************************\/ */
/* void H1seminorminterp(REAL *norm,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order)  */
/* { */
	
/*   /\* Computes the H1 Semi-Norm of an approximation using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				nelm			Number of Elements */
/*    *				ndof			Number of DOF */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature Points per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    * */
/*    * Output:		norm				H1 semi-norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
	
/*   REAL el = 0.0; */
/*   REAL elx = 0.0;			/\* Error *\/ */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0; */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* uh = calloc(mydim,sizeof(REAL));			/\* Grad of approximation interpolated at Quadrature Nodes *\/ */
	
/*   if (mydim==2 ) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
				
/* 	H1Grad_interpolation(uh,u,x,y,z,ipn,xn,yn,zn,0,element_order,mydim,1); */
				
/* 	elx = elx + uh[0]*uh[0]*w; */
/* 	ely = ely + uh[1]*uh[1]*w; */
/*       } */
/*     } */
		
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
				
/* 	H1Grad_interpolation(uh,u,x,y,z,ipn,xn,yn,zn,0,element_order,mydim,1); */
				
/* 	elx = elx + uh[0]*uh[0]*w; */
/* 	ely = ely + uh[1]*uh[1]*w; */
/* 	elz = elz + uh[2]*uh[2]*w; */
/*       } */
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3.  What is this the Twighlight Zone??"); */
/*     exit(0); */
/*   } */
	
	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *norm = el; */
	
/*   if (ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\*******************************************************************************************************************************************************\/ */
/* void L2normNed(REAL *norm,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT *iel_ed,INT *jel_ed,INT *ied_n,INT *jed_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT edge_order)  */
/* { */
	
/*   /\* Computes the L2 Norm of an approximation defined on Nedelec elements using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				iel_ed,jel_ed	Element to Edge Map */
/*    *				ied_n,jed_n		Edge to Node Map */
/*    *				nelm			Number of Elements */
/*    *				nve				Number of vertices per element */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature Points per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    *				edge_order		Number of Edges per Element */
/*    * */
/*    * Output:		norm			L2 Norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
	
/*   REAL el = 0.0; */
/*   REAL elx = 0.0;			/\* Error *\/ */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0;	 */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
/*   INT* ipe = calloc(edge_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* uh = calloc(mydim,sizeof(REAL));			/\* Numerical Solution interpolated at Quadrature Nodes *\/ */
	
/*   if (mydim==2) { */
/*     for (i=0; i<nelm; i++) { */
		
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
		
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
		
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
			
/* 	Ned_interpolation(uh,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
			
/* 	elx = elx + uh[0]*uh[0]*w; */
/* 	ely = ely + uh[1]*uh[1]*w; */
/*       } */
		
/*     }	 */
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
				
/* 	Ned_interpolation(uh,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
				
/* 	elx = elx + uh[0]*uh[0]*w; */
/* 	ely = ely + uh[1]*uh[1]*w; */
/* 	elz = elz + uh[2]*uh[2]*w; */
/*       } */
			
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3\n"); */
/*     exit(0); */
/*   } */

	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *norm = el; */
	
/*   if (ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(ipe) free(ipe); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\*******************************************************************************************************************************************************\/ */
/* void L2normRT(REAL *norm,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT *iel_f,INT *jel_f,INT *if_n,INT *jf_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT face_order,REAL* f_area,REAL* f_norm)  */
/* { */
	
/*   /\* Computes the L2 Norm of an approximation defined on Nedelec elements using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				iel_ed,jel_ed	Element to Edge Map */
/*    *				ied_n,jed_n		Edge to Node Map */
/*    *				nelm			Number of Elements */
/*    *				nve				Number of vertices per element */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature Points per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    *				edge_order		Number of Edges per Element */
/*    * */
/*    * Output:		norm			L2 Norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
	
/*   REAL el = 0.0; */
/*   REAL elx = 0.0;			/\* Error *\/ */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0;	 */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
/*   INT* ipf = calloc(face_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* uh = calloc(mydim,sizeof(REAL));			/\* Numerical Solution interpolated at Quadrature Nodes *\/ */
	
/*   if (mydim==2) { */
/*     for (i=0; i<nelm; i++) { */
		
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_f[i]-1; */
/*       for (j=0; j<face_order; j++) { */
/* 	ipf[j] = jel_f[rowa+j]; */
/*       } */
		
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
		
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
			
/* 	RT_interpolation(uh,u,x,y,z,ipn,ipf,if_n,jf_n,xn,yn,zn,element_order,face_order,mydim,f_area,f_norm); */
			
/* 	elx = elx + uh[0]*uh[0]*w; */
/* 	ely = ely + uh[1]*uh[1]*w; */
/*       } */
		
/*     }	 */
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_f[i]-1; */
/*       for (j=0; j<face_order; j++) { */
/* 	ipf[j] = jel_f[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
				
/* 	RT_interpolation(uh,u,x,y,z,ipn,ipf,if_n,jf_n,xn,yn,zn,element_order,face_order,mydim,f_area,f_norm); */
				
/* 	elx = elx + uh[0]*uh[0]*w; */
/* 	ely = ely + uh[1]*uh[1]*w; */
/* 	elz = elz + uh[2]*uh[2]*w; */
/*       } */
			
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3\n"); */
/*     exit(0); */
/*   } */

	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *norm = el; */
	
/*   if (ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(ipf) free(ipf); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\*******************************************************************************************************************************************************\/ */
/* void NedsemiCurlnorm(REAL *norm,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT *iel_ed,INT *jel_ed,INT *ied_n,INT *jed_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT edge_order)  */
/* { */
	
/*   /\* Computes the H curl Semi-Norm of an approximation using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				iel_ed,jel_ed	Element to Edge Map */
/*    *				ied_n,jed_n		Edge to Node Map */
/*    *				nelm			Number of Elements */
/*    *				nve				Number of vertices per element */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature Points per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    *				edge_order		Number of Edges per Element */
/*    * */
/*    * Output:		norm			L2 Norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
	
/*   REAL el = 0.0; */
/*   REAL elx = 0.0;			/\* Error *\/ */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0; */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
/*   INT* ipe = calloc(edge_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* uh = calloc(mydim,sizeof(REAL));			/\* Curl of approximation interpolated at Quadrature Nodes *\/ */
	
/*   if (mydim==2 ) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
				
/* 	NedCurl_interpolation(uh,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
				
/* 	elx = elx + uh[0]*uh[0]*w; */
/*       } */
/*     } */
		
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
				
/* 	NedCurl_interpolation(uh,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
				
/* 	elx = elx + uh[0]*uh[0]*w; */
/* 	ely = ely + uh[1]*uh[1]*w; */
/* 	elz = elz + uh[2]*uh[2]*w; */
/*       } */
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3.  What is this the Twighlight Zone??"); */
/*     exit(0); */
/*   } */
	
	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *norm = el; */
	
/*   if (ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(ipe) free(ipe); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */
