/*
 *  assemble_global.c   
 *  
 *  Created by James Adler and Xiaozhe Hu on 4/22/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

/* This code will build global stiffness matrices for various PDE systems */

#include "hazmat.h"

/******************************************************************************************************/
void assemble_global_withBC(dCSRmat* A,dvector *b,void (*local_assembly)(REAL *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),void (*bc)(REAL *,REAL *,REAL),void (*coeff)(REAL *,REAL *,REAL),REAL time) 
{
	
  /* Computes the global stiffness matrix and rhs for any a(u,v) = <f,v> bilinear form using various element types
   * (eg. P1, P2, Nedelec, and Raviart-Thomas).
   * Also takes care of Dirichlet boundary conditions.  If the node is a boundary the row will be
   * zeroed out except for the diagonal entry being 1.  The corresponding column will also be 0 and 
   * the right-hand side adjusted.
   *
   * For this problem we compute:
   *
   * Lu = f  ---->   a(u,v) = <f,v>
   *
   * which gives Ax = b,
   *
   * A_ij = a( phi_j, phi_i)
   * b_i  = <f,phi_i>
   * 
   *    INPUT:
   *            local_assembly      Function to local assembly routine
   *            fespace		    Finite-Element Space Struct
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights
   *            rhs                 Function that gives right-hand side at given coordinates rhs(&val,x,time)
   *            bc                  Function that gives boundary condition at given coordinates bc(&val,x,time)
   *            coeff               Function that gives coefficient (for now assume constant)
   *            time                Time if time dependent
   *
   *    OUTPUT:
   *	       	A                   Global Stiffness Matrix (CSR Format)
   *            b	       	    Global Right hand side vector
   */

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT i,j;

  // Allocate Row Array
  A->row = FE->ndof;
  A->col = FE->ndof;
  A->IA = (INT *) calloc(FE->ndof+1,sizeof(INT));
  if(rhs!=NULL) {
    b->row = FE->ndof;
    b->val = (REAL *) calloc(b->row,sizeof(REAL));
  }

  // Get Sparsity Structure First
  // Non-zeros of A and IA (ignores cancellations, so maybe more than necessary)
  stiffG_nnzBC(A,FE);
	
  // Columns of A -> JA
  A->JA = (INT *) calloc(A->nnz,sizeof(INT));
  stiffG_colsBC(A,FE);
  
  // Set values
  A->val = (REAL *) calloc(A->nnz,sizeof(REAL));
  for (i=0; i<A->nnz; i++) {
    A->val[i] = 0;
  }
  	
  // Now Build Global Matrix entries
  // First deal with boundary rows
  for (i=0; i<FE->ndof; i++) {
    if (FE->dof_bdry[i]==1) { /* This is a boundary row.  Just make identity and fix right hand side */
      A->val[A->IA[i]-1] = 1.0;
      if(rhs!=NULL)
	b->val[i] = FE_Evaluate_DOF(bc,FE,mesh,time,i);
    }
  }

  // Now adjust other rows
  /* Loop over all Elements and build local matrix and rhs */
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* ALoc = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* bLoc=NULL;
  if(rhs!=NULL)
    bLoc = (REAL *) calloc(dof_per_elm,sizeof(REAL));

  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT rowa,rowb,jcntr;
  for (i=0; i<FE->nelm; i++) {	
    // Zero out local matrices
    for (j=0; j<local_size; j++) {
      ALoc[j]=0;
    }
    if(rhs!=NULL) {
      for (j=0; j<dof_per_elm; j++) {
	bLoc[j]=0;
      }
    }
		
    // Find DOF for given Element
    rowa = FE->el_dof->IA[i]-1;
    rowb = FE->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm[jcntr] = FE->el_dof->JA[j];
      jcntr++;
    }

    // Find vertices for given Element
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }
		
    // Compute Local Stiffness Matrix for given Element
    (*local_assembly)(ALoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,coeff,time);
    if(rhs!=NULL)
      FEM_RHS_Local(bLoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,rhs,time);

    // Loop over DOF and place in appropriate slot globally
    LocaltoGlobalBC(dof_on_elm,FE,b,A,ALoc,bLoc);
  }	

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(ALoc) free(ALoc);
  if(bLoc) free(bLoc);
		
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void assemble_global(dCSRmat* A,dvector *b,void (*local_assembly)(REAL *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),void (*coeff)(REAL *,REAL *,REAL),REAL time) 
{
	
  /* Computes the global stiffness matrix and rhs for any a(u,v) = <f,v> bilinear form using various element types
   * (eg. P1, P2, Nedelec, and Raviart-Thomas).
   * DOES NOT take care of Dirichlet boundary conditions.  A separate routine will eliminate them later
   * This allows for several matrices to be assembled then added or concatenated together.
   *
   * For this problem we compute:
   *
   * Lu = f  ---->   a(u,v) = <f,v>
   *
   * which gives Ax = b,
   *
   * A_ij = a( phi_j, phi_i)
   * b_i  = <f,phi_i>
   * 
   *    INPUT:
   *            local_assembly      Function to local assembly routine
   *            fespace		    Finite-Element Space Struct
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights
   *            rhs                 Function that gives right-hand side at given coordinates rhs(&val,x,time)
   *            coeff               Function that gives coefficient (for now assume constant)
   *            time                Time if time dependent
   *
   *    OUTPUT:
   *	       	A                   Global Stiffness Matrix (CSR Format)
   *            b	       	    Global Right hand side vector
   */

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT i,j;

  // Allocate Row Array
  A->row = FE->ndof;
  A->col = FE->ndof;
  A->IA = (INT *) calloc(FE->ndof+1,sizeof(INT));
  if(rhs!=NULL) {
    b->row = FE->ndof;
    b->val = (REAL *) calloc(b->row,sizeof(REAL));
  }

  // Get Sparsity Structure First
  // Non-zeros of A and IA (ignores cancellations, so maybe more than necessary)
  stiffG_nnz(A,FE);
	
  // Columns of A -> JA
  A->JA = (INT *) calloc(A->nnz,sizeof(INT));
  stiffG_cols(A,FE);
  
  // Set values
  A->val = (REAL *) calloc(A->nnz,sizeof(REAL));
  for (i=0; i<A->nnz; i++) {
    A->val[i] = 0;
  }
  	
  // Now Build Global Matrix entries
 
  /* Loop over all Elements and build local matrix and rhs */
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* ALoc = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* bLoc=NULL;
  if(rhs!=NULL) 
    bLoc = (REAL *) calloc(dof_per_elm,sizeof(REAL));

  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT rowa,rowb,jcntr;
  for (i=0; i<FE->nelm; i++) {	
    // Zero out local matrices
    for (j=0; j<local_size; j++) {
      ALoc[j]=0;
    }
    if(rhs!=NULL) {
      for (j=0; j<dof_per_elm; j++) {
	bLoc[j]=0;
      }
    }
		
    // Find DOF for given Element
    rowa = FE->el_dof->IA[i]-1;
    rowb = FE->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm[jcntr] = FE->el_dof->JA[j];
      jcntr++;
    }

    // Find vertices for given Element
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }
		
    // Compute Local Stiffness Matrix for given Element
    (*local_assembly)(ALoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,coeff,time);
    if(rhs!=NULL)
      FEM_RHS_Local(bLoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,rhs,time);

    // Loop over DOF and place in appropriate slot globally
    LocaltoGlobal(dof_on_elm,FE,b,A,ALoc,bLoc);
  }	

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(ALoc) free(ALoc);
  if(bLoc) free(bLoc);
		
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void assemble_global_RHS(dvector *b,fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),REAL time) 
{
	
  /* Computes the global rhs for any a(u,v) = <f,v> bilinear form using various element types
   * (eg. P1, P2, Nedelec, and Raviart-Thomas).
   * DOES NOT take care of Dirichlet boundary conditions.  A separate routine will eliminate them later
   *
   * For this problem we compute RHS of:
   *
   * Lu = f  ---->   a(u,v) = <f,v>
   *
   * which gives Ax = b,
   *
   * A_ij = a( phi_j, phi_i)
   * b_i  = <f,phi_i>
   * 
   *    INPUT:
   *            fespace		    Finite-Element Space Struct
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights
   *            rhs                 Function that gives right-hand side at given coordinates rhs(&val,x,time)
   *            time                Time if time dependent
   *
   *    OUTPUT:
   *            b	       	    Global Right hand side vector
   */

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT i,j,row;

  // Allocate Arrays
  b->row = FE->ndof;
  b->val = (REAL *) calloc(b->row,sizeof(REAL));
 
  /* Loop over all Elements and build local rhs */
  REAL* bLoc= (REAL *) calloc(dof_per_elm,sizeof(REAL));

  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT rowa,rowb,jcntr;
  for (i=0; i<FE->nelm; i++) {	
    
    for (j=0; j<dof_per_elm; j++) {
      bLoc[j]=0;
    }
		
    // Find DOF for given Element
    rowa = FE->el_dof->IA[i]-1;
    rowb = FE->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm[jcntr] = FE->el_dof->JA[j];
      jcntr++;
    }

    // Find vertices for given Element
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }
		
    // Compute Local RHS for given Element
    FEM_RHS_Local(bLoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,rhs,time);

    // Loop over DOF and place in appropriate slot globally
    for (j=0; j<dof_per_elm; j++) { 
      row = dof_on_elm[j]-1;	
      b->val[row] = b->val[row] + bLoc[j];
    }
  }	

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(bLoc) free(bLoc);
		
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void assemble_global_face(dCSRmat* A,void (*local_assembly_face)(REAL *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT *,INT,INT,void (*)(REAL *,REAL *,REAL),REAL),fespace *FE,trimesh *mesh,qcoordinates *cq,void (*coeff)(REAL *,REAL *,REAL),REAL time,INT flag) 
{
	
  /* Computes the global stiffness matrix for any "boundary" bilinear form using various element types
   * (eg. P1, P2, Nedelec, and Raviart-Thomas).
   * This does integration over a surface or boundary (i.e., faces of your domain: faces in 3D, edges in 2D)
   * a(u,v)_i, where i denotes a set of faces (or edges) with in a boundary region marked with flag
   * DOES NOT take care of Dirichlet boundary conditions.  A separate routine will eliminate them later
   * This allows for several matrices to be assembled then added or concatenated together.
   *
   * For this problem we compute:
   *
   * Lu = f  ---->   a(u,v) = <f,v>
   *
   * which gives Ax = b,
   *
   * A_ij = a( phi_j, phi_i)
   * b_i  = <f,phi_i>
   * 
   *    INPUT:
   * local_assembly_face        Function to local assembly routine on face
   *             fespace        Finite-Element Space Struct
   *	      	    mesh        Mesh Struct
   *                  cq        Quadrature Coordinates and Weights
   *               coeff        Function that gives coefficient (for now assume constant)
   *                time        Time if time dependent
   *                flag        Marker for which faces are included in boundary integration
   *
   *    OUTPUT:
   *	         	A       Global Stiffness Matrix (CSR Format)
   */

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT f_per_elm = mesh->f_per_elm;
  INT dof_per_face = 0;
  INT dim = mesh->dim;
  if(FE->FEtype==1) { // P1 Elements
    dof_per_face = dim;
  } else if(FE->FEtype==-1) { // Nedelec Elements
    dof_per_face = 2*dim - 3; // 3 edges per face in 3D; face is edge in 2D
  } else if(FE->FEtype==-2) { // Raviart-Thomas Elements
    dof_per_elm = 1;
  } else {
    printf("Face integration isn't set up for the FEM space you chose\n");
    exit(0);
  }

  INT i,j,elm;

  // Allocate Row Array
  A->row = FE->ndof;
  A->col = FE->ndof;
  A->IA = (INT *) calloc(FE->ndof+1,sizeof(INT));

  // Get Sparsity Structure First
  // Non-zeros of A and IA (ignores cancellations, so maybe more than necessary)
  stiffG_nnz_subset(A,FE,flag);
	
  // Columns of A -> JA
  A->JA = (INT *) calloc(A->nnz,sizeof(INT));
  stiffG_cols_subset(A,FE,flag);
  
  // Set values
  A->val = (REAL *) calloc(A->nnz,sizeof(REAL));
  for (i=0; i<A->nnz; i++) {
    A->val[i] = 0;
  }
  	
  // Now Build Global Matrix entries
 
  /* Loop over all Faces and build local matrix */
  INT local_size = dof_per_face*dof_per_face;
  REAL* ALoc = (REAL *) calloc(local_size,sizeof(REAL));

  // Get mappings for given element and face
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT* dof_on_f = (INT *) calloc(dof_per_face,sizeof(INT));
  INT rowa,rowb,jcntr;

  // We will need the face to element map
  iCSRmat f_el;
  icsr_trans_1(mesh->el_f,&f_el);
  
  // Loop over boundary faces
  for (i=0; i<mesh->nface; i++) {	
    // Only grab the faces on the flagged boundary
    if(mesh->f_bdry[i]==flag) {
      // Zero out local matrices
      for (j=0; j<local_size; j++) {
	ALoc[j]=0;
      }

      // Find DOF for given Face
      rowa = FE->f_dof->IA[i]-1;
      rowb = FE->f_dof->IA[i+1]-1;
      jcntr = 0;
      for (j=rowa; j<rowb; j++) {
	dof_on_f[jcntr] = FE->f_dof->JA[j];
	jcntr++;
      }
      if(rowb-rowa!=dof_per_face) {
	printf("Face Integration has error.  Do you know how many DOF's per Face you have??\n");
	exit(0);
      }

      // Find the corresponding element associated with the face
      // Assume only 1 matters (if on boundary)
      rowa = f_el.IA[i]-1;
      elm = f_el.JA[rowa]-1;

      // Find DOF on that element
      rowa = FE->el_dof->IA[elm]-1;
      rowb = FE->el_dof->IA[elm+1]-1;
      jcntr = 0;
      for (j=rowa; j<rowb; j++) {
	dof_on_elm[jcntr] = FE->el_dof->JA[j];
	jcntr++;
      }
      
      // Find vertices for given Element
      rowa = mesh->el_v->IA[elm]-1;
      rowb = mesh->el_v->IA[elm+1]-1;
      jcntr = 0;
      for (j=rowa; j<rowb; j++) {
	v_on_elm[jcntr] = mesh->el_v->JA[j];
	jcntr++;
      }
		
      // Compute Local Stiffness Matrix for given Element
      (*local_assembly_face)(ALoc,FE,mesh,cq,dof_on_f,dof_on_elm,v_on_elm,i,elm,coeff,time);
      
      // Loop over DOF and place in appropriate slot globally
      LocaltoGlobal_face(dof_on_f,dof_per_face,FE,A,ALoc,flag); 
    }	
  }

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(dof_on_f) free(dof_on_f);
  if(ALoc) free(ALoc);
  icsr_free(&f_el);
  
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void assemble_global_Ned_GradH1_RHS(dvector *b,fespace *FE_H1,fespace *FE_Ned,trimesh *mesh,qcoordinates *cq,dvector* u) 
{
	
  /* Computes the global rhs for any <u,grad v>, where u is a Nedelec FE function and v is in P1
   * DOES NOT take care of Dirichlet boundary conditions.  A separate routine will eliminate them later
   *
   * 
   *    INPUT:
   *            FE_H1		    Finite-Element Space Struct for P1
   *            FE_Ned              Finite_element Space Struct for Nedelec
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights
   *            u                   Nedelec function for RHS
   *
   *    OUTPUT:
   *            b	       	    Global Right hand side vector
   */

  INT ed_per_elm = FE_Ned->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT i,j,row;

  // Allocate Arrays
  b->row = FE_H1->ndof;
  b->val = (REAL *) calloc(b->row,sizeof(REAL));
 
  /* Loop over all Elements and build local rhs */
  REAL* bLoc= (REAL *) calloc(v_per_elm,sizeof(REAL));

  INT* ed_on_elm = (INT *) calloc(ed_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT rowa,rowb,jcntr;
  for (i=0; i<FE_H1->nelm; i++) {	
    
    for (j=0; j<v_per_elm; j++) {
      bLoc[j]=0;
    }
		
    // Find Edges for given Element
    rowa = FE_Ned->el_dof->IA[i]-1;
    rowb = FE_Ned->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      ed_on_elm[jcntr] = FE_Ned->el_dof->JA[j];
      jcntr++;
    }

    // Find vertices for given Element
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }
		
    // Compute Local RHS for given Element
    Ned_GradH1_RHS_local(bLoc,FE_H1,FE_Ned,mesh,cq,ed_on_elm,v_on_elm,i,u);

    // Loop over DOF and place in appropriate slot globally
    for (j=0; j<v_per_elm; j++) { 
      row = v_on_elm[j]-1;	
      b->val[row] += bLoc[j];
    }
  }	

  if(ed_on_elm) free(ed_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(bLoc) free(bLoc);
		
  return;
}
/******************************************************************************************************/

