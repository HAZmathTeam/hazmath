/*!
 *
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program solves Biot's PDE for poroelasticity using finite elements
 * Locking-Free enrichecd Galerkin is used for the mechanics
 * Locally-conservative enriched Galerkin in used for the pressure
 *
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "hazmath.h"

#include "eg_stokes_params.h"
#include "eg_stokes_error.h"
#include "eg_stokes_system.h"
/****************************************************************/
static void csr2matlab_code(FILE *fp,dCSRmat *a, const char *varname)
{
  char *s;
  if((varname!="") && (varname))
    s=strndup(varname,8);
  else
    s=strdup("a");
  //
  fprintf(fp,"\n%%\n%s=[",s);
  csr_print_matlab(fp,a);
  fprintf(fp,"\n];\n%s=sparse(%s(:,1),%s(:,2),%s(:,3),%d,%d);",s,s,s,s,a->row,a->col);
  fprintf(fp,"\n%%\n");
  free(s);
  return;
}
/*********************************************************/
static void bdcsr_extend(block_dCSRmat *ba,				\
			 REAL *vcol, REAL *vrow,			\
			 const INT iblock,				\
			 const REAL c1, const REAL c2)
{
  /*
     if vrow and vcol are NULL, then adds one row and one column of
     zeros at block (iblock,iblock) also reallocs the corresponding
     blocks->IA, blocks->JA, and blocks->val. They are not reallocated
     (except blocks->IA) if vrow and vcol are null because we add
     column and row of zeros. --ltz 20200405
  */

  //
  INT k,j,n,nnz;
  REAL *vals=NULL;
  dCSRmat *atmp=NULL,*ctmp=NULL,*add=NULL;
  INT ptype=-1;
  FILE *fp;///=stdout;
  if(vcol==NULL && vrow==NULL){
    for(k=0;k<ba->brow;k++){
      // just add rows and columns of zeros to the matrix
      atmp=ba->blocks[iblock*ba->bcol+k];
      atmp->row++;
      atmp->IA=realloc(atmp->IA,sizeof(INT)*(atmp->row+1));
      atmp->IA[atmp->row]=atmp->nnz;
      atmp=ba->blocks[k*ba->bcol+iblock];
      atmp->col++;
    }
    return;
  }
  if(vcol==NULL){
    vals=vrow;
    ptype=2;
  } else if(vrow==NULL){
    vals=vcol;
    ptype=1;
  } else
    ptype=0;
  fprintf(stdout,"\n**** ptype=%d\n",ptype);
  if(ptype){
    // block rows
    atmp=ba->blocks[iblock*ba->bcol+iblock];
    n=atmp->row;
    nnz=n;
    // new matrix to add to the existing one
    add=dcsr_create_p(n,n,nnz);
    add->IA[0]=0;
    for(k=0;k<n;k++){
      add->JA[add->IA[k]]=k;
      add->val[add->IA[k]]=vals[k];
      add->IA[k+1]=add->IA[k]+1;
    }
  } else {
    for(k=0;k<ba->brow;k++){
      // block pressure rows
      atmp=ba->blocks[iblock*ba->bcol+k];
      atmp->row++;
      atmp->IA=realloc(atmp->IA,sizeof(INT)*(atmp->row+1));
      atmp->IA[atmp->row]=atmp->nnz;
      atmp=ba->blocks[k*ba->bcol+iblock];
      atmp->col++;
    }
    atmp=ba->blocks[iblock*ba->bcol+iblock];
    n=atmp->row;
    nnz=2*(n-1);
    /* fp=fopen("debug/x123.m","w"); */
    /* csr2matlab_code(fp,atmp,"xtmp"); */
    /* fclose(fp); */
    /* fprintf(stdout,"\nnnz=%d;ncol=%d,nrow=%d\n\n",atmp->nnz,atmp->col,atmp->row);       */
    /* fprintf(stdout,"\nn=%d;n=%d,nnz=%d\n\n",n,n,nnz); */
    /* fflush(stdout); */
    add=dcsr_create_p(n,n,nnz);
    add->IA[0]=0;
    add->val[0]=vcol[0];
    for(k=0;k<(n-1);k++){
      add->JA[add->IA[k]]=(n-1);
      add->val[add->IA[k]]=vcol[k];
      add->IA[k+1]=add->IA[k]+1;
    }
    add->IA[n]=nnz;
    //    fprintf(stdout,"\nn=%d,add->IA[n,n-1]=(%d,%d); (c1,c2)=(%e,%e)\n\n",n,add->IA[n],add->IA[n-1],c1,c2);
    j=0;
    for(k=add->IA[n-1];k<(add->IA[n]);k++){
      add->JA[k]=j;
      add->val[k]=vrow[j];
      j++;
    }
  }
  // now we add to atmp and then copy back to atmp;
  ctmp=malloc(1*sizeof(dCSRmat));
  // ctmp=c1*atmp+c2*add;
  dcsr_add(atmp,c1,add,c2,ctmp);
  free(add);
  atmp->IA=(INT *)realloc(atmp->IA,sizeof(INT)*(ctmp->row+1));
  atmp->JA=(INT *)realloc(atmp->JA,sizeof(INT)*ctmp->nnz);
  atmp->val=(REAL *)realloc(atmp->val,sizeof(REAL)*ctmp->nnz);
  // copy atmp=ctmp
  dcsr_cp(ctmp,atmp);
  /////////////////// print to debug/*
  /* fp=fopen("debug/app.m","w"); */
  /* csr2matlab_code(fp,atmp,"app"); */
  /* fclose(fp); */
  /* ///////////////////// */
  /* fprintf(stdout,"\nnnz-c=%d;nnz-a=%d\n",ctmp->nnz,atmp->nnz); fflush(stdout); */
  /* fp=fopen("debug/a.data","w"); */
  /* bdcsr_print_matlab(fp,ba); */
  /* fclose(fp); */
  //////////////////end print for debug
  // free the working memory.
  free(ctmp);
  return;
}
/**************************************************************/
void local_assembly_Elasticity_FACE(block_dCSRmat* A, block_fespace *FE, \
				    mesh_struct *mesh, qcoordinates *cq, \
				    INT face, iCSRmat *f_el,REAL time,	\
				    REAL timestep)
{

  //printf("Local Assemble Face - start\n");
  bool bEG = BOOL_EG_MECHANICS;
  bool bEG_Pressure = BOOL_EG_PRESSURE;
  REAL alpha = 1.;

  // INDEX
  INT i,j,k;
  INT dim = mesh->dim;
  INT rowa,rowb;

  // Mesh Stuff
  INT nspaces = FE->nspaces;

  INT v_per_elm = mesh->v_per_elm;
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT* v_on_elm_neighbor = (INT *) calloc(v_per_elm,sizeof(INT));

  INT test, trial;
  INT local_row_index, local_col_index;
  REAL kij = 0.0;

  INT dof_per_elm = 0;
  for(i=0;i<FE->nspaces;i++) {
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  }

  INT local_size = dof_per_elm*dof_per_elm;

  //Local Matrix for Each Cases on Faces
  REAL* ALoc = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* ALoc_u_vneighbor = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* ALoc_u_v = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* ALoc_uneighbor_v = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* ALoc_uneighbor_vneighbor = (REAL *) calloc(local_size,sizeof(REAL));

  for (j=0; j<local_size; j++) {
    ALoc[j]=0;
    ALoc_u_vneighbor[j]=0;
    ALoc_u_v[j]=0;
    ALoc_uneighbor_v[j]=0;
    ALoc_uneighbor_vneighbor[j]=0;
  }

  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* dof_on_elm_neighbor = (INT *) calloc(dof_per_elm,sizeof(INT));

  double lambda = LAME_LAMBDA_GLOBAL;

  INT quad,quad_face;
  INT jk, jkl,jcntr,ed;

  INT rowa_neighbor,rowb_neighbor,jcntr_neighbor;
  INT k_neighbor, j_neighbor;
  double fiarea;

  // For Locking-Free EG
  coordinates *barycenter = allocatecoords(1,dim);
  coordinates *barycenter_neighbor = allocatecoords(1,dim);


  // printf("===== FACE = %d  ===== \n",face);
  INT neighbor_index[2];
  neighbor_index[0] = -1;
  neighbor_index[1] = -1;
  INT counter = 0;

  //printf("=============================\n");
  //printf("** FACE = %d \n", face );
  //printf("=============================\n");
  int pq,nbr0;
  for(pq=f_el->IA[face];pq<f_el->IA[face+1];pq++){
    //printf("-- pq = %d, f_el->IA[face] = %d, f_el->IA[face+1] = %d \n", pq, f_el->IA[face], f_el->IA[face+1]);
    nbr0=f_el->JA[pq];
    //printf("-- nbr0 = %d  \n", nbr0);
    neighbor_index[counter] = nbr0;
    counter++;

  }
  // Sanity Check
  /*
    for(pq=0;pq<2;++pq)
    {
    if(counter == 2)
    {
    //printf("neighbor_index[%d]= %d || counter  = %d\n", pq, neighbor_index[pq],counter);
    }
    else if(counter == 1){
    if(pq == 0)
    {
    //printf("neighbor_index[%d]= %d || counter  = %d\n", pq, neighbor_index[pq],counter);
    }
    else{
    //printf("-\n");
    }
    }
    }
  */

  // Saniyy Check
  // counter == 1 means the BD
  if(counter == 1 && mesh->f_flag[face] == 0)
    {printf("check-error FALSE\n"); exit(0);}

  if(neighbor_index[0] == -1)
    {
      printf("Something is Wrong! \n");
      exit(0);
    }

  // Find DOF on element
  jcntr = 0;
  for(i=0;i<nspaces;i++) {
    rowa = FE->var_spaces[i]->el_dof->IA[neighbor_index[0]];
    rowb = FE->var_spaces[i]->el_dof->IA[neighbor_index[0]+1];
    for (j=rowa; j<rowb; j++) {
      dof_on_elm[jcntr] = FE->var_spaces[i]->el_dof->JA[j];
      jcntr++;
    }
  }

  // Find vertices for given Element
  get_incidence_row(neighbor_index[0],mesh->el_v,v_on_elm);

  // Find barycenter for given element
  barycenter->x[0] = mesh->el_mid[neighbor_index[0]*dim];
  barycenter->y[0] = mesh->el_mid[neighbor_index[0]*dim + 1];


  // Find DOF on the neighbor Element
  if(counter == 2)// in the interface
    {
      if(neighbor_index[1] == -1)
	{
	  printf("!#!@$#@!$ \n");
	  exit(0);
	}

      jcntr_neighbor = 0;
      for(k_neighbor=0;k_neighbor<FE->nspaces;k_neighbor++) {
	rowa_neighbor = FE->var_spaces[k_neighbor]->el_dof->IA[neighbor_index[1]];
	rowb_neighbor = FE->var_spaces[k_neighbor]->el_dof->IA[neighbor_index[1]+1];
	for (j_neighbor=rowa_neighbor; j_neighbor<rowb_neighbor; j_neighbor++) {
	  dof_on_elm_neighbor[jcntr_neighbor] = FE->var_spaces[k_neighbor]->el_dof->JA[j_neighbor];
	  jcntr_neighbor++;
	}
      }
      // Find vertices for given Element
      get_incidence_row(neighbor_index[1],mesh->el_v,v_on_elm_neighbor);

      // Find barycenter for given element
      barycenter_neighbor->x[0] = mesh->el_mid[neighbor_index[1]*dim];
      barycenter_neighbor->y[0] = mesh->el_mid[neighbor_index[1]*dim + 1];
    }

  // Sanity Check
  /*
    for(INT zz=0; zz<v_per_elm; ++zz)
    {
    printf("V_on_elm[%d] = %d  ||  dof_on_elm[%d] =%d \n", zz,v_on_elm[zz],zz,dof_on_elm[zz]);
    }

    if(counter ==2) // if it has neighbor
    for(INT zz=0; zz<v_per_elm; ++zz)
    {
    printf("V_on_elm_neighbor[%d] = %d  ||  dof_on_elm_neighbor[%d] =%d \n", zz,v_on_elm_neighbor[zz],zz,dof_on_elm_neighbor[zz]);
    }

    printf("FACE = %d, ELM = %d, (X,Y) = (%f,%f) \n", face, neighbor_index[0], barycenter->x[0], barycenter->y[0]);
    if(counter ==2) // if it has neighbor
    printf("FACE = %d, ELM = %d, (X,Y) = (%f,%f) \n", face, neighbor_index[1], barycenter_neighbor->x[0], barycenter_neighbor->y[0]);
  */

  INT* local_dof_on_elm_face = NULL;
  INT *local_dof_on_elm_face_interface = NULL;
  INT *local_dof_on_elm_neighbor = NULL;

  REAL *data_face=calloc(dim+dim*dim, sizeof(REAL)); //coords of normal and vertices of a face.
  REAL *xfi=data_face; //coords of vertices on face i.
  REAL *finrm=xfi+dim*dim; //coords of normal vector on face i
  //REAL *data_face_end=finrm + dim; //end

  // Get normal vector values.
  finrm[0]=mesh->f_norm[face*dim+0];
  finrm[1]=mesh->f_norm[face*dim+1];

  //Sanity Check; Print the normal vectors..
  //printf("FACE = %d, ELM = %d, (Nx,Ny) = (%f,%f) \n", face, neighbor_index[0], finrm[0], finrm[1]);

  INT nq1d_face=3;
  qcoordinates *cq_face = allocateqcoords_bdry(nq1d_face,1,dim,2);

  INT maxdim=2;
  REAL qx_face[maxdim];
  REAL w_face;
  INT nquad_face= cq_face->nq_per_elm; //quad<cq_face->nq_per_elm;


  // Get xfi points..
  for(jkl=mesh->f_v->IA[face];jkl<mesh->f_v->IA[face+1];jkl++){
    //printf("M** jkl = %d, mesh->f_v->IA[face] = %d, mesh->f_v->IA[face+1] = %d \n", jkl, mesh->f_v->IA[face], mesh->f_v->IA[face+1]);
    j=jkl-mesh->f_v->IA[face];
    k=mesh->f_v->JA[jkl];

    xfi[j*dim+0]=mesh->cv->x[k];
    xfi[j*dim+1]=mesh->cv->y[k];
    //printf("M** xfi[j*dim+0] = %f,  xfi[j*dim+1] = %f \n",  xfi[j*dim+0] ,  xfi[j*dim+1]);
  }

  //define the length of the face-
  fiarea=mesh->f_area[face];
  //printf("FACE = %d, ELM = %d, fiarea = %f \n", face, neighbor_index[0], fiarea);

  // get the quadrautre
  zquad_face(cq_face,nq1d_face,dim,xfi,fiarea);

  // get the penalty (FACE ASSEMBLE)
  double penalty_term = PENALTY_PARAMETER_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));
  double BC_penalty_term = BC_PENALTY_PARAMETER_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));

  double penalty_term_pressure = PENALTY_PARAMETER_PRESSURE_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));
  double BC_penalty_term_pressure = BC_PENALTY_PARAMETER_PRESSURE_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));


  // DEBUG mode
  // Get the BD values
  //REAL* val_true_face = (REAL *) calloc(nquad_face,sizeof(REAL));

  for (quad_face=0;quad_face<nquad_face;quad_face++) {

    qx_face[0] = cq_face->x[quad_face];
    qx_face[1] = cq_face->y[quad_face];

    if(dim==3) qx_face[2] = cq_face->z[quad_face];

    w_face = cq_face->w[quad_face];

    //DEBUG mode
    //(*exact_sol2D)(val_true_face,qx_face,time,
    //	     &(mesh->el_flag[neighbor_index[0]])); //DUMMY VALUE at the end // ???


    ///////////////////////////////////////////////////
    if(counter == 1){ //only at boundary
      //printf("== BD = %d \n", face);
      //Sanity Check
      if(mesh->f_flag[face] == 0)
	{printf("check-error0\n"); exit(0);}

      bool bWeakBC_FACE = BOOL_WEAKLY_IMPOSED_BC;

      if(bWeakBC_FACE){

	get_FEM_basis(FE->var_spaces[0]->phi,
		      FE->var_spaces[0]->dphi,
		      qx_face,
		      v_on_elm,
		      dof_on_elm,
		      mesh,
		      FE->var_spaces[0]);

	local_dof_on_elm_face = dof_on_elm + FE->var_spaces[0]->dof_per_elm;

	get_FEM_basis(FE->var_spaces[1]->phi,
		      FE->var_spaces[1]->dphi,
		      qx_face,
		      v_on_elm,
		      local_dof_on_elm_face,
		      mesh,
		      FE->var_spaces[1]);

	local_dof_on_elm_face += FE->var_spaces[1]->dof_per_elm;

	//for the EG part
	local_dof_on_elm_face += FE->var_spaces[2]->dof_per_elm;


	get_FEM_basis(FE->var_spaces[3]->phi,
		      FE->var_spaces[3]->dphi,
		      qx_face,
		      v_on_elm,
		      local_dof_on_elm_face,
		      mesh,
		      FE->var_spaces[3]);


	//for the P part
	local_dof_on_elm_face += FE->var_spaces[3]->dof_per_elm;

	//---- #2 BOUNDARY ---------------------------- START //
	local_row_index = 0;
	local_col_index = 0;
	// u0-v0 block:
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	    kij = -2.0 * FE->var_spaces[0]->dphi[trial*dim]* finrm[0] * FE->var_spaces[0]->phi[test];
	    kij -= FE->var_spaces[0]->dphi[trial*dim+1] * finrm[1] *    FE->var_spaces[0]->phi[test];

	    kij -= lambda * FE->var_spaces[0]->dphi[trial*dim] * finrm[0] * FE->var_spaces[0]->phi[test];

	    //penalty term
	    kij += BC_penalty_term * FE->var_spaces[0]->phi[trial] * FE->var_spaces[0]->phi[test];

	    ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face *kij;
	  }
	}

	// u0-v1 block <dy(u1),dx(v2)>
	local_row_index = FE->var_spaces[0]->dof_per_elm;
	local_col_index = 0;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	    kij = -FE->var_spaces[0]->dphi[trial*dim+1] * finrm[0] * FE->var_spaces[1]->phi[test];

	    kij -= lambda * FE->var_spaces[0]->dphi[trial*dim] * finrm[1] * FE->var_spaces[1]->phi[test];

	    ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// u1-v0 block :
	local_row_index = 0;
	local_col_index = FE->var_spaces[0]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	    kij = - FE->var_spaces[1]->dphi[trial*dim] * finrm[1] * FE->var_spaces[0]->phi[test];

	    kij -= lambda * FE->var_spaces[1]->dphi[trial*dim+1] * finrm[0] * FE->var_spaces[0]->phi[test];

	    ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// u1-v1 block
	local_row_index = FE->var_spaces[0]->dof_per_elm;
	local_col_index = FE->var_spaces[0]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	    kij = -2.0 * FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * FE->var_spaces[1]->phi[test] ;
	    kij -= 1.0 * FE->var_spaces[1]->dphi[trial*dim] * finrm[0]   * FE->var_spaces[1]->phi[test];

	    // NEW SLEE : elasticity div part
	    // u2-v2 block : <dy(u2),dx(v2)>
	    kij -= lambda * FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * FE->var_spaces[1]->phi[test];

	    //penalty term
	    kij += BC_penalty_term *  FE->var_spaces[1]->phi[trial]  * FE->var_spaces[1]->phi[test];

	    ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	//q-q block:
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){

	    kij = -2.0 * 1. * finrm[0] * (qx_face[0] - barycenter->x[0]);
	    kij -= 2.0 * 1. * finrm[1] * (qx_face[1] - barycenter->y[0]);

	    // NEW SLEE : elasticity div part
	    // u2-v1 block : <dy(u2),dx(v1)>
	    kij -= lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);
	    kij -= lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);

	    kij -= lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);
	    kij -= lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);


	    //penalty term
	    kij += BC_penalty_term * (qx_face[0]- barycenter->x[0]) * (qx_face[0]- barycenter->x[0]);
	    kij += BC_penalty_term * (qx_face[1]- barycenter->y[0]) * (qx_face[1]- barycenter->y[0]);


	    ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Face Boundary
	if(bEG){
	  //EG PART
	  //NEW SLEE :
	  // u0- q block
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
	  local_col_index = 0;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	      //u0-q0
	      //u0-q1
	      kij = -2.0*FE->var_spaces[0]->dphi[trial*dim] * finrm[0]  * (qx_face[0]- barycenter->x[0]); //FE->var_spaces[dim]->phi[test];
	      kij -= FE->var_spaces[0]->dphi[trial*dim+1] * finrm[1] * (qx_face[0]- barycenter->x[0]);

	      kij -= FE->var_spaces[0]->dphi[trial*dim+1] * finrm[0]  * (qx_face[1]- barycenter->y[0]); //FE->var_spaces[dim]->phi[test];

	      // Divergence
	      // u1-q block : <dx(u1),dx(v1)>
	      kij -= lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[0] * (qx_face[0]- barycenter->x[0]);//FE->var_spaces[dim]->dphi[test*dim];
	      kij -= lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[1] * (qx_face[1]- barycenter->y[0]);

	      //penalty term
	      kij += penalty_term * FE->var_spaces[0]->phi[trial] * (qx_face[0]- barycenter->x[0]);

	      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	  // u1-q block
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
	  local_col_index = FE->var_spaces[0]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	      //u1-q0
	      //u1-q1
	      kij = -2.0*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * (qx_face[1]- barycenter->y[0]);
	      kij -= FE->var_spaces[1]->dphi[trial*dim] * finrm[0] * (qx_face[1]- barycenter->y[0]);

	      kij -= FE->var_spaces[1]->dphi[trial*dim] * finrm[1] * (qx_face[0]- barycenter->x[0]);

	      kij -= lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[0] *  (qx_face[0]- barycenter->x[0]);
	      kij -= lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] *  (qx_face[1]- barycenter->y[0]);


	      //penalty term
	      kij += penalty_term * FE->var_spaces[1]->phi[trial] * (qx_face[1]- barycenter->y[0]);

	      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	  // q-v0 block:
	  local_row_index = 0;
	  local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){

	      //q0 - v0
	      //q1 - v0
	      kij = -2. * 1. * finrm[0] * FE->var_spaces[0]->phi[test];

	      // Divergence
	      // u1-q block : <dx(u1),dx(v1)>
	      kij -= lambda* 1. * finrm[0] * FE->var_spaces[0]->phi[test];//FE->var_spaces[dim]->dphi[test*dim];
	      kij -= lambda* 1. * finrm[0] * FE->var_spaces[0]->phi[test];


	      //penalty term
	      kij += penalty_term * FE->var_spaces[0]->phi[test] * (qx_face[0]- barycenter->x[0]);

	      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }


	  // q-v1 block:
	  local_row_index = FE->var_spaces[0]->dof_per_elm;
	  local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){

	      //q0-v1
	      //q1-v1
	      kij = -2.0* 1. * finrm[1] * FE->var_spaces[1]->phi[test];

	      // NEW SLEE : elasticity div part
	      // u2-v1 block : <dy(u2),dx(v1)>
	      kij -= lambda * 1. * finrm[1] * FE->var_spaces[1]->phi[test];
	      kij -= lambda * 1. * finrm[1] * FE->var_spaces[1]->phi[test];

	      //penalty term
	      kij += penalty_term * FE->var_spaces[1]->phi[test] * (qx_face[1]- barycenter->y[0]);

 	      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	} // bEG
	//----- #2 BOUNDARY ------------------------------------------- END //


	//---- Poro P2 Boundary -------------------------------------- START //
	///////////////////////////////////////////////////////////////////////////
	// porop2
	// On BOUNDARY
	// PoroElasticity   + ({p},  [u] dot n)
	// Mechanics coupeld with pressure
	//////////////////////////////////////////////////////////////////////////
	//      ({p},  [u] dot n)
	//       CG  CG u0
	//       trial (column) 3  test (row) 0
	// test - u0 cg
	local_row_index = 0;
	// trial - p cg
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * FE->var_spaces[3]->phi[trial] *  FE->var_spaces[0]->phi[test]*finrm[0];

	    ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}
	//      ({p},  [u] dot n)
	//       CG  CG u1
	//       trial (column) 3  test (row) 1
	// test - u0 cg
	local_row_index = FE->var_spaces[0]->dof_per_elm;
	// trial - p cg
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * FE->var_spaces[3]->phi[trial] *  FE->var_spaces[1]->phi[test]*finrm[1];

	    ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	if(bEG){

	  //      ({p},  [u] dot n)
	  //       CG  EG
	  //       trial (column) 3  test (row) 2
	  // test - u0 eg
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  // trial - p cg
	  local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	      kij = alpha * FE->var_spaces[3]->phi[trial] *
		(finrm[0] *  (qx_face[0]- barycenter->x[0]) + finrm[1] *  (qx_face[1]- barycenter->y[0]));

	      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	}

	//---- Poro P2 Boundary -------------------------------------- END //


	//---- Stokes P2 Boundary -------------------------------------- START //
	///////////////////////////////////////////////////////////////////////////
	// On BOUNDARY
	// STokes  + ({p},  [u] dot n)
	//////////////////////////////////////////////////////////////////////////
	//      ({p},  [u] dot n)
	//       CG  CG u0
	//       trial (column) 3  test (row) 0
	// test - p cg
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// trial - u0 cg
	local_col_index = 0;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	    kij = alpha * FE->var_spaces[3]->phi[test] *  FE->var_spaces[0]->phi[trial]*finrm[0];

	    ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}
	//      ({p},  [u] dot n)
	//       CG  CG u1
	//       trial (column) 3  test (row) 1
	// test - p cg
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// trial - u0 cg
	local_col_index = FE->var_spaces[0]->dof_per_elm;

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	    kij = alpha * FE->var_spaces[3]->phi[test] *  FE->var_spaces[1]->phi[trial]*finrm[1];

	    ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	if(bEG){

	  //      ({p},  [u] dot n)
	  //       CG  EG
	  //       trial (column) 3  test (row) 2
	  // test - p cg
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	  // trial - u0 cg
	  local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){

	      kij = alpha * FE->var_spaces[3]->phi[test] *
		(finrm[0] *  (qx_face[0]- barycenter->x[0]) + finrm[1] *  (qx_face[1]- barycenter->y[0]));

	      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	}

	//---- Stokes P2 Boundary -------------------------------------- END //

      }// Weakly Face

    }
    else if(counter ==2)
      { //on interface

	local_dof_on_elm_face_interface = NULL;
	local_dof_on_elm_neighbor = NULL;
	// HERE IN THE INTERFACE, WE HAVE FOUR DIFFERENT LOCAL MATRICES
	// CASE I   u v
	// CASE II  u v_neighbor
	// CASE III u_neighbor v
	// CASE IV  u_neighbor v_neighbor
	//Sanity Cehck
	if(mesh->f_flag[face] != 0)
	  {printf("well? something is wrong \n"); exit(0);}


	//printf("counter 2 \n");

	get_FEM_basis(FE->var_spaces[0]->phi,
		      FE->var_spaces[0]->dphi,
		      qx_face,
		      v_on_elm,
		      dof_on_elm,
		      mesh,
		      FE->var_spaces[0]);

	local_dof_on_elm_face_interface = dof_on_elm + FE->var_spaces[0]->dof_per_elm;


	get_FEM_basis(FE->var_spaces[1]->phi,
		      FE->var_spaces[1]->dphi,
		      qx_face,
		      v_on_elm,
		      local_dof_on_elm_face_interface,
		      mesh,
		      FE->var_spaces[1]);

	local_dof_on_elm_face_interface += FE->var_spaces[1]->dof_per_elm;

	//FOR EG
	local_dof_on_elm_face_interface += FE->var_spaces[2]->dof_per_elm;

	get_FEM_basis(FE->var_spaces[3]->phi,
		      FE->var_spaces[3]->dphi,
		      qx_face,
		      v_on_elm,
		      local_dof_on_elm_face_interface,
		      mesh,
		      FE->var_spaces[3]);

	//FOR P
	local_dof_on_elm_face_interface += FE->var_spaces[3]->dof_per_elm;

	// Now, get the basis for the neighbor
	REAL* neighbor_basis_0_phi  = (REAL *) calloc(3,sizeof(REAL));
	REAL* neighbor_basis_0_dphi = (REAL *) calloc(3*mesh->dim,sizeof(REAL));

	REAL* neighbor_basis_1_phi  = (REAL *) calloc(3,sizeof(REAL));
	REAL* neighbor_basis_1_dphi = (REAL *) calloc(3*mesh->dim,sizeof(REAL));

	REAL* neighbor_basis_3_phi  = (REAL *) calloc(3,sizeof(REAL));
	REAL* neighbor_basis_3_dphi = (REAL *) calloc(3*mesh->dim,sizeof(REAL));


	get_FEM_basis( neighbor_basis_0_phi ,
		       neighbor_basis_0_dphi ,
		       qx_face,
		       v_on_elm_neighbor,
		       dof_on_elm_neighbor,
		       mesh,
		       FE->var_spaces[0]);

	local_dof_on_elm_neighbor = dof_on_elm_neighbor + FE->var_spaces[0]->dof_per_elm;

	//printf("counter 111 \n");

	get_FEM_basis( neighbor_basis_1_phi ,
		       neighbor_basis_1_dphi ,
		       qx_face,
		       v_on_elm_neighbor,
		       local_dof_on_elm_neighbor,
		       mesh,
		       FE->var_spaces[1]);

	// u2
	local_dof_on_elm_neighbor += FE->var_spaces[1]->dof_per_elm;

	//printf("counter 222 \n");
	// u eg
	local_dof_on_elm_neighbor += FE->var_spaces[2]->dof_per_elm;


	get_FEM_basis( neighbor_basis_3_phi ,
		       neighbor_basis_3_dphi ,
		       qx_face,
		       v_on_elm_neighbor,
		       local_dof_on_elm_neighbor,
		       mesh,
		       FE->var_spaces[3]);

	// p
	local_dof_on_elm_neighbor += FE->var_spaces[3]->dof_per_elm;


	// ------ #2 INTERFACE -------------------------- START //
	//q-q block:
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){

	    kij = -0.5* 2.0 * 1. * finrm[0] * (qx_face[0]- barycenter->x[0]);
	    kij -= 0.5* 2.0 * 1. * finrm[1] * (qx_face[1]- barycenter->y[0]);

	    // NEW SLEE : elasticity div part
	    // u2-v1 block : <dy(u2),dx(v1)>
	    kij -= 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);
	    kij -= 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);

	    kij -= 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);
	    kij -= 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);

	    //penalty term
	    kij += penalty_term * (qx_face[0]- barycenter->x[0]) * (qx_face[0]- barycenter->x[0]);
	    kij += penalty_term * (qx_face[1]- barycenter->y[0]) * (qx_face[1]- barycenter->y[0]);

	    ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	//q-neighbor q block: where q is from neighbor
	// column (trial is from neighbor)
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){
	    if(test == 1)
	      {
		printf("oops !? \n");
		exit(0);
	      }

	    kij = 0.5* -2.0 * 1. * finrm[0] * (qx_face[0]- barycenter->x[0]);
	    kij -= 0.5*2.0 * 1. * finrm[1] * (qx_face[1]- barycenter->y[0]);

	    // NEW SLEE : elasticity div part
	    // u2-v1 block : <dy(u2),dx(v1)>
	    kij -= 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);
	    kij -= 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);

	    kij -= 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);
	    kij -= 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);

	    //penalty_term
	    kij += -penalty_term * (qx_face[0]- barycenter_neighbor->x[0]) * (qx_face[0]- barycenter->x[0]) ;
	    kij += -penalty_term * (qx_face[1]- barycenter_neighbor->y[0]) * (qx_face[1]- barycenter->y[0]) ;



	    ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}


	//q q-neighbor block: where q (Test) is from neighbor => we need minus sign
	// row (test is from neighbor)
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){
	    if(test == 1)
	      {
		printf("oops !? \n");
		exit(0);
	      }

	    kij = 0.5* 2.0 * 1. * finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);
	    kij += 0.5*2.0 * 1. * finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);

	    // NEW SLEE : elasticity div part
	    // u2-v1 block : <dy(u2),dx(v1)>
	    kij += 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);
	    kij += 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);

	    kij += 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);
	    kij += 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);

	    //penalty_term
	    kij += -penalty_term * (qx_face[0]- barycenter_neighbor->x[0]) * (qx_face[0]- barycenter->x[0]) ;
	    kij += -penalty_term * (qx_face[1]- barycenter_neighbor->y[0]) * (qx_face[1]- barycenter->y[0]) ;


	    ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	//q-neighobr q-neighbor block: where q (Test) is from neighbor => we need minus sign
	// column (trial  is from neighbor)
	// row (test is from neighbor)
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){
	    if(test == 1)
	      {
		printf("oops !? \n");
		exit(0);
	      }

	    kij = 0.5* 2.0 * 1. * finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);
	    kij += 0.5*2.0 * 1. * finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);

	    // NEW SLEE : elasticity div part
	    // u2-v1 block : <dy(u2),dx(v1)>
	    kij += 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);
	    kij += 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);

	    kij += 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);
	    kij += 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);

	    //penalty_term
	    kij += penalty_term * (qx_face[0]- barycenter_neighbor->x[0]) * (qx_face[0]- barycenter_neighbor->x[0]) ;
	    kij += penalty_term * (qx_face[1]- barycenter_neighbor->y[0]) * (qx_face[1]- barycenter_neighbor->y[0]) ;


	    ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}


	if(bEG){
	  // u0 - q block
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = 0;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	      //CG-DG
	      // u0 - q0
	      // u0 - q1
	      kij = -0.5 * 2.0*FE->var_spaces[0]->dphi[trial*dim] * finrm[0]  * (qx_face[0]- barycenter->x[0]); //FE->var_spaces[dim]->phi[test];
	      kij -= 0.5 * FE->var_spaces[0]->dphi[trial*dim+1] * finrm[1]    * (qx_face[0]- barycenter->x[0]);
	      kij -= 0.5 * FE->var_spaces[0]->dphi[trial*dim+1] * finrm[0]    * (qx_face[1]- barycenter->y[0]); //FE->var_spaces[dim]->phi[test];

	      // Divergence
	      // u1-q block :
	      kij -= 0.5 * lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[0] * (qx_face[0]- barycenter->x[0]);//FE->var_spaces[dim]->dphi[test*dim];
	      kij -= 0.5 * lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[1] * (qx_face[1]- barycenter->y[0]);

	      ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	  // u1-q block
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = FE->var_spaces[0]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	      kij = -0.5 * 2.0*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * (qx_face[1]- barycenter->y[0]);
	      kij -= 0.5 * FE->var_spaces[1]->dphi[trial*dim] * finrm[0]       * (qx_face[1]- barycenter->y[0]);
	      kij -= 0.5 * FE->var_spaces[1]->dphi[trial*dim] * finrm[1]       * (qx_face[0]- barycenter->x[0]);

	      // NEW SLEE : elasticity div part
	      // u2-v1 block : <dy(u2),dx(v1)>
	      kij -= 0.5 *lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[0] *  (qx_face[0]- barycenter->x[0]);
	      kij -= 0.5 *lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] *  (qx_face[1]- barycenter->y[0]);

	      ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }


	  // u0 - q (where u0 is from neighbor)
	  // column (trial is from neighbor)
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = 0;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	      //DEBUG
	      kij = -0.5 * 2.0*neighbor_basis_0_dphi[trial*dim] * finrm[0]  * (qx_face[0]- barycenter->x[0]); //FE->var_spaces[dim]->phi[test];
	      kij -= 0.5 * neighbor_basis_0_dphi[trial*dim+1] * finrm[1]    * (qx_face[0]- barycenter->x[0]);
	      kij -= 0.5 * neighbor_basis_0_dphi[trial*dim+1] * finrm[0]    * (qx_face[1]- barycenter->y[0]); //FE->var_spaces[dim]->phi[test];

	      // Divergence
	      // u1-q block :
	      kij -= 0.5 * lambda*neighbor_basis_0_dphi[trial*dim]* finrm[0] * (qx_face[0]- barycenter->x[0]);//FE->var_spaces[dim]->dphi[test*dim];
	      kij -= 0.5 * lambda*neighbor_basis_0_dphi[trial*dim]* finrm[1] * (qx_face[1]- barycenter->y[0]);

	      ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	  // u1-q (where u1 is from neighbor)
	  // column (trial is from neighbor)
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = FE->var_spaces[0]->dof_per_elm;

	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
	      //DEBUG
	      kij = -0.5 * 2.0*neighbor_basis_1_dphi[trial*dim+1] * finrm[1] * (qx_face[1]- barycenter->y[0]);
	      kij -= 0.5 * neighbor_basis_1_dphi[trial*dim]       * finrm[0] * (qx_face[1]- barycenter->y[0]);
	      kij -= 0.5 * neighbor_basis_1_dphi[trial*dim]       * finrm[1] * (qx_face[0]- barycenter->x[0]);

	      // NEW SLEE : elasticity div part
	      // u2-v1 block : <dy(u2),dx(v1)>
	      kij -= 0.5 *lambda*neighbor_basis_1_dphi[trial*dim+1] * finrm[0] *  (qx_face[0]- barycenter->x[0]);
	      kij -= 0.5 *lambda*neighbor_basis_1_dphi[trial*dim+1] * finrm[1] *  (qx_face[1]- barycenter->y[0]);

	      ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	  // u0 - q_neighbor block
	  // (q has minus sign)
	  // row (test is from neighbor)
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = 0;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	      //CG-DG
	      // u0 - q0
	      // u0 - q1
	      kij = 0.5 * 2.0*FE->var_spaces[0]->dphi[trial*dim] * finrm[0]  * (qx_face[0]- barycenter_neighbor->x[0]); //FE->var_spaces[dim]->phi[test];
	      kij += 0.5 * FE->var_spaces[0]->dphi[trial*dim+1] * finrm[1] * (qx_face[0]- barycenter_neighbor->x[0]);
	      kij += 0.5 * FE->var_spaces[0]->dphi[trial*dim+1] * finrm[0]  * (qx_face[1]- barycenter_neighbor->y[0]); //FE->var_spaces[dim]->phi[test];

	      // Divergence
	      // u1-q block :
	      kij += 0.5 * lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);//FE->var_spaces[dim]->dphi[test*dim];
	      kij += 0.5 * lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);

	      ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	  // u1-q_neighbor block
	  // (q has minus sign)
	  // row (test is from neighbor)
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = FE->var_spaces[0]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	      kij = 0.5 * 2.0*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);
	      kij += 0.5 * FE->var_spaces[1]->dphi[trial*dim] * finrm[0] * (qx_face[1]- barycenter_neighbor->y[0]);
	      kij += 0.5 * FE->var_spaces[1]->dphi[trial*dim] * finrm[1] * (qx_face[0]- barycenter_neighbor->x[0]);

	      // NEW SLEE : elasticity div part
	      // u2-v1 block : <dy(u2),dx(v1)>
	      kij += 0.5 *lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[0] *  (qx_face[0]- barycenter_neighbor->x[0]);
	      kij += 0.5 *lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] *  (qx_face[1]- barycenter_neighbor->y[0]);


	      ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	  // u0-neighbor - q_neighbor block
	  // (q has minus sign)
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = 0;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	      //CG-DG
	      // u0 - q0
	      // u0 - q1
	      kij = 0.5 * 2.0*neighbor_basis_0_dphi[trial*dim] * finrm[0]  * (qx_face[0]- barycenter_neighbor->x[0]); //FE->var_spaces[dim]->phi[test];
	      kij += 0.5 * neighbor_basis_0_dphi[trial*dim+1] * finrm[1] * (qx_face[0]- barycenter_neighbor->x[0]);

	      kij += 0.5 * neighbor_basis_0_dphi[trial*dim+1] * finrm[0]  * (qx_face[1]- barycenter_neighbor->y[0]); //FE->var_spaces[dim]->phi[test];

	      // Divergence
	      // u1-q block :
	      kij += 0.5 * lambda*neighbor_basis_0_dphi[trial*dim]* finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);//FE->var_spaces[dim]->dphi[test*dim];
	      kij += 0.5 * lambda*neighbor_basis_0_dphi[trial*dim]* finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);

	      ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	  // u1-neighbor - q_neighbor block
	  // (q has minus sign)
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = FE->var_spaces[0]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	      kij = 0.5 * 2.0*neighbor_basis_1_dphi[trial*dim+1] * finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);
	      kij += 0.5 * neighbor_basis_1_dphi[trial*dim] * finrm[0] * (qx_face[1]- barycenter_neighbor->y[0]);
	      kij += 0.5 * neighbor_basis_0_dphi[trial*dim] * finrm[1] * (qx_face[0]- barycenter_neighbor->x[0]);

	      // NEW SLEE : elasticity div part
	      // u2-v1 block : <dy(u2),dx(v1)>
	      kij += 0.5 *lambda*neighbor_basis_0_dphi[trial*dim+1] * finrm[0] *  (qx_face[0]- barycenter_neighbor->x[0]);
	      kij += 0.5 *lambda*neighbor_basis_0_dphi[trial*dim+1] * finrm[1] *  (qx_face[1]- barycenter_neighbor->y[0]);


	      ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }
	}// bEG TERMS.
	// ------ #2 INTERFACE -------------------------- END //


	// ------ PORO P1 INTERFACE -------------------------- START //
	///////////////////////////////////////////////////////////////////////////
	//porop1
	//PoroElasticity   + ({p},  [u] dot n)
	//////////////////////////////////////////////////////////////////////////
	//      ({p},  [u] dot n)
	//       CG  CG
	//       trial (column) 3  test (row) 0
	// test - u0 cg

	local_row_index = 0;
	// trial - p cg
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * FE->var_spaces[3]->phi[trial] *  FE->var_spaces[0]->phi[test]*finrm[0];

	    ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * FE->var_spaces[3]->phi[trial] *  neighbor_basis_0_phi[test]*finrm[0];

	    ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * neighbor_basis_3_phi[trial] *  FE->var_spaces[0]->phi[test]*finrm[0];

	    ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * neighbor_basis_3_phi[trial] *  neighbor_basis_0_phi[test]*finrm[0];

	    ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}


	//      ({p},  [u] dot n)
	//       CG  CG u1
	//       trial (column) 3  test (row) 1
	// test - u1 cg
	local_row_index = FE->var_spaces[0]->dof_per_elm;
	// trial - p cg
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * FE->var_spaces[3]->phi[trial] *  FE->var_spaces[1]->phi[test]*finrm[1];

	    ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * FE->var_spaces[3]->phi[trial] *  neighbor_basis_1_phi[test]*finrm[1];

	    ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * neighbor_basis_3_phi[trial] *  FE->var_spaces[1]->phi[test]*finrm[1];

	    ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * neighbor_basis_3_phi[trial] *  neighbor_basis_1_phi[test]*finrm[1];

	    ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}




	//porop1
	//      ({p},  [u] dot n)
	//       CG  DG
	//       trial (column) 3  test (row) 2
	// test - u eg
	local_row_index = FE->var_spaces[0]->dof_per_elm +FE->var_spaces[1]->dof_per_elm;
	// trial - p cg
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * FE->var_spaces[3]->phi[trial] *
	      (finrm[0] *  (qx_face[0]- barycenter->x[0]) + finrm[1] *  (qx_face[1]- barycenter->y[0]));
	    ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * FE->var_spaces[3]->phi[trial] *
	      (finrm[0] *  (qx_face[0]- barycenter_neighbor->x[0]) + finrm[1] *  (qx_face[1]- barycenter_neighbor->y[0]));

	    ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * neighbor_basis_3_phi[trial] *
	      (finrm[0] *  (qx_face[0]- barycenter->x[0]) + finrm[1] *  (qx_face[1]- barycenter->y[0]));

	    ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * neighbor_basis_3_phi[trial] *
	      (finrm[0] *  (qx_face[0]- barycenter_neighbor->x[0]) + finrm[1] *  (qx_face[1]- barycenter_neighbor->y[0]));

	    ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}



	// ------ PORO P1 INTERFACE -------------------------- END//


	// ------ STOKES P1 INTERFACE -------------------------- START //
	///////////////////////////////////////////////////////////////////////////
	// Stokes   + ({p},  [u] dot n)
	//////////////////////////////////////////////////////////////////////////
	//      ({p},  [u] dot n)
	//       CG  CG
	//       trial (column) 3  test (row) 0
	// trial - u0 cg
	local_col_index = 0;
	// test - p cg
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * FE->var_spaces[3]->phi[test] *  FE->var_spaces[0]->phi[trial]*finrm[0];

	    ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * FE->var_spaces[3]->phi[test] *  neighbor_basis_0_phi[trial]*finrm[0];

	    ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * neighbor_basis_3_phi[test] *  FE->var_spaces[0]->phi[trial]*finrm[0];

	    ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * neighbor_basis_3_phi[test] *  neighbor_basis_0_phi[trial]*finrm[0];

	    ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}


	//      ({p},  [u] dot n)
	//       CG  CG u1
	//       trial (column) 3  test (row) 1
	// trial - u1 cg
	local_col_index = FE->var_spaces[0]->dof_per_elm;
	// test - p cg
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * FE->var_spaces[3]->phi[test] *  FE->var_spaces[1]->phi[trial]*finrm[1];

	    ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * FE->var_spaces[3]->phi[test] *  neighbor_basis_1_phi[trial]*finrm[1];

	    ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * neighbor_basis_3_phi[test] *  FE->var_spaces[1]->phi[trial]*finrm[1];

	    ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * neighbor_basis_3_phi[test] *  neighbor_basis_1_phi[trial]*finrm[1];

	    ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}


	//Stokes1
	//      ({p},  [u] dot n)
	//       CG  DG
	//       test (column) 3  trial (row) 2
	// trial - u eg
	local_col_index = FE->var_spaces[0]->dof_per_elm +FE->var_spaces[1]->dof_per_elm;
	// test - p cg
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * FE->var_spaces[3]->phi[test] *
	      (finrm[0] *  (qx_face[0]- barycenter->x[0]) + finrm[1] *  (qx_face[1]- barycenter->y[0]));
	    ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * FE->var_spaces[3]->phi[test] *
	      (finrm[0] *  (qx_face[0]- barycenter_neighbor->x[0]) + finrm[1] *  (qx_face[1]- barycenter_neighbor->y[0]));

	    ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * neighbor_basis_3_phi[test] *
	      (finrm[0] *  (qx_face[0]- barycenter->x[0]) + finrm[1] *  (qx_face[1]- barycenter->y[0]));

	    ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * neighbor_basis_3_phi[test] *
	      (finrm[0] *  (qx_face[0]- barycenter_neighbor->x[0]) + finrm[1] *  (qx_face[1]- barycenter_neighbor->y[0]));

	    ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}


	//Stokes1


	// Stokes Interafec -- end //



	free(neighbor_basis_0_phi);
	free(neighbor_basis_0_dphi);

	free(neighbor_basis_1_phi);
	free(neighbor_basis_1_dphi);

	free(neighbor_basis_3_phi);
	free(neighbor_basis_3_dphi);
      }
  }// face q loop


  // printf("ok\n");
  if(counter == 1)
    block_LocaltoGlobal_neighbor(dof_on_elm,dof_on_elm,FE,A,ALoc);
  else if(counter == 2)
    {
      //printf("=============== 111111111111 =================== \n");
      block_LocaltoGlobal_neighbor(dof_on_elm,dof_on_elm,FE,A,ALoc_u_v);
      // row is neighbor

      //printf("=============== 222222222222 =================== \n");
      block_LocaltoGlobal_neighbor(dof_on_elm_neighbor,dof_on_elm,FE,A,ALoc_u_vneighbor);
      //block_LocaltoGlobal_neighbor(dof_on_elm,dof_on_elm_neighbor,FE,A,ALoc_u_vneighbor);

      //printf("=============== 333333333333 =================== \n");
      block_LocaltoGlobal_neighbor(dof_on_elm,dof_on_elm_neighbor,FE,A,ALoc_uneighbor_v);
      //block_LocaltoGlobal_neighbor(dof_on_elm_neighbor,dof_on_elm,FE,A,ALoc_uneighbor_v);

      //printf("=============== 444444444444 =================== \n");
      block_LocaltoGlobal_neighbor(dof_on_elm_neighbor,dof_on_elm_neighbor,FE,A,ALoc_uneighbor_vneighbor);
    }
  else{
    printf("!@#$\n");exit(0);
  }





  //printf("Local Assemble Face - end\n");
}

/******************************************************************************************************/
/*!
 * \fn void FEM_Block_RHS_Local(REAL* bLoc,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
 *
 * \brief Computes the local Right hand side vector for a block FEM system
 *        b_i  = <f,phi_i>
 *
 * \param FE            Block FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param dof_on_elm    Specific DOF on element
 * \param v_on_elm      Specific Vertices on element
 * \param elm           Current element
 * \param rhs           Function that gives RHS (in FEM block ordering
 * \param time          Physical Time if time dependent
 *
 \return bLoc         Local RHS Vector
 *
 *
 */

// ISSUE : this doesn't allow to get the exact solution out...
void FEM_Block_RHS_Local_Elasticity(dvector *b,REAL* bLoc, REAL *solution, \
				    block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,
				    void (*rhs)(REAL *,REAL *,REAL,void *),REAL time, REAL timestep, \
				    void (*truesol)(REAL *,REAL *,REAL,void *),
				    void (*D_truesol)(REAL *,REAL *,REAL,void *),
				    void (*truesol_dt)(REAL *,REAL *,REAL,void *),
				    REAL *err, INT* switch_on_face)
{

  //printf("RHS ASSEMBLE  \n");
  //fflush(stdout);

  // Loop Indices
  INT i,quad,test;
  REAL C_0 = 0.;
  bool bEG = BOOL_EG_MECHANICS;
  bool bEG_Pressure = BOOL_EG_PRESSURE;

  // Mesh and FE data
  INT dim = mesh->dim;
  INT dof_per_elm = 0;
  INT nun=0;

  REAL alpha = 1.;

  for(i=0;i<FE->nspaces;i++) {
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
    if(FE->var_spaces[i]->FEtype<20) /* Scalar Element */
      nun++;
    else /* Vector Element */
      nun += dim;
  }

  INT* local_dof_on_elm = NULL;
  INT local_row_index=0;
  INT unknown_index=0;

  // Quadrature Weights and Nodes
  REAL w;
  INT maxdim=nun;
  REAL qx[maxdim];

  // Right-hand side function at Quadrature Nodes
  REAL rhs_val[nun];

  // For initial condition
  REAL* val_true = (REAL *) calloc(nun,sizeof(REAL));
  REAL* val_true_D = (REAL *) calloc(nun*2,sizeof(REAL));

  // for the preivous time step value
  //REAL* val_sol = (REAL *) calloc(nun,sizeof(REAL));


  coordinates *barycenter = allocatecoords(1,dim);
  barycenter->x[0] = mesh->el_mid[elm*dim];
  barycenter->y[0] = mesh->el_mid[elm*dim + 1];


  //  Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];

    (*rhs)(rhs_val,qx,time,&(mesh->el_flag[elm]));

    //if(time == timestep){

    //printf("********\n");
    //printf("********\n");
    //printf("***** time = %f ***\n", time);
    (*truesol)(val_true,qx,time - timestep,&(mesh->el_flag[elm]));
    (*D_truesol)(val_true_D,qx,time - timestep,&(mesh->el_flag[elm]));


    //if(time != timestep){
    //printf("qx[0] = %f, qx[1] = %f, time = %f, val_true_D[0] = %f, val_true_D[3] = %f \n",
    //	   qx[0], qx[1], time-timestep, val_true_D[0], val_true_D[3]);
    //exit(0);
    //}

    local_row_index=0;
    unknown_index=0;
    local_dof_on_elm=dof_on_elm;

    REAL* u_comp = solution;
    REAL u0_value_at_q = 0.;
    REAL u1_value_at_q = 0.;
    REAL val[dim];
    val[0] = 0.0;
    val[1] = 0.0;

    REAL p_value_at_q = 0.;
    REAL p_eg_value_at_q = 0.;
    int j,dof;

    REAL grad_val0[dim];
    grad_val0[0]=0;
    grad_val0[1]=0;
    REAL grad_val1[dim];
    grad_val1[0]=0;
    grad_val1[1]=0;

    // U  bEG
    REAL grad_val2[dim];
    grad_val2[0] = 0.0;
    grad_val2[1] = 0.0;

    REAL grad_val3[dim];
    grad_val3[0] = 0.0;
    grad_val3[1] = 0.0;

    // Pressure
    REAL grad_val4[dim];
    grad_val4[0]=0;
    grad_val4[1]=0;


    //printf("time = %f, timestep = %f, val_true[3] = %f, val_true[4] = %f \n", time, timestep, val_true[3], val_true[4]);

    for(i=0;i<FE->nspaces;i++) {

      // Basis Functions and its derivatives if necessary
      get_FEM_basis(FE->var_spaces[i]->phi,
		    FE->var_spaces[i]->dphi,
		    qx,
		    v_on_elm,
		    local_dof_on_elm,
		    mesh,
		    FE->var_spaces[i]);




      for (test=0; test<FE->var_spaces[i]->dof_per_elm;test++) {

	//SLEE
	if(i == 0 || i == 1)
	  {
	    //This is for  u0 and u1 (CG part)
	    //printf(" i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
	    bLoc[(local_row_index+test)] += w*rhs_val[i]*FE->var_spaces[i]->phi[test];
	  }
	else if(i == 2)
	  {
	    //SLEE
	    // This is for  u2:: (EG part)
	    //Note that new DOF is \phi^3 = [x ; y]
	    //printf(" i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
	    bLoc[(local_row_index+test)] += w*(rhs_val[0]*  (qx[0]-barycenter->x[0])  +rhs_val[1]* (qx[1]-barycenter->y[0]));
	  }
	else if(i == 3)
	  {

	    bLoc[(local_row_index+test)] += w*rhs_val[3]*FE->var_spaces[3]->phi[test];

	    //printf(" FE->var_spaces[3]->phi[test] = %f\n", FE->var_spaces[3]->phi[test]);


	  }



      }

      local_row_index  += FE->var_spaces[i]->dof_per_elm;



    }
  }// end quad


  INT v_per_elm = mesh->v_per_elm;
  //Sanity Check // SLEE
  //for(i=0;i<FE->nspaces;i++) {
  //printf("dof_per_face = %d,   dof_per_face_blk[%d] = %d \n", dof_per_face, i, dof_per_face_blk[i]);
  //}
  //exit(0);
  REAL *data_face=calloc(dim+dim*dim, sizeof(REAL)); //coords of normal and vertices of a face.
  REAL *xfi=data_face; //coords of vertices on face i.
  REAL *finrm=xfi+dim*dim; //coords of normal vector on face i

  //REAL *data_face_end=finrm + dim; //end
  INT nq1d_face=3;
  qcoordinates *cq_face = allocateqcoords_bdry(nq1d_face,1,dim,2);
  //REAL fiarea=0e0;
  INT jk,k,face,quad_face,rowa, rowb, jcntr,ed, jkl, j;
  REAL qx_face[maxdim];
  REAL w_face;
  INT nquad_face= cq_face->nq_per_elm; //quad<cq_face->nq_per_elm;

  // SLEE

  bool bWeakBC_RHS = BOOL_WEAKLY_IMPOSED_BC;

  if(bWeakBC_RHS){

    //printf("=============================\n");
    //printf("** ELEMENT = %d \n", elm );
    //printf("=============================\n");
    for(jk=mesh->el_f->IA[elm];jk<mesh->el_f->IA[elm+1];jk++){

      //printf("jk = %d, mesh->el_f->IA[element] = %d, mesh->el_f->IA[element+1] = %d  \n",
      //     jk, mesh->el_f->IA[elm], mesh->el_f->IA[elm+1]);
      //  j = face0, face1, face2
      j=jk - mesh->el_f->IA[elm];
      // face is the global number
      face=mesh->el_f->JA[jk];

      // Get normal vector values.
      finrm[0]=mesh->f_norm[face*dim+0];
      if(dim>1)
	finrm[1]=mesh->f_norm[face*dim+1];
      if(dim>2)
	finrm[2]=mesh->f_norm[face*dim+2];

      for(jkl=mesh->f_v->IA[face];jkl<mesh->f_v->IA[face+1];jkl++){

	//printf("** jkl = %d, mesh->f_v->IA[face] = %d, mesh->f_v->IA[face+1] = %d \n", jkl, mesh->f_v->IA[face], mesh->f_v->IA[face+1]);

	j=jkl-mesh->f_v->IA[face];
	k=mesh->f_v->JA[jkl];

	//printf("** j = %d, k = %d \n", j, k );

	xfi[j*dim+0]=mesh->cv->x[k];
	if(dim>1)
	  xfi[j*dim+1]=mesh->cv->y[k];
	if(dim>2)
	  xfi[j*dim+2]=mesh->cv->z[k];

	//printf("** xfi[j*dim+0] = %f,  xfi[j*dim+1] = %f \n",  xfi[j*dim+0] ,  xfi[j*dim+1]);

      }

      // Get the BD values
      REAL* val_true_face = (REAL *) calloc(nun,sizeof(REAL));
      REAL* val_true_face_n = (REAL *) calloc(nun,sizeof(REAL));
      REAL* val_true_face_n_neighbor = (REAL *) calloc(nun,sizeof(REAL));
      REAL* val_true_dt_face = (REAL *) calloc(nun,sizeof(REAL));


      // FOR NEIGHBOR..
      INT* v_on_elm_neighbor = (INT *) calloc(v_per_elm,sizeof(INT));
      INT* dof_on_elm_neighbor = (INT *) calloc(dof_per_elm,sizeof(INT));
      INT *local_dof_on_elm_face_interface = NULL;
      INT *local_dof_on_elm_neighbor = NULL;


      // NOW FOR FACES (at BOUNDARY)
      INT* local_dof_on_elm_face = NULL;
      //Neighbor
      INT neighbor_index[2];
      neighbor_index[0] = -1;
      neighbor_index[1] = -1;
      INT counter = 0;

      iCSRmat *f_el=NULL;
      f_el=(iCSRmat *)malloc(1*sizeof(iCSRmat)); // face_to_element;
      icsr_trans(mesh->el_f,f_el); // f_el=transpose(el_f);

      // printf("=============================\n");
      //printf("** ELM  = %d   FACE = %d \n",  elm, face );
      //printf("=============================\n");
      int pq,nbr0;
      for(pq=f_el->IA[face];pq<f_el->IA[face+1];pq++){

	//printf("-- pq = %d, f_el->IA[face] = %d, f_el->IA[face+1] = %d \n", pq, f_el->IA[face], f_el->IA[face+1]);
	nbr0=f_el->JA[pq];
	//printf("-- nbr0 = %d  \n", nbr0);

	neighbor_index[counter] = nbr0;
	counter++;

      }

      //Sanity Check
      /* print out
	 for(pq=0;pq<2;++pq)
	 {
	 if(counter == 2)
	 {
	 //printf("neighbor_index[%d]= %d || counter  = %d\n", pq, neighbor_index[pq],counter);
	 }
	 else if(counter == 1){
	 if(pq == 0)
	 {
	 //printf("neighbor_index[%d]= %d || counter  = %d\n", pq, neighbor_index[pq],counter);
	 }
	 else{
	 //printf("-\n");
	 }
	 }
	 }
      */

      double fiarea=mesh->f_area[face];
      double lambda = LAME_LAMBDA_GLOBAL;//000000.;//000000.;//000000.; //000000.0;
      double penalty_term = PENALTY_PARAMETER_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));
      double BC_penalty_term = BC_PENALTY_PARAMETER_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));

      //penalty_term*=lambda;

      double penalty_term_pressure =  PENALTY_PARAMETER_PRESSURE_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));
      double BC_penalty_term_pressure =  BC_PENALTY_PARAMETER_PRESSURE_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));

      //nq1d_face == 3, 3 quad points.
      zquad_face(cq_face,nq1d_face,dim,xfi,fiarea);

      REAL edge_length = mesh->ed_len[face];

      if(mesh->f_flag[face]>0) {


	for (quad_face=0;quad_face<nquad_face;quad_face++) {

	  qx_face[0] = cq_face->x[quad_face];
	  qx_face[1] = cq_face->y[quad_face];

	  if(dim==3) qx_face[2] = cq_face->z[quad_face];

	  w_face = cq_face->w[quad_face];

	  (*exact_sol2D)(val_true_face,qx_face,time,
			 &(mesh->el_flag[elm]));  // ???


	  (*exact_sol2D)(val_true_face_n,qx_face,time-timestep,
			 &(mesh->el_flag[elm]));  // ???


	  //true solution is in n+1 time
	  (*truesol_dt)(val_true_dt_face,
			qx_face,
			time,
			&(mesh->el_flag[elm]));

	  local_row_index=0;
	  unknown_index=0;
	  local_dof_on_elm_face = dof_on_elm;

	  for(i=0;i<FE->nspaces;i++) {

	    get_FEM_basis(FE->var_spaces[i]->phi,FE->var_spaces[i]->dphi,
			  qx_face,
			  v_on_elm,
			  local_dof_on_elm_face,
			  mesh,FE->var_spaces[i]);

	    for (test=0; test<FE->var_spaces[i]->dof_per_elm;test++) {
	      //SLEE
	      if(i==0 || i ==1)
		{
		  bLoc[(local_row_index+test)] += BC_penalty_term * w_face*
		    (val_true_face[unknown_index]*  FE->var_spaces[i]->phi[test]);

		  //DEBUG100
		  //bLoc[(local_row_index+test)] += penalty_term * w_face*
		  //(val_true_face[3]*  FE->var_spaces[i]->phi[test] *finrm[i]);

		}

	      else if(i == 2)
		{
		  //SLEE
		  // Note that new DOF is \phi^3 = [x ; y]
		  //printf("i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
		  bLoc[(local_row_index+test)] +=
		    BC_penalty_term *  w_face * (val_true_face[0] * (qx_face[0] - barycenter->x[0])
						 + val_true_face[1]*  (qx_face[1] - barycenter->y[0]));

		}
	      else if(i == 3)
		{

		  //Stokes 3
		  bLoc[(local_row_index+test)] +=
		    w_face * (val_true_face[0] * finrm[0] + val_true_face[1] * finrm[1])
		    * FE->var_spaces[3]->phi[test];

		}


	      //sanity check
	      if(unknown_index != i)
		{
		  printf("unknown index != i \n");
		  exit(0);
		}
	      //SLEE

	    }
	    unknown_index++;

	    local_dof_on_elm_face += FE->var_spaces[i]->dof_per_elm;
	    local_row_index += FE->var_spaces[i]->dof_per_elm;

	  } // i = 0

	} //face

      }// else if


      //printf("FEREE \n");
      icsr_free(f_el);
      //printf("FEREE DONE\n");
    }// for each face


  }


  block_LocaltoGlobal_RHS(dof_on_elm,FE,b,bLoc);

  //if(val_true) free(val_true);
  //if(val_sol) free(val_sol);

  //printf("RHS ASSEMBLE END \n");
  return;
}





void local_assembly_Elasticity(block_dCSRmat* A,dvector *b,REAL* ALoc, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time, REAL timestep,  INT* switch_on_face)
//void local_assembly_Elasticity(REAL* ALoc, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)
{

  //printf("Local Assembly -- start\n");

  bool bEG = BOOL_EG_MECHANICS;
  bool bEG_Pressure = BOOL_EG_PRESSURE;

  REAL alpha = 1.;

  REAL C_0 = 0.001;
  // Loop indices
  INT i,j,idim,quad,test,trial;
  INT dim = mesh->dim;
  // Mesh and FE data
  INT dof_per_elm = 0;
  for (i=0; i<FE->nspaces;i++)
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;

  //printf("CHECK POINT: FE->nspaces = %d, dof_per_elem = %d, FE->var_spaces[0]->dof_per_elm = %d,  FE->var_spaces[1]->dof_per_elm = %d ,  FE->var_spaces[1]->dof_per_elm = %d \n",
  //	 FE->nspaces, dof_per_elm, FE->var_spaces[0]->dof_per_elm,FE->var_spaces[1]->dof_per_elm, FE->var_spaces[2]->dof_per_elm);
  //exit(0);

  INT *local_dof_on_elm = NULL;
  INT *local_dof_on_elm_face = NULL;
  INT *local_dof_on_elm_face_interface = NULL;
  INT *local_dof_on_elm_neighbor = NULL;
  //printf("dof per elm = %d" , dof_per_elm);
  INT* dof_on_elm_neighbor = (INT *) calloc(dof_per_elm,sizeof(INT));

  //SLEE
  coordinates *barycenter = allocatecoords(1,dim);
  coordinates *barycenter_neighbor = allocatecoords(1,dim);
  barycenter->x[0] = mesh->el_mid[elm*dim];
  barycenter->y[0] = mesh->el_mid[elm*dim + 1];

  //printf("ELEMENT = %D, x = %f , y = %f \n", elm,  barycenter->x[0],   barycenter->y[0]);
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* ALoc_neighbor = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* bLoc=NULL;

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim];

  // Stiffness Matrix Entry
  REAL kij = 0.0;

  // Keep track of local indexing
  INT local_row_index, local_col_index;

  double lambda = LAME_LAMBDA_GLOBAL ;

  // Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {

    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];

    if(mesh->dim==3)
      qx[2] = cq->z[elm*cq->nq_per_elm+quad];

    w = cq->w[elm*cq->nq_per_elm+quad];

    //  Get the Basis Functions at each quadrature node
    // u = (u1,u2,u3) and v = (v1,v2,v3)
    get_FEM_basis(FE->var_spaces[0]->phi,
		  FE->var_spaces[0]->dphi,
		  qx,
		  v_on_elm,
		  dof_on_elm,
		  mesh,
		  FE->var_spaces[0]);
    /////////

    // u2
    local_dof_on_elm = dof_on_elm + FE->var_spaces[0]->dof_per_elm;

    get_FEM_basis(FE->var_spaces[1]->phi,
		  FE->var_spaces[1]->dphi,
		  qx,
		  v_on_elm,
		  local_dof_on_elm,
		  mesh,
		  FE->var_spaces[1]);

    // u_eg
    local_dof_on_elm += FE->var_spaces[1]->dof_per_elm;
    get_FEM_basis(FE->var_spaces[2]->phi,
		  FE->var_spaces[2]->dphi,
		  qx,
		  v_on_elm,
		  local_dof_on_elm,
		  mesh,
		  FE->var_spaces[2]);


    // p
    local_dof_on_elm += FE->var_spaces[2]->dof_per_elm;

    get_FEM_basis(FE->var_spaces[3]->phi,
		  FE->var_spaces[3]->dphi,
		  qx,
		  v_on_elm,
		  local_dof_on_elm,
		  mesh,
		  FE->var_spaces[3]);

    // p_eg
    local_dof_on_elm += FE->var_spaces[3]->dof_per_elm;



    //------#1---------------------------------------------------- START//
    // u1-v1 block: 2*<dx(u1),dx(v1)> + <dy(u1),dy(v1)>
    local_row_index = 0;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	//printf("quad = %d ::   local_row_index = %d,  qx[0] = %f, qx[1] = %f,  test = %d \n",
	//     quad,  local_row_index,  qx[0], qx[1] ,   test);

	kij = 2.0*FE->var_spaces[0]->dphi[test*dim]*FE->var_spaces[0]->dphi[trial*dim];
	kij += 1.0*FE->var_spaces[0]->dphi[test*dim+1]*FE->var_spaces[0]->dphi[trial*dim+1];

	// NEW SLEE : elasticity div part
	// u1-v1 block : <dx(u1),dx(v1)>
	kij += lambda*FE->var_spaces[0]->dphi[test*dim]*FE->var_spaces[0]->dphi[trial*dim];

        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // u1-v2 block <dy(u1),dx(v2)>
    local_row_index = FE->var_spaces[0]->dof_per_elm;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
        kij = 1.0*FE->var_spaces[1]->dphi[test*dim+0]*FE->var_spaces[0]->dphi[trial*dim+1];

	// NEW SLEE : elasticity div part
	// u1-v2 block : <dx(u1),dx(v2)>
	kij += lambda*FE->var_spaces[1]->dphi[test*dim+1]*FE->var_spaces[0]->dphi[trial*dim];

        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // u2-v1 block : <dx(u2),dy(v1)>>
    local_row_index = 0;
    local_col_index = FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
        kij = 1.0*FE->var_spaces[0]->dphi[test*dim+1]*FE->var_spaces[1]->dphi[trial*dim+0];

	// NEW SLEE : elasticity div part
	// u2-v1 block : <dy(u2),dx(v1)>
	kij += lambda*FE->var_spaces[0]->dphi[test*dim]*FE->var_spaces[1]->dphi[trial*dim+1];

        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // u2-v2 block <dx(u2),dx(v2)> + 2*<dy(u2),dy(v2)>
    local_row_index = FE->var_spaces[0]->dof_per_elm;
    local_col_index = FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	kij = 2.0*FE->var_spaces[1]->dphi[test*dim+1]*FE->var_spaces[1]->dphi[trial*dim+1];
	kij += 1.0*FE->var_spaces[1]->dphi[test*dim]*FE->var_spaces[1]->dphi[trial*dim];

	// NEW SLEE : elasticity div part
	// u2-v2 block : <dy(u2),dx(v2)>
   	kij += lambda*FE->var_spaces[1]->dphi[test*dim+1]*FE->var_spaces[1]->dphi[trial*dim+1];

        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    //q-q block:
    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
    local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
    if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){

	kij = 2.0 * 1. * 1.;
	kij += 2.0 * 1. * 1.;

	// NEW SLEE : elasticity div part
	// u2-v1 block : <dy(u2),dx(v1)>
	kij += lambda * 1. * 1.;
	kij += lambda * 1. * 1.;

	kij += lambda * 1. * 1.;
	kij += lambda * 1. * 1.;

	ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    //Local Matrix
    if(bEG){
      //EG PART
      // u0- q block < 2mu e(u1) * e(q) > + lambda <div u1 : e(q) >
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
      local_col_index = 0;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	  //u0 - q0
	  kij  = 2.0*FE->var_spaces[0]->dphi[trial*dim] * 1.; //FE->var_spaces[dim]->dphi[test*dim];
	  //u0 - q1
	  kij += 0.;

	  // Divergence
	  // u0-q0 block : <dx(u1),dx(v1)>
	  kij += lambda*FE->var_spaces[0]->dphi[trial*dim]* 1.;//FE->var_spaces[dim]->dphi[test*dim];

	  // u0-q1
	  kij += lambda*FE->var_spaces[0]->dphi[trial*dim]* 1.;

	  ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
	}
      }

      // u1-q block
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	  //u1 - q0
	  // = 0
	  //u1 - q1
	  kij = 2.0*FE->var_spaces[1]->dphi[trial*dim+1] * 1.;

	  // NEW SLEE : elasticity div part
	  // u2-v1 block : <dy(u2),dx(v1)>
	  kij += lambda*FE->var_spaces[1]->dphi[trial*dim+1] * 1.;
	  kij += lambda*FE->var_spaces[1]->dphi[trial*dim+1] * 1.;

	  ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
	}
      }


      // q-v0 block:
      local_row_index = 0;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){

	  //q0 - v0
	  kij = 2.0* 1. * FE->var_spaces[0]->dphi[test*dim]; //FE->var_spaces[dim]->dphi[test*dim];
	  //q1 - v0
	  // = 0.;

	  // Divergence
	  // u1-q block : <dx(u1),dx(v1)>
	  kij += lambda * 1. * FE->var_spaces[0]->dphi[test*dim];//FE->var_spaces[dim]->dphi[test*dim];
	  kij += lambda * 1. * FE->var_spaces[0]->dphi[test*dim];

	  ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
	}
      }


      // q-v1 block:
      local_row_index = FE->var_spaces[0]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){

	  //q0 - v1
	  //q1 - v1
	  kij = 2.0* 1. * FE->var_spaces[1]->dphi[test*dim+1];

	  // NEW SLEE : elasticity div part
	  // u2-v1 block : <dy(u2),dx(v1)>
	  kij += lambda*FE->var_spaces[1]->dphi[test*dim+1] * 1.;
	  kij += lambda*FE->var_spaces[1]->dphi[test*dim+1] * 1.;

	  ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
	}
      }



    } // bEG
    //------#1---------------------------------------------------- END//


    //------PORO #0---------------------------------------------------- START//
    //////////////////////////////////////////////////////////////////////////
    // porop0
    // PoroElasticity   -alpha(p, div v)
    // Add Pressure to Mechanics
    //////////////////////////////////////////////////////////////////////////
    //   -alpha(p, div v)
    //          CG  CG
    //          trial (column) -3  test (row) -0
    // test - u0 cg
    local_row_index = 0.;
    // trial - p
    local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	kij = -alpha * FE->var_spaces[3]->phi[trial] *  (FE->var_spaces[0]->dphi[test*dim]);

	ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // test - u1 cg
    local_row_index = FE->var_spaces[0]->dof_per_elm;
    // trial - p
    local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	kij = -alpha * FE->var_spaces[3]->phi[trial] *  (FE->var_spaces[1]->dphi[test*dim+1]);

	ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    if(bEG){
      //   -alpha(p, div v)
      //          CG  EG
      //        trial (column) 3  test (row) 2
      // test - u eg
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      // trial - p
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	  kij = -alpha * FE->var_spaces[3]->phi[trial] * 2.;

	  ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
	}
      }
    }


    //------PORO #0---------------------------------------------------- END//


    //------STOKES  #0---------------------------------------------------- START//
    //////////////////////////////////////////////////////////////////////////
    // Transpose of Stokes #0
    // PoroElasticity   -(div v, p)
    // Add Pressure to Mechanics
    //////////////////////////////////////////////////////////////////////////
    //   -alpha(div v, p)
    //          CG  CG
    //          trial (column) -3  test (row) -0
    // test - p cg
    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
    // trial - u0 cg
    local_col_index = 0.;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	kij = -alpha * FE->var_spaces[3]->phi[test] *  (FE->var_spaces[0]->dphi[trial*dim]);

	ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // test - p cg
    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
    // trial - u1 cg
    local_col_index = FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	kij = -alpha * FE->var_spaces[3]->phi[test] *  (FE->var_spaces[1]->dphi[trial*dim+1]);

	ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    if(bEG){
      //   -alpha(div v, p)
      //          EG  CG
      //        test (column) 3  trial (row) 2
      // test - p cg
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
      // trial - u eg
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){

	  kij = -alpha * FE->var_spaces[3]->phi[test] * 2.;

	  ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
	}
      }

    }

    //------STOKES #0---------------------------------------------------- END//



    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    // PRESSURE BLOCK
    // p-p block:
    /*
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

      kij =  FE->var_spaces[3]->phi[test] *  FE->var_spaces[3]->phi[trial];

      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
      }
    */



  }//QUAD


  block_LocaltoGlobal_neighbor(dof_on_elm,dof_on_elm,FE,A,ALoc);

  //printf("*********** {DONE} @ END  ************************ \n");


  //printf("Local Assembly --end\n");

  return;
}
/*********************************************************************/
/****** MAIN DRIVER **************************************************************/
INT main(int argc, char* argv[])
{

  printf("\n===========================================================================\n");
  printf("Beginning Program to solve Elasticity Equation.\n");
  printf("===========================================================================\n");
  // Aug.3.2020 SLEE
  // Define variables forthe error convergence test
  int total_num_cycle = TOTAL_NUM_CYCLES_GLOBAL ;
  // SLEE initialize the vectors to save the errors for each cycle
  double L2_error_per_cycle[total_num_cycle];
  double H1_error_per_cycle[total_num_cycle];
  double H1_stress_error_per_cycle[total_num_cycle];

  double L2_EG_error_per_cycle[total_num_cycle];
  double H1_EG_error_per_cycle[total_num_cycle];
  double H1_stress_EG_error_per_cycle[total_num_cycle];
  double H1_energy_EG_error_per_cycle[total_num_cycle];


  // For Pressure
  double L2_error_p_per_cycle[total_num_cycle];
  double H1_error_p_per_cycle[total_num_cycle];

  //double L2_p_EG_error_per_cycle[total_num_cycle];
  //double H1_p_EG_error_per_cycle[total_num_cycle];


  // SLEE vector to save the DOF
  int dof_per_cycle_CG[total_num_cycle];
  int dof_per_cycle_EG[total_num_cycle];

  // For Pressure
  int dof_per_cycle_CG_p[total_num_cycle];
  int dof_per_cycle_EG_p[total_num_cycle];

  double mesh_size_per_cycle[total_num_cycle];

  // SLEE vector to save the convergence rate
  double L2_conv_rate_per_cycle[total_num_cycle];
  double H1_conv_rate_per_cycle[total_num_cycle];
  double H1_stress_conv_rate_per_cycle[total_num_cycle];


  double L2_EG_conv_rate_per_cycle[total_num_cycle];
  double H1_EG_conv_rate_per_cycle[total_num_cycle];
  double H1_stress_EG_conv_rate_per_cycle[total_num_cycle];
  double H1_energy_EG_conv_rate_per_cycle[total_num_cycle];

  // For Pressure

  double L2_p_conv_rate_per_cycle[total_num_cycle];
  double H1_p_conv_rate_per_cycle[total_num_cycle];

  //double L2_p_EG_conv_rate_per_cycle[total_num_cycle];
  //double H1_p_EG_conv_rate_per_cycle[total_num_cycle];


  int global_dim_space = 0;




  // ALL TIME STEPPING ALGORITHMS /////
  REAL time = 0.;
  REAL timestep = 0.01;
  INT timestep_number = 0;
  INT total_timestep = 10;

  for(int cycle=0; cycle<total_num_cycle; ++cycle) {
      //Aug.3.2020 SLEE initilize
      L2_error_per_cycle[cycle] = 0.;
      H1_error_per_cycle[cycle] = 0.;
      H1_stress_error_per_cycle[cycle] = 0.;

      L2_EG_error_per_cycle[cycle] = 0.;
      H1_EG_error_per_cycle[cycle] = 0.;
      H1_stress_EG_error_per_cycle[cycle] = 0.;
      H1_energy_EG_error_per_cycle[cycle] = 0.;

      // For Pressure
      L2_error_p_per_cycle[cycle] = 0.;
      H1_error_p_per_cycle[cycle] = 0.;


      dof_per_cycle_CG[cycle]=0;
      dof_per_cycle_EG[cycle]=0;

      // For Pressure
      dof_per_cycle_CG_p[cycle]=0;
      dof_per_cycle_EG_p[cycle]=0;

      mesh_size_per_cycle[cycle]=0.;

      L2_conv_rate_per_cycle[cycle]=0.;
      H1_conv_rate_per_cycle[cycle]=0.;
      H1_stress_conv_rate_per_cycle[cycle]=0.;

      L2_EG_conv_rate_per_cycle[cycle]=0.;
      H1_EG_conv_rate_per_cycle[cycle]=0.;
      H1_stress_EG_conv_rate_per_cycle[cycle]=0.;
      H1_energy_EG_conv_rate_per_cycle[cycle]=0.;

      // For Pressure
      L2_p_conv_rate_per_cycle[cycle]=0.;
      H1_p_conv_rate_per_cycle[cycle]=0.;


      printf("************ CYCLE   %d  /   %d  ************** \n", cycle, total_num_cycle);

      /****** INITIALIZE PARAMETERS **************************************************/
      // loops
      INT i;
      // Overall CPU Timing
      clock_t clk_overall_start = clock();

      // Set Parameters from Reading in Input File
      input_param inparam;
      param_input_init(&inparam);
      param_input("./input.dat", &inparam);

      // Open gridfile for reading
      printf("\nCreating mesh and FEM spaces:\n");
      //FILE* gfid = HAZ_fopen(inparam.gridfile,"r");
      //SLEE
      FILE* gfid;

      //Jul.10.2020 SLEE: setup the code to read the different mesh files for each cycle
      //gfid = HAZ_fopen(inparam.gridfile,"r");
      char filename_per_cycle[512]={'\0'};
      //sprintf(filename_per_cycle, "%s%d.haz", inparam.gridfile,cycle);

      //DEBUG SIMPLE MESH
      sprintf(filename_per_cycle, "%s%d.haz", inparam.gridfile,cycle+1);
      gfid = HAZ_fopen(filename_per_cycle,"r");
      if(gfid == NULL){
	perror("Could not find and open the file !!!! ");
      }

      // Create the mesh (now we assume triangles in 2D or tetrahedra in 3D)
      // File types possible are 0 - old format; 1 - vtk format
      INT mesh_type = 0;
      clock_t clk_mesh_start = clock(); // Time mesh generation FE setup
      mesh_struct mesh;

      //printf(" --> loading grid from file: %s\n",inparam.gridfile);
      //Jul.10. 2020 SLEE
      printf(" --> loading grid from file: %s\n",filename_per_cycle);

      initialize_mesh(&mesh);   // Q1. Why only here?
      creategrid_fread(gfid,mesh_type,&mesh);
      fclose(gfid);

      INT dim = mesh.dim;
      // Jul.10.2020 SLEE : for further use in the convergence test
      global_dim_space = dim;

      // Get Quadrature Nodes for the Mesh
      INT nq1d = inparam.nquad; /* Quadrature points per dimension */
      qcoordinates *cq = get_quadrature(&mesh,nq1d);

      // Get info for and create FEM spaces
      // Order of elements: 0 - P0; 1 - P1; 2 - P2; 20 - Nedlec; 30 - Raviart-Thomas
      INT order_u = 1;
      INT order_u_eg = 0;

      // Need Spaces for each component of the Mechanics plus pressure
      fespace FE_ux; // Mechanics in x direction
      create_fespace(&FE_ux,&mesh,order_u);
      fespace FE_uy; // Mechanics in y direction
      create_fespace(&FE_uy,&mesh,order_u);
      fespace FE_uz; // Mechanics in z direction
      if(dim==3) create_fespace(&FE_uz,&mesh,order_u);
      fespace FE_u_eg; // Mechanics EG
      create_fespace(&FE_u_eg,&mesh,order_u_eg);

      INT order_p = 0;
      //INT order_p_eg = 0;
      fespace FE_p; // Pressuer
      create_fespace(&FE_p,&mesh,order_p);
      //fespace FE_p_eg; // Pressuer
      //create_fespace(&FE_p_eg,&mesh,order_p_eg);
      // Set Dirichlet Boundaries

      set_dirichlet_bdry(&FE_ux,&mesh,1,1);
      set_dirichlet_bdry(&FE_uy,&mesh,1,1);
      //if(dim==3) set_dirichlet_bdry(&FE_uz,&mesh,1,1);
      set_dirichlet_bdry(&FE_u_eg,&mesh,1,1);

      set_dirichlet_bdry(&FE_p,&mesh,1,1);
      //set_dirichlet_bdry(&FE_p_eg,&mesh,1,1);

      if(BOOL_WEAKLY_IMPOSED_BC){
        for(i=0;i<FE_u_eg.ndof;i++) {
          FE_u_eg.dirichlet[i] = 0;
        }
        for(i=0;i<FE_ux.ndof;i++) {
          FE_ux.dirichlet[i] = 0;
        }
        for(i=0;i<FE_uy.ndof;i++) {
          FE_uy.dirichlet[i] = 0;
        }

        for(i=0;i<FE_p.ndof;i++) {
          FE_p.dirichlet[i] = 0;
        }
      }

      // Create Block System with ordering (u,p)
      INT u_ndof = FE_ux.ndof + FE_uy.ndof + FE_u_eg.ndof;
      INT p_ndof = FE_p.ndof;

      //p_debug
      INT ndof = u_ndof + p_ndof;
      if(dim==3) ndof += FE_uz.ndof;

      // Get Global FE Space
      block_fespace FE;
      FE.nun = dim+1 +1;  // p_debug
      FE.ndof = ndof;
      FE.nbdof = FE_ux.nbdof + FE_uy.nbdof + FE_u_eg.nbdof +FE_p.nbdof; //+FE_p_eg.nbdof;
      //if(dim==3) FE.nbdof += FE_uz.nbdof;
      FE.nspaces = dim+1 +1; // SLEE?
      FE.var_spaces = (fespace **) calloc(dim+1 +1,sizeof(fespace *));

      /*
	FE.nun = dim+1;  // p_debug
	FE.ndof = ndof;
	FE.nbdof = FE_ux.nbdof + FE_uy.nbdof + FE_u_eg.nbdof;// +FE_p.nbdof+FE_p_eg.nbdof;
	FE.nspaces = dim+1;
	FE.var_spaces = (fespace **) calloc(dim+1,sizeof(fespace *));
      */

      FE.var_spaces[0] = &FE_ux;
      FE.var_spaces[1] = &FE_uy;
      if(dim==3) FE.var_spaces[2] = &FE_uz;
      FE.var_spaces[dim] = &FE_u_eg;

      //p_debug
      FE.var_spaces[dim+1]   = &FE_p;
      //FE.var_spaces[dim+1+1] = &FE_p_eg;

      // Set Dirichlet Boundaries
      if(!BOOL_WEAKLY_IMPOSED_BC)
	set_dirichlet_bdry_block(&FE,&mesh);


      clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
      printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
	     (REAL) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
      /*******************************************************************************/

      printf("***********************************************************************************\n");
      printf("Number of Elements = %d\tOrder of Quadrature = %d\n",mesh.nelm,2*nq1d-1);
      printf("\n\t--- Element Type ---\n");
      printf("Mechanics Element Type = %d\t Meachnics EG  Type = %d\n",order_u,order_u_eg);
      printf("\n\t--- Degrees of Freedom ---\n");
      printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nv,mesh.nedge,mesh.nface);
      printf("\t--> DOF: %d\n",FE.ndof);
      printf("\n\t--- Boundaries ---\n");
      printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nbv,mesh.nbedge,mesh.nbface);
      printf("\t--> Boundary DOF: %d\n",FE.nbdof);
      printf("***********************************************************************************\n\n");

      //Jul.10.2020 SLEE: insert the total number of DOF for convergnece computation
      dof_per_cycle_CG[cycle]  = FE_ux.ndof + FE_uy.ndof;////FE.ndof + FE.nbdof;
      dof_per_cycle_EG[cycle] =  FE_ux.ndof + FE_uy.ndof + FE_u_eg.ndof;////FE.ndof + FE.nbdof;

      dof_per_cycle_CG_p[cycle]  = FE_p.ndof;
      dof_per_cycle_EG_p[cycle] =  FE_p.ndof;// + FE_p_eg.ndof;


      printf("FE.ndof = %d | u_ndof = %d | FE_ux.ndof = %d | FE_uy.ndof = %d |  FE_u_eg.ndof = %d \n",
	     FE.ndof , u_ndof, FE_ux.ndof  , FE_uy.ndof , FE_u_eg.ndof );
      printf("FE.ndof = %d | p_ndof = %d | FE_p.ndof = %d  \n",
	     FE.ndof , p_ndof, FE_p.ndof);


      printf("FE.nbdof = %d \n", FE.nbdof);

      printf("###########\n");
      printf("CG DOF for Mechanics = %d\n",   dof_per_cycle_CG[cycle] );
      printf("EG DOF for Mechanics = %d\n",   dof_per_cycle_EG[cycle] );
      printf("###########\n");

      printf("###########\n");
      printf("CG DOF for Pressure = %d\n",   dof_per_cycle_CG_p[cycle] );
      printf("EG DOF for Pressure = %d\n",   dof_per_cycle_EG_p[cycle] );
      printf("###########\n");


      // Get the minimum mesh size
      int zz=0;
      double min_mesh_size = 10000000.;
      double tmp_size =0;
      mesh_struct *mesh_2 = &mesh;

      for (zz=0; zz<mesh_2->nface; zz++) {

	tmp_size=mesh_2->f_area[zz];

	if(tmp_size < min_mesh_size)
	  min_mesh_size = tmp_size;
      }

      mesh_size_per_cycle[cycle] = min_mesh_size;
      if(time == 0)
	{
	  //Set Initial Condition
	  //Inside of the assemble for p0 ...
	}

      dvector sol = dvec_create(ndof);
      dvector old_timestep_sol = dvec_create(ndof);

      //for(time = timestep; timestep_number < total_timestep; time += timestep){
      dvec_cp(&sol, &old_timestep_sol);
      timestep_number = timestep_number+1;
      printf(" ---  CYCLE = %d  ----------- TIME STEPPING TIME = %f  (timestep # = %d | %d) ------------------- \n", \
	     cycle,time,timestep_number,total_timestep);
      //printf(" [Time Step = %f]  \n", timestep);
      //fflush(stdout);

      clock_t clk_assembly_start = clock();

      // Allocate the right-hand and declare the csr matrix
      dvector b;
      // Put into Block Form
      block_dCSRmat A;
      bdcsr_alloc(dim+1+1,dim+1+1,&A);
      /*** Assemble the matrix and right hand side *******************************/
      printf("Assembling the matrix and right-hand side:\n");fflush(stdout);


      assemble_global_block_neighbor(&A,&b,old_timestep_sol.val,	\
				     local_assembly_Elasticity_FACE,	\
				     local_assembly_Elasticity,		\
				     FEM_Block_RHS_Local_Elasticity,	\
				     &FE,				\
				     &mesh,				\
				     cq,source2D,exact_sol2D,Dexact_sol2D, exact_sol2D_dt, time,timestep);
      printf("\n------ Assemble done: \n");fflush(stdout);
      //printf("cycle = %d -- total = %d \n", cycle, total_num_cycle);
      /*
	if((cycle) == total_num_cycle-1)
	{

	FILE* fid;
	fid = fopen("matrix.txt","w");

	dCSRmat Amerge = bdcsr_2_dcsr(&A);
	csr_print_matlab(fid,&Amerge);
	dcsr_free(&Amerge);
	exit(0);
	}
      */
      // Eliminate boundary conditions in matrix and rhs
      if(!BOOL_WEAKLY_IMPOSED_BC)
	eliminate_DirichletBC_blockFE_blockA(bc2D,&FE,&mesh,&b,&A,time);
      /**************************************************/
      //  Apply Pressure "BCs" (removes singularity)
      /*
	REAL pressureval =0.;
	INT pressureloc = 0;

	clock_t clk_assembly_end = clock();
	printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL)
	(clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);

	// Prepare diagonal blocks
	dCSRmat *A_diag;
	A_diag = (dCSRmat *)calloc(dim+1, sizeof(dCSRmat));

	for(i=0;i<dim;i++){ // copy block diagonal to A_diag
	dcsr_alloc(A.blocks[i*(dim+2)]->row, A.blocks[i*(dim+2)]->col, A.blocks[i*(dim+2)]->nnz, &A_diag[i]);
	dcsr_cp(A.blocks[i*(dim+2)], &A_diag[i]);
	}

	// Get Mass Matrix for p
	dCSRmat Mp;
	assemble_global(&Mp,NULL,assemble_mass_local,&FE_p,&mesh,cq,NULL,one_coeff_scal,0.0);
	dcsr_alloc(Mp.row, Mp.col, Mp.nnz, &A_diag[dim]);
	dcsr_cp(&Mp, &A_diag[dim]);

	printf("Solving the System:\n");
	clock_t clk_solve_start = clock();

	INT solver_flag = -20;

	// Allocate solution
	dvector sol = dvec_create(ndof);
	dvector v_ux = dvec_create(FE_ux.ndof);
	dvector v_uy = dvec_create(FE_uy.ndof);
	dvector v_uz;
	if(dim==3) v_uz = dvec_create(FE_uz.ndof);
	dvector v_p = dvec_create(FE_p.ndof);

	// Set initial guess to be all zero
	dvec_set(sol.row, &sol, 0.0);
	// Set initial guess for pressure to match the known "boundary condition" for pressure
	if(dim==2) sol.val[FE_ux.ndof + FE_uy.ndof + pressureloc]  = pressureval;
	if(dim==3) sol.val[FE_ux.ndof + FE_uy.ndof + FE_uz.ndof + pressureloc]  = pressureval;

	// Set parameters for linear iterative methods
	linear_itsolver_param linear_itparam;
	param_linear_solver_init(&linear_itparam);
	param_linear_solver_set(&linear_itparam,&inparam);
	INT solver_type = linear_itparam.linear_itsolver_type;
	INT solver_printlevel = linear_itparam.linear_print_level;

	// Solve
	if(solver_type==0) { // Direct Solver
	solver_flag = block_directsolve_UMF(&A,&b,&sol,solver_printlevel);
	} else { // Iterative Solver
	if (linear_itparam.linear_precond_type == PREC_NULL) {
	solver_flag = linear_solver_bdcsr_krylov(&A, &b, &sol, &linear_itparam);
	} else {
	if(dim==2) solver_flag = linear_solver_bdcsr_krylov_block_3(&A, &b, &sol, &linear_itparam, NULL, A_diag);
	if(dim==3) solver_flag = linear_solver_bdcsr_krylov_block_4(&A, &b, &sol, &linear_itparam, NULL, A_diag);
	}
	}

	// Error Check
	if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n",solver_flag);

	clock_t clk_solve_end = clock();
	printf("Elapsed CPU Time for Solve = %f seconds.\n\n",
	(REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);

      */

      /**************************************************/
      //  Apply Pressure "BCs" (removes singularity)

      //	  REAL pressureval =0.;
      //	  INT pressureloc = 0;

      clock_t clk_assembly_end = clock();
      printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL)
	     (clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);

      printf("Solving the System:\n");fflush(stdout);
      clock_t clk_solve_start = clock();

      //*****************************************
      //    SOLVER-SOLVER AND average P	=0
      //*****************************************
      INT jj=-10,solver_flag=-20;
      REAL pmin=1e20,pmax=-1e20,sum=-1e20;
      void *numeric=NULL;
      sol.row++;
      sol.val=realloc(sol.val,sol.row*sizeof(REAL));
      dvec_set(sol.row, &sol, 0.0);
      b.row++;
      b.val=realloc(b.val,b.row*sizeof(REAL));
      b.val[b.row-1]=0e0;
      //sol.val[FE_ux.ndof + FE_uy.ndof + pressureloc]  = pressureval;
      //
      //	  add row and column for the pressure block
      // which is the pressure block?Ans: dim+1;
      // extend

      //-----------------------
      // set paramters for linear solver
      //-----------------------
      linear_itsolver_param linear_itparam;
      param_linear_solver_init(&linear_itparam);
      param_linear_solver_set(&linear_itparam,&inparam);
      INT solver_type = linear_itparam.linear_itsolver_type;
      INT solver_printlevel = linear_itparam.linear_print_level;

      // Set parameters for AMG
      AMG_param amgparam;
      param_amg_init(&amgparam);
      param_amg_set(&amgparam, &inparam);
      param_amg_print(&amgparam);

      // Prepare diagonal blocks for block preconditioner
    	dCSRmat *A_diag;
    	A_diag = (dCSRmat *)calloc(dim+2, sizeof(dCSRmat)); // number of blocks = dim+2

      // get the blocks
      // for ux, uy, and u_eg, grab the blocks directly
      for(i=0;i<dim+1;i++){ // copy block diagonal to A_diag
        dcsr_alloc(A.blocks[i*(dim+3)]->row, A.blocks[i*(dim+3)]->col, A.blocks[i*(dim+3)]->nnz, &A_diag[i]);
    	  dcsr_cp(A.blocks[i*(dim+3)], &A_diag[i]);
    	}
      // for pressure, use the mass matrix
      dcsr_alloc(p_ndof, p_ndof, p_ndof, &A_diag[dim+1]);
      for (i=0; i<=A_diag[dim+1].row; i++)
        A_diag[dim+1].IA[i] = i;
      for (i=0; i<A_diag[dim+1].row; i++)
        A_diag[dim+1].JA[i] = i;
      for (i=0; i<A_diag[dim+1].row; i++)
        A_diag[dim+1].val[i] = mesh.el_vol[i];

      // Linear Solver
      if(solver_type==0) { // Direct Solver
        bdcsr_extend(&A,mesh.el_vol,mesh.el_vol,(dim+1),1.0,1.0);
        printf("nblocks = %d\n", A.brow);
        getchar();
        numeric=block_factorize_UMF(&A,0);//inparam.print_level);
        // solve
        solver_flag=(INT )block_solve_UMF(&A,
  					&b, // rhs.
  					&sol,  // solution.
  					numeric,
  					0);//     inparam.print_level);
        free(numeric);
      } else { // Iterative Solver
        if (linear_itparam.linear_precond_type == PREC_NULL) {
          solver_flag = linear_solver_bdcsr_krylov(&A, &b, &sol, &linear_itparam);
        } else {
          //if(dim==2) solver_flag = linear_solver_bdcsr_krylov_block_3(&A, &b, &sol, &linear_itparam, NULL, A_diag);
          solver_flag = linear_solver_bdcsr_krylov_block_4(&A, &b, &sol, &linear_itparam, &amgparam, A_diag);
        }
      }



      //// printing for debug
      /* FILE *fptmp=NULL; */
      /* fptmp=fopen("debug/u.data","w"); dvector_print(fptmp,&sol); fclose(fptmp); */
      /* fptmp=fopen("debug/b.data","w"); dvector_print(fptmp,&b); fclose(fptmp); */
      //// end printing for debug
      jj=sol.row - mesh.nelm - 1; // begin index of pressure block
      // after extension
      pmin=sol.val[jj];
      pmax=sol.val[jj];
      sum=0e0;
      for(i=0;i<mesh.nelm;i++){
  	       //	    fprintf(stdout,"\nnel:%7d, dof=%7d: sol=%e",	\
  	       //		    i,jj,sol.val[jj]);
  	        sum+=mesh.el_vol[i]*sol.val[jj];
  	         if(sol.val[jj]<pmin) pmin=sol.val[jj];
  	          if(sol.val[jj]>pmax) pmax=sol.val[jj];
  	           jj++;
      }
      fprintf(stdout,"\nINTEGRAL(p)=%13.6e, min(p)=%11.4e, max(p)=%11.4e\n",sum,pmin,pmax);
      // Error Check
      if (solver_flag < 0) fprintf(stdout,"### ERROR: Solver does not converge with error code = %d!\n", solver_flag);
      b.row--;
      b.val=realloc(b.val,b.row*sizeof(REAL));
      sol.row--;
      sol.val=realloc(sol.val,sol.row*sizeof(REAL));
      //*****************************************
      //    END-SOLVER-SOLVER AND average P=0
      //*****************************************


      /*
      // Prepare diagonal blocks
      dCSRmat *A_diag;
      A_diag = (dCSRmat *)calloc(dim+1, sizeof(dCSRmat));

      for(i=0;i<dim;i++){ // copy block diagonal to A_diag
      dcsr_alloc(A.blocks[i*(dim+2)]->row, A.blocks[i*(dim+2)]->col, A.blocks[i*(dim+2)]->nnz, &A_diag[i]);
      dcsr_cp(A.blocks[i*(dim+2)], &A_diag[i]);
      }
      // Get Mass Matrix for p
      dCSRmat Mp;
      assemble_global(&Mp,NULL,assemble_mass_local,&FE_p,&mesh,cq,NULL,one_coeff_scal,0.0);
      dcsr_alloc(Mp.row, Mp.col, Mp.nnz, &A_diag[dim]);
      dcsr_cp(&Mp, &A_diag[dim]);

      // Set parameters for linear iterative methods
      linear_itsolver_param linear_itparam;
      param_linear_solver_init(&linear_itparam);
      param_linear_solver_set(&linear_itparam,&inparam);
      INT solver_type = linear_itparam.linear_itsolver_type;
      INT solver_printlevel = linear_itparam.linear_print_level;

      // Solve
      if(solver_type==0) { // Direct Solver
      solver_flag = block_directsolve_UMF(&A,&b,&sol,solver_printlevel);
      } else { // Iterative Solver
      if (linear_itparam.linear_precond_type == PREC_NULL) {
      solver_flag = linear_solver_bdcsr_krylov(&A, &b, &sol, &linear_itparam);
      } else {
      if(dim==2) solver_flag = linear_solver_bdcsr_krylov_block_3(&A, &b, &sol, &linear_itparam, NULL, A_diag);
      if(dim==3) solver_flag = linear_solver_bdcsr_krylov_block_4(&A, &b, &sol, &linear_itparam, NULL, A_diag);
      }
      }
      */
      /////////////////////////////////////////////////////////////////
      // Error Check
      if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n",solver_flag);



      clock_t clk_solve_end = clock();
      printf("Elapsed CPU Time for Solve = %f seconds.\n\n",
	     (REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);


      //if(timestep_number == total_timestep){
      {
	printf("Compute Error at Time = %f \n", time);

	//////////////////////////////////////////////////////////////////////////////////////////////
	/********************* Compute Errors if you have exact solution ****************************/
	clock_t clk_error_start = clock();

	REAL* solerrL2 = (REAL *) calloc(dim+1+1, sizeof(REAL));
	REAL* solerrH1 = (REAL *) calloc(dim+1+1, sizeof(REAL)); // Note: No H1 error for P0 elements
	REAL* solerr_stress = (REAL *) calloc(dim+1+1, sizeof(REAL)); // Note: No H1 error for P0 elements

	L2error_block(solerrL2, sol.val, exact_sol2D, &FE, &mesh, cq, time);
	//HDerror_block(solerrH1, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
	HDsemierror_block(solerrH1, sol.val, Dexact_sol2D, &FE, &mesh, cq, time);
	//NEW SLEE Aug 17 2020
	HDsemierror_block_Stress(solerr_stress, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, time);

	REAL uerrL2 = 0;
	REAL uerrH1 = 0;
	REAL uerr_stressH1 = 0;

	//For Pressure
	REAL uerrL2_p = 0;
	REAL uerrH1_p = 0;

	for(i=0;i<dim;i++)
	  uerrL2 += solerrL2[i]*solerrL2[i];
	for(i=0;i<dim;i++)
	  uerrH1 += solerrH1[i]*solerrH1[i];
	for(i=0;i<dim;i++)
	  uerr_stressH1 += solerr_stress[i]*solerr_stress[i];

	uerrL2 = sqrt(uerrL2);
	uerrH1 = sqrt(uerrH1);
	uerr_stressH1 = sqrt(uerr_stressH1);

	//For Pressure
	//uerrL2_p += solerrL2[3]*solerrL2[3];
	//uerrH1_p += solerrL2[4]*solerrL2[4];
	//uerrL2_p = sqrt(uerrL2_p);
	//uerrH1_p = sqrt(uerrH1_p);

	uerrL2_p = sqrt(solerrL2[3]);
	uerrH1_p = sqrt(solerrH1[3]);


	//REAL perrL2 = solerrL2[dim];
	//REAL perrH1 = solerrH1[dim];

	printf("************* MECHANCIS   *****************************\n");
	printf("[CG] L2 Norm of u error    = %26.13e\n",uerrL2);
	printf("[CG] H1 Norm of u error    = %26.13e\n",uerrH1);
	printf("[CG] Stress Norm of u error    = %26.13e\n",uerrH1);
	printf("*******************************************************\n\n");

	//Jul. 10. 2020 SLEE save the errors for convergence computation
	L2_error_per_cycle[cycle] = uerrL2;
	H1_error_per_cycle[cycle] = uerrH1;
	H1_stress_error_per_cycle[cycle] = uerr_stressH1;

	//For Pressure
	L2_error_p_per_cycle[cycle] = uerrL2_p;
	H1_error_p_per_cycle[cycle] = uerrH1_p;

	printf("************* PRESSURE    *****************************\n");
	printf("[CG] L2 Norm of u error    = %26.13e\n",uerrL2_p);
	printf("[CG] H1 Norm of u error    = %26.13e\n",uerrH1_p);
	//printf("[CG] Stress Norm of u error    = %26.13e\n",uerrH1);
	printf("*******************************************************\n\n");


	//L2_error_p_per_cycle[cycle] = perrL2;
	//NEW SLEE Aug 23 2020
	//NEW ERROR FOR EG

	REAL* solerrL2_EG = (REAL *) calloc(dim+1+1, sizeof(REAL));
	REAL* solerrH1_EG = (REAL *) calloc(dim+1+1, sizeof(REAL)); // Note: No H1 error for P0 elements
	REAL* solerr_stress_EG = (REAL *) calloc(dim+1+1, sizeof(REAL)); // Note: No H1 error for P0 elements
	REAL* solerr_energy_EG = (REAL *) calloc(dim+1+1, sizeof(REAL)); // Note: No H1 error for P0 elements

	L2error_block_EG(solerrL2_EG, sol.val, exact_sol2D, &FE, &mesh, cq, time);
	HDerror_block_EG(solerrH1_EG, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, time);
	HDsemierror_block_Stress_EG(solerr_stress_EG, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, time);
	//HDsemierror_block_EnergyNorm_EG(solerr_energy_EG, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
	HDsemierror_block_EnergyNorm_EG_FaceLoop(solerr_energy_EG, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, time);

	//NEW SLEE Aug 17 2020
	//NEW ERROR FOR STRESS, \mu \epsilon(u) + \lambda \div u
	REAL uerrL2_EG = 0;
	REAL uerrH1_EG = 0;
	REAL uerr_stress_EG = 0;
	REAL uerr_energy_EG = 0;

	// For Pressure
	REAL uerrL2_EG_p = 0;
	REAL uerrH1_EG_p = 0;


	for(i=0;i<dim;i++)
	  uerrL2_EG += solerrL2_EG[i]*solerrL2_EG[i];
	for(i=0;i<dim;i++)
	  uerrH1_EG += solerrH1_EG[i]*solerrH1_EG[i];
	for(i=0;i<dim;i++)
	  uerr_stress_EG += solerr_stress_EG[i]*solerr_stress_EG[i];
	//for(i=0;i<dim;i++)
	//uerr_energy_EG += solerr_energy_EG[i]*solerr_energy_EG[i];

	uerr_energy_EG = solerr_energy_EG[0]+solerr_energy_EG[1];
	//DEBUG
	uerr_energy_EG += uerrH1_EG;

	// For Pressure
	uerrL2_EG_p = solerrL2_EG[3];//*solerrL2_EG[3];
	uerrH1_EG_p = solerrH1_EG[3];//*solerrL2_EG[4];

	// DEUBG L2 - H1
	uerrL2_EG = sqrt(uerrL2_EG);
	uerrH1_EG = sqrt(uerrH1_EG);
	uerr_stress_EG = sqrt(uerr_stress_EG);
	uerr_energy_EG = sqrt(uerr_energy_EG);

	// For Pressure
	uerrL2_EG_p = sqrt(uerrL2_EG_p);
	uerrH1_EG_p = sqrt(uerrH1_EG_p);

	L2_EG_error_per_cycle[cycle] = uerrL2_EG;
	H1_EG_error_per_cycle[cycle] = uerrH1_EG;
	H1_stress_EG_error_per_cycle[cycle] = uerr_stress_EG;
	H1_energy_EG_error_per_cycle[cycle] = uerr_energy_EG;

	// For Pressure
	//L2_p_EG_error_per_cycle[cycle] = uerrL2_EG_p;
	//H1_p_EG_error_per_cycle[cycle] = uerrH1_EG_p;

	printf("# of elements = %d \n", mesh.nelm);
	printf("*************     MECHANCIS   **************************\n");
	printf("L2 Norm of u (EG) error    = %26.13e\n",uerrL2_EG);
	printf("H1 Norm of u (EG) error    = %26.13e\n",uerrH1_EG);
	printf("Stress Norm of u (EG) error    = %26.13e\n",uerr_stress_EG);
	printf("Energy Norm of u (EG) error    = %26.13e\n",uerr_energy_EG);
	printf("*******************************************************\n");
	printf("*************     PRESSURE   **************************\n");
	printf("L2 Norm of u (EG) error    = %26.13e\n",uerrL2_EG_p);
	printf("H1 Norm of u (EG) error    = %26.13e\n",uerrH1_EG_p);


	/////////////////////////////////////
	/*fprintf(stdout,"\n%d,%d\n\n",A.brow,A.bcol);
	  INT j;
	  for(i=0;i<A.brow;i++){
	  for(j=0;j<A.brow;j++){
	  fprintf(stdout,"\n(%d,%d):::%d,%d,%d\n\n",i,j,
	  A.blocks[j*A.bcol+i]->row,
	  A.blocks[j*A.bcol+i]->col,
	  A.blocks[j*A.bcol+i]->nnz);
	  }
	  }

	  fflush(stdout);
	*/
	clock_t clk_error_end = clock();
	printf("Elapsed CPU time for getting errors = %lf seconds.\n\n",(REAL)
	       (clk_error_end-clk_error_start)/CLOCKS_PER_SEC);

	/*******************************************************************************************/
	/// Plotting
	// Allocate solution
	dvector v_ux = dvec_create(FE_ux.ndof);
	dvector v_uy = dvec_create(FE_uy.ndof);
	dvector v_uz;
	if(dim==3) v_uz = dvec_create(FE_uz.ndof);
	dvector v_u_eg = dvec_create(FE_u_eg.ndof);

	dvector v_p = dvec_create(FE_p.ndof);
	//dvector v_p_eg = dvec_create(FE_p_eg.ndof);


	get_unknown_component(&v_ux,&sol,&FE,0);
	get_unknown_component(&v_uy,&sol,&FE,1);
	if(dim==3) get_unknown_component(&v_uz,&sol,&FE,2);
	get_unknown_component(&v_u_eg,&sol,&FE,dim);

	get_unknown_component(&v_p,&sol,&FE,3);
	//get_unknown_component(&v_p_eg,&sol,&FE,4);


	char** varname;

	//printf("OUTPUT?\n");


	if(inparam.print_level > 3){

	  char output_filename_per_cycle[512]={'\0'};
	  sprintf( output_filename_per_cycle, "output/solution_%d_%d.vtu", cycle,timestep_number);
	  char* soldump = output_filename_per_cycle;//"output/solution.vtu";
	  //char* soldump = "output/solution.vtu";

	  varname = malloc(5*FE.nspaces*sizeof(char *));
	  varname[0] = "ux";
	  varname[1] = "uy";
	  if(dim==3) varname[2] = "uz";
	  varname[dim] = "u_eg ";

	  varname[3] = "p";
	  varname[4] = "p_eg";

	  dump_blocksol_vtk(soldump,varname,&mesh,&FE,sol.val);

	  // Print in Matlab format to show vector field in nice way.
	  //if(dim==3)
	  //print_matlab_vector_field(&v_ux,&v_uy,&v_uz,&FE_ux);
	}

	if(solerrL2) free(solerrL2);
	if(solerrL2_EG) free(solerrL2_EG);
	if( solerr_stress ) free( solerr_stress);
	if(solerrH1) free(solerrH1);
	dvec_free( &v_ux );
	dvec_free( &v_uy );
	if(dim==3) dvec_free( &v_uz );
	dvec_free( &v_u_eg );
	dvec_free( &v_p );
	//dvec_free( &v_p_eg );


      }//if timestep == 10

      //
      bdcsr_free( &A );
      dvec_free( &b );
      // Quadrature


      //	}//Time Loop

      //RESET TIME
      time = 0.;
      timestep_number = 0;
      //timestep = timestep/2.;
      //total_timestep = total_timestep *2;
      /************ Free All the Arrays ***********************************************************/
      // CSR
      /*
	bdcsr_free( &A );
	if(solerrL2) free(solerrL2);
	if(solerrL2_EG) free(solerrL2_EG);
	if( solerr_stress ) free( solerr_stress);
	if(solerrH1) free(solerrH1);

	dvec_free( &b );
	dvec_free( &sol );
	dvec_free( &v_ux );
	dvec_free( &v_uy );
	if(dim==3) dvec_free( &v_uz );
	dvec_free( &v_u_eg );
	dvec_free( &v_p );
	dvec_free( &v_p_eg );
      */

      dvec_free( &sol );
      dvec_free( &old_timestep_sol );
      // FE Spaces
      free_fespace(&FE_ux);
      free_fespace(&FE_uy);
      if(dim==3) free_fespace(&FE_uz);
      free_fespace(&FE_u_eg);
      free_fespace(&FE_p);
      //free_fespace(&FE_p_eg);

      //free_blockfespace(&FE);

      // Quadrature
      if(cq){
	free_qcoords(cq);
	free(cq);
	cq=NULL;
      }

      // Mesh
      free_mesh(&mesh);
      //*/
      // Strings

      //if(inparam.print_level > 3){
      //if(varname) free(varname);
      //}

      /*******************************************************************************************/
      clock_t clk_overall_end = clock();
      printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",
	     (REAL) (clk_overall_end-clk_overall_start)/CLOCKS_PER_SEC);


    }//SLEE :: cycle loop end


  // Jul.10.2020 SLEE: Error Computation
  //SLEE compute the convergence rate and print
  INT tmp;
  for(tmp=0; tmp<total_num_cycle; ++tmp)
    {
      if(tmp == 0){
	L2_conv_rate_per_cycle[tmp] = 0;
	H1_conv_rate_per_cycle[tmp] = 0;
	H1_stress_conv_rate_per_cycle[tmp] = 0;

	L2_EG_conv_rate_per_cycle[tmp] = 0;
	H1_EG_conv_rate_per_cycle[tmp] = 0;
	H1_stress_EG_conv_rate_per_cycle[tmp] = 0;
	H1_energy_EG_conv_rate_per_cycle[tmp] = 0;

	//For Pressure
	L2_p_conv_rate_per_cycle[tmp] = 0;
	H1_p_conv_rate_per_cycle[tmp] = 0;

	//L2_p_EG_conv_rate_per_cycle[tmp] = 0;
	//H1_p_EG_conv_rate_per_cycle[tmp] = 0;

      }
      else{
	// multiplied dim since we use DOFs not h here.
	L2_conv_rate_per_cycle[tmp] = global_dim_space * (log(L2_error_per_cycle[tmp]) -log(L2_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_CG[tmp-1]) -log(dof_per_cycle_CG[tmp]) );
	H1_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_error_per_cycle[tmp]) -log(H1_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_CG[tmp-1]) -log(dof_per_cycle_CG[tmp]) );
	H1_stress_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_stress_error_per_cycle[tmp]) -log(H1_stress_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_CG[tmp-1]) -log(dof_per_cycle_CG[tmp]) );

	L2_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(L2_EG_error_per_cycle[tmp]) -log(L2_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG[tmp-1]) -log(dof_per_cycle_EG[tmp]) );
	H1_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_EG_error_per_cycle[tmp]) -log(H1_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG[tmp-1]) -log(dof_per_cycle_EG[tmp]) );

	H1_stress_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_stress_EG_error_per_cycle[tmp]) -log(H1_stress_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG[tmp-1]) -log(dof_per_cycle_EG[tmp]) );
	H1_energy_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_energy_EG_error_per_cycle[tmp]) -log(H1_energy_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG[tmp-1]) -log(dof_per_cycle_EG[tmp]) );

	//For Pressure
	L2_p_conv_rate_per_cycle[tmp] = global_dim_space * (log(L2_error_p_per_cycle[tmp]) -log(L2_error_p_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_CG_p[tmp-1]) -log(dof_per_cycle_CG_p[tmp]) );
	H1_p_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_error_p_per_cycle[tmp]) -log(H1_error_p_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_CG_p[tmp-1]) -log(dof_per_cycle_CG_p[tmp]) );

	//L2_p_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(L2_p_EG_error_per_cycle[tmp]) -log(L2_p_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG_p[tmp-1]) -log(dof_per_cycle_EG_p[tmp]) );
	//H1_p_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_p_EG_error_per_cycle[tmp]) -log(H1_p_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG_p[tmp-1]) -log(dof_per_cycle_EG_p[tmp]) );


      }


      printf("****** CYCLE = %d ****** \n", tmp);
      printf("----- MECHANCICS   ----------------------------------\n");

      printf("L2 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, L2_error_per_cycle[tmp], dof_per_cycle_CG[tmp],L2_conv_rate_per_cycle[tmp]);
      printf("H1 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_error_per_cycle[tmp], dof_per_cycle_CG[tmp],H1_conv_rate_per_cycle[tmp]);
      printf("H1 Stress Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_stress_error_per_cycle[tmp], dof_per_cycle_CG[tmp],H1_stress_conv_rate_per_cycle[tmp]);

      printf("EG - L2 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, L2_EG_error_per_cycle[tmp], dof_per_cycle_EG[tmp],
	     L2_EG_conv_rate_per_cycle[tmp]);
      printf("EG - H1 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_EG_error_per_cycle[tmp], dof_per_cycle_EG[tmp],
	     H1_EG_conv_rate_per_cycle[tmp]);

      printf("EG - Stress_Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_stress_EG_error_per_cycle[tmp], dof_per_cycle_EG[tmp],
	     H1_stress_EG_conv_rate_per_cycle[tmp]);

      printf("EG - EnergyNorm_Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_energy_EG_error_per_cycle[tmp], dof_per_cycle_EG[tmp],
	     H1_energy_EG_conv_rate_per_cycle[tmp]);

      printf("----------------------------------------------------\n");


      printf("-----  PRESSURE   ----------------------------------\n");

      printf("L2 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, L2_error_p_per_cycle[tmp], dof_per_cycle_CG_p[tmp],L2_p_conv_rate_per_cycle[tmp]);
      printf("H1 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_error_p_per_cycle[tmp], dof_per_cycle_CG_p[tmp],H1_p_conv_rate_per_cycle[tmp]);

      //printf("EG - L2 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, L2_p_EG_error_per_cycle[tmp], dof_per_cycle_EG_p[tmp],L2_p_EG_conv_rate_per_cycle[tmp]);
      //printf("EG - H1 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_p_EG_error_per_cycle[tmp], dof_per_cycle_EG_p[tmp],H1_p_EG_conv_rate_per_cycle[tmp]);

      printf("----------------------------------------------------\n");

    }

  /*
  //for Latex Print Table
  printf("** LATEX TABLE CG MECHANICS ** \n");
  for(int tmp=0; tmp<total_num_cycle; ++tmp)
  {
  printf("%d & %f & %f &  %f &  %f   &  %f &  %f \\\\ \\hline \n", dof_per_cycle_CG[tmp], L2_error_per_cycle[tmp], L2_conv_rate_per_cycle[tmp],H1_error_per_cycle[tmp], H1_conv_rate_per_cycle[tmp],H1_stress_error_per_cycle[tmp], H1_stress_conv_rate_per_cycle[tmp] );
  }

  printf("** LATEX TABLE CG PRESSURE ** \n");
  for(int tmp=0; tmp<total_num_cycle; ++tmp)
  {
  printf("%d & %f & %f &  %f &  %f   \\\\ \\hline \n", dof_per_cycle_CG[tmp], L2_error_p_per_cycle[tmp], L2_p_conv_rate_per_cycle[tmp],H1_error_p_per_cycle[tmp], H1_p_conv_rate_per_cycle[tmp]);
  }

  printf("** LATEX TABLE ALL CG PRESSURE ** \n");
  for(int tmp=0; tmp<total_num_cycle; ++tmp)
  {
  printf("%d & %f & %f &  %f &  %f   &  %f &  %f \\\\ \\hline \n", dof_per_cycle_CG[tmp], L2_error_per_cycle[tmp], L2_conv_rate_per_cycle[tmp],H1_error_per_cycle[tmp], H1_conv_rate_per_cycle[tmp],H1_error_p_per_cycle[tmp], H1_p_conv_rate_per_cycle[tmp] );
  }


  for(int tmp=0; tmp<total_num_cycle; ++tmp)
  {
  printf("%f \n",  mesh_size_per_cycle[tmp]);
  }
  */

  /*
    printf("** LATEX TABLE EG MECHANICS ** \n");
    for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
    printf("%d & %f & %f &  %f &  %f   &  %f &  %f  & %f & %f \\\\ \\hline \n",  dof_per_cycle_EG[tmp],L2_EG_error_per_cycle[tmp], L2_EG_conv_rate_per_cycle[tmp],H1_EG_error_per_cycle[tmp], H1_EG_conv_rate_per_cycle[tmp],
    H1_stress_EG_error_per_cycle[tmp], H1_stress_EG_conv_rate_per_cycle[tmp], H1_energy_EG_error_per_cycle[tmp], H1_energy_EG_conv_rate_per_cycle[tmp] );
    }

    printf("** LATEX TABLE EG PRESSURE ** \n");
    for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
    printf("%d & %f & %f &  %f &  %f  \\\\ \\hline \n", dof_per_cycle_EG[tmp], L2_p_EG_error_per_cycle[tmp], L2_p_EG_conv_rate_per_cycle[tmp],H1_p_EG_error_per_cycle[tmp], H1_p_EG_conv_rate_per_cycle[tmp]);
    }
  */
  /*
    printf("** LATEX TABLE ALL CG MECHANICS && EG PRESSURE ** \n");
    for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
    printf("%f & %d & %f & %f &  %f &  %f  & %d  & %f  & %f  \\\\ \\hline \n",mesh_size_per_cycle[tmp],  dof_per_cycle_CG[tmp],L2_error_per_cycle[tmp], L2_conv_rate_per_cycle[tmp],H1_error_per_cycle[tmp], H1_conv_rate_per_cycle[tmp],
    dof_per_cycle_EG[tmp], H1_p_EG_error_per_cycle[tmp], H1_p_EG_conv_rate_per_cycle[tmp]);
    }


    printf("** LATEX TABLE ALL EG MECHANICS && EG PRESSURE ** \n");
    for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
    printf("%f & %d & %f & %f &  %f &  %f  & %f  & %f  \\\\ \\hline \n", mesh_size_per_cycle[tmp], dof_per_cycle_EG[tmp],L2_EG_error_per_cycle[tmp], L2_EG_conv_rate_per_cycle[tmp],H1_energy_EG_error_per_cycle[tmp], H1_energy_EG_conv_rate_per_cycle[tmp],
    H1_p_EG_error_per_cycle[tmp], H1_p_EG_conv_rate_per_cycle[tmp]);
    }

    printf("** NO L2 CG ** \n");
    for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
    printf("%f & %d & %f &  %f  & %d  & %f  & %f  \\\\ \\hline \n", mesh_size_per_cycle[tmp], dof_per_cycle_CG[tmp],H1_error_per_cycle[tmp], H1_conv_rate_per_cycle[tmp],
    dof_per_cycle_EG[tmp], H1_p_EG_error_per_cycle[tmp], H1_p_EG_conv_rate_per_cycle[tmp]);
    }

    printf("** NO L2 EG ** \n");
    for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
    printf("%f & %d &  %f &  %f & %d  & %f  & %f  \\\\ \\hline \n", mesh_size_per_cycle[tmp], dof_per_cycle_EG[tmp],H1_energy_EG_error_per_cycle[tmp], H1_energy_EG_conv_rate_per_cycle[tmp],
    dof_per_cycle_EG[tmp], H1_p_EG_error_per_cycle[tmp], H1_p_EG_conv_rate_per_cycle[tmp]);
    }
  */

  return 0;


}  /* End of Program */
/*********************************************************************************************/
