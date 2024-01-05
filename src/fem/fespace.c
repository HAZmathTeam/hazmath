/*! \file src/fem/fespace.c
*
* \brief Creates and destroys the structres for the finite-element spaces.
*
*  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 2/17/15.
*  Copyright 2015__HAZMATH__. All rights reserved.
*
* \note modified by James Adler 11/14/2016
* \note Updated on 10/3/2018 for 0-1 fix.
*/

#include "hazmath.h"

/*!
* \fn void initialize_fespace(fespace *FE)
*
* \brief Initializes all components of the structure FE.
*
* \return FE     Struct for FE space
*
*/
void initialize_fespace(fespace *FE)
{
  FE->FEtype = -666;
  FE->scal_or_vec = -666;
  FE->dof_form = -666;
  FE->nelm = -666;
  FE->cdof = NULL;
  FE->ndof = -666;
  FE->nbdof = -666;
  FE->dof_per_elm = -666;
  FE->dof_per_edge = -666;
  FE->dof_per_face = -666;
  FE->el_dof = NULL;
  FE->ed_dof = NULL;
  FE->f_dof = NULL;
  FE->dirichlet = NULL;
  FE->dof_flag = NULL;
  FE->periodic = NULL;
  FE->phi = NULL;
  FE->dphi = NULL;
  FE->ddphi = NULL;

  return;
}
/******************************************************************************/

/*!
* \fn void create_fespace(fespace *FE,mesh_struct* mesh,INT FEtype)
*
* \brief Allocates memory and properties of fespace struct.
*
* \param mesh      Mesh struc
* \param FEtype    Element Type:
*                  Implemented:  0-2: PX | 20: Ned-0 | 30: RT-0 | 62: P1 + facebubble | 99: single DoF for constraints
*                  TBI: -9--1: DGX | 10-19: QX | 21-29: Ned-X | 31:39: RT-X | 40-59: Ned/RT-X on quads | 61: MINI element | 70-79: Vector of PX scalars
*
* \return FE       Struct for FE space
*
*/
void create_fespace(fespace *FE,mesh_struct* mesh,INT FEtype)
{
  // Flag for errors
  SHORT status;
  INT i;

  // Initialize First
  initialize_fespace(FE);

  // Set parameters
  FE->FEtype = FEtype;
  FE->nelm = mesh->nelm;
  INT dim = mesh->dim;

  INT* dirichlet;
  INT* dof_flag;
  REAL* phi;
  REAL* dphi;

  /**/
  iCSRmat* temp;
  temp = (iCSRmat *) malloc(sizeof(struct iCSRmat));
  iCSRmat ed_el;
  iCSRmat f_el;
  iCSRmat ed_f;
  switch (FEtype)
  {
    // Contants - P0
    case 0:
    FE->scal_or_vec = 0;
    FE->ndof = mesh->nelm;
    coordinates *barycenter = allocatecoords(mesh->nelm,dim);
    for (i=0; i<mesh->nelm; i++) {
      barycenter->x[i] = mesh->el_mid[i*dim];
      if(mesh->dim>1)
      barycenter->y[i] = mesh->el_mid[i*dim + 1];
      if (mesh->dim>2)
      barycenter->z[i] = mesh->el_mid[i*dim + 2];
    }
    FE->cdof = barycenter;
    FE->nbdof = 0;
    FE->dof_per_elm = 1;
    FE->dof_per_face = 1; // This is wrong but needs to be here
    FE->el_dof = malloc(sizeof(struct iCSRmat));
    *(FE->el_dof) = icsr_create_identity(mesh->nelm, 0);
    if(dim>1) {
      icsr_trans(mesh->el_ed,&ed_el);
      FE->ed_dof = malloc(sizeof(struct iCSRmat));
      *(FE->ed_dof) = ed_el;
      icsr_trans(mesh->el_f,&f_el);
      FE->f_dof = malloc(sizeof(struct iCSRmat));
      *(FE->f_dof) = f_el;
    }
    dirichlet = (INT *) calloc(mesh->nelm,sizeof(INT));
    dof_flag = (INT *) calloc(mesh->nelm,sizeof(INT));
    for(i=0;i<mesh->nelm;i++) {
      dirichlet[i] = 0;
      dof_flag[i] = mesh->el_flag[i];
    }
    FE->dirichlet = dirichlet;
    FE->dof_flag = dof_flag;
    phi = (REAL *) calloc(FE->dof_per_elm,sizeof(REAL));
    FE->phi = phi;
    break;

    // Linears - P1
    case 1:
    FE->scal_or_vec = 0;
    FE->cdof = mesh->cv;
    FE->ndof = mesh->nv;
    FE->nbdof = mesh->nbv;
    FE->dof_per_elm = mesh->v_per_elm;
    FE->dof_per_face = mesh->dim;
    FE->dof_per_edge = 2;
    FE->el_dof = mesh->el_v;
    if(dim>1) {
      FE->ed_dof = mesh->ed_v;
      FE->f_dof = mesh->f_v;
    }
    dirichlet = (INT *) calloc(FE->ndof,sizeof(INT));
    dof_flag = (INT *) calloc(FE->ndof,sizeof(INT));
    for(i=0;i<FE->ndof;i++) {
      dirichlet[i] = mesh->v_flag[i];
      dof_flag[i] = mesh->v_flag[i];
    }
    FE->dirichlet = dirichlet;
    FE->dof_flag = dof_flag;
    phi = (REAL *) calloc(FE->dof_per_elm,sizeof(REAL));
    FE->phi = phi;
    dphi = (REAL *) calloc(FE->dof_per_elm*mesh->dim,sizeof(REAL));
    FE->dphi = dphi;
    break;

    // Quadratics - P2
    case 2:
    FE->scal_or_vec = 0;
    FE->ndof = mesh->nv + mesh->nelm; // In 1D
    FE->nbdof = mesh->nbv; // In 1D
    FE->dof_per_elm = mesh->v_per_elm + 1; // In 1D
    FE->dof_per_face = dim + (2*mesh->dim-3);
    FE->el_dof = malloc(sizeof(struct iCSRmat));
    if(mesh->dim>1) { // Not in 1D
      FE->ed_dof = malloc(sizeof(struct iCSRmat));
      FE->f_dof = malloc(sizeof(struct iCSRmat));
      FE->ndof = mesh->nv + mesh->nedge;
      FE->nbdof = mesh->nbv + mesh->nbedge;
      FE->dof_per_elm = mesh->v_per_elm + mesh->ed_per_elm;
    }
    get_P2(FE,mesh);
    phi = (REAL *) calloc(FE->dof_per_elm,sizeof(REAL));
    FE->phi = phi;
    dphi = (REAL *) calloc(FE->dof_per_elm*mesh->dim,sizeof(REAL));
    FE->dphi = dphi;
    break;

    // Nedelec Elements
    case 20:
    FE->scal_or_vec = 1;
    FE->cdof = array_2_coord ( mesh->ed_mid, mesh->nedge, mesh->dim);
    FE->ndof = mesh->nedge;
    FE->nbdof = mesh->nbedge;
    FE->dof_per_elm = mesh->ed_per_elm;
    FE->dof_per_face = (2*mesh->dim-3);
    FE->el_dof = mesh->el_ed;
    FE->ed_dof = malloc(sizeof(struct iCSRmat));
    *(FE->ed_dof) = icsr_create_identity(mesh->nedge, 0);
    FE->f_dof = mesh->f_ed;
    dirichlet = (INT *) calloc(FE->ndof,sizeof(INT));
    dof_flag = (INT *) calloc(FE->ndof,sizeof(INT));
    for(i=0;i<FE->ndof;i++) {
      dirichlet[i] = mesh->ed_flag[i];
      dof_flag[i] = mesh->ed_flag[i];
    }
    FE->dirichlet = dirichlet;
    FE->dof_flag = dof_flag;
    phi = (REAL *) calloc(FE->dof_per_elm*mesh->dim,sizeof(REAL));
    FE->phi = phi;
    if(mesh->dim==2) // Scalar Curl
    dphi = (REAL *) calloc(FE->dof_per_elm,sizeof(REAL));
    else // Vector Curl
    dphi = (REAL *) calloc(FE->dof_per_elm*mesh->dim,sizeof(REAL));
    FE->dphi = dphi;
    break;

    // Raviart-Thomas Elements
    case 30:
    FE->scal_or_vec = 1;
    FE->cdof = array_2_coord ( mesh->f_mid, mesh->nface, mesh->dim);
    FE->ndof = mesh->nface;
    FE->nbdof = mesh->nbface;
    FE->dof_per_elm = mesh->f_per_elm;
    FE->dof_per_face = 1;
    FE->el_dof = mesh->el_f;
    icsr_trans(mesh->f_ed,&ed_f);
    FE->ed_dof = malloc(sizeof(struct iCSRmat));
    *(FE->ed_dof) = ed_f;
    FE->f_dof = malloc(sizeof(struct iCSRmat));
    *(FE->f_dof) = icsr_create_identity(mesh->nface, 0);
    dirichlet = (INT *) calloc(FE->ndof,sizeof(INT));
    dof_flag = (INT *) calloc(FE->ndof,sizeof(INT));
    for(i=0;i<FE->ndof;i++) {
      dirichlet[i] = mesh->f_flag[i];
      dof_flag[i] = mesh->f_flag[i];
    }
    FE->dirichlet = dirichlet;
    FE->dof_flag = dof_flag;
    phi = (REAL *) calloc(FE->dof_per_elm*mesh->dim,sizeof(REAL));
    FE->phi = phi;
    dphi = (REAL *) calloc(FE->dof_per_elm,sizeof(REAL));
    FE->dphi = dphi;
    break;

    // Bubble routines.
    // Numbering will be X0Y: X = polynomial space order; Y = order of bubbles
    //                        Ex: 103 is P1 plus cubic bubble -> MINI element
    // MINI element
    // Note: numbering is as follows: 0, 1, ..., nv-1, nv, nv+1, ..., nv+nelm-1
    case 103:
    FE->scal_or_vec = 0; // Scalar
    FE->ndof = mesh->nv + mesh->nelm;
    FE->nbdof = mesh->nbv;
    FE->dof_per_elm = mesh->v_per_elm + 1;
    FE->dof_per_face = mesh->dim;
    FE->dof_per_edge = 2;
    FE->el_dof = malloc(sizeof(struct iCSRmat));

    // Get coordinates and adjust el_dof maps to account for bubble
    get_MINI(FE,mesh);

    // Edge to DoF and Face to Dof is same as P1, since bubble is at element center
    if(dim>1) {
      FE->ed_dof = mesh->ed_v;
      FE->f_dof = mesh->f_v;
    }

    // Allocate space for basis functions
    phi = (REAL *) calloc(FE->dof_per_elm,sizeof(REAL));
    FE->phi = phi;
    dphi = (REAL *) calloc(FE->dof_per_elm*mesh->dim,sizeof(REAL));
    FE->dphi = dphi;
    break;

    // Vector Velocity TODO: This is not complete yet (Seems to be only for linears though...)
    case 60:
    FE->scal_or_vec = 1;
    FE->cdof = NULL;
    FE->ndof = mesh->nv * dim;
    FE->nbdof = mesh->nbv * dim;
    FE->dof_per_elm = mesh->v_per_elm * mesh->dim;
    FE->dof_per_face = mesh->dim*mesh->dim;
    FE->el_dof = malloc(sizeof(struct iCSRmat));
    //printf("Before Concat\n");
    if(dim==1){
      FE->el_dof = mesh->el_v;
    } else if(dim==2){
      icsr_concat(mesh->el_v, mesh->el_v, FE->el_dof);
    } else if(mesh->dim==3){
      icsr_concat(mesh->el_v,mesh->el_v,temp);
      icsr_concat(temp,mesh->el_v,FE->el_dof);
      icsr_free(temp);
    }
    if(dim>1) {
      FE->ed_dof = malloc(sizeof(struct iCSRmat));
      if(mesh->dim==2){
        icsr_concat(mesh->ed_v, mesh->ed_v, FE->ed_dof);
      } else if(mesh->dim==3) {
        icsr_concat(mesh->ed_v, mesh->ed_v, temp);
        icsr_concat(temp, mesh->ed_v, FE->ed_dof);
        icsr_free(temp);
      }
      FE->f_dof = malloc(sizeof(struct iCSRmat));
      if(mesh->dim==2){
        icsr_concat(mesh->f_v, mesh->f_v, FE->f_dof);
      } else if(mesh->dim==3){
        icsr_concat(mesh->f_v, mesh->f_v, temp);
        icsr_concat(temp, mesh->f_v, FE->f_dof);
        icsr_free(temp);
      }
    }
    dirichlet = (INT *) calloc(FE->ndof,sizeof(INT));
    dof_flag = (INT *) calloc(FE->ndof,sizeof(INT));
    for(i=0;i<FE->ndof;i++) {
      dirichlet[i] = mesh->v_flag[i % mesh->nv];
      dof_flag[i] = mesh->v_flag[i % mesh->nv];
    }
    FE->dirichlet = dirichlet;
    FE->dof_flag = dof_flag;
    phi = (REAL *) calloc(FE->dof_per_elm,sizeof(REAL)); // This is basis functions
    FE->phi = phi;
    dphi = (REAL *) calloc(FE->dof_per_elm*mesh->dim,sizeof(REAL));
    FE->dphi = dphi;
    break;

    // Face bubbles // TODO Make one vector space plus linears
    case 61:
    FE->scal_or_vec = 1;
    FE->cdof = array_2_coord ( mesh->f_mid, mesh->nface, mesh->dim);
    FE->ndof = mesh->nface;
    FE->nbdof = mesh->nbface;
    FE->dof_per_elm = mesh->f_per_elm;
    FE->dof_per_face = 1;
    FE->el_dof = mesh->el_f;
    icsr_trans(mesh->f_ed,&ed_f);
    FE->ed_dof = malloc(sizeof(struct iCSRmat));
    *(FE->ed_dof) = ed_f;
    FE->f_dof = malloc(sizeof(struct iCSRmat));
    *(FE->f_dof) = icsr_create_identity(mesh->nface, 0);
    dirichlet = (INT *) calloc(FE->ndof,sizeof(INT));
    dof_flag = (INT *) calloc(FE->ndof,sizeof(INT));
    for(i=0;i<FE->ndof;i++) {
      dirichlet[i] = mesh->f_flag[i];
      dof_flag[i] = mesh->f_flag[i];
    }
    FE->dirichlet = dirichlet;
    FE->dof_flag = dof_flag;
    phi = (REAL *) calloc(FE->dof_per_elm*mesh->dim,sizeof(REAL));
    FE->phi = phi;
    dphi = (REAL *) calloc(FE->dof_per_elm*mesh->dim*mesh->dim,sizeof(REAL));
    FE->dphi = dphi;
    break;

    // Constraint space
    case 99: // 1 DOF FE Space (i.e., for an integral constraint)
    FE->scal_or_vec = 0;
    FE->ndof = 1;
    coordinates *domain = allocatecoords(1,dim);
    domain->x[0] = 0.0;
    if(mesh->dim>1) domain->y[0] = 0.0;
    if(mesh->dim>2) domain->z[0] = 0.0;
    FE->cdof = domain;
    FE->nbdof = 0;
    FE->dof_per_elm = 1;
    FE->dof_per_face = 1;
    FE->el_dof = malloc(sizeof(struct iCSRmat));
    *(FE->el_dof) = icsr_create(mesh->nelm,1,mesh->nelm);
    for(i=0;i<mesh->nelm;i++) {
      FE->el_dof->IA[i] = i;
      FE->el_dof->JA[i] = 0;
      FE->el_dof->val[i] = 1;
    }
    FE->el_dof->IA[mesh->nelm] = mesh->nelm;
    FE->f_dof = malloc(sizeof(struct iCSRmat));
    *(FE->f_dof) = icsr_create(mesh->nface,1,mesh->nface);
    for(i=0;i<mesh->nface;i++) {
      FE->f_dof->IA[i] = i;
      FE->f_dof->JA[i] = 0;
      FE->f_dof->val[i] = 1;
    }
    FE->f_dof->IA[mesh->nface] = mesh->nface;
    FE->ed_dof = malloc(sizeof(struct iCSRmat));
    *(FE->ed_dof) = icsr_create(mesh->nedge,1,mesh->nedge);
    for(i=0;i<mesh->nedge;i++) {
      FE->ed_dof->IA[i] = i;
      FE->ed_dof->JA[i] = 0;
      FE->ed_dof->val[i] = 1;
    }
    FE->ed_dof->IA[mesh->nedge] = mesh->nedge;
    dirichlet = (INT *) calloc(1,sizeof(INT));
    dof_flag = (INT *) calloc(1,sizeof(INT));
    dirichlet[0] = 0;
    dof_flag[0] = 0;
    FE->dirichlet = dirichlet;
    FE->dof_flag = dof_flag;
    phi = (REAL *) calloc(FE->dof_per_elm,sizeof(REAL));
    phi[0] = 1;
    FE->phi = phi;
    break;

    // Not a valid space
    default:
    status = ERROR_FE_TYPE;
    check_error(status, __FUNCTION__);
  }

  // Assume all DOF are not periodic to start
  FE->periodic = (INT *) calloc(FE->ndof,sizeof(INT));
  for(i=0;i<FE->ndof;i++) {
    FE->periodic[i] = -1;
  }

  // clean temp
  if (temp) {
    free(temp);
    temp = NULL;
  }

  return;
}
/******************************************************************************/

/*!
* \fn void free_fespace(fespace* FE)
*
* \brief Frees memory of arrays of fespace struct
*
* \return FE       Struct for FE space to be freed
*
*/
void free_fespace(fespace* FE)
{
  if(FE->cdof && (FE->FEtype!=1)) { // If Linears, free_mesh will destroy coordinate struct
    free_coords(FE->cdof);
    free(FE->cdof);
    FE->cdof = NULL;
  }

  if(FE->el_dof && (FE->FEtype==2 || FE->FEtype==0 || FE->FEtype==99 || FE->FEtype>=101)) { // If not P2 or P0 or bubbles, free_mesh will destroy el_dof struct
    icsr_free(FE->el_dof);
    free(FE->el_dof);
    FE->el_dof = NULL;
  }

  if(FE->ed_dof && FE->FEtype!=1 && FE->FEtype!=103) { // If Linears, free_mesh will destroy ed_dof struct
    icsr_free(FE->ed_dof);
    free(FE->ed_dof);
    FE->ed_dof = NULL;
  }

  if(FE->f_dof && FE->FEtype!=1 && FE->FEtype!=20 && FE->FEtype!=103) { // If Linears or Nedelec, free_mesh will destroy f_dof
    icsr_free(FE->f_dof);
    free(FE->f_dof);
    FE->f_dof = NULL;
  }

  if(FE->dirichlet) {
    free(FE->dirichlet);
    FE->dirichlet = NULL;
  }

  if(FE->dof_flag) {
    free(FE->dof_flag);
    FE->dof_flag = NULL;
  }

  if(FE->periodic) {
    free(FE->periodic);
    FE->periodic = NULL;
  }

  if(FE->phi) {
    free(FE->phi);
    FE->phi = NULL;
  }

  if(FE->dphi) {
    free(FE->dphi);
    FE->dphi = NULL;
  }

  if(FE->ddphi) {
    free(FE->ddphi);
    FE->ddphi = NULL;
  }

  return;
}
/******************************************************************************/

/*!
* \fn void initialize_fesystem(block_fespace *FE)
*
* \brief Initializes all components of the structure block FE system
*
* \param nspaces  Number of FE spaces
* \param nun      Number of unknowns (includes components of vector functions)
* \param ndof     Total number of DoF in FE system
*
* \return FE     Struct for block FE system
*
*/
void initialize_fesystem(block_fespace *FE,INT nspaces,INT nun,INT ndof,INT nelm)
{
  FE->nspaces = nspaces;
  FE->nun = nun;
  FE->ndof = ndof;
  FE->nelm = nelm;
  FE->var_spaces = (fespace **) calloc(nspaces,sizeof(fespace *));
  FE->dirichlet = NULL;
  FE->dof_flag = NULL;
  FE->simplex_data = (simplex_local_data *) calloc(1,sizeof(simplex_local_data));
  FE->fe_data = (fe_local_data *) calloc(1,sizeof(fe_local_data));

  return;
}
/******************************************************************************/

/*!
* \fn void initialize_localdata_elm(simplex_local_data *simplex_data,fe_local_data *fe_data,mesh_struct* mesh,block_fespace *FE,INT nq1d)
*
* \brief Initializes some of the components of simplex_local_data and
*        fe_local_data (things that come from the global data and are fixed).
*        Assumes data is on elements
*
* \param mesh     mesh data
* \param nq1d     number of quadrature points in a 1D direction
* \param FE       Struct for block FE system
*
* \return simplex_data Local Mesh data on elm
* \return fe_data      Local FE data on elm
*
*/
void initialize_localdata_elm(simplex_local_data *simplex_data,fe_local_data *fe_data,mesh_struct* mesh,block_fespace *FE,INT nq1d)
{
  // counters
  INT i;

  INT dim = mesh->dim;
  INT nspaces = FE->nspaces;
  INT nun = FE->nun;
  INT v_per_elm = mesh->v_per_elm;
  INT ed_per_elm = mesh->ed_per_elm;
  INT f_per_elm = mesh->f_per_elm;
  INT dof_per_elm_tot = 0;

  // Simplex Data first that is fixed
  simplex_data->dim = dim;
  simplex_data->n_v = v_per_elm;
  simplex_data->local_v = (INT *) calloc(v_per_elm,sizeof(INT));
  simplex_data->xv = (REAL *) calloc(v_per_elm*dim,sizeof(REAL));
  simplex_data->n_ed = ed_per_elm;
  simplex_data->local_ed = (INT *) calloc(ed_per_elm,sizeof(INT));
  simplex_data->v_on_ed = (INT *) calloc(2*ed_per_elm,sizeof(INT));
  simplex_data->ed_len = (REAL *) calloc(ed_per_elm,sizeof(REAL));
  simplex_data->ed_tau = (REAL *) calloc(ed_per_elm*dim,sizeof(REAL));
  simplex_data->ed_mid = (REAL *) calloc(ed_per_elm*dim,sizeof(REAL));
  simplex_data->n_f = f_per_elm;
  simplex_data->local_f = (INT *) calloc(f_per_elm,sizeof(INT));
  simplex_data->v_on_f = (INT *) calloc(dim*f_per_elm,sizeof(INT));
  simplex_data->f_area = (REAL *) calloc(f_per_elm,sizeof(REAL));
  simplex_data->f_norm = (REAL *) calloc(f_per_elm*dim,sizeof(REAL));
  simplex_data->f_mid = (REAL *) calloc(f_per_elm*dim,sizeof(REAL));

  // Quadrature and reference mappings
  simplex_data->quad_local = alloc_quadrature(nq1d,1,dim);
  // Fill in with reference element for now to get lambdas (overwite later for elm)
  quad_refelm(simplex_data->quad_local,nq1d,dim);
  INT nq = simplex_data->quad_local->nq;
  simplex_data->lams = (REAL *) calloc(nq*v_per_elm,sizeof(REAL));
  // Lams are ordered for each quadrature point, lam0, lam1, etc.
  // We fill in gradlams, knowing that this will be overwritten later anyway
  REAL* lam_at_q = simplex_data->lams;
  REAL* dlam_at_q = simplex_data->gradlams;
  REAL* qx = simplex_data->quad_local->x;
  for(i=0;i<nq;i++) {
    P1_basis_ref(lam_at_q,dlam_at_q,qx,dim);
    lam_at_q += dim+1;
    dlam_at_q += (dim+1)*dim;
    qx += dim;
  }
  simplex_data->ref_map = (REAL *) calloc(dim*dim,sizeof(REAL));
  simplex_data->gradlams = (REAL *) calloc(nq*v_per_elm*dim,sizeof(REAL));


  // FE local data that is fixed
  fe_data->nspaces = nspaces;
  fe_data->nun = nun;
  fe_data->fe_types = (INT *) calloc(nspaces,sizeof(INT));
  fe_data->scal_or_vec = (INT *) calloc(nspaces,sizeof(INT));
  // Count of DoF on element
  fe_data->n_dof_per_space = (INT *) calloc(nspaces,sizeof(INT));
  // Basis functions
  fe_data->phi = (REAL **) calloc(nspaces,sizeof(REAL *));
  fe_data->dphi = (REAL **) calloc(nspaces,sizeof(REAL *));
  fe_data->ddphi = (REAL **) calloc(nspaces,sizeof(REAL *));
  // Fill in data for each fe space in system
  for(i=0;i<nspaces;i++) {
    fe_data->fe_types[i] = FE->var_spaces[i]->FEtype;
    fe_data->scal_or_vec[i] = FE->var_spaces[i]->scal_or_vec;
    dof_per_elm_tot += FE->var_spaces[i]->dof_per_elm;
    fe_data->n_dof_per_space[i] = FE->var_spaces[i]->dof_per_elm;
    // Point to correct basis function for space
    fe_data->phi[i] = FE->var_spaces[i]->phi;
    fe_data->dphi[i] = FE->var_spaces[i]->dphi;
    fe_data->ddphi[i] = FE->var_spaces[i]->ddphi;
  }
  fe_data->n_dof = dof_per_elm_tot;
  fe_data->local_dof = (INT *) calloc(dof_per_elm_tot,sizeof(INT));
  fe_data->local_dof_flags = (INT *) calloc(dof_per_elm_tot,sizeof(INT));
  fe_data->u_local = (REAL *) calloc(dof_per_elm_tot,sizeof(REAL));

  return;
}
/******************************************************************************/

/*!
* \fn void get_elmlocaldata(simplex_local_data* elm_data,mesh_struct* mesh,INT elm)
*
* \brief On a given element, grab the local mesh data.  We assume fixed values
*        for the mesh have been initialized already, and just filling in element
*        specific data here.
*
* \param mesh     mesh data
* \param elm      Element we are considering
*
* \return elm_data Struct for local simplex data on elm
*
*/
void get_elmlocaldata(simplex_local_data* elm_data,mesh_struct* mesh,INT elm)
{
  // counters
  INT i,j,vind;

  INT dim = mesh->dim;
  INT v_per_elm = mesh->v_per_elm;
  INT ed_per_elm = mesh->ed_per_elm;
  INT f_per_elm = mesh->f_per_elm;

  // Generic elm data
  elm_data->sindex = elm;
  elm_data->vol = mesh->el_vol[elm];
  elm_data->flag = mesh->el_flag[elm];

  // Get vertices on element
  get_incidence_row(elm,mesh->el_v,elm_data->local_v);
  // Temporary fix until coordinates is made dimensionless
  for(i=0;i<v_per_elm;i++) {
    vind = elm_data->local_v[i];
    for(j=0;j<dim;j++) {
      elm_data->xv[i*dim+j] = mesh->cv->x[vind*dim+j];
    }
  }

  // Edge data
  INT v_per_ed = 2;
  INT edge;
  get_incidence_row(elm,mesh->el_ed,elm_data->local_ed);
  for(i=0;i<ed_per_elm;i++) {
    edge = elm_data->local_ed[i];
    get_incidence_row(edge,mesh->ed_v,elm_data->v_on_ed+i*v_per_ed);
    elm_data->ed_len[i] = mesh->ed_len[edge];
    for(j=0;j<dim;j++) {
      elm_data->ed_tau[i*dim+j] = mesh->ed_tau[edge*dim+j];
      elm_data->ed_mid[i*dim+j] = mesh->ed_mid[edge*dim+j];
    }
  }

  // Face data
  INT v_per_f = dim;
  INT face;
  get_incidence_row(elm,mesh->el_f,elm_data->local_f);
  for(i=0;i<f_per_elm;i++) {
    face = elm_data->local_f[i];
    get_incidence_row(face,mesh->f_v,elm_data->v_on_f+i*v_per_f);
    elm_data->f_area[i] = mesh->f_area[face];
    for(j=0;j<dim;j++) {
      elm_data->f_norm[i*dim+j] = mesh->f_norm[face*dim+j];
      elm_data->f_mid[i*dim+j] = mesh->f_mid[face*dim+j];
    }
  }

  // Quadrature and mappings (recall lamda at ref element already built)
  quad_elm_local(elm_data->quad_local,elm_data,elm_data->quad_local->nq1d);
  compute_refelm_mapping(elm_data->ref_map,elm_data->gradlams,elm_data->xv,dim);

  return;
}
/******************************************************************************/

/*!
* \fn void initialize_localdata_face(simplex_local_data *simplex_data,fe_local_data *fe_data,mesh_struct* mesh,block_fespace *FE,INT nq1d)
*
* \brief Initializes some of the components of simplex_local_data and
*        fe_local_data (things that come from the global data and are fixed).
*        Assumes data is on faces
*
* \param mesh     mesh data
* \param nq1d     number of quadrature points in a 1D direction
* \param FE       Struct for block FE system
*
* \return simplex_data Local Mesh data on face
* \return fe_data      Local FE data on face
*
*/
void initialize_localdata_face(simplex_local_data *simplex_data,fe_local_data *fe_data,mesh_struct* mesh,block_fespace *FE,INT nq1d)
{
  // counters
  INT i;

  INT dim = mesh->dim;
  INT nspaces = FE->nspaces;
  INT nun = FE->nun;
  INT v_per_f = dim;
  INT ed_per_f = 2*dim-3;
  INT dof_per_face_tot = 0;

  // Simplex Data first that is fixed
  simplex_data->dim = dim;
  simplex_data->n_v = v_per_f;
  simplex_data->local_v = (INT *) calloc(v_per_f,sizeof(INT));
  simplex_data->xv = (REAL *) calloc(v_per_f*dim,sizeof(REAL));
  simplex_data->n_ed = ed_per_f;
  simplex_data->local_ed = (INT *) calloc(ed_per_f,sizeof(INT));
  simplex_data->v_on_ed = (INT *) calloc(2*ed_per_f,sizeof(INT));
  simplex_data->ed_len = (REAL *) calloc(ed_per_f,sizeof(REAL));
  simplex_data->ed_tau = (REAL *) calloc(ed_per_f*dim,sizeof(REAL));
  simplex_data->ed_mid = (REAL *) calloc(ed_per_f*dim,sizeof(REAL));
  simplex_data->n_f = 1;
  simplex_data->f_area = (REAL *) calloc(1,sizeof(REAL));
  simplex_data->f_norm = (REAL *) calloc(dim,sizeof(REAL));
  simplex_data->f_mid = (REAL *) calloc(dim,sizeof(REAL));

  // Quadrature (will ignore ref maps for now)
  simplex_data->quad_local = alloc_quadrature_bdry(nq1d,1,dim,2);

  // FE local data that is fixed
  fe_data->nspaces = nspaces;
  fe_data->nun = nun;
  fe_data->fe_types = (INT *) calloc(nspaces,sizeof(INT));
  fe_data->scal_or_vec = (INT *) calloc(nspaces,sizeof(INT));
  fe_data->n_dof_per_space = (INT *) calloc(nspaces,sizeof(INT));
  // Fill in data for each fe space in system
  for(i=0;i<nspaces;i++) {
    fe_data->fe_types[i] = FE->var_spaces[i]->FEtype;
    fe_data->scal_or_vec[i] = FE->var_spaces[i]->scal_or_vec;
    dof_per_face_tot += FE->var_spaces[i]->dof_per_face;
    fe_data->n_dof_per_space[i] = FE->var_spaces[i]->dof_per_face;
  }
  fe_data->n_dof = dof_per_face_tot;
  fe_data->local_dof = (INT *) calloc(dof_per_face_tot,sizeof(INT));
  fe_data->local_dof_flags = (INT *) calloc(dof_per_face_tot,sizeof(INT));
  fe_data->u_local = (REAL *) calloc(dof_per_face_tot,sizeof(REAL));

  // Don't construct basis functions on faces (yet...)

  return;
}
/****************************************************************************************/

/****************************************************************************************/
/*!
* \fn void get_facelocaldata(simplex_local_data* face_data,mesh_struct* mesh,INT face)
*
* \brief On a given face, grab the local mesh data.  We assume fixed values
*        for the mesh have been initialized already, and just filling in face
*        specific data here.
*
* \param mesh     mesh data
* \param face     face we are considering
*
* \return face_data Struct for local simplex data on face
*
*/
void get_facelocaldata(simplex_local_data* face_data,mesh_struct* mesh,INT face)
{
  // counters
  INT i,j,vind;

  INT dim = mesh->dim;
  INT v_per_f = dim;
  INT ed_per_f = 2*dim-3;

  // Generic elm data
  face_data->sindex = face;
  face_data->vol = mesh->f_area[face];
  face_data->flag = mesh->f_flag[face];

  // Get vertices on element
  get_incidence_row(face,mesh->f_v,face_data->local_v);
  // Temporary fix until coordinates is made dimensionless
  for(i=0;i<v_per_f;i++) {
    vind = face_data->local_v[i];
    for(j=0;j<dim;j++) {
      face_data->xv[i*dim+j] = mesh->cv->x[vind*dim+j];
    }
  }

  // Edge data
  INT v_per_ed = 2;
  INT edge;
  get_incidence_row(face,mesh->f_ed,face_data->local_ed);
  for(i=0;i<ed_per_f;i++) {
    edge = face_data->local_ed[i];
    get_incidence_row(edge,mesh->ed_v,face_data->v_on_ed+i*v_per_ed);
    face_data->ed_len[i] = mesh->ed_len[edge];
    for(j=0;j<dim;j++) {
      face_data->ed_tau[i*dim+j] = mesh->ed_tau[edge*dim+j];
      face_data->ed_mid[i*dim+j] = mesh->ed_mid[edge*dim+j];
    }
  }

  // Face data
  face_data->f_area[0] = mesh->f_area[face];
  for(j=0;j<dim;j++) {
    face_data->f_norm[j] = mesh->f_norm[face*dim+j];
    face_data->f_mid[j] = mesh->f_mid[face*dim+j];
  }

  // Quadrature (will ignore reference mapping for now)
  quad_face_local(face_data->quad_local,face_data,face_data->quad_local->nq1d,face);

  return;
}
/******************************************************************************/

/*!
* \fn void initialize_localdata_edge(simplex_local_data *simplex_data,fe_local_data *fe_data,mesh_struct* mesh,block_fespace *FE,INT nq1d)
*
* \brief Initializes some of the components of simplex_local_data and
*        fe_local_data (things that come from the global data and are fixed).
*        Assumes data is on edges
*
* \param mesh     mesh data
* \param nq1d     number of quadrature points in a 1D direction
* \param FE       Struct for block FE system
*
* \return simplex_data Local Mesh data on edge
* \return fe_data      Local FE data on edge
*
*/
void initialize_localdata_edge(simplex_local_data *simplex_data,fe_local_data *fe_data,mesh_struct* mesh,block_fespace *FE,INT nq1d)
{
  // counters
  INT i;

  INT dim = mesh->dim;
  INT nspaces = FE->nspaces;
  INT nun = FE->nun;
  INT v_per_ed = 2;
  INT dof_per_edge_tot = 0;

  // Simplex Data first that is fixed
  simplex_data->dim = dim;
  simplex_data->n_v = v_per_ed;
  simplex_data->local_v = (INT *) calloc(v_per_ed,sizeof(INT));
  simplex_data->xv = (REAL *) calloc(v_per_ed*dim,sizeof(REAL));
  simplex_data->n_ed = 1;
  simplex_data->ed_len = (REAL *) calloc(1,sizeof(REAL));
  simplex_data->ed_tau = (REAL *) calloc(dim,sizeof(REAL));
  simplex_data->ed_mid = (REAL *) calloc(dim,sizeof(REAL));

  // Quadrature (will ignore ref maps for now)
  simplex_data->quad_local = alloc_quadrature_bdry(nq1d,1,dim,1);

  // FE local data that is fixed
  fe_data->nspaces = nspaces;
  fe_data->nun = nun;
  fe_data->fe_types = (INT *) calloc(nspaces,sizeof(INT));
  fe_data->fe_types = (INT *) calloc(nspaces,sizeof(INT));
  fe_data->n_dof_per_space = (INT *) calloc(nspaces,sizeof(INT));
  // Fill in data for each fe space in system
  for(i=0;i<nspaces;i++) {
    fe_data->fe_types[i] = FE->var_spaces[i]->FEtype;
    fe_data->scal_or_vec[i] = FE->var_spaces[i]->scal_or_vec;
    dof_per_edge_tot += FE->var_spaces[i]->dof_per_edge;
    fe_data->n_dof_per_space[i] = FE->var_spaces[i]->dof_per_edge;
  }
  fe_data->n_dof = dof_per_edge_tot;
  fe_data->local_dof = (INT *) calloc(dof_per_edge_tot,sizeof(INT));
  fe_data->local_dof_flags = (INT *) calloc(dof_per_edge_tot,sizeof(INT));
  fe_data->u_local = (REAL *) calloc(dof_per_edge_tot,sizeof(REAL));

  // Don't construct basis functions on edges (yet...)

  return;
}
/****************************************************************************************/

/****************************************************************************************/
/*!
* \fn void get_edgelocaldata(simplex_local_data* edge_data,mesh_struct* mesh,INT edge)
*
* \brief On a given edge, grab the local mesh data.  We assume fixed values
*        for the mesh have been initialized already, and just filling in edge
*        specific data here.
*
* \param mesh     mesh data
* \param edge     edge we are considering
*
* \return edge_data Struct for local simplex data on edge
*
*/
void get_edgelocaldata(simplex_local_data* edge_data,mesh_struct* mesh,INT edge)
{
  // counters
  INT i,j,vind;

  INT dim = mesh->dim;
  INT v_per_ed = 2;

  // Generic elm data
  edge_data->sindex = edge;
  edge_data->vol = mesh->ed_len[edge];
  edge_data->flag = mesh->ed_flag[edge];

  // Get vertices on edge
  get_incidence_row(edge,mesh->ed_v,edge_data->local_v);
  // Temporary fix until coordinates is made dimensionless
  for(i=0;i<v_per_ed;i++) {
    vind = edge_data->local_v[i];
    for(j=0;j<dim;j++) {
      edge_data->xv[i*dim+j] = mesh->cv->x[vind*dim+j];
    }
  }

  // Edge data
  edge_data->ed_len[0] = mesh->ed_len[edge];
  for(j=0;j<dim;j++) {
    edge_data->ed_tau[j] = mesh->ed_tau[edge*dim+j];
    edge_data->ed_mid[j] = mesh->ed_mid[edge*dim+j];
  }

  // Quadrature (will ignore reference mapping for now)
  quad_edge_local(edge_data->quad_local,edge_data,edge_data->quad_local->nq1d,edge);

  return;
}
/******************************************************************************/

/*!
* \fn void free_blockfespace(block_fespace* FE)
*
* \brief Frees memory of arrays of block_fespace struct
*
* \return FE       Struct for BLOCK FE space to be freed
*
*/
void free_blockfespace(block_fespace* FE)
{
  if (FE == NULL) return; // Nothing need to be freed!

  INT i;
  INT num_spaces = (FE->nspaces);

  for ( i=0; i<num_spaces; i++ ) {
    free_fespace(FE->var_spaces[i]);
    FE->var_spaces[i] = NULL;
  }
  free(FE->var_spaces);
  FE->var_spaces = NULL;

  if(FE->dirichlet) free(FE->dirichlet);
  if(FE->dof_flag) free(FE->dof_flag);
  if(FE->simplex_data) {
    free_simplexlocaldata(FE->simplex_data);
    free(FE->simplex_data);
    FE->simplex_data=NULL;
  }
  if(FE->fe_data) {
    free_felocaldata(FE->fe_data);
    free(FE->fe_data);
    FE->fe_data=NULL;
  }
  return;
}

/*!
* \fn void free_felocaldata(fe_local_data* loc_data)
*
* \brief Frees memory of arrays of fe_local_data structure
*
* \return loc_data       Struct to be freed
*
*/
void free_felocaldata(fe_local_data* loc_data)
{
  if (loc_data == NULL) return; // Nothing need to be freed!

  if(loc_data->fe_types) free(loc_data->fe_types);
  if(loc_data->scal_or_vec) free(loc_data->scal_or_vec);
  if(loc_data->n_dof_per_space) free(loc_data->n_dof_per_space);
  if(loc_data->local_dof) free(loc_data->local_dof);
  if(loc_data->local_dof_flags) free(loc_data->local_dof_flags);
  if(loc_data->u_local) free(loc_data->u_local);

  // Actual basis functions are freed by free_fespace
  if(loc_data->phi) free(loc_data->phi);
  loc_data->phi = NULL;
  if(loc_data->dphi) free(loc_data->dphi);
  loc_data->dphi = NULL;
  if(loc_data->ddphi) free(loc_data->ddphi);
  loc_data->ddphi = NULL;

  return;
}

/*!
* \fn void free_simplexlocaldata(simplex_local_data* loc_data)
*
* \brief Frees memory of arrays of simplex_local_data structure
*
* \return loc_data       Struct to be freed
*
*/
void free_simplexlocaldata(simplex_local_data* loc_data)
{
  if (loc_data == NULL) return; // Nothing need to be freed!

  if(loc_data->local_v) free(loc_data->local_v);
  if(loc_data->xv) free(loc_data->xv);
  if(loc_data->local_ed) free(loc_data->local_ed);
  if(loc_data->v_on_ed) free(loc_data->v_on_ed);
  if(loc_data->ed_len) free(loc_data->ed_len);
  if(loc_data->ed_tau) free(loc_data->ed_tau);
  if(loc_data->ed_mid) free(loc_data->ed_mid);
  if(loc_data->local_f) free(loc_data->local_f);
  if(loc_data->v_on_f) free(loc_data->v_on_f);
  if(loc_data->f_area) free(loc_data->f_area);
  if(loc_data->f_norm) free(loc_data->f_norm);
  if(loc_data->f_mid) free(loc_data->f_mid);
  if(loc_data->ref_map) free(loc_data->ref_map);
  if(loc_data->lams) free(loc_data->lams);
  if(loc_data->gradlams) free(loc_data->gradlams);
  if(loc_data->dwork) free(loc_data->dwork);
  if(loc_data->iwork) free(loc_data->iwork);
  if(loc_data->quad_local){
    free_quadrature(loc_data->quad_local);
    free(loc_data->quad_local);
    loc_data->quad_local=NULL;
  }

  return;
}

/*!
* \fn void get_P2(fespace* FE,mesh_struct* mesh)
*
* \brief Converts mesh data to account for P2 elements
*
* \param mesh      Mesh struct
*
* \return FE       Struct for P2 FE space
*
*/
void get_P2(fespace* FE,mesh_struct* mesh)
{
  // Loop indices
  INT i,j,k,s,jcntr;
  INT n1,n2,n3,va,vb,ea,eb,v1,v2,mye,ed,na,ed1,ed2,face;

  INT ndof = FE->ndof;
  INT dof_per_elm = FE->dof_per_elm;
  INT dim = mesh->dim;
  INT nv = mesh->nv;
  INT nedge = mesh->nedge;
  INT nface = mesh->nface;
  INT nelm = mesh->nelm;
  INT v_per_elm = mesh->v_per_elm;
  iCSRmat *el_v = mesh->el_v;
  iCSRmat *el_ed = NULL;
  iCSRmat *ed_v = NULL;
  iCSRmat *f_v = NULL;
  iCSRmat *f_ed = NULL;

  if(dim>1) {
    el_ed = mesh->el_ed;
    ed_v = mesh->ed_v;
    f_v = mesh->f_v;
    f_ed = mesh->f_ed;
  }

  INT* ipv = (INT *) calloc(mesh->v_per_elm,sizeof(INT));

  // Get Coordinates
  coordinates *cdof = allocatecoords(ndof,dim);

  // First go through all vertices.
  for (i=0; i<nv; i++) {
    cdof->x[i] = mesh->cv->x[i];
    if(dim>1)
    cdof->y[i] = mesh->cv->y[i];
    if (dim>2)
    cdof->z[i] = mesh->cv->z[i];
  }
  // Now, go through and add extra nodes
  s = nv;

  // In 1D this is just the midpoint of the elements
  if(dim==1) {
    for(i=0;i<nelm;i++) {
      cdof->x[s] = mesh->el_mid[i];
      s++;
    }
  } else { // In 2D or 3D, these are simply the midpoints of the edges
    for (i=0; i<nedge; i++) {
      ed = ed_v->IA[i];
      n1 = ed_v->JA[ed];
      n2 = ed_v->JA[ed+1];
      cdof->x[s] = 0.5*(mesh->cv->x[n1]+mesh->cv->x[n2]);
      cdof->y[s] = 0.5*(mesh->cv->y[n1]+mesh->cv->y[n2]);
      if (dim>2) { cdof->z[s] = 0.5*(mesh->cv->z[n1]+mesh->cv->z[n2]); }
      s++;
    }
  }
  // Get Element to Node map
  iCSRmat el_n = icsr_create(nelm,ndof,dof_per_elm*nelm);
  // Rows of Element to Node map
  for(i=0;i<nelm+1;i++) {
    el_n.IA[i] = dof_per_elm*i;
  }
  // Columns
  // In 1D just add the midpoint of elements
  if(dim==1) {
    for(i=0;i<nelm;i++) {
      va = el_v->IA[i];
      vb = el_v->IA[i+1];
      na = el_n.IA[i];
      jcntr = 0;
      for(j=va;j<vb;j++) {
        el_n.JA[na + jcntr] = el_v->JA[j];
        jcntr++;
      }
      el_n.JA[na+jcntr] = nv+i;
    }
  } else {
    // For 2D or 3D we add the midpoint of the edges
    for(i=0;i<nelm;i++) {
      va = el_v->IA[i];
      vb = el_v->IA[i+1];
      ea = el_ed->IA[i];
      eb = el_ed->IA[i+1];
      na = el_n.IA[i];
      jcntr = 0;
      for(j=va;j<vb;j++) {
        ipv[jcntr] = el_v->JA[j];
        el_n.JA[na+jcntr] = ipv[jcntr];
        jcntr++;
      }
      for (j=0;j<v_per_elm;j++) {
        for (k=j+1;k<v_per_elm;k++) {
          v1 = ipv[j];
          v2 = ipv[k];
          for(s=ea;s<eb;s++) {
            mye = el_ed->JA[s];
            n1 = ed_v->JA[ed_v->IA[mye]];
            n2 = ed_v->JA[ed_v->IA[mye]+1];
            if ((n1==v1 && n2==v2) || (n2==v1 && n1==v2)) {
              el_n.JA[v_per_elm+na] = mye+nv;
              na++;
            }
          }
        }
      }
    }
  }

  // Fix Boundaries
  INT* dirichlet = (INT *) calloc(ndof,sizeof(INT));
  INT* dof_flag = (INT *) calloc(ndof,sizeof(INT));
  // First set of nodes are vertices
  for (i=0; i<nv; i++) {
    dirichlet[i] = mesh->v_flag[i];
    dof_flag[i] = mesh->v_flag[i];
  }
  // In 1D rest are interior
  if(dim==1) {
    for(i=0;i<nelm;i++) {
      dirichlet[nv+i] = 0;
      dof_flag[nv+i] = 0;
    }
  } else {
    // In 2D or 3D, rest are edges
    for(i=0;i<nedge;i++) {
      dirichlet[nv+i] = mesh->ed_flag[i];
      dof_flag[nv+i] = mesh->ed_flag[i];
    }
  }

  // Get edge to DOF map
  if(dim>1) {
    iCSRmat ed_n = icsr_create(nedge,ndof,3*nedge);
    s = nv;
    for (i=0; i<nedge; i++) {
      ed = ed_v->IA[i];
      ed_n.IA[i] = ed+i;
      n1 = ed_v->JA[ed];
      n2 = ed_v->JA[ed+1];
      ed_n.JA[ed+i] = n1;
      ed_n.JA[ed+1+i] = n2;
      ed_n.JA[ed+2+i] = s;
      s++;
    }
    ed_n.IA[nedge] = 3*nedge;

    // Get face to DOF map
    INT n_per_f = 3*(dim-1);
    INT extra_n = 2*dim-3;
    iCSRmat f_n = icsr_create(nface,ndof,n_per_f*nface);

    for (i=0; i<nface; i++) {
      face = f_v->IA[i];
      f_n.IA[i] = face+(2*dim-3)*i;
      n1 = f_v->JA[face];
      n2 = f_v->JA[face+1];
      if(dim==3) n3 = f_v->JA[face+2];
      f_n.JA[face+extra_n*i] = n1;
      f_n.JA[face+extra_n*i+1] = n2;
      if(dim==3) f_n.JA[face+extra_n*i+2] = n3;
      // Fill in rest with edge numbers
      ed1 = f_ed->IA[i];
      ed2 = f_ed->IA[i+1];
      jcntr = 0;
      for(j=ed1;j<ed2;j++) {
        f_n.JA[face+extra_n*i+dim+jcntr] = f_ed->JA[j]+nv;
        jcntr++;
      }
    }
    f_n.IA[nface] = n_per_f*nedge;
    *(FE->ed_dof) = ed_n;
    *(FE->f_dof) = f_n;
  }

  FE->dirichlet = dirichlet;
  FE->dof_flag = dof_flag;
  *(FE->el_dof) = el_n;
  FE->cdof = cdof;

  if(ipv) free(ipv);

  return;
}
/******************************************************************************************/

/*!
* \fn void get_MINI(fespace* FE,mesh_struct* mesh)
*
* \brief Converts mesh data to account for Mini elements
*
* \param mesh      Mesh struct
*
* \return FE       Struct for Mini FE space (P1 + cubic bubble)
*
*/
void get_MINI(fespace* FE,mesh_struct* mesh)
{
  // Loop indices
  INT i,j,s,jcntr;
  INT va,vb,dofa;

  INT ndof = FE->ndof;
  INT dof_per_elm = FE->dof_per_elm;
  INT dim = mesh->dim;
  INT nv = mesh->nv;
  INT nelm = mesh->nelm;
  iCSRmat *el_v = mesh->el_v;

  // Get Coordinates
  coordinates *cdof = allocatecoords(ndof,dim);
  // First go through all vertices.
  for (i=0; i<nv; i++) {
    for (j=0;j<dim;j++) {
      cdof->x[j*nv+i] = mesh->cv->x[j*nv+i];
    }
  }
  // Now, go through and add extra nodes for the midpoint of element
  s = nv;
  // In 1D this is just the midpoint of the elements
  for(i=0;i<nelm;i++) {
    for (j=0;j<dim;j++) {
      cdof->x[j*nv+s] = mesh->el_mid[i*dim+j];
    }
    s++;
  }

  // Get Element to DoF map
  iCSRmat el_dof = icsr_create(nelm,ndof,dof_per_elm*nelm);
  // Rows of Element to DOF map
  for(i=0;i<nelm+1;i++) {
    el_dof.IA[i] = dof_per_elm*i;
  }
  // Columns
  // Just add the midpoint of elements
  for(i=0;i<nelm;i++) {
    va = el_v->IA[i];
    vb = el_v->IA[i+1];
    dofa = el_dof.IA[i];
    jcntr = 0;
    for(j=va;j<vb;j++) {
      el_dof.JA[dofa + jcntr] = el_v->JA[j];
      jcntr++;
    }
    el_dof.JA[dofa+jcntr] = nv+i;
  }

  // Fix Boundaries
  INT* dirichlet = (INT *) calloc(ndof,sizeof(INT));
  INT* dof_flag = (INT *) calloc(ndof,sizeof(INT));
  // First set of nodes are vertices
  for (i=0; i<nv; i++) {
    dirichlet[i] = mesh->v_flag[i];
    dof_flag[i] = mesh->v_flag[i];
  }
  // Rest are interior element centers
  for(i=0;i<nelm;i++) {
    dirichlet[nv+i] = 0;
    dof_flag[nv+i] = mesh->el_flag[i];
  }

  FE->dirichlet = dirichlet;
  FE->dof_flag = dof_flag;
  *(FE->el_dof) = el_dof;
  FE->cdof = cdof;

  return;
}

/****************************************************************************************/
/*!
* \fn void dump_el_dof(FILE* fid,iCSRmat *el_dof)
*
* \brief Dump the element to DOF map to file for plotting purposes.
*
* \param el_dof      Element to DOF map
* \param fid         Output FILE ID
*
* \return el_dof.dat Output file with el_dof data
*
*/
void dump_el_dof(FILE* fid,iCSRmat *el_dof)
{
  // Loop indices
  INT i,j,acol,bcol;

  for (i=0; i<el_dof->row; i++) {
    acol = el_dof->IA[i];
    bcol = el_dof->IA[i+1];
    for (j=acol; j<bcol; j++) {
      fprintf(fid,"%lld\t",(long long )el_dof->JA[j]);
    }
    fprintf(fid,"\n");
  }

  return;
}
/****************************************************************************************/

/****************************************************************************************/
/*!
* \fn void dump_fespace(fespace *FE,char *varname,char *dir)
*
* \brief Dump the FE space data to file for plotting purposes
*
* \param FE               FE space to dump
* \param varname          Output file name
* \param dir              Directory to dump data
*
* \return dir/varname.dat Output file with FE data
*
*/
void dump_fespace(fespace *FE,char *varname,char *dir)
{
  INT i;
  INT totdof = FE->ndof;
  char eldofname[20];
  char bdryname[20];

  sprintf(eldofname,"%s/el_dof_%s.dat",dir,varname);
  sprintf(bdryname,"%s/bdry_%s.dat",dir,varname);
  FILE* fid1 = fopen(eldofname,"w");
  FILE* fid2 = fopen(bdryname,"w");

  if(fid1==NULL || fid2==NULL) {
    printf("\n\nFilenames: %s\t%s are incorrect or do not exist.  No files dumped.\n\n",eldofname,bdryname);
  } else {
    // Dump Element to DOF map
    dump_el_dof(fid1,FE->el_dof);

    // Dump boundary data
    for(i=0;i<totdof;i++) {
      fprintf(fid2,"%lld\n",(long long )FE->dirichlet[i]);
    }

    fclose(fid1);
    fclose(fid2);
  }
  return;
}
/****************************************************************************************/

/****************************************************************************************/
/*!
* \fn void set_dirichlet_bdry(fespace* FE,mesh_struct* mesh, const INT flag0, const INT flag1)
*
* \brief Determine which boundary DOF are Dirichlet.  Determined by the FE space type
*        and by the given flag from the mesh file.
*
* \param FE               FE space struct
* \param mesh             Mesh struct
* \param flag0            min flag value for Dirichlet DOF
*                         (e.g. in fem.h: MARKER_DIRICHLET)
* \param flag1            max flag value for Dirichlet DOF
*                         e.g. in fem.h (MARKER_NEUMANN - 1)
*
* \return FE.dirichlet    Binary boundary array for DOF
*
*/
void set_dirichlet_bdry(fespace* FE,mesh_struct* mesh, const INT flag0, const INT flag1)
{
  INT i;
  for(i=0;i<FE->ndof;i++) {
    if((FE->dof_flag[i])>=flag0 && (FE->dof_flag[i]<=flag1)) {
      FE->dirichlet[i] = 1;
    } else {
      FE->dirichlet[i] = 0;
    }
  }
  return;
}
/****************************************************************************************/

/****************************************************************************************/
/*!
* \fn void set_dirichlet_bdry_block(fespace* FE,mesh_struct* mesh)
*
* \brief Determine which boundary DOF are Dirichlet.  Determined by the BLOCK FE space type
*        and by the given flag from the mesh file.
*
* \param FE               BLOCK FE space struct
* \param mesh             Mesh struct
*
* \return FE.dirichlet    Binary boundary array for DOF
* \return FE.dof_flag     Also set DOF flags based on each FE space
*
*/
void set_dirichlet_bdry_block(block_fespace* FE,mesh_struct* mesh)
{
  INT i,j,ndof,cnt;

  INT* isdirichlet = (INT *) calloc(FE->ndof,sizeof(INT));
  INT* dof_flags = (INT *) calloc(FE->ndof,sizeof(INT));

  cnt = 0;
  for(i=0;i<FE->nspaces;i++) {
    ndof = FE->var_spaces[i]->ndof;
    for(j=0;j<ndof;j++) {
      isdirichlet[cnt+j] = FE->var_spaces[i]->dirichlet[j];
      dof_flags[cnt+j] = FE->var_spaces[i]->dof_flag[j];
    }
    cnt += ndof;
  }

  FE->dirichlet = isdirichlet;
  FE->dof_flag = dof_flags;
  return;
}
/****************************************************************************************/

/****************************************************************************************/
/*!
* \fn void set_periodic_bdry(fespace* FE,mesh_struct* mesh, const REAL minx,const REAL maxx,const REAL miny,const REAL maxy,const REAL minz,const REAL maxz)
*
* \brief Determine which boundary DOF are periodic, and for each DOF which is
*        the corresponding periodic DOF.
*
* \param FE                  FE space struct
* \param mesh                Mesh struct
* \param minx, miny, minz    Min values of x, y and z
* \param maxx, maxy, maxz    Max values of x, y and z
*
* \return FE.periodic     Array for each DOF, which contains the corresponding
*                         periodic DOF on the "other side".
*
* \note If Max and Min values are the same, the user indicates that that
*       dimension is not periodic.
*
*/
void set_periodic_bdry(fespace* FE,mesh_struct* mesh,const REAL minx,const REAL maxx,const REAL miny,const REAL maxy,const REAL minz,const REAL maxz)
{
  // Flag for errors
  SHORT status;

  // Loop indices
  INT i,j;

  // Indicator for periodicity
  INT xper=1,yper=1,zper=1;

  // Check which dimensions we are periodic in
  if(minx==maxx) xper=0;
  if(miny==maxy) yper=0;
  if(minz==maxz) zper=0;

  // Array to hold values of other dimensions
  REAL* xvals = (REAL *) calloc(FE->ndof,sizeof(REAL));
  REAL* yvals = (REAL *) calloc(FE->ndof,sizeof(REAL));
  REAL* zvals=NULL;
  if(mesh->dim==3) zvals = (REAL *) calloc(FE->ndof,sizeof(REAL));
  for(i=0;i<FE->ndof;i++) {
    xvals[i] = -666.66;
    yvals[i] = -666.66;
    if(mesh->dim==3) zvals[i] = -666.66;
  }

  // Tolerance for finding difference in real number values
  REAL mytol = 1e-13;

  // Cycle through DOF finding corresponding values.
  // This depends on the FE, but we assume FE->cdof contains the Coordinates
  if(mesh->dim==1) {
    for(i=0;i<FE->ndof;i++) {
      if(xper && fabs(FE->cdof->x[i]-minx)<=mytol) {
        for(j=0;j<FE->ndof;j++) {
          if(fabs(FE->cdof->x[j]-maxx)<=mytol) {
            FE->periodic[i]=j;
            FE->periodic[j]=i;
          }
        }
      }
    }
  } else if(mesh->dim==2) {
    for(i=0;i<FE->ndof;i++) {
      if(xper && fabs(FE->cdof->x[i]-minx)<=mytol) {
        yvals[i] = FE->cdof->y[i];
        for(j=0;j<FE->ndof;j++) {
          if(fabs(FE->cdof->x[j]-maxx)<=mytol && fabs(FE->cdof->y[j]-yvals[i])<=mytol) {
            FE->periodic[i]=j;
            FE->periodic[j]=i;
          }
        }
      }
      if(yper && fabs(FE->cdof->y[i]-miny)<=mytol) {
        xvals[i] = FE->cdof->x[i];
        for(j=0;j<FE->ndof;j++) {
          if(fabs(FE->cdof->y[j]-maxy)<=mytol && fabs(FE->cdof->x[j]-xvals[i])<=mytol) {
            FE->periodic[i]=j;
            FE->periodic[j]=i;
          }
        }
      }
    }
  } else if(mesh->dim==3) {
    for(i=0;i<FE->ndof;i++) {
      if(xper && fabs(FE->cdof->x[i]-minx)<=mytol) {
        yvals[i] = FE->cdof->y[i];
        zvals[i] = FE->cdof->z[i];
        for(j=0;j<FE->ndof;j++) {
          if(fabs(FE->cdof->x[j]-maxx)<=mytol && fabs(FE->cdof->y[j]-yvals[i])<=mytol && fabs(FE->cdof->z[j]-zvals[i])<=mytol) {
            FE->periodic[i]=j;
            FE->periodic[j]=i;
          }
        }
      }
      if(yper && fabs(FE->cdof->y[i]-miny)<=mytol) {
        xvals[i] = FE->cdof->x[i];
        zvals[i] = FE->cdof->z[i];
        for(j=0;j<FE->ndof;j++) {
          if(fabs(FE->cdof->y[j]-maxy)<=mytol && fabs(FE->cdof->x[j]-xvals[i])<=mytol && fabs(FE->cdof->z[j]-zvals[i])<=mytol) {
            FE->periodic[i]=j;
            FE->periodic[j]=i;
          }
        }
      }
      if(zper && fabs(FE->cdof->z[i]-minz)<=mytol) {
        xvals[i] = FE->cdof->x[i];
        yvals[i] = FE->cdof->y[i];
        for(j=0;j<FE->ndof;j++) {
          if(fabs(FE->cdof->z[j]-maxz)<=mytol && fabs(FE->cdof->x[j]-xvals[i])<=mytol && fabs(FE->cdof->y[j]-yvals[i])<=mytol) {
            FE->periodic[i]=j;
            FE->periodic[j]=i;
          }
        }
      }
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  if(xvals) free(xvals);
  if(yvals) free(yvals);
  if(zvals) free(zvals);
  return;
}
/****************************************************************************************/

/******************************************************************************/
/*!
* \fn void get_incidence_row(INT row,iCSRmat *fem_map,INT* thisrow)
*
* \brief Gets single row of an incidence map (i.e., gets vertices of given element from el_v)
*
* \param row       Row to grab (indexed at 0)
* \param fem_map   Incidence matrix to grab
*
* \return thisrow  Given row of the incidence matrix
*
*/
void get_incidence_row(INT row,iCSRmat *fem_map,INT* thisrow)
{
  INT j;
  INT rowa = fem_map->IA[row];
  INT rowb = fem_map->IA[row+1];
  INT jcntr = 0;
  for (j=rowa; j<rowb; j++) {
    thisrow[jcntr] = fem_map->JA[j];
    jcntr++;
  }

  return;
}
/******************************************************************************/

/******************************************************************************/
/*!
* \fn void get_coords(REAL* x,INT dof,coordinates* cdof,INT dim)
*
* \brief Gets coordinates of a particular DOF
*
* \param dof       DOF to grab (indexed at 0)
* \param cdof      Coordinate struct to grab from
* \param dim
*
* \return x        (x,y,z) coordinates
*
*/
void get_coords(REAL* x,INT dof,coordinates* cdof,INT dim)
{
  x[0] = cdof->x[dof];
  if(dim==2 || dim==3)
  x[1] = cdof->y[dof];
  if(dim==3)
  x[2] = cdof->z[dof];

  return;
}
/******************************************************************************/
