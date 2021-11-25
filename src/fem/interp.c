/*! \file src/fem/interp.c
 *
 *  \brief This code contains functions for interpolating FE approximations
 *         and their derivatives using FE basis functions, as well as routines
 *         for evaluating expressions on the FE spaces (projections)
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 2/1/15.
 *  Copyright 2016__HAZMATH__. All rights reserved.
 *
 *  \note modified by Casey Cavanaugh 04/18/2019
 *  \note modified by James Adler     04/19/2019
 *  \note Updated on 11/3/2018 for 0-1 fix.
 *  \note Updated on 03/31/2021 for efficiency
 *
 */

#include "hazmath.h"

/********************************************************************/
/**************** Interpolation Routines ****************************/
/********************************************************************/

/*!
 * \fn void fe_interpolation_to_x(REAL* val,REAL* u,REAL* x,fe_local_data* fe_data,simplex_local_data* elm_data,INT space_index)
 *
 * \brief Interpolate a finite-element approximation to any other point in the
 *        given element using the given type of elements.
 *
 * \note We assume u is given and properly ordered on the given element
 *
 * \param u           Values of FEM approximation at DoF of space on element only
 * \param x           Coordinates where to compute value
 * \param fe_data     Local FE space data
 * \param elm_data    Local Mesh Data
 * \param space_index Index of fe space inside fe system
 *
 * \return val        Pointer to value of approximation at given x
 *
 */
void fe_interpolation_to_x(REAL* val,REAL* u,REAL* x,fe_local_data* fe_data,simplex_local_data* elm_data,INT space_index)
{
  INT i,j;

  // Get FE and Mesh data
  INT dof_per_elm = fe_data->n_dof_per_space[space_index];
  INT scal_or_vec = fe_data->scal_or_vec[space_index];
  INT dim = elm_data->dim;

  // Spot to store data
  REAL coef[dim];

  // Grab basis functions on element at point x
  get_FEM_basis_at_x(elm_data,fe_data,space_index,x);
  REAL* phi = fe_data->phi[space_index];

  if(scal_or_vec==0) { // Scalar Element
    coef[0] = 0.0;
    for(j=0; j<dof_per_elm; j++) {
      coef[0] += u[j]*phi[j];
    }
    val[0] = coef[0];
  } else { // Vector Element
    for(i=0;i<dim;i++) {
      coef[i] = 0.0;
      for(j=0; j<dof_per_elm; j++) {
        coef[i] += u[j]*phi[j*dim+i];
      }
      val[i] = coef[i];
    }
  }

  return;
}
/******************************************************************************/

/*!
 * \fn void fe_interpolation_to_quadpt(REAL* val,fe_local_data* fe_data,simplex_local_data* elm_data,INT space_index,INT quadpt)
 *
 * \brief Interpolate a finite-element approximation to a quadrature point in the
 *        given element using the given type of elements.
 *
 * \note We assume u is contained in fe_data
 *
 * \param fe_data     Local FE space data
 * \param elm_data    Local Mesh Data
 * \param space_index Index of fe space inside fe system
 * \param quadpt      Index of quadrature point in elment to interpolate to
 *
 * \return val        Pointer to value of approximation at given quadpt
 *
 */
void fe_interpolation_to_quadpt(REAL* val,fe_local_data* fe_data,simplex_local_data* elm_data,INT space_index,INT quadpt)
{
  INT i,j;

  // Get FE and Mesh data
  INT dof_per_elm = fe_data->n_dof_per_space[space_index];
  INT scal_or_vec = fe_data->scal_or_vec[space_index];
  INT dim = elm_data->dim;
  REAL* u = fe_data->u_local;
  for(i=0;i<space_index;i++) u += fe_data->n_dof_per_space[space_index];

  // Spot to store data
  REAL coef[dim];

  // Grab basis functions on element at point x
  get_FEM_basis_at_quadpt(elm_data,fe_data,space_index,quadpt);
  REAL* phi = fe_data->phi[space_index];

  if(scal_or_vec==0) { // Scalar Element
    coef[0] = 0.0;
    for(j=0; j<dof_per_elm; j++) {
      coef[0] += u[j]*phi[j];
    }
    val[0] = coef[0];
  } else { // Vector Element
    for(i=0;i<dim;i++) {
      coef[i] = 0.0;
      for(j=0; j<dof_per_elm; j++) {
        coef[i] += u[j]*phi[j*dim+i];
      }
      val[i] = coef[i];
    }
  }

  return;
}
/******************************************************************************/

/*!
 * \fn void fe_dinterp_to_x(REAL* val,REAL* u,REAL* x,fe_local_data* fe_data,simplex_local_data* elm_data,INT space_index)
 *
 * \brief Interpolate the derivative of a finite-element approximation to any other point in the
 *        given element using the given type of elements.
 *
 * \note We assume u is given and properly ordered on the given element
 * \note Since the "type" of derivative depends on the FE type, we will have to
 *       check for various cases
 *
 *
 * \param u           Values of FEM approximation at DoF of space on element only
 * \param x           Coordinates where to compute value
 * \param fe_data     Local FE space data
 * \param elm_data    Local Mesh Data
 * \param space_index Index of fe space inside fe system
 *
 * \return val        Pointer to value of derivative of approximation at given x
 *
 */
 void fe_dinterp_to_x(REAL* val,REAL* u,REAL* x,fe_local_data* fe_data,simplex_local_data* elm_data,INT space_index)
 {
   INT i,j,k;

   // Get FE and Mesh data
   INT dof_per_elm = fe_data->n_dof_per_space[space_index];
   INT dim = elm_data->dim;
   INT fe_type = fe_data->fe_types[space_index];

   // Spot to store data
   REAL coef[dim*dim];

   // Grab basis functions on element at point x
   get_FEM_basis_at_x(elm_data,fe_data,space_index,x);
   REAL* dphi = fe_data->dphi[space_index];

   // Determine what type of space it is and interpolate appropriately
   INT grad_type=0;
   // Don't compute gradients for these Spaces
   if(fe_type==0 || fe_type==99) grad_type=0;
   // Derivatives are scalar-valued
   if((fe_type==20 && dim==2) || fe_type==30) grad_type=1;
   // Derivatives are vector-valued
   if((fe_type==20 && dim==3) || (fe_type>0 && fe_type<20)) grad_type=2;
   // Derivates are matrix-valued
   if(fe_type==60) grad_type=3;

   // Build interpolant
   if(grad_type==0) {
     // Don't compute derivative
     for(i=0;i<dim;i++) val[i] = 0.0;
   } else if(grad_type==1) {
     // Scalar-valued derivatives
     coef[0] = 0.0;
     for (j=0; j<dof_per_elm; j++) {
       coef[0] += u[j]*dphi[j];
     }
     val[0] = coef[0];
   } else if(grad_type==2) {
     // Vector-valued derivatives
     for(i=0;i<dim;i++) {
       coef[i] = 0.0;
       for(j=0; j<dof_per_elm; j++) {
         coef[i] += u[j]*dphi[j*dim+i];
       }
       val[i] = coef[i];
     }
   } else if(grad_type==3) {
     // Matrix-valued derivatives
     for(k=0;k<dim;k++) {
       for(i=0;i<dim;i++) {
         coef[k*dim + i] = 0.0;
         for(j=0; j<dof_per_elm; j++) {
           coef[k*dim + i] += u[j]*dphi[j*dim*dim + k*dim + i];
         }
         val[k*dim + i] = coef[k*dim + i];
       }
     }
   } else {
     check_error(ERROR_FE_TYPE,__FUNCTION__);
   }

   return;
 }
/******************************************************************************/

/*!
 * \fn void fe_dinterp_to_quadpt(REAL* val,fe_local_data* fe_data,simplex_local_data* elm_data,INT space_index,INT quadpt)
 *
 * \brief Interpolate the derivative of a finite-element approximation to a quadrature point in the
 *        given element using the given type of elements.
 *
 * \note Since the "type" of derivative depends on the FE type, we will have to
 *       check for various cases
 *
 * \note We assume u is contained in fe_data
 *
 * \param fe_data     Local FE space data
 * \param elm_data    Local Mesh Data
 * \param space_index Index of fe space inside fe system
 * \param quadpt      Index of quadrature point in elment to interpolate to
 *
 * \return val        Pointer to value of derivative of approximation at given quad point
 *
 */
 void fe_dinterp_to_quadpt(REAL* val,fe_local_data* fe_data,simplex_local_data* elm_data,INT space_index,INT quadpt)
 {
   INT i,j,k;

   // Get FE and Mesh data
   INT dof_per_elm = fe_data->n_dof_per_space[space_index];
   INT dim = elm_data->dim;
   INT fe_type = fe_data->fe_types[space_index];
   REAL* u = fe_data->u_local;
   for(i=0;i<space_index;i++) u += fe_data->n_dof_per_space[space_index];

   // Spot to store data
   REAL coef[dim*dim];

   // Grab basis functions on element at point x
   get_FEM_basis_at_quadpt(elm_data,fe_data,space_index,quadpt);
   REAL* dphi = fe_data->dphi[space_index];

   // Determine what type of space it is and interpolate appropriately
   INT grad_type=0;
   // Don't compute gradients for these Spaces
   if(fe_type==0 || fe_type==99) grad_type=0;
   // Derivatives are scalar-valued
   if((fe_type==20 && dim==2) || fe_type==30) grad_type=1;
   // Derivatives are vector-valued
   if((fe_type==20 && dim==3) || (fe_type>0 && fe_type<20)) grad_type=2;
   // Derivates are matrix-valued
   if(fe_type==60) grad_type=3;

   // Build interpolant
   if(grad_type==0) {
     // Don't compute derivative
     for(i=0;i<dim;i++) val[i] = 0.0;
   } else if(grad_type==1) {
     // Scalar-valued derivatives
     coef[0] = 0.0;
     for (j=0; j<dof_per_elm; j++) {
       coef[0] += u[j]*dphi[j];
     }
     val[0] = coef[0];
   } else if(grad_type==2) {
     // Vector-valued derivatives
     for(i=0;i<dim;i++) {
       coef[i] = 0.0;
       for(j=0; j<dof_per_elm; j++) {
         coef[i] += u[j]*dphi[j*dim+i];
       }
       val[i] = coef[i];
     }
   } else if(grad_type==3) {
     // Matrix-valued derivatives
     for(k=0;k<dim;k++) {
       for(i=0;i<dim;i++) {
         coef[k*dim + i] = 0.0;
         for(j=0; j<dof_per_elm; j++) {
           coef[k*dim + i] += u[j]*dphi[j*dim*dim + k*dim + i];
         }
         val[k*dim + i] = coef[k*dim + i];
       }
     }
   } else {
     check_error(ERROR_FE_TYPE,__FUNCTION__);
   }

   return;
 }
/******************************************************************************/

/*!
 * \fn void blockfe_interpolation_to_x(REAL* val,REAL* u,REAL* x,fe_local_data* fe_data,simplex_local_data* elm_data)
 *
 * \brief Interpolate a finite-element approximation to any other point in the
 *        given element for ALL FE spaces in system
 *
 * \note We assume u is given and properly ordered on the given element for each space
 *
 * \param u           Values of FEM approximation at DoF of all spaces on element only
 * \param x           Coordinates where to compute value
 * \param fe_data     Local FE space data
 * \param elm_data    Local Mesh Data
 *
 * \return val        Pointer to value of approximation at given x for all spaces
 *
 */
void blockfe_interpolation_to_x(REAL* val,REAL* u,REAL* x,fe_local_data* fe_data,simplex_local_data* elm_data)
{

  INT i;
  // Get FE and Mesh data
  INT nspaces = fe_data->nspaces;
  INT* dof_per_elm = fe_data->n_dof_per_space;
  INT* scal_or_vec = fe_data->scal_or_vec;
  INT dim = elm_data->dim;

  // Spot to store data per space
  REAL* val_space = val;
  REAL* u_space = u;

  // Loop over spaces and perform individual interpolation
  for(i=0;i<nspaces;i++) {
    fe_interpolation_to_x(val_space,u_space,x,fe_data,elm_data,i);
    if(scal_or_vec[i] == 0) {
      // Scalar quantity
      val_space++;
    } else {
      // Vector quantity
      val_space += dim;
    }
    u_space += dof_per_elm[i];
  }

  return;
}
/******************************************************************************/

/*!
 * \fn void blockfe_interpolation_to_quadpt(REAL* val,fe_local_data* fe_data,simplex_local_data* elm_data,INT quadpt)
 *
 * \brief Interpolate a finite-element approximation to a quadrature point in the
 *        given element using the given type of elements for ALL FE spaces in system
 *
 * \note We assume u is contained in fe_data for each space
 *
 * \param fe_data     Local FE space data
 * \param elm_data    Local Mesh Data
 * \param quadpt      Index of quadrature point in elment to interpolate to
 *
 * \return val        Pointer to value of approximation for all spaces at given quadpt
 *
 */
void blockfe_interpolation_to_quadpt(REAL* val,fe_local_data* fe_data,simplex_local_data* elm_data,INT quadpt)
{

  INT i;

  // Get FE and Mesh data
  INT nspaces = fe_data->nspaces;
  INT* scal_or_vec = fe_data->scal_or_vec;
  INT dim = elm_data->dim;

  // Spot to store data per space
  REAL* val_space = val;

  // Loop over spaces and perform individual interpolation
  for(i=0;i<nspaces;i++) {
    fe_interpolation_to_quadpt(val_space,fe_data,elm_data,i,quadpt);
    if(scal_or_vec[i] == 0) {
      // Scalar quantity
      val_space++;
    } else {
      // Vector quantity
      val_space += dim;
    }
  }

  return;
}
/******************************************************************************/

/*!
 * \fn void blockfe_dinterp_to_x(REAL* val,REAL* u,REAL* x,fe_local_data* fe_data,simplex_local_data* elm_data)
 *
 * \brief Interpolate derivative of a finite-element approximation to any other point in the
 *        given element for ALL FE spaces in system
 *
 * \note We assume u is given and properly ordered on the given element for each space
 *
 * \param u           Values of FEM approximation at DoF of all spaces on element only
 * \param x           Coordinates where to compute value
 * \param fe_data     Local FE space data
 * \param elm_data    Local Mesh Data
 *
 * \return val        Pointer to value of derivative of approximation at given x for all spaces
 *
 */
 void blockfe_dinterp_to_x(REAL* val,REAL* u,REAL* x,fe_local_data* fe_data,simplex_local_data* elm_data)
 {
   INT i;

   // Get FE and Mesh data
   INT nspaces = fe_data->nspaces;
   INT fe_type;
   INT* dof_per_elm = fe_data->n_dof_per_space;
   INT dim = elm_data->dim;

   // Spot to store data per space
   REAL* val_space = val;
   REAL* u_space = u;

   // Loop over spaces and call individual interpolation routines
   for(i=0;i<nspaces;i++) {
     fe_type = fe_data->fe_types[i];
     fe_dinterp_to_x(val_space,u_space,x,fe_data,elm_data,i);
     // Derivatives are scalar-valued
     if((fe_type==20 && dim==2) || fe_type==30) val_space++;
     // Derivatives are vector-valued
     if((fe_type==20 && dim==3) || (fe_type>0 && fe_type<20)) val_space += dim;
     // Derivates are matrix-valued
     if(fe_type==60) val_space += dim*dim;

     // Update approximation pointer
     u_space += dof_per_elm[i];
   }

   return;
 }
/******************************************************************************/

/*!
 * \fn void blockfe_dinterp_to_quadpt(REAL* val,fe_local_data* fe_data,simplex_local_data* elm_data,INT quadpt)
 *
 * \brief Interpolate derivative of a finite-element approximation to a quadrature point in the
 *        given element using the given type of elements for ALL FE spaces in system
 *
 * \note We assume u is contained in fe_data for each space
 *
 * \param fe_data     Local FE space data
 * \param elm_data    Local Mesh Data
 * \param quadpt      Index of quadrature point in elment to interpolate to
 *
 * \return val        Pointer to derivative of approximation for all spaces at given quadpt
 *
 */
void blockfe_dinterp_to_quadpt(REAL* val,fe_local_data* fe_data,simplex_local_data* elm_data,INT quadpt)
{

  INT i;

  // Get FE and Mesh data
  INT fe_type;
  INT nspaces = fe_data->nspaces;
  INT dim = elm_data->dim;

  // Spot to store data per space
  REAL* val_space = val;

  // Loop over spaces and call individual interpolation routines
  for(i=0;i<nspaces;i++) {
    fe_type = fe_data->fe_types[i];
    fe_dinterp_to_quadpt(val_space,fe_data,elm_data,i,quadpt);
    // Derivatives are scalar-valued
    if((fe_type==20 && dim==2) || fe_type==30) val_space++;
    // Derivatives are vector-valued
    if((fe_type==20 && dim==3) || (fe_type>0 && fe_type<20)) val_space += dim;
    // Derivates are matrix-valued
    if(fe_type==60) val_space += dim*dim;
  }

  return;
}
/******************************************************************************/

/**************** Projection Routines *******************************/
/********************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn REAL FE_Evaluate_DOF(void (*expr)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,REAL time,INT DOF)
 *
 * \brief Evaluate a given analytical function on the specific degree of freedom of the finite-element space given
 *
 * \param expr 	    Function call to analytical expression
 * \param FE        FE Space
 * \param mesh      Mesh Data
 * \param time      Physical time to evaluate function at
 * \param DOF       DOF index to evaluate (start at 0)
 *
 * \return val      FE approximation of function on fespace
 *
 */
REAL FE_Evaluate_DOF(void (*expr)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,REAL time,INT DOF)
{
  INT j,m;
  INT maxdim=4; // Used to allocate coordinate arrays
  REAL x[maxdim];
  INT dim = mesh->dim;
  INT FEtype = FE->FEtype;
  INT nq1d = 3; // If quadrature not given, fix order.
  REAL val=-666e+00;
  REAL* valx=NULL;

  INT* face_vertex = (INT *) calloc(dim,sizeof(INT));

  // P0 elements u[dof] = 1/elvol \int_el u
  if(FEtype==0) {
    val = (1.0/mesh->el_vol[DOF])*integrate_elm(expr,1,0,nq1d,NULL,mesh,time,DOF);

  // Lagrange Elements u[dof] = u[x_i]
  } else if(FEtype>0 && FEtype<10) {
    x[0] = FE->cdof->x[DOF];
    if(dim==2 || dim==3)
      x[1] = FE->cdof->y[DOF];
    if(dim==3)
      x[2] = FE->cdof->z[DOF];
    (*expr)(&val,x,time,&(FE->dof_flag[DOF]));

  // Nedelec u[dof] = (1/elen) \int_edge u*t_edge
  } else if (FEtype==20) {
    val = (1.0/mesh->ed_len[DOF])*integrate_edge_vector_tangent(expr,dim,0,nq1d,NULL,mesh,time,DOF);

  // Raviart-Thomas u[dof] = 1/farea \int_face u*n_face
  } else if (FEtype==30) {
    val = (1.0/mesh->f_area[DOF])*integrate_face_vector_normal(expr,dim,0,nq1d,NULL,mesh,time,DOF);

  // Face Bubbles
  // For now, we assume face bubbles only and only in 2D or 3D.
  // Here we compute the integral over the face, but then subtract off the
  // linear components:
  // u[dof] = (int_f v*n_f - |f|/dim * sum(vertices on face) v(i)*n_f) / int_f phi_f*n_f
  // Note that since phi_f = 2^d lam_i lam_j ..., we have
  // int_f phi_f*n_f = 2^d (d-1)! / (2d-1)! |f|
  } else if (FEtype>=60 && FEtype<70) {
    valx = (REAL *) calloc(dim,sizeof(REAL));
    val = integrate_face_vector_normal(expr,dim,0,nq1d,NULL,mesh,time,DOF);
    get_incidence_row(DOF,mesh->f_v,face_vertex);
    for (m=0;m<dim;m++) {
      x[0] = mesh->cv->x[face_vertex[m]];
      x[1] = mesh->cv->y[face_vertex[m]];
      if(dim==3) x[2] = mesh->cv->z[face_vertex[m]];
      // The following only works for 2D
      (*expr)(valx,x,time,&(FE->dof_flag[DOF]));
      for(j=0;j<dim;j++) val+= -(mesh->f_area[DOF]/dim)*mesh->f_norm[DOF*dim+j]*valx[j];
    }
    // Divide by int_f \phi_f*n_f
    if(dim==2) { // In 2D, int_f \phi_f*n_f = 2/3 |f|
      val = 1.5*val/(mesh->f_area[DOF]);
    } else if(dim==3) { // In 3D, int_f \phi_f*n_f = 2/15 |f|
      val = 7.5*val/(mesh->f_area[DOF]);
    } else {
      check_error(ERROR_DIM,__FUNCTION__);
    }
  // Constraint Element just grab constant value
  } else if(FEtype==99) {
    x[0] =0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    (*expr)(&val,x,time,&(FE->dof_flag[DOF]));
  // No other FEM types implemented
  } else {
    check_error(ERROR_FE_TYPE,__FUNCTION__);
  }

  if (face_vertex) free(face_vertex);
  //if (valx) free(valx);

  return val;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void FE_Evaluate(REAL* val,void (*expr)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,REAL time)
 *
 * \brief Evaluate (Project) a given analytical function on the finite-element space given.
 *
 * \param expr 	      Function call to analytical expression
 * \param FE          FE Space
 * \param mesh        Mesh Data
 * \param time        Physical time to evaluate function at
 *
 * \return val        FE approximation of function on fespace
 *
 * \note Just calls FE_Evalue_DOF in a loop over all DOF
 *
 */
void FE_Evaluate(REAL* val,void (*expr)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,REAL time)
{
  INT i;

  for(i=0;i<FE->ndof;i++) {
    val[i] = FE_Evaluate_DOF(expr,FE,mesh,time,i);
  }

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn REAL blockFE_Evaluate_DOF(void (*expr)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,REAL time,INT comp,INT DOF)
 *
 * \brief Evaluate a given analytical function on the specific degree of freedom of the block finite-element space given
 *
 * \param expr 	    Function call to analytical expression
 * \param FE        FE Space
 * \param mesh      Mesh Data
 * \param time      Physical time to evaluate function at
 * \param comp      FE block to consider (indexes at 0)
 * \param DOF       DOF index to evaluate (start at 0)
 *
 * \return val         FE approximation of function on fespace
 *
 * \note For all "integrals" we use the midpoint rule
 *
 */
REAL blockFE_Evaluate_DOF(void (*expr)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,REAL time,INT comp,INT DOF)
{
  INT i,j,m;
  INT nq1d = 3; // Number of quadrature points in 1D if not given.
  INT dim = mesh->dim;
  INT maxdim=4;
  REAL x[maxdim];
  REAL valx[dim*FE->nspaces];
  INT local_dim = 0;
  REAL val=-666e+00;
  INT face_vertex[maxdim];

  for(i=0;i<comp;i++) {
    if((FE->var_spaces[i]->FEtype>=0 && FE->var_spaces[i]->FEtype<10) || FE->var_spaces[i]->FEtype==99) { // Scalar Element
      local_dim += 1;
    } else if(FE->var_spaces[i]->FEtype == 61) { // Bubble Element
      local_dim += 0; // Assuming bubbles are written such that linears follow it immediately
    } else { // Vector Element
      local_dim += dim;
    }
  }

  // P0 elements u[dof] = 1/elvol \int_el u
  if(FE->var_spaces[comp]->FEtype==0) {
    val = (1.0/mesh->el_vol[DOF])*integrate_elm(expr,FE->nun,local_dim,nq1d,NULL,mesh,time,DOF);

  // Lagrange Elements u[dof] = u[x_i]
  } else if(FE->var_spaces[comp]->FEtype>0 && FE->var_spaces[comp]->FEtype<10) {
    x[0] = FE->var_spaces[comp]->cdof->x[DOF];
    if(dim==2 || dim==3)
      x[1] = FE->var_spaces[comp]->cdof->y[DOF];
    if(dim==3)
      x[2] = FE->var_spaces[comp]->cdof->z[DOF];
    (*expr)(valx,x,time,&(FE->var_spaces[comp]->dof_flag[DOF]));
    val = valx[local_dim];

  // Nedelec u[dof] = (1/elen) \int_edge u*t_edge
  } else if (FE->var_spaces[comp]->FEtype==20) {
    val = (1.0/mesh->ed_len[DOF])*integrate_edge_vector_tangent(expr,FE->nun,local_dim,nq1d,NULL,mesh,time,DOF);

  // Raviart-Thomas u[dof] = 1/farea \int_face u*n_face
  } else if (FE->var_spaces[comp]->FEtype==30) {
    val = (1.0/mesh->f_area[DOF])*integrate_face_vector_normal(expr,FE->nun,local_dim,nq1d,NULL,mesh,time,DOF);

  // Face Bubbles
  // For now, we assume face bubbles only and only in 2D or 3D
  // Here we compute the integral over the face, but then subtract off the
  // linear components:
  // u[dof] = (int_f v*n_f - |f|/dim * sum(vertices on face) v(i)*n_f) / int_f phi_f*n_f
  // Note that since phi_f = 2^d lam_i lam_j ..., we have
  // int_f phi_f*n_f = 2^d (d-1)! / (2d-1)! |f|
  } else if (FE->var_spaces[comp]->FEtype>=60 && FE->var_spaces[comp]->FEtype<70) {
    val = integrate_face_vector_normal(expr,FE->nun,local_dim,nq1d,NULL,mesh,time,DOF);
    get_incidence_row(DOF,mesh->f_v,face_vertex);
    for (m=0;m<dim;m++) {
      x[0] = mesh->cv->x[face_vertex[m]];
      x[1] = mesh->cv->y[face_vertex[m]];
      if(dim==3) x[2] = mesh->cv->z[face_vertex[m]];
      (*expr)(valx,x,time,&(FE->var_spaces[comp]->dof_flag[DOF]));
      for(j=0;j<dim;j++) val+= -(mesh->f_area[DOF]/dim)*mesh->f_norm[DOF*dim+j]*valx[local_dim + j];
    }
    // Divide by int_f \phi_f*n_f
    if(dim==2) { // In 2D, int_f \phi_f*n_f = 2/3 |f|
      val = 1.5*val/(mesh->f_area[DOF]);
    } else if(dim==3) { // In 3D, int_f \phi_f*n_f = 2/15 |f|
      val = 7.5*val/(mesh->f_area[DOF]);
    } else {
      check_error(ERROR_DIM,__FUNCTION__);
    }
    // Constraint Element just grab constant value
  } else if(FE->var_spaces[comp]->FEtype==99) {
    x[0] =0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    (*expr)(valx,x,time,&(FE->var_spaces[comp]->dof_flag[DOF]));
    val = valx[local_dim];
    // Not a FEM implemented
  } else {
    check_error(ERROR_FE_TYPE,__FUNCTION__);
  }

  return val;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void blockFE_Evaluate(REAL* val,void (*expr)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,REAL time)
 *
 * \brief Evaluate (Project) a given analytical function on the block finite-element space given.
 *
 * \param expr 	    Function call to analytical expression
 * \param FE        block FE Space for multiple variables
 * \param mesh      Mesh Data
 * \param time      Physical time to evaluate function at
 * \param val       FE approximation of function on fespace
 *
 */
void blockFE_Evaluate(REAL* val,void (*expr)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,REAL time)
{
  int i,k;
  INT entry = 0;

  for(k=0;k<FE->nspaces;k++) {
    for(i=0;i<FE->var_spaces[k]->ndof;i++) {
      val[entry] = blockFE_Evaluate_DOF(expr,FE,mesh,time,k,i);
      entry++;
    }
  }

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void blockFE_EvaluateBoundary(REAL* val,void (*expr)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,REAL time)
 *
 * \brief Evaluate (Project) a given analytical function on the block
 *        finite-element space given only at Dirichlet DoF
 *
 * \param expr 	    Function call to analytical expression
 * \param FE        block FE Space for multiple variables
 * \param mesh      Mesh Data
 * \param time      Physical time to evaluate function at
 * \param val       FE approximation of function on fespace
 *
 */
void blockFE_EvaluateBoundary(REAL* val,void (*expr)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,REAL time)
{
  int i,k;
  INT entry = 0;
  INT nspaces = FE->nspaces;
  INT ndof = 0;

  for(k=0;k<nspaces;k++) {
    ndof = FE->var_spaces[k]->ndof;
    for(i=0;i<ndof;i++) {
      if(FE->var_spaces[k]->dirichlet[i]==1) {
        val[entry] = blockFE_Evaluate_DOF(expr,FE,mesh,time,k,i);
      }
      entry++;
    }
  }

  return;
}
/****************************************************************************************************************************/

/***********************************************************************************************/
/*!
 * \fn void Project_to_Vertices(REAL* u_on_V,REAL *u,fespace *FE,mesh_struct *mesh)
 *
 * \brief Interpolate a finite-element approximation to the vertices of a mesh
 *
 * \param u_on_V  Solution on vertices of mesh
 * \param u       Approximation to Interpolate
 * \param FE      FE space
 * \param mesh    Mesh Data
 *
 */
void Project_to_Vertices(REAL* u_on_V,REAL *u,fespace *FE,mesh_struct *mesh)
{
  INT i,k,j;
  INT dim = mesh->dim;
  INT maxdim = 4;
  REAL* x = (REAL *) calloc(maxdim,sizeof(REAL));
  REAL val[maxdim];

  // Get FE and Mesh data
  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  //INT FEtype = FE->FEtype;
  INT scal_or_vec = FE->scal_or_vec;
  INT nelm = mesh->nelm;
  INT nv = mesh->nv;

  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  // Loop over Elements
  for(i=0;i<nelm;i++) {

    // Find DOF for given Element
    get_incidence_row(i,FE->el_dof,dof_on_elm);

    // Find vertices for given Element
    get_incidence_row(i,mesh->el_v,v_on_elm);

    // Interpolate FE approximation to vertices
    for(j=0;j<v_per_elm;j++) {
      get_coords(x,v_on_elm[j],mesh->cv,dim);
      FE_Interpolation(val,u,x,dof_on_elm,v_on_elm,FE,mesh);

      if(scal_or_vec==0) { // Scalar Element
        u_on_V[v_on_elm[j]] = val[0];
      } else { // Vector Element
        for(k=0;k<dim;k++) {
          u_on_V[k*nv + v_on_elm[j]] = val[k];
        }
      }
    }
  }

  if(v_on_elm) free(v_on_elm);
  if(dof_on_elm) free(dof_on_elm);
  if(x) free(x);
  return;
}
/***********************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void get_unknown_component(dvector* u,dvector* ublock,block_fespace *FE,INT comp)
 *
 * \brief Grab a component of a block FE function.
 *
 * \param u	         FE approximation of function on specific block of fespace
 * \param ublock       FE vector in blocks
 * \param FE           block FE Space struct for multiple variables
 * \param comp         Which block to grab (starts count at 0)
 *
 */
void get_unknown_component(dvector* u,dvector* ublock,block_fespace *FE,INT comp)
{
  INT i=0;
  INT entry = 0;

  for(i=0;i<comp;i++) {
    entry += FE->var_spaces[i]->ndof;
  }

  for(i=0;i<FE->var_spaces[comp]->ndof;i++) {
    u->val[i] = ublock->val[entry + i];
  }

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void set_unknown_component(dvector* u,dvector* ublock,block_fespace *FE,INT comp)
 *
 * \brief Set a component of a block FE function
 *
 * \param u	         FE approximation of function on specific block of fespace
 * \param ublock       FE vector in blocks
 * \param FE           block FE Space struct for multiple variables
 * \param comp         Which block to grab (starts count at 0)
 *
 */
void set_unknown_component(dvector* u,dvector* ublock,block_fespace *FE,INT comp)
{
  INT i=0;
  INT entry = 0;

  for(i=0;i<comp;i++) {
    entry += FE->var_spaces[i]->ndof;
  }

  for(i=0;i<FE->var_spaces[comp]->ndof;i++) {
    ublock->val[entry + i] = u->val[i];
  }

  return;
}
/****************************************************************************************************************************/

/***********************************************************************************************/
/*!
 * \fn void get_grad_H1toNed(dCSRmat* Grad,mesh_struct* mesh)
 *
 * \brief Computes Gradient operator.
 *        Applying the resulting matrix computes the gradient of an H1 approximation.
 *        Takes Nodal (H1) DOF vector to Edge (Nedelec) DOF vector
 *        Ordering determined by edge to node map: bigger node is +1 smaller is -1
 *
 * \param Grad	     dCSRmat matrix that takes the gradient of an H1 approximation
 * \param mesh       Mesh data
 *
 * \note This only makes sense for P1 elements and for dimensions 2 or 3.
 *
 */
void get_grad_H1toNed(dCSRmat* Grad,mesh_struct* mesh)
{
  INT i,j,k,rowa;
  INT nedge = mesh->nedge;
  REAL oneoverlen;
  dCSRmat Gtmp;

  Gtmp.row = mesh->ed_v->row;
  Gtmp.col = mesh->ed_v->col;
  Gtmp.nnz = mesh->ed_v->nnz;
  Gtmp.IA = (INT *) calloc(Gtmp.row+1,sizeof(INT));
  Gtmp.JA = (INT *) calloc(Gtmp.nnz,sizeof(INT));
  Gtmp.val = (REAL *) calloc(Gtmp.nnz,sizeof(REAL));

  for(i=0;i<=nedge;i++) {
    Gtmp.IA[i] = mesh->ed_v->IA[i];
  }
  for(i=0;i<Gtmp.nnz;i++) {
    Gtmp.JA[i] = mesh->ed_v->JA[i];
  }

  for (i=0; i<nedge; i++) {
    oneoverlen = 1.0/(mesh->ed_len[i]);
    rowa = mesh->ed_v->IA[i];
    j = mesh->ed_v->JA[rowa];
    k = mesh->ed_v->JA[rowa+1];
    if(j>k) {
      Gtmp.val[rowa] = oneoverlen;
      Gtmp.val[rowa+1] = -oneoverlen;
    } else {
      Gtmp.val[rowa+1] = oneoverlen;
      Gtmp.val[rowa] = -oneoverlen;
    }
  }

  *Grad = Gtmp;

  return;
}
/***********************************************************************************************/

/***********************************************************************************************/
/*!
 * \fn void get_curl_NedtoRT(dCSRmat* Curl,mesh_struct* mesh)
 *
 * \brief Computes Curl operator in 3D ONLY.
 *        Applying the resulting matrix computes the Curl of a Nedelec approximation.
 *        Takes Edge (Nedelec) DOF vector to Face (RT) DOF vector
 *
 * \param Curl 	     dCSRmat matrix that takes the curl of a Nedelec approximation
 * \param mesh         Mesh data
 *
 * \note This only makes sense in 3D and assuming lowest order elements
 *
 */
void get_curl_NedtoRT(dCSRmat* Curl,mesh_struct* mesh)
{
  INT i,j,k,s,col_a,nd1,nd2,nd3,rowa,rowb,jcntr;
  INT ndpf = mesh->dim;
  INT mydim = mesh->dim;
  INT nface = mesh->nface;
  INT* inf = (INT *) calloc(ndpf,sizeof(INT));
  REAL vec1[3],vec2[3],vec3[3];
  REAL mydet;

  dCSRmat Ktmp;

  Ktmp.row = mesh->f_ed->row;
  Ktmp.col = mesh->f_ed->col;
  Ktmp.nnz = mesh->f_ed->nnz;
  Ktmp.IA = (INT *) calloc(Ktmp.row+1,sizeof(INT));
  Ktmp.JA = (INT *) calloc(Ktmp.nnz,sizeof(INT));
  Ktmp.val = (REAL *) calloc(Ktmp.nnz,sizeof(REAL));

  for(i=0;i<=Ktmp.row;i++) {
    Ktmp.IA[i] = mesh->f_ed->IA[i];
  }
  for(i=0;i<Ktmp.nnz;i++) {
    Ktmp.JA[i] = mesh->f_ed->JA[i];
  }

  // Get Kcurl -> if_ed,jf_ed, sign = sign(det(xi-xk,xj-xk,n_fijk))
  for (i=0;i<nface;i++) {
    // Get Normal Vector
    vec3[0] = mesh->f_norm[i*mydim];
    vec3[1] = mesh->f_norm[i*mydim+1];
    vec3[2] = mesh->f_norm[i*mydim+2];
    // Get nodes of Face
    rowa = mesh->f_v->IA[i];
    rowb = mesh->f_v->IA[i+1];
    jcntr=0;
    for(j=rowa;j<rowb;j++) {
      inf[jcntr] = mesh->f_v->JA[j];
      jcntr++;
    }
    // Get edges of face
    rowa = mesh->f_ed->IA[i];
    rowb = mesh->f_ed->IA[i+1];
    for(j=rowa;j<rowb;j++) {
      k = mesh->f_ed->JA[j];
      // Get nodes of edge
      col_a = mesh->ed_v->IA[k];
      nd1 = mesh->ed_v->JA[col_a];
      nd2 = mesh->ed_v->JA[col_a+1];
      // Determine what other node on face is
      for(s=0;s<ndpf;s++) {
        if(inf[s]!=nd1 && inf[s]!=nd2) {
          nd3 = inf[s];
        }
      }
      vec1[0] = mesh->cv->x[nd1]-mesh->cv->x[nd3];
      vec2[0] = mesh->cv->x[nd2]-mesh->cv->x[nd3];
      vec1[1] = mesh->cv->y[nd1]-mesh->cv->y[nd3];
      vec2[1] = mesh->cv->y[nd2]-mesh->cv->y[nd3];
      vec1[2] = mesh->cv->z[nd1]-mesh->cv->z[nd3];
      vec2[2] = mesh->cv->z[nd2]-mesh->cv->z[nd3];
      if(nd1>nd2) {
        det3D(&mydet,vec2,vec1,vec3);
      } else {
        det3D(&mydet,vec1,vec2,vec3);
      }
      if(mydet>0) {
        Ktmp.val[j]=1;
      } else {
        Ktmp.val[j]=-1;
      }
    }
  }

  // K -> Df^(-1) K De
  for(i=0;i<nface;i++) {
    rowa = mesh->f_ed->IA[i];
    rowb = mesh->f_ed->IA[i+1];
    for(j=rowa;j<rowb;j++) {
      k = mesh->f_ed->JA[j];
      Ktmp.val[j] = (1.0/(mesh->f_area[i]))*( (REAL) Ktmp.val[j])*(mesh->ed_len[k]);
    }
  }

  *Curl = Ktmp;

  if(inf) free(inf);

  return;
}
/***********************************************************************************************/

/***********************************************************************************************/
/*!
 * \fn void get_div_RTtoL2(dCSRmat* Div,mesh_struct* mesh)
 *
 * \brief Computes Divergence operator.
 *        Applying the resulting matrix computes the divergence of an RT approximation.
 *        Takes Face (RT) DOF vector to Volume (L2) DOF vector.
 *        Ordering determined by element to face map: orientation determined by normal vector of face.
 *        Normal vectors point from lower number element to higher one or outward from external boundary
 *        Entry gets positive 1 if normal is outward of element and -1 if inward of element
 *
 * \param Div 	     dCSRmat matrix that takes the div of a RT approximation
 * \param mesh         Mesh data
 *
 * \note This only makes sense for lowest order elements...I think...and 2D or 3D
 *
 */
void get_div_RTtoL2(dCSRmat* Div,mesh_struct* mesh)
{
  INT i,j,rowa,rowb,rowc,rowd,face,elm1,elm2,elm_big,n_felm;
  INT nelm = mesh->nelm;
  REAL oneovervol=0.0;
  REAL farea=0.0;
  dCSRmat Dtmp;

  // We will need the face to element map as well
  // Each face will have two elements (or one if on boundary).
  iCSRmat f_el;
  icsr_trans(mesh->el_f,&f_el);

  Dtmp.row = mesh->el_f->row;
  Dtmp.col = mesh->el_f->col;
  Dtmp.nnz = mesh->el_f->nnz;
  Dtmp.IA = (INT *) calloc(Dtmp.row+1,sizeof(INT));
  Dtmp.JA = (INT *) calloc(Dtmp.nnz,sizeof(INT));
  Dtmp.val = (REAL *) calloc(Dtmp.nnz,sizeof(REAL));

  for(i=0;i<=nelm;i++) {
    Dtmp.IA[i] = mesh->el_f->IA[i];
  }
  for(i=0;i<Dtmp.nnz;i++) {
    Dtmp.JA[i] = mesh->el_f->JA[i];
  }

  for (i=0; i<nelm; i++) {
    oneovervol = 1.0/(mesh->el_vol[i]);
    // Get faces of element
    rowa = mesh->el_f->IA[i];
    rowb = mesh->el_f->IA[i+1];
    for(j=rowa;j<rowb;j++) {
      face = mesh->el_f->JA[j];
      farea = mesh->f_area[face];
      // Get elements of face
      rowc = f_el.IA[face];
      rowd = f_el.IA[face+1];
      n_felm = rowd-rowc;
      if(n_felm==1)
        Dtmp.val[j] = farea*oneovervol;
      else if(n_felm==2) {
        elm1 = f_el.JA[rowc];
        elm2 = f_el.JA[rowc+1];
        elm_big = MAX(elm1,elm2);
        if(i==elm_big)
          Dtmp.val[j] = -farea*oneovervol;
        else
          Dtmp.val[j] = farea*oneovervol;
      } else {
        printf("something is wrong in the divergence operator\n");
        exit(-50);
      }
    }
  }

  *Div = Dtmp;

  return;
}
/***********************************************************************************************/

/***********************************************************************************************/
/*!
 * \fn void get_Pigrad_H1toNed(dCSRmat* Pgrad,mesh_struct* mesh)
 *
 * \brief Computes Gradient operator into scalar components for HX preconditioner.
 *        Ordering determined by edge to node map: bigger node is +1 smaller is -1
 *
 * \param Pgrad 	     dCSRmat matrix that takes the the \Pi for HX preconditioner.
 * \param mesh         Mesh data
 *
 * \note This only makes sense in 2D or 3D and P1 elements.
 *       Also assumes shuffled ordering.
 *
 */
void get_Pigrad_H1toNed(dCSRmat* Pgrad,mesh_struct* mesh)
{
  INT i,j,k,rowa,cola;
  INT nedge = mesh->nedge;
  INT dim = mesh->dim;
  REAL oneoverlen,xL,yL,zL;
  dCSRmat Ptmp;

  Ptmp.row = mesh->ed_v->row;
  Ptmp.col = (mesh->ed_v->col)*dim;
  Ptmp.nnz = (mesh->ed_v->nnz)*dim;
  Ptmp.IA = (INT *) calloc(Ptmp.row+1,sizeof(INT));
  Ptmp.JA = (INT *) calloc(Ptmp.nnz,sizeof(INT));
  Ptmp.val = (REAL *) calloc(Ptmp.nnz,sizeof(REAL));

  for (i=0; i<nedge; i++) {
    oneoverlen = 1.0/(mesh->ed_len[i]);
    rowa = mesh->ed_v->IA[i];
    Ptmp.IA[i] = rowa + i*(dim-1)*2;
    cola = Ptmp.IA[i];
    j = mesh->ed_v->JA[rowa];
    k = mesh->ed_v->JA[rowa+1];
    if(j>k) {
      xL = 0.5*oneoverlen*((mesh->cv->x[j])-(mesh->cv->x[k]));
      yL = 0.5*oneoverlen*((mesh->cv->y[j])-(mesh->cv->y[k]));
      if(dim==3)
        zL = 0.5*oneoverlen*((mesh->cv->z[j])-(mesh->cv->z[k]));
    } else {
      xL = 0.5*oneoverlen*((mesh->cv->x[k])-(mesh->cv->x[j]));
      yL = 0.5*oneoverlen*((mesh->cv->y[k])-(mesh->cv->y[j]));
      if(dim==3)
        zL = 0.5*oneoverlen*((mesh->cv->z[k])-(mesh->cv->z[j]));
    }
    Ptmp.JA[cola] = j*dim;
    Ptmp.val[cola] = xL;
    Ptmp.JA[cola+dim] = k*dim;
    Ptmp.val[cola+dim] = xL;
    Ptmp.JA[cola+1] = j*dim+1;
    Ptmp.val[cola+1] = yL;
    Ptmp.JA[cola+dim+1] = k*dim+1;
    Ptmp.val[cola+dim+1] = yL;
    if(dim==3) {
      Ptmp.JA[cola+2] = j*dim+2;
      Ptmp.val[cola+2] = zL;
      Ptmp.JA[cola+dim+2] = k*dim+2;
      Ptmp.val[cola+dim+2] = zL;
    }
  }

  *Pgrad = Ptmp;

  return;
}
/***********************************************************************************************/

/***********************************************************************************************/
/*!
 * \fn void get_Pigrad_H1toRT( dCSRmat* Pdiv, dCSRmat* Pcurl, dCSRmat* Curl, mesh_struct* mesh)
 *
 * \brief Computes Raviart-Thomas interpolation operator for HX preconditioner
 *
 * \param Pgrad 	     dCSRmat matrix that takes the the \Pi for HX preconditioner.
 * \param mesh         Mesh data
 *
 * \note This only makes sense in 2D or 3D and P1 elements.
 *       Also assumes shuffled ordering.
 *
 */
/***********************************************************************************************/
void get_Pigrad_H1toRT( dCSRmat* Pdiv, dCSRmat* Pcurl, dCSRmat* Curl, mesh_struct* mesh)
{
  INT i,j,rowa,cola;
  INT v1,v2,v3;
  INT begin_row,end_row;
  INT nface = mesh->nface;
  INT nedge = mesh->nedge;
  INT dim = mesh->dim;

  REAL temp1;
  REAL temp2;
  REAL temp3;

  dCSRmat Ptmp;

  Ptmp.row = mesh->f_v->row;
  Ptmp.col = (mesh->f_v->col)*dim;
  Ptmp.nnz = (mesh->f_v->nnz)*dim;
  Ptmp.IA = (INT *) calloc(Ptmp.row+1,sizeof(INT));
  Ptmp.JA = (INT *) calloc(Ptmp.nnz,sizeof(INT));
  Ptmp.val = (REAL *) calloc(Ptmp.nnz,sizeof(REAL));

  // Need to get the f-th component of
  //    -Curl*Pcurl*(y,0,0)^T
  REAL* PcurlY = (REAL*)calloc(nedge,sizeof(REAL));
  REAL* PcurlZ = (REAL*)calloc(nedge,sizeof(REAL));
  REAL* PcurlX = (REAL*)calloc(nedge,sizeof(REAL));

  for(i=0; i<nedge; i++) {
    rowa = mesh->ed_v->IA[i];
    v1 = mesh->ed_v->JA[rowa];
    v2 = mesh->ed_v->JA[rowa+1];
    cola = Pcurl->IA[i];
    PcurlY[i] = (mesh->cv->y[v1])*(Pcurl->val[cola])   + (mesh->cv->y[v2])*(Pcurl->val[cola+dim]);
    PcurlZ[i] = (mesh->cv->z[v1])*(Pcurl->val[cola+1]) + (mesh->cv->z[v2])*(Pcurl->val[cola+dim+1]);
    PcurlX[i] = (mesh->cv->x[v1])*(Pcurl->val[cola+2]) + (mesh->cv->x[v2])*(Pcurl->val[cola+dim+2]);

//    // PiV1 * y
//    cola = Pcurl->IA[i]-1;
//    j = Pcurl->JA[cola+1]-1;
//    k = Pcurl->JA[cola+dim+1]-1;
//    PcurlY[i] = (mesh->cv->y[v1])*(Pcurl->val[j]) + (mesh->cv->y[v2])*(Pcurl->val[k]);
//    // PiV2 * z
//    cola = Pcurl->IA[i]-1;
//    j = Pcurl->JA[cola+2]-1;
//    k = Pcurl->JA[cola+dim+2]-1;
//    PcurlZ[i] = (mesh->cv->z[v1])*(Pcurl->val[j]) + (mesh->cv->z[v2])*(Pcurl->val[k]);
//    // PiV3 * x
//    cola = Pcurl->IA[i]-1;
//    j = Pcurl->JA[cola]-1;
//    k = Pcurl->JA[cola+dim]-1;
//    PcurlX[i] = (mesh->cv->x[v1])*(Pcurl->val[j]) + (mesh->cv->x[v2])*(Pcurl->val[k]);
  }


  for(i=0; i<nface; i++) {
    temp1 = 0.0;
    temp2 = 0.0;
    temp3 = 0.0;
    //TODO: check if Curl and PcurlX is indexed correctly to do this.
    begin_row = Curl->IA[i];
    end_row = Curl->IA[i+1];
    for(j=begin_row; j<end_row; j++){
      temp1 += -(Curl->val[j])*PcurlZ[Curl->JA[j]];
      temp2 += -(Curl->val[j])*PcurlX[Curl->JA[j]];
      temp3 += -(Curl->val[j])*PcurlY[Curl->JA[j]];
    }
    // divide by number of vertices per face
    temp1 = temp1 / 3;
    temp2 = temp2 / 3;
    temp3 = temp3 / 3;
//    temp1 = temp1 / mesh->f_area[i];
//    temp2 = temp2 / mesh->f_area[i];
//    temp3 = temp3 / mesh->f_area[i];

    //printf("F_norm(%f,%f,%f)\ttemp(%f,%f,%f)\n",mesh->f_norm[i*dim],mesh->f_norm[i*dim+1],mesh->f_norm[i*dim+2],temp1,temp2,temp3);

    //temp1 = mesh->f_area[i]*mesh->f_norm[i*dim+0]/3;
    //temp2 = mesh->f_area[i]*mesh->f_norm[i*dim+1]/3;
    //temp3 = mesh->f_area[i]*mesh->f_norm[i*dim+2]/3;
    temp1 = mesh->f_norm[i*dim+0]/3;
    temp2 = mesh->f_norm[i*dim+1]/3;
    temp3 = mesh->f_norm[i*dim+2]/3;

    // Build Pdiv
    rowa = mesh->f_v->IA[i];
    Ptmp.IA[i] = rowa + i*3*(dim-1);
    cola = Ptmp.IA[i];

    v1 = mesh->f_v->JA[rowa];
    v2 = mesh->f_v->JA[rowa+1];
    v3 = mesh->f_v->JA[rowa+2];

    Ptmp.JA[cola]       = v1*dim;
    Ptmp.val[cola]      = temp1;
    Ptmp.JA[cola+dim]   = v2*dim;
    Ptmp.val[cola+dim]  = temp1;
    Ptmp.JA[cola+2*dim] = v3*dim;
    Ptmp.val[cola+2*dim]= temp1;

    Ptmp.JA[cola+1]       = v1*dim+1;
    Ptmp.val[cola+1]      = temp2;
    Ptmp.JA[cola+dim+1]   = v2*dim+1;
    Ptmp.val[cola+dim+1]  = temp2;
    Ptmp.JA[cola+2*dim+1] = v3*dim+1;
    Ptmp.val[cola+2*dim+1]= temp2;

    Ptmp.JA[cola+2]       = v1*dim+2;
    Ptmp.val[cola+2]      = temp3;
    Ptmp.JA[cola+dim+2]   = v2*dim+2;
    Ptmp.val[cola+dim+2]  = temp3;
    Ptmp.JA[cola+2*dim+2] = v3*dim+2;
    Ptmp.val[cola+2*dim+2]= temp3;
  }
  Ptmp.IA[nface] = Ptmp.nnz;

  *Pdiv = Ptmp;

  free(PcurlY);
  free(PcurlZ);
  free(PcurlX);
}
/***********************************************************************************************/

/*!
 * \fn void ProjectOut_Grad(dvector* u,fespace* FE_H1,fespace* FE_Ned,mesh_struct* mesh,qcoordinates* cq,dCSRmat* G)
 *
 * \brief Takes an approximation in H(curl) using lowest-order Nedelec FE space
 *        and projects out the gradient from it's decomposition:
 *
 *        u = grad p + curl q
 *
 *        Thus, we remove the grad p, so that div u = 0 (at least weakly).
 *        Here p is in H0^1 and we assume P1 elements.
 *          u <-- u - grad p
 *        p is found by solving discrete Laplacian: <grad p, grad v> = <E, grad v>
 *
 * \param u       Nedelec Approximation with gradient kernel removed
 * \param FE_H1   P1 FE space
 * \param FE_Ned  Nedelec FE space
 * \param mesh    Mesh Data
 * \param cq      Quadrature points
 * \param G       Gradient matrix in dCSRmat format
 *
 * \note Assumes 2D or 3D.
 *
 */
void ProjectOut_Grad(dvector* u,fespace* FE_H1,fespace* FE_Ned,mesh_struct* mesh,qcoordinates* cq,dCSRmat* G)
{
  // Construct the Laplacian using P1 elements: <grad p, grad v>
  dCSRmat Alap;
  assemble_global(&Alap,NULL,assemble_DuDv_local,FE_H1,mesh,cq,NULL,NULL,0.0);

  // Construct RHS vector: <E,grad v>
  dvector b;
  assemble_global_Ned_GradH1_RHS(&b,FE_H1,FE_Ned,mesh,cq,u);

  // Eliminate Boundary Conditions (p=0 on boundary)
  eliminate_DirichletBC(zero_coeff_scal,FE_H1,mesh,&b,&Alap,0.0); //

  // Solve Laplacian System to get p: <grad p, grad v> = <E, grad v> -> use CG
  // Allocate the solution and set initial guess to be all zero (except at boundaries)
  dvector p = dvec_create(mesh->nv);
  dvec_set(p.row, &p, 0.0);

  // Solve for p
  dcsr_pcg(&Alap,&b,&p,NULL,1e-15,50000,1,0);

  // Multiply by Gradient
  dvector gradp = dvec_create(mesh->nedge);
  dcsr_mxv(G,p.val,gradp.val);

  // Update E
  dvec_axpy(-1.0,&gradp,u);

  // Free Vectors
  dcsr_free(&Alap);
  dvec_free(&b);
  dvec_free(&p);
  dvec_free(&gradp);

  return;
}
/******************************************************************************/


/******************************************************************************/
/**********************Deprecated versions*************************************/
/******************************************************************************/

/*!
 * \fn void FE_Interpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,fespace *FE,mesh_struct *mesh)
 *
 * \brief Interpolate a finite-element approximation to any other point in the given element using the given type of elements.
 *
 * \param u 	      Approximation to interpolate
 * \param x           Coordinates where to compute value
 * \param dof_on_elm  DOF belonging to particular element
 * \param v_on_elm    Vertices belonging to particular element
 * \param FE          FE Space
 * \param mesh        Mesh Data
 * \param val         Pointer to value of approximation at given values
 *
 * \note This will be deprecated soon.
 *
 */
void FE_Interpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,fespace *FE,mesh_struct *mesh)
{
  INT i,j,dof;

  // Get FE and Mesh data
  INT dof_per_elm = FE->dof_per_elm;
  //INT FEtype = FE->FEtype;
  INT scal_or_vec = FE->scal_or_vec;
  INT dim = mesh->dim;

  REAL coef[dim];

  get_FEM_basis(FE->phi,FE->dphi,x,v_on_elm,dof_on_elm,mesh,FE);

  if(scal_or_vec==0) { // Scalar Element
    coef[0] = 0.0;
    for(j=0; j<dof_per_elm; j++) {
      dof = dof_on_elm[j];
      coef[0] += u[dof]*FE->phi[j];
    }
    val[0] = coef[0];
  } else { // Vector Element
    for(i=0;i<dim;i++) {
      coef[i] = 0.0;
      for(j=0; j<dof_per_elm; j++) {
        dof = dof_on_elm[j];
        coef[i] += u[dof]*FE->phi[j*dim+i];
      }
      val[i] = coef[i];
    }
  }

  return;
}
/******************************************************************************/

/*!
 * \fn void FE_DerivativeInterpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,fespace *FE,mesh_struct *mesh)
 *
 * \brief Interpolate the "derivative" of a finite-element approximation to any other point in the given element using the given type of elements.
 *        Note that for Lagrange Elements this means the Gradient, grad u, for Nedelec it means the Curl, curl u, and for RT it is the Divergence, div u.
 *
 * \param u 	      Approximation to interpolate
 * \param x           Coordinates where to compute value
 * \param dof_on_elm  DOF belonging to particular element
 * \param v_on_elm    Vertices belonging to particular element
 * \param FE          FE Space
 * \param mesh        Mesh Data
 * \param val         Pointer to value of approximation at given values
 *
 * \note This will be deprecated soon.
 *
 */
void FE_DerivativeInterpolation(REAL* val,REAL *u,REAL *x,INT *dof_on_elm,INT *v_on_elm,fespace *FE,mesh_struct *mesh)
{
  INT dof,j,k,i;

  // Get FE and Mesh data
  INT dof_per_elm = FE->dof_per_elm;
  INT FEtype = FE->FEtype;
  INT dim = mesh->dim;

  // Basis Functions and its derivatives if necessary
  //REAL coef[dim];
  REAL coef[dim*dim];

  get_FEM_basis(FE->phi,FE->dphi,x,v_on_elm,dof_on_elm,mesh,FE);

  if(FEtype==0 || FEtype==99) { // Don't compute derivatives of P0 elements (set to 0)
    for(j=0;j<dim;j++) {
      val[j] = 0.0;
    }
  } else if(FEtype<20 && FEtype!=0) { // Scalar Element
    for(j=0;j<dim;j++) {
      coef[j] = 0.0;
      for(k=0; k<dof_per_elm; k++) {
        dof = dof_on_elm[k];
        coef[j] += u[dof]*FE->dphi[k*dim+j];
      }
      val[j] = coef[j];
    }
  } else if (FEtype==20) { // Nedelec
    for(j=0;j<dim;j++)
      coef[j] = 0.0;
    if (dim==2) { // Curl is scalar
      for (j=0; j<dof_per_elm; j++) {
        dof = dof_on_elm[j];
        coef[0] += u[dof]*FE->dphi[j];
      }
      val[0] = coef[0];
    } else if (dim==3) { // Curl is vector
      for (j=0; j<dof_per_elm; j++) {
        dof = dof_on_elm[j];
        coef[0] += u[dof]*FE->dphi[j*dim+0];
        coef[1] += u[dof]*FE->dphi[j*dim+1];
        coef[2] += u[dof]*FE->dphi[j*dim+2];
      }
      val[0] = coef[0];
      val[1] = coef[1];
      val[2] = coef[2];
    }
  } else if (FEtype==30) { // Raviart-Thomas (div is scalar)
    coef[0] = 0.0;

    for (j=0; j<dof_per_elm; j++) {
      dof = dof_on_elm[j];
      coef[0] += u[dof]*FE->dphi[j];
    }
    val[0] = coef[0];
  } else if (FEtype>=60) { // Vector Lagrange Functions (i.e. bubbles) -> gradients of vectors -> tensor
    for(j=0;j<dim;j++) {
      for(i=0;i<dim;i++) {
        coef[j*dim + i] = 0.0;
        for(k=0; k<dof_per_elm; k++) {
          dof = dof_on_elm[k];
          coef[j*dim + i] += u[dof]*FE->dphi[k*dim*dim + j*dim + i];
        }
        val[j*dim + i] = coef[j*dim + i];
      }
    }

  } else {
    check_error(ERROR_FE_TYPE,__FUNCTION__);
  }

  return;
}
/******************************************************************************/

/*!
 * \fn void blockFE_Interpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,block_fespace *FE,mesh_struct *mesh)
 *
 * \brief Interpolate a block finite-element approximation to any other point in the given element using the given type of elements.
 *
 * \param u 	        Approximation to interpolate in block form
 * \param x           Coordinates where to compute value
 * \param dof_on_elm  DOF belonging to particular element
 * \param v_on_elm    Vertices belonging to particular element
 * \param FE          Block FE Space
 * \param mesh        Mesh Data
 * \param nun         Number of unknowns in u (1 is a scalar)
 * \param val         Pointer to value of approximation at given values
 *
 */
void blockFE_Interpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,block_fespace *FE,mesh_struct *mesh)
{
  INT k;
  INT dim = mesh->dim;

  INT* local_dof_on_elm = dof_on_elm;
  REAL* val_sol = val;
  REAL* u_comp = u;

  for(k=0;k<FE->nspaces;k++) {
    FE_Interpolation(val_sol,u_comp,x,local_dof_on_elm,v_on_elm,FE->var_spaces[k],mesh);
    if(FE->var_spaces[k]->scal_or_vec==0) { // Scalar
      val_sol++;
    } else { // Vector
      val_sol += dim;
    }
    u_comp += FE->var_spaces[k]->ndof;
    local_dof_on_elm += FE->var_spaces[k]->dof_per_elm;
  }

  return;
}
/******************************************************************************/

/*!
 * \fn void blockFE_DerivativeInterpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,block_fespace *FE,mesh_struct *mesh)
 *
 * \brief Interpolate the "derivative" of a block finite-element approximation to any other point in the given
 *        element using the given type of elements.  Note that for Lagrange Elements this means the Gradient, grad u,
 *        for Nedelec it means the Curl, curl u, and for RT it is the Divergence, div u.
 *
 *
 * \param u 	        Approximation to interpolate in block form
 * \param x           Coordinates where to compute value
 * \param dof_on_elm  DOF belonging to particular element
 * \param v_on_elm    Vertices belonging to particular element
 * \param FE          Block FE Space
 * \param mesh        Mesh Data
 * \param nun         Number of unknowns in u (1 is a scalar)
 * \param val         Pointer to value of approximation at given values
 *
 */
void blockFE_DerivativeInterpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,block_fespace *FE,mesh_struct *mesh)
{
  INT k;
  INT dim = mesh->dim;

  INT* local_dof_on_elm = dof_on_elm;
  REAL* val_sol = val;
  REAL* u_comp = u;

  for(k=0;k<FE->nspaces;k++) {
    FE_DerivativeInterpolation(val_sol,u_comp,x,local_dof_on_elm,v_on_elm,FE->var_spaces[k],mesh);
    if(FE->var_spaces[k]->FEtype<20 || FE->var_spaces[k]->FEtype==99) { // Scalar
      val_sol += dim;
    } else if(FE->var_spaces[k]->FEtype==20 && dim==2) { // Curl in 2D is Scalar
      val_sol++;
    } else if(FE->var_spaces[k]->FEtype==20 && dim==3) { // Curl in 3D is Vector
      val_sol+=dim;
    } else if(FE->var_spaces[k]->FEtype==30) { // Div is Scalar
      val_sol++;
    } else if(FE->var_spaces[k]->FEtype>=60) { // VectorLagrange (i.e., bubbles) Gradient is matrix.
      val_sol+=dim*dim;
    } else {
      check_error(ERROR_FE_TYPE,__FUNCTION__);
    }
    u_comp += FE->var_spaces[k]->ndof;
    local_dof_on_elm += FE->var_spaces[k]->dof_per_elm;
  }

  return;
}
/******************************************************************************/
