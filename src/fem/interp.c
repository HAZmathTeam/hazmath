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
 *
 */

#include "hazmath.h"

/**************** Interpolation Routines ****************************/
/********************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void FE_Interpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,fespace *FE,trimesh *mesh)
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
 */
void FE_Interpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,fespace *FE,trimesh *mesh)
{
  INT i,j,dof;

  // Get FE and Mesh data
  INT dof_per_elm = FE->dof_per_elm;
  INT FEtype = FE->FEtype;
  INT dim = mesh->dim;

  REAL coef[dim];

  get_FEM_basis(FE->phi,FE->dphi,x,v_on_elm,dof_on_elm,mesh,FE);

  if(FEtype<20) { // Scalar Element
    coef[0] = 0.0;
    for(j=0; j<dof_per_elm; j++) {
      dof = dof_on_elm[j] - 1;
      coef[0] += u[dof]*FE->phi[j];
    }
    val[0] = coef[0];
  } else { // Vector Element
    for(i=0;i<dim;i++) {
      coef[i] = 0.0;
      for(j=0; j<dof_per_elm; j++) {
        dof = dof_on_elm[j] - 1;
        coef[i] += u[dof]*FE->phi[j*dim+i];
      }
      val[i] = coef[i];
    }
  }

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void FE_DerivativeInterpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,fespace *FE,trimesh *mesh)
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
 */
void FE_DerivativeInterpolation(REAL* val,REAL *u,REAL *x,INT *dof_on_elm,INT *v_on_elm,fespace *FE,trimesh *mesh)
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

  if(FEtype==0) { // Don't compute derivatives of P0 elements (set to 0)
    for(j=0;j<dim;j++) {
      val[j] = 0.0;
    }
  } else if(FEtype<20 && FEtype!=0) { // Scalar Element
    for(j=0;j<dim;j++) {
      coef[j] = 0.0;
      for(k=0; k<dof_per_elm; k++) {
        dof = dof_on_elm[k] - 1;
        coef[j] += u[dof]*FE->dphi[k*dim+j];
      }
      val[j] = coef[j];
    }
  } else if (FEtype==20) { // Nedelec
    for(j=0;j<dim;j++)
      coef[j] = 0.0;
    if (dim==2) { // Curl is scalar
      for (j=0; j<dof_per_elm; j++) {
        dof = dof_on_elm[j]-1;
        coef[0] += u[dof]*FE->dphi[j];
      }
      val[0] = coef[0];
    } else if (dim==3) { // Curl is vector
      for (j=0; j<dof_per_elm; j++) {
        dof = dof_on_elm[j]-1;
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
      dof = dof_on_elm[j]-1;
      coef[0] += u[dof]*FE->dphi[j];
    }
    val[0] = coef[0];
  } else if (FEtype==61) { // bubble
    for(j=0;j<dim;j++) {
      for(i=0;i<dim;i++) {
        coef[j*dim + i] = 0.0;
        for(k=0; k<dof_per_elm; k++) {
          dof = dof_on_elm[k] - 1;
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
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void mult_FE_Interpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,fespace *FE,trimesh *mesh,INT nun)
 *
 * \brief Interpolate a finite-element approximation to any other point in the given element using the given type of elements.
 *
 * \note This allows you to interpolate multiple approximations at once.  OUTDATED
 *
 * \param u 	      Approximation to interpolate
 * \param x           Coordinates where to compute value
 * \param dof_on_elm  DOF belonging to particular element
 * \param v_on_elm    Vertices belonging to particular element
 * \param FE          FE Space
 * \param mesh        Mesh Data
 * \param nun         Number of unknowns in u (1 is a scalar)
 * \param val         Pointer to value of approximation at given values
 *
 */
void mult_FE_Interpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,fespace *FE,trimesh *mesh,INT nun)
{
  INT i,j,k,dof;

  // Get FE and Mesh data
  INT dof_per_elm = FE->dof_per_elm;
  INT FEtype = FE->FEtype;
  INT dim = mesh->dim;
  INT ndof = FE->ndof;

  REAL coef[dim];

  get_FEM_basis(FE->phi,FE->dphi,x,v_on_elm,dof_on_elm,mesh,FE);

  if(FEtype<20) { // Scalar Element
    for(i=0; i<nun; i++) {
      coef[0] = 0.0;
      for(j=0; j<dof_per_elm; j++) {
        dof = i*ndof + dof_on_elm[j] - 1;
        coef[0] += u[dof]*FE->phi[j];
      }
      val[i] = coef[0];
    }
  } else { // Vector Element
    for(k=0; k<nun; k++) {
      for(i=0;i<dim;i++) {
        coef[i] = 0.0;
        for(j=0; j<dof_per_elm; j++) {
          dof = dof_on_elm[j] - 1;
          coef[i] += u[dof]*FE->phi[j*dim+i];
        }
        val[k*dim + i] = coef[i];
      }
    }
  }

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void mult_FE_DerivativeInterpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,fespace *FE,trimesh *mesh,INT nun)
 *
 * \brief Interpolate the "derivative" of a finite-element approximation to any other point in the given element using the given type of elements.
 *        Note that for Lagrange Elements this means the Gradient, grad u, for Nedelec it means the Curl, curl u, and for RT it is the Divergence, div u.
 *
 * \note This routine allows you to interpolate multiple approximations at once. OUTDATED
 *
 * \param u 	        Approximation to interpolate
 * \param x           Coordinates where to compute value
 * \param dof_on_elm  DOF belonging to particular element
 * \param v_on_elm    Vertices belonging to particular element
 * \param FE          FE Space
 * \param mesh        Mesh Data
 * \param nun         Number of unknowns in u (1 is a scalar)
 * \param val         Pointer to value of approximation at given values
 *
 */
void mult_FE_DerivativeInterpolation(REAL* val,REAL *u,REAL *x,INT *dof_on_elm,INT *v_on_elm,fespace *FE,trimesh *mesh,INT nun)
{
  INT i,dof,j,k;

  // Get FE and Mesh data
  INT dof_per_elm = FE->dof_per_elm;
  INT FEtype = FE->FEtype;
  INT dim = mesh->dim;
  INT ndof = FE->ndof;

  // Basis Functions and its derivatives if necessary
  REAL coef[dim];

  get_FEM_basis(FE->phi,FE->dphi,x,v_on_elm,dof_on_elm,mesh,FE);

  if(FEtype==0) { // Don't compute derivatives of P0 elements (set to 0)
    for(i=0; i<nun; i++) {
      for(j=0;j<dim;j++) {
        val[i*dim+j] = 0.0;
      }
    }
  } else if(FEtype<20 && FEtype!=0) { // Scalar Element
    for(i=0; i<nun; i++) {
      for(j=0;j<dim;j++) {
        coef[j] = 0.0;
        for(k=0; k<dof_per_elm; k++) {
          dof = i*ndof + dof_on_elm[k] - 1;
          coef[j] += u[dof]*FE->dphi[k*dim+j];
        }
        val[i*dim+j] = coef[j];
      }
    }
  } else if (FEtype==20) { // Nedelec
    for(j=0;j<dim;j++)
      coef[j] = 0.0;
    if (dim==2) { // Curl is scalar
      for (j=0; j<dof_per_elm; j++) {
        dof = dof_on_elm[j]-1;
        coef[0] += u[dof]*FE->dphi[j];
      }
      val[0] = coef[0];
    } else if (dim==3) { // Curl is vector
      for (j=0; j<dof_per_elm; j++) {
        dof = dof_on_elm[j]-1;
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
      dof = dof_on_elm[j]-1;
      coef[0] += u[dof]*FE->dphi[j];
    }
    val[0] = coef[0];
  } else {
    check_error(ERROR_FE_TYPE,__FUNCTION__);
  }

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void blockFE_Interpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,block_fespace *FE,trimesh *mesh)
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
void blockFE_Interpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,block_fespace *FE,trimesh *mesh)
{
  INT k;
  INT dim = mesh->dim;

  INT* local_dof_on_elm = dof_on_elm;
  REAL* val_sol = val;
  REAL* u_comp = u;

  for(k=0;k<FE->nspaces;k++) {
    FE_Interpolation(val_sol,u_comp,x,local_dof_on_elm,v_on_elm,FE->var_spaces[k],mesh);
    if(FE->var_spaces[k]->FEtype<20) { // Scalar
      val_sol++;
    } else { // Vector
      val_sol += dim;
    }
    u_comp += FE->var_spaces[k]->ndof;
    local_dof_on_elm += FE->var_spaces[k]->dof_per_elm;
  }

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void blockFE_DerivativeInterpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,block_fespace *FE,trimesh *mesh)
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
void blockFE_DerivativeInterpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT *v_on_elm,block_fespace *FE,trimesh *mesh)
{
  INT k;
  INT dim = mesh->dim;

  INT* local_dof_on_elm = dof_on_elm;
  REAL* val_sol = val;
  REAL* u_comp = u;

  for(k=0;k<FE->nspaces;k++) {
    FE_DerivativeInterpolation(val_sol,u_comp,x,local_dof_on_elm,v_on_elm,FE->var_spaces[k],mesh);
    if(FE->var_spaces[k]->FEtype<20) { // Scalar
      val_sol += dim;
    } else if(FE->var_spaces[k]->FEtype==20 && dim==2) { // Curl in 2D is Scalar
      val_sol++;
    } else if(FE->var_spaces[k]->FEtype==20 && dim==3) { // Curl in 3D is Vector
      val_sol+=dim;
    } else if(FE->var_spaces[k]->FEtype==30) { // Div is Scalar
      val_sol++;
    } else if(FE->var_spaces[k]->FEtype==61) { // bubble
      val_sol+=dim*dim;
    } else {
      check_error(ERROR_FE_TYPE,__FUNCTION__);
    }
    u_comp += FE->var_spaces[k]->ndof;
    local_dof_on_elm += FE->var_spaces[k]->dof_per_elm;
  }

  return;
}
/****************************************************************************************************************************/

/**************** Projection Routines *******************************/
/********************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn REAL FE_Evaluate_DOF(void (*expr)(REAL *,REAL *,REAL,void *),fespace *FE,trimesh *mesh,REAL time,INT DOF)
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
REAL FE_Evaluate_DOF(void (*expr)(REAL *,REAL *,REAL,void *),fespace *FE,trimesh *mesh,REAL time,INT DOF)
{
  REAL* x = (REAL *) calloc(mesh->dim,sizeof(REAL));
  INT dim = mesh->dim;
  INT FEtype = FE->FEtype;
  INT nq1d = 3; // If quadrature not given, fix order.
  REAL val=-666e+00;

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
    val = (1.0/mesh->ed_len[DOF])*integrate_edge_vector_tangent(expr,1,0,nq1d,NULL,mesh,time,DOF);

  // Raviart-Thomas u[dof] = 1/farea \int_face u*n_face
  } else if (FEtype==30) {
    val = (1.0/mesh->f_area[DOF])*integrate_face_vector_normal(expr,1,0,nq1d,NULL,mesh,time,DOF);

  // Bubbles
  } else if (FEtype>=60 && FEtype<70) {
    val = (1.0/mesh->f_area[DOF])*integrate_face_vector_normal(expr,1,0,nq1d,NULL,mesh,time,DOF);

  // No other FEM types implemented
  } else {
    check_error(ERROR_FE_TYPE,__FUNCTION__);
  }

  if (x) free(x);

  return val;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void FE_Evaluate(REAL* val,void (*expr)(REAL *,REAL *,REAL,void *),fespace *FE,trimesh *mesh,REAL time)
 *
 * \brief Evaluate a given analytical function on the finite-element space given.
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
void FE_Evaluate(REAL* val,void (*expr)(REAL *,REAL *,REAL,void *),fespace *FE,trimesh *mesh,REAL time)
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
 * \fn REAL blockFE_Evaluate_DOF(void (*expr)(REAL *,REAL *,REAL,void *),block_fespace *FE,trimesh *mesh,REAL time,INT comp,INT DOF)
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
REAL blockFE_Evaluate_DOF(void (*expr)(REAL *,REAL *,REAL,void *),block_fespace *FE,trimesh *mesh,REAL time,INT comp,INT DOF)
{
  int i,j,m;
  INT nq1d = 3; // Number of quadrature points in 1D if not given.
  REAL* x = (REAL *) calloc(mesh->dim,sizeof(REAL));
  REAL* valx = (REAL *) calloc(FE->nun,sizeof(REAL));
  INT dim = mesh->dim;
  INT local_dim = 0;
  REAL val=-666e+00;

  INT* face_vertex = (INT *) calloc(dim,sizeof(INT));

  for(i=0;i<comp;i++) {
    if(FE->var_spaces[i]->FEtype>=0 && FE->var_spaces[i]->FEtype<10) { // Scalar Element
      local_dim += 1;
    } else if(FE->var_spaces[i]->FEtype == 61) { // Bubble Element
      local_dim += 0;
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

  // Bubbles (PETER HELP??)
  } else if (FE->var_spaces[comp]->FEtype>=60 && FE->var_spaces[comp]->FEtype<70) {
    x[0] = mesh->f_mid[DOF*dim];
    x[1] = mesh->f_mid[DOF*dim+1];
    if(dim==3) x[2] = mesh->f_mid[DOF*dim+2];
    (*expr)(valx,x,time,&(FE->var_spaces[comp]->dof_flag[DOF]));
    val = 0.0;
    for(j=0;j<dim;j++) val+=mesh->f_norm[DOF*dim+j]*valx[local_dim + j];
    get_incidence_row(DOF,mesh->f_v,face_vertex);
    for (m=0;m<dim;m++) {
      x[0] = mesh->cv->x[face_vertex[m]-1];
      x[1] = mesh->cv->y[face_vertex[m]-1];
      if(dim==3) x[2] = mesh->cv->z[face_vertex[m]-1];
      // The following only works for 2D
      (*expr)(valx,x,time,&(FE->var_spaces[0]->dof_flag[DOF]));
      //for(j=0;j<dim;j++) val+= -(1.0/dim)*mesh->f_area[DOF]*mesh->f_norm[DOF*dim+j]*valx[local_dim + j];
      for(j=0;j<dim;j++) val+= -(1.0/dim)*mesh->f_norm[DOF*dim+j]*valx[local_dim + j];
    }

  // Not a FEM implemented
  } else {
    check_error(ERROR_FE_TYPE,__FUNCTION__);
  }

  if (x) free(x);
  if(valx) free(valx);
  if (face_vertex) free(face_vertex);

  return val;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void blockFE_Evaluate(REAL* val,void (*expr)(REAL *,REAL *,REAL,void *),block_fespace *FE,trimesh *mesh,REAL time)
 *
 * \brief Evaluate a given analytical function on the block finite-element space given.
 *
 * \param expr 	    Function call to analytical expression
 * \param FE        block FE Space for multiple variables
 * \param mesh      Mesh Data
 * \param time      Physical time to evaluate function at
 * \param val       FE approximation of function on fespace
 *
 */
void blockFE_Evaluate(REAL* val,void (*expr)(REAL *,REAL *,REAL,void *),block_fespace *FE,trimesh *mesh,REAL time)
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

/***********************************************************************************************/
/*!
 * \fn void Project_to_Vertices(REAL* u_on_V,REAL *u,fespace *FE,trimesh *mesh)
 *
 * \brief Interpolate a finite-element approximation to the vertices of a mesh
 *
 * \param u_on_V  Solution on vertices of mesh
 * \param u       Approximation to Interpolate
 * \param FE      FE space
 * \param mesh    Mesh Data
 *
 */
void Project_to_Vertices(REAL* u_on_V,REAL *u,fespace *FE,trimesh *mesh)
{
  INT i,k,j;
  REAL* val = NULL;
  INT dim = mesh->dim;
  REAL* x = (REAL *) calloc(dim,sizeof(REAL));

  // Get FE and Mesh data
  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT FEtype = FE->FEtype;
  INT nelm = mesh->nelm;
  INT nv = mesh->nv;

  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  if(FEtype>=0 && FEtype<20) { // Scalar Element
    val = (REAL *) calloc(1,sizeof(REAL));
  } else { // Vector Element
    val = (REAL *) calloc(dim,sizeof(REAL));
  }

  // Loop over Elements
  for(i=0;i<nelm;i++) {

    // Find DOF for given Element
    get_incidence_row(i,FE->el_dof,dof_on_elm);

    // Find vertices for given Element
    get_incidence_row(i,mesh->el_v,v_on_elm);

    // Interpolate FE approximation to vertices
    for(j=0;j<v_per_elm;j++) {
      get_coords(x,v_on_elm[j]-1,mesh->cv,dim);
      FE_Interpolation(val,u,x,dof_on_elm,v_on_elm,FE,mesh);

      if(FEtype>=0 && FEtype<20) { // Scalar Element
        u_on_V[v_on_elm[j]-1] = val[0];
      } else { // Vector Element
        for(k=0;k<dim;k++) {
          u_on_V[k*nv + v_on_elm[j]-1] = val[k];
        }
      }
    }
  }

  if (x) free(x);
  if(val) free(val);
  if(v_on_elm) free(v_on_elm);
  if(dof_on_elm) free(dof_on_elm);
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
  int i;
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
  int i;
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
 * \fn void get_grad_H1toNed(dCSRmat* Grad,trimesh* mesh)
 *
 * \brief Computes Gradient operator.
 *        Applying the resulting matrix computes the gradient of an H1 approximation.
 *        Takes Nodal (H1) DOF vector to Edge (Nedelec) DOF vector
 *        Ordering determined by edge to node map: bigger node is +1 smaller is -1
 *
 * \param Grad	     dCSRmat matrix that takes the gradient of an H1 approximation
 * \param mesh         Mesh data
 *
 * \note This only makes sense for P1 elements and for dimensions 2 or 3.
 *
 */
void get_grad_H1toNed(dCSRmat* Grad,trimesh* mesh)
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
    rowa = mesh->ed_v->IA[i]-1;
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
 * \fn void get_curl_NedtoRT(dCSRmat* Curl,trimesh* mesh)
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
void get_curl_NedtoRT(dCSRmat* Curl,trimesh* mesh)
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
    rowa = mesh->f_v->IA[i]-1;
    rowb = mesh->f_v->IA[i+1]-1;
    jcntr=0;
    for(j=rowa;j<rowb;j++) {
      inf[jcntr] = mesh->f_v->JA[j];
      jcntr++;
    }
    // Get edges of face
    rowa = mesh->f_ed->IA[i]-1;
    rowb = mesh->f_ed->IA[i+1]-1;
    for(j=rowa;j<rowb;j++) {
      k = mesh->f_ed->JA[j];
      // Get nodes of edge
      col_a = mesh->ed_v->IA[k-1]-1;
      nd1 = mesh->ed_v->JA[col_a];
      nd2 = mesh->ed_v->JA[col_a+1];
      // Determine what other node on face is
      for(s=0;s<ndpf;s++) {
        if(inf[s]!=nd1 && inf[s]!=nd2) {
          nd3 = inf[s];
        }
      }
      vec1[0] = mesh->cv->x[nd1-1]-mesh->cv->x[nd3-1];
      vec2[0] = mesh->cv->x[nd2-1]-mesh->cv->x[nd3-1];
      vec1[1] = mesh->cv->y[nd1-1]-mesh->cv->y[nd3-1];
      vec2[1] = mesh->cv->y[nd2-1]-mesh->cv->y[nd3-1];
      vec1[2] = mesh->cv->z[nd1-1]-mesh->cv->z[nd3-1];
      vec2[2] = mesh->cv->z[nd2-1]-mesh->cv->z[nd3-1];
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
    rowa = mesh->f_ed->IA[i]-1;
    rowb = mesh->f_ed->IA[i+1]-1;
    for(j=rowa;j<rowb;j++) {
      k = mesh->f_ed->JA[j]-1;
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
 * \fn void get_div_RTtoL2(dCSRmat* Div,trimesh* mesh)
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
void get_div_RTtoL2(dCSRmat* Div,trimesh* mesh)
{
  INT i,j,rowa,rowb,rowc,rowd,face,elm1,elm2,elm_big,n_felm;
  INT nelm = mesh->nelm;
  REAL oneovervol=0.0;
  REAL farea=0.0;
  dCSRmat Dtmp;

  // We will need the face to element map as well
  // Each face will have two elements (or one if on boundary).
  iCSRmat f_el;
  icsr_trans_1(mesh->el_f,&f_el);

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
    rowa = mesh->el_f->IA[i]-1;
    rowb = mesh->el_f->IA[i+1]-1;
    for(j=rowa;j<rowb;j++) {
      face = mesh->el_f->JA[j]-1;
      farea = mesh->f_area[face];
      // Get elements of face
      rowc = f_el.IA[face]-1;
      rowd = f_el.IA[face+1]-1;
      n_felm = rowd-rowc;
      if(n_felm==1)
        Dtmp.val[j] = farea*oneovervol;
      else if(n_felm==2) {
        elm1 = f_el.JA[rowc]-1;
        elm2 = f_el.JA[rowc+1]-1;
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
 * \fn void get_Pigrad_H1toNed(dCSRmat* Pgrad,trimesh* mesh)
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
void get_Pigrad_H1toNed(dCSRmat* Pgrad,trimesh* mesh)
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
    rowa = mesh->ed_v->IA[i]-1;
    Ptmp.IA[i] = rowa+1 + i*(dim-1)*2;
    cola = Ptmp.IA[i]-1;
    j = mesh->ed_v->JA[rowa];
    k = mesh->ed_v->JA[rowa+1];
    if(j>k) {
      xL = 0.5*oneoverlen*((mesh->cv->x[j-1])-(mesh->cv->x[k-1]));
      yL = 0.5*oneoverlen*((mesh->cv->y[j-1])-(mesh->cv->y[k-1]));
      if(dim==3)
        zL = 0.5*oneoverlen*((mesh->cv->z[j-1])-(mesh->cv->z[k-1]));
    } else {
      xL = 0.5*oneoverlen*((mesh->cv->x[k-1])-(mesh->cv->x[j-1]));
      yL = 0.5*oneoverlen*((mesh->cv->y[k-1])-(mesh->cv->y[j-1]));
      if(dim==3)
        zL = 0.5*oneoverlen*((mesh->cv->z[k-1])-(mesh->cv->z[j-1]));
    }
    Ptmp.JA[cola] = (j-1)*dim+1;
    Ptmp.val[cola] = xL;
    Ptmp.JA[cola+dim] = (k-1)*dim+1;
    Ptmp.val[cola+dim] = xL;
    Ptmp.JA[cola+1] = (j-1)*dim+2;
    Ptmp.val[cola+1] = yL;
    Ptmp.JA[cola+dim+1] = (k-1)*dim+2;
    Ptmp.val[cola+dim+1] = yL;
    if(dim==3) {
      Ptmp.JA[cola+2] = (j-1)*dim+3;
      Ptmp.val[cola+2] = zL;
      Ptmp.JA[cola+dim+2] = (k-1)*dim+3;
      Ptmp.val[cola+dim+2] = zL;
    }
  }

  *Pgrad = Ptmp;

  return;
}
/***********************************************************************************************/

/***********************************************************************************************/
/*!
 * \fn void get_Pigrad_H1toRT( dCSRmat* Pdiv, dCSRmat* Pcurl, dCSRmat* Curl, trimesh* mesh)
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
void get_Pigrad_H1toRT( dCSRmat* Pdiv, dCSRmat* Pcurl, dCSRmat* Curl, trimesh* mesh)
{
  INT i,j,rowa,cola;
  INT v1,v2,v3;
  INT begin_row,end_row;
  INT nface = mesh->nface;
  INT nedge = mesh->nedge;
  //INT nvert = mesh->nv;
  INT dim = mesh->dim;

  REAL temp1;
  REAL temp2;
  REAL temp3;

  dCSRmat Ptmp;

  Ptmp.row = mesh->f_v->row;
  Ptmp.col = (mesh->f_v->col)*dim;
  Ptmp.nnz = (mesh->f_v->nnz)*dim;// ((dim+1)*nface)*dim
  Ptmp.IA = (INT *) calloc(Ptmp.row+1,sizeof(INT));
  Ptmp.JA = (INT *) calloc(Ptmp.nnz,sizeof(INT));
  Ptmp.val = (REAL *) calloc(Ptmp.nnz,sizeof(REAL));

  // Need to get the f-th component of
  //    -Curl*Pcurl*(y,0,0)^T
  REAL* PcurlY = (REAL*)calloc(nedge,sizeof(REAL));
  REAL* PcurlZ = (REAL*)calloc(nedge,sizeof(REAL));
  REAL* PcurlX = (REAL*)calloc(nedge,sizeof(REAL));

  for(i=0; i<nedge; i++) {
    rowa = mesh->ed_v->IA[i]-1;
    v1 = mesh->ed_v->JA[rowa]-1;
    v2 = mesh->ed_v->JA[rowa+1]-1;
    cola = Pcurl->IA[i]-1;
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
    begin_row = Curl->IA[i]-1;
    end_row = Curl->IA[i+1]-1;
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
    rowa = mesh->f_v->IA[i]-1;
    Ptmp.IA[i] = rowa+1 + i*3*(dim-1);
    cola = Ptmp.IA[i]-1;

    v1 = mesh->f_v->JA[rowa];
    v2 = mesh->f_v->JA[rowa+1];
    v3 = mesh->f_v->JA[rowa+2];

    Ptmp.JA[cola]       = (v1-1)*dim+1;
    Ptmp.val[cola]      = temp1;
    Ptmp.JA[cola+dim]   = (v2-1)*dim+1;
    Ptmp.val[cola+dim]  = temp1;
    Ptmp.JA[cola+2*dim] = (v3-1)*dim+1;
    Ptmp.val[cola+2*dim]= temp1;

    Ptmp.JA[cola+1]       = (v1-1)*dim+2;
    Ptmp.val[cola+1]      = temp2;
    Ptmp.JA[cola+dim+1]   = (v2-1)*dim+2;
    Ptmp.val[cola+dim+1]  = temp2;
    Ptmp.JA[cola+2*dim+1] = (v3-1)*dim+2;
    Ptmp.val[cola+2*dim+1]= temp2;

    Ptmp.JA[cola+2]       = (v1-1)*dim+3;
    Ptmp.val[cola+2]      = temp3;
    Ptmp.JA[cola+dim+2]   = (v2-1)*dim+3;
    Ptmp.val[cola+dim+2]  = temp3;
    Ptmp.JA[cola+2*dim+2] = (v3-1)*dim+3;
    Ptmp.val[cola+2*dim+2]= temp3;
  }
  Ptmp.IA[nface] = Ptmp.nnz+1;

  *Pdiv = Ptmp;

  free(PcurlY);
  free(PcurlZ);
  free(PcurlX);
}
/***********************************************************************************************/

/*!
 * \fn void ProjectOut_Grad(dvector* u,fespace* FE_H1,fespace* FE_Ned,trimesh* mesh,qcoordinates* cq,dCSRmat* G)
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
void ProjectOut_Grad(dvector* u,fespace* FE_H1,fespace* FE_Ned,trimesh* mesh,qcoordinates* cq,dCSRmat* G)
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
  dcsr_mxv_1(G,p.val,gradp.val);

  // Update E
  dvec_axpy(-1.0,&gradp,u);

  // Free Vectors
  dcsr_free(&Alap);
  dvec_free(&b);
  dvec_free(&p);
  dvec_free(&gradp);

  return;
}
/***********************************************************************************************/
