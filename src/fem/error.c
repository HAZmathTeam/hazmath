/*! \file src/fem/error.c
 *
 * \brief This code contains functions for computing errors and norms
 *        mostly using FEM assembly routines.
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 5/14/15.
 *  Copyright 2016__HAZMATH__. All rights reserved.
 *
 *  \note modified by James Adler 11/1/2016
 *  \note Updated on 9/26/2018 for 0-1 fix.
 *
 */

#include "hazmath.h"

/***************************************************************************/
/*!
 * \fn REAL L2norm(REAL *u,fespace *FE,mesh_struct *mesh,qcoordinates *cq)
 *
 * \brief Computes the L2 Norm of a FE approximation using the mass matrix
 *        assembly for any type of element.
 *
 * \param u 	    Numerical Solution at DOF
 * \param FE      FE Space
 * \param mesh    Mesh Data
 * \param cq      Quadrature Nodes
 *
 * \return norm   L2 Norm
 *
 */
REAL L2norm(REAL *u,fespace *FE,mesh_struct *mesh,qcoordinates *cq)
{

  INT i,j,k;
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

    // Find DOF for given Element
    get_incidence_row(i,FE->el_dof,dof_on_elm);

    //Find Nodes for given Element if not H1 elements
    get_incidence_row(i,mesh->el_v,v_on_elm);

    // Compute Local Stiffness Matrix for given Element
    assemble_mass_local(MLoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,constant_coeff_scal,1.0);

    for(j=0;j<dof_per_elm;j++) {
      for(k=0;k<dof_per_elm;k++) {
        sum+=u[dof_on_elm[j]]*MLoc[j*dof_per_elm+k]*u[dof_on_elm[k]];
      }
    }
  }

  if(MLoc) free(MLoc);
  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);

  return sqrt(sum);
}
/***************************************************************************/

/***************************************************************************/
/*!
 * \fn void L2norm_block(REAL *norm,REAL *u,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq)
 *
 * \brief Computes the L2 Norm of a block FE approximation using the mass matrix
 *        assembly for any type of element.
 *
 * \param u 	    Numerical Solution at DOF
 * \param FE      Block FE Space
 * \param mesh    Mesh Data
 * \param cq      Quadrature Nodes
 *
 * \return norm   L2 Norm
 *
 */
void L2norm_block(REAL *norm,REAL *u,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq)
{
  INT i;
  REAL* udof = u;

  for(i=0;i<FE->nspaces;i++) {
    norm[i] = L2norm(udof,FE->var_spaces[i],mesh,cq);
    udof += FE->var_spaces[i]->ndof;
  }

  return;
}
/***************************************************************************/

/***************************************************************************/
/*!
 * \fn REAL L2_InnerProduct(REAL *u,REAL *v,fespace *FE,mesh_struct *mesh,qcoordinates *cq)
 *
 * \brief Computes the L2 inner product of two FE approximations using the mass matrix
 *        assembly for any type of element.
 *
 * \param u 	    Numerical Solution 1 at DOF
 * \param v       Numerical Solution 2 at DOF
 * \param FE      FE Space
 * \param mesh    Mesh Data
 * \param cq      Quadrature Nodes
 *
 * \return product  L2 Inner Product of u and v, <u,v>
 *
 */
REAL L2_InnerProduct(REAL *u,REAL *v,fespace *FE,mesh_struct *mesh,qcoordinates *cq)
{
  INT i,j,k;
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

    // Find DOF for given Element
    get_incidence_row(i,FE->el_dof,dof_on_elm);

    //Find Nodes for given Element if not H1 elements
    get_incidence_row(i,mesh->el_v,v_on_elm);

    // Compute Local Stiffness Matrix for given Element
    assemble_mass_local(MLoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,constant_coeff_scal,1.0);

    for(j=0;j<dof_per_elm;j++) {
      for(k=0;k<dof_per_elm;k++) {
        sum+=v[dof_on_elm[j]]*MLoc[j*dof_per_elm+k]*u[dof_on_elm[k]];
      }
    }
  }

  if(MLoc) free(MLoc);
  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);

  return sum;
}
/***************************************************************************/

/***************************************************************************/
/*!
 * \fn void L2_InnerProduct_block(REAL *prod,REAL *u,REAL *v,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq)
 *
 * \brief Computes the L2 inner product of two block FE approximations using the mass matrix
 *        assembly for any type of element.
 *
 * \param u 	    Numerical Solution 1 at DOF
 * \param v       Numerical Solution 2 at DOF
 * \param FE      Block FE Space
 * \param mesh    Mesh Data
 * \param cq      Quadrature Nodes
 *
 * \return product  L2 Inner Product of u and v, <u,v> (for each component of block FE space)
 *
 */
void L2_InnerProduct_block(REAL *prod,REAL *u,REAL *v,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq)
{
  INT i;
  REAL* udof = u;
  REAL* vdof = v;

  for(i=0;i<FE->nspaces;i++) {
    prod[i] = L2_InnerProduct(udof,vdof,FE->var_spaces[i],mesh,cq);
    udof += FE->var_spaces[i]->ndof;
    vdof += FE->var_spaces[i]->ndof;
  }

  return;
}
/***************************************************************************/

/***************************************************************************/
/*!
 * \fn REAL L2error(REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
 *
 * \brief Computes the L2 Norm of the error of a FE approximation and a true
 *        solution given by a function using quadrature for any type of element.
 *
 * \param u 	          Numerical Solution at DOF
 * \param truesol       Function to get true solution at a given point
 * \param FE            FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param time          Physical time to compute solution at
 *
 * \return error        L2 Error
 *
 */
REAL L2error(REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
{
  INT dim = mesh->dim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(mesh->dim,sizeof(REAL));

  // FE Stuff
  //INT FEtype = FE->FEtype;
  INT scal_or_vec = FE->scal_or_vec;
  INT elm,quad,j;
  REAL sum = 0.0;
  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  // True Solution and FE Solution at Quadrature Nodes
  INT ncomp = 0;
  if(scal_or_vec==0) { // Lagrange Elements (i.e., Scalar elements)
    ncomp = 1;
  } else { // Vector Elements
    ncomp = dim;
  }
  REAL* val_true = (REAL *) calloc(ncomp,sizeof(REAL));
  REAL* val_sol = (REAL *) calloc(ncomp,sizeof(REAL));

  /* Loop over all Elements */
  for (elm=0; elm<FE->nelm; elm++) {

    // Find DOF for given Element
    get_incidence_row(elm,FE->el_dof,dof_on_elm);

    // Find Vertices for given Element if not H1 elements
    get_incidence_row(elm,mesh->el_v,v_on_elm);

    // Loop over quadrature nodes on element
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(mesh->dim==2 || mesh->dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];

      // Get True Solution at Quadrature Nodes
      (*truesol)(val_true,qx,time,&(mesh->el_flag[elm]));

      // Interpolate FE solution to quadrature point
      FE_Interpolation(val_sol,u,qx,dof_on_elm,v_on_elm,FE,mesh);

      // Compute Error on Element
      for(j=0;j<ncomp;j++) {
        sum+=w*(ABS(val_sol[j] - val_true[j]))*(ABS(val_sol[j] - val_true[j]));
      }
    }
  }


  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(qx) free(qx);
  if(val_true) free(val_true);
  if(val_sol) free(val_sol);

  return sqrt(sum);
}
/****************************************************************************/

/***************************************************************************/
/*!
 * \fn void L2error_block(REAL *err,REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
 *
 * \brief Computes the L2 Norm of the error of a block FE approximation and a true
 *        solution given by a function using quadrature for any type of element.
 *
 * \param u 	          Numerical Solution at DOF
 * \param truesol       Function to get true solution at a given point
 * \param FE            FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param time          Physical time to compute solution at
 *
 * \return err          L2 Error
 *
 */
void L2error_block(REAL *err,REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
{
  // Loop Indices
  INT i,elm,quad,j,rowa,rowb,jcntr;

  // Mesh Stuff
  INT dim = mesh->dim;
  INT v_per_elm = mesh->v_per_elm;
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));

  // FEM Stuff
  INT nspaces = FE->nspaces;
  INT dof_per_elm = 0;
  INT* ncomp = (INT *) calloc(FE->nspaces,sizeof(INT));
  INT nun=0;
  for(i=0;i<FE->nspaces;i++) {
    err[i]=0;
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
    if(FE->var_spaces[i]->scal_or_vec==0) /* Scalar Element */
      ncomp[i]=1;
    else /* Vector Element */
      ncomp[i] = dim;
    nun += ncomp[i];
  }
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  REAL* val_true = (REAL *) calloc(nun,sizeof(REAL));
  REAL* val_sol = (REAL *) calloc(nun,sizeof(REAL));

  /* Loop over all Elements */
  for (elm=0; elm<mesh->nelm; elm++) {

    // Find DOF for given Element
    // Note this is "local" ordering for the given FE space of the block
    // Not global ordering of all DOF
    jcntr = 0;
    for(i=0;i<nspaces;i++) {
      rowa = FE->var_spaces[i]->el_dof->IA[elm];
      rowb = FE->var_spaces[i]->el_dof->IA[elm+1];
      for (j=rowa; j<rowb; j++) {
        dof_on_elm[jcntr] = FE->var_spaces[i]->el_dof->JA[j];
        jcntr++;
      }
    }
    // Find vertices for given Element
    get_incidence_row(elm,mesh->el_v,v_on_elm);

    // Loop over quadrature nodes on element
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(mesh->dim==2 || mesh->dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];

      // Get True Solution at Quadrature Nodes
      (*truesol)(val_true,qx,time,&(mesh->el_flag[elm]));

      // Interpolate FE solution to quadrature point
      blockFE_Interpolation(val_sol,u,qx,dof_on_elm,v_on_elm,FE,mesh);

      // Compute Square of Error on Element for each component of FE space
      jcntr=0;
      for(i=0;i<nspaces;i++) {
        for(j=0;j<ncomp[i];j++) {
          err[i]+=w*(ABS(val_sol[jcntr+j] - val_true[jcntr+j]))*(ABS(val_sol[jcntr+j] - val_true[jcntr+j]));
        }
        jcntr+=ncomp[i];
      }
    }
  }

  for(i=0;i<nspaces;i++) {
    err[i] = sqrt(err[i]);
  }

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(qx) free(qx);
  if(val_true) free(val_true);
  if(val_sol) free(val_sol);
  if(ncomp) free(ncomp);

  return;
}
/************************************************************************************************/

/***************************************************************************/
/*!
   * \fn L2error_mass(REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
   *
   * \brief Computes the L2 Norm of the error of a FE approximation and a true
   *        solution given by a function using mass matrix assembly for any type of element.
   *
   * \param u 	          Numerical Solution at DOF
   * \param truesol       Function to get true solution at a given point
   * \param FE            FE Space
   * \param mesh          Mesh Data
   * \param cq            Quadrature Nodes
   * \param time          Physical time to compute solution at
   *
   * \return error        L2 Error
   *
   */
REAL L2error_mass(REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
{
  INT i,j,k;
  REAL sum = 0.0;
  REAL utk,utj,erk,erj;

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

    // Find DOF for given Element
    get_incidence_row(i,FE->el_dof,dof_on_elm);

    //Find Nodes for given Element if not H1 elements
    get_incidence_row(i,mesh->el_v,v_on_elm);

    // Compute Local Stiffness Matrix for given Element
    assemble_mass_local(MLoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,constant_coeff_scal,1.0);

    for(j=0;j<dof_per_elm;j++) {
      for(k=0;k<dof_per_elm;k++) {
        utj = FE_Evaluate_DOF(truesol,FE,mesh,time,dof_on_elm[j]);
        utk = FE_Evaluate_DOF(truesol,FE,mesh,time,dof_on_elm[k]);
        erj = (utj - u[dof_on_elm[j]]);
        erk = (utk - u[dof_on_elm[k]]);
        sum+=erj*MLoc[j*dof_per_elm+k]*erk;
      }
    }
  }

  if(MLoc) free(MLoc);
  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);

  return sqrt(ABS(sum));
}
/*******************************************************************************************************************************************************/

/***************************************************************************/
/*!
 * \fn void L2error_block_mass(REAL *err, REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
 *
 * \brief Computes the L2 Norm of the error of a block FE approximation and a true
 *        solution given by a function using mass matrix assembly for any type of element.
 *
 * \param u 	          Numerical Solution at DOF in blocks
 * \param truesol       Function to get true solution at a given point
 * \param FE            Block FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param time          Physical time to compute solution at
 *
 * \return err          L2 Error
 *
 */
void L2error_block_mass(REAL *err, REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
{
  INT i,j,k,elm;
  REAL utk,utj,erk,erj;

  INT v_per_elm = mesh->v_per_elm;
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  INT local_size=0;
  INT u_dof = 0;
  INT nspaces = FE->nspaces;

  for(i=0;i<nspaces;i++) {
    err[i] = 0.0;
  }

  /* Loop over all Elements */
  for (elm=0; elm<mesh->nelm; elm++) {

    // Find vertices for given Element
    get_incidence_row(elm,mesh->el_v,v_on_elm);

    // Get DOF and error on DOF for given element for each FE space
    u_dof=0;
    for(i=0;i<nspaces;i++) {
      local_size = FE->var_spaces[i]->dof_per_elm;
      INT* dof_on_elm = (INT *) calloc(local_size,sizeof(INT));
      get_incidence_row(elm,FE->var_spaces[i]->el_dof,dof_on_elm);

      REAL* MLoc = calloc(local_size*local_size,sizeof(REAL));
      for (j=0; j<local_size*local_size; j++) MLoc[j] = 0.0;

      // Compute Local Stiffness Matrix for given Element
      assemble_mass_local(MLoc,FE->var_spaces[i],mesh,cq,dof_on_elm,v_on_elm,elm,constant_coeff_scal,1.0);

      for(j=0;j<local_size;j++) {
        for(k=0;k<local_size;k++) {
          utj = blockFE_Evaluate_DOF(truesol,FE,mesh,time,i,dof_on_elm[j]);
          utk = blockFE_Evaluate_DOF(truesol,FE,mesh,time,i,dof_on_elm[k]);
          erj = (utj - u[u_dof + dof_on_elm[j]]);
          erk = (utk - u[u_dof + dof_on_elm[k]]);
          err[i]+=erj*MLoc[j*local_size+k]*erk;
        }
      }
      u_dof += FE->var_spaces[i]->ndof;
      if(MLoc) free(MLoc);
      if(dof_on_elm) free(dof_on_elm);
    }
  }

  for(i=0;i<nspaces;i++) {
    err[i] = sqrt(ABS(err[i]));
  }

  if(v_on_elm) free(v_on_elm);

  return;
}
/*******************************************************************************************************************************************************/

/***************************************************************************/
/*!
 * \fn REAL HDseminorm(REAL *u,fespace *FE,mesh_struct *mesh,qcoordinates *cq)
 *
 * \brief Computes the H(D) semi-Norm of a FE approximation using the
 *        <Du,Dv> matrix assembly for any type of element.
 *          Nodal   - <grad u, grad u>  -> |u|_1
 *          RT      - <div u, div u>    -> |u|_(H(div))
 *          Nedelec - <curl u, curl u>  -> |u|_(H(curl))
 *
 * \param u 	    Numerical Solution at DOF
 * \param FE      FE Space
 * \param mesh    Mesh Data
 * \param cq      Quadrature Nodes
 *
 * \return norm   HD Semi Norm
 *
 * \note Return 0.0 if P0 elements are used.
 *
 */
REAL HDseminorm(REAL *u,fespace *FE,mesh_struct *mesh,qcoordinates *cq)
{
  INT i,j,k;
  REAL sum = 0.0;
  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* ALoc = calloc(local_size,sizeof(REAL));
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  if(FE->FEtype==0 || FE->FEtype==99) {
    sum=0.0;
  } else {
    /* Loop over all Elements */
    for (i=0; i<FE->nelm; i++) {

      // Zero out local matrices
      for (j=0; j<local_size; j++) ALoc[j] = 0.0;

      // Find DOF for given Element
      get_incidence_row(i,FE->el_dof,dof_on_elm);

      //Find Nodes for given Element if not H1 elements
      get_incidence_row(i,mesh->el_v,v_on_elm);

      // Compute Local Stiffness Matrix for given Element
      assemble_DuDv_local(ALoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,constant_coeff_scal,1.0);

      for(j=0;j<dof_per_elm;j++) {
        for(k=0;k<dof_per_elm;k++) {
          sum+=u[dof_on_elm[j]]*ALoc[j*dof_per_elm+k]*u[dof_on_elm[k]];
        }
      }
    }
  }

  if(ALoc) free(ALoc);
  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);

  if(sum<0.0) {
    printf("Your H1 Semi Norm Squared is negative (%25.16e)!  Taking ABS before squarerooting itself\n",sum);
    return sqrt(ABS(sum));
  } else {
    return sqrt(sum);
  }

}
/*******************************************************************************************************************************************************/

/***************************************************************************/
/*!
 * \fn void HDseminorm_block(REAL *norm, REAL *u,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq)
 *
 * \brief Computes the H(D) semi-Norm of a block FE approximation using the
 *        <Du,Dv> matrix assembly for any type of element.
 *          Nodal   - <grad u, grad u>  -> |u|_1
 *          RT      - <div u, div u>    -> |u|_(H(div))
 *          Nedelec - <curl u, curl u>  -> |u|_(H(curl))
 *
 * \param u 	    Numerical Solution at DOF
 * \param FE      Block FE Space
 * \param mesh    Mesh Data
 * \param cq      Quadrature Nodes
 *
 * \return norm   HD Semi Norm
 *
 * \note Return 0.0 if P0 elements are used for that block.
 *
 */
void HDseminorm_block(REAL *norm,REAL *u,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq)
{
  INT i;
  REAL* udof = u;

  for(i=0;i<FE->nspaces;i++) {
    norm[i] = HDseminorm(udof,FE->var_spaces[i],mesh,cq);
    udof += FE->var_spaces[i]->ndof;
  }

  return;
}
/*******************************************************************************************************************************************************/

/***************************************************************************/
/*!
 * \fn REAL HDsemierror(REAL *u,void (*D_truesol)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
 *
 * \brief Computes the H(D) semi-norm of the error of a FE approximation and a true
 *        solution given by a function using quadrature for any type of element.
 *          Nodal   - <grad u, grad u>  -> |u|_1
 *          RT      - <div u, div u>    -> |u|_(H(div))
 *          Nedelec - <curl u, curl u>  -> |u|_(H(curl))
 * \param u 	          Numerical Solution at DOF
 * \param D_truesol     Function to get derivative of true solution at a given point
 * \param FE            FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param time          Physical time to compute solution at
 *
 * \return error        Semi-Norm Error
 *
 */
REAL HDsemierror(REAL *u,void (*D_truesol)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
{
  INT dim = mesh->dim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(mesh->dim,sizeof(REAL));

  // FE Stuff
  INT FEtype = FE->FEtype;
  INT elm,quad,j;
  REAL sum = 0.0;
  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  // Derivative of True Solution and FE Solution at Quadrature Nodes
  INT ncomp = 0;
  if(FEtype<20) { // Lagrange Elements -> Grad
    ncomp = dim;
  } else if(FEtype>=20 && FEtype<30 && dim==3) { // Nedelec Elements in 3D -> Curl is a Vector
    ncomp = dim;
  } else { // 2D Nedelec -> Curl is a scalar or RT -> Div is a scalar
    ncomp = 1;
  }
  REAL* val_true = (REAL *) calloc(ncomp,sizeof(REAL));
  REAL* val_sol = (REAL *) calloc(ncomp,sizeof(REAL));

  /* Loop over all Elements */
  for (elm=0; elm<FE->nelm; elm++) {

    // Find DOF for given Element
    get_incidence_row(elm,FE->el_dof,dof_on_elm);

    //Find Vertices for given Element if not H1 elements
    get_incidence_row(elm,mesh->el_v,v_on_elm);

    // Loop over quadrature nodes on element

    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(mesh->dim==2 || mesh->dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];

      // Get True Solution at Quadrature Nodes
      (*D_truesol)(val_true,qx,time,&(mesh->el_flag[elm]));

      // Interpolate FE solution to quadrature point
      FE_DerivativeInterpolation(val_sol,u,qx,dof_on_elm,v_on_elm,FE,mesh);

      // Compute Error on Element
      for(j=0;j<ncomp;j++) {
        sum+=w*(ABS(val_sol[j] - val_true[j]))*(ABS(val_sol[j] - val_true[j]));
      }
    }
  }


  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(qx) free(qx);
  if(val_true) free(val_true);
  if(val_sol) free(val_sol);

  return sqrt(sum);
}
/*******************************************************************************************************************************************************/

/***************************************************************************/
/*!
 * \fn void HDsemierror_block(REAL *err,REAL *u,void (*D_truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
 *
 * \brief Computes the H(D) semi-norm of the error of a block FE approximation and a true
 *        solution given by a function using quadrature for any type of element.
 *          Nodal   - <grad u, grad u>  -> |u|_1
 *          RT      - <div u, div u>    -> |u|_(H(div))
 *          Nedelec - <curl u, curl u>  -> |u|_(H(curl))
 * \param u 	          Numerical Solution at DOF
 * \param D_truesol     Function to get derivative of true solution at a given point
 * \param FE            Block FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param time          Physical time to compute solution at
 *
 * \return err          Semi-Norm Error
 *
 */
void HDsemierror_block(REAL *err,REAL *u,void (*D_truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
{
  // Loop Indices
  INT i,elm,quad,j,rowa,rowb,jcntr;

  // Mesh Stuff
  INT dim = mesh->dim;
  INT v_per_elm = mesh->v_per_elm;
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));

  // FEM Stuff
  INT nspaces = FE->nspaces;
  INT dof_per_elm = 0;
  INT* ncomp = (INT *) calloc(FE->nspaces,sizeof(INT));
  INT nun=0;
  for(i=0;i<FE->nspaces;i++) {
    err[i] = 0.0;
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
    if(FE->var_spaces[i]->FEtype<20 || FE->var_spaces[i]->FEtype==103) /* Scalar Gradient */
      ncomp[i]=dim;
    else if(FE->var_spaces[i]->FEtype==20 && dim==2) /* Curl in 2D is scalar */
      ncomp[i] = 1;
    else if(FE->var_spaces[i]->FEtype==20 && dim==3) /* Curl in 3D is vector */
      ncomp[i] = dim;
    else /* Div is scalar */
      ncomp[i] = 1;
    nun += ncomp[i];
  }
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  REAL* val_true = (REAL *) calloc(nun,sizeof(REAL));
  REAL* val_sol = (REAL *) calloc(nun,sizeof(REAL));

  /* Loop over all Elements */
  for (elm=0; elm<mesh->nelm; elm++) {

    // Find DOF for given Element
    // Note this is "local" ordering for the given FE space of the block
    // Not global ordering of all DOF
    jcntr = 0;
    for(i=0;i<nspaces;i++) {
      rowa = FE->var_spaces[i]->el_dof->IA[elm];
      rowb = FE->var_spaces[i]->el_dof->IA[elm+1];
      for (j=rowa; j<rowb; j++) {
        dof_on_elm[jcntr] = FE->var_spaces[i]->el_dof->JA[j];
        jcntr++;
      }
    }
    // Find vertices for given Element
    get_incidence_row(elm,mesh->el_v,v_on_elm);

    // Loop over quadrature nodes on element
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(mesh->dim==2 || mesh->dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];

      // Get True Solution at Quadrature Nodes
      (*D_truesol)(val_true,qx,time,&(mesh->el_flag[elm]));

      // Interpolate FE solution to quadrature point
      blockFE_DerivativeInterpolation(val_sol,u,qx,dof_on_elm,v_on_elm,FE,mesh);

      // Compute Square of Error on Element for each component of FE space
      jcntr=0;
      for(i=0;i<nspaces;i++) {
        for(j=0;j<ncomp[i];j++) {
          err[i]+=w*(ABS(val_sol[jcntr+j] - val_true[jcntr+j]))*(ABS(val_sol[jcntr+j] - val_true[jcntr+j]));
        }
        jcntr+=ncomp[i];
      }
    }
  }

  for(i=0;i<nspaces;i++) {
    if(FE->var_spaces[i]->FEtype==0) { // Don't get H1 norm for P0 Elements
      err[i] = 0.0;
    } else {
      err[i] = sqrt(err[i]);
    }
  }

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(qx) free(qx);
  if(val_true) free(val_true);
  if(val_sol) free(val_sol);
  if(ncomp) free(ncomp);

  return;
}
/*******************************************************************************************************************************************************/

/***************************************************************************/
/*!
 * \fn REAL HDsemierror_stiff(REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
 *
 * \brief Computes the H(D) semi-norm of the error of a FE approximation and a true
 *        solution given by a function using the <Du,Dv> matrix assembly for any type of element.
 *          Nodal   - <grad u, grad u>  -> |u|_1
 *          RT      - <div u, div u>    -> |u|_(H(div))
 *          Nedelec - <curl u, curl u>  -> |u|_(H(curl))
 * \param u 	          Numerical Solution at DOF
 * \param D_truesol     Function to get derivative of true solution at a given point
 * \param FE            FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param time          Physical time to compute solution at
 *
 * \return error        Semi-Norm Error
 *
 */
REAL HDsemierror_stiff(REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
{
  INT i,j,k;
  REAL sum = 0.0;
  REAL utk,utj,erk,erj;

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* ALoc = calloc(local_size,sizeof(REAL));
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  /* Loop over all Elements */
  for (i=0; i<FE->nelm; i++) {

    // Zero out local matrices
    for (j=0; j<local_size; j++) ALoc[j] = 0.0;

    // Find DOF for given Element
    get_incidence_row(i,FE->el_dof,dof_on_elm);

    //Find Nodes for given Element if not H1 elements
    get_incidence_row(i,mesh->el_v,v_on_elm);

    // Compute Local Stiffness Matrix for given Element
    assemble_DuDv_local(ALoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,constant_coeff_scal,1.0);

    for(j=0;j<dof_per_elm;j++) {
      for(k=0;k<dof_per_elm;k++) {
        utj = FE_Evaluate_DOF(truesol,FE,mesh,time,dof_on_elm[j]);
        utk = FE_Evaluate_DOF(truesol,FE,mesh,time,dof_on_elm[k]);
        erj = (utj - u[dof_on_elm[j]]);
        erk = (utk - u[dof_on_elm[k]]);
        sum+=erj*ALoc[j*dof_per_elm+k]*erk;
      }
    }
  }

  if(ALoc) free(ALoc);
  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);

  return sqrt(ABS(sum));
}
/*******************************************************************************************************************************************************/

/***************************************************************************/
/*!
 * \fn void HDsemierror_block_stiff(REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
 *
 * \brief Computes the H(D) semi-norm of the error of a block FE approximation and a true
 *        solution given by a function using the <Du,Dv> matrix assembly for any type of element.
 *          Nodal   - <grad u, grad u>  -> |u|_1
 *          RT      - <div u, div u>    -> |u|_(H(div))
 *          Nedelec - <curl u, curl u>  -> |u|_(H(curl))
 * \param u 	          Numerical Solution at DOF
 * \param D_truesol     Function to get derivative of true solution at a given point
 * \param FE            Block FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param time          Physical time to compute solution at
 *
 * \return err          Semi-Norm Error
 *
 */
void HDsemierror_block_stiff(REAL *err, REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
{
  INT i,j,k,elm;
  REAL utk,utj,erk,erj;

  INT v_per_elm = mesh->v_per_elm;
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  INT local_size=0;
  INT u_dof = 0;
  INT nspaces = FE->nspaces;

  for(i=0;i<nspaces;i++) {
    err[i] = 0.0;
  }

  /* Loop over all Elements */
  for (elm=0; elm<mesh->nelm; elm++) {

    // Find vertices for given Element
    get_incidence_row(elm,mesh->el_v,v_on_elm);

    // Get DOF and error on DOF for given element for each FE space
    u_dof=0;
    for(i=0;i<nspaces;i++) {
      local_size = FE->var_spaces[i]->dof_per_elm;
      INT* dof_on_elm = (INT *) calloc(local_size,sizeof(INT));
      get_incidence_row(elm,FE->var_spaces[i]->el_dof,dof_on_elm);

      REAL* ALoc = calloc(local_size*local_size,sizeof(REAL));
      for (j=0; j<local_size*local_size; j++) ALoc[j] = 0.0;

      // Compute Local Stiffness Matrix for given Element
      assemble_DuDv_local(ALoc,FE->var_spaces[i],mesh,cq,dof_on_elm,v_on_elm,elm,constant_coeff_scal,1.0);

      for(j=0;j<local_size;j++) {
        for(k=0;k<local_size;k++) {
          utj = blockFE_Evaluate_DOF(truesol,FE,mesh,time,i,dof_on_elm[j]);
          utk = blockFE_Evaluate_DOF(truesol,FE,mesh,time,i,dof_on_elm[k]);
          erj = (utj - u[u_dof + dof_on_elm[j]]);
          erk = (utk - u[u_dof + dof_on_elm[k]]);
          err[i]+=erj*ALoc[j*local_size+k]*erk;
        }
      }
      u_dof += FE->var_spaces[i]->ndof;
      if(ALoc) free(ALoc);
      if(dof_on_elm) free(dof_on_elm);
    }
  }

  if(v_on_elm) free(v_on_elm);

  for(i=0;i<nspaces;i++) {
    err[i] = sqrt(ABS(err[i]));
  }

  return;
}
/*******************************************************************************************************************************************************/

/***************************************************************************/
/*!
 * \fn REAL HDnorm(REAL *u,fespace *FE,mesh_struct *mesh,qcoordinates *cq)
 *
 * \brief Computes the H(D) norm of a FE approximation using the
 *        <Du,Dv> matrix assembly for any type of element.
 *          Nodal   - <u,u> + <grad u, grad u>  -> ||u||_1
 *          RT      - <u,u> + <div u, div u>    -> ||u||_(H(div))
 *          Nedelec - <u,u> + <curl u, curl u>  -> ||u||_(H(curl))
 *
 * \param u 	    Numerical Solution at DOF
 * \param FE      FE Space
 * \param mesh    Mesh Data
 * \param cq      Quadrature Nodes
 *
 * \return norm   HD Norm
 *
 * \note Return just L2 norm if P0 elements are used.
 *
 */
REAL HDnorm(REAL *u,fespace *FE,mesh_struct *mesh,qcoordinates *cq)
{
  REAL sumL2 = L2norm(u,FE,mesh,cq);
  REAL sumSemi = HDseminorm(u,FE,mesh,cq);

  return sqrt(sumL2*sumL2 + sumSemi*sumSemi);

}
/*******************************************************************************************************************************************************/

/***************************************************************************/
/*!
 * \fn void HDnorm(REAL *norm, REAL *u,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq)
 *
 * \brief Computes the H(D) norm of a block FE approximation using the
 *        <Du,Dv> matrix assembly for any type of element.
 *          Nodal   - <u,u> + <grad u, grad u>  -> ||u||_1
 *          RT      - <u,u> + <div u, div u>    -> ||u||_(H(div))
 *          Nedelec - <u,u> + <curl u, curl u>  -> ||u||_(H(curl))
 *
 * \param u 	    Numerical Solution at DOF
 * \param FE      Block FE Space
 * \param mesh    Mesh Data
 * \param cq      Quadrature Nodes
 *
 * \return norm   HD Norm
 *
 * \note Return just L2 norm if P0 elements are used for that block.
 *
 */
void HDnorm_block(REAL *norm,REAL *u,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq)
{
  INT i;
  REAL* sumL2 = (REAL *) calloc(FE->nspaces,sizeof(REAL));
  REAL* sumSemi = (REAL *) calloc(FE->nspaces,sizeof(REAL));

  L2norm_block(sumL2,u,FE,mesh,cq);
  HDseminorm_block(sumSemi,u,FE,mesh,cq);

  for(i=0;i<FE->nspaces;i++)
    norm[i] = sqrt(sumL2[i]*sumL2[i] + sumSemi[i]*sumSemi[i]);

  if(sumL2) free(sumL2);
  if(sumSemi) free(sumSemi);

  return;
}
/***************************************************************************/

/***************************************************************************/
/*!
 * \fn REAL HDerror(REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),void (*D_truesol)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
 *
 * \brief Computes the H(D) norm of the error of a FE approximation
 *        and a true solution given by a function using the
 *        <u,v> and <Du,Dv> matrix assembly for any type of element.
 *          Nodal   - <u,u> + <grad u, grad u>  -> ||u||_1
 *          RT      - <u,u> + <div u, div u>    -> ||u||_(H(div))
 *          Nedelec - <u,u> + <curl u, curl u>  -> ||u||_(H(curl))
 *
 * \param u 	          Numerical Solution at DOF
 * \param truesol       Function to get true solution at given point
 * \param D_truesol     Function to get derivative of true solution at a given point
 * \param FE            FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param time          Physical time to compute solution at
 *
 * \return norm         HD Norm of Error
 *
 */
REAL HDerror(REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),void (*D_truesol)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
{
  REAL sumL2 = L2error(u,truesol,FE,mesh,cq,time);
  REAL sumSemi = HDsemierror(u,D_truesol,FE,mesh,cq,time);

  return sqrt(sumL2*sumL2 + sumSemi*sumSemi);

}
/*******************************************************************************************************************************************************/

/***************************************************************************/
/*!
 * \fn void HDerror_block(REAL *err,REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),void (*D_truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
 *
 * \brief Computes the H(D) norm of the error of a block FE approximation
 *        and a true solution given by a function using the
 *        <u,v> and <Du,Dv> matrix assembly for any type of element.
 *          Nodal   - <u,u> + <grad u, grad u>  -> ||u||_1
 *          RT      - <u,u> + <div u, div u>    -> ||u||_(H(div))
 *          Nedelec - <u,u> + <curl u, curl u>  -> ||u||_(H(curl))
 *
 * \param u 	          Numerical Solution at DOF
 * \param truesol       Function to get true solution at given point
 * \param D_truesol     Function to get derivative of true solution at a given point
 * \param FE            Block FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param time          Physical time to compute solution at
 *
 * \return norm         HD Norm of Error
 *
 */
void HDerror_block(REAL *err,REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),void (*D_truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
{
  INT i;
  REAL* sumL2 = (REAL *) calloc(FE->nspaces,sizeof(REAL));
  REAL* sumSemi = (REAL *) calloc(FE->nspaces,sizeof(REAL));

  L2error_block(sumL2,u,truesol,FE,mesh,cq,time);
  HDsemierror_block(sumSemi,u,D_truesol,FE,mesh,cq,time);
  for(i=0;i<FE->nspaces;i++)  {
    err[i] = sqrt(sumL2[i]*sumL2[i] + sumSemi[i]*sumSemi[i]);
  }
  if(sumL2) free(sumL2);
  if(sumSemi) free(sumSemi);

  return;
}
/*******************************************************************************************************************************************************/
/***************************************************************************/
/*!
 * \fn REAL energynorm_discrete(REAL *u,fespace *FE,mesh_struct *mesh,qcoordinates *cq, void (*coeff)(REAL *,REAL *,REAL,void *),void (*local_assembly_routine)(REAL *,fespace *,mesh_struct *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL,void *),REAL), REAL param)
 *
 * \brief Computes the matrix energy norm of a FE approximation using the local matrix
 *        assembly for any type of element.
 *
 * \param u 	    				Numerical Solution at DOF
 * \param FE      				FE Space
 * \param mesh    				Mesh Data
 * \param cq     					Quadrature Nodes
 * \param coeff   				PDE coefficients
 * \param local_assembly_routine	local assembly routine for LHS matrix (stiffness/mass matrix)
 * \param param						extra real param needed by local assemble (could be time)
 *
 * \return norm   				energy norm
 *
 */
REAL energynorm_discrete(REAL *u,fespace *FE,mesh_struct *mesh,qcoordinates *cq, void (*coeff)(REAL *,REAL *,REAL,void *),void (*local_assembly_routine)(REAL *,fespace *,mesh_struct *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL,void *),REAL), REAL param)
{

  INT i,j,k;
  REAL sum = 0.0;

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* ALoc = calloc(local_size,sizeof(REAL));
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  // Loop over all Elements
  for (i=0; i<FE->nelm; i++) {

    // Zero out local matrices
    for (j=0; j<local_size; j++) ALoc[j] = 0.0;

    // Find DOF for given Element
    get_incidence_row(i,FE->el_dof,dof_on_elm);

    // Find Nodes for given Element if not H1 elements
    get_incidence_row(i,mesh->el_v,v_on_elm);

    // Compute Local Stiffness Matrix(operator associated with energy norm) for given Element
    (*local_assembly_routine)(ALoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,coeff,param);

	// Loop over DOF
    for(j=0;j<dof_per_elm;j++) {
      for(k=0;k<dof_per_elm;k++) { //compute u^T*A*u
        sum+=u[dof_on_elm[j]]*ALoc[j*dof_per_elm+k]*u[dof_on_elm[k]];
      }
    }
  }

  if(ALoc) free(ALoc);
  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);

  return sqrt(sum);
}

/***************************************************************************/
