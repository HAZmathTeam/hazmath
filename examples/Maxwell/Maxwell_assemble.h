/*! \file Maxwell_assemble.h
 *
 *  Created by Adler, Hu, Zikatanov on 10/1/17.
 *  Copyright 2016_HAZMATH__. All rights reserved.
 *
 * \brief This contains all the local assembly routine for the impedance boundary routine
 *        for the Maxwell example.
 *
 */

/******************************************************************************************************/
/*!
 * \fn void impedancebdry_local(REAL* ZLoc,dvector *old_sol,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *ed_on_f, \
                       INT *ed_on_elm,INT *v_on_elm,INT face,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time)
 *
 * \brief Computes the local weak formulation of the Impedance boundary condition for Maxwell's Equations
 *         Uses midpoint rule to integrate on edges of boundary face
 *         For this problem we compute the left-hand side of:
 *
 *         <n x E,n x F>_bdryobstacle    for all F in H_imp(curl) (Nedelec)
 *
 * \param old_sol       Solution at previous step (not needed in this function)
 * \param FE            FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param ed_on_f       Specific edges on the given face
 * \param ed_on_elm     Specific edges on the given element
 * \param v_on_elm      Specific vertices on the given element
 * \param face          Current face
 * \param elm           Current element
 * \param coeff         Function that gives coefficient (for now assume constant)
 * \param time          Physical Time if time dependent
 *
 * \return ZLoc         Local Boundary Matrix (Full Matrix)
 *
 * \note                ASSUMING 3D ONLY
 *
 */
void impedancebdry_local(REAL* ZLoc,dvector *old_sol,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *ed_on_f, \
                         INT *ed_on_elm,INT *v_on_elm,INT face,INT elm,void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
{
  // Mesh and FE data
  INT ed_per_elm = FE->dof_per_elm;
  INT dim = mesh->dim;

  // Loop Indices
  INT j,quad,test,trial,ed,edt,edb;

  // Quadrature Weights and Nodes
  INT nq = 2*dim-3; // = ed_per_face
  REAL* qx = (REAL *) calloc(nq,sizeof(REAL));
  // 3D: Using triangle midpoint rule, so qx is midpoint of edges and w is |F|/3
  REAL w = mesh->f_area[face]/3.0;

  // Get normal vector components on face
  REAL nx = mesh->f_norm[face*dim];
  REAL ny = mesh->f_norm[face*dim+1];
  REAL nz = mesh->f_norm[face*dim+2];

  // Stiffness Matrix Entry
  REAL kij,kij1,kij2,kij3,kij4,kij5,kij6;

  // Basis Functions and its curl
  REAL* phi= (REAL *) calloc(ed_per_elm*dim,sizeof(REAL));
  REAL* cphi = (REAL *) calloc(ed_per_elm*dim,sizeof(REAL));

  // Coefficient Value at Quadrature Nodes
  REAL coeff_val=0.0;

  //  Sum over midpoints of edges
  for (quad=0;quad<nq;quad++) {
    ed = ed_on_f[quad]-1;
    qx[0] = mesh->ed_mid[ed*dim];
    qx[1] = mesh->ed_mid[ed*dim+1];
    qx[2] = mesh->ed_mid[ed*dim+2];

    if(coeff!=NULL) {
      (*coeff)(&coeff_val,qx,time,&(mesh->ed_flag[ed]));
    } else {
      coeff_val = 1.0;
    }

    //  Get the Basis Functions at each quadrature node
    ned_basis(phi,cphi,qx,v_on_elm,ed_on_elm,mesh);

    // Loop over Test Functions (Rows - edges)
    for (test=0; test<nq;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<nq; trial++) {
        // Make sure ordering for global matrix is right
        for(j=0;j<mesh->ed_per_elm;j++) {
          if(ed_on_f[test]==ed_on_elm[j]) {
            edt = j;
          }
          if(ed_on_f[trial]==ed_on_elm[j]) {
            edb = j;
          }
        }
        kij1 = phi[edb*dim+1]*nz - phi[edb*dim+2]*ny;
        kij2 = phi[edb*dim+2]*nx - phi[edb*dim]*nz;
        kij3 = phi[edb*dim]*ny - phi[edb*dim+1]*nx;
        kij4 = phi[edt*dim+1]*nz - phi[edt*dim+2]*ny;
        kij5 = phi[edt*dim+2]*nx - phi[edt*dim]*nz;
        kij6 = phi[edt*dim]*ny - phi[edt*dim+1]*nx;
        kij = coeff_val*(kij1*kij4+kij2*kij5+kij3*kij6);
        ZLoc[test*nq+trial]+=w*kij;
      }
    }
  }

  if (phi) free(phi);
  if(cphi) free(cphi);
  if(qx) free(qx);
  return;
}
/******************************************************************************************************/
