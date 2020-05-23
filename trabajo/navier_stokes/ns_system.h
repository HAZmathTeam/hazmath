/*! \file examples/navier_stokes/ns_system.h
 *
 *  Created by James Adler on 05/20/20.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This contains all the local assembly routines
 *        for the Navier-Stokes example.
 */

 /*!
 * \fn void local_assembly_NavierStokes(REAL *ALoc,REAL* bLoc, dvector *old_sol, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
 *
 * \brief Computes the local stiffness matrix (Jacobian) for the linearized
 *        steady-state portion of the Navier-Stokes system.  It also includes
 *        the nonlinear residual as the right-hand side.
 *
 *
 *        <2 eps(u), eps(v)> + Re*<u_old*grad u + u*grad u_old, v> - <p, div v> = L1(u_old,v)
 *                   - <div u, q> = <div u_old, q>
 *
 *        where eps(u) = (grad u + (grad u)^T)/2 is the symmetric gradient,
 *              L1(u_old, v) = <f,v> - <2 eps(u_old), eps(v)> - Re*<u_old*grad u_old, v> + <p_old, div v>
 *
 * \note Note the scaling of Reynolds Number here.  We just put this on the
 *       nonlinear inertia term.  We assume the pressure, p, contains the appropriate
 *       scaling.  This allows for the case of Re=0 equaling Stokes' Equations.
 *
 * \param old_sol       FE approximation of previous Newton step
 * \param FE            Block FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param dof_on_elm    Specific DOF on element
 * \param v_on_elm      Specific vertices on element
 * \param elm           Current element
 * \param rhs           Function for rhs (NULL here)
 * \param time          Physical Time if time dependent
 *
 * \return ALoc         Local Stiffness Matrix (Full Matrix) ordered (u1,u2,p)
 * \return bLoc         Local RHS ordered (u1,u2,p)
 *
 * \note Assumes 2D or 3D only
 *
 * \note We provide the computation for the 0-0 block:
 *       <2 eps(u), eps(v)> =
 *                <2 dx(u1),dx(v1)> + <dy(u1),dy(v1)>  +     <dx(u2),dy(v1)>
 *                <dy(u1),dx(v2)>                      +     <dx(u2),dx(v2)> + <2 dy(u2),dy(v2)>
 *
 *       for Dirichlet boundary conditions and when div u = 0 exactly, then
 *       <2 eps(u), eps(v)> = <grad u, grad v>
 *
 *
 */
 void local_assembly_NavierStokes(REAL *ALoc,REAL* bLoc, dvector *old_sol, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
 {

   // Loop indices
   INT i,j,quad,test,trial;

   // Mesh and FE data
   INT dim = mesh->dim;
   INT dof_per_elm = 0;
   for (i=0; i<FE->nspaces;i++)
     dof_per_elm += FE->var_spaces[i]->dof_per_elm;
   INT* local_dof_on_elm;
   // Total DoF for each unknwon and DoF per element for each unknown
   INT u1dof = FE->var_spaces[0]->ndof;
   INT u1dofpelm = FE->var_spaces[0]->dof_per_elm;
   INT u2dof = FE->var_spaces[1]->ndof;
   INT u2dofpelm = FE->var_spaces[1]->dof_per_elm;
   INT u3dof,u3dofpelm;
   if(dim==3) {
     INT u3dof = FE->var_spaces[2]->ndof;
     INT u3dofpelm = FE->var_spaces[2]->dof_per_elm;
   }
   INT pdof = FE->var_spaces[dim]->ndof;
   INT pdofpelm = FE->var_spaces[dim]->dof_per_elm;

   // Quadrature Weights and Nodes
   REAL w;
   REAL qx[dim];

   // Stiffness Matrix Entry
   REAL kij = 0.0;
   REAL rij = 0.0;

   // Keep track of local indexing
   INT local_row_index, local_col_index;

   // Stuff for previous solution
   REAL* local_uprev = NULL;
   REAL u1_prev;
   REAL u2_prev;
   REAL u3_prev;
   REAL p_prev;
   REAL* du1_prev = (REAL *) calloc( dim, sizeof(REAL));
   REAL* du2_prev = (REAL *) calloc( dim, sizeof(REAL));
   REAL* du3_prev = NULL;
   if(dim==3) du3_prev = (REAL *) calloc( dim, sizeof(REAL));
   // Divergence of u_prev
   REAL divuk;
   // Components of eps(u_prev) (symmetry means we only need 6 pieces)
   REAL epsu11,epsu21,epsu31,epsu22,epsu32,epsu33;
   // Components of u_prev*grad(u_prev)
   REAL ukgraduk1,ukgraduk2,ukgraduk3;

   // Trial Functions
   REAL u1,u2,u3,p,u1x,u1y,u1z,u2x,u2y,u2z,u3x,u3y,u3z;
   // Test Functions
   REAL v1,v2,v3,v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z,q;

   // Get Constants (assumes Reynolds number is constant for now)
   REAL Re = 0.0;
   get_reynolds_number(&Re);

   // Source Term Values if any
   REAL rhs_val[dim+1];

   // Sum over quadrature points
   for (quad=0;quad<cq->nq_per_elm;quad++) {
     qx[0] = cq->x[elm*cq->nq_per_elm+quad];
     qx[1] = cq->y[elm*cq->nq_per_elm+quad];
     if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
     w = cq->w[elm*cq->nq_per_elm+quad];
     (*rhs)(rhs_val,qx,time,&(mesh->el_flag[elm]));

     //  Get the Basis Functions and previous solutions at each quadrature node
     // u = (u1,u2) and v = (v1,v2)
     local_dof_on_elm = dof_on_elm;
     local_uprev = old_sol->val;
     FE_Interpolation(&u1_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[0],mesh);
     FE_DerivativeInterpolation(du1_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[0],mesh);
     get_FEM_basis(FE->var_spaces[0]->phi,FE->var_spaces[0]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[0]);

     local_dof_on_elm += u1dofpelm;
     local_uprev += u1dof;
     FE_Interpolation(&u2_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh);
     FE_DerivativeInterpolation(du2_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh);
     get_FEM_basis(FE->var_spaces[1]->phi,FE->var_spaces[1]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[1]);

     local_dof_on_elm += u2dofpelm;
     local_uprev += u2dof;
     if(dim==3) {
       FE_Interpolation(&u3_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[2],mesh);
       FE_DerivativeInterpolation(du3_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[2],mesh);
       get_FEM_basis(FE->var_spaces[2]->phi,FE->var_spaces[2]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[2]);
       local_dof_on_elm += u3dofpelm;
       local_uprev += u3dof;
     }

     // p
     FE_Interpolation(&p_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[dim],mesh);
     get_FEM_basis(FE->var_spaces[dim]->phi,FE->var_spaces[dim]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[dim]);

     // Some precomputations
     divuk = du1_prev[0] + du2_prev[1];
     epsu11 = du1_prev[0];
     epsu21 = 0.5*(du1_prev[1] + du2_prev[0]);
     epsu22 = du2_prev[1];
     ukgraduk1 = u1_prev*du1_prev[0] + u2_prev*du1_prev[1];
     ukgraduk2 = u1_prev*du2_prev[0] + u2_prev*du2_prev[1];
     if(dim==3) {
       divuk += du3_prev[2];
       epsu31 = 0.5*(du1_prev[2] + du3_prev[0]);
       epsu32 = 0.5*(du2_prev[2] + du3_prev[1]);
       ukgraduk1 += u3_prev*du1_prev[2];
       ukgraduk2 += u3_prev*du2_prev[2];
       ukgraduk3 = u1_prev*du3_prev[0] + u2_prev*du3_prev[1] + u3_prev*du3_prev[2];
     }

     // v1 block row
     local_row_index = 0;
     // Loop over Test Functions (Rows)
     for (test=0; test<u1dofpelm;test++){
       v1=FE->var_spaces[0]->phi[test];
       v1x=FE->var_spaces[0]->dphi[test*dim];
       v1y=FE->var_spaces[0]->dphi[test*dim+1];
       if(dim==3) v1z=FE->var_spaces[0]->dphi[test*dim+2];

       // u1-v1 block:
       // 2*<dx(u1),dx(v1)> + <dy(u1),dy(v1)> + <dz(u1),dz(v1)> + Re*<u1k*dx(u1)+u2k*dy(u1) + u1*dx(u1k),v>
       local_col_index = 0;
       // Loop over Trial Functions (Columns)
       for (trial=0; trial<u1dofpelm;trial++){
         u1=FE->var_spaces[0]->phi[trial];
         u1x=FE->var_spaces[0]->dphi[trial*dim];
         u1y=FE->var_spaces[0]->dphi[trial*dim+1];
         kij = 2*u1x*v1x + u1y*v1y + Re*(u1_prev*u1x+u2_prev*u1y+u1*du1_prev[0])*v1;
         if(dim==3) {
           u1z=FE->var_spaces[0]->dphi[trial*dim+2];
           kij+=u1z*v1z + Re*(u3_prev*u1z)*v1;
         }
         ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
       }

       // u2-v1 block:
       // <dx(u2),dy(v1)> + Re*<u2*dy(u1k),v1>
       local_col_index += u1dofpelm;
       // Loop over Trial Functions (Columns)
       for (trial=0; trial<u2dofpelm;trial++){
         u2=FE->var_spaces[1]->phi[trial];
         u2x=FE->var_spaces[1]->dphi[trial*dim];
         kij = u2x*v1y + Re*u2*du1_prev[1]*v1;
         ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
       }

       // u3-v1 block <dx(u3),dz(v1)> + Re*<u3*dz(u1k),v1>
       if(dim==3) {
         local_col_index += u2dofpelm;
         for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){
           u3=FE->var_spaces[2]->phi[trial];
           u3x=FE->var_spaces[2]->dphi[trial*dim];
           kij = u3x*v1z + Re*u3*du1_prev[2]*v1;
           ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
         }
       }

       // p-v1 block: -<p, dx(v1)>
       local_col_index += FE->var_spaces[dim-1]->dof_per_elm;
       // Loop over Trial Functions (Columns)
       for (trial=0; trial<pdofpelm;trial++){
         p = FE->var_spaces[dim]->phi[trial];
         kij = -p*v1x;
         ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
       }

       // v1 part of RHS
       rij = rhs_val[0]*v1 - 2*du1_prev[0]*v1x - 2*epsu21*v1y - Re*ukgraduk1*v1 + p_prev*v1x;
       if(dim==3) rij+= -2*epsu31*v1z;
       bLoc[(local_row_index+test)] += w*rij;
     }

     // v2 block row
     local_row_index += u1dofpelm;
     // Loop over Test Functions (Rows)
     for (test=0; test<u2dofpelm;test++){
       v2=FE->var_spaces[1]->phi[test];
       v2x=FE->var_spaces[1]->dphi[test*dim];
       v2y=FE->var_spaces[1]->dphi[test*dim+1];
       if(dim==3) v2z=FE->var_spaces[1]->dphi[test*dim+2];

       // u1-v2 block
       // <dy(u1),dx(v2)> + Re*<u1*dx(u2k),v2>
       local_col_index = 0;
       // Loop over Trial Functions (Columns)
       for (trial=0; trial<u1dofpelm;trial++){
         u1=FE->var_spaces[0]->phi[trial];
         u1y=FE->var_spaces[0]->dphi[trial*dim+1];
         kij = u1y*v2x + Re*u1*du2_prev[0]*v2;
         ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
       }

       // u2-v2 block:
       // <dx(u2),dx(v2)> + 2*<dy(u2),dy(v2)> + <dz(u2),dz(v2)> + Re*<u_prev*grad(u2) + u2*dy(u2k),v2>
       local_col_index += u1dofpelm;
       // Loop over Trial Functions (Columns)
       for (trial=0; trial<u2dofpelm;trial++){
         u2=FE->var_spaces[1]->phi[trial];
         u2x=FE->var_spaces[1]->dphi[trial*dim];
         u2y=FE->var_spaces[1]->dphi[trial*dim+1];
         kij = u2x*v2x + 2*u2y*v2y + Re*(u1_prev*u2x + u2_prev*u2y + u2*du2_prev[1])*v2;
         if(dim==3) {
           u1z=FE->var_spaces[1]->dphi[trial*dim+2];
           kij+=u2z*v2z + Re*(u3_prev*u2z)*v2;
         }
         ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
       }

       // u3-v2 block
       // <dy(u3),dz(v2)> + Re*<u3*dz(u2k),v2>
       if(dim==3) {
         local_col_index += u2dofpelm;
         for (trial=0; trial<u2dofpelm;trial++){
           u3=FE->var_spaces[2]->phi[trial];
           u3y=FE->var_spaces[2]->dphi[trial*dim+1];
           kij = u3y*v2z + Re*u3*du2_prev[2]*v2;
           ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
         }
       }

       // p-v2 block: -<p, dy(v2)>
       local_col_index += FE->var_spaces[dim-1]->dof_per_elm;
       // Loop over Trial Functions (Columns)
       for (trial=0; trial<pdofpelm;trial++){
         p = FE->var_spaces[dim]->phi[trial];
         kij = -p*v2y;
         ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
       }

       // v2 part of RHS
       rij = rhs_val[1]*v2 - 2*du2_prev[1]*v2y - 2*epsu21*v2x - Re*ukgraduk2*v2 + p_prev*v2y;
       if(dim==3) rij+= -2*epsu32*v2z;
       bLoc[(local_row_index+test)] += w*rij;
     }

     if(dim==3) {
       // v3 block row
       local_row_index += u2dofpelm;
       // Loop over Test Functions (Rows)
       for (test=0; test<u3dofpelm;test++){
         v3=FE->var_spaces[2]->phi[test];
         v3x=FE->var_spaces[2]->dphi[test*dim];
         v3y=FE->var_spaces[2]->dphi[test*dim+1];
         v3z=FE->var_spaces[2]->dphi[test*dim+2];

         // u1-v3 block
         // <dz(u1),dx(v3)> + Re*<u1*dx(u3k),v3>
         local_col_index = 0;
         // Loop over Trial Functions (Columns)
         for (trial=0; trial<u1dofpelm;trial++){
           u1=FE->var_spaces[0]->phi[trial];
           u1z=FE->var_spaces[0]->dphi[trial*dim+2];
           kij = u1z*v3x + Re*u1*du3_prev[0]*v3;
           ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
         }

         // u2-v3 block:
         // <dz(u2),dy(v3)> + Re*<u2*dy(u3k),v3>
         local_col_index += u1dofpelm;
         // Loop over Trial Functions (Columns)
         for (trial=0; trial<u2dofpelm;trial++){
           u2=FE->var_spaces[1]->phi[trial];
           u2z=FE->var_spaces[1]->dphi[trial*dim+2];
           kij = u2z*v3y + Re*u2*du3_prev[1]*v3;
           ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
         }

         // u3-v3 block
         // <dx(u3),dx(v3)> + <dy(u3),dy(v3)> + 2*<dz(u3),dz(v3)> + Re*<u_prev*grad(u3) + u3*dz(u3k),v3>
         local_col_index += u2dofpelm;
         for (trial=0; trial<u3dofpelm;trial++){
           u3=FE->var_spaces[2]->phi[trial];
           u3x=FE->var_spaces[2]->dphi[trial*dim];
           u3y=FE->var_spaces[2]->dphi[trial*dim+1];
           u3z=FE->var_spaces[2]->dphi[trial*dim+2];
           kij = u3x*v3x + u3y*v3y + 2*u3z*v3z + Re*(u1_prev*u3x + u2_prev*u3y + u3_prev*u3z + u3*du3_prev[2])*v3;
           ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
         }

         // p-v3 block: -<p, dz(v3)>
         local_col_index += u3dofpelm;
         // Loop over Trial Functions (Columns)
         for (trial=0; trial<pdofpelm;trial++){
           p = FE->var_spaces[dim]->phi[trial];
           kij = -p*v3z;
           ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
         }

         // v3 part of RHS
         rij = rhs_val[2]*v3 - 2*du3_prev[2]*v3z - 2*epsu31*v3x - 2*epsu32*v3y - Re*ukgraduk3*v3 + p_prev*v3z;
         bLoc[(local_row_index+test)] += w*rij;
       }
     }

     // q block row
     local_row_index += FE->var_spaces[dim-1]->dof_per_elm;
     // Loop over Test Functions (Rows)
     for (test=0; test<pdofpelm;test++){
       q=FE->var_spaces[dim]->phi[test];

       // u1-q block:
       // -<dx(u1), q>
       local_col_index = 0;
       // Loop over Trial Functions (Columns)
       for (trial=0; trial<u1dofpelm;trial++){
         u1x=FE->var_spaces[0]->dphi[trial*dim];
         kij = -u1x*q;
         ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
       }

       // u2-q block:
       // -<dy(u2), q>
       local_col_index += u1dofpelm;
       // Loop over Trial Functions (Columns)
       for (trial=0; trial<u2dofpelm;trial++){
         u2y=FE->var_spaces[1]->dphi[trial*dim+1];
         kij = -u2y*q;
         ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
       }

       if(dim==3) {
         // u3-q block
         // -<dz(u3), q>
         local_col_index += u2dofpelm;
         for (trial=0; trial<u3dofpelm;trial++){
           u3z=FE->var_spaces[2]->dphi[trial*dim+2];
           kij = -u3z*q;
           ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
         }
       }

       // p-q block: NULL

       // q part of RHS
       rij = divuk*q;
       bLoc[(local_row_index+test)] += w*rij;
     }
   }

   return;
 }
 /*****************************************************************************************************/
