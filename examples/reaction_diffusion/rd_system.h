/*! \file examples/reaction_diffusion/rd_system.h
*
*
* \brief This contains the local assembly routines for the reaction diffusion example:
*        -div A(x) grad u(x) + c(x) u(x) = f(x)
*
*/

/*!\fn void assemble_local_diffusion(REAL* ALoc,fespace *FE,mesh_struct *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
*
* \brief Computes the local stiffness matrix for the reaction-diffusion system.
*        For this problem we compute the local diffusion matrix
*
*	          ALoc = <a(x) grad phi_i, grad phi_j>
*
*	      of the problem
*
*			-div(a(x)grad(u)) + c(x)u = f(x)
*
*        where a(x) is a dim x dim matrix of scalar functions in x for dim=1,2,3
*
* \param FE            FE Space
* \param mesh          Mesh Data
* \param cq            Quadrature Nodes
* \param dof_on_elm    Specific DOF on element
* \param elm           Current element
* \param coeff         Function that gives coefficient
* \param time          Physical Time if time dependent
*
* \return ALoc         Local Stiffness Matrix (Full Matrix)
*
*
*
*/
void assemble_local_diffusion(REAL* ALoc,fespace *FE,mesh_struct *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
{
  INT dim = mesh->dim;

  // Loop Indices
  INT quad,test,trial,idim, jdim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[3];
  // Stiffness Matrix Entry
  REAL kij;
  REAL aij;
  // Coefficient Value at Quadrature Nodes
  REAL coeff_val[dim*dim+1];

  // Vector Derivatives: Gradients (PX)
  //  Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    if(dim==2 || dim==3)
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(dim==3)
    qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];

    //fill in coeff matrix a(x)
    if(coeff!=NULL) {
      (*coeff)(coeff_val,qx,time,&dim);
    } else { //make coeff the identity matrix
      coeff_val[dim*dim] = 1.0;
      for (int j = 0; j < dim*dim; j ++){
        if (j%(dim+1) == 0){
          coeff_val[j] = 1.0;
        }
        else{
          coeff_val[j] = 0.0;
        }
      }
    }

    // Basis Functions and its derivatives if necessary
    get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,mesh,FE);

    // Loop over Test Functions (Rows)
    for (test=0; test<FE->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->dof_per_elm; trial++) {
        kij=0.0;
        //do the mat-vec multiplication
        for(idim=0;idim<dim;idim++){
          aij = 0.0;
          for (jdim = 0; jdim < dim; jdim++){
            // mat-vec multipliation
            aij +=  coeff_val[dim*idim+jdim]*(FE->dphi[test*dim+jdim]);
          }
          kij += aij*(FE->dphi[trial*dim+idim]);		//diffusion matrix
        }
        kij+= coeff_val[dim*dim]*(FE->phi[trial])*(FE->phi[test]); //add mass matrix
        ALoc[test*FE->dof_per_elm+trial] += w*kij;
      }
    }
  }

  return;
}

/***************************************************************************/
/*!
* \fn REAL energyerror(REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *), void (*D_truesol)(REAL *,REAL *,REAL,void *),void (*coeff)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
*
* \brief Computes the Energy Norm of the error of a FE approximation and a true
*        solution given by a function using quadrature for any type of element.
*
* \param u 	           Numerical Solution at DOF
* \param truesol       Function to get true solution at a given point
* \param D_truesol     Function to get derivative of true solution at a given point
* \param coeff   	  	 PDE coefficients
* \param FE            FE Space
* \param mesh          Mesh Data
* \param cq            Quadrature Nodes
* \param time          Physical time to compute solution at
*
* \return error        Energy Error
*
*/
REAL energyerror(REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *), void (*D_truesol)(REAL *,REAL *,REAL,void *),void (*coeff)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
{
  INT dim = mesh->dim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(mesh->dim,sizeof(REAL));

  // FE Stuff
  INT FEtype = FE->FEtype;
  INT elm,quad,j, i;
  REAL sum = 0.0;
  REAL d_sum = 0.0;
  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  REAL* coeff_val = (REAL *) calloc(dim*dim+1,sizeof(REAL));


  REAL* val_true = (REAL *) calloc(1,sizeof(REAL));	   //true solution
  REAL* val_sol = (REAL *) calloc(1,sizeof(REAL));     //numerical solution
  REAL* d_val_true = (REAL *) calloc(dim,sizeof(REAL));  //gradient of true solution
  REAL* d_val_sol = (REAL *) calloc(dim,sizeof(REAL));   //gradient of numerical solution

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

      //define the coeff values on the quadrature nodes
      if(coeff!=NULL) {
        (*coeff)(coeff_val,qx,time,&dim);

      } else { //make coeff the identity matrix
        coeff_val[dim*dim] = 1.0;
        for (j = 0; j < dim*dim; j ++){
          if (j%(dim+1) == 0){
            coeff_val[j] = 1.0;
          }
          else{
            coeff_val[j] = 0.0;
          }
        }
      }

      // Get True Solution at Quadrature Nodes
      (*truesol)(val_true,qx,time,&(mesh->el_flag[elm]));

      // Interpolate FE solution to quadrature point
      FE_Interpolation(val_sol,u,qx,dof_on_elm,v_on_elm,FE,mesh);

      // Get True Derivative Solution at Quadrature Nodes
      (*D_truesol)(d_val_true,qx,time,&(mesh->el_flag[elm]));

      // Interpolate FE gradient solution to quadrature point
      FE_DerivativeInterpolation(d_val_sol,u,qx,dof_on_elm,v_on_elm,FE,mesh);

      // Compute Error on Element
      for(j=0;j<1;j++) { //weight * c(x) * (u - u_h)^T * (u-u_h)
        sum+=w*coeff_val[dim*dim]*(ABS(val_sol[j] - val_true[j]))*(ABS(val_sol[j] - val_true[j]));
      }


      for(i=0; i<dim; i++){
        for(j=0; j<dim; j++){ //weight * (grad u - grad u_h)^T * A(x) * (grad u - grad u_h)
          d_sum+=w*ABS(coeff_val[i*dim +j]*(d_val_sol[i] - d_val_true[i])*(d_val_sol[j] - d_val_true[j]));
        }
      }
    }
  }


  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(qx) free(qx);
  if(val_true) free(val_true);
  if(val_sol) free(val_sol);
  if(d_val_true) free(d_val_true);
  if(d_val_sol) free(d_val_sol);

  return sqrt(sum + d_sum);
}
