/*! \file Maxwell_data.h
*
*  Created by Adler, Hu, Zikatanov on 09/29/2017.
*  Edited by Casey Cavanaugh 5/20/19
*  Copyright 2016_HAZMATH__. All rights reserved.
*
* \brief This contains all the Data parameters and coefficients
*        for the Maxwell example.  This includes exact solutions,
*        RHS functions, coefficients, and boundary conditions.
*
*/

// PDE Coefficients
void permitivity(REAL *val,REAL* x,REAL time,void *param) {
  *val = 1.0;
}
void permeability(REAL *val,REAL* x,REAL time,void *param) {
  *val = 1.0;
}
void oneovermu(REAL *val,REAL* x,REAL time,void *param) {
  REAL mu = -6.66;
  permeability(&mu,x,time,param);
  *val = 1.0/mu;
}


// True Solution (if you have one)
void truesol(REAL *val,REAL* x,REAL time,void *param) {

  REAL a = 1.0;
  REAL myexp = exp(-a*time);

  //ordering (E1,E2,E3,B1,B2,B3,p)

  //Real test problem
  val[0] = -cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]) *myexp/M_PI;
  val[1] = sin(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]) *myexp/M_PI;
  val[2] = 0.0;
  val[3] = -sin(M_PI*x[0])*cos(M_PI*x[1])*cos(M_PI*x[2]) *myexp;
  val[4] = -cos(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*x[2]) *myexp;
  val[5] = 2*cos(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]) *myexp;
  val[6] = 0.0;

  /* val[0] = myexp*(x[0]*x[0] - x[0]*x[1] + 2*x[1]*x[2]);
  val[1] = myexp*(2*x[2]*x[2] - 2*x[0]*x[1]);
  val[2] = myexp*(x[1]*x[2] - 3*x[1]*x[1]);
  val[3] = myexp*(-6*x[1] - 3*x[2]);
  val[4] = myexp*2*x[1];
  val[5] = myexp*(x[0] - 2*x[1] - 2*x[2]);
  val[6] = 0.0;
  */

}
void Etrue(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[7];
  truesol(myu,x,time,param);
  val[0] = myu[0];
  val[1] = myu[1];
  val[2] = myu[2];
}
void Btrue(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[7];
  truesol(myu,x,time,param);
  val[0] = myu[3];
  val[1] = myu[4];
  val[2] = myu[5];
}
void ptrue(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[7];
  truesol(myu,x,time,param);
  *val = myu[6];
}

// Right-hand Side
//NEED TO INPUT -j as RHS
void current_density(REAL *val,REAL* x,REAL time,void *param) {
  REAL a = 1.0;
  REAL myexp = exp(-a*time);

  val[0] = myexp*(3*M_PI + 1/M_PI)*cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
  val[1] = myexp*(3*M_PI + 1/M_PI)*-sin(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]);
  val[2] = 0.0;

  /* val[0] = -myexp*(x[0]*x[0] - x[0]*x[1] + 2*x[1]*x[2] - 2);
  val[1] = -myexp*(2*x[2]*x[2] - 2*x[0]*x[1] - 4);
  val[2] = -myexp*(x[1]*x[2] - 3*x[1]*x[1] + 6);
  */
}


// Boundary Conditions
void bc(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[7];
  truesol(myu,x,time,param);
  val[0] = myu[0];
  val[1] = myu[1];
  val[2] = myu[2];
  val[3] = myu[3];
  val[4] = myu[4];
  val[5] = myu[5];
  val[6] = myu[6];
}
void bc_test(REAL *val,REAL* x,REAL time,void *param) {
  *val = 0.0;
}



/***********************************************************************************************/
/*!
   * \fn dCSRmat dcsr_invert_diagonal_matrix(dvector *diag)
   *
   * \brief create inverse of diagonal matrix using diag as diagonal entres
   *
   * \param diag  Pointer to the diagonal as a dvector
   *
   * \return D    Pointer to dCSRmat CSR inverse diagonal matrix
   *
   */
dCSRmat dcsr_invert_diagonal_matrix(dvector *diag)
{
  //local variable
  INT n = diag->row;
  INT i;

  // form the diaongal matrix
  dCSRmat D = dcsr_create(n,n,n);
  for (i=0;i<n;i++)
  {
    D.IA[i] = i;
    D.JA[i] = i;
    D.val[i]   = 1.0/diag->val[i];
  }
  D.IA[n] = n;

  // return
  return D;
}

/***********************************************************************************************/
/*!
   * \fn void build_precond_maxwell(REAL dt, dCSRmat *G, dCSRmat *Gt, dCSRmat *K , dCSRmat *Kt,
					dCSRmat *Mf, dCSRmat *Me, dCSRmat *Mv, dCSRmat *A_diag, dCSRmat *Gtb,
							dCSRmat *Ktb, dCSRmat *Kb, dCSRmat *Gb)
   *
   * \brief build the matrices for the Maxwell block preconditioner to send into the Maxwell
   *					Krylov solver
   *
   * \param dt  timestep size
   * \param G   Pointer to the discrete gradient (gradD for MFD)
   * \param G   Pointer to the discrete divergence (divD for MFD)
   * \param K   Pointer to the discrete curl (curlD for MFD)
   * \param G   Pointer to the discrete curl* (curlV for MFD)
   * \param Mf  Pointer to the face mass matrix (idenity for MFD)
   * \param Me  Pointer to the edge mass matrix (identity for MFD)
   * \param Mv  Pointer to the vertex mass matrix (identity for MFD)

   * \param A_diag    pointer to diag in LDU decomp (schur complements)
   * \param Gtb    	pointer to discrete div with BC eliminated (divD for MFD)
   * \param Ktb    	pointer to discrete curl with BCs applied (curlV for MFD)
   * \param Gb   		pointer to discrete gradient with BCs applied (divD for MFD)
   * \param Kb		pointer to discrete curl with BCs applied (curlD for MFD)
   *
   */

void build_precond_maxwell(mesh_struct* mesh, block_fespace *FE, REAL dt, dCSRmat *G, dCSRmat *Gt, dCSRmat *K , dCSRmat *Kt,
					dCSRmat *Mf, dCSRmat *Me, dCSRmat *Mv, dCSRmat *A_diag,
					dCSRmat *Gtb, dCSRmat *Ktb, dCSRmat *Kb, dCSRmat *Gb)
{




//-------------------------------
  // prepare block preconditioner

  //-------------------------------------
  // lower block triangular for LU solve
  //-------------------------------------
  block_dCSRmat Lb;
  Lb.brow = 3; Lb.bcol = 3;
  Lb.blocks = (dCSRmat **) calloc(9,sizeof(dCSRmat *));


  dCSRmat IB =  dcsr_create_identity_matrix(mesh->nface,0);
  dCSRmat IE =  dcsr_create_identity_matrix(mesh->nedge,0);
  dCSRmat Ip =  dcsr_create_identity_matrix(mesh->nv,0);

  Lb.blocks[0] = &IE;
  Lb.blocks[1] = Kt;
  Lb.blocks[2] = G;
  Lb.blocks[3] = K;
  Lb.blocks[4] = &IB;
  Lb.blocks[5] = NULL;
  Lb.blocks[6] = Gt;
  Lb.blocks[7] = NULL;
  Lb.blocks[8] = &Ip;

 block_eliminate_DirichletBC(NULL,FE,mesh,NULL, &Lb, 1,0.0);

	dcsr_cp(Lb.blocks[1], Ktb);
	dcsr_cp(Lb.blocks[3], Kb);
	dcsr_cp(Lb.blocks[2], Gb);
	dcsr_cp(Lb.blocks[6], Gtb);



  // scale offdiagonal blocks
  //we need the negatives, right?
  dcsr_axm(Kb,  dt/2.0);
  dcsr_axm(Gb,  dt/2.0);
  dcsr_axm(Ktb, -dt/2.0);
  dcsr_axm(Gtb, -dt/2.0);

  // clean up
  bdcsr_free(&Lb);
  dcsr_free(&IB);
  dcsr_free(&IE);
  dcsr_free(&Ip);
  //-------------------------------------

  // prepare diagonal blocks

  // first diagonal block: 2/dt * M_f
  dcsr_alloc(Mf->row, Mf->col, Mf->nnz, &A_diag[0]);
  dcsr_cp(Mf, &A_diag[0]);
  dcsr_axm(&A_diag[0], 2.0/dt);


  eliminate_DirichletBC(NULL,FE->var_spaces[1],mesh,NULL,&A_diag[0],0.0);


  // second diagonal block: dt/2 * K^tMfK + 2/dt*Me
  dCSRmat KTMfK;
  dcsr_rap(Kt, Mf, K, &KTMfK);
  dcsr_add(&KTMfK, dt/2.0, Me, 2.0/dt, &A_diag[1]);
  dcsr_free(&KTMfK);
  dcsr_alloc(A_diag[1].row, A_diag[1].col, A_diag[1].nnz, &KTMfK);
  dcsr_cp(&A_diag[1],&KTMfK);
  dcsr_free(&A_diag[1]);
  dcsr_free(&KTMfK);


  eliminate_DirichletBC(NULL,FE->var_spaces[0],mesh,NULL,&A_diag[1],0.0);


  // third diagonal block: dt/2 G^tMeG + 2/dt Mv
  dCSRmat GTMeG;
  dcsr_rap(Gt, Me, G, &GTMeG);
  dcsr_add(&GTMeG, dt/2.0, Mv, 2.0/dt, &A_diag[2]);
  dcsr_free(&GTMeG);

  eliminate_DirichletBC(NULL,FE->var_spaces[2],mesh,NULL,&A_diag[2],0.0);

}
