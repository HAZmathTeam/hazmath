/*! \file examples/stokes/stokes_precond.h
*
*  Created by Xiaozhe Hu on 08/24/2021.
*  Copyright 2021_HAZMATH__. All rights reserved.
*
* \brief This contains all preconditioners for Stokes example.
*
*/

/*!
* \fn precond_block_data *get_precond_block_data_stokes(block_dCSRmat *Ab, dCSRmat *Mp, linear_itsolver_param *itparam,
AMG_param *amgparam)
*
* \brief get data for block preconditioner for solving the Stokes problem
*
* \param Ab                Point to a block_dCSRmat matrix
* \param Mp                Point to the pressue mass matrix in the dCSR format
* \param itparam           Parameters of iterative methods
* \param amgparam          Parameters of AMG methods
*
* \author Xiaozhe Hu 08/24/2021
*
* \note this is a special function only for Stokes -- Xiaozhe
*
*/
static precond_block_data *get_precond_block_data_stokes(block_dCSRmat *Ab,dCSRmat *Mp,linear_itsolver_param *itparam,AMG_param *amgparam)
{
  // return variable
  precond_block_data *precdata=(precond_block_data *)malloc(1*sizeof(precond_block_data));

  // local variables
  const SHORT prtlvl = itparam->linear_print_level;
  INT brow=Ab->brow, bcol=Ab->bcol;
  INT n1, n2;
  INT dim = brow-1; // 2D: 3x3 block matrix | 3D: 4x4 block matrix
  INT status = SUCCESS;

  // initialize the block preconditioner
  precond_block_data_null(precdata);

  // store the big block
  precdata->Abcsr=Ab;

  // data for diagonal blocks
  precdata->A_diag = (dCSRmat *)calloc(2, sizeof(dCSRmat)); // two blocks, one for velocity and one for pressure

  //-----------------------------------------
  // solver data for the velocity part
  //-----------------------------------------
  // grab the velocity block
  n1=0; n2=dim;
  bdcsr_getdiagblk_dcsr(Ab,n1,n2, &(precdata->A_diag[0]));

  /*
  // get the LU factorization of this block
  precdata->LU_diag=malloc(sizeof(void *));
  precdata->LU_diag[0]=factorize_UMF(&(precdata->A_diag[0]),0);
  */

  // get AMG data for this block
  AMG_data **mgl = (AMG_data **)calloc(1, sizeof(AMG_data *));

  SHORT max_levels;
  if (amgparam) max_levels = amgparam->max_levels;

  mgl[0] = amg_data_create(max_levels);
  dcsr_alloc(precdata->A_diag[0].row, precdata->A_diag[0].row, precdata->A_diag[0].nnz, &mgl[0][0].A);
  dcsr_cp(&(precdata->A_diag[0]), &mgl[0][0].A);
  mgl[0][0].b=dvec_create(precdata->A_diag[0].row);
  mgl[0][0].x=dvec_create(precdata->A_diag[0].row);

  switch (amgparam->AMG_type) {

    case UA_AMG: // Unsmoothed Aggregation AMG
    if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
    status = amg_setup_ua(mgl[0], amgparam);
    break;

    case SA_AMG: // Smoothed Aggregation AMG
    if ( prtlvl > PRINT_NONE ) printf("\n Calling SA AMG ...\n");
    status = amg_setup_sa(mgl[0], amgparam);
    break;

    default: // UA AMG
    if ( prtlvl > PRINT_NONE ) printf("\n Calling UA AMG ...\n");
    status = amg_setup_ua(mgl[0], amgparam);
    break;
  }

  // assign mgl
  if (amgparam) precdata->amgparam = amgparam;
  precdata->mgl = mgl;

  //-----------------------------------------
  // solver data for the pressure part
  //-----------------------------------------
  // grab the pressure block
  dcsr_alloc(Mp->row, Mp->row, Mp->nnz, &(precdata->A_diag[1]));
  dcsr_cp(Mp, &(precdata->A_diag[1]));

  // get the diagonal of the pressure block
  precdata->diag = (dvector **)calloc(1, sizeof(dvector *));
  precdata->diag[0] = (dvector *)calloc(1, sizeof(dvector));
  dcsr_getdiag(0, Mp, precdata->diag[0]);

  //-----------------------------------------
  // allocate work space
  //-----------------------------------------
  INT N = precdata->A_diag[0].row+precdata->A_diag[1].row; // total degrees of freedoms
  precdata->r = dvec_create(2*N);

  // return
  return precdata;
}

/*!
* \fn void precond_block_data_free_stokes(precond_block_data *precdata)
*
* \brief Free precond_block_data structure for Stokes problem (set values to 0 and pointers to NULL)
*
* \param precdata      Pointer to the precond_block_data structure (OUTPUT)
*
* \author Xiaozhe Hu 08/25/2021
*
*/
static void precond_block_data_free_stokes(precond_block_data *precdata)
{
  // free velocity block
  if (&(precdata->A_diag[0])) dcsr_free(&(precdata->A_diag[0]));
  // free pressure block
  if (&(precdata->A_diag[1])) dcsr_free(&(precdata->A_diag[1]));
  if (precdata->A_diag) free(precdata->A_diag);

  /*
  if (precdata->LU_diag[0]) free(precdata->LU_diag[0]);
  if (precdata->LU_diag) free(precdata->LU_diag);
  */

  // free AMG data for velocity block
  amg_data_free(precdata->mgl[0], &precdata->amgparam[0]);
  free(precdata->mgl[0]);
  if(precdata->mgl) free(precdata->mgl);

  // free diagonal of the pressure block
  if (precdata->diag[0]) dvec_free(precdata->diag[0]);
  if (precdata->diag[0]) free(precdata->diag[0]);
  if (precdata->diag) free(precdata->diag);

  // free work space
  if (&(precdata->r)) dvec_free(&(precdata->r));

  // free the whole precond data
  if (precdata) free(precdata);

  return;
}

/**
* \fn void precond_block_diag_stokes(REAL *r, REAL *z, void *data)
* \brief block diagonal preconditioning for Stokes problem
*        AMG for velocity block
*        Jacobi for pressure block
*
* \param r     Pointer to the vector needs preconditioning
* \param z     Pointer to preconditioned vector
* \param data  Pointer to precondition data
*
* \author Xiaozhe Hu
* \date   08/24/2021
*/
static void precond_block_diag_stokes(REAL *r, REAL *z,void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  //-----------------------------------------
  // solver data for the whole block matrix
  //-----------------------------------------
  block_dCSRmat *A = precdata->Abcsr;

  //-----------------------------------------
  // solver data for the velocity part
  //-----------------------------------------
  // data for the vecolity block
  dCSRmat *Au = &(precdata->A_diag[0]);
  const INT Nu = Au->row;

  //void **LU_diag = precdata->LU_diag;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;

  //-----------------------------------------
  // solver data for the pressure part
  //-----------------------------------------
  dCSRmat *Mp = &(precdata->A_diag[1]);
  const INT Np = Mp->row;

  //-----------------------------------------
  // local variabl
  //-----------------------------------------
  INT i;
  // total size
  const INT N = Nu + Np;

  //-----------------------------------------
  // back up r, setup z;
  //-----------------------------------------
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);

  //-----------------------------------------
  // prepare
  //-----------------------------------------
  dvector ru, rp, zu, zp;

  ru.row = Nu; zu.row = Nu;
  rp.row = Np; zp.row = Np;

  ru.val = r, rp.val = &(r[Nu]);
  zu.val = z, zp.val = &(z[Nu]);

  //-----------------------------------------
  // main part of the preconditioner
  //-----------------------------------------
  // Preconditioning the velocity block
  //-----------------------------------------
  // use direct solver
  //solve_UMF(Au, &ru, &zu, LU_diag[0], 0);

  // use AMG solver
  mgl[0]->b.row=Nu; array_cp(Nu, ru.val, mgl[0]->b.val); // residual is an input
  mgl[0]->x.row=Nu; dvec_set(Nu, &mgl[0]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
  array_cp(Nu, mgl[0]->x.val, zu.val);
  //-----------------------------------------

  //-----------------------------------------
  // Preconditioning the pressure block
  //-----------------------------------------
  // Simple Jacobi iteration
  smoother_dcsr_jacobi(&zp, 0, Np, 1, Mp, &rp, 1);
  //-----------------------------------------

  // restore r
  array_cp(N, tempr->val, r);

}

/**
* \fn void precond_block_diag_stokes_krylov(REAL *r, REAL *z, void *data)
* \brief block diagonal preconditioning for Stokes problem
*        AMG+krylov for velocity block
*        Jacobi+Krylov for pressure block
*
* \param r     Pointer to the vector needs preconditioning
* \param z     Pointer to preconditioned vector
* \param data  Pointer to precondition data
*
* \author Xiaozhe Hu
* \date   08/24/2021
*/
static void precond_block_diag_stokes_krylov(REAL *r,REAL *z,void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  //-----------------------------------------
  // solver data for the whole block matrix
  //-----------------------------------------
  block_dCSRmat *A = precdata->Abcsr;

  //-----------------------------------------
  // solver data for the velocity part
  //-----------------------------------------
  // data for the vecolity block
  dCSRmat *Au = &(precdata->A_diag[0]);
  const INT Nu = Au->row;

  //void **LU_diag = precdata->LU_diag;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;

  //-----------------------------------------
  // solver data for the pressure part
  //-----------------------------------------
  dCSRmat *Mp = &(precdata->A_diag[1]);
  const INT Np = Mp->row;

  //-----------------------------------------
  // local variabl
  //-----------------------------------------
  INT i;
  // total size
  const INT N = Nu + Np;

  //-----------------------------------------
  // back up r, setup z;
  //-----------------------------------------
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);

  //-----------------------------------------
  // prepare
  //-----------------------------------------
  dvector ru, rp, zu, zp;

  ru.row = Nu; zu.row = Nu;
  rp.row = Np; zp.row = Np;

  ru.val = r, rp.val = &(r[Nu]);
  zu.val = z, zp.val = &(z[Nu]);

  //-----------------------------------------
  // main part of the preconditioner
  //-----------------------------------------
  // Preconditioning the velocity block
  //-----------------------------------------
  // AMG preconditioned Krylov method
  precond_data pcdata_u;
  param_amg_to_prec(&pcdata_u,amgparam);
  precond pc_u;
  pc_u.fct = precond_amg;

  pcdata_u.max_levels = mgl[0][0].num_levels;
  pcdata_u.mgl_data = mgl[0];

  pc_u.data = &pcdata_u;

  dcsr_pvfgmres(Au, &ru, &zu, &pc_u, 1e-3, 100, 100, 1, 1);
  //-----------------------------------------

  //-----------------------------------------
  // Preconditioning the pressure block
  //-----------------------------------------
  // Jacobi preconditioned Krylov method
  precond pc_p;
  pc_p.data = precdata->diag[0];
  pc_p.fct  = precond_diag;

  dcsr_pvfgmres(Mp, &rp, &zp, &pc_p, 1e-3, 100, 100, 1, 1);
  //-----------------------------------------

  // restore r
  array_cp(N, tempr->val, r);

}

/**
* \fn void precond_block_lower_stokes(REAL *r, REAL *z, void *data)
* \brief block lower trianglar preconditioning for Stokes problem
*        AMG for velocity block
*        Jacobi for pressure block
*
* \param r     Pointer to the vector needs preconditioning
* \param z     Pointer to preconditioned vector
* \param data  Pointer to precondition data
*
* \author Xiaozhe Hu
* \date   08/24/2021
*/
static void precond_block_lower_stokes(REAL *r,REAL *z,void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  //-----------------------------------------
  // solver data for the whole block matrix
  //-----------------------------------------
  block_dCSRmat *A = precdata->Abcsr;

  //-----------------------------------------
  // solver data for the velocity part
  //-----------------------------------------
  // data for the vecolity block
  dCSRmat *Au = &(precdata->A_diag[0]);
  const INT Nu = Au->row;

  //void **LU_diag = precdata->LU_diag;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;

  //-----------------------------------------
  // solver data for the pressure part
  //-----------------------------------------
  dCSRmat *Mp = &(precdata->A_diag[1]);
  const INT Np = Mp->row;

  //-----------------------------------------
  // local variabl
  //-----------------------------------------
  INT i;
  // total size
  const INT N = Nu + Np;

  //-----------------------------------------
  // back up r, setup z;
  //-----------------------------------------
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);

  //-----------------------------------------
  // prepare
  //-----------------------------------------
  dvector ru, rp, zu, zp;

  ru.row = Nu; zu.row = Nu;
  rp.row = Np; zp.row = Np;

  ru.val = r, rp.val = &(r[Nu]);
  zu.val = z, zp.val = &(z[Nu]);

  //-----------------------------------------
  // main part of the preconditioner
  //-----------------------------------------
  // Preconditioning the velocity block
  //-----------------------------------------
  // use direct solver
  //solve_UMF(Au, &ru, &zu, LU_diag[0], 0);

  // use AMG solver
  mgl[0]->b.row=Nu; array_cp(Nu, ru.val, mgl[0]->b.val); // residual is an input
  mgl[0]->x.row=Nu; dvec_set(Nu, &mgl[0]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
  array_cp(Nu, mgl[0]->x.val, zu.val);
  //-----------------------------------------

  //-----------------------------------------
  // update the residual
  //-----------------------------------------
  //array_cp(N, tempr->val, r);
  bdcsr_aAxpy(-1.0, A, z, r);

  //-----------------------------------------
  // Preconditioning the pressure block
  //-----------------------------------------
  // Simple Jacobi iteration
  smoother_dcsr_jacobi(&zp, 0, Np, 1, Mp, &rp, 1);
  //-----------------------------------------

  // restore r
  array_cp(N, tempr->val, r);

}

/**
* \fn void precond_block_lower_stokes_krylov(REAL *r, REAL *z, void *data)
* \brief block lower trianglar preconditioning for Stokes problem
*        AMG+Krylov for velocity block
*        Jacobi+Krylov for pressure block
*
* \param r     Pointer to the vector needs preconditioning
* \param z     Pointer to preconditioned vector
* \param data  Pointer to precondition data
*
* \author Xiaozhe Hu
* \date   08/24/2021
*/
static void precond_block_lower_stokes_krylov(REAL *r,REAL *z,void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  //-----------------------------------------
  // solver data for the whole block matrix
  //-----------------------------------------
  block_dCSRmat *A = precdata->Abcsr;

  //-----------------------------------------
  // solver data for the velocity part
  //-----------------------------------------
  // data for the vecolity block
  dCSRmat *Au = &(precdata->A_diag[0]);
  const INT Nu = Au->row;

  //void **LU_diag = precdata->LU_diag;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;

  //-----------------------------------------
  // solver data for the pressure part
  //-----------------------------------------
  dCSRmat *Mp = &(precdata->A_diag[1]);
  const INT Np = Mp->row;

  //-----------------------------------------
  // local variabl
  //-----------------------------------------
  INT i;
  // total size
  const INT N = Nu + Np;

  //-----------------------------------------
  // back up r, setup z;
  //-----------------------------------------
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);

  //-----------------------------------------
  // prepare
  //-----------------------------------------
  dvector ru, rp, zu, zp;

  ru.row = Nu; zu.row = Nu;
  rp.row = Np; zp.row = Np;

  ru.val = r, rp.val = &(r[Nu]);
  zu.val = z, zp.val = &(z[Nu]);

  //-----------------------------------------
  // main part of the preconditioner
  //-----------------------------------------
  // Preconditioning the velocity block
  //-----------------------------------------
  // AMG preconditioned Krylov method
  precond_data pcdata_u;
  param_amg_to_prec(&pcdata_u,amgparam);
  precond pc_u;
  pc_u.fct = precond_amg;

  pcdata_u.max_levels = mgl[0][0].num_levels;
  pcdata_u.mgl_data = mgl[0];

  pc_u.data = &pcdata_u;

  dcsr_pvfgmres(Au, &ru, &zu, &pc_u, 1e-3, 100, 100, 1, 1);
  //-----------------------------------------

  //-----------------------------------------
  // update the residual
  //-----------------------------------------
  //array_cp(N, tempr->val, r);
  bdcsr_aAxpy(-1.0, A, z, r);

  //-----------------------------------------
  // Preconditioning the pressure block
  //-----------------------------------------
  // Jacobi preconditioned Krylov method
  precond pc_p;
  pc_p.data = precdata->diag[0];
  pc_p.fct  = precond_diag;

  dcsr_pvfgmres(Mp, &rp, &zp, &pc_p, 1e-3, 100, 100, 1, 1);
  //-----------------------------------------

  // restore r
  array_cp(N, tempr->val, r);

}

/**
* \fn void precond_block_upper_stokes(REAL *r, REAL *z, void *data)
* \brief block upper triangle preconditioning for Stokes problem
*        AMG for velocity block
*        Jacobi for pressure block
*
* \param r     Pointer to the vector needs preconditioning
* \param z     Pointer to preconditioned vector
* \param data  Pointer to precondition data
*
* \author Xiaozhe Hu
* \date   08/24/2021
*/
static void precond_block_upper_stokes(REAL *r,REAL *z,void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  //-----------------------------------------
  // solver data for the whole block matrix
  //-----------------------------------------
  block_dCSRmat *A = precdata->Abcsr;

  //-----------------------------------------
  // solver data for the velocity part
  //-----------------------------------------
  // data for the vecolity block
  dCSRmat *Au = &(precdata->A_diag[0]);
  const INT Nu = Au->row;

  //void **LU_diag = precdata->LU_diag;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;

  //-----------------------------------------
  // solver data for the pressure part
  //-----------------------------------------
  dCSRmat *Mp = &(precdata->A_diag[1]);
  const INT Np = Mp->row;

  //-----------------------------------------
  // local variabl
  //-----------------------------------------
  INT i;
  // total size
  const INT N = Nu + Np;

  //-----------------------------------------
  // back up r, setup z;
  //-----------------------------------------
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);

  //-----------------------------------------
  // prepare
  //-----------------------------------------
  dvector ru, rp, zu, zp;

  ru.row = Nu; zu.row = Nu;
  rp.row = Np; zp.row = Np;

  ru.val = r, rp.val = &(r[Nu]);
  zu.val = z, zp.val = &(z[Nu]);

  //-----------------------------------------
  // main part of the preconditioner
  //-----------------------------------------
  // Preconditioning the pressure block
  //-----------------------------------------
  // Simple Jacobi iteration
  smoother_dcsr_jacobi(&zp, 0, Np, 1, Mp, &rp, 1);
  //-----------------------------------------

  //-----------------------------------------
  // update the residual
  //-----------------------------------------
  //array_cp(N, tempr->val, r);
  bdcsr_aAxpy(-1.0, A, z, r);
  //-----------------------------------------

  //-----------------------------------------
  // Preconditioning the velocity block
  //-----------------------------------------
  // use direct solver
  //solve_UMF(Au, &ru, &zu, LU_diag[0], 0);

  // use AMG solver
  mgl[0]->b.row=Nu; array_cp(Nu, ru.val, mgl[0]->b.val); // residual is an input
  mgl[0]->x.row=Nu; dvec_set(Nu, &mgl[0]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
  array_cp(Nu, mgl[0]->x.val, zu.val);
  //-----------------------------------------

  // restore r
  array_cp(N, tempr->val, r);

}

/**
* \fn void precond_block_upper_stokes_krylov(REAL *r, REAL *z, void *data)
* \brief block upper triangle preconditioning for Stokes problem
*        AMG+krylov for velocity block
*        Jacobi+krylov for pressure block
*
* \param r     Pointer to the vector needs preconditioning
* \param z     Pointer to preconditioned vector
* \param data  Pointer to precondition data
*
* \author Xiaozhe Hu
* \date   08/24/2021
*/
static void precond_block_upper_stokes_krylov(REAL *r,REAL *z,void *data)
{
  precond_block_data *precdata=(precond_block_data *)data;
  dvector *tempr = &(precdata->r);

  //-----------------------------------------
  // solver data for the whole block matrix
  //-----------------------------------------
  block_dCSRmat *A = precdata->Abcsr;

  //-----------------------------------------
  // solver data for the velocity part
  //-----------------------------------------
  // data for the vecolity block
  dCSRmat *Au = &(precdata->A_diag[0]);
  const INT Nu = Au->row;

  //void **LU_diag = precdata->LU_diag;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;

  //-----------------------------------------
  // solver data for the pressure part
  //-----------------------------------------
  dCSRmat *Mp = &(precdata->A_diag[1]);
  const INT Np = Mp->row;

  //-----------------------------------------
  // local variabl
  //-----------------------------------------
  INT i;
  // total size
  const INT N = Nu + Np;

  //-----------------------------------------
  // back up r, setup z;
  //-----------------------------------------
  array_cp(N, r, tempr->val);
  array_set(N, z, 0.0);

  //-----------------------------------------
  // prepare
  //-----------------------------------------
  dvector ru, rp, zu, zp;

  ru.row = Nu; zu.row = Nu;
  rp.row = Np; zp.row = Np;

  ru.val = r, rp.val = &(r[Nu]);
  zu.val = z, zp.val = &(z[Nu]);

  //-----------------------------------------
  // main part of the preconditioner
  //-----------------------------------------
  // Preconditioning the pressure block
  //-----------------------------------------
  // Jacobi preconditioned Krylov method
  precond pc_p;
  pc_p.data = precdata->diag[0];
  pc_p.fct  = precond_diag;

  dcsr_pvfgmres(Mp, &rp, &zp, &pc_p, 1e-3, 100, 100, 1, 1);
  //-----------------------------------------

  //-----------------------------------------
  // update the residual
  //-----------------------------------------
  //array_cp(N, tempr->val, r);
  bdcsr_aAxpy(-1.0, A, z, r);
  //-----------------------------------------

  //-----------------------------------------
  // Preconditioning the velocity block
  //-----------------------------------------
  // AMG preconditioned Krylov method
  precond_data pcdata_u;
  param_amg_to_prec(&pcdata_u,amgparam);
  precond pc_u;
  pc_u.fct = precond_amg;

  pcdata_u.max_levels = mgl[0][0].num_levels;
  pcdata_u.mgl_data = mgl[0];

  pc_u.data = &pcdata_u;

  dcsr_pvfgmres(Au, &ru, &zu, &pc_u, 1e-3, 100, 100, 1, 1);
  //-----------------------------------------

  // restore r
  array_cp(N, tempr->val, r);

}

/*!
* \fn static void set_precond_type_stokes(precond *prec,linear_itsolver_param *itparam)
*
* \brief set the type of block preconditioner for solving the Stokes problem
*
* \param prec              preconditioner struct
* \param itparam           Parameters of iterative methods
*
* \author Xiaozhe Hu 08/24/2021
*
* \note this is a special function only for Stokes -- Xiaozhe
*
*/
static void set_precond_type_stokes(precond *prec,linear_itsolver_param *itparam)
{

  // Input file will set the appropriate preconditioner type
  switch (itparam->linear_precond_type) {
    // Block Diagonal with AMG for velocity and Jacobi for pressure
    case 40:
    prec->fct = precond_block_diag_stokes;
    printf(" --> Using block diagonal preconditioner for Stokes (AMG for velocity, Jacobi for pressure)\n");
    break;

    // Block Lower Triangular with AMG for velocity and Jacobi for pressure
    case 41:
    prec->fct = precond_block_lower_stokes;
    printf(" --> Using block lower triangular preconditioner for Stokes (AMG for velocity, Jacobi for pressure)\n");
    break;

    // Block Upper Traingular with AMG for velocity and Jacobi for pressure
    case 42:
    prec->fct = precond_block_upper_stokes;
    printf(" --> Using block upper triangular preconditioner for Stokes (AMG for velocity, Jacobi for pressure)\n");
    break;

    // Block Diagonal with AMG+krylov for velocity and Jacobi+krylov for pressure
    case 50:
    prec->fct = precond_block_diag_stokes_krylov;
    printf(" --> Using block diagonal preconditioner for Stokes (AMG+krylov for velocity, Jacobi+krylov for pressure)\n");
    break;

    // Block Lower Triangular with AMG+krylov for velocity and Jacobi+krylov for pressure
    case 51:
    prec->fct = precond_block_lower_stokes_krylov;
    printf(" --> Using block lower triangular preconditioner for Stokes (AMG+krylov for velocity, Jacobi+krylov for pressure)\n");
    break;

    // Block Upper Triangular with AMG+krylov for velocity and Jacobi+krylov for pressure
    case 52:
    prec->fct = precond_block_upper_stokes_krylov;
    printf(" --> Using block upper triangular preconditioner for Stokes (AMG+krylov for velocity, Jacobi+krylov for pressure)\n");
    break;

    // Default is Block Lower Triangular with AMG+krylov for velocity and Jacobi+krylov for pressure
    default:
    prec->fct = precond_block_lower_stokes_krylov;
    printf(" --> Using block lower triangular preconditioner for Stokes (AMG+krylov for velocity, Jacobi+krylov for pressure)\n");
    break;
  }

  return;
}
