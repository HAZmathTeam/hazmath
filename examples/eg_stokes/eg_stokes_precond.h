//-----------------------------------------------------------
// subroutines for preconditioner
//-----------------------------------------------------------
/*********************************************************************************/
/*!
 * \fn dvector *get_diag_bdcsr(block_dCSRmat *Ab, const INT n1, const INT n2)
 *
 * \brief   extracts the diagonal of Ab(n1:n2,n1:n2) and stored them in a dvector
 *
 * \param Ab    Point to a block_dCSRmat matrix
 * \param n1    staring index for the blocks
 * \param n2    ending index for the blocks
 *
 * \note index starts with 0
 *
 */
static dvector *get_diag_bdcsr(block_dCSRmat *Ab,
                               const INT n1,
                               const INT n2)
{

  // local variables
  INT i;
  INT total_size = 0;
  INT nb = Ab->brow;

  // loop 1: get size
  for (i=n1; i<n2; i++)
  {
    total_size = total_size + Ab->blocks[i*(nb+1)]->row;
  }

  // loop 2: get diagonals
  dvector *diag_A = dvec_create_p(total_size); // allocate
  dvector temp_vec;
  INT temp_n;

  // reset total_size
  total_size = 0;

  for (i=n1; i<n2; i++)
  {
      printf("i=%d\n",i);
     // get size for current block
     temp_n = Ab->blocks[i*(nb+1)]->row;

     // get the diagonal entries of the current block
     dcsr_getdiag(temp_n, Ab->blocks[i*(nb+1)], &temp_vec);

     // copy diagonal entry to the correct place
     array_cp(temp_n, temp_vec.val, diag_A->val+total_size);

     // update total size
     total_size = total_size + temp_n;

     // free temp_vec
     dvec_free(&temp_vec);

  }

  return diag_A;

}

/*********************************************************************************/
/*!
 * \fn dCSRmat *get_diag_blocks(block_dCSRmat *Ab, const INT n10, const INT n20)
 *
 * \brief   get the diagonal blocks Ab(n1:n2,n1:n2) and store them in a dCSR matrix;
 *
 * \param Ab    Point to a block_dCSRmat matrix
 * \param n10    staring index for the blocks
 * \param n20    ending index for the blocks
 *
 * \note    Memory space for the dCSRmat matrix is allocated inside this function! -- Xiaozhe Hu
 * \note    modeled on bdcsr_2_dcsr from utilities/format.c -- Ludmil
 *
 */
static INT get_diag_blocks(block_dCSRmat *Ab,
                           const INT n10,
                           const INT n20,
                           dCSRmat *A)
{
  // local variables
  INT m=0,n=0,nnz=0;
  const INT mb=Ab->brow, nb=Ab->bcol, n_blocks=mb*nb;
  dCSRmat **blockptr=Ab->blocks, *blockptrij;
  INT i,j,ij,ir,i1,length,ilength,start,irmrow,irmrowp1;
  INT *row, *col;
  INT n1=n10,n2=n20;
  if(n10<0) n1 = 0;
  if(n20>mb) n2=mb;
  if(n2<n1) {j=n2;n2=n1;n1=j;}
  INT nblk=n2-n1+1; // number of blocks
  // flag for errors
  SHORT status = SUCCESS;
  row = (INT *)calloc(mb+1,sizeof(INT));
  col = (INT *)calloc(nb+1,sizeof(INT));
  // get the size of A
  row[0]=0; col[0]=0;

  // count number of rows
  for (i=n1;i<n2;++i) {
    status = ERROR_BLKMAT_ZERO;
    for (j=n1; j<n2; ++j){
      if (blockptr[i*nb+j]) {
	m+=blockptr[i*nb+j]->row;
	row[i+1]=m;
	status = SUCCESS;
	break;
      }
    }
    // check error
    if (status < SUCCESS) check_error(ERROR_BLKMAT_ZERO, __FUNCTION__);
  }

  // count number of columns
  for (i=n1;i<n2;++i) {
    status = ERROR_BLKMAT_ZERO;
    for (j=n1;j<n2;++j){
      if (blockptr[j*mb+i]) {
	n+=blockptr[j*mb+i]->col;
	col[i+1]=n;
	status = SUCCESS;
	break;
      }
    }
    // check error
    if (status < SUCCESS) check_error(ERROR_BLKMAT_ZERO, __FUNCTION__);
  }
  // count number of nonzeros
  for (i=n1;i<n2;++i) {
    for (j=n1;j<n2;++j){
      if (blockptr[i*mb+j]) {
	nnz+=blockptr[i*mb+j]->nnz;
      }
    }
  }
  // memory space allocation
  //A = dcsr_create_p(m,n,nnz);
  dcsr_alloc(m,n,nnz,A);
  // set dCSRmat for A
  A->IA[0]=0;
  for (i=n1;i<n2;++i) {
    for (ir=row[i];ir<row[i+1];ir++) {
      for (length=j=n1;j<n2;++j) {
	ij=i*nb+j;
	blockptrij=blockptr[ij];
	if (blockptrij && blockptrij->nnz>0) {
	  start=A->IA[ir]+length;
	  irmrow=ir-row[i];irmrowp1=irmrow+1;
	  ilength=blockptrij->IA[irmrowp1]-blockptrij->IA[irmrow];
	  if (ilength>0) {
	    memcpy((A->val+start),(blockptrij->val+blockptrij->IA[irmrow]),ilength*sizeof(REAL));
	    memcpy((A->JA+start),(blockptrij->JA+blockptrij->IA[irmrow]), ilength*sizeof(INT));
	    // shift column index
	    for (i1=0;i1<ilength;i1++) A->JA[start+i1]+=col[j];
	    length+=ilength;
	  }
	}
      } // end for j
      A->IA[ir+1]=A->IA[ir]+length;
    } // end for ir
  } // end for i
  A->nnz=A->IA[row[n2]];
  /* for(i=n1;i<=n2;i++){ */
  /*   fprintf(stdout,"\nblk=%d,row=%d",i,row[i]); */
  /* }   */
  /* for(i=n1;i<=n2;i++){ */
  /*   fprintf(stdout,"\nblk=%d,row=%d",i,col[i]); */
  /* } */
  /* fprintf(stdout,"\n*** IAend=%d\n",A->IA[row[n2]]); */
  /* fprintf(stdout,"\nA11 data:(%d,%d,%d):rows:(%d,%d)\n",A->row,A->col,A->nnz,row[n2-1],row[n2]); */
  free(row);
  free(col);
  return 0;
}
/*********************************************************************************/

/**************************************************************************************/
/*!
 * \fn precond_block_data *get_precond_block_data_eg_stokes(block_dCSRmat *Ab, const INT p_ndof, REAL *el_vol)
 *
 * \brief get data for block preconditioner for solving the EG stokes
 *
 * \param Ab    Point to a block_dCSRmat matrix
 *
 * \note this is a special function only for eg stokes -- Xiaozhe
 *
 */
static precond_block_data *get_precond_block_data_eg_stokes(block_dCSRmat *Ab,
                                                            const INT p_ndof,
                                                            REAL *el_vol,
                                                            linear_itsolver_param *itparam,
                                                            AMG_param *amgparam)
{
  // return variable
  precond_block_data *precdata=(precond_block_data *)malloc(1*sizeof(precond_block_data));

  // local variables
  const SHORT prtlvl = itparam->linear_print_level;
  INT brow=Ab->brow, bcol=Ab->bcol;
  INT n1, n2;
  INT dim = brow-2;
  INT status = SUCCESS;
  //nblk,iblk,n1,n2,j,k,l,m,n,iaa,iab;

  // initialize the block preconditioner
  precond_block_data_null(precdata);

  // store the big block
  precdata->Abcsr=Ab;

  precdata->A_diag = (dCSRmat *)calloc(2, sizeof(dCSRmat));

  //-----------------------------------------
  // solver data for the velocity part
  //-----------------------------------------
  // grab the velocity block without the eg part
  n1=0; n2=dim;
  get_diag_blocks(Ab,n1,n2, &(precdata->A_diag[0]));

  // get the LU factorization of this block
  precdata->LU_diag=malloc(sizeof(void *));
  precdata->LU_diag[0]=factorize_UMF(&(precdata->A_diag[0]),0);

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

  // grab the velocity block including the eg part
  n1=0; n2 = dim+1;  // this time we need to include the EG part
  get_diag_blocks(Ab,n1,n2, &(precdata->A_diag[1]));

  // grab the diagonal of velocity block including EG part
  precdata->diag = malloc(sizeof(dvector *));
  n1=0; n2 = dim+1;  // this time we need to include the EG part
  precdata->diag[0] = get_diag_bdcsr(Ab, n1, n2);

  //-----------------------------------------
  // solver data for the pressure part
  //-----------------------------------------
  precdata->el_vol = dvec_create_p(p_ndof);
  array_cp(p_ndof, el_vol, precdata->el_vol->val);

  //-----------------------------------------
  // allocate work space
  //-----------------------------------------
  INT N = precdata->diag[0]->row+precdata->el_vol->row; // total degrees of freedoms
  precdata->r = dvec_create(N);

  // return
  return precdata;
}

/**************************************************************************************/
/*!
 * \fn void precond_block_data_free_eg_stokes(precond_block_data *precdata)
 *
 * \brief Free precond_block_data structure (set values to 0 and pointers to NULL)
 *
 * \param precdata      Pointer to the precond_block_data structure (OUTPUT)
 *
 */
void precond_block_data_free_eg_stokes(precond_block_data *precdata)
{

    if (&(precdata->A_diag[0])) dcsr_free(&(precdata->A_diag[0]));
    if (&(precdata->A_diag[1])) dcsr_free(&(precdata->A_diag[1]));
    if (precdata->A_diag) free(precdata->A_diag);

    if (precdata->LU_diag[0]) free(precdata->LU_diag[0]);
    if (precdata->LU_diag) free(precdata->LU_diag);

    amg_data_free(precdata->mgl[0], &precdata->amgparam[0]);
    free(precdata->mgl[0]);
    if(precdata->mgl) free(precdata->mgl);

    // the next three free looks strange to me...  --XH
    if (precdata->diag[0]) free(precdata->diag[0]);
    if (precdata->diag) free(precdata->diag);

    if (precdata->el_vol) free(precdata->el_vol);

    if (&(precdata->r)) dvec_free(&(precdata->r));

    return;
}


/***********************************************************************************************/
/**
 * \fn void precond_block_diag_eg_stokes_additive (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning for EG Stokes problem
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/04/2021
 */
void precond_block_diag_eg_stokes_additive(REAL *r,
                                           REAL *z,
                                           void *data)
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
  // data for the vecolity block without eg part
  dCSRmat *A_diag = precdata->A_diag;
  void **LU_diag = precdata->LU_diag;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;

  // data for the velocity block including eg part
  dvector **velocity_diag = precdata->diag;

  //-----------------------------------------
  // solver data for the pressure part
  //-----------------------------------------
  dvector *el_vol = precdata->el_vol;

  //-----------------------------------------
  // local variabl
  //-----------------------------------------
  INT i;
  INT brow = A->brow;

  //-----------------------------------------
  // get different sizes
  //-----------------------------------------
  // get size of the velocity block without the eg part
  const INT Nu_wo_eg = A_diag[0].row;

  // get size of the veclosity block including the eg part
  const INT Nu = velocity_diag[0]->row;

  // get size of the pressure block
  const INT Np = el_vol->row;

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
  dvector ru_wo_eg, ru, rp, zu_wo_eg, zu, zp;

  ru_wo_eg.row = Nu_wo_eg; zu_wo_eg.row = Nu_wo_eg;
  ru.row = Nu; zu.row = Nu;
  rp.row = Np; zp.row = Np;

  ru_wo_eg.val = r; ru.val = r, rp.val = &(r[Nu]);
  zu_wo_eg.val = z; zu.val = z, zp.val = &(z[Nu]);

  //-----------------------------------------
  // main part of the preconditioner
  //-----------------------------------------
  // Preconditioning the velocity block without the eg part
  // use direct solver
  solve_UMF(&(A_diag[0]), &ru_wo_eg, &zu_wo_eg, LU_diag[0], 0);

  /*
  mgl[0]->b.row=Nu_wo_eg; array_cp(Nu_wo_eg, ru_wo_eg.val, mgl[0]->b.val); // residual is an input
  mgl[0]->x.row=Nu_wo_eg; dvec_set(Nu_wo_eg, &mgl[0]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
  array_cp(Nu_wo_eg, mgl[0]->x.val, zu_wo_eg.val);
  */

  // Preconditioning the velocity block including the eg part
  for(i=0; i<Nu; i++){
    zu.val[i] = zu.val[i] + ru.val[i]/velocity_diag[0]->val[i];
  }

  //getchar();

  // Preconditioning the pressure block
  // Diagonal matrix for P0
  for(i=0;i<Np;i++){
    zp.val[i] = rp.val[i]/el_vol->val[i];
  }

  // restore r
  array_cp(N, tempr->val, r);

}


/***********************************************************************************************/
/**
 * \fn void precond_block_diag_eg_stokes_multiplicative (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning for EG Stokes problem
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   03/04/2021
 */
void precond_block_diag_eg_stokes_multiplicative(REAL *r,
                                                 REAL *z,
                                                 void *data)
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
  // data for the vecolity block without eg part
  dCSRmat *A_diag = precdata->A_diag;
  void **LU_diag = precdata->LU_diag;
  AMG_param *amgparam = precdata->amgparam;
  AMG_data **mgl = precdata->mgl;

  // data for the velocity block including eg part
  dvector **velocity_diag = precdata->diag;

  //-----------------------------------------
  // solver data for the pressure part
  //-----------------------------------------
  dvector *el_vol = precdata->el_vol;

  //-----------------------------------------
  // local variables
  //-----------------------------------------
  INT i;
  INT brow = A->brow;

  //-----------------------------------------
  // get different sizes
  //-----------------------------------------
  // get size of the velocity block without the eg part
  const INT Nu_wo_eg = A_diag[0].row;

  // get size of the veclosity block including the eg part
  const INT Nu = velocity_diag[0]->row;

  // get size of the pressure block
  const INT Np = el_vol->row;

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
  dvector ru_wo_eg, ru, rp, zu_wo_eg, zu, zp;

  ru_wo_eg.row = Nu_wo_eg; zu_wo_eg.row = Nu_wo_eg;
  ru.row = Nu; zu.row = Nu;
  rp.row = Np; zp.row = Np;

  ru_wo_eg.val = r; ru.val = r, rp.val = &(r[Nu]);
  zu_wo_eg.val = z; zu.val = z, zp.val = &(z[Nu]);

  //-----------------------------------------
  // main part of the preconditioner
  //-----------------------------------------
  // Preconditioning the velocity block without the eg part
  // use direct solver
  solve_UMF(&(A_diag[0]), &ru_wo_eg, &zu_wo_eg, LU_diag[0], 0);

  /*
  // use AMG solver
  mgl[0]->b.row=Nu_wo_eg; array_cp(Nu_wo_eg, ru_wo_eg.val, mgl[0]->b.val); // residual is an input
  mgl[0]->x.row=Nu_wo_eg; dvec_set(Nu_wo_eg, &mgl[0]->x,0.0);

  for(i=0;i<amgparam->maxit;++i) mgcycle(mgl[0], amgparam);
  array_cp(Nu_wo_eg, mgl[0]->x.val, zu_wo_eg.val);
  */

  // update residual
  array_cp(N, tempr->val, r);
  bdcsr_aAxpy(-1.0, A, z, r);

  // Preconditioning the velocity block including the eg part
  for(i=0; i<Nu; i++){
    zu.val[i] = zu.val[i] + ru.val[i]/velocity_diag[0]->val[i];
  }


  // update residual
  array_cp(N, tempr->val, r);
  bdcsr_aAxpy(-1.0, A, z, r);

  // Preconditioning the pressure block
  // Diagonal matrix for P0
  for(i=0;i<Np;i++){
    zp.val[i] = rp.val[i]/el_vol->val[i];
  }

  // restore r
  array_cp(N, tempr->val, r);

}

/* End of Program */
/*********************************************************************************************/
