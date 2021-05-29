//-----------------------------------------------------------
// subroutines for preconditioner
//-----------------------------------------------------------
/*********************************************************************************/
static dvector *get_diag_bdcsr(block_dCSRmat *Ab, const REAL omega)
{
/*!
 * \fn dvector *get_diag_bdcsr(block_dCSRmat *Ab)
 *
 * \brief   extracts the diagonal of a block matrix in a dvector 
 *
 * \param   
 *
 * \note    
 *
 */
  // 
  INT nzd,i,j,k,ii,nn,iaa,iab,iblk,brow=Ab->brow,bcol=Ab->bcol;
  //  dvector **d_a=malloc(brow*sizeof(dvector *));
  //  for(iblk=0;iblk<brow;iblk++){
  //    j=Ab->blocks[iblk*bcol+iblk]->row;
  //    d_a[iblk]=dvec_create_p(j);
  //}
    //    d_a[iblk]=dvec_create_p(j);
  nn=0;// here: total number of rows 
  for(iblk=0;iblk<brow;iblk++){
    nn+=Ab->blocks[iblk*bcol+iblk]->row;
  }
  dvector *d_a = dvec_create_p(nn);
  nzd=0;
  for(iblk=0;iblk<brow;iblk++){
    ii=iblk*bcol+iblk;
    for(i=0;i<Ab->blocks[ii]->row;i++){
      iaa=Ab->blocks[ii]->IA[i];
      iab=Ab->blocks[ii]->IA[i+1];
      for(k=iaa;k<iab;k++){
	j=Ab->blocks[ii]->JA[k];
	if(i!=j) continue;
	d_a->val[nzd]=Ab->blocks[ii]->val[k];//aii
	nzd++;
      }
    }
  }
  for(i=0;i<nn;i++) d_a->val[i]*=omega;
  fprintf(stdout,"\nentries on the diagonal: %d (%d); scalling: omega=%.5e\n",nzd,nn,omega);  
  return d_a;
}
/*********************************************************************************/
static dCSRmat *get_diag_blocks(block_dCSRmat *Ab, const INT n10, const INT n20)
/*********************************************************************************/
/*!
 * \fn precond_block_data *get_pblock_data(block_dCSRmat *Ab, INT n2)
 *
 * \brief   get the diagonal blocks Ab(n1:n2,n1:n2) and store them in a dCSR matrix;
 *
 * \param   Ab   Pointer to a block_dCSRmat matrix
 *
 * \note    Memory space for the dCSRmat matrix is allocated inside this function! -- Xiaozhe Hu
 * \note    modeled on bdcsr_2_dcsr from utilities/format.c -- Ludmil
 *
 */
{
  // local variables
  INT m=0,n=0,nnz=0;
  const INT mb=Ab->brow, nb=Ab->bcol, n_blocks=mb*nb;
  dCSRmat **blockptr=Ab->blocks, *blockptrij, *A;
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
  A = dcsr_create_p(m,n,nnz);
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
  return A;
}
/**************************************************************************************/
static precond_block_data *get_precond_block_data(block_dCSRmat *Ab)
{
  precond_block_data *pblk=(precond_block_data *)malloc(1*sizeof(precond_block_data));
  //   block_dCSRmat *
  precond_block_data_null(pblk);
  pblk->Abcsr=Ab;
  INT brow=Ab->brow, bcol=Ab->bcol;
  INT nblk,iblk,n1,n2,j,k,l,m,n,iaa,iab;
  // blocks outside of [n1,n2) are ignored
  n1=0;
  n2=brow-1;  
  pblk->A_diag=get_diag_blocks(Ab,n1,n2);
  //dvector ** the diagonal of the stiffness matrix:
  pblk->diag=malloc(sizeof(dvector *));       
  pblk->diag[0] = get_diag_bdcsr(Ab,1.0);
  /*-------------------------------------------------*/  
  /* FILE *fptmp; */
  /* fptmp=fopen("a.dat","w"); */
  /* bdcsr_print_matlab(fptmp,Ab); */
  /* fclose(fptmp); */
  /* //////////////////////////// */
  /* fptmp=fopen("a11.dat","w"); */
  /* csr_print_matlab(fptmp,pblk->A_diag); */
  /* fclose(fptmp); */
  /* fptmp=fopen("d.dat","w"); */
  /* dvector_print(fptmp,pblk->diag[0]); */
  /* fclose(fptmp); */
  /*----------------------------------------------------------------*/
  pblk->LU_diag=malloc(sizeof(void *));
  pblk->LU_diag[0]=factorize_UMF(pblk->A_diag,0);
  return pblk;
}
/**************************************************************************************/
static void free_precond_block_data(precond_block_data *pblk)
{
  free(pblk->A_diag);
  free(pblk->diag[0]);
  free(pblk->diag);
  free(pblk->LU_diag[0]);  
  free(pblk->LU_diag);
  free(pblk);
  return;
}
