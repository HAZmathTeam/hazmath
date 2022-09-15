/*! \file src/utilities/format.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 02/17/16.
 *  Copyright 2016__HAZMATH__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/28/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/28/2016 & 08/28/2021
 *
 */

#include "hazmath.h"

/***********************************************************************************************/
/*!
 * \fn dCSRmat bdcsr_2_dcsr (block_dCSRmat *Ab)
 *
 * \brief   Transform a block_dCSRmat matrix to a dCSRmat matrix
 *
 * \param   Ab   Pointer to a block_dCSRmat matrix
 *
 * \return  A    dCSRmat matrix if succeed, NULL if fail
 *
 * \note    Memory space for the dCSRmat matrix is allocated inside this functions! -- Xiaozhe Hu
 *
 */
dCSRmat bdcsr_2_dcsr (block_dCSRmat *Ab)
{
    // local variables
    INT m=0,n=0,nnz=0;
    const INT mb=Ab->brow, nb=Ab->bcol, n_blocks=mb*nb;
    dCSRmat **blockptr=Ab->blocks, *blockptrij, A;
    INT i,j,ij,ir,i1,length,ilength,start,irmrow,irmrowp1;
    INT *row, *col;

    // flag for errors
    SHORT status = SUCCESS;

    row = (INT *)calloc(mb+1,sizeof(INT));
    col = (INT *)calloc(nb+1,sizeof(INT));

    // get the size of A
    row[0]=0; col[0]=0;

    // count number of rows
    for (i=0;i<mb;++i) {

        status = ERROR_BLKMAT_ZERO;

        for (j=0; j<nb; ++j){
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
    for (i=0;i<nb;++i) {

        status = ERROR_BLKMAT_ZERO;

        for (j=0;j<mb;++j){
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
    for (i=0;i<n_blocks;++i) {
        if (blockptr[i]) {
            nnz+=blockptr[i]->nnz;
        }
    }

    // memory space allocation
    A = dcsr_create(m,n,nnz);

    // set dCSRmat for A
    A.IA[0]=0;
    for (i=0;i<mb;++i) {

        for (ir=row[i];ir<row[i+1];ir++) {

            for (length=j=0;j<nb;++j) {

                ij=i*nb+j;
                blockptrij=blockptr[ij];

                if (blockptrij && blockptrij->nnz>0) {

                    start=A.IA[ir]+length;
                    irmrow=ir-row[i];irmrowp1=irmrow+1;
                    ilength=blockptrij->IA[irmrowp1]-blockptrij->IA[irmrow];

                    if (ilength>0) {
                        memcpy(&(A.val[start]),&(blockptrij->val[blockptrij->IA[irmrow]]),ilength*sizeof(REAL));
                        memcpy(&(A.JA[start]), &(blockptrij->JA[blockptrij->IA[irmrow]]), ilength*sizeof(INT));
                        for (i1=0;i1<ilength;i1++) A.JA[start+i1]+=col[j];
                        length+=ilength;
                    }

                 }

            } // end for j

            A.IA[ir+1]=A.IA[ir]+length;

        } // end for ir

    } // end for i

    free(row);
    free(col);

    return(A);
}

/***********************************************************************************************/
/*!
 * \fn dCSRmat bdcsr_subblk_2_dcsr (block_dCSRmat *Ab)
 *
 * \brief   Transform a block_dCSRmat matrix to a dCSRmat matrix
 *
 * \param   Ab   Pointer to a block_dCSRmat matrix
 *
 * \return  A    dCSRmat matrix if succeed, NULL if fail
 *
 * \note    Memory space for the dCSRmat matrix is allocated inside this functions! -- Xiaozhe Hu
 *
 */
dCSRmat bdcsr_subblk_2_dcsr (block_dCSRmat *Ab, INT brow_start, INT brow_end, INT bcol_start, INT bcol_end)
{
    INT i, j;
    INT ii, jj;
    block_dCSRmat Asub;
    bdcsr_alloc_minimal( brow_end-brow_start+1, bcol_end-bcol_start+1, &Asub);

    for(j = bcol_start; j<bcol_end+1; j++){
      jj = j-bcol_start;
      for(i = brow_start; i<brow_end+1; i++){
        ii = i-brow_start;
        Asub.blocks[ii*Asub.bcol + jj] = Ab->blocks[i*Ab->bcol + j];
      }
    }

    dCSRmat A = bdcsr_2_dcsr( &Asub );
//    bdcsr_free_minimal( &Asub );
    return(A);
}

/***********************************************************************************************/
/*!
 * \fn dBSRmat bdcsr_2_dbsr (block_dCSRmat *Ab)
 *
 * \brief   Convert a block_dCSRmat matrix to a dCSRmat matrix
 *
 * \param   Ab   Pointer to a block_dCSRmat matrix
 *
 * \return  A    dBSRmat matrix if succeed, NULL if fail
 *
 *
 */
dBSRmat bdcsr_2_dbsr(block_dCSRmat *Ab)
{
    // local variables
    dCSRmat A0 = bdcsr_2_dcsr(Ab);

    INT i,j;
    INT nb = Ab->brow;

    // generate the ordering
    ivector order = ivec_create(A0.row);
    for (i=0; i<order.row/nb; i++){
        for (j=0; j<nb; j++){
            order.val[i*nb+j] = j*(order.row/nb)+i;
        }
    }

    // reorder dCOOmat
    dCSRmat A1 = dcsr_reorder(&A0, order.val);
    dcsr_free(&A0);

    // convert to dBSRmat
    dBSRmat B = dcsr_2_dbsr(&A1, nb);

    // clean
    dcsr_free(&A1);

    // return
    return B;

}

/***********************************************************************************************/
/*!
 * \fn block_dCSRmat dcsr_2_bdcsr (dCSRmat *A, INT bnum, INT *bsize)
 *
 * \brief
 *
 * \param A         Pointer to dCSRmat matrix
 * \param bnum      Block size of the block_dCSRmat matrix
 * \param bsize     Pointer to the size of each diagonal block
 *
 * \return Ab       block_dCSRmat matrix if succeed, NULL if fail
 *
 * \note Assume dCSRmat A has block form already and just need to explicitly form the blocks without reordering! -- Xiaozhe Hu
 * \note SUM(bsize) = A->row = A->col, i.e., size has to be consistent!!  -- Xiaozhe Hu
 *
 */
block_dCSRmat dcsr_2_bdcsr (dCSRmat *A,
                            INT bnum,
                            INT *bsize)
{
    // local variable
    INT i, j;
    SHORT status = SUCCESS;

    // allocate block dCSRmat
    block_dCSRmat Ab;
    Ab.brow = bnum; Ab.bcol = bnum;

    Ab.blocks = (dCSRmat **)calloc(bnum*bnum, sizeof(dCSRmat));
    for (i=0; i<bnum*bnum; i++) Ab.blocks[i] = (dCSRmat *) calloc(1, sizeof(dCSRmat));


    // allocate
    INT *idx_row = (INT *)calloc(A->row, sizeof(INT));
    INT *idx_col = (INT *)calloc(A->col, sizeof(INT));

    for (i=0; i<A->row; i++) idx_row[i] = i;  // generate row index
    for (i=0; i<A->col; i++) idx_col[i] = i;  // generate col index

    INT *idx_start = (INT *)calloc(bnum, sizeof(INT));
    idx_start[0] = 0;
    for (i=0; i<bnum-1; i++) idx_start[i+1] = idx_start[i] + bsize[i];

    for (i=0; i<bnum; i++){

        for (j=0; j<bnum; j++) {

            status = dcsr_getblk(A, idx_row+idx_start[i], idx_col+idx_start[j], bsize[i], bsize[j], Ab.blocks[i*bnum+j]);
            if (status < SUCCESS) check_error(status, __FUNCTION__);

        }

    }


    // free memory
    free(idx_row);
    free(idx_col);
    free(idx_start);

    return Ab;
}

/***********************************************************************************************/
/*!
 * \fn SHORT dcoo_2_dcsr (dCOOmat *A, dCSRmat *B)
 *
 * \brief Transform a dCOOmat matrix to a dCSRmat format.
 *
 * \param A   Pointer to dCOOmat matrix
 * \param B   Pointer to dCSRmat matrix
 *
 * \return    SUCCESS if successed; otherwise, error information.
 *
 */
SHORT dcoo_2_dcsr (dCOOmat *A,
                   dCSRmat *B)
{
    // get size
    const INT m=A->row, n=A->col, nnz=A->nnz;

    // allocate
    dcsr_alloc(m,n,nnz,B);

    // local variable
    INT *ia = B->IA;
    INT *ja = B->JA;
    REAL *Bval = B->val;
    INT *row_idx = A->rowind;
    INT *col_idx = A->colind;
    REAL *Aval = A->val;
    INT i, iind, jind;

    INT *ind = (INT *) calloc(m+1,sizeof(INT));

    // initialize
    memset(ind, 0, sizeof(INT)*(m+1));

    // count number of nonzeros in each row
    for (i=0; i<nnz; ++i) ind[row_idx[i]+1]++;

    // set row pointer
    ia[0] = 0;
    for (i=1; i<=m; ++i) {
        ia[i] = ia[i-1]+ind[i];
        ind[i] = ia[i];
    }

    // set column index and values
    for (i=0; i<nnz; ++i) {
        iind = row_idx[i];
        jind = ind[iind];
        ja[jind] = col_idx[i];
        Bval[jind] = Aval[i];
        ind[iind] = ++jind;
    }

    if (ind) free(ind);

    return SUCCESS;
}

/***********************************************************************************************/
/*!
 * \fn SHORT dcsr_2_dcoo (dCSRmat *A, dCOOmat *B)
 *
 * \brief Transform a dCSRmat matrix to a dCOOmat format.
 *
 * \param A   Pointer to dCSRmat matrix
 * \param B   Pointer to dCOOmat matrix
 *
 * \return    SUCCESS if successed; otherwise, error information.
 *
 */
SHORT dcsr_2_dcoo (dCSRmat *A,
		   dCOOmat *B)
{
  //  dcoo_alloc(m,n,nnz,B);
  B->rowind=calloc(A->nnz,sizeof(INT));
  B->colind=calloc(A->nnz,sizeof(INT));
  B->val=calloc(A->nnz,sizeof(REAL));
  INT i,ij;
  B->row=A->row;
  B->col=A->col;
  B->nnz=A->nnz;
  for(i=0;i<A->row;++i){
    for(ij=A->IA[i];ij<A->IA[i+1];ij++){
      B->rowind[ij]=i;
      B->colind[ij]=A->JA[ij];
      B->val[ij]=A->val[ij];
    }
  }
  //  memcpy(B->colind,A->JA,nnz*sizeof(INT));
  //  memcpy(B->val,A->val,nnz*sizeof(REAL));
  return 0;
}
/***********************************************************************************************/
/*!
 * \fn SHORT dcsr_unique (dCSRmat *A)
 *
 * \brief Removes repetitions from column indices in a dCSRmat A.
 *
 * \param A   Pointer to dCSRmat matrix
 *
 * \return SUCCESS if successful;
 *
 * \note the input matrix A is overwritten with the new CSR matrix
 *       which has no repetitions in the column indices and the value
 *       corresponding to a column index is the sum of all column
 *       indices in a row that have the same column index.
 *
 *  Ludmil 20210530.
 */
SHORT dcsr_unique (dCSRmat *A)
{
    // get size
    INT m=A->row, n=A->col, nnz=A->nnz;
    // copy pointers for easier reference
    INT *ia = A->IA;
    INT *ja = A->JA;
    REAL *a = A->val;
    INT i, ij,j,ih,iai,k;
    SHORT norepeat;
    INT maxdeg=ia[1]-ia[0];

    INT *ind = (INT *) calloc(n,sizeof(INT));
    for (i=0; i<n; ++i) ind[i]=-1;
    // clean up because there might be some repeated indices
    // compute max degree of all vertices (for memory allocation):
    for (i=1;i<m;++i){
      ih=ia[i+1]-ia[i];
      if(maxdeg<ih) maxdeg=ih;
    }
    REAL *atmp=calloc(maxdeg,sizeof(REAL));
    INT *jatmp=calloc(maxdeg,sizeof(INT));
    nnz=0;
    for (i=0;i<m;++i){
      // loop over each row. first find the length of the row:
      ih=ia[i+1]-ia[i];
      // copy the indices in tmp arrays.
      memcpy(jatmp,(ja+ia[i]),ih*sizeof(INT));
      memcpy(atmp,(a+ia[i]),ih*sizeof(REAL));
      norepeat=1;
      for(ij=0;ij<ih;++ij){
	j=jatmp[ij];
	if(ind[j]<0){
	  ind[j]=ij;
	} else {
	  norepeat=0; // we have a repeated index.
	  atmp[ind[j]]+=atmp[ij];
	  jatmp[ij]=-abs(jatmp[ij]+1);
	}
      }
      for(k=0;k<ih;++k)
      for(ij=0;ij<ih;++ij){
	j=jatmp[ij];
	if(j<0) continue;// if j is negative, this has repeated somewhere. do nothing
	ind[j]=-1;// make, for all all visited j, ind[j]=-1;
      }
      if(norepeat) continue; // do nothing if no indices repeat.
      // put everything back, but now we have negative column indices
      // on the repeated column indices and we have accumulated the
      // values in the first position of j on every row.
      memcpy(&ja[ia[i]],jatmp,ih*sizeof(INT));
      memcpy(&a[ia[i]],atmp,ih*sizeof(REAL));
    }
    if (ind) free(ind);
    if(atmp) free(atmp);
    if(jatmp) free(jatmp);
    // run over the matrix and remove all negative column indices.
    iai=ia[0];
    nnz=0;
    for (i=0;i<m;++i){
      for(ij=iai;ij<ia[i+1];++ij){
	j=ja[ij];
	if(j<0) continue;
	ja[nnz]=ja[ij];
	a[nnz]=a[ij];
	++nnz;
      }
      iai=ia[i+1];
      ia[i+1]=nnz;
    }
    A->nnz=nnz;
    A->val=realloc(A->val,A->nnz*sizeof(REAL));
    A->JA=realloc(A->JA,A->nnz*sizeof(INT));
    return SUCCESS;
}

/***********************************************************************************************/
/*!
 * \fn coordinates* array_2_coord (REAL *xyz, INT ndof, INT dim)
 *
 * \brief Transform an array of coordinates from [x0,y0,z0,x1,y1,z1,...] to a coordinates structure
 *
 * \param xyz  Array with coordinate values
 * \param ndof Number of degrees of freedom
 * \param dim  Dimension
 *
 * \return cv Coordinate struct
 *
 */
coordinates* array_2_coord ( REAL* xyz, INT ndof, INT dim)
{
  // Local Variables
  INT dof;//,i,j;
  // Allocate
  coordinates* cv = allocatecoords(ndof, dim);
  // Fill
  switch (dim)
  {
  case 1:
    for( dof=0; dof<ndof; dof++){
      cv->x[dof] = xyz[dof*dim+0];
    }
    break;
  case 2:
    for( dof=0; dof<ndof; dof++){
      cv->x[dof] = xyz[dof*dim+0];
      cv->y[dof] = xyz[dof*dim+1];
    }
    break;
  case 3:
    for( dof=0; dof<ndof; dof++){
      cv->x[dof] = xyz[dof*dim+0];
      cv->y[dof] = xyz[dof*dim+1];
      cv->z[dof] = xyz[dof*dim+2];
    }
    break;
  default:
    check_error( ERROR_DIM, __FUNCTION__);
  }
  return cv;
}


/*!
 * \fn dCSRmat dbsr_2_dcsr (const dBSRmat *B)
 *
 * \brief Transfer a 'dBSRmat' type matrix into a dCSRmat.
 *
 * \param B   Pointer to dBSRmat matrix
 *
 * \return    dCSRmat matrix
 *
 *
 * \note Works for general nb (Xiaozhe)
 */
dCSRmat dbsr_2_dcsr (const dBSRmat *B)
{
    dCSRmat A;

    /* members of B */
    INT     ROW = B->ROW;
    INT     COL = B->COL;
    INT     NNZ = B->NNZ;
    INT     nb  = B->nb;
    INT    *IA  = B->IA;
    INT    *JA  = B->JA;
    REAL   *val = B->val;

    INT     storage_manner = B->storage_manner;

    INT jump = nb*nb;
    INT rowA = ROW*nb;
    INT colA = COL*nb;
    INT nzA  = NNZ*jump;

    INT     *ia = NULL;
    INT     *ja = NULL;
    REAL    *a  = NULL;

    INT i,j,k;
    INT mr,mc;
    INT rowstart0,rowstart,colstart0,colstart;
    INT colblock,nzperrow;

    REAL  *vp = NULL;
    REAL  *ap = NULL;
    INT  *jap = NULL;

    //--------------------------------------------------------
    // Create a CSR Matrix
    //--------------------------------------------------------
    A  = dcsr_create(rowA, colA, nzA);
    ia = A.IA;
    ja = A.JA;
    a  = A.val;

    //--------------------------------------------------------------------------
    // Compute the number of nonzeros per row, and after this loop,
    // ia[i],i=1:rowA, will be the number of nonzeros of the (i-1)-th row.
    //--------------------------------------------------------------------------
    for (i = 0; i < ROW; ++i)
    {
        rowstart = i*nb + 1;
        colblock = IA[i+1] - IA[i];
        nzperrow = colblock*nb;
        for (j = 0; j < nb; ++j)
        {
            ia[rowstart+j] = nzperrow;
        }
    }

    //-----------------------------------------------------
    // Generate the real 'ia' for CSR of A
    //-----------------------------------------------------
    ia[0] = 0;
    for (i = 1; i <= rowA; ++i)
    {
        ia[i] += ia[i-1];
    }

    //-----------------------------------------------------
    // Generate 'ja' and 'a' for CSR of A
    //-----------------------------------------------------
    switch (storage_manner)
    {
        case 0: // each non-zero block elements are stored in row-major order
        {
            for (i = 0; i < ROW; ++i)
            {
                for (k = IA[i]; k < IA[i+1]; ++k)
                {
                    j = JA[k];
                    rowstart = i*nb;
                    colstart = j*nb;
                    vp = &val[k*jump];
                    for (mr = 0; mr < nb; mr ++)
                    {
                        ap  = &a[ia[rowstart]];
                        jap = &ja[ia[rowstart]];
                        for (mc = 0; mc < nb; mc ++)
                        {
                            *ap = *vp;
                            *jap = colstart + mc;
                            vp ++; ap ++; jap ++;
                        }
                        ia[rowstart] += nb;
                        rowstart ++;
                    }
                }
            }
        }
            break;

        case 1: // each non-zero block elements are stored in column-major order
        {
            for (i = 0; i < ROW; ++i)
            {
                for (k = IA[i]; k < IA[i+1]; ++k)
                {
                    j = JA[k];
                    rowstart0 = i*nb;
                    colstart0 = j*nb;
                    vp = &val[k*jump];
                    for (mc = 0; mc < nb; mc ++)
                    {
                        rowstart = rowstart0;
                        colstart = colstart0 + mc;
                        for (mr = 0; mr < nb; mr ++)
                        {
                            a[ia[rowstart]] = *vp;
                            ja[ia[rowstart]] = colstart;
                            vp ++; ia[rowstart]++; rowstart++;
                        }
                    }
                }
            }
        }
            break;
    }

    //-----------------------------------------------------
    // Map back the real 'ia' for CSR of A
    //-----------------------------------------------------
    for (i = rowA; i > 0; i --) {
        ia[i] = ia[i-1];
    }
    ia[0] = 0;

    return (A);
}

/*!
 * \fn dBSRmat dcsr_2_dbsr ( const dCSRmat *A, const INT nb )
 *
 * \brief Transfer a dCSRmat type matrix into a dBSRmat.
 *
 * \param A   Pointer to the dCSRmat type matrix
 * \param nb  size of each block
 *
 * \return    dBSRmat matrix
 *
 */
dBSRmat dcsr_2_dbsr(const dCSRmat  *A,
                    const INT       nb)
{
    INT i, j, k, ii, jj, kk, l, mod, nnz;
    INT row   = A->row/nb;
    INT col   = A->col/nb;
    INT nb2   = nb*nb;
    INT *IA   = A->IA;
    INT *JA   = A->JA;
    REAL *val = A->val;

    dBSRmat B;	// Safe-guard check
	INT *col_flag, *ia, *ja;
	REAL *bval;

    if ((A->row)%nb!=0) {
        printf("### ERROR: A.row=%lld is not a multiplication of nb=%lld!\n",
                 (long long )A->row,   (long long )nb);
        check_error(ERROR_MAT_SIZE, __FUNCTION__);
    }

    if ((A->col)%nb!=0) {
        printf("### ERROR: A.col=%lld is not a multiplication of nb=%lld!\n",
                 (long long )A->col,   (long long )nb);
        check_error(ERROR_MAT_SIZE, __FUNCTION__);
    }

    B.ROW = row;
    B.COL = col;
    B.nb  = nb;
    B.storage_manner = 0;

    // allocate memory for B
	col_flag = (INT *)calloc(col, sizeof(INT));
    ia = (INT *)calloc(row+1, sizeof(INT));

    iarray_set(col, col_flag, -1);

    // Get ia for BSR format
    nnz = 0;
	for (i=0; i<row; ++i) {
        ii = nb*i;
        for (j=0; j<nb; ++j) {
            jj = ii+j;
            for (k=IA[jj]; k<IA[jj+1]; ++k) {
                kk = JA[k]/nb;
                if (col_flag[kk]!=0) {
                    col_flag[kk] = 0;
                    //ja[nnz] = kk;
                    nnz ++;
                }
			}
		}
        ia[i+1] = nnz;
        iarray_set(col, col_flag, -1);
	}

    // set NNZ
    B.NNZ = nnz;

    // allocate ja and bval
    ja = (INT*)calloc(nnz, sizeof(INT));
    bval = (REAL*)calloc(nnz*nb2, sizeof(REAL));

    // Get ja for BSR format
    nnz = 0;
    for (i=0; i<row; ++i) {
        ii = nb*i;
        for(j=0; j<nb; ++j) {
            jj = ii+j;
            for(k=IA[jj]; k<IA[jj+1]; ++k) {
                kk = JA[k]/nb;
                if (col_flag[kk]!=0) {
                    col_flag[kk] = 0;
                    ja[nnz] = kk;
                    nnz ++;
                }
			}
		}
        ia[i+1] = nnz;
        iarray_set(col, col_flag, -1);
	}

    // Get non-zeros of BSR
	for (i=0; i<row; ++i) {
		ii = nb*i;
        for (j=0; j<nb; ++j) {
			jj = ii+j;
			for (k=IA[jj]; k<IA[jj+1]; ++k) {
				for (l=ia[i]; l<ia[i+1]; ++l) {
					if (JA[k]/nb ==ja[l]) {
                        mod = JA[k]%nb;
                        bval[l*nb2+j*nb+mod] = val[k];
                        break;
                    }
				}
			}
		}
	}

    B.IA = ia;
    B.JA = ja;
    B.val = bval;

    free(col_flag); col_flag = NULL;

    return B;
}

/***********************************************************************************************/
/*!
 * \fn dCSRmat *dcoo_2_dcsr_p(dCOOmat *A)
 *
 * \brief Transform a dCOOmat matrix to a dCSRmat format.
 *
 * \param A   Pointer to dCOOmat matrix
 *
 * \return B   Pointer to dCSRmat matrix
 *
 *
 */
dCSRmat *dcoo_2_dcsr_p (dCOOmat *A)
{
    // get size
    const INT m=A->row, n=A->col, nnz=A->nnz;
    // allocate
    dCSRmat *B=dcsr_create_p(m,n,nnz);
    // local variable
    INT *ia = B->IA;
    INT *ja = B->JA;
    REAL *Bval = B->val;
    INT *row_idx = A->rowind;
    INT *col_idx = A->colind;
    REAL *Aval = A->val;
    INT i, iind, jind;

    INT *ind = (INT *) calloc(m+1,sizeof(INT));

    // initialize
    memset(ind, 0, sizeof(INT)*(m+1));

    // count number of nonzeros in each row
    for (i=0; i<nnz; ++i) ind[row_idx[i]+1]++;

    // set row pointer
    ia[0] = 0;
    for (i=1; i<=m; ++i) {
        ia[i] = ia[i-1]+ind[i];
        ind[i] = ia[i];
    }

    // set column index and values
    for (i=0; i<nnz; ++i) {
        iind = row_idx[i];
        jind = ind[iind];
        ja[jind] = col_idx[i];
        Bval[jind] = Aval[i];
        ind[iind] = ++jind;
    }

    if (ind) free(ind);

    return B;
}
/********************************  END  ********************************************************/
