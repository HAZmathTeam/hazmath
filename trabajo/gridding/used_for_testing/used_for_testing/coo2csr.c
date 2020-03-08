/* /\*---------------------------------------------------------------------*\/ */
/* void coo2csr(INT nrow,INT ncol,INT nnz,					\ */
/* 	     INT *row_idx,INT *col_idx, void *aval,			\ */
/* 	     INT *ia,INT *ja, void *bval,				\ */
/* 	     size_t elsize) */
/* { */
/*   // converts (i,j,value) to (ia,ja,bval) for matrices matrices whose */
/*   // elements have elsize bytes */
/*   // */
/*   // get dimensions of the matrix */
/*   // */
/*   const INT m=nrow, n=ncol; */
/*   INT i, iind, jind;   */
/*   INT *ind = (INT *) calloc(m+1,sizeof(INT)); */
/*   // initialize */
/*   memset(ind, 0, sizeof(INT)*(m+1));     */
/*   // count number of nonzeros in each row */
/*   for (i=0; i<nnz; ++i) ind[row_idx[i]+1]++;     */
/*   // set row pointer */
/*   ia[0] = 0; */
/*   for (i=1; i<=m; ++i) { */
/*     ia[i] = ia[i-1]+ind[i]; */
/*     ind[i] = ia[i]; */
/*   }     */
/*   // set column index and values */
/*   for (i=0; i<nnz; ++i) { */
/*     iind = row_idx[i]; */
/*     jind = ind[iind]; */
/*     ja[jind] = col_idx[i]; */
/*     memcpy((bval+jind*elsize),(aval+i*elsize),elsize); */
/*     ind[iind] = ++jind; */
/*   }     */
/*   if (ind) free(ind);     */
/*   return; */
/* } */
/*---------------------------------------------------------------------*/
